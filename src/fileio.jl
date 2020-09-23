module Files

using HDF5
using FileIO
using Requires
using SparseArrays
using SparseArrays: getcolptr
import UnixMmap

export read_obsmat, write_obsmat

"""
    read_obsmat(filename::String; keywords...)

Read a sparse observing matrix and the corresponding pixelization descriptors from
`filename`. The return value is a `NamedTuple` with the following fields:

- `R`: the sparse observing matrix ``R``.
- `pixr`: either `missing` or a description of the "right-hand side" pixelization format
  of the observing matrix — i.e. the pixelization of a map vector ``v`` for which the
  matrix-vector multiplication ``R v`` is defined.
- `pixl`: either `missing` or a description of the "left-hand side" pixelication format
  of the observing matrix — i.e. the pixelization of a map vector ``w`` which results
  from the matrix-vector multiplication ``w = R v``.

The pixel descriptions `pixr` and `pixl` are loaded from the datasets named by the keyword
arguments of the same name, if available. If the named dataset exists, then it is loaded
and returned. (The data type is format- and situation-specific.) It is not an error for the
named dataset to not exist, and the field will be returned with the value `missing`.

**See also:** [`write_obsmat`](@ref)

# Extended help

## Keywords

The following keywords identify the specific data structures or fields within the
given data file that correspond to the observing matrix and pixelization data structures
or fields, and they are available on all backends.
Additional keywords may be supported on a backend-specific basis.

- `name`: Name of the observing matrix ``R``. Defaults to `"R"`.
- `pixr`: Name of the data structure or field describing the right-hand side
  pixelization format. Defaults to `"pixels_right"`.
- `pixl`: Name of the data structure or field describing the "left-hand side"
  pixelization format. Defaults to `"pixels_right"`.

In all cases, the field names may be `nothing` to indicate the corresponding data should
not be loaded, and the corresponding field in the returned named tuple will have the
value `missing`. It is an error for `name` to point to a non-existent field, whereas the
pixel descriptions will silently ignore an invalid name and just return `missing` instead.

## Backends

The `HDF5.jl` storage backend is always loaded with `CMB.jl` and is considered the native
storage format.
See [`write_obsmat`](@ref) for writing a native HDF5 file to disk.

Importing observing matrices from the following additional data formats is supported
via `Requires.jl`, which requires the user to first load the extra backend of choice.

- `JLD` and `JLD2`-flavored HDF5 files with `JLD.jl` and `JLD2.jl`, respectively.
- MATLAB v5, v6, v7, and v7.3 save files with `MAT.jl`.
- `scipy.sparse` CSC and CSR matrices saved to HDF5 files with `h5sparse`. This case
  is supported without needing to load any extra packages.

The pixelization descriptions are imported as data format specific types. For instance,
all formats support loading strings and simple numerical scalars or arrays. Additionally,
JLD and JLD2 formats can load named datasets as arbitrary Julia types, MATLAB structs
are deserialized as `Dict`s, and named HDF5 groups are read as (nested) `Dict`s.
"""
function read_obsmat(filename::String; kws...)
    file = query(filename)
    return read_obsmat(file; kws...)
end

include("fileio_hdf5.jl")

function __init__()
    if Sys.islinux()
        READ_OBSMAT_MMAP_FLAGS[] |= UnixMmap.MAP_POPULATE
    end

    @require JLD = "4138dd39-2aa7-5051-a626-17a0bb65d9c8" begin
        using .JLD
        """
            read_obsmat(file::FileIO.File{format"JLD"}; keywords...)

        Reads an observing matrix and pixelization descriptors from a JLD-formatted HDF5 file.
        This function is conditionally included via `Requires.jl` and requires the user to
        `import JLD` or `using JLD` first.
        """
        function read_obsmat(file::File{format"JLD"};
                             name::Union{String,Nothing} = "R",
                             pixr::Union{String,Nothing} = "pixels_right",
                             pixl::Union{String,Nothing} = "pixels_left")
            hid = JLD.jldopen(FileIO.filename(file))
            try
                R = name !== nothing ? read(hid, name) : missing
                pixelsr = (pixr !== nothing && exists(hid, pixr)) ?
                        read(hid, pixr) : missing
                pixelsl = (pixr !== nothing && exists(hid, pixl)) ?
                        read(hid, pixl) : missing
                return (; R = R, pixr = pixelsr, pixl = pixelsl)
            finally
                close(hid)
            end
        end
    end
    @require JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819" begin
        using .JLD2
        """
            read_obsmat(file::FileIO.File{format"JLD2"}; keywords...)

        Reads an observing matrix and pixelization descriptors from a JLD2-formatted HDF5 file.
        This function is conditionally included via `Requires.jl` and requires the user to
        `import JLD2` or `using JLD2` first.
        """
        function read_obsmat(file::File{format"JLD2"};
                             name::Union{String,Nothing} = "R",
                             pixr::Union{String,Nothing} = "pixels_right",
                             pixl::Union{String,Nothing} = "pixels_left")
            hid = JLD2.jldopen(FileIO.filename(file))
            try
                R = name !== nothing ? read(hid, name) : missing
                pixelsr = (pixr !== nothing && haskey(hid, pixr)) ?
                        read(hid, pixr) : missing
                pixelsl = (pixr !== nothing && haskey(hid, pixl)) ?
                        read(hid, pixl) : missing
                return (; R = R, pixr = pixelsr, pixl = pixelsl)
            finally
                close(hid)
            end
        end
    end

    @require MAT = "23992714-dd62-5051-b70f-ba57cb901cac" begin
        # FileIO doesn't define a MAT format, so do it ourselves
        add_format(format"MAT", FileIO.detecthdf5, [".mat"], [:MAT])
        using .MAT
        """
            read_obsmat(file::FileIO.File{format"MAT"}; keywords...)

        Reads an observing matrix and pixelization descriptors from a MATLAB save file in
        v5, v6, v7, or v7.3 format.
        This function is conditionally included via `Requires.jl` and requires the user to
        `import MAT` or `using MAT` first.
        """
        function read_obsmat(file::File{format"MAT"};
                             name::Union{String,Nothing} = "R",
                             pixr::Union{String,Nothing} = "pixels_right",
                             pixl::Union{String,Nothing} = "pixels_left")
            hid = matopen(FileIO.filename(file))
            try
                R = name !== nothing ? read(hid, name) : missing
                pixelsr = (pixr !== nothing && exists(hid, pixr)) ?
                        read(hid, pixr) : missing
                pixelsl = (pixr !== nothing && exists(hid, pixl)) ?
                        read(hid, pixl) : missing
                return (; R = R, pixr = pixelsr, pixl = pixelsl)
            finally
                close(hid)
            end
        end
    end
end

end # module Files
