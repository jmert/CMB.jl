module Files

using HDF5
using FileIO
using Requires
using SparseArrays
using SparseArrays: getcolptr
import UnixMmap

export read_obsmat, write_obsmat

"""
    R, metadata = read_obsmat(filename::String; keywords...)

Read a sparse observing matrix `R` and the corresponding metadata `metadata` from
`filename`. The metadata data will include the following fields:

- `fields`: a description of the Stokes fields for which the observing matrix applies —
  e.g. the value "QU" signals that `R` is the block matrix `[R_QQ R_QU; R_UQ R_UU]`.
- `pixels_right`: a description of the "right-hand side" pixelization format of the
  observing matrix — i.e. the pixelization of a map vector ``v`` for which the matrix-vector
  multiplication ``R v`` is defined.
- `pixels_left`: a description of the "left-hand side" pixelication format of the observing
  matrix — i.e. the pixelization of a map vector ``w`` which results from the matrix-vector
  multiplication ``w = R v``.

The metadata fields are loaded from datasets named by keyword arguments of the same name.
Any keyword set to `nothing` indicates that the corresponding data should not be loaded.
The data types are format- and situation-specific. It is not an error for the named dataset
to not exist, and if the dataset does not exist, the field of the metadata named tuple will
be filled with the value `missing`.

The observing matrix dataset within the data file is given by the keyword `name` and
defaults to `"R"`. It is an error for the named dataset to not exist (if not `nothing`).

**See also:** [`write_obsmat`](@ref)

# Extended help

## Backends

The `HDF5.jl` storage backend is always loaded with `CMB.jl` and is considered the native
storage format.
See [`write_obsmat`](@ref) for writing a native HDF5 file to disk.

Importing observing matrices from the following additional data formats is supported
via `Requires.jl`, which requires the user to first load the extra backend of choice.

- `JLD2`-flavored HDF5 files with `JLD2.jl`.
- MATLAB v5, v6, v7, and v7.3 save files with `MAT.jl`.
- `scipy.sparse` CSC and CSR matrices saved to HDF5 files with `h5sparse`. This case
  is supported without needing to load any extra packages.

The pixelization descriptions are imported as data format specific types. For instance,
all formats support loading strings and simple numerical scalars or arrays. Additionally,
the JLD2 format can load named datasets as arbitrary Julia types, MATLAB structs
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
                             fields::Union{String,Nothing} = "fields",
                             pixels_right::Union{String,Nothing} = "pixels_right",
                             pixels_left::Union{String,Nothing} = "pixels_left")
            hid = JLD2.jldopen(FileIO.filename(file))
            @inline function _missing_read(name)
                name === nothing && return missing
                !haskey(hid, name) && return missing
                return hid[name]
            end
            try
                R = name !== nothing ? hid[name] : missing
                metadata = (;
                            fields = _missing_read(fields),
                            pixels_right = _missing_read(pixels_right),
                            pixels_left = _missing_read(pixels_left)
                           )
                return R, metadata
            finally
                close(hid)
            end
        end
    end

    @require MAT = "23992714-dd62-5051-b70f-ba57cb901cac" begin
        # FileIO doesn't define a MAT format, so do it ourselves
        add_format(format"MAT", FileIO.detecthdf5, [".mat"], [:MAT])
        using .MAT

        # Because Matlab doesn't have vectors, we're forced to try to infer a vector
        # from a matrix based on its size. ಠ_ಠ Gahh!!!
        _mat2vec(x::Any) = x
        function _mat2vec(x::AbstractDict)
            for (k, v) in x
                x[k] = _mat2vec(v)
            end
            return x
        end
        function _mat2vec(x::AbstractMatrix)
            if size(x, 2) == 1
                return vec(x)
            else
                return x
            end
        end

        """
            read_obsmat(file::FileIO.File{format"MAT"}; keywords...)

        Reads an observing matrix and pixelization descriptors from a MATLAB save file in
        v5, v6, v7, or v7.3 format.
        This function is conditionally included via `Requires.jl` and requires the user to
        `import MAT` or `using MAT` first.
        """
        function read_obsmat(file::File{format"MAT"};
                             name::Union{String,Nothing} = "R",
                             fields::Union{String,Nothing} = "fields",
                             pixels_right::Union{String,Nothing} = "pixels_right",
                             pixels_left::Union{String,Nothing} = "pixels_left")
            hid = matopen(FileIO.filename(file))
            @inline function _missing_read(name)
                name === nothing && return missing
                !exists(hid, name) && return missing
                return read(hid, name)
            end
            try
                if name === nothing
                    R = missing
                else
                    !exists(hid, name) && error("Error reading /", name)
                    R = read(hid, name)
                end
                fields = _missing_read(fields)
                fields = fields isa Char ? join(fields) : fields
                metadata = (;
                            fields = fields,
                            pixels_right = _mat2vec(_missing_read(pixels_right)),
                            pixels_left = _mat2vec(_missing_read(pixels_left))
                           )
                return R, metadata
            finally
                close(hid)
            end
        end
    end
end

end # module Files
