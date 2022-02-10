import Base: datatype_alignment

# Attributes on the data group which must be present to interpret the rest of the
# data stream.
const FORMAT_ATTR_NAMES = ["h5sparse_format", "format"]
const SHAPE_ATTR_NAMES  = ["h5sparse_shape", "shape"]

# Recognizable names of the formats which are supported
const CSR_NAMES = ["csr"]
const CSC_NAMES = ["csc"]

# Recognized dataset names for the CSC and CSR formats
const CSX_INDICES_NAMES = ["indices"]
const CSX_POINTER_NAMES = ["indptr"]
const CSX_VALUES_NAMES  = ["data"]

"""
    READ_OBSMAT_MMAP = Ref{Bool}(true)

Controls the default memory mapping behavior of [`read_obsmat`](@ref read_obsmat).
Defaults to `true`.
"""
const READ_OBSMAT_MMAP = Ref{Bool}(true)

"""
    READ_OBSMAT_MMAP_FLAGS = Ref{UnixMmap.MmapFlags}(MAP_SHARED)

Controls the flags passed to `mmap` when [`read_obsmat`](@ref read_obsmat) memory maps
the observing matrix arrays. Defaults to `MAP_SHARED` on all systems, with Linux also
including `MAP_POPULATE` by default.
"""
const READ_OBSMAT_MMAP_FLAGS = Ref{UnixMmap.MmapFlags}(UnixMmap.MAP_SHARED)

"""
    read_obsmat(file::FileIO.File{format"HDF5"};
                mmap::Union{Bool,Val{true},Val{false}} = READ_OBSMAT_MMAP[],
                keywords...)

Reads an observing matrix and pixelization descriptors from an HDF5 file. The native
format as written by [`write_obsmat`](@ref) as well as `h5sparse`-formatted data
structures are supported.

With HDF5 files in the native format, the additional keyword `mmap` argument controls
whether the observing matrix arrays are loaded via memory-mapping or not.
If memory-mapped on Unix systems, the `MADV_WILLNEED` advice will be applied to the
arrays. Additionally, on Linux systems, the memory mapping includes the `MAP_POPULATE`
flag to pre-fault the data pages into active RAM. The default behavior is controlled
by the global `Ref` [`READ_OBSMAT_MMAP`](@ref).
"""
function read_obsmat(file::File{format"HDF5"};
                     name::Union{String,Nothing} = "R",
                     fields::Union{String,Nothing} = "fields",
                     pixels_right::Union{String,Nothing} = "pixels_right",
                     pixels_left::Union{String,Nothing} = "pixels_left",
                     mmap::Union{Bool,Val{true},Val{false}} = Val(READ_OBSMAT_MMAP[]))
    hfile = h5open(FileIO.filename(file), "r")
    try
        if name !== nothing
            if !haskey(hfile, name)
                error(name, " does not exist in ", HDF5.file(hfile))
            end
            group = hfile[name]
            if !isa(group, HDF5.Group)
                error(HDF5.name(group), " is not a group in ", HDF5.file(hfile))
            end

            mmap = (mmap isa Bool) ? Val(mmap) : mmap
            R = _read_obsmat(group, mmap)
        else
            R = missing
        end

        @inline function _missing_read(name)
            name === nothing && return missing
            !haskey(hfile, name) && return missing
            return _read(hfile[name])
        end
        metadata = (;
                    fields = _missing_read(fields),
                    pixels_right = _missing_read(pixels_right),
                    pixels_left = _missing_read(pixels_left)
                   )
        return R, metadata
    finally
        close(hfile)
    end
end

function _read_obsmat(group::HDF5.Group,
                      mmap::Union{Val{true},Val{false}} = Val(READ_OBSMAT_MMAP[]))
    @noinline throw_msg(obj, msg) = error(HDF5.name(obj), " ", msg, " in ", HDF5.file(obj))
    @inline function matchname(names, choices, errmsg)
        for nn in choices
            nn in names && return nn
        end
        throw_msg(group, errmsg)
    end

    usemmap = mmap === Val{true}()

    meta = attributes(group)
    attribs = map(lowercase, keys(meta))
    dsetnames = map(lowercase, keys(group))

    formatfld = matchname(FORMAT_ATTR_NAMES, attribs, "does not have a recognized sparse format attribute")
    shapefld  = matchname(SHAPE_ATTR_NAMES,  attribs, "does not have a recognized matrix shape attribute")

    shape = convert(Vector{Int}, read(meta[shapefld]))::Vector{Int}
    if length(shape) != 2
        throw_msg(group, "has an unexpected shape attribute")
    end

    format = read(meta[formatfld])::String
    if format != "csc" && format != "csr"
        throw_msg(group, "is not a recognized sparse format")
    end
    indicesfld  = matchname(CSX_INDICES_NAMES, dsetnames, "does not have a recognized CSC/CSR indices array")
    pointersfld = matchname(CSX_POINTER_NAMES, dsetnames, "does not have a recognized CSC/CSR pointer array")
    valuesfld   = matchname(CSX_VALUES_NAMES,  dsetnames, "does not have a recognized CSC/CSR values array")

    if usemmap
        if format != "csc"
            throw_msg(group, "cannot be memory mapped; not a CSC sparse matrix")
        end

        indices_dset  = group[indicesfld]
        pointers_dset = group[pointersfld]
        if eltype(indices_dset) != eltype(pointers_dset)
            throw_msg(group, "cannot be memory mapped; indptr and indices array types differ")
        end

        indices  = read_mmap(indices_dset)
        pointers = read_mmap(pointers_dset)
        values   = read_mmap(group[valuesfld])
    else
        indices  = read(group[indicesfld])
        pointers = read(group[pointersfld])
        values   = read(group[valuesfld])
        indices, pointers = promote(indices, pointers)
    end

    if length(pointers) > 0 && iszero(pointers[1])
        if usemmap
            throw_msg(group, "is 0-indexed, but 1-indexing required if memory-mapping")
        else
            # adjust for 0-based indexing
            indices  .+= one(eltype(indices))
            pointers .+= one(eltype(pointers))
        end
    end

    if format == "csr"
        R = SparseMatrixCSC(shape[2], shape[1], pointers, indices, values)
        R = permutedims(R)
    else
        R = SparseMatrixCSC(shape[1], shape[2], pointers, indices, values)
    end
    return R
end

function _read(d::Union{HDF5.Group,HDF5.Dataset})
    if d isa HDF5.Dataset
        return read(d)
    else
        return Dict{String,Any}(k => _read(d[k]) for k in keys(d))
    end
end

function _write(parent::HDF5.Group, name::String, d::Any)
    if d isa Dict
        hobj = create_group(parent, name)
        for (k, v) in d
            _write(hobj, k, v)
        end
    else
        dtype = datatype(d)
        hobj = create_dataset(parent, name, dtype, dataspace(d))
        write_dataset(hobj, dtype, d)
    end
    return hobj
end

@static if Sys.isunix()
    # To gain access to the Unix mmap flags, we make our own version `HDF5.readmmap`.
    function read_mmap(dset::HDF5.Dataset)
        HDF5.iscontiguous(dset) || error("Cannot mmap discontiguous dataset")

        AT = HDF5.get_jl_type(dset)
        HDF5.ismmappable(AT) || error("Cannot mmap datasets of type ", AT)
        offset = HDF5.API.h5d_get_offset(dset)

        io = open(HDF5.filename(dset), read = true)

        flags = READ_OBSMAT_MMAP_FLAGS[]
        data = UnixMmap.mmap(io, Array{AT}, size(dset); offset = offset,
                             flags = flags, prot = UnixMmap.PROT_READ)
        UnixMmap.madvise!(data, UnixMmap.MADV_WILLNEED)
        return data
    end
else
    const read_mmap = HDF5.readmmap
end

"""
    write_obsmat(filename::String, obsmat::SparseMatrixCSC;
                 obsmat_name::String = "R",
                 fields = missing, fields_name::String = "fields",
                 pixels_right = missing, pixels_right_name::String = "pixels_right",
                 pixels_left  = missing, pixels_left_name::String  = "pixels_left")

Write an observing matrix `obsmat` to the HDF5 file `filename`. If the file exists, it will
be overwritten completely (to ensure proper data alignment required for memory mapping),
otherwise it will be created.

Additional metadata — pixelization descriptions `pixels_right` and `pixels_left` and an
annotation of the Stokes field the observing matrix applies to `fields` —  may also be
provided; for a `missing` value, the description will not be written to disk.

The HDF5 dataset and group names of the observing matrix and metadata can
be changed from their defaults by setting the associated `*_name` keywords which effect.

**See also:** [`read_obsmat`](@ref)
"""
function write_obsmat(filename, obsmat::SparseMatrixCSC;
                      obsmat_name::String = "R",
                      fields = missing, fields_name::String = "fields",
                      pixels_right = missing, pixels_right_name::String = "pixels_right",
                      pixels_left  = missing, pixels_left_name::String  = "pixels_left")
    indptr  = getcolptr(obsmat)
    indices = rowvals(obsmat)
    data    = nonzeros(obsmat)

    TI = eltype(indptr)
    TV = eltype(data)
    align = max(datatype_alignment.((TI, TV))...)

    h5open(filename, "w"; alignment = (0, align)) do hfile
        # Write the sparse matrix
        g = create_group(hfile, obsmat_name)
        g["indptr",  alloc_time = :early] = indptr
        g["indices", alloc_time = :early] = indices
        g["data",    alloc_time = :early] = data
        a = attributes(g)
        a["format"] = "csc"
        a["shape"] = [size(obsmat)...]
        a["description"] = "sparse matrix stored in CSC format"
        close(g)

        groot = HDF5.root(hfile)
        a = attributes(hfile)
        a["creator"] = "CMB.jl (https://github.com/jmert/CMB.jl)"
        a["description"] = "sparse observing matrix"

        # Write the (optional) metadata
        if !ismissing(fields) && !isnothing(fields)
            o_fields = _write(groot, fields_name, fields)
            a = attributes(o_fields)
            a["description"] = "the Stokes fields to which the observing matrix applies"
            close(o_fields)
        end
        if !ismissing(pixels_right) && !isnothing(pixels_right)
            o_pixels_right = _write(groot, pixels_right_name, pixels_right)
            a = attributes(o_pixels_right)
            a["description",] = "pixelization of acted-upon map vector (right-hand side of matrix-vector multiplication)"
            close(o_pixels_right)
        end
        if !ismissing(pixels_left) && !isnothing(pixels_left)
            o_pixels_left = _write(groot, pixels_left_name, pixels_left)
            a = attributes(o_pixels_left)
            a["description"] = "pixelization of resultant map vector (left-hand side of matrix-vector multiplication)"
            close(o_pixels_left)
        end
        close(groot)
    end
    return nothing
end
