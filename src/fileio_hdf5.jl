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
                     pixr::Union{String,Nothing} = "pixels_right",
                     pixl::Union{String,Nothing} = "pixels_left",
                     mmap::Union{Bool,Val{true},Val{false}} = Val(READ_OBSMAT_MMAP[]))
    hfile = h5open(FileIO.filename(file), "r")
    try
        if name !== nothing
            if !exists(hfile, name)
                error(name, " does not exist in ", HDF5.file(hfile))
            end
            group = hfile[name]
            if !isa(group, HDF5Group)
                error(HDF5.name(group), " is not a group in ", HDF5.file(hfile))
            end

            mmap = (mmap isa Bool) ? Val(mmap) : mmap
            R = _read_obsmat(group, mmap)
        else
            R = missing
        end

        pixelsr = (pixr !== nothing && exists(hfile, pixr)) ?
                _read(hfile[pixr]) : missing
        pixelsl = (pixl !== nothing && exists(hfile, pixl)) ?
                _read(hfile[pixl]) : missing

        return (; R = R, pixr = pixelsr, pixl = pixelsl)
    finally
        close(hfile)
    end
end

function _read_obsmat(group::HDF5Group,
                      mmap::Union{Val{true},Val{false}} = Val(READ_OBSMAT_MMAP[]))
    @noinline throw_msg(obj, msg) = error(HDF5.name(obj), " ", msg, " in ", HDF5.file(obj))
    @inline function matchname(names, choices, errmsg)
        for nn in choices
            nn in names && return nn
        end
        throw_msg(group, errmsg)
    end

    usemmap = mmap === Val{true}()

    meta = attrs(group)
    attribs = map(lowercase, names(meta))
    dsetnames = map(lowercase, names(group))

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

function _read(d::Union{HDF5Group,HDF5Dataset})
    if d isa HDF5Dataset
        return read(d)
    else
        return Dict{String,Any}(k => _read(d[k]) for k in names(d))
    end
end

function _write(parent::HDF5Group, name::String, d::Any)
    if d isa Dict
        hobj = g_create(parent, name)
        for (k, v) in d
            _write(hobj, k, v)
        end
    else
        hobj = d_create(parent, name, datatype(d), dataspace(d))
        write(hobj, d)
    end
    return hobj
end

@static if Sys.isunix()
    # To gain access to the Unix mmap flags, we make our own version `HDF5.readmmap`.
    function read_mmap(dset::HDF5Dataset)
        iscontiguous(dset) || error("Cannot mmap discontiguous dataset")

        AT = HDF5.hdf5_to_julia(dset)
        ismmappable(AT) || error("Cannot mmap datasets of type ", AT)
        offset = HDF5.h5d_get_offset(dset)

        io = open(HDF5.filename(dset), read = true)

        flags = READ_OBSMAT_MMAP_FLAGS[]
        data = UnixMmap.mmap(io, AT, size(dset); offset = offset,
                             flags = flags, prot = UnixMmap.PROT_READ)
        UnixMmap.madvise!(data, UnixMmap.MADV_WILLNEED)
        return data
    end
else
    const read_mmap = HDF5.readmmap
end

"""
    write_obsmat(filename::String, obsmat::SparseMatrixCSC, pixr = missing, pixl = missing;
                 obsmat_name::String = "R", pixr_name::String = "pixels_right",
                 pixl_name::String = "pixels_left")

Write an observing matrix `obsmat` to the HDF5 file `filename`. If the file exists, it will
be overwritten completely (to ensure proper data alignment required for memory mapping),
otherwise it will be created. Pixelization descriptions `pixr` and `pixl` may also be
provided; for a `missing` value, the description will not be written to disk.

The HDF5 dataset and group names of the observing matrix and pixelization descriptions can
be changed from their defaults by setting the `obsmat_name`, `pixr_name`, and `pixl_name`
keywords which effect the `obsmat`, `pixr`, and `pixl` storage locations, respectively.

**See also:** [`read_obsmat`](@ref)
"""
function write_obsmat(filename, obsmat::SparseMatrixCSC, pixr = missing, pixl = missing;
                      obsmat_name::String = "R",
                      pixr_name::String = "pixels_right",
                      pixl_name::String = "pixels_left")
    indptr  = getcolptr(obsmat)
    indices = rowvals(obsmat)
    data    = nonzeros(obsmat)

    TI = eltype(indptr)
    TV = eltype(data)
    align = max(datatype_alignment.((TI, TV))...)

    h5open(filename, "w", "alignment", (0, align)) do hfile
        # Write the sparse matrix
        g = g_create(hfile, obsmat_name)
        EARLY = HDF5.H5D_ALLOC_TIME_EARLY
        d_indptr  = d_create(g, "indptr",  datatype(TI), dataspace(size(indptr)),  "alloc_time", EARLY)
        d_indices = d_create(g, "indices", datatype(TI), dataspace(size(indices)), "alloc_time", EARLY)
        d_data    = d_create(g, "data",    datatype(TV), dataspace(size(data)),    "alloc_time", EARLY)
        write(d_indptr, indptr)
        write(d_indices, indices)
        write(d_data, data)
        a_write(g, "format", "csc")
        a_write(g, "shape", [size(obsmat)...])
        a_write(g, "description", "sparse matrix stored in CSC format")

        groot = root(hfile)
        a_write(hfile, "creator", "CMB.jl (https://github.com/jmert/CMB.jl)")
        a_write(hfile, "description", "sparse observing matrix")

        # Write the (optional) pixelization information
        if !ismissing(pixr) && !isnothing(pixr)
            o_pixr = _write(groot, pixr_name, pixr)
            a_write(o_pixr, "description", "pixelization of acted-upon map vector (right-hand side of matrix-vector multiplication)")
        end
        if !ismissing(pixl) && !isnothing(pixl)
            o_pixl = _write(groot, pixl_name, pixl)
            a_write(o_pixl, "description", "pixelization of resultant map vector (left-hand side of matrix-vector multiplication)")
        end
        close(groot)
    end
    return nothing
end