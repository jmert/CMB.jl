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

function read_obsmat(file::File{format"HDF5"}, name;
                     kws...)
    hfile = h5open(FileIO.filename(file), "r")
    try
        return read_obsmat(hfile, name; kws...)
    finally
        close(hfile)
    end
end

function read_obsmat(hfile::HDF5File, name::String;
                     mmap::Union{Val{true},Val{false}} = Val(READ_OBSMAT_MMAP[]))
    @noinline throw_msg(hfile, name, msg) = error("$name $msg in $(file(hfile))")

    @inline function matchname(names, choices, errmsg)
        for nn in choices
            nn in names && return nn
        end
        throw_msg(hfile, name, errmsg)
    end

    usemmap = mmap === Val{true}()
    group = hfile[name]
    if !isa(group, HDF5Group)
        throw_msg(hfile, name, "is not a group")
    end

    meta = attrs(group)
    attribs = map(lowercase, names(meta))
    dsetnames = map(lowercase, names(group))

    formatfld = matchname(FORMAT_ATTR_NAMES, attribs, "does not have a recognized sparse format attribute")
    shapefld  = matchname(SHAPE_ATTR_NAMES,  attribs, "does not have a recognized matrix shape attribute")

    shape = convert(Vector{Int}, read(meta[shapefld]))::Vector{Int}
    if length(shape) != 2
        throw_msg(hfile, name, "has an unexpected shape attribute")
    end

    format = read(meta[formatfld])::String
    if format != "csc" && format != "csr"
        throw_msg(hfile, name, "is not a recognized sparse format")
    end
    indicesfld  = matchname(CSX_INDICES_NAMES, dsetnames, "does not have a recognized CSC/CSR indices array")
    pointersfld = matchname(CSX_POINTER_NAMES, dsetnames, "does not have a recognized CSC/CSR pointer array")
    valuesfld   = matchname(CSX_VALUES_NAMES,  dsetnames, "does not have a recognized CSC/CSR values array")

    if usemmap
        if format != "csc"
            throw_msg(hfile, name, "cannot be memory mapped; not a CSC sparse matrix")
        end

        indices_dset  = group[indicesfld]
        pointers_dset = group[pointersfld]
        if eltype(indices_dset) != eltype(pointers_dset)
            throw_msg(hfile, name, "cannot be memory mapped; indptr and indices array types differ")
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
            throw_msg(hfile, name, "is 0-indexed, but 1-indexing required if memory-mapping")
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

function write_obsmat(filename, obsmat::SparseMatrixCSC)
    indptr  = getcolptr(obsmat)
    indices = rowvals(obsmat)
    data    = nonzeros(obsmat)

    TI = eltype(indptr)
    TV = eltype(data)
    align = max(datatype_alignment.((TI, TV))...)

    h5open(filename, "w", "alignment", (0, align)) do hfile
        g = g_create(hfile, "R")

        EARLY = HDF5.H5D_ALLOC_TIME_EARLY
        d_indptr  = d_create(g, "indptr",  datatype(TI), dataspace(size(indptr)),  "alloc_time", EARLY)
        d_indices = d_create(g, "indices", datatype(TI), dataspace(size(indices)), "alloc_time", EARLY)
        d_data    = d_create(g, "data",    datatype(TV), dataspace(size(data)),    "alloc_time", EARLY)

        # Now actually write out the data to disk
        write(d_indptr, indptr)
        write(d_indices, indices)
        write(d_data, data)

        # Finally, fill in the attributes required to make this a valid observing matrix
        # to be read back in with `read_obsmat`.
        attrs(g)["format"] = "csc"
        attrs(g)["shape"] = [size(obsmat)...]
    end
    return nothing
end
