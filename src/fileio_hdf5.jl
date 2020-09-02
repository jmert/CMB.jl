# Attributes on the data group which must be present to interpret the rest of the
# data stream.
const FORMAT_ATTR_NAMES = ["h5sparse_format"]
const SHAPE_ATTR_NAMES  = ["h5sparse_shape"]

# Recognizable names of the formats which are supported
const CSR_NAMES = ["csr"]
const CSC_NAMES = ["csc"]

# Recognized dataset names for the CSC and CSR formats
const CSX_INDICES_NAMES = ["indices"]
const CSX_POINTER_NAMES = ["indptr"]
const CSX_VALUES_NAMES  = ["data"]

function read_obsmat(file::File{format"HDF5"}, name)
    h5open(FileIO.filename(file), "r") do hfile
        return read_obsmat(hfile, name)
    end
end

function read_obsmat(hfile::HDF5File, name::String, indtype::Type=Int, valtype=Float64)
    @noinline throw_msg(hfile, name, msg) = error("$name $msg in $(file(hfile))")
    @inline function matchname(names, choices, errmsg)
        for nn in choices
            nn in names && return nn
        end
        throw_msg(hfile, name, errmsg)
    end

    group = hfile[name]
    if !isa(group, HDF5Group)
        throw_msg(hfile, name, "is not a group")
    end

    meta = attrs(group)
    attribs = map(lowercase, names(meta))
    dsetnames = map(lowercase, names(group))

    formatfld = matchname(FORMAT_ATTR_NAMES, attribs, "does not have a recognized sparse format attribute")
    shapefld  = matchname(SHAPE_ATTR_NAMES,  attribs, "does not have a recognized matrix shape attribute")

    shape = read(meta[shapefld], Array{Int})::Vector{Int}
    if length(shape) != 2
        throw_msg(hfile, name, "has an unexpected shape attribute")
    end

    format = read(meta[formatfld], String)
    if format != "csc" && format != "csr"
        throw_msg(hfile, name, "is not a recognized sparse format")
    end
    indicesfld  = matchname(CSX_INDICES_NAMES, dsetnames, "does not have a recognized CSC/CSR indices array")
    pointersfld = matchname(CSX_POINTER_NAMES, dsetnames, "does not have a recognized CSC/CSR pointer array")
    valuesfld   = matchname(CSX_VALUES_NAMES,  dsetnames, "does not have a recognized CSC/CSR values array")
    indices  = read(group[indicesfld], Array{indtype})::Vector{indtype}
    pointers = read(group[pointersfld], Array{indtype})::Vector{indtype}
    values   = read(group[valuesfld], Array{valtype})::Vector{valtype}
    indices  .+= one(indtype)
    pointers .+= one(indtype)
    if format == "csr"
        R = SparseMatrixCSC(shape[2], shape[1], pointers, indices, values)
        R = permutedims(R)
    else
        R = SparseMatrixCSC(shape[1], shape[2], pointers, indices, values)
    end

    return R
end

import Base: datatype_alignment

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
