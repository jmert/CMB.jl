module Files

using HDF5
using FileIO
using Requires
using SparseArrays
using SparseArrays: getcolptr

export read_obsmat, write_obsmat

"""
    read_obsmat(filename::String, name; kws...)

Read a sparse observing matrix from `filename`, where the data structure is identified
by `name` (which is file type specific).
"""
function read_obsmat(filename::String, name; kws...)
    file = query(filename)
    return read_obsmat(file, name; kws...)
end

"""
    write_obsmat(filename::String, obsmat::SparseMatrixCSC)
"""
function write_obsmat end

include("fileio_hdf5.jl")

function __init__()
    @require JLD = "4138dd39-2aa7-5051-a626-17a0bb65d9c8" begin
        using .JLD
        read_obsmat(file::File{format"JLD"}, name::String) = JLD.load(FileIO.filename(file), name)::SparseMatrixCSC
    end
    @require JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819" begin
        using .JLD2
        read_obsmat(file::File{format"JLD2"}, name::String) = JLD.load(FileIO.filename(file), name)::SparseMatrixCSC
    end

    @require MAT = "23992714-dd62-5051-b70f-ba57cb901cac" begin
        # FileIO doesn't define a MAT format, so do it ourselves
        add_format(format"MAT", FileIO.detecthdf5, [".mat"], [:MAT])
        using .MAT
        function read_obsmat(file::File{format"MAT"}, name::String)
            matopen(FileIO.filename(file), "r") do mat
                return read(mat, name)::SparseMatrixCSC
            end
        end
    end
end

end # module Files
