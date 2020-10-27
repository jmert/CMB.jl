import HDF5, MAT, JLD, JLD2
using SparseArrays

const obsmat_ref = sparse(
        [
            0.1 0.2 0.0 0.0 0.0 0.0
            0.0 0.0 0.5 0.5 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 1.0
        ])
const pixr_ref = Dict{String,Any}(
        "type" => "dummy_pixels",
        "index" => collect(0:5),
        "sub" => Dict{String,Any}("extra" => Int8(1))
        )
const pixl_ref = collect(1:4)
const meta_ref = (;
                  fields = "T",
                  pixels_right = pixr_ref,
                  pixels_left = pixl_ref
                 )
const meta_nothing = map(_ -> nothing, meta_ref)
const meta_nonexist = map(_ -> "doesntexist", meta_ref)

const pathbase = joinpath(@__DIR__, "testdata")

CMB.Files.READ_OBSMAT_MMAP[] = false

@testset "Reading observing matrices" begin
    @testset "Julia JLD" begin
        !isdefined(Main, :SparseArrays) && @eval Main using SparseArrays
        !isdefined(Main, :JLD) && @eval Main using JLD
        source = joinpath(pathbase, "obsmat_sparse.jld")

        R, meta = read_obsmat(source)
        @test R isa SparseMatrixCSC
        @test R == obsmat_ref
        @test meta == meta_ref

        # pixel descriptors are optional so ignore non-existent names
        R, meta = read_obsmat(source; meta_nonexist...)
        @test all(ismissing, meta)
        # can also explicitly say they should not be read
        R, meta = read_obsmat(source; meta_nothing...)
        @test all(ismissing, meta)

        # Also valid to explicitly specify not loading the actual matrix (so that the
        # pixel descriptions can be loaded first for pre-processing, for example)
        R, meta = read_obsmat(source, name = nothing)
        @test ismissing(R)
        @test meta == meta_ref
        # but it's an error to request a matrix name which doesn't exist or is invalid
        @test_throws ErrorException read_obsmat(source, name = "R_dne")
    end

    @testset "Julia JLD2" begin
        source = joinpath(pathbase, "obsmat_sparse.jld2")

        R, meta = read_obsmat(source)
        @test R isa SparseMatrixCSC
        @test R == obsmat_ref
        @test meta == meta_ref

        # pixel descriptors are optional so ignore non-existent names
        R, meta = read_obsmat(source; meta_nonexist...)
        @test all(ismissing, meta)
        # can also explicitly say they should not be read
        R, meta = read_obsmat(source; meta_nothing...)
        @test all(ismissing, meta)

        # Also valid to explicitly specify not loading the actual matrix (so that the
        # pixel descriptions can be loaded first for pre-processing, for example)
        R, meta = read_obsmat(source, name = nothing)
        @test ismissing(R)
        @test meta == meta_ref
        # but it's an error to request a matrix name which doesn't exist or is invalid
        @test_throws KeyError read_obsmat(source, name = "R_dne")
    end

    @testset "Python h5sparse" begin
        cscpath = joinpath(pathbase, "obsmat_sparse_pycsc.h5")
        csrpath = joinpath(pathbase, "obsmat_sparse_pycsr.h5")

        R, meta = read_obsmat(cscpath)
        @test R isa SparseMatrixCSC
        @test R == obsmat_ref
        @test meta == meta_ref

        # pixel descriptors are optional so ignore non-existent names
        R, meta = read_obsmat(cscpath; meta_nonexist...)
        @test all(ismissing, meta)
        # can also explicitly say they should not be read
        R, meta = read_obsmat(cscpath; meta_nothing...)
        @test all(ismissing, meta)

        # Also valid to explicitly specify not loading the actual matrix (so that the
        # pixel descriptions can be loaded first for pre-processing, for example)
        R, meta = read_obsmat(cscpath, name = nothing)
        @test ismissing(R)
        @test meta == meta_ref
        # but it's an error to request a matrix name which doesn't exist or is invalid
        @test_throws ErrorException(match"does not exist") read_obsmat(cscpath, name = "R_dne")
        @test_throws ErrorException(match"is not a group") read_obsmat(cscpath, name = "R/indices")

        R, meta = read_obsmat(csrpath)
        @test R isa SparseMatrixCSC
        @test R == obsmat_ref
        @test meta == meta_ref

        # Construct variations of a valid matrix to test error conditions
        function testerror(prepfn, testmsg)
            dest = tempname()
            cp(cscpath, dest; force = true, follow_symlinks=true)
            HDF5.h5open(dest, "r+") do hfile
                prepfn(hfile)
            end
            res = @test_throws ErrorException(testmsg) read_obsmat(dest)
            rm(dest)
            return res
        end

        testerror(match"R does not exist") do hfile
            HDF5.o_delete(hfile["R"])
        end
        testerror(match"does not have a recognized sparse format attribute") do hfile
            HDF5.a_delete(hfile["R"], "h5sparse_format")
        end
        testerror(match"does not have a recognized matrix shape attribute") do hfile
            HDF5.a_delete(hfile["R"], "h5sparse_shape")
        end
        testerror(match"is not a recognized sparse format") do hfile
            HDF5.a_delete(hfile["R"], "h5sparse_format")
            HDF5.attrs(hfile["R"])["h5sparse_format"] = "coo"
        end
        testerror(match"has an unexpected shape attribute") do hfile
            HDF5.a_delete(hfile["R"], "h5sparse_shape")
            HDF5.attrs(hfile["R"])["h5sparse_shape"] = Int[2, 3, 4]
        end
        testerror(match"does not have a recognized CSC/CSR indices array") do hfile
            HDF5.o_delete(hfile["R/indices"])
        end
        testerror(match"does not have a recognized CSC/CSR pointer array") do hfile
            HDF5.o_delete(hfile["R/indptr"])
        end
        testerror(match"does not have a recognized CSC/CSR values array") do hfile
            HDF5.o_delete(hfile["R/data"])
        end
    end

    @testset "Matlab" begin
        source = joinpath(pathbase, "obsmat_sparse.mat")

        R, meta = read_obsmat(source)
        @test R isa SparseMatrixCSC
        @test R == obsmat_ref
        @test meta == meta_ref

        # pixel descriptors are optional so ignore non-existent names
        R, meta = read_obsmat(source; meta_nonexist...)
        @test all(ismissing, meta)
        # can also explicitly say they should not be read
        R, meta = read_obsmat(source; meta_nothing...)
        @test all(ismissing, meta)

        # Also valid to explicitly specify not loading the actual matrix (so that the
        # pixel descriptions can be loaded first for pre-processing, for example)
        R, meta = read_obsmat(source, name = nothing)
        @test ismissing(R)
        @test meta == meta_ref
        # but it's an error to request a matrix name which doesn't exist or is invalid
        @test_throws ErrorException read_obsmat(source, name = "R_dne")
    end
end

@testset "Writing observing matrices" begin
    mktemp() do path, io
        write_obsmat(path, obsmat_ref; meta_ref...)
        flush(io)
        @test HDF5.ishdf5(path)
        R, meta = read_obsmat(path)
        @test R == obsmat_ref
        @test meta == meta_ref
    end

    mktemp() do path, io
        # write observing matrix file without metadata
        write_obsmat(path, obsmat_ref)
        flush(io)
        @test HDF5.ishdf5(path)
        R, meta = read_obsmat(path)
        @test R == obsmat_ref
        @test all(ismissing, meta)
    end
end

@testset "Memory-mapping" begin
    # Check that the reference matrix is memory mappable. The matrix requires no
    # padding (for either Int32 or Int64 indices).
    mktemp() do path, io
        write_obsmat(path, obsmat_ref)
        flush(io)
        R, meta = read_obsmat(path, mmap = Val(true))
        @test R == obsmat_ref
    end

    # Construct matrix which requires padding --- Int16 pointers for a 6x4 matrix with
    # 5 nonzeros has (5 + 5) * 2 = 20 bytes of indices, which is not a multiple of
    # 8 bytes (alignment of Float64 data array).
    mktemp() do path, io
        R′ = convert(SparseMatrixCSC{Float64,Int16}, obsmat_ref')
        @test sizeof(rowvals(R′)) + sizeof(SparseArrays.getcolptr(R′)) == 20
        write_obsmat(path, R′)
        flush(io)
        R, meta = read_obsmat(path, mmap = Val(true))
        @test R == R′
    end

    # Cannot mmap a non-CSC matrix
    @test_throws ErrorException(match"cannot be memory mapped; not a CSC sparse matrix"
                               ) read_obsmat(joinpath(pathbase, "obsmat_sparse_pycsr.h5"), mmap = Val(true))

    # Cannot mmap a CSC matrix where indptr and indices arrays are different integer
    # types.
    @test_throws ErrorException(match"cannot be memory mapped; indptr and indices array types differ"
                               ) read_obsmat(joinpath(pathbase, "obsmat_sparse_pycsc_i64i32f64.h5"), mmap = Val(true))

    # Python's sparse saves use 0-index index arrays which cannot be used if we want
    # to memory map them into Julia data structures.
    @test_throws ErrorException(match"is 0-indexed, but 1-indexing required if memory-mapping"
                               ) read_obsmat(joinpath(pathbase, "obsmat_sparse_pycsc_i32i32f64.h5"), mmap = Val(true))
end
