import HDF5, MAT, JLD, JLD2
using SparseArrays

const obsmat_ref = sparse(
        [
            0.1 0.2 0.0 0.0 0.0 0.0
            0.0 0.0 0.5 0.5 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 1.0
        ])
const pathbase = joinpath(@__DIR__, "testdata")

@testset "Reading observing matrices" begin
    @testset "Julia JLD" begin
        !isdefined(Main, :SparseArrays) && @eval Main using SparseArrays
        R = read_obsmat(joinpath(pathbase, "obsmat_sparse.jld"), "R")
        @test R isa SparseMatrixCSC
        @test R == obsmat_ref
    end

    @testset "Julia JLD2" begin
        R = read_obsmat(joinpath(pathbase, "obsmat_sparse.jld2"), "R")
        @test R isa SparseMatrixCSC
        @test R == obsmat_ref
    end

    @testset "Python h5sparse" begin
        cscpath = joinpath(pathbase, "obsmat_sparse_pycsc.h5")
        csrpath = joinpath(pathbase, "obsmat_sparse_pycsr.h5")

        Rcsc = read_obsmat(cscpath, "R")
        @test Rcsc isa SparseMatrixCSC
        @test Rcsc == obsmat_ref

        Rcsr = read_obsmat(csrpath, "R")
        @test Rcsr isa SparseMatrixCSC
        @test Rcsr == obsmat_ref

        # Construct variations of a valid matrix to test error conditions
        function testerror(prepfn, testmsg)
            dest = tempname()
            cp(cscpath, dest; force = true)
            HDF5.h5open(dest, "r+") do hfile
                prepfn(hfile)
            end
            res = @test_throws ErrorException(testmsg) read_obsmat(dest, "R")
            rm(dest)
            return res
        end

        @test_throws ErrorException(match"is not a group") read_obsmat(cscpath, "R/indices")
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
        R = read_obsmat(joinpath(pathbase, "obsmat_sparse.mat"), "R")
        @test R isa SparseMatrixCSC
        @test R == obsmat_ref
    end
end
