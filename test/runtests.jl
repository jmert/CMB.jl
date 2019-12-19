using Test
const NumTypes = (Float32, Float64, BigFloat)

const TESTLIST = Dict(
    "sphere" => "Sphere",
    "legendre" => "Legendre",
    "healpix" => "HEALPix",
    "sparse" => "Sparse"
   )

@testset "CMB" begin
    @testset "$desc" for (id,desc) in TESTLIST
        modpath = joinpath(dirname(@__FILE__), "$(id).jl")
        # Include the file and have it processed within this module
        print("running $desc tests... ")
        t0 = time_ns()
        include(modpath)
        t1 = time_ns()
        println( (t1-t0)/1e9, " seconds")
    end
end
