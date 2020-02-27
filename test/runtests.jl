using Test, TestSetExtensions
using CMB # For doctests, include CMB as a binding in Main
const NumTypes = (Float32, Float64, BigFloat)

function prettytime(t)
    v, u = t < 1e3 ? (t, "ns") :
           t < 1e6 ? (t/1e3, "μs") :
           t < 1e9 ? (t/1e6, "ms") :
                     (t/1e9, "s")
    return string(round(v, digits=3), " ", u)
end
macro include(file, desc)
    mod = gensym(first(splitext(file)))
    quote
        print($desc, ": ")
        t0 = time_ns()
        @testset $desc begin
            @eval module $mod
                using Test, CMB
                import ..NumTypes
                include($file)
            end
        end
        printstyled("  ", prettytime(time_ns() - t0), "\n", color=:light_black)
    end
end

@testset ExtendedTestSet "CMB" begin
    @include "sphere.jl" "Spherical functions"
    @include "healpix.jl" "HEALPix functions"
    @include "sphericalharmonics.jl" "Spherical Harmonics"
    @include "util.jl" "Utility functions"
    @include "doctests.jl" "Doctests"
end
