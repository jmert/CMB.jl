using Test, TestSetExtensions
using CMB # For doctests, include CMB as a binding in Main
const NumTypes = (Float32, Float64, BigFloat)

struct IsEqualMatch{T} <: AbstractString
    msg::T
end
Base.isequal(m::IsEqualMatch, o) = occursin(m.msg, o)
function Base.show(io::IO, m::IsEqualMatch)
    print(io, "match")
    if m isa IsEqualMatch{Regex}
        Base.print_quoted(io, m.msg.pattern)
        print(io, "r")
    else
        Base.print_quoted(io, m.msg)
    end
end
macro match_str(msg, flags="")
    if flags == "r"
        return IsEqualMatch(Regex(msg))
    elseif flags == ""
        return IsEqualMatch(msg)
    else
        error("Unrecognized flag(s) `$flags`.")
    end
end

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
                import ..NumTypes, ..@match_str
                include($file)
            end
        end
        printstyled("  ", prettytime(time_ns() - t0), "\n", color=:light_black)
    end
end

@testset ExtendedTestSet "CMB" begin
    @include "numerics.jl" "Numerics utilities"
    @include "conventions.jl" "Conventions and Definitions"
    @include "sphere.jl" "Spherical functions"
    @include "healpix.jl" "HEALPix functions"
    @include "sphericalharmonics.jl" "Spherical Harmonics"
    @include "pixelizations.jl" "Pixelizations"
    @include "pixelcovariance.jl" "Pixel Covariance functions"
    @include "fileio.jl" "File I/O"
    @include "doctests.jl" "Doctests"
end
