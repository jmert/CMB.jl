"""
Miscellaneous utility functions.
"""
module Util

include("sparse.jl")

"""
    @relerr fncall(args...)

Takes the function call `fncall(args...)` and rewrites the expression to calculate the
relative deviation (in ulps) between the regular call and one with all arguments promoted to
`BigFloat` or `BigInt`.

# Example
```jldoctest
julia> CMB.Util.@relerr sin(π-1e-5)
0.24759425226749643
```

# See Also
[`@abserr`](@ref), [`@absrelerr`](@ref)
"""
macro relerr(fncall)
    (fncall isa Expr && fncall.head == :call) ||
            error("Expected a function call, found $(fncall)")
    fn = esc(fncall.args[1])
    args = esc.(fncall.args[2:end])
    bigargs = Any[:(big($a)) for a in args]

    quote
        begin
            r₁ = $fn($(args...))
            r₂ = $fn($(bigargs...))
            u = (r₁-r₂) / eps(abs(r₁))
            convert(Float64, u)
        end
    end
end

"""
    @abserr fncall(args...)

Takes the function call `fncall(args...)` and rewrites the expression to calculate the
absolute deviation between the regular call and one with all arguments promoted to
`BigFloat` or `BigInt`.

# Example
```jldoctest
julia> CMB.Util.@abserr sin(π-1e-5)
4.1944097844272447e-22
```

# See Also
[`@relerr`](@ref), [`@absrelerr`](@ref)
"""
macro abserr(fncall)
    (fncall isa Expr && fncall.head == :call) ||
            error("Expected a function call, found $(fncall)")
    fn = esc(fncall.args[1])
    args = esc.(fncall.args[2:end])
    bigargs = Any[:(big($a)) for a in args]

    quote
        begin
            r₁ = $fn($(args...))
            r₂ = $fn($(bigargs...))
            convert(Float64, r₁-r₂)
        end
    end
end

"""
    @absrelerr fncall(args...)

Takes the function call `fncall(args...)` and rewrites the expression to calculate the
absolute and relative deviation between the regular call and one with all arguments promoted
to `BigFloat` or `BigInt`, returning both as a tuple pair.

# Example
```jldoctest
julia> CMB.Util.@absrelerr sin(π-1e-5)
(4.1944097844272447e-22, 0.24759425226749643)
```

# See Also
[`@abserr`](@ref), [`@relerr`](@ref)
"""
macro absrelerr(fncall)
    (fncall isa Expr && fncall.head == :call) ||
            error("Expected a function call, found $(fncall)")
    fn = esc(fncall.args[1])
    args = esc.(fncall.args[2:end])
    bigargs = Any[:(big($a)) for a in args]

    quote
        begin
            r₁ = $fn($(args...))
            r₂ = $fn($(bigargs...))
            a = convert(Float64, r₁-r₂)
            u = convert(Float64, a / eps(abs(r₁)))
            (a,u)
        end
    end
end

end
