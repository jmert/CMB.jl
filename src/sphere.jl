"""
Collection of routines for working with coordinates on the sphere.
"""
module Sphere

export
    bearing, distance, cosdistance

# Fast-path the isapprox tolerance for Float32 and Float64. The computation appears to be
# too complex to constant-fold automatically, so we do it manually for these two most
# important case.
@eval rtepspi(::Type{Float32}) = $(sqrt(eps(convert(Float32,π))))
@eval rtepspi(::Type{Float64}) = $(sqrt(eps(convert(Float64,π))))
# Fallback for all other floating point values.
rtepspi(::Type{T}) where {T<:AbstractFloat} = sqrt(eps(convert(T,π)))

# Use a local isapprox function instead of Base.isapprox. We get far fewer instructions with
# this implementation. (Probably related to the keyword-argument penalty?)
local ≈(x::T, y::T) where {T} = @fastmath x==y || abs(x-y) < rtepspi(T)

"""
    bearing(θ₁, ϕ₁, θ₂, ϕ₂) -> α

Calculates at the coordinate ``(θ_1, ϕ_1)`` the bearing angle ``α`` between North and the great
circle connecting ``(θ_1, ϕ_1)`` to ``(θ_2,ϕ_2)``. Coordinates are to be given in radians where
``θ`` is the colatitude angle from the North Pole and ``ϕ`` is the azimuthal angle.
"""
function bearing end

bearing(θ₁, ϕ₁, θ₂, ϕ₂) = bearing(promote(θ₁, ϕ₁, θ₂, ϕ₂)...)

function bearing(θ₁::T, ϕ₁::T, θ₂::T, ϕ₂::T) where T<:Number
    local π = convert(T, Base.π)
    # For a coordinate at a pole, force both longitudes to 0. Since the pole is degenerate,
    # only the latitude makes any difference, and setting to zero makes the parallel checks
    # below also function properly.
    if (θ₁ ≈ zero(T) || θ₁ ≈ π) || (θ₂ ≈ zero(T) || θ₂ ≈ π)
        ϕ₁ = zero(T)
        ϕ₂ = zero(T)
    end

    # For two parallel vectors, return 0 explicitly to avoid hitting any ambiguity in the
    # rotation angle.
    if (θ₁ ≈ θ₂ && ϕ₁ ≈ ϕ₂) || (θ₁ ≈ -θ₂ && ϕ₁ ≈ -ϕ₂)
        return zero(T)
    end

    ϕ₂₁ = ϕ₂ - ϕ₁
    α = @fastmath atan2(-sin(ϕ₂₁), cos(θ₁)*cos(ϕ₂₁) - sin(θ₁)*cot(θ₂))
    return @fastmath mod(α, π)
end


"""
    distance(θ₁, ϕ₁, θ₂, ϕ₂) -> σ

Calculates the inner angle ``σ`` between unit vectors pointing from the center of the sphere
to the two points ``(θ_1,ϕ_1)`` and ``(θ_2,ϕ_2)`` on its surface, where ``θ`` is the colatitude
angle from the North Pole and ``ϕ`` is the azimuthal angle, both in radians.

See also [`cosdistance`](@ref)
"""
function distance end

distance(θ₁, ϕ₁, θ₂, ϕ₂) = distance(promote(θ₁, ϕ₁, θ₂, ϕ₂)...)

function distance(θ₁::T, ϕ₁::T, θ₂::T, ϕ₂::T) where T<:Number
    local π = convert(T, Base.π)
    # For a coordinate at a pole, force both longitudes to 0. Since the pole is degenerate,
    # only the latitude makes any difference, and setting to zero makes the parallel checks
    # below also function properly.
    if (θ₁ ≈ zero(T) || θ₁ ≈ π) || (θ₂ ≈ zero(T) || θ₂ ≈ π)
        ϕ₁ = zero(T)
        ϕ₂ = zero(T)
    end

    # For two parallel vectors, return 0 or π explicitly to avoid numerical precision
    # limits.
    if θ₁ ≈ θ₂ && ϕ₁ ≈ ϕ₂
        return zero(T)
    elseif θ₁ ≈ -θ₂ && ϕ₁ ≈ -ϕ₂
        return π
    end

    return @fastmath 2.0*asin(sqrt(sin(0.5*(θ₂-θ₁))^2 + sin(θ₁)*sin(θ₂)*sin(0.5*(ϕ₂-ϕ₁))^2))
end

"""
    cosdistance(θ₁, ϕ₁, θ₂, ϕ₂) -> z

Calculates the cosine of the inner angle ``z`` between unit vectors pointing from the center
of the sphere to the two points ``(θ_1,ϕ_1)`` and ``(θ_2,ϕ_2)`` on its surface, where ``θ`` is the
colatitude angle from the North Pole and ``ϕ`` is the azimuthal angle, both in radians.

See also [`distance`](@ref)
"""
function cosdistance end

cosdistance(θ₁, ϕ₁, θ₂, ϕ₂) = cosdistance(promote(θ₁, ϕ₁, θ₂, ϕ₂)...)

function cosdistance(θ₁::T, ϕ₁::T, θ₂::T, ϕ₂::T) where T<:Number
    local π = convert(T, Base.π)
    # For a coordinate at a pole, force both longitudes to 0. Since the pole is degenerate,
    # only the latitude makes any difference, and setting to zero makes the parallel checks
    # below also function properly.
    if (θ₁ ≈ zero(T) || θ₁ ≈ π) || (θ₂ ≈ zero(T) || θ₂ ≈ π)
        ϕ₁ = zero(T)
        ϕ₂ = zero(T)
    end

    # For two parallel vectors, return 0 or π explicitly to avoid numerical precision
    # limits.
    if θ₁ ≈ θ₂ && ϕ₁ ≈ ϕ₂
        return one(T)
    elseif θ₁ ≈ -θ₂ && ϕ₁ ≈ -ϕ₂
        return -one(T)
    end

    return @fastmath one(T) - 2*(sin(0.5*(θ₂-θ₁))^2 + sin(θ₁)*sin(θ₂)*sin(0.5*(ϕ₂-ϕ₁))^2)
end

end
