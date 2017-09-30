"""
Collection of routines for working with coordinates on the sphere.
"""
module Sphere

export
    bearing, distance, cosdistance

using StaticArrays
import Base: @propagate_inbounds

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
Calculates the bearing angle (``α``), defined as the angle between the meridian (at the
first coordinate) and the great circle connection the first coordinate to the second. Angles
are measured clockwise and will be in the range ``[0,π)``.
"""
function bearing end

"""
    bearing(θ₁, ϕ₁, θ₂, ϕ₂) -> α

Points on the sphere are given as coordinate pairs ``(θ₁,ϕ₁)`` and ``(θ₂,ϕ₂)`` where ``θ``
is the colatitude angle from the North Pole and ``ϕ`` is the azimuthal angle, both in
radians.
"""
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

# Helper function to get the v × ẑ value, with a special case for SVectors so that the
# compiler can inline all the operations without needing to actually allocate anything.
cross_ẑ(v::AbstractVector) = [v[2], -v[1], zero(eltype(v))]
cross_ẑ(v::SVector) = SVector(v[2], -v[1], zero(eltype(v)))

"""
    bearing(r₁, r₂) -> α

Points on the sphere are given as unit vectors ``r₁`` and ``r₂``.
"""
@propagate_inbounds function bearing(r₁::AbstractVector, r₂::AbstractVector)
    r₁₂ = r₁ × r₂
    r₁′ = r₁ |> cross_ẑ
    num = (r₁₂ × r₁′) ⋅ r₁
    den = r₁₂ ⋅ r₁′
    # Flip signs of both to move from quadrants 3 and 4 back into 1 and 2 iff the numerator
    # is negative.
    den = flipsign(den, num)
    num = flipsign(num, num)
    return @fastmath atan2(num, den)
end

"""
Calculates the inner angle (``σ``) between unit vectors pointing from the center of the
sphere to two points on its surface.
"""
function distance end

"""
    distance(θ₁, ϕ₁, θ₂, ϕ₂) -> σ

Points on the sphere are given as coordinate pairs ``(θ₁,ϕ₁)`` and ``(θ₂,ϕ₂)`` where ``θ``
is the colatitude angle from the North Pole and ``ϕ`` is the azimuthal angle, both in
radians.
"""
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
    distance(r₁, r₂) -> σ

Points on the sphere are given as unit vectors ``r₁`` and ``r₂``.
"""
@propagate_inbounds function distance(r₁::AbstractVector, r₂::AbstractVector)
    return acos(cosdistance(r₁, r₂))
end

"""
Calculates the cosine of the inner angle (``z``) between unit vectors pointing from the
center of the sphere to two points on its surface.
"""
function cosdistance end

"""
    cosdistance(θ₁, ϕ₁, θ₂, ϕ₂) -> z

Points on the sphere are given as coordinate pairs ``(θ₁,ϕ₁)`` and ``(θ₂,ϕ₂)`` where ``θ``
is the colatitude angle from the North Pole and ``ϕ`` is the azimuthal angle, both in
radians.
"""
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

"""
    cosdistance(r₁, r₂) -> z

Points on the sphere are given as unit vectors ``r₁`` and ``r₂``.
"""
@propagate_inbounds function cosdistance(r₁::AbstractVector, r₂::AbstractVector)
    return r₁[1]*r₂[1] + r₁[2]*r₂[2] + r₁[3]*r₂[3]
end



end
