"""
Collection of routines for working with coordinates on the sphere.
"""
module Sphere

export
    bearing, bearing2, distance, cosdistance

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
@inline simpleapprox(x::T, y::T) where {T} = @fastmath x==y || abs(x-y) < rtepspi(T)

# Make a couple of functions which will make vector math easier for us

"""
    ∥(u, v)

Test whether vector ``u`` is parallel to vector ``v``. Assumes that both are unit
normalized. See also [`⟂`](@ref).
"""
@generated function ∥(u, v)
    T = promote_type(eltype(u), eltype(v))
    return :( ($(one(T)) - abs(u⋅v)) < $(sqrt(eps(one(T)))) )
end

"""
    ⟂(u, v)

Test whether vector ``u`` is perpendicular to vector ``v``. Assumes that both are unit
normalized. See also [`∥`](@ref).
"""
@generated function ⟂(u, v)
    T = promote_type(eltype(u), eltype(v))
    return :( abs(u⋅v) < $(sqrt(eps(one(T)))) )
end

"""
    const x̂ = SVector(1, 0, 0)
    const ŷ = SVector(0, 1, 0)
    const ẑ = SVector(0, 0, 1)

Constant unit vectors in the Cartesian directions.
"""
x̂, ŷ, ẑ

const x̂ = SVector(1, 0, 0)
const ŷ = SVector(0, 1, 0)
const ẑ = SVector(0, 0, 1)

"""
Calculates the bearing angle (``α``), defined as the angle between the meridian (at the
first coordinate) and the great circle connecting the first coordinate to the second. Angles
are measured clockwise and will be in the range ``[0,π)``. See also [`bearing2`](@ref).
"""
function bearing end

"""
    α = bearing(θ₁, ϕ₁, θ₂, ϕ₂)

Points on the sphere are given as coordinate pairs ``(θ₁,ϕ₁)`` and ``(θ₂,ϕ₂)`` where ``θ``
is the colatitude angle from the North Pole and ``ϕ`` is the azimuthal angle, both in
radians.

# Examples
```jldoctest
julia> bearing(π/2, 0.0, π/4, π/4)
0.6154797086703871
```
"""
bearing(θ₁, ϕ₁, θ₂, ϕ₂) = bearing(promote(θ₁, ϕ₁, θ₂, ϕ₂)...)

function bearing(θ₁::T, ϕ₁::T, θ₂::T, ϕ₂::T) where T<:Number
    local π = convert(T, Base.π)
    local ≈ = simpleapprox
    # For a coordinate at a pole, force both longitudes to 0. Since the pole is degenerate,
    # only the latitude makes any difference, and setting to zero makes the parallel checks
    # below also function properly.
    if (θ₁ ≈ zero(T) || θ₁ ≈ π) || (θ₂ ≈ zero(T) || θ₂ ≈ π)
        ϕ₁ = zero(T)
        ϕ₂ = zero(T)
    end

    # For two parallel vectors, return 0 explicitly to avoid hitting any ambiguity in the
    # rotation angle.
    if any(abs(θ₁-θ₂) .≈ (zero(T),π)) || any(abs(ϕ₁-ϕ₂) .≈ (zero(T),2π))
        return zero(T)
    end

    ϕ₂₁ = ϕ₂ - ϕ₁
    num = -sin(ϕ₂₁)
    den = cos(θ₁)*cos(ϕ₂₁) - sin(θ₁)*cot(θ₂)
    # Flip signs of both to move from quadrants 3 and 4 back into 1 and 2 iff the numerator
    # is negative.
    den = flipsign(den, num)
    num = flipsign(num, num)
    return atan2(num, den)
end

"""
    α = bearing(r₁, r₂)

Points on the sphere are given as unit vectors ``r₁`` and ``r₂``.

# Examples
```jldoctest
julia> bearing([1.0, 0.0, 0.0], [0.5, 0.5, sqrt(2)/2])
0.6154797086703873
```
"""
@propagate_inbounds function bearing(r₁::AbstractVector, r₂::AbstractVector)
    r₁ ∥ r₂ && return zero(eltype(r₁))
    r₁₂ = r₁ × r₂
    r₁′ = r₁ ∥ ẑ ? ẑ × r₂ : r₁ × ẑ
    num = (r₁₂ × r₁′) ⋅ r₁
    den = r₁₂ ⋅ r₁′
    # Flip signs of both to move from quadrants 3 and 4 back into 1 and 2 iff the numerator
    # is negative.
    den = flipsign(den, num)
    num = flipsign(num, num)
    return atan2(num, den)
end

"""
Calculates the latitude/longitude vector components of the bearing angle (i.e. ``δθ =
\\cos(α), δϕ = \\sin(α)``), defined as the angle between the meridian (at the first
coordinate) and the great circle connecting the first coordinate to the second. See also
[`bearing`](@ref).
"""
function bearing2 end

"""
    (δθ, δϕ) = bearing2(θ₁, ϕ₁, θ₂, ϕ₂)

Points on the sphere are given as coordinate pairs ``(θ₁,ϕ₁)`` and ``(θ₂,ϕ₂)`` where ``θ``
is the colatitude angle from the North Pole and ``ϕ`` is the azimuthal angle, both in
radians.

# Examples
```jldoctest
julia> bearing2(π/2, 0.0, π/4, π/4)
(0.8164965809277261, 0.5773502691896256)
```
"""
bearing2(θ₁, ϕ₁, θ₂, ϕ₂) = bearing2(promote(θ₁, ϕ₁, θ₂, ϕ₂)...)

function bearing2(θ₁::T, ϕ₁::T, θ₂::T, ϕ₂::T) where T<:Number
    local π = convert(T, Base.π)
    local ≈ = simpleapprox
    # For a coordinate at a pole, force both longitudes to 0. Since the pole is degenerate,
    # only the latitude makes any difference, and setting to zero makes the parallel checks
    # below also function properly.
    if (θ₁ ≈ zero(T) || θ₁ ≈ π) || (θ₂ ≈ zero(T) || θ₂ ≈ π)
        ϕ₁ = zero(T)
        ϕ₂ = zero(T)
    end

    # For two parallel vectors, return [1,0] explicitly to avoid hitting any ambiguity in
    # the rotation angle.
    if any(abs(θ₁-θ₂) .≈ (zero(T),π)) || any(abs(ϕ₁-ϕ₂) .≈ (zero(T),2π))
        return (one(T), zero(T))
    end

    ϕ₂₁ = ϕ₂ - ϕ₁
    num = -sin(ϕ₂₁)
    den = cos(θ₁)*cos(ϕ₂₁) - sin(θ₁)*cot(θ₂)
    # Flip signs of both to move from quadrants 3 and 4 back into 1 and 2 iff the numerator
    # is negative.
    den = flipsign(den, num)
    num = flipsign(num, num)
    # Normalize the 2-vector [den, num]
    len = sqrt(den*den + num*num)
    return (den / len, num / len)
end

"""
    (δθ, δϕ) = bearing2(r₁, r₂)

Points on the sphere are given as unit vectors ``r₁`` and ``r₂``.

# Examples
```jldoctest
julia> bearing2([1.0, 0.0, 0.0], [0.5, 0.5, sqrt(2)/2])
(0.816496580927726, 0.5773502691896257)
```
"""
@propagate_inbounds function bearing2(r₁::AbstractVector, r₂::AbstractVector)
    r₁ ∥ r₂ && return (one(eltype(r₁)), zero(eltype(r₁)))
    r₁₂ = normalize(r₁ × r₂)
    r₁′ = normalize(r₁ ∥ ẑ ? ẑ × r₂ : r₁ × ẑ)
    num = (r₁₂ × r₁′) ⋅ r₁
    den = r₁₂ ⋅ r₁′
    # Flip signs of both to move from quadrants 3 and 4 back into 1 and 2 iff the numerator
    # is negative.
    den = flipsign(den, num)
    num = flipsign(num, num)
    return (den, num)
end

"""
Calculates the inner angle (``σ``) between unit vectors pointing from the center of the
sphere to two points on its surface.
"""
function distance end

"""
    σ = distance(θ₁, ϕ₁, θ₂, ϕ₂)

Points on the sphere are given as coordinate pairs ``(θ₁,ϕ₁)`` and ``(θ₂,ϕ₂)`` where ``θ``
is the colatitude angle from the North Pole and ``ϕ`` is the azimuthal angle, both in
radians.

# Examples
```jldoctest
julia> distance(π/2, 0.0, π/4, π/4)
1.0471975511965979
```
"""
distance(θ₁, ϕ₁, θ₂, ϕ₂) = distance(promote(θ₁, ϕ₁, θ₂, ϕ₂)...)

function distance(θ₁::T, ϕ₁::T, θ₂::T, ϕ₂::T) where T<:Number
    local π = convert(T, Base.π)
    local ≈ = simpleapprox
    # For a coordinate at a pole, force both longitudes to 0. Since the pole is degenerate,
    # only the latitude makes any difference, and setting to zero makes the parallel checks
    # below also function properly.
    if (θ₁ ≈ zero(T) || θ₁ ≈ π) || (θ₂ ≈ zero(T) || θ₂ ≈ π)
        ϕ₁ = zero(T)
        ϕ₂ = zero(T)
    end

    # For two parallel vectors, return 0 or π explicitly to avoid numerical precision
    # limits.
    if abs(θ₁ - θ₂) ≈ zero(T) && abs(ϕ₁ - ϕ₂) ≈ zero(T)
        return zero(T)
    elseif abs(θ₁ - θ₂) ≈ π && abs(ϕ₁ - ϕ₂) ≈ 2π
        return π
    end

    return @fastmath 2.0*asin(sqrt(sin(0.5*(θ₂-θ₁))^2 + sin(θ₁)*sin(θ₂)*sin(0.5*(ϕ₂-ϕ₁))^2))
end

"""
    σ = distance(r₁, r₂)

Points on the sphere are given as unit vectors ``r₁`` and ``r₂``.

# Examples
```jldoctest
julia> distance([1.,0.,0.], [0.5,0.5,sqrt(2)/2])
1.0471975511965979
```
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
    z = cosdistance(θ₁, ϕ₁, θ₂, ϕ₂)

Points on the sphere are given as coordinate pairs ``(θ₁,ϕ₁)`` and ``(θ₂,ϕ₂)`` where ``θ``
is the colatitude angle from the North Pole and ``ϕ`` is the azimuthal angle, both in
radians.

# Examples
```jldoctest
julia> cosdistance(π/2, 0.0, π/4, π/4)
0.5
```
"""
cosdistance(θ₁, ϕ₁, θ₂, ϕ₂) = cosdistance(promote(θ₁, ϕ₁, θ₂, ϕ₂)...)

function cosdistance(θ₁::T, ϕ₁::T, θ₂::T, ϕ₂::T) where T<:Number
    local π = convert(T, Base.π)
    local ≈ = simpleapprox
    # For a coordinate at a pole, force both longitudes to 0. Since the pole is degenerate,
    # only the latitude makes any difference, and setting to zero makes the parallel checks
    # below also function properly.
    if (θ₁ ≈ zero(T) || θ₁ ≈ π) || (θ₂ ≈ zero(T) || θ₂ ≈ π)
        ϕ₁ = zero(T)
        ϕ₂ = zero(T)
    end

    # For two parallel vectors, return 0 or π explicitly to avoid numerical precision
    # limits.
    if abs(θ₁ - θ₂) ≈ zero(T) && abs(ϕ₁ - ϕ₂) ≈ zero(T)
        return one(T)
    elseif abs(θ₁ - θ₂) ≈ π && abs(ϕ₁ - ϕ₂) ≈ 2π
        return -one(T)
    end

    return @fastmath one(T) - 2*(sin(0.5*(θ₂-θ₁))^2 + sin(θ₁)*sin(θ₂)*sin(0.5*(ϕ₂-ϕ₁))^2)
end

"""
    z = cosdistance(r₁, r₂)

Points on the sphere are given as unit vectors ``r₁`` and ``r₂``.

# Examples
```jldoctest
julia> cosdistance([1.,0.,0.], [0.5,0.5,sqrt(2)/2])
0.5
```
"""
@propagate_inbounds function cosdistance(r₁::AbstractVector, r₂::AbstractVector)
    return clamp(r₁ ⋅ r₂, -one(eltype(r₁)), one(eltype(r₁)))
end

end
