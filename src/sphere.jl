"""
Collection of routines for working with coordinates on the sphere.
"""
module Sphere

export
    bearing, bearing2, distance, cosdistance

using StaticArrays
import Base: @propagate_inbounds
import LinearAlgebra: ⋅, ×, normalize

@eval rtepsone(::Type{Float32}) = $(sqrt(eps(one(Float32))))
@eval rtepsone(::Type{Float64}) = $(sqrt(eps(one(Float64))))
rtepsone(::Type{T}) where T = sqrt(eps(one(T)))

"""
    ∥(u, v)

Test whether vector ``u`` is parallel to vector ``v``. Assumes that both are unit
normalized.
"""
function ∥(u, v)
    T = promote_type(eltype(u), eltype(v))
    return (one(T) - abs(u⋅v)) < rtepsone(T)
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
    cart(θ, ϕ)

Converts the colatitude-azimuth pair to a Cartesian unit vector.
"""
function cart(θ, ϕ)
    sθ, cθ = sincos(θ)
    sϕ, cϕ = sincos(ϕ)
    return SVector(sθ * cϕ, sθ * sϕ, cθ)
end

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
0.6154797086703873
```
"""
bearing(θ₁, ϕ₁, θ₂, ϕ₂) = bearing(promote(θ₁, ϕ₁, θ₂, ϕ₂)...)

function bearing(θ₁::T, ϕ₁::T, θ₂::T, ϕ₂::T) where T<:Number
    r₁ = cart(θ₁, ϕ₁)
    r₂ = cart(θ₂, ϕ₂)
    return bearing(r₁, r₂)
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
    T = eltype(r₁)
    r₁ ∥ r₂ && return r₁ ⋅ r₂ < zero(T) ? convert(T, π) : zero(T)
    r₁ ∥ ẑ  && return r₁[3] > zero(T) ? convert(T, π) : zero(T)
    r₂ ∥ ẑ  && return r₂[3] > zero(T) ? zero(T) : convert(T, π)
    r₁₂ = r₁ × r₂
    r₁′ = r₁ × ẑ
    num = (r₁₂ × r₁′) ⋅ r₁
    den = r₁₂ ⋅ r₁′
    # Flip signs of both to move from quadrants 3 and 4 back into 1 and 2 iff the numerator
    # is negative.
    den = flipsign(den, num)
    num = flipsign(num, num)
    return atan(num, den)
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
(0.8164965809277261, 0.5773502691896257)
```
"""
bearing2(θ₁, ϕ₁, θ₂, ϕ₂) = bearing2(promote(θ₁, ϕ₁, θ₂, ϕ₂)...)

function bearing2(θ₁::T, ϕ₁::T, θ₂::T, ϕ₂::T) where T<:Number
    r₁ = cart(θ₁, ϕ₁)
    r₂ = cart(θ₂, ϕ₂)
    return bearing2(r₁, r₂)
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
    T = eltype(r₁)
    r₁ ∥ r₂ && return (copysign(one(T), r₁ ⋅ r₂), zero(T))
    r₁ ∥ ẑ  && return (copysign(one(T), -r₁[3]), zero(T))
    r₂ ∥ ẑ  && return (copysign(one(T),  r₂[3]), zero(T))
    r₁₂ = normalize(r₁ × r₂)
    r₁′ = normalize(r₁ × ẑ)
    num = clamp((r₁₂ × r₁′) ⋅ r₁, -one(T), one(T))
    den = clamp(r₁₂ ⋅ r₁′, -one(T), one(T))
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
    r₁ = cart(θ₁, ϕ₁)
    r₂ = cart(θ₂, ϕ₂)
    return distance(r₁, r₂)
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
    r₁ = cart(θ₁, ϕ₁)
    r₂ = cart(θ₂, ϕ₂)
    return cosdistance(r₁, r₂)
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
