"""
Collection of routines for working with coordinates on the sphere.
"""
module Sphere

export
    bearing, bearing2, distance, cosdistance, reckon

using StaticArrays
import LinearAlgebra: ⋅, ×, normalize

rtepsone(::Type{T}) where T = sqrt(eps(one(T)))

function _promote_sphere_type(r₁, r₂)
    T₁, T₂ = float(eltype(r₁)), float(eltype(r₂))
    promote_type(T₁, T₂, typeof(zero(T₁)*zero(T₂) + zero(T₁)*zero(T₂)))
end

"""
    ∥(u, v) -> Bool

Test whether vector ``u`` is parallel to vector ``v``. Assumes that both are unit
normalized.
"""
function ∥(u, v)
    T = _promote_sphere_type(u, v)
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
    r = cartvec(θ, ϕ)
    r = cartvec((θ, ϕ))

Converts the colatitude-azimuth pair ``(θ, ϕ)`` to a Cartesian unit vector ``r``.
"""
function cartvec(θ, ϕ)
    sθ, cθ = sincos(float(θ))
    sϕ, cϕ = sincos(float(ϕ))
    return SVector(sθ * cϕ, sθ * sϕ, cθ)
end
@inline cartvec((θ, ϕ)::Tuple{Any,Any}) = cartvec(θ, ϕ)

"""
    θ, ϕ = colataz(r)

Converts the Cartesian unit vector ``r`` to a colatitude-azimuth pair ``(θ, ϕ)`` in radians.
"""
function colataz(r::SVector{3})
    θ = acos(r[3])
    ϕ = rem2pi(atan(r[2], r[1]), RoundDown)
    (θ, ϕ)
end

"""
    θ, ϕ = colataz(δ, λ)
    θ, ϕ = colataz((δ, λ))

Converts the latitude-longitude pair ``(δ, λ)`` in degrees to colatitude-azimuth ``(θ, ϕ)``
in radians.
"""
function colataz(δ, λ)
    θ, ϕ = unsafe_colataz(δ, λ)
    return (θ, rem2pi(ϕ, RoundDown))
end
@inline colataz((δ, λ)::Tuple{Any,Any}) = colataz(δ, λ)

function unsafe_colataz(δ, λ)
    δ′, λ′ = promote(float(δ), float(λ))
    θ = oftype(δ′, π)/2 - deg2rad(δ′)
    ϕ = deg2rad(λ′)
    return (θ, ϕ)
end
@inline unsafe_colataz((δ, λ)::Tuple{Any,Any}) = unsafe_colataz(δ, λ)

"""
    δ, λ = latlon(r)

Converts the Cartesian unit vector ``r`` to a latitude-longitude pair ``(δ, λ)`` in degrees.
"""
function latlon(r::SVector{3})
    δ = asind(r[3])
    λ = atand(r[2], r[1])
    (δ, λ)
end

"""
    δ, λ = latlon(θ, ϕ)
    δ, λ = latlon((θ, ϕ))

Converts the colatitude-azimuth pair ``(θ, ϕ)`` in radians to latitude-longitude ``(δ, λ)``
in radians.
"""
function latlon(θ, ϕ)
    θ′, ϕ′ = promote(float(θ), float(ϕ))
    local pi = oftype(θ′, π)
    δ = rad2deg(pi/2 - θ′)
    λ = rad2deg(rem2pi(ϕ′, RoundNearest))
    return (δ, λ)
end
@inline latlon((θ, ϕ)::Tuple{Any,Any}) = latlon(θ, ϕ)

"""
Calculates the bearing angle (``α``), defined as the angle between the meridian (at the
first coordinate) and the great circle connecting the first coordinate to the second. Angles
are measured eastward of north and will be in the range ``[-π,π]``. See also
[`bearing2`](@ref).
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
function bearing(θ₁, ϕ₁, θ₂, ϕ₂)
    r₁ = cartvec(θ₁, ϕ₁)
    r₂ = cartvec(θ₂, ϕ₂)
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
function bearing(r₁::AbstractVector, r₂::AbstractVector)
    T = _promote_sphere_type(r₁, r₂)
    let r₁ = SVector{3,T}(r₁), r₂ = SVector{3,T}(r₂)
        r₁ ∥ r₂ && return r₁ ⋅ r₂ < zero(T) ? convert(T, π) : zero(T)
        r₁ ∥ ẑ  && return r₁[3] > zero(T) ? convert(T, π) : zero(T)
        r₂ ∥ ẑ  && return r₂[3] > zero(T) ? zero(T) : convert(T, π)
        r₁₂ = r₁ × r₂
        r₁′ = r₁ × ẑ
        num = (r₁₂ × r₁′) ⋅ r₁
        den = r₁₂ ⋅ r₁′
        return atan(num, den)
    end
end

"""
Calculates the north/east vector components of the bearing angle (i.e.
``δn = \\cos(α), δe = \\sin(α)``), defined as the angle between the meridian (at the first
coordinate) and the great circle connecting the first coordinate to the second. See also
[`bearing`](@ref).
"""
function bearing2 end

"""
    (δn, δe) = bearing2(θ₁, ϕ₁, θ₂, ϕ₂)

Points on the sphere are given as coordinate pairs ``(θ₁,ϕ₁)`` and ``(θ₂,ϕ₂)`` where ``θ``
is the colatitude angle from the North Pole and ``ϕ`` is the azimuthal angle, both in
radians.

# Examples
```jldoctest
julia> bearing2(π/2, 0.0, π/4, π/4)
(0.8164965809277261, 0.5773502691896257)
```
"""
function bearing2(θ₁, ϕ₁, θ₂, ϕ₂)
    r₁ = cartvec(θ₁, ϕ₁)
    r₂ = cartvec(θ₂, ϕ₂)
    return bearing2(r₁, r₂)
end

"""
    (δn, δe) = bearing2(r₁, r₂)

Points on the sphere are given as unit vectors ``r₁`` and ``r₂``.

# Examples
```jldoctest
julia> bearing2([1.0, 0.0, 0.0], [0.5, 0.5, sqrt(2)/2])
(0.816496580927726, 0.5773502691896257)
```
"""
function bearing2(r₁::AbstractVector, r₂::AbstractVector)
    T = _promote_sphere_type(r₁, r₂)
    let r₁ = SVector{3,T}(r₁), r₂ = SVector{3,T}(r₂)
        r₁ ∥ r₂ && return (copysign(one(T), r₁ ⋅ r₂), zero(T))
        r₁ ∥ ẑ  && return (copysign(one(T), -r₁[3]), zero(T))
        r₂ ∥ ẑ  && return (copysign(one(T),  r₂[3]), zero(T))
        r₁₂ = normalize(r₁ × r₂)
        r₁′ = normalize(r₁ × ẑ)
        num = clamp((r₁₂ × r₁′) ⋅ r₁, -one(T), one(T))
        den = clamp(r₁₂ ⋅ r₁′, -one(T), one(T))
        return (den, num)
    end
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
function distance(θ₁, ϕ₁, θ₂, ϕ₂)
    r₁ = cartvec(θ₁, ϕ₁)
    r₂ = cartvec(θ₂, ϕ₂)
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
function distance(r₁::AbstractVector, r₂::AbstractVector)
    T = _promote_sphere_type(r₁, r₂)
    let r₁ = SVector{3,T}(r₁), r₂ = SVector{3,T}(r₂)
        return acos(cosdistance(r₁, r₂))
    end
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
function cosdistance(θ₁, ϕ₁, θ₂, ϕ₂)
    r₁ = cartvec(θ₁, ϕ₁)
    r₂ = cartvec(θ₂, ϕ₂)
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
function cosdistance(r₁::AbstractVector, r₂::AbstractVector)
    T = _promote_sphere_type(r₁, r₂)
    let r₁ = SVector{3,T}(r₁), r₂ = SVector{3,T}(r₂)
        return clamp(r₁ ⋅ r₂, -one(T), one(T))
    end
end

"""
Calculates a position on the sphere a given [`distance`](@ref) (``σ``, in radians) and
relative [`bearing`](@ref) angle (``α``, in radians) away from a given point (measuring the
eastward-of-north orientation of the great circle connecting the source and destination
points with respect to the merdian passing through the source).
"""
function reckon end

"""
    r′ = reckon(r::AbstractVector, σ, α)

The point on the sphere is given as a unit vector ``r``.

!!! note
    When ``r`` points to either the north or south pole, the meridian is defined to be
    **prime meridian** and the bearing angle ``α`` is oriented with respect to it.

    For example, moving a distance ``π/2`` with no bearing goes to the negative ``x``
    axis (i.e. 0° N, 180° W):
    ```jldoctest; setup = (import StaticArrays.SVector)
    julia> reckon([0.0, 0.0, 1.0], π/2, 0.0)
    3-element SVector{3, Float64} with indices SOneTo(3):
     -1.0
      0.0
      6.123233995736766e-17
    ```
"""
function reckon(r::AbstractVector, σ, α)
    sσ, cσ = sincos(σ)
    sα, cα = sincos(α)
    T = typeof(zero(eltype(r)) * sσ * sα)
    r̂ = SVector{3,T}(r)
    if r ∥ ẑ
        e = ŷ
        n = sign(r[3]) * -x̂
    else
        e = ẑ × r̂
        n = ẑ - (r̂⋅ẑ)r̂ # == r × (ẑ × r) == r × e by vector triple product
    end
    d = normalize(n * cα + e * sα)
    r′ = r̂ * cσ + d * sσ
    return r′
end

"""
    (θ′, ϕ′) = reckon(θ, ϕ, σ, α)

The point on the sphere is given by the colatitude-azimuth pair ``(θ, ϕ)``, both given
in radians.

!!! note
    When ``r`` points to either the north or south pole, the meridian is defined to be
    **``θ`` meridian**, and the bearing angle ``α`` is oriented with respect to it.

    For example, moving a distance ``π/2`` with no bearing goes to the equator, with the
    longitude dependent on the input longitude:
    ```jldoctest
    julia> reckon(0.0, 0.0, π/2, 0.0)
    (1.5707963267948966, 3.141592653589793)

    julia> reckon(0.0, π/2, π/2, 0.0)
    (1.5707963267948966, 4.71238898038469)
    ```
"""
function reckon(θ, ϕ, σ, α)
    r = cartvec(θ, ϕ)
    θ′, ϕ′ = colataz(reckon(r, σ, α))
    if r ∥ ẑ
        ϕ′ = mod2pi(ϕ + ϕ′)
    end
    return (θ′, ϕ′)
end

end
