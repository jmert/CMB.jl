"""
Collection of routines for working with coordinates on the sphere.
"""
module Sphere

export
    bearing, distance, cosdistance

doc"""
    bearing(θ₁, ϕ₁, θ₂, ϕ₂) -> α

Calculates at the coordinate $(\theta_1, \phi_1)$ the bearing angle $\alpha$ between North
and the great circle connecting $(\theta_1, \phi_1)$ to $(\theta_2,\phi_2)$. Coordinates are
to be given in radians where $\theta$ is the colatitude angle from the North Pole and $\phi$
is the azimuthal angle.
"""
@inline function bearing{T<:Number}(θ₁::T, ϕ₁::T, θ₂::T, ϕ₂::T)
    local π = convert(T, Base.π)
    # Defining a local isapprox seems to produce far fewer instruction than Base.isapprox.
    local ≈(x::T, y::T) = @fastmath (x==y || abs(x-y) < sqrt(eps(π)))
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
    α = atan2(-sin(ϕ₂₁), cos(θ₁)*cos(ϕ₂₁) - sin(θ₁)*cot(θ₂))
    return mod(α, π)
end

doc"""
    distance(θ₁, ϕ₁, θ₂, ϕ₂) -> σ

Calculates the inner angle $\sigma$ between unit vectors pointing from the center of the
sphere to the two points $(\theta_1,\phi_1)$ and $(\theta_2,\phi_2)$ on its surface, where
$\theta$ is the colatitude angle from the North Pole and $\phi$ is the azimuthal angle, both
in radians.

See also [`cosdistance`](@ref)
"""
@inline function distance{T<:Number}(θ₁::T, ϕ₁::T, θ₂::T, ϕ₂::T)
    const π = convert(T, Base.π)
    # Defining a local isapprox seems to produce far fewer instruction than Base.isapprox.
    const ≈(x::T, y::T) = @fastmath (x==y || abs(x-y) < sqrt(eps(π)))
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

    return 2.0*asin(sqrt(sin(0.5*(θ₂-θ₁))^2 + sin(θ₁)*sin(θ₂)*sin(0.5*(ϕ₂-ϕ₁))^2))
end

doc"""
    cosdistance(θ₁, ϕ₁, θ₂, ϕ₂) -> z

Calculates the cosine of the inner angle $z$ between unit vectors pointing from the center
of the sphere to the two points $(\theta_1,\phi_1)$ and $(\theta_2,\phi_2)$ on its surface,
where $\theta$ is the colatitude angle from the North Pole and $\phi$ is the azimuthal
angle, both in radians.

See also [`distance`](@ref)
"""
@inline function cosdistance{T<:Number}(θ₁::T, ϕ₁::T, θ₂::T, ϕ₂::T)
    local π = convert(T, Base.π)
    # Defining a local isapprox seems to produce far fewer instruction than Base.isapprox.
    local ≈(x::T, y::T) = @fastmath (x==y || abs(x-y) < sqrt(eps(π)))
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

    return one(T) - 2*(sin(0.5*(θ₂-θ₁))^2 + sin(θ₁)*sin(θ₂)*sin(0.5*(ϕ₂-ϕ₁))^2)
end

end
