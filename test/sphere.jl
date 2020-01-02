module Sphere
using Test
using LinearAlgebra, Random, StaticArrays
using CMB.Sphere
import ..NumTypes
import CMB.Sphere: cart, ẑ

@testset "Antipodes ($T)" for T in NumTypes
    npole = SVector{3,T}(ẑ)
    spole = SVector{3,T}(-ẑ)

    # North and South Pole antipodes
    @test @inferred(bearing(    npole, spole)) == T(π)
    @test @inferred(bearing(    spole, npole)) == T(π)
    @test @inferred(bearing2(   npole, spole)) == (T(-1.0), T(0.0))
    @test @inferred(bearing2(   spole, npole)) == (T(-1.0), T(0.0))
    @test @inferred(distance(   npole, spole)) == T(π)
    @test @inferred(cosdistance(npole, spole)) == T(-1.0)

    # North and North Pole
    @test @inferred(bearing(    npole, npole)) == T(0.0)
    @test @inferred(bearing2(   npole, npole)) == (T(1.0), T(0.0))
    @test @inferred(distance(   npole, npole)) == T(0.0)
    @test @inferred(cosdistance(npole, npole)) == T(1.0)

    # South and South Pole
    @test @inferred(bearing(    spole, spole)) == T(0.0)
    @test @inferred(bearing2(   spole, spole)) == (T(1.0), T(0.0))
    @test @inferred(distance(   spole, spole)) == T(0.0)
    @test @inferred(cosdistance(spole, spole)) == T(1.0)
end

# Pole with arbitrary points around equator. Distance at the equator is easy to know
# analytically, so vary longitudes to not nice multiples of π.
@testset "Pole to Equator ($T)" for T in NumTypes
    local pi = T===BigFloat ? T(π) : 1.0π
    atol = max(eps(1π), eps(T(π)))

    # Use colatitude-azimuth pairs to test angle interfaces
    npole = (T(0.0), T(0.0))
    for ii=1:7
        eqpnt = (T(pi/2), T(2pi*ii/7))
        @test @inferred(bearing(    npole..., eqpnt...)) == T(pi)
        @test @inferred(bearing(    eqpnt..., npole...)) == T(0.0)
        @test @inferred(bearing2(   npole..., eqpnt...)) == (T(-1.0), T(0.0))
        @test @inferred(bearing2(   eqpnt..., npole...)) == (T( 1.0), T(0.0))
        @test @inferred(distance(   npole..., eqpnt...)) ≈ T(pi/2) atol=atol
        @test @inferred(cosdistance(npole..., eqpnt...)) ≈ T(0.0)  atol=atol
    end
end

# Pairs of points on the equator will always have bearing angles of 90°, and
# the angular separation is equal to the difference in azimuths.
@testset "Points around equator ($T)" for T in NumTypes
    local pi = T===BigFloat ? T(π) : 1.0π
    atol = max(eps(1π), eps(T(π)))

    Random.seed!(1928)
    for nn in 1:4
        θ  = T(pi / 2)
        ϕ₁ = T(pi) * rand(T)
        for Δϕ in T.((pi/8, pi/4, pi/2, 3pi/4, 7pi/8))
            ϕ₂ = ϕ₁ + Δϕ
            @test all((≈).(bearing2(θ, ϕ₁, θ, ϕ₂), T.((0, 1)), atol=atol))
            @test all((≈).(bearing2(θ, ϕ₂, θ, ϕ₁), T.((0, 1)), atol=atol))
            @test bearing( θ, ϕ₁, θ, ϕ₂) ≈ T(pi / 2) atol=atol
            @test bearing( θ, ϕ₂, θ, ϕ₁) ≈ T(pi / 2) atol=atol
            @test distance(θ, ϕ₁, θ, ϕ₂) ≈ Δϕ atol=atol
            @test distance(θ, ϕ₂, θ, ϕ₁) ≈ Δϕ atol=atol
        end
    end
end

# For a pair of points on the same latitude and separated in azimuth by 180°, the
# angular separation will be twice the colatitude angle.
@testset "Isolatitudes ($T)" for T in NumTypes
    local pi = T===BigFloat ? T(π) : 1.0π
    # TODO: Figure out if there's a way to improve the situation so that this can be
    #       dropped down to just max(eps) like the other test sets.
    atol = sqrt(max(T(eps(1π)), eps(T(π))))

    Random.seed!(4759)
    for nn in 1:4
        θ  = T(pi) * rand(T)
        ϕ₁ = T(pi) * rand(T)
        ϕ₂ = ϕ₁ + T(π)
        σ  = θ > pi/2 ? T(2pi) - 2θ : 2θ

        @test distance(θ, ϕ₁, θ, ϕ₂) ≈ σ atol=atol
        @test distance(θ, ϕ₂, θ, ϕ₁) ≈ σ atol=atol
    end
end

@testset "Clamping to ±1" begin
    # Specific case where
    #   norm(a) == 1.0 && dot(a, a) > 1.0
    a = [0.6468729683743681,-0.18340815863750012,-0.7402140299479167]
    @test norm(a) == 1.0
    @test dot(a, a) > 1.0
    @test cosdistance(a, a) == 1.0

    # Specific case where internal dot product returns value > 1.0
    r₁ = [0.0,0.0,1.0]
    r₂ = [cos(2π*3/7), sin(2π*3/7), 0.0]
    @test bearing2(r₁, r₂) == (-1.0, 0.0)
end

# Test failed with world age problems when generated functions were incorrectly used;
# numerical derivative via newer world's Dual numbers tests for this error.
@testset "Derivatives via dual numbers" begin
    using ForwardDiff: gradient
    p₁ = @SVector[π/4, 0.0]
    p₂ = @SVector[π/2, π/4]
    pts = SVector{4,Float64}([p₁..., p₂...])

    # Analytic derivatives for cosdistance with respect to θ,ϕ
    function anal_deriv1(pts)
        r₁ = Sphere.cart(pts[1:2]...)
        r₂ = Sphere.cart(pts[3:4]...)
        @inline function dcosσ_dθ(r₁, r₂)
            x₁, y₁, z₁ = r₁
            x₂, y₂, z₂ = r₂
            ζ = @fastmath sqrt(1 - z₁^2)
            return (x₁*x₂ + y₁*y₂) * z₁/ζ - ζ*z₂
        end
        @inline function dcosσ_dϕ(r₁, r₂)
            x₁, y₁, z₁ = r₁
            x₂, y₂, z₂ = r₂
            return -y₁*x₂ + x₁*y₂
        end
        return SVector{4}([
            dcosσ_dθ(r₁, r₂),
            dcosσ_dϕ(r₁, r₂),
            dcosσ_dθ(r₂, r₁),
            dcosσ_dϕ(r₂, r₁)
           ])
    end
    # Numerical derivative using Dual numbers
    dual_deriv1(pts) = gradient(p -> cosdistance(p...), pts)

    @test anal_deriv1(pts) ≈ @inferred dual_deriv1(pts)
end

end
