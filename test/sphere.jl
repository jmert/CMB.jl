module Sphere
using Test
using LinearAlgebra, StaticArrays
using ..CMBTests: NumTypes
using CMB.Sphere
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

end

