module Sphere
using Compat.Test
using ..CMBTests: NumTypes
using CMB.Sphere, StaticArrays

@testset "Antipodes (angles, $T)" for T in NumTypes
    npole = (T(0.0),  T(1.0))
    spole = (T(π), T(2.0))

    # North and South Pole antipodes
    @test @inferred bearing(    npole..., spole...) == T(0.0)
    @test @inferred bearing(    spole..., npole...) == T(0.0)
    @test @inferred bearing2(   npole..., spole...) == (T(1.0), T(0.0))
    @test @inferred bearing2(   spole..., npole...) == (T(1.0), T(0.0))
    @test @inferred distance(   npole..., spole...) == T(π)
    @test @inferred cosdistance(npole..., spole...) == T(-1.0)

    # North and North Pole
    @test @inferred bearing(    npole..., npole...) == T(0.0)
    @test @inferred bearing2(   npole..., npole...) == (T(1.0), T(0.0))
    @test @inferred distance(   npole..., npole...) == T(0.0)
    @test @inferred cosdistance(npole..., npole...) == T(1.0)

    # South and South Pole
    @test @inferred bearing(    spole..., spole...) == T(0.0)
    @test @inferred bearing2(   spole..., spole...) == (T(1.0), T(0.0))
    @test @inferred distance(   spole..., spole...) == T(0.0)
    @test @inferred cosdistance(spole..., spole...) == T(1.0)
end

@testset "Antipodes (vectors, $T)" for T in NumTypes
    npole = @SVector T[0.0, 0.0,  1.0]
    spole = @SVector T[0.0, 0.0, -1.0]

    # North and South Pole antipodes
    @test @inferred bearing(    npole, spole) == T(0.0)
    @test @inferred bearing(    spole, npole) == T(0.0)
    @test @inferred bearing2(   npole, spole) == (T(1.0), T(0.0))
    @test @inferred bearing2(   spole, npole) == (T(1.0), T(0.0))
    @test @inferred distance(   npole, spole) == T(π)
    @test @inferred cosdistance(npole, spole) == T(-1.0)

    # North and North Pole
    @test @inferred bearing(    npole, npole) == T(0.0)
    @test @inferred bearing2(   npole, npole) == (T(1.0), T(0.0))
    @test @inferred distance(   npole, npole) == T(0.0)
    @test @inferred cosdistance(npole, npole) == T(1.0)

    # South and South Pole
    @test @inferred bearing(    spole, spole) == T(0.0)
    @test @inferred bearing2(   spole, spole) == (T(1.0), T(0.0))
    @test @inferred distance(   spole, spole) == T(0.0)
    @test @inferred cosdistance(spole, spole) == T(1.0)
end

@testset "Pole to Equator (angles, $T)" for T in NumTypes
    # Pole with arbitrary points around equator. Distance at the equator is easy to know
    # analytically.
    npole = (T(0.0), T(1.0))
    atol = max(eps(1π), eps(T(π)))
    for ii=1:7
        local pi = T===BigFloat ? T(π) : 1.0π
        eqpnt = (T(pi/2), T(2pi*ii/7))
        @test @inferred bearing(    npole..., eqpnt...) == T(0.0)
        @test @inferred bearing2(   npole..., eqpnt...) == (T(1.0), T(0.0))
        @test @inferred(distance(   npole..., eqpnt...)) ≈ T(π/2) atol=atol
        @test @inferred(cosdistance(npole..., eqpnt...)) ≈ T(0.0) atol=atol
    end
end

@testset "Pole to Equator (vectors, $T)" for T in NumTypes
    # Pole with arbitrary points around equator. Distance at the equator is easy to know
    # analytically, so vary longitudes to not nice multiples of π.
    npole = @SVector T[0.0, 0.0, 1.0]
    atol = max(eps(1π), eps(T(π)))
    for ii=1:7
        local pi = T===BigFloat ? T(π) : 1.0π
        eqpnt = @SVector T[cos(2pi*ii/7), sin(2pi*ii/7), 0.0]
        @test @inferred bearing(    npole, eqpnt) == T(0.0)
        @test @inferred bearing2(   npole, eqpnt) == (T(1.0), T(0.0))
        @test @inferred(distance(   npole, eqpnt)) ≈ T(π/2) atol=atol
        @test @inferred(cosdistance(npole, eqpnt)) ≈ T(0.0) atol=atol
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
    @test bearing2(r₁, r₂) == (1.0, 0.0)
end

end

