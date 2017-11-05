module Sphere
using Compat.Test
using CMB.Sphere, StaticArrays

@testset "Sphere antipodes (angles)" begin
    npole = (0.0,  1.0)
    spole = (1.0π, 2.0)

    # North and South Pole antipodes
    @test bearing(    npole..., spole...) == 0.0
    @test bearing(    spole..., npole...) == 0.0
    @test bearing2(   npole..., spole...) == (1.0, 0.0)
    @test bearing2(   spole..., npole...) == (1.0, 0.0)
    @test distance(   npole..., spole...) == 1π
    @test cosdistance(npole..., spole...) == -1.0

    # North and North Pole
    @test bearing(    npole..., npole...) == 0.0
    @test bearing2(   npole..., npole...) == (1.0, 0.0)
    @test distance(   npole..., npole...) == 0.0
    @test cosdistance(npole..., npole...) == 1.0

    # South and South Pole
    @test bearing(    spole..., spole...) == 0.0
    @test bearing2(   spole..., spole...) == (1.0, 0.0)
    @test distance(   spole..., spole...) == 0.0
    @test cosdistance(spole..., spole...) == 1.0
end

@testset "Sphere antipodes (vectors)" begin
    npole = @SVector [0.0, 0.0,  1.0]
    spole = @SVector [0.0, 0.0, -1.0]

    # North and South Pole antipodes
    @test bearing(    npole, spole) == 0.0
    @test bearing(    spole, npole) == 0.0
    @test bearing2(   npole, spole) == (1.0, 0.0)
    @test bearing2(   spole, npole) == (1.0, 0.0)
    @test distance(   npole, spole) == 1π
    @test cosdistance(npole, spole) == -1.0

    # North and North Pole
    @test bearing(    npole, npole) == 0.0
    @test bearing2(   npole, npole) == (1.0, 0.0)
    @test distance(   npole, npole) == 0.0
    @test cosdistance(npole, npole) == 1.0

    # South and South Pole
    @test bearing(    spole, spole) == 0.0
    @test bearing2(   npole, npole) == (1.0, 0.0)
    @test distance(   spole, spole) == 0.0
    @test cosdistance(spole, spole) == 1.0
end

@testset "Sphere Pole to Equator (angles)" begin
    # Pole with arbitrary points around equator. Distance at the equator is easy to know
    # analytically.
    npole = (0.0, 1.0)
    for ii=1:7
        eqpnt = (π/2, 2π*ii/7)
        @test bearing(    npole..., eqpnt...) == 0.0
        @test bearing2(   npole..., eqpnt...) == (1.0, 0.0)
        @test distance(   npole..., eqpnt...) ≈ π/2
        @test cosdistance(npole..., eqpnt...) ≈ 0.0 atol=eps(1π)
    end
end

@testset "Sphere Pole to Equator (vectors)" begin
    # Pole with arbitrary points around equator. Distance at the equator is easy to know
    # analytically, so vary longitudes to not nice multiples of π.
    npole = @SVector [0.0, 0.0, 1.0]
    for ii=1:7
        eqpnt = @SVector [cos(2π*ii/7), sin(2π*ii/7), 0.0]
        @test bearing(    npole, eqpnt) == 0.0
        @test bearing2(   npole, eqpnt) == (1.0, 0.0)
        @test distance(   npole, eqpnt) ≈ π/2
        @test cosdistance(npole, eqpnt) ≈ 0.0 atol=eps(1.0)
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

