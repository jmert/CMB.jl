module Sphere

using CMB.Sphere
using Base.Test

const npole = (0.0,  1.0)
const spole = (1.0π, 2.0)

# North and South Pole antipodes
@test bearing(    npole..., spole...) == 0.0
@test bearing(    spole..., npole...) == 0.0
@test distance(   npole..., spole...) == 1π
@test cosdistance(npole..., spole...) == -1.0

# North and North Pole
@test bearing(    npole..., npole...) == 0.0
@test distance(   npole..., npole...) == 0.0
@test cosdistance(npole..., npole...) == 1.0

# South and South Pole
@test bearing(    spole..., spole...) == 0.0
@test distance(   spole..., spole...) == 0.0
@test cosdistance(spole..., spole...) == 1.0

# Pole with arbitrary points around equator. Distance at the equator is easy to know
# analytically, so vary longitudes to not nice multiples of π.
for ii=1:7
  @test bearing(    npole..., π/2, 2π*ii/7) == 0.0
  @test distance(   npole..., π/2, 2π*ii/7) ≈ π/2
  @test cosdistance(npole..., π/2, 2π*ii/7) ≈ 0.0 atol=eps(1π)
end

end

