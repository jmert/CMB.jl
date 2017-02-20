module Healpix
    using CMB.Healpix
    using Base.Test

    # Read off ring beginning/ending for Nside=4 map from Figure 5, Gorski
    # et al (2004).
    p_start_nside4 = [0, 4,12,24,40,56,72, 88,104,120,136,152,168,180,188]
    p_end_nside4   = [3,11,23,39,55,71,87,103,119,135,151,167,179,187,191]

    # Known ring start and ends for each ring
    @test all(pix2ring_ring.(4, p_start_nside4) .== 1:15)
    @test all(pix2ring_ring.(4, p_end_nside4)   .== 1:15)
end

