module Healpix
    using CMB.Healpix
    using Base.Test

    # Read off ring beginning/ending for Nside=4 map from Figure 5, Gorski
    # et al (2004).
    p_start_nside4 = [0; 4;12;24;40;56;72; 88;104;120;136;152;168;180;188]
    p_end_nside4   = [3;11;23;39;55;71;87;103;119;135;151;167;179;187;191]

    # Known ring start and ends for each ring
    @test all(pix2ring_ring.(4, p_start_nside4) .== 1:15)
    @test all(pix2ring_ring.(4, p_end_nside4)   .== 1:15)

    # Start of each ring is always 1
    @test all(pix2ringidx_ring.(4, p_start_nside4) .== 1)
    # End of rings vary based on whether in a cap or the equatorial belt
    @test all(pix2ringidx_ring.(4, p_end_nside4) .== [4;8;12;16ones(Int,9);12;8;4])

    # Make sure pixels at starts and ends of rings agree on their z value
    @test all(pix2z_ring.(4,p_start_nside4) .== pix2z_ring.(4,p_end_nside4))

    # Type inference doesn't handle the nside -> Val{nside} transformation
    # without losing type stability, so these currently fail.
    @test_broken @inferred pix2ring_ring(4, 0)
    @test_broken @inferred pix2ringidx(4, 0)
    @test_broken @inferred pix2z_ring(4, 0)
end

