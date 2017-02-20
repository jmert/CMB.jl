module Healpix
    using CMB.Healpix
    using Base.Test

    # Known ring numbers for start- and end-of-ring pixels in Nside=4 map.
    # Read of off Figure 5, Gorski et al (2004).
    function pix2ring_ring_nside4()
        nside = 4
        p_start = [0, 4,12,24,40,56,72, 88,104,120,136,152,168,180,188]
        p_end   = [3,11,23,39,55,71,87,103,119,135,151,167,179,187,191]
        @test all(pix2ring_ring.(nside, p_start) .== 1:15)
        @test all(pix2ring_ring.(nside, p_end)   .== 1:15)
    end

    function runtests()
        pix2ring_ring_nside4()
    end
end

