using SnoopCompile
using CMB
using CMB.Sphere: x̂, ŷ, cartvec, colataz, latlon

# Really fast functions: no minimum time
inf_timing = @snoopi tmin = 0.0 begin
    ###
    ### Spherical functions

    #   - one angle pair
    for f in (cartvec, colataz, latlon)
        f((0.1, 0.2))
        f(0.1, 0.2)
    end
    #   - one pointing vector
    for f in (colataz, latlon)
        f(x̂)
    end
    #   - two points as angle pairs
    for f in (bearing, bearing2, distance, cosdistance)
        f(0.1, 0.1, 0.2, 0.2)
    end
    #   - two points as vectors
    for f in (bearing, bearing2, distance, cosdistance)
        f(x̂, ŷ)
    end
    reckon(x̂, 0.1, 0.1)
    reckon(0.1, 0.1, 0.1, 0.1)

    ###
    ### HEALPix pixel indexing:

    #   - Single parameter functions
    for f in (nside2npix, nside2npixcap, nside2npixequ, nside2pixarea,
              nside2nring, nring2nside, checkhealpix)
        f(4)
    end
    #   - nside + pixel index functions
    for f in (pix2ring, pix2ringidx, pix2z, pix2theta, pix2phi, pix2ang, pix2vec)
        f(4, 1)
    end
    ang2pix(4, 1.0, 1.0)
    vec2pix(4, CMB.Sphere.x̂)

end

# More expensive functions
next_inf = @snoopi tmin = 0.005 begin
    ###
    ### Pixel-pixel covariance
    Fweights!(LegendreUnitNorm(), zeros(5, 4), 4, 0.1)
end

append!(inf_timing, next_inf)
pc = SnoopCompile.parcel(inf_timing)
pc = filter!(p -> p.first == :CMB, pc)
sort!(pc[:CMB])

outfile = joinpath(tempdir(), "precompile")
println("Saving precompile scripts to $outfile")
SnoopCompile.write(outfile, pc)
