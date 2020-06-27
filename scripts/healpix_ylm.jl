using Base.Threads
using CMB
using CMB.SphericalHarmonics: synthesize_ring
using PyPlot, PyCall, Random

const healpy = pyimport("healpy")

nside = 512
lmax = 1000 #3nside - 1

Random.seed!(2)
alms = zeros(ComplexF64, lmax, lmax)
for l in 1:size(alms,1)
    alms[l,1] = randn() ./ l
    alms[l,2:end] .= complex.(randn.(), randn.()) ./ (l*sqrt(2))
end
llim = 100 #min(2nside, lmax)
alms = alms[1:llim, 1:llim]

function synthesize_healpix(nside, alms)
    syn = zeros(nside2npix(nside))
    @threads for rn in 1:2nside
        # rn counts rings from the north pole; rs counts from the south pole
        rs = 4nside - rn
        if rn < nside # polar caps
            npix = 4rn
            pixfirstn = nside2npixcap(rn)
            pixfirsts = nside2npix(nside) - pixfirstn - npix
        else # equatorial belt
            npix = 4nside
            pixfirstn = nside2npixcap(nside) + 4nside * (rn - nside)
            pixfirsts = pixfirstn + 4nside * (rs - rn)
        end
        θ, ϕ = pix2ang(nside, pixfirstn)

        rings = synthesize_ring(alms, θ, ϕ, npix, Val(true))
        syn[pixfirstn .+ (1:npix)] .= rings[1]
        syn[pixfirsts .+ (1:npix)] .= rings[2]
    end
    return syn
end
syn = synthesize_healpix(nside, alms)

# To match cartopy.Orthographic(30, 30), require arguments to be
#     rot = (30-180, 30), flip = "geo", half_sky = true
healpy.orthview(syn, rot = (30-180, 30), flip = "geo",
                cmap = "RdBu_r", fig = 1)
