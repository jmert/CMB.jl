using FITSIO
using CMB

# For plotting, we're going to use Julia's interface to pyplot (to make sure the environment
# is setup correctly) and Python's healpy (to actually interact with HEALPix maps).
using PyPlot
using PyCall
@pyimport healpy

"""
    plot_map(hmap; fig=nothing, cmap="hot", kwds...)

Makes a [relatively stupid] plot of the given HEALPix map, assuming that the map covers the
BICEP/Keck field.

This *should* be updated in the future to match BK plotting styles more, but right now, this
is simplest since it reuses existing HEALPix plotting machinery (available in Python).
"""
function plot_map(hmapin; fig=nothing, cmap="hot", kwds...)
    # healpy complains about the NaN values, so overwrite NaN with the HEALPix sentinel value
    # UNSEEN
    hmap = copy(hmapin)
    hmap[isnan.(hmap)] .= healpy.UNSEEN
    if fig === nothing
        # Default to something relatively appropriately sized
        hfig = PyPlot.figure(figsize=[8,4])
    end
    # Plot the figure
    hax = healpy.gnomview(hmap, fig=hfig[:number], rot=[0, -60],
                          reso=15, xsize=236, ysize=120,
                          cmap=cmap, kwds...)
    return hax
end

# Download a sample mask, if not already exists
if ~isfile("bk14_mask_cel_n0512.fits")
    run(`curl -kLO http://bicepkeck.org/BK14_datarelease/bk14_mask_cel_n0512.fits`);
end

# Read the apodization mask in to memory and plot it for visualization
hmask = FITS("bk14_mask_cel_n0512.fits");
apmask = read(hmask[2], "AP_MASK");
if false
    plot_map(apmask)
end

# Figure out which pixels must be computed upon. These are all the non-NaN pixels
obsmask = (~).(isnan.(apmask));
if false
    plot_map(float(obsmask), cmap="gray")
end

# Get the complete pixel list of observed pixels
obspix = find(obsmask)

# Choose a pixel at the center of the map
#   Currently retreived using healpy since I haven't coded up
#   the angle -> index functions yet.
pix = healpy.ang2pix(512, deg2rad(90-(-57.5)), 0.0)

# Get the distance between this pixel and every other pixel in the map
(θ₀,ϕ₀) = (pix2theta(512, pix), pix2phi(512, pix))
θ = pix2theta.(512, obspix);
ϕ = pix2phi.(512, obspix);
σ = cosdistance.(θ₀, ϕ₀, θ, ϕ);

# For plotting, make a dummy array which we can just keep filling.
tmp = Vector{Float64}(nside2npix(512));
tmp .= NaN;

if false
    # Plot the distance between pixels
    tmp[obspix] .= σ;
    plot_map(tmp, cmap="hot")
end

# To compute the polarized covariance, we need to know the angle of the great
# circles connecting each pair of points with respect to the meridian.
αij = bearing.(θ₀, ϕ₀, θ, ϕ);
αji = bearing.(θ, ϕ, θ₀, ϕ₀);
# For computation, we actually use cosines/sines of twice the angle (spin
# two field)
cij = cos.(2.*αij);
sij = sin.(2.*αij);
cji = cos.(2.*αji);
sji = sin.(2.*αji);

if false
    # Show samples of the angles
    tmp[obspix] .= sji;
    plot_map(tmp, cmap="jet")
end

# Arbitrarily choose a spectrum. Choose to have the spectrum be 1μK² at ℓ==80
# in ℓ(ℓ+1)C_ℓ/2π
#
# (Optionally) Then normalized away the covariances' normalization
# factor as well.
ℓ = collect(0:700);
Cl = @. (2π / (ℓ*(ℓ+1))) # * (4π / (2ℓ+1));
Cl[1] = 0.0; # division by zero, so get rid of the NaN

if false
    semilogy(ℓ, @. ℓ*(ℓ+1)*Cl/(2π) )
end

# Setup to comput the TT covariance
coeff = PixelCovariance.PixelCovarianceCoeff{Float64}(700);
distP = Vector{Float64}(701);
if false
    LegendreP!(coeff.λ, distP, 700, 0, σ[1])
    plot(distP)
end

# Loop over all pixels, summing the Legendre polynomials as we go
covTT = similar(σ);
for (i,z) in enumerate(σ)
    LegendreP!(coeff.λ, distP, 700, 0, z)
    # Sum over all ℓ in this pixel
    val = zero(eltype(distP))
    #for ll in length(distP):-1:1
    for ll in collect(55:90)
        val += coeff.η[ll] * Cl[ll] * distP[ll]
    end
    covTT[i] = val
end

if true
    # Plot the pixel-pixel covariance
    tmp[obspix] .= covTT
    plot_map(tmp, cmap="magma")
end


