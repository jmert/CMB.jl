using FITSIO
include("CMB.jl"); using CMB

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
function plot_map(hmap; fig=nothing, cmap="hot", kwds...)
    # healpy complains about the NaN values, so overwrite NaN with the HEALPix sentinel value
    # UNSEEN
    hmap[isnan.(hmap)] .= healpy.UNSEEN
    if fig === nothing
        # Default to something relatively appropriately sized
        hfig = PyPlot.figure(figsize=[8,4])
    end
    # Plot the figure
    hax = healpy.gnomview(hmap, fig=hfig[:number], rot=[0, -60],
                           reso=15, xsize=236, ysize=120,
                           min=0, max=1, cmap=cmap, kwds...)
    return hax
end

# Download a sample mask, if not already exists
if ~isfile("bk14_mask_cel_n0512.fits")
    run(`curl -kLO http://bicepkeck.org/BK14_datarelease/bk14_mask_cel_n0512.fits`);
end

# Read the apodization mask in to memory and plot it for visualization
hmask = FITS("bk14_mask_cel_n0512.fits");
apmask = read(hmask[2], "AP_MASK");
plot_map(apmask)

# Figure out which pixels must be computed upon. These are all the non-NaN pixels
obsmask = ~.(isnan.(apmask));
plot_map(float(obsmask), cmap="gray")

# Get the complete pixel list of observed pixels
obspix = find(obsmask)


