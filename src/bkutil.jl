module BKUtil

export plot_healpix_map

using ..Healpix: UNSEEN
using ..Mapping
using PyPlot
using PyCall

const healpy = PyCall.PyNULL()
const m_bicep = BicepMapDefn()

function __init__()
    copy!(healpy, pyimport("healpy"))
end

function plot_healpix_map(m::ECPMapPatchDefn, hmapin::Vector; kwds...)
    # healpy complains about the NaN values, so overwrite NaN with the HEALPix sentinel value
    # UNSEEN
    hmap = copy(hmapin)
    hmap[isnan.(hmap)] .= UNSEEN

    # Working in Julia, we may be given a PyPlot.Figure object instead of the figure number
    # which healpy is actually expecting. Automatically convert these if found.
    _conv_fig(p) = ( (k,v)=p; (v isa PyPlot.Figure) ? (k,v[:number]) : (k,v) )
    map!(_conv_fig, kwds, kwds)

    # Plot the figure
    healpy[:cartview](hmap; lonra=[m.lx,m.hx], latra=[m.ly,m.hy],
                    kwds...)
    ax = gca()
    # Fix the aspect ratio so that it looks much more BK like
    ax[:set_aspect](m.xdos / m.ydos)
end

plot_healpix_map(hmapin::Vector; kwds...) = plot_healpix_map(m_bicep, hmapin; kwds...)

end
