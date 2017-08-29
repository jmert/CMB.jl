```@setup guide
using CMB
using FITSIO
using PyPlot
using PyCall
@pyimport healpy
```

# Sample Computation

```@example guide
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
                          cmap=cmap; kwds...)
    return hax
end
nothing # hide
```

```@setup guide
# Download a sample mask, if not already exists
cd() do
    if ~isfile("bk14_mask_cel_n0512.fits")
        run(`curl -kLO http://bicepkeck.org/BK14_datarelease/bk14_mask_cel_n0512.fits`)
    end
end
# Get FITS module compiled to avoid dumping deprecation warnings to the
# rendered output
hmask = FITS("~/bk14_mask_cel_n0512.fits")
apmask = read(hmask[2], "AP_MASK")
```
```@example guide
hmask = FITS("~/bk14_mask_cel_n0512.fits")
apmask = read(hmask[2], "AP_MASK")
plot_map(apmask; min=0)
savefig("bk14_apmask.png") # hide
nothing # hide
```
![](bk14_apmask.png)

```@example guide
obsmask = (~).(isnan.(apmask))
plot_map(float(obsmask); cmap="gray")
savefig("bk14_binary_mask.png") # hide
nothing # hide
```
![](bk14_binary_mask.png)

```@example guide
obspix = find(obsmask);
nothing # hide
```

```@example guide
pix = healpy.ang2pix(512, deg2rad(90-(-57.5)), 0.0)
```

```@example guide
(θ₀,ϕ₀) = (pix2theta(512, pix), pix2phi(512, pix))
θ = pix2theta.(512, obspix)
ϕ = pix2phi.(512, obspix)
σ = cosdistance.(θ₀, ϕ₀, θ, ϕ)
```

```@example guide
tmp = Vector{Float64}(nside2npix(512));
tmp .= NaN;
nothing # hide
```

```@example guide
tmp[obspix] .= σ
plot_map(tmp, cmap="hot")
savefig("center_pix_distance.png") # hide
```
![](center_pix_distance.png)

```@example guide
αij = bearing.(θ₀, ϕ₀, θ, ϕ)
αji = bearing.(θ, ϕ, θ₀, ϕ₀)
cij = cos.(2.*αij)
sij = sin.(2.*αij)
cji = cos.(2.*αji)
sji = sin.(2.*αji)
tmp[obspix] .= sji
plot_map(tmp, cmap="jet")
savefig("center_pix_azimuths.png") # hide
```
![](center_pix_azimuths.png)

```@example guide
ℓ = collect(0:700)
Cl = @. (2π / (ℓ*(ℓ+1)))
Cl[1] = 0.0 # division by zero, so get rid of the NaN
nothing # hide
```

```@example guide
# Setup to comput the TT covariance
coeff = PixelCovariance.PixelCovarianceCoeff{Float64}(700)
distP = Vector{Float64}(701)
nothing # hide
```

```@example guide
covTT = similar(σ)
for (i,z) in enumerate(σ)
    LegendreP!(coeff.λ, distP, 700, 0, z)
    # Sum over all ℓ in this pixel
    val = zero(eltype(distP))
    for ll in length(distP):-1:1
        val += coeff.η[ll] * Cl[ll] * distP[ll]
    end
    covTT[i] = val
end

tmp[obspix] .= covTT
plot_map(tmp, cmap="magma")
savefig("center_pix_covTT.png") # hide
nothing # hide
```
![](center_pix_covTT.png)

