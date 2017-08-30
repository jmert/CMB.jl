# Sample Computation

```@example guide
using CMB
using CMB.BKUtils: plot_healpix_map
using FITSIO
using PyPlot
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
plot_healpix_map(apmask; min=0, cmap="hot", title="BK14 Apodization Mask")
savefig("bk14_apmask.png"); close(gcf()) # hide
nothing # hide
```
![](bk14_apmask.png)

```@example guide
obsmask = (~).(isnan.(apmask))
plot_healpix_map(float(obsmask); cmap="gray", title="Observed Pixel Mask")
savefig("bk14_binary_mask.png"); close(gcf()) # hide
nothing # hide
```
![](bk14_binary_mask.png)

```@example guide
obspix = find(obsmask);
@eval function expand_obspix(v)
    f = Vector{Float64}($(length(obsmask)))
    fill!(f, Healpix.UNSEEN)
    f[$obspix] .= v
    return f
end
nothing # hide
```

```@example guide
using PyCall
@pyimport healpy
pix = healpy.ang2pix(512, deg2rad(90-(-57.5)), 0.0)
```

```@example guide
(θ₀,ϕ₀) = (pix2theta(512, pix), pix2phi(512, pix))
θ = pix2theta.(512, obspix)
ϕ = pix2phi.(512, obspix)
σ = cosdistance.(θ₀, ϕ₀, θ, ϕ)
```

```@example guide
plot_healpix_map(expand_obspix(σ), cmap="magma", title="Cosine of angular pixel separation")
savefig("center_pix_distance.png"); close(gcf()) # hide
```
![](center_pix_distance.png)

```@example guide
αij = bearing.(θ₀, ϕ₀, θ, ϕ)
αji = bearing.(θ, ϕ, θ₀, ϕ₀)
cij = cos.(2.*αij)
sij = sin.(2.*αij)
cji = cos.(2.*αji)
sji = sin.(2.*αji)
fig = figure(figsize=(12,8))
plot_healpix_map(expand_obspix(cij); fig=fig, cmap="magma", sub=(2,2,1), title="cos(2α_ij)", cbar=false)
plot_healpix_map(expand_obspix(sij); fig=fig, cmap="magma", sub=(2,2,2), title="sin(2α_ij)", cbar=false)
plot_healpix_map(expand_obspix(cji); fig=fig, cmap="magma", sub=(2,2,3), title="cos(2α_ji)")
plot_healpix_map(expand_obspix(sji); fig=fig, cmap="magma", sub=(2,2,4), title="sin(2α_ji)")
savefig("center_pix_azimuths.png"); close(gcf()) # hide
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

plot_healpix_map(expand_obspix(covTT); cmap="magma", title="TT pixel-pixel covariance")
savefig("center_pix_covTT.png"); close(gcf()) # hide
nothing # hide
```
![](center_pix_covTT.png)

