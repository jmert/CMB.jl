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
pix = first(pix) # healpy returns a 0-dimensional array
pixind = first(find(obspix .== pix))
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

Make the reddened spectrum...
```@example guide
ℓ = collect(0:700)
Cl = @. (2π / (ℓ*(ℓ+1)))
Cl[1] = 0.0 # division by zero, so get rid of the NaN
nothing # hide
```
Make the input spectrum be a TT + BB (no EE) spectrum:
```@example guide
spec = hcat(Cl, zeros(size(Cl)), Cl)
```

The [`pixelcovariance`](@ref) function makes use of the angular separation and
bearing angles to compute the pixel-pixel covariance terms according to
Tegmark, et. al. and sum over the input spectrum. It then returns terms
appropriate for a single column of each of the pixel-pixel covariance block
submatrices.

Each of these block columns can be visualized as a map:
```@example guide
covcol = PixelCovariance.pixelcovariance(512, obspix, pixind, spec)
fig = figure(figsize=(16,16))
plot_healpix_map(expand_obspix(covcol[:TT]); sub=(3,3,1), cmap="magma", title="TT pixel-pixel covariance")
plot_healpix_map(expand_obspix(covcol[:QQ]); sub=(3,3,5), cmap="magma", title="QQ pixel-pixel covariance")
plot_healpix_map(expand_obspix(covcol[:QU]); sub=(3,3,6), cmap="magma", title="QU pixel-pixel covariance")
plot_healpix_map(expand_obspix(covcol[:UQ]); sub=(3,3,8), cmap="magma", title="UQ pixel-pixel covariance")
plot_healpix_map(expand_obspix(covcol[:UU]); sub=(3,3,9), cmap="magma", title="UU pixel-pixel covariance")

savefig("center_pix_cov.png"); close(gcf()) # hide
nothing # hide
```
![](center_pix_cov.png)

