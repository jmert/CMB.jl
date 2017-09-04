# Sample Computation

```@example guide
using CMB
using CMB.BKUtil: plot_healpix_map
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

The sines/cosines are used...
```math
\newcommand{\covF}[1]{F_\ell^{#1}(x)}
\newcommand{\TT}{\langle T_i T_j\rangle}
\newcommand{\TQ}{\langle T_i Q_j\rangle}
\newcommand{\TU}{\langle T_i U_j\rangle}
\newcommand{\QQ}{\langle Q_i Q_j\rangle}
\newcommand{\QU}{\langle Q_i U_j\rangle}
\newcommand{\UU}{\langle U_i U_j\rangle}
\begin{align}
    \TT &\equiv \sum_\ell C_\ell^{TT} \covF{00}
    \\
    \TQ &\equiv -\sum_\ell C_\ell^{TE} \covF{10}
    \\
    \TU &\equiv -\sum_\ell C_\ell^{BT} \covF{10}
    \\
    \QQ &\equiv \sum_\ell
        \left[ C_\ell^{EE} \covF{12} - C_\ell^{BB} \covF{22} \right]
    \\
    \UU &\equiv \sum_\ell
        \left[ C_\ell^{BB} \covF{12} - C_\ell^{EE} \covF{22} \right]
    \\
    \QU &\equiv \sum_\ell
        C_\ell^{EB} \left[ \covF{12} + \covF{22} \right]
\end{align}
```

```math
\begin{align}
    R(\alpha) &= \begin{bmatrix}
        1 & 0 & 0 \\
        0 &  \cos(2\alpha) & \sin(2\alpha) \\
        0 & -\sin(2\alpha) & \cos(2\alpha)
    \end{bmatrix}
    &
    M(z_{ij}) &= \begin{bmatrix}
        \TT & \TQ & \TU \\
        \TQ & \QQ & \QU \\
        \TU & \QU & \UU
    \end{bmatrix}
\end{align}
```

```math
\newcommand{\cij}{c_{ij}}
\newcommand{\sij}{s_{ij}}
\newcommand{\cji}{c_{ji}}
\newcommand{\sji}{s_{ji}}
\begin{align}
    R(\alpha_{ij}) &= \begin{bmatrix}
        1 & 0 & 0 \\
        0 &  \cij & \sij \\
        0 & -\sij & \cij
    \end{bmatrix}
    &
    R(\alpha_{ji})^T &= \begin{bmatrix}
        1 & 0 & 0 \\
        0 & \cji & -\sji \\
        0 & \sji &  \cji
    \end{bmatrix}
\end{align}
```

```math
\begin{align*}
    \begin{aligned}
    C_{ij} = R(\alpha_{ij}) M(z_{ij}) R(\alpha_{ji})^T =
        % Show the first two columns together first:
        &\left[ \begin{matrix}
            \TT \\
            \TQ\cij + \TU\sij \\
            -\TQ\sij + \TU\cij
        \end{matrix} \right.
    \\
        \qquad\ldots\qquad
        &\begin{matrix}
            \TQ\cji + \TU\sji \\
            \QQ\cij\cji + \QU(\cij\sji + \sij\cji) + \UU\sij\sji \\
            -\QQ\sij\cji + \QU(\cij\cji - \sij\sji) + \UU\cij\sji \\
        \end{matrix}
    \\
        \qquad\ldots\qquad
        &\left. \begin{matrix}
            -\TQ\sji + \TU\cji \\
            -\QQ\cij\sji + \QU(\cij\cji - \sij\sji) + \UU\sij\cji \\
            \QQ\sij\sji - \QU(\cij\sji + \sij\cji) + \UU\cij\cji
        \end{matrix} \right]
    \end{aligned}
\end{align*}
```

Make the reddened spectrum...
```@example guide
ℓ = collect(0:700)
Cl = @. (2π / (ℓ*(ℓ+1)))
Cl[1] = 0.0 # division by zero, so get rid of the NaN
nothing # hide
```
Make the input spectrum be a TT + BB (no EE) spectrum:
```@example guide
spec = zeros(Float64, 701, 6)
spec[:,1] .= Cl # TT
spec[:,3] .= Cl # EE
```

The [`pixelcovariance`](@ref) function makes use of the angular separation and
bearing angles to compute the pixel-pixel covariance terms according to
Tegmark, et. al. and sum over the input spectrum. It then returns terms
appropriate for a single column of each of the pixel-pixel covariance block
submatrices.

Each of these block columns can be visualized as a map:
```@example guide
cache = PixelCovarianceCache(512, 700, obspix, [:TT,:QQ,:QU,:UQ,:UU])
covmat = Matrix{Float64}(length(obspix), 9)
updatespectra!(cache, spec)
selectpixel!(cache, pixind)

pixelcovariance!(cache, covmat)
fig = figure(figsize=(16,16))
plot_healpix_map(expand_obspix(covmat[:,1]); sub=(3,3,1), cmap="magma", title="TT pixel-pixel covariance")
plot_healpix_map(expand_obspix(covmat[:,5]); sub=(3,3,5), cmap="magma", title="QQ pixel-pixel covariance")
plot_healpix_map(expand_obspix(covmat[:,8]); sub=(3,3,6), cmap="magma", title="QU pixel-pixel covariance")
plot_healpix_map(expand_obspix(covmat[:,6]); sub=(3,3,8), cmap="magma", title="UQ pixel-pixel covariance")
plot_healpix_map(expand_obspix(covmat[:,9]); sub=(3,3,9), cmap="magma", title="UU pixel-pixel covariance")

savefig("center_pix_cov.png"); close(gcf()) # hide
nothing # hide
```
![](center_pix_cov.png)

