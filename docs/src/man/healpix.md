# `HEALPix` Pixelization

```@contents
Pages = ["healpix.md"]
Depth = 2
```

The [`Healpix`](@ref) module implements a selection of functions for interacting with
the **H**ierarchical **E**qual-**A**rea Iso-**L**atitude **Pix**elization (`HEALPix`) as
described in [Górski et al. (2005)](@ref bib-healpix).
The emphasis has been on working with pixel index and spherical coordinate conversions
for the RING-ordered scheme only.
If a complete `HEALPix` implementation is required, try
[LibHealpix.jl](https://github.com/mweastwood/LibHealpix.jl) which provides Julia bindings
to the official C++ library or try using the Python
[healpy](https://github.com/healpy/healpy) interface via
[PyCall.jl](https://github.com/JuliaPy/PyCall.jl).
Additionally, the [Healpix.jl](https://github.com/ziotom78/Healpix.jl) package aims for a
complete native Julia reimplementation.
``\newcommand{\Nside}{N_\mathrm{side}}``

## Basic pixelization properties

A core defining attribute of the `HEALPix` map is the ``\Nside`` parameter.
The minimum valid value is ``\Nside = 1``, representing the coarsest pixelization of the
sphere which can be represented by `HEALPix`.
The ``\Nside`` then increases by factors of two—denoting ever finer resolutions—with
the number of pixels covering the sphere being a function of ``\Nside`` alone.
[`CMB.Healpix.nside2npix`](@ref) returns the number of pixels in a given map:
```jldoctest healpix
julia> using CMB.Healpix

julia> nside = 4;

julia> nside2npix(nside)
192
```
Because all pixels are of equal area and the number of pixels is derived from only the
``\Nside``, the pixel surface area must be as well.
For convenience, [`CMB.Healpix.nside2pixarea`](@ref) is provided and is equivalent to the
simple from-scratch calculation (in steradians).
```jldoctest healpix
julia> 4π / nside2npix(nside)
0.06544984694978735

julia> nside2pixarea(nside)
0.06544984694978735
```
The number of iso-latitude rings is also a function of only ``\Nside`` and calculated by
[`CMB.Healpix.nside2nring`](@ref):
```jldoctest healpix
julia> nside2nring(nside)
15
```

## Working with pixel indices

The pixels are enumerated as 0-indexed integers from west to east along the iso-latitude
rings, from north to south.
For example, pixel 0 is the first pixel in the first ring, and pixel 103 is the sixteenth
pixel in the eighth ring for an ``\Nside = 4`` map:
```jldoctest healpix
julia> pix = (0, 103);

julia> pix2ring.(nside, pix)    # Ring
(1, 8)

julia> pix2ringidx.(nside, pix) # Index in ring
(1, 16)
```

!!! note
    Be careful to note that pixels are 0-indexed, but the rings and indices within a ring
    are 1-indexed.

The `HEALPix` grid is symmetric about equator, with the equatorial ring considered part of
the northern hemisphere by convention.
Membership as part of the northern or southern hemisphere can be tested with the
[`CMB.Healpix.isnorth`](@ref) and [`CMB.Healpix.issouth`](@ref) functions, respectively.
Pixel 103 is actually the last pixel in the northern hemisphere, so
```jldoctest healpix
julia> isnorth(nside, 103)
true

julia> isnorth(nside, 104)
false

julia> issouth(nside, 104)
true
```
In fact, each hemisphere is further composed of a so-called *polar cap* and *equatorial
belt* region of pixels (a property derived from the mathematical details of the `HEALPix`
grid's definition).
According to the ring-ordered definition, pixel 0 should be in the polar cap (tested via
[`CMB.Healpix.iscap`](@ref)), while pixel 103 in the equatorial ring is expected to be
part of the equitorial belt (tested via [`CMB.Healpix.isequbelt`](@ref)).
```jldoctest healpix
julia> iscap.(nside, pix)
(true, false)

julia> isequbelt.(nside, pix)
(false, true)
```
Membership in a particular hemisphere's polar cap or equatorial belt is accomplished with
variants inserting `north` and `south` into the function names, i.e. polar caps are
distinguished by [`CMB.Healpix.isnorthcap`](@ref) and [`CMB.Healpix.issouthcap`](@ref),
and the halves of the equatorial belt are distinguished by
[`CMB.Healpix.isnorthequbelt`](@ref) and [`CMB.Healpix.issouthequbelt`](@ref).
```jldoctest healpix
julia> pix = (0, 103, 104, 191);

julia> isnorthcap.(nside, pix)
(true, false, false, false)

julia> isnorthequbelt.(nside, pix)
(false, true, false, false)

julia> issouthequbelt.(nside, pix)
(false, false, true, false)

julia> issouthcap.(nside, pix)
(false, false, false, true)
```

## Working with spherical coordinates

Up to now, all the features shown have concerned working with properties of the
pixelization scheme, but the utility of the `HEALPix` grid is its ability to describe the
surface of a sphere.
Using spherical coordinates is more useful and more natural for more algorithms than the
`HEALPix`-specific indexing scheme.

The first method of describing the location of a particular `HEALPix` pixel is as a
colatitude/azimuth pair of angles on the surface of the sphere identifying the pixel
center.
Colatitude measures the angle (in radians) south of the North Pole, and azimuth measures
the angle (in radians) east of the Prime Meridian.
To get the colatitude, use [`CMB.Healpix.pix2theta`](@ref),
```jldoctest healpix
julia> pix2theta(nside, 103)
1.5707963267948966
```
and to get the azimuth, use [`CMB.Healpix.pix2phi`](@ref)
```jldoctest healpix
julia> pix2phi(nside, 103)
6.086835766330224
```
(both named to follow the mathematical convention that colatitude/azimuth pairs in
spherical coordinates are the variable pair ``(θ, ϕ)``).
When the coordinate pair is required, the method [`CMB.Healpix.pix2ang`](@ref) returns a
2-tuple with the coordinates:
```jldoctest healpix
julia> pix2ang(nside, 103)
(1.5707963267948966, 6.086835766330224)

julia> pix2ang(nside, 103) .|> rad2deg
(90.0, 348.75)
```

The other common way to represent coordinates on the sphere is via unit vectors.
The corresponding vector for a given pixel is retrieved with
[`CMB.Healpix.pix2vec`](@ref).
```jldoctest healpix
julia> pix2vec(nside, 103)
3-element SVector{3,Float64}:
  0.980785
 -0.19509
  0.0
```
where the elements correspond to the typical ``(x, y, z)`` right-handed coordinates with
the positive ``z``-axis passing through the North Pole and the positive ``x``-axis
passing through the Prime Meridian.

In reverse, converting an arbitrary spherical coordinate to a pixel index... *...TO BE
IMPLEMENTED...*
