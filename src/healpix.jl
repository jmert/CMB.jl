"""
A module of functions implementing function which interact with the HEALPix pixel
definitions. In most cases, only the RING ordering functions are being provided.

See "HEALPix: A Framework for High-Resolution Discretization and Fast Analysis of Data
Distributed on the Sphere" Górski, Hivon, & Banday et al (2005) ApJ 622:759–771
arXiv: astro-ph/0409513
"""
module Healpix

export
    npix2nside, nring2nside,
    nside2npix, nside2nring, nside2npixcap, nside2npixequ, nside2npixarea,
    isnorth, issouth,
    isnorthcap, issouthcap, iscap,
    isnorthequbelt, issouthequbelt, isequbelt,
    pix2ring, pix2ringidx, pix2z, pix2theta, pix2phi, pix2ang, pix2vec,
    UNSEEN

using StaticArrays

"""
    const UNSEEN = -1.6375e+30

Special value recognized by the libhealpix/healpy routines as an unobserved/masked
pixel.
"""
const UNSEEN = -1.6375e+30

Base.@pure npix2nside(npix::Integer)   = trunc(typeof(npix), sqrt(npix/12))
Base.@pure nring2nside(nring::Integer) = (nring + one(nring)) ÷ 4

Base.@pure nside2npix(nside::Integer)    = 12*nside*nside
Base.@pure nside2nring(nside::Integer)   =  4*nside - one(nside)
Base.@pure nside2npixcap(nside::Integer) =  2*nside*(nside - one(nside))
Base.@pure nside2npixequ(nside::Integer) =  2*nside*(3*nside + one(nside))
Base.@pure nside2pixarea(nside::Integer) = 4π / nside2npix(nside)

Base.@pure isnorth(nside, p) = p < nside2npixequ(nside)
Base.@pure issouth(nside, p) = p ≥ nside2npixequ(nside)

Base.@pure isnorthcap(nside, p) = p < nside2npixcap(nside)
Base.@pure issouthcap(nside, p) = p ≥ nside2npix(nside) - nside2npixcap(nside)
Base.@pure iscap(nside, p)      = isnorthcap(nside, p) | issouthcap(nside, p)

Base.@pure isnorthequbelt(nside, p) = isnorth(nside, p) & ~isnorthcap(nside, p)
Base.@pure issouthequbelt(nside, p) = issouth(nside, p) & ~issouthcap(nside, p)
Base.@pure isequbelt(nside, p)      = ~iscap(n, p)

# The function definitions are short enough that I wanted to keep them visually together
# without having the documentation break them up. Add documentation now

"""
    npix2nside(npix) -> Nside

Returns the equivalent `Nside` corresponding to the number of pixels `npix`. Note that no
validation is performed, so non-conformant values of `npix` will give non-conformant
Nside values.
""" npix2nside(npix)

"""
    nring2nside(nring) -> Nside

Returns the equivalent `Nside` corresponding to the number of iso-latitude rings `nring`.
Note that no validation is performed, so non-conformant values of `nring` will give
non-conformant Nside values.
""" nring2nside(npix)

"""
    nside2npix(nside) -> N

Returns the total number of pixels `N` in an `nside` HEALPix map.

N.B.: HEALPix pixel indexing is 0-based, so valid pixel values are in the range
0:(N-1).
""" nside2npix

"""
    nside2nring(nside) -> N

Returns the number of iso-latitude rings `N` in the `nside` HEALPix map.
""" nside2nring

"""
    nside2npixcap(nside) -> N

Returns the number of pixels `N` in the polar caps for the given `nside` HEALPix map.
""" nside2npixcap

"""
    nside2npixequ(nside) -> N

Returns the number of pixels `N` in the northern hemisphere, including the equatorial
ring, for the given `nside` HEALPix map.
""" nside2npixequ

"""
    nside2pixarea(nside) -> σ

Returns the surface area `σ` (in steradiands) of each pixel in the given `nside` HEALPix
map.
""" nside2pixarea

"""
    isnorthcap(nside, pix) -> TF

Test for whether the given pixel `pix` is in the northern polar cap for an `nside`
HEALPix map.
""" isnorthcap

"""
    issouthcap(nside, pix) -> TF

Test for whether the given pixel `pix` is in the southern polar cap for an `nside`
HEALPix map.
""" issouthcap

"""
    iscap(nside, pix) -> TF

Test for whether the given pixel `pix` is in either polar cap for an `nside` HEALPix
map.
""" iscap

"""
    isnorthequbelt(nside, pix) -> TF

Test for whether the given pixel `pix` is in the northern equatorial belt (including the
equatorial ring) for an `nside` HEALPix map.
""" isnorthequbelt

"""
    issouthequbelt(nside, pix) -> TF

Test for whether the given pixel `pix` is in the southern equatorial belt (excluding the
equatorial ring) for an `nside` HEALPix map.
""" issouthequbelt

"""
    isequbelt(nside, pix) -> TF

Test for whether the given pixel `pix` is in the equatorial belt for an `nside` HEALPix
map.
""" isequbelt

"""
    isnorth(nside, pix) -> TF

Test for whether the given pixel `pix` is in the northern hemisphere (including the
equatorial ring) for an `nside` HEALPix map.
""" isnorth

"""
    issouth(nside, pix) -> TF

Test for whether the given pixel `pix` is in the southern hemisphere (excluding the
equatorial ring) for an `nside` HEALPix map.
""" issouth

####

"""
    pix2ring(nside, p) -> i

Computes the ring index `i` for the given pixel `p`. `nside` is the Nside resolution
factor.
"""
@fastmath function pix2ring(nside::I, p::I) where I<:Integer
    F = float(I)
    p′ = isnorth(nside, p) ? p : (nside2npix(nside)-one(I)) - p
    if isnorthcap(nside, p′)
        i′ = trunc(I, sqrt((p′+one(I))/2 - sqrt(convert(F,(p′+one(I))>>1)))) + one(I)
    else
        i′ = div(p′-nside2npixcap(nside), 4nside) + nside
    end
    i = isnorth(nside, p) ? i′ : (nside2nring(nside)+one(I)) - i′
    return i
end
pix2ring(nside, p) = pix2ring(promote(nside, p)...)

"""
    pix2ringidx(nside, p) -> j

Computes the index `j` within the ring for the given pixel `p`. `nside` is the Nside
resolution factor.
"""
@fastmath function pix2ringidx(nside::I, p::I) where I<:Integer
    F = float(I)
    p′ = isnorth(nside, p) ? p : (nside2npix(nside)-1) - p
    if isnorthcap(nside, p′)
        i′ = trunc(I, sqrt((p′+one(I))/2 - sqrt(convert(F,(p′+one(I))>>1)))) + one(I)
        j′ = p′ + one(I) - nside2npixcap(i′)
        l′ = 4i′
    else
        j′ = rem(p′-nside2npixcap(nside), 4nside) + one(I)
        l′ = 4nside
    end
    j = isnorth(nside, p) ? j′ : (l′ - j′ + one(I))
    return j
end
pix2ringidx(nside, p) = pix2ringidx(promote(nside, p)...)

"""
    pix2z(nside, p) -> z

Computes the cosine of the colatitude `z` for the given pixel `p`. `nside` is the Nside
resolution factor.
"""
@fastmath function pix2z(nside::I, p::I) where I<:Integer
    F = float(I)
    p′ = isnorth(nside, p) ? p : (nside2npix(nside)-1) - p
    if isnorthcap(nside, p′)
        i′ = trunc(I, sqrt((p′+one(I))/2 - sqrt(convert(F,(p′+one(I))>>1)))) + one(I)
        z′ = one(F) - i′^2 / (3nside^2)
    else
        i′ = div(p′-nside2npixcap(nside), 4nside) + nside
        z′ = 4/3 - 2i′ / (3nside)
    end
    z = p < nside2npixequ(nside) ? z′ : -z′
    return z
end
pix2z(nside, p) = pix2z(promote(nside, p)...)

"""
    pix2theta(nside, p) -> θ

Computes the colatitude `θ` for the given pixel `p`. `nside` is the Nside resolution
factor.
"""
pix2theta(nside, p) = @fastmath acos(pix2z(promote(nside, p)...))

"""
    pix2phi(nside, p) -> ϕ

Computes the azimuth `ϕ` for the given pixel `p`. `nside` is the Nside resolution
factor.
"""
@fastmath function pix2phi(nside::I, p::I) where I<:Integer
    F = float(I)
    p′ = isnorth(nside, p) ? p : (nside2npix(nside)-one(I)) - p
    if isnorthcap(nside, p′)
        i′ = trunc(I, sqrt((p′+one(I))/2 - sqrt(convert(F,(p′+one(I))>>1)))) + one(I)
        j′ = p′ + one(I) - nside2npixcap(i′)
        ϕ′ = (π/2)/i′ * (j′ - 0.5)
    else
        i′,j′ = divrem(p′-nside2npixcap(nside), 4nside)
        i′ += nside
        j′ += one(I)
        δ′ = rem(i′ - nside, 2)
        ϕ′ = (π/2/nside) * (j′ - 0.5*(one(I) + δ′))
        if issouth(nside, p) && δ′ ≠ 0
            ϕ′ += π/2/nside
        end
    end
    ϕ = isnorth(nside, p) ? ϕ′ : (2π - ϕ′)
    return ϕ
end
pix2phi(nside, p) = pix2phi(promote(nside, p)...)

"""
    pix2ang(nside, p) -> (θ,ϕ)

Computes the colatitude and azimuth pair `(θ,ϕ)` for the given pixel `p`.
"""
pix2ang(nside, p) = (pix2theta(nside, p), pix2phi(nside, p))

"""
    pix2vec(nside, p) -> SVector{3}(x, y, z)

Computes the unit vector pointing to the pixel center of the given pixel `p`.
"""
@fastmath function pix2vec(nside::I, p::I) where I<:Integer
    z = pix2z(nside, p)
    ϕ = pix2phi(nside, p)
    # If z ≈ ±1, this form should cause less "catastrophic cancellation" than the simpler
    # invocation `one(z) - z*z` (I think...).
    sinθ = sqrt((one(z)-z)*(one(z)+z))
    # TODO: switch to sincos() when julia v0.7 is base version
    x = sinθ * cos(ϕ)
    y = sinθ * sin(ϕ)
    return SVector{3}(x, y, z)
end
pix2vec(nside, p) = pix2vec(promote(nside, p)...)

end # module Healpix

