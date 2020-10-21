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
    nside2npix, nside2nring, nside2npixcap, nside2npixequ, nside2pixarea,
    isnorth, issouth,
    isnorthcap, issouthcap, iscap,
    isnorthequbelt, issouthequbelt, isequbelt,
    pix2ring, pix2ringidx, pix2z, pix2theta, pix2phi, pix2ang, pix2vec,
    ang2pix, vec2pix,
    UNSEEN, ishealpixok, checkhealpix, InvalidNside, InvalidPixel

using Base: @propagate_inbounds
using StaticArrays
import ..sincospi, ..unchecked_sqrt

"""
    const UNSEEN = -1.6375e+30

Special value recognized by the libhealpix/healpy routines as an unobserved/masked
pixel.
"""
const UNSEEN = -1.6375e+30

"""
    const MAX_NSIDE = 2^29

Maximum valid ``N_\\mathrm{side}`` parameter value.
"""
const MAX_NSIDE = 2^29

"""
    InvalidNside(nside)

An invalid `nside` value was provided.
"""
struct InvalidNside <: Exception
    nside::Int
end
Base.showerror(io::IO, e::InvalidNside) =
        print(io, "$(e.nside) is not a valid Nside parameter (must be power of 2)")

"""
    InvalidPixel(nside, pix)

An invalid pixel index `pix` was provided for the given `nside`.
"""
struct InvalidPixel <: Exception
    nside::Int
    pix::Int
end
Base.showerror(io::IO, e::InvalidPixel) =
        print(io, "$(e.pix) is not a valid pixel index for Nside = $(e.nside) " *
              "(must be from 0 to $(nside2npix(e.nside)-1))")

isvalidnside(nside) = (1 ≤ nside ≤ MAX_NSIDE) && ispow2(nside)
isvalidpixel(nside, pix) = 0 ≤ pix < nside2npix(nside)

"""
    ishealpixok(nside)

Returns `true` if `nside` is a power of two in the range `1` to
`2^$(Int(log2(MAX_NSIDE)))`, otherwise `false`.
"""
ishealpixok(nside) = isvalidnside(nside)

"""
    isheapixok(nside, pix)

Returns `true` if `nside` is valid and `pix` is in the range `0` to `nside2npix(nside) - 1`,
otherwise `false`.
"""
ishealpixok(nside, pix) = isvalidnside(nside) && isvalidpixel(nside, pix)

"""
    checkhealpix(nside)

Throws an [`InvalidNside`](@ref) exception if `nside` is not a valid value.
"""
@noinline function checkhealpix(nside)
    isvalidnside(nside) || throw(InvalidNside(nside))
    return nothing
end

"""
    checkhealpix(nside, pix)

Throws an [`InvalidNside`](@ref) exception if `nside` is not a valid value or an
[`InvalidPixel`](@ref) exception if `pix` is out of range for the given
``N_\\mathrm{side}``.
"""
@noinline function checkhealpix(nside, pix)
    isvalidnside(nside) || throw(InvalidNside(nside))
    isvalidpixel(nside, pix) || throw(InvalidPixel(nside, pix))
    return nothing
end

# As relatively low-level functions, we don't validate nside or pixel values

npix2nside(npix::T) where {T<:Integer}   = trunc(T, sqrt(npix / T(12)))
nring2nside(nring::T) where {T<:Integer} = div(nring + one(nring), T(4))

nside2npix(nside::T) where {T<:Integer}    = T(12)*nside*nside
nside2nring(nside::T) where {T<:Integer}   =  T(4)*nside - one(nside)
nside2npixcap(nside::T) where {T<:Integer} =  T(2)*nside*(nside - one(nside))
nside2npixequ(nside::T) where {T<:Integer} =  T(2)*nside*(T(3)*nside + one(nside))
nside2pixarea(nside::T) where {T<:Integer} = 4convert(float(T), π) / nside2npix(nside)

isnorth(nside, p) = p < nside2npixequ(nside)
issouth(nside, p) = p ≥ nside2npixequ(nside)

isnorthcap(nside, p) = p < nside2npixcap(nside)
issouthcap(nside, p) = p ≥ nside2npix(nside) - nside2npixcap(nside)
iscap(nside, p)      = isnorthcap(nside, p) | issouthcap(nside, p)

isnorthequbelt(nside, p) = isnorth(nside, p) & ~isnorthcap(nside, p)
issouthequbelt(nside, p) = issouth(nside, p) & ~issouthcap(nside, p)
isequbelt(nside, p)      = ~iscap(nside, p)

# The function definitions are short enough that I wanted to keep them visually together
# without having the documentation break them up. Add documentation now

"""
    Nisde = npix2nside(npix)

Returns the equivalent `Nside` corresponding to the number of pixels `npix`.
""" npix2nside(npix)

"""
    Nside = nring2nside(nring)

Returns the equivalent `Nside` corresponding to the number of iso-latitude rings `nring`.
""" nring2nside(npix)

"""
    Npix = nside2npix(nside)

Returns the total number of pixels `Npix` in an `nside` HEALPix map. Note that `HEALPix`
pixel indexing is 0-based, so valid pixel values are in the range `0` to `Npix - 1`.
""" nside2npix

"""
    Nring = nside2nring(nside)

Returns the number of iso-latitude rings `Nring` in the `nside` HEALPix map.
""" nside2nring

"""
    Npix = nside2npixcap(nside)

Returns the number of pixels `Npix` in the polar caps for the given `nside` HEALPix map.
""" nside2npixcap

"""
    Npix = nside2npixequ(nside)

Returns the number of pixels `Npix` in the northern hemisphere, including the equatorial
ring, for the given `nside` HEALPix map.
""" nside2npixequ

"""
    σ = nside2pixarea(nside)

Returns the surface area `σ` (in steradians) of each pixel in the given `nside` HEALPix
map.
""" nside2pixarea

"""
    isnorthcap(nside, pix)

Test for whether the given pixel `pix` is in the northern polar cap for an `nside`
HEALPix map.
""" isnorthcap

"""
    issouthcap(nside, pix)

Test for whether the given pixel `pix` is in the southern polar cap for an `nside`
HEALPix map.
""" issouthcap

"""
    iscap(nside, pix)

Test for whether the given pixel `pix` is in either polar cap for an `nside` HEALPix
map.
""" iscap

"""
    isnorthequbelt(nside, pix)

Test for whether the given pixel `pix` is in the northern equatorial belt (including the
equatorial ring) for an `nside` HEALPix map.
""" isnorthequbelt

"""
    issouthequbelt(nside, pix)

Test for whether the given pixel `pix` is in the southern equatorial belt (excluding the
equatorial ring) for an `nside` HEALPix map.
""" issouthequbelt

"""
    isequbelt(nside, pix)

Test for whether the given pixel `pix` is in the equatorial belt for an `nside` HEALPix
map.
""" isequbelt

"""
    isnorth(nside, pix)

Test for whether the given pixel `pix` is in the northern hemisphere (including the
equatorial ring) for an `nside` HEALPix map.
""" isnorth

"""
    issouth(nside, pix)

Test for whether the given pixel `pix` is in the southern hemisphere (excluding the
equatorial ring) for an `nside` HEALPix map.
""" issouth

####

# internal method - pixel ring index for polar cap pixel
@inline function _capring(p::I) where {I <: Integer}
    F = float(I)
    pp1 = p + one(p)
    ring0 = unchecked_sqrt(muladd(F(pp1), F(0.5), -unchecked_sqrt(pp1 >> 1)))
    return unsafe_trunc(I, ring0) + one(I)
end

"""
    i = pix2ring(nside, p)

Computes the ring index `i` for the given pixel `p`, where `nside` is the Nside resolution
factor.
"""
function pix2ring(nside::I, p::I) where I<:Integer
    F = float(I)
    p′ = isnorth(nside, p) ? p : (nside2npix(nside)-one(I)) - p
    if isnorthcap(nside, p′)
        i′ = _capring(p′)
    else
        i′ = div(p′-nside2npixcap(nside), 4nside) + nside
    end
    i = isnorth(nside, p) ? i′ : (nside2nring(nside)+one(I)) - i′
    return i
end
pix2ring(nside, p) = pix2ring(promote(nside, p)...)

"""
    j = pix2ringidx(nside, p)

Computes the index `j` within the ring for the given pixel `p`, where `nside` is the Nside
resolution factor.
"""
function pix2ringidx(nside::I, p::I) where I<:Integer
    F = float(I)
    p′ = isnorth(nside, p) ? p : (nside2npix(nside)-1) - p
    if isnorthcap(nside, p′)
        i′ = _capring(p′)
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
    z = pix2z(nside, p)

Computes the cosine of the colatitude `z` for the given pixel `p`, where `nside` is the
Nside resolution factor.
"""
function pix2z(nside::I, p::I) where I<:Integer
    checkhealpix(nside, p)
    return unsafe_pix2z(nside, p)
end
pix2z(nside, p) = pix2z(promote(nside, p)...)

"""
    (θ,ϕ) = unsafe_pix2z(nside, p)

Like [`pix2z`](@ref) but does not call [`checkhealpix`](@ref) to check `nside` and pixel
index validity.
"""
function unsafe_pix2z(nside::I, p::I) where I<:Integer
    z, iscap = _pix2z_capcompl(nside, p)
    z = iscap ? copysign(one(z), z) - z : z
    return z
end

# Internal helper that helps preserve precision for further internal use.
# The return value is location dependent --- if the pixel `p` is within the equatorial
# belt, then the return value is `(z, iscap == false)`, while for the polar caps the
# return value instead gives the "complement" `(±(1 - z), iscap == true)`.
# The advantage is that as `|z| → 1` approaches the poles, there is loss of precision
# while for the complement `|1 - z| → 0` the density of floating point values is
# increasing.
#
# The later definition of `pix2vec` can make use of this to eek out a bit more precision
# because for z = ±(1 - 𝓏),
#   sin θ == sin[acos(z)] == √(1 - z²)
#                         == √(1 - (1 - 𝓏)²)
#                         == √(2𝓏 - 𝓏²)
#                         == sqrt(fma(-𝓏, 𝓏, 2𝓏))
@inline function _pix2z_capcompl(nside::I, p::I) where I<:Integer
    F = float(I)
    p′ = isnorth(nside, p) ? p : (nside2npix(nside)-1) - p
    if isnorthcap(nside, p′)
        i′ = _capring(p′)
        z′ = i′^2 / (F(3) * nside^2)
        iscap = true
    else
        i′ = div(p′-nside2npixcap(nside), 4nside) + nside
        z′ = (4nside - 2i′) / (F(3) * nside)
        iscap = false
    end
    z = p < nside2npixequ(nside) ? z′ : -z′
    return z, iscap
end

"""
    θ = pix2theta(nside, p)

Computes the colatitude `θ` for the given pixel `p`, where `nside` is the Nside resolution
factor.
"""
pix2theta(nside, p) = acos(pix2z(promote(nside, p)...))

"""
    (θ,ϕ) = unsafe_pix2theta(nside, p)

Like [`pix2theta`](@ref) but does not call [`checkhealpix`](@ref) to check `nside` and pixel
index validity.
"""
unsafe_pix2theta(nside, p) = acos(unsafe_pix2z(promote(nside, p)...))

"""
    ϕ = pix2phi(nside, p)

Computes the azimuth `ϕ` for the given pixel `p`, where `nside` is the Nside resolution
factor.
"""
function pix2phi(nside::I, p::I) where I<:Integer
    checkhealpix(nside, p)
    return unsafe_pix2phi(nside, p)
end
pix2phi(nside, p) = pix2phi(promote(nside, p)...)

"""
    (θ,ϕ) = unsafe_pix2phi(nside, p)

Like [`pix2phi`](@ref) but does not call [`checkhealpix`](@ref) to check `nside` and pixel
index validity.
"""
unsafe_pix2phi(nside::I, p::I) where {I<:Integer} = convert(float(I), π) * _pix2phibypi(nside, p)

# Internal helper that helps preserve precision for further internal use.
# Since π is trancendental, truncating it to any given precision necessarily loses
# information. We can eliminate the rounding in some cases by instead returning the scaled
# angle ``ϕ/π`` which is naturally calculated by the geometric construction of the grid.
#
# This is used later in `pix2vec` to eek out a little bit more precision for all azimuthal
# angles and is a good complement to the earlier optimization of `_pix2z_capcompl`.
@inline function _pix2phibypi(nside::I, p::I) where I<:Integer
    F = float(I)
    p′ = isnorth(nside, p) ? p : (nside2npix(nside)-one(I)) - p
    if isnorthcap(nside, p′)
        i′ = _capring(p′)
        j′ = p′ + one(I) - nside2npixcap(i′)
        ϕ′ = (j′ - F(0.5)) / 2i′
    else
        i′,j′ = divrem(p′-nside2npixcap(nside), 4nside)
        i′ += nside
        j′ += one(I)
        δ′ = rem(i′ - nside, 2)
        ϕ′ = (j′ - F(0.5)*(one(I) + δ′)) / 2nside
        if issouth(nside, p) && δ′ ≠ 0
            ϕ′ += one(F) / 2nside
        end
    end
    ϕ = isnorth(nside, p) ? ϕ′ : (F(2) - ϕ′)
    return ϕ
end

"""
    (θ,ϕ) = pix2ang(nside, p)

Computes the colatitude and azimuth pair `(θ,ϕ)` for the given pixel `p`, where
`nside` is the Nside resolution factor.
"""
pix2ang(nside, p) = (pix2theta(nside, p), pix2phi(nside, p))

"""
    (θ,ϕ) = unsafe_pix2ang(nside, p)

Like [`pix2ang`](@ref) but does not call [`checkhealpix`](@ref) to check `nside` and pixel
index validity.
"""
unsafe_pix2ang(nside, p) = (unsafe_pix2theta(nside, p), unsafe_pix2phi(nside, p))

"""
    r = pix2vec(nside, p)

Computes the unit vector `r` pointing to the pixel center of the given pixel `p`, where
`nside` is the Nside resolution factor.
"""
function pix2vec(nside::I, p::I) where I<:Integer
    checkhealpix(nside, p)
    return unsafe_pix2vec(nside, p)
end
pix2vec(nside, p) = pix2vec(promote(nside, p)...)

"""
    (θ,ϕ) = unsafe_pix2vec(nside, p)

Like [`pix2vec`](@ref) but does not call [`checkhealpix`](@ref) to check `nside` and pixel
index validity.
"""
function unsafe_pix2vec(nside::I, p::I) where I<:Integer
    z, iscap = _pix2z_capcompl(nside, p)
    φbypi = _pix2phibypi(nside, p)
    if iscap
        sinθ = unchecked_sqrt(fma(-z, z, 2abs(z)))
        z = copysign(one(z), z) - z
    else
        sinθ = unchecked_sqrt(fma(-z, z, one(z)))
    end
    y, x = sinθ .* sincospi(φbypi)
    return SVector{3}(x, y, z)
end

"""
    p = ang2pix(nside, θ, ϕ)

Computes the HEALPix pixel index `p` which contains the point ``(θ,ϕ)`` given by the
colatitude `θ` and azimuth `ϕ`, where `nside` is the Nside resolution factor.
"""
function ang2pix(nside, θ, ϕ)
    checkhealpix(nside)
    zero(θ) ≤ θ ≤ oftype(θ, π) || throw(DomainError("θ must be in [0,π], but got $θ"))
    ϕ = mod2pi(ϕ)
    return unsafe_ang2pix(nside, θ, ϕ)
end

"""
    p = vec2pix(nside, r)

Computes the HEALPix pixel index `p` which contains the point at the end of the unit
vector `r`, where `nside` is the Nside resolution factor.
"""
function vec2pix(nside, r)
    checkhealpix(nside)
    length(r) == 3 || throw(DimensionMismatch("r must be a 3-vector"))
    return @inbounds unsafe_vec2pix(nside, r)
end

"""
    p = unsafe_ang2pix(nside, θ, ϕ)

Like [`ang2pix`](@ref) but neither calls [`checkhealpix`](@ref) to check the validity of
`nside` nor checks the domain of the spherical coordinates `θ` and `ϕ`.
"""
function unsafe_ang2pix(nside, θ, ϕ)
    z = cos(θ)
    return unsafe_zphi2pix(nside, z, ϕ)
end

"""
    p = unsafe_vec2pix(nside, r)

Like [`vec2pix`](@ref) but does not check the validity of the `nside` or length of `r`.
"""
@propagate_inbounds function unsafe_vec2pix(nside, r)
    z = r[3]
    ϕ = atan(r[2], r[1])
    ϕ += ifelse(ϕ < zero(ϕ), 2oftype(ϕ, π), zero(ϕ))
    return unsafe_zphi2pix(nside, z, ϕ)
end

"""
    p = unsafe_zphi2pix(nside, z, ϕ)

Like [`unsafe_ang2pix`](@ref) but uses the value ``z = \\cos(θ)`` instead.
"""
function unsafe_zphi2pix(nside, z, ϕ)
    z′ = abs(z)
    α = 2ϕ / π      # scaled distance around ring in [0,4)
                    # later advantage is that mod(ϕ, π/2) becomes modf(α)

    if z′ > 2/3
        αt,_ = modf(α)  # fraction across first quadrant

        σ = nside * sqrt(3 * (1 - z′))
        kp = unsafe_trunc(Int, σ * αt)          # NW pixel boundary
        km = unsafe_trunc(Int, σ * (1 - αt))    # SW pixel boundary

        i = kp + km + 1     # intersection of (kp, km+1) or (kp+1, km) lines
        j = unsafe_trunc(Int, i * α) + 1

        if z > 0
            # counting from north pole
            p = nside2npixcap(i) + j - 1
        else
            # counting from south pole
            p = nside2npix(nside) - nside2npixcap(i+1) + j - 1
        end
    else
        tmp = 3z / 4
        kp = unsafe_trunc(Int, nside * (1/2 - tmp + α)) # NW pixel boundary
        km = unsafe_trunc(Int, nside * (1/2 + tmp + α)) # SW pixel boundary
        i′ = nside + kp - km  # ring offset w.r.t. northernmost equatorial belt ring

        s = mod(i′, 2) + 1
        j = (km + kp + s - nside) >> 1
        # rings with first pixel center at ϕ == 0 have region where ϕ == 2π - δ "wrap
        # around" back to first pixel
        if j == 4nside
            j = 0
        end

        p = nside2npixcap(nside) + 4nside * i′ + j
    end
    return p
end

end # module Healpix
