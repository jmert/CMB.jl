"""
Collection of functions which compute the pixel-pixel covariance of the CMB
sky.

Based on equations given in Tegmark and de Oliveira-Costa (2001) *“How to
measure CMB polarization power spectra without losing information”*
arXiv:astro-ph/0012120v3
"""
module PixelCovariance

export Fweights!,
    PixelCovarianceCache,
    selectpixel!,
    pixelcovariance, pixelcovariance!

using BitFlags
using StaticArrays
using Legendre
using Legendre: unsafe_legendre!

# For pixelcovariance() wrapper function
import ..Sphere: bearing2, cosdistance
import ..Healpix: pix2vec
import ..unchecked_sqrt

import Base: @propagate_inbounds, checkindex, checkbounds_indices, OneTo, Slice

struct FweightsWork{T,N,V<:AbstractArray{T}}
    legwork::Legendre.Work{T,N,V}
    x::V
    y::V
    xy::V
end
function FweightsWork(norm::AbstractLegendreNorm, F, z)
    T = promote_type(eltype(norm), eltype(F), eltype(z))
    legwork = Legendre.Work(norm, F, z)
    x = legwork.z
    y = similar(x)
    xy = similar(x)
    return FweightsWork(legwork, x, y, xy)
end

@inline function coeff_η(::Type{T}, l::Integer) where T
    return T(2l + 1)
end
@inline function coeff_χ(::Type{T}, l::Integer) where T
    lT = convert(T, l)
    fac1 = 4 * T(2lT + one(T))^2 * lT
    fac2 = @evalpoly(lT, -2, -1, 2, 1)
    return unchecked_sqrt(fac1 / fac2)
end
@inline function coeff_γ(::Type{T}, l::Integer) where T
    lT = convert(T, l)
    fac1 = 2 * T(lT + one(T))
    fac2 = @evalpoly(lT, 0, -2, -1, 2, 1)
    return fac1 / fac2
end

function _chkdomain(lmax)
    @noinline _chkdomain_throw_lmax(l) = throw(DomainError(l, "degree lmax must be non-negative"))
    0 ≤ lmax || _chkdomain_throw_lmax(lmax)
    nothing
end

function _chkbounds(F, lmax, x)
    @noinline _chkbounds_throw_dims(M, N) = throw(DimensionMismatch(
            "Output has $M dimensions, expected $(N+2)"))
    @noinline _chkbounds_throw_axes(F, x) = begin
        throw(DimensionMismatch(
            "Output has leading axes $(ntuple(i -> axes(F,i), ndims(x))), expected $(axes(x))"))
    end
    @noinline _chkbounds_throw_lmax() = throw(DimensionMismatch(
            "lmax incompatible with output array axes"))
    @noinline _chkbounds_throw_poldim() = throw(DimensionMismatch(
            "incompatible output array axes; expected last dimension with axes 1:4"))

    M = ndims(F)
    N = ndims(x)
    M == N + 2 || _chkbounds_throw_dims(M, N)
    # Leading dimensions of F are storage for the same dimensions as x
    axesF = axes(F)
    if N > 0
        axF = ntuple(i -> axesF[i], N)
        axx = axes(x)
        checkbounds_indices(Bool, axF, axx) || _chkbounds_throw_axes(F, x)
    end
    # Trailing dimensions of F are storage for range of ell and 4-pol types
    checkindex(Bool, axesF[N+1], OneTo(lmax+1)) || _chkbounds_throw_lmax()
    checkindex(Bool, axesF[N+2], OneTo(4)) || _chkbounds_throw_poldim()
    nothing
end
@noinline function _chkbounds_F(F, lmax)
    (0 ≤ lmax) || throw(DomainError())
    (size(F,1)≥lmax+1 && size(F,2)≥4) || throw(DimensionMismatch())
end

function unsafe_Fweights!(workornorm::Union{AbstractLegendreNorm,FweightsWork}, F, lmax, x)
    if ndims(x) > 1
        M = ndims(F)
        N = ndims(x)
        S = prod(size(x))
        x′ = reshape(x, S)
        F′ = reshape(F, S, ntuple(i->size(F,i+N), M-N)...)
    else
        x′ = x
        F′ = F
    end
    if workornorm isa AbstractLegendreNorm
        work = FweightsWork(LegendreUnitNorm(), F′, x′)
        @inbounds _Fweights_impl!(work, F′, lmax, x′)
    else
        @inbounds _Fweights_impl!(workornorm, F′, lmax, x′)
    end
    return F
end

@propagate_inbounds function _Fweights_impl!(work::FweightsWork, F, lmax, z)
    # Use a local isapprox function instead of Base.isapprox. We get far fewer instructions with
    # this implementation. (Probably related to the keyword-argument penalty?)
    local @inline ≈(x, y) = @fastmath x==y || abs(x-y) < eps(one(y))

    x  = work.x
    y  = work.y
    xy = work.xy
    legwork = work.legwork
    T = eltype(x)
    half = one(T) / 2

    Is = map(Slice, axes(x))
    I = CartesianIndices(Is)
    P = @view F[Is..., :, 1]

    @simd for ii in I
        x[ii] = x′ = convert(T, z[ii])
        y′ = inv(fma(-x′, x′, one(T)))
        y[ii], xy[ii] = y′, x′ * y′
    end

    # Calculations below implicitly assume P has been zeroed in domain where ℓ < m, so
    # zero first two ℓ entries before filling the P^2_ℓ(x) terms
    fill!(@view(F[Is...,1:2,:]), zero(T))

    # Fill with the P^2_ℓ(x) terms initially
    unsafe_legendre!(legwork, P, lmax, 2, x)
    # Calculate the F12 and F22 terms using P^2_ℓ
    for ll in 2:lmax
        η = coeff_η(T, ll)
        γ = coeff_γ(T, ll)
        lT = convert(T, ll)
        lp2 = lT + 2one(T)
        lm1 = lT - one(T)
        lm4 = lT - 4one(T)

        # Closure over pre-computed shared terms
        @inline function F12F22(Plm1, Pl, x, y, xy)
            if abs(x) ≈ one(T)
            # Case where two points are antipodes
                # x == +1, both F12 and F22 are constants (before applying the normalization)
                # x == -1, the F12 and F22 terms flip signs at each ℓ (modulo normalization)
                shalf = iseven(ll) ? half : -half
                F12, F22 = !signbit(x) ?
                        (η *  half, η * -half) :
                        (η * shalf, η * shalf)
            else # abs(x) ≈ one(T)
            # Case where two points are not antipodes
                F12 =  γ * (lp2 * xy * Plm1 - (lm4 * y + half * lT * lm1) * Pl)
                F22 = 2γ * (lp2 *  y * Plm1 - lm1 * xy * Pl)
            end # abs(x) ≈ one(T)
            return (F12, F22)
        end

        @simd for ii in I
            x′, y′, xy′ = x[ii], y[ii], xy[ii]
            Plm1, Pl = P[ii,ll], P[ii,ll+1]
            F[ii,ll+1,3], F[ii,ll+1,4] = F12F22(Plm1, Pl, x′, y′, xy′)
        end
    end
    # Replace with P^0_ℓ(x) terms
    unsafe_legendre!(legwork, P, lmax, 0, x)
    # Then calculate the F10 terms
    for ll in 2:lmax
        lm1 = convert(T, ll) - one(T)
        χ = coeff_χ(T, ll)

        @inline function F10(Plm1, Pl, x, y, xy)
            if abs(x) ≈ one(T)
            # Case where two points are antipodes
                return zero(T)
            else
            # Case where two points are not antipodes
                return χ * (xy * Plm1 - (y + half * lm1) * Pl)
            end
        end

        @simd for ii in I
            x′, y′, xy′ = x[ii], y[ii], xy[ii]
            Plm1, Pl = P[ii,ll], P[ii,ll+1]
            F[ii,ll+1,2] = F10(Plm1, Pl, x′, y′, xy′)
        end
    end

    # Now finally apply the normalization to the P^0_ℓ(x) function to make F00
    for ll in 0:lmax
        η = coeff_η(T, ll)
        @simd for ii in I
            P[ii,ll+1] *= η
        end
    end

    return F
end

function Fweights!(norm::AbstractLegendreNorm, F, lmax::Integer, x)
    _chkdomain(lmax)
    _chkbounds(F, lmax, x)
    return unsafe_Fweights!(norm, F, lmax, x)
end

"""
    @bitflag CovarianceFields

A bitfield for identifying subblocks of the pixel-pixel covariance matrix. There are 9
subblocks, named as the Cartesian product of elements T, Q, and U:

    TT  TQ  TU
    QT  QQ  QU
    UT  UQ  UU
"""
@bitflag CovarianceFields TT TQ TU QT QQ QU UT UQ UU NO_FIELD=0
const TPol = TQ | TU | UT | QT
const Pol  = QQ | UU | QU | UQ

"""
    struct PixelCovarianceCache

Data structure which contains all the information and buffers required to compute the
pixel-pixel covariance terms for a given pixel with respect to all other pixels.
"""
struct PixelCovarianceCache
    nside::Int
    lmax::Int
    pixels::Vector{Int}
    fields::CovarianceFields

    pixind::Ref{Int}
    r::Vector{SVector{3,Float64}}
    z::Vector{Float64}
    cij::Vector{Float64}
    sij::Vector{Float64}
    cji::Vector{Float64}
    sji::Vector{Float64}

    F::Matrix{Float64}

    """
        PixelCovarianceCache(nside, lmax, pixels, fields::CovarianceFields = Pol)
    """
    function PixelCovarianceCache(nside, lmax, pixels,
                                  fields::CovarianceFields = Pol)
        N = length(pixels)

        r = zeros(SVector{3,Float64}, N)
        z = zeros(Float64, N)
        cij = zeros(Float64, N)
        sij = zeros(Float64, N)
        cji = zeros(Float64, N)
        sji = zeros(Float64, N)

        F = zeros(Float64, lmax+1, 4)

        # The pixel coordinates can all be precomputed just once
        r .= pix2vec.(nside, pixels)

        return new(nside, lmax, pixels, bitfields,
                   1, r, z, cij, sij, cji, sji, F)
    end
end

# Improve printing somewhat
Base.show(io::IO, norm::PixelCovarianceCache) = print(io, PixelCovarianceCache)
function Base.show(io::IO, ::MIME"text/plain", C::PixelCovarianceCache)
    pix = C.pixind[]
    println(io, C)
    println(io, "    HEALPix nside: $(C.nside)")
    println(io, "    number of pixels: $(length(C.pixels))")
    println(io, "    maximum ℓ mode: $(C.lmax)")
    println(io, "    selected pixel: ", pix, ", at r = $(C.r[pix])")
    println(io, "    covariance blocks: $(C.fields)")
end

function selectpixel!(cache, pixind)
    cache.pixind[] = pixind
    @inbounds for ii in 1:length(cache.z)
        cache.z[ii] = cosdistance(cache.r[pixind], cache.r[ii])
        c,s = bearing2(cache.r[ii], cache.r[pixind])
        cache.sij[ii] = 2*c*s
        cache.cij[ii] = c*c - s*s
        c,s = bearing2(cache.r[pixind], cache.r[ii])
        cache.sji[ii] = 2*c*s
        cache.cji[ii] = c*c - s*s
    end

    return cache
end

function pixelcovariance(nside, pixels, pixind)
    lmax = size(spec,1)
    cache = PixelCovarianceCache(nside, lmax, pixels)
    selectpixel!(cache, pixind)
    return cache
end

function pixelcovariance!(cache::PixelCovarianceCache, cov::AbstractMatrix, Cl::AbstractMatrix)
    T = eltype(cache.F)
    fourpi = 4 * convert(T, π)
    # N.B. In general, the spectrum causes Cl*F to decrease rapidly as ℓ → ∞; reverse
    #      the order of summation to accumulate from smallest to largest magnitude values.
    # TODO: Look into using and compare this assumption against Julia's built-in sum() which
    #       uses a divide-and-conquered summation to increase numerical precision in general
    #       cases.
    R = reverse(axes(cache.F, 1))

    @inbounds for (i,z) in enumerate(cache.z)
        Fweights!(LegendreUnitNorm(), cache.F, cache.lmax, z)

        # TT
        if cache.fields & TT != NO_FIELD
            tt = zero(T)
            @simd for ll in R
                ClTT = Cl[ll,1]
                F00 = cache.F[ll, 1]
                tt = muladd(ClTT, F00, tt) # tt + ClTT*F00
            end
            tt /= fourpi
            cov[i,1] = tt # TT
        end

        # TQ and TU
        if cache.fields & TPol != NO_FIELD
            tq = zero(T)
            tu = zero(T)
            @simd for ll in R
                ClTE = Cl[ll,4]
                ClTB = Cl[ll,5]
                F10 = F[ll,2]
                tq = muladd(-ClTE, F10, tq) # tq - ClTE*F10
                tu = muladd(-ClTB, F10, tu) # tu - ClTB*F10
            end
            tq /= fourpi
            tu /= fourpi
            cov[i,2] =  tq*cache.cij[i] + tu*cache.sij[i] # QT
            cov[i,3] = -tq*cache.sij[i] + tu*cache.cij[i] # QU
            cov[i,4] =  tq*cache.cji[i] + tu*cache.sji[i] # TQ
            cov[i,7] = -tq*cache.sji[i] + tu*cache.cji[i] # TU
        end

        # QQ, QU, and UU
        if cache.fields & Pol != NO_FIELD
            qq = zero(T)
            qu = zero(T)
            uu = zero(T)
            @simd for ll in R
                ClEE = Cl[ll,2]
                ClBB = Cl[ll,3]
                ClEB = Cl[ll,6]
                F12 = cache.F[ll,3]
                F22 = cache.F[ll,4]

                qq = muladd(ClEE, F12, muladd(-ClBB, F22, qq)) # qq + ClEE*F12 - ClBB*F22
                uu = muladd(ClBB, F12, muladd(-ClEE, F22, uu)) # uu + ClBB*F12 - ClEE*F22
                qu = muladd(ClEB, F12 + F22, qu)               # qu + ClEB*(F12 + F22)
            end
            qq /= fourpi
            uu /= fourpi
            qu /= fourpi
            cij = cache.cij[i]
            cji = cache.cji[i]
            sij = cache.sij[i]
            sji = cache.sji[i]
            cov[i,5] =  qq*cij*cji + qu*(cij*sji+sij*cji) + uu*sij*sji # QQ
            cov[i,6] = -qq*sij*cji + qu*(cij*cji-sij*sji) + uu*cij*sji # UQ
            cov[i,8] = -qq*cij*sji + qu*(cij*cji-sij*sji) + uu*sij*cji # QU
            cov[i,9] =  qq*sij*sji - qu*(cij*sji+sij*cji) + uu*cij*cji # UU
        end
    end

    return cov
end

end # module PixelCovariance

