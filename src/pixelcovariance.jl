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
    copyspectra!, makespectra!, applybeam!, selectpixel!,
    pixelcovariance, pixelcovariance!

using StaticArrays

# For computing the Legendre terms
import Legendre: AbstractLegendreNorm, legendre!
# also reuse some utilities developed in Legendre
import Legendre: _similar

# For pixelcovariance() wrapper function
import ..Sphere: bearing2, cosdistance
import ..Healpix: pix2vec

import Base.@boundscheck, Base.@propagate_inbounds

@inline function Fweights_η(::Type{T}, l::Integer) where {T}
    return convert(T, 2l+1) / 4convert(T, π)
end

@inline function Fweights_αβ(::Type{T}, l::Integer) where {T}
    η = Fweights_η(T, l)
    lT = convert(T, l)
    invql = inv( (lT-one(T)) * lT * (lT+one(T)) * (lT+2one(T)) )

    α = 2η * lT * sqrt(invql)
    β = 2η * invql
    return (α, β)
end

@noinline function _chkbounds_F(F, lmax)
    (0 ≤ lmax) || throw(DomainError())
    (size(F,1)≥lmax+1 && size(F,2)≥4) || throw(DimensionMismatch())
end

function Fweights!(Λ::AbstractLegendreNorm, F, lmax::Integer, z)
    # Use a local isapprox function instead of Base.isapprox. We get far fewer instructions with
    # this implementation. (Probably related to the keyword-argument penalty?)
    local @inline ≈(x, y) = @fastmath x==y || abs(x-y) < eps(one(y))

    #@boundscheck _chkbounds_F(F, lmax)
    T = promote_type(eltype(Λ), eltype(F), eltype(z))
    half = one(T) / 2

    x  = _similar(z)
    y  = _similar(z)
    xy = _similar(z)
    @inbounds @simd for I in eachindex(x)
        x[I] = x′ = convert(T, z[I])
        y′ = -inv(fma(x′, x′, -one(T)))
        y[I], xy[I] = y′, x′ * y′

        # Clear out the ℓ == 0 and ℓ == 1 terms since they are undefined (denominators
        # of α and β are singular).
        F[I,1,1] = F[I,2,1] = zero(T)
        F[I,1,2] = F[I,2,2] = zero(T)
        F[I,1,3] = F[I,2,3] = zero(T)
        F[I,1,4] = F[I,2,4] = zero(T)
    end

    Iv = ntuple(_ -> :, ndims(x) + 1)
    P = view(F, Iv..., 1)

    # Fill with the P^2_ℓ(x) terms initially
    legendre!(Λ, P, lmax, 2, x)
    # Calculate the F12 and F22 terms using P^2_ℓ
    for ll in 2:lmax
        η = Fweights_η(T, ll)
        _, β = Fweights_αβ(T, ll)
        lT = convert(T, ll)
        lp2 = lT + convert(T, 2)
        lm1 = lT - one(T)
        lm4 = lT - convert(T, 4)

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
                F12 =  β * (lp2 * xy * Plm1 - (lm4 * y + half * lT * lm1) * Pl)
                F22 = 2β * (lp2 *  y * Plm1 - lm1 * xy * Pl)
            end # abs(x) ≈ one(T)
            return (F12, F22)
        end

        @inbounds @simd for I in eachindex(x)
            x′, y′, xy′ = x[I], y[I], xy[I]
            Plm1, Pl = P[I,ll], P[I,ll+1]
            F[I,ll+1,3], F[I,ll+1,4] = F12F22(Plm1, Pl, x′, y′, xy′)
        end
    end
    # Replace with P^0_ℓ(x) terms
    legendre!(Λ, P, lmax, 0, x)
    # Then calculate the F10 terms
    for ll in 2:lmax
        lm1 = convert(T, ll) - one(T)
        α, _ = Fweights_αβ(T, ll)

        @inline function F10(Plm1, Pl, x, y, xy)
            if abs(x) ≈ one(T)
            # Case where two points are antipodes
                return zero(T)
            else
            # Case where two points are not antipodes
                return α * (xy * Plm1 - (y + half * lm1) * Pl)
            end
        end

        @inbounds @simd for I in eachindex(x)
            x′, y′, xy′ = x[I], y[I], xy[I]
            Plm1, Pl = P[I,ll], P[I,ll+1]
            F[I,ll+1,2] = F10(Plm1, Pl, x′, y′, xy′)
        end
    end

    # Now finally apply the normalization to the P^0_ℓ(x) function to make F00
    for ll in 0:lmax
        η = Fweights_η(T, ll)
        @inbounds @simd for I in eachindex(x)
            P[I,ll+1] *= η
        end
    end

    return F
end

"""
    const FIELDMAP

A symbol array which states the canonical ordering of the block matrices within a
pixel-pixel covariance matrix.
"""
const FIELDMAP = [:TT :TQ :TU;
                  :QT :QQ :QU;
                  :UT :UQ :UU]

"""
    const SPECTRAMAP

A symbol array which states the canonical ordering of spectra.
"""
const SPECTRAMAP = [:TT :EE :BB :TE :TB :EB]

"""
    struct PixelCovarianceCache

Data structure which contains all the information and buffers required to compute the
pixel-pixel covariance terms for a given pixel with respect to all other pixels.
"""
struct PixelCovarianceCache
    nside::Int
    lmax::Int
    pixels::Vector{Int}
    fields::BitArray{2}

    spectra::Matrix{Float64}

    pixind::Ref{Int}
    r::Vector{SVector{3,Float64}}
    z::Vector{Float64}
    cij::Vector{Float64}
    sij::Vector{Float64}
    cji::Vector{Float64}
    sji::Vector{Float64}

    F::Matrix{Float64}

    """
        PixelCovarianceCache(nside, lmax, pixels, fields=Symbol[:QQ,:QU,:UQ,:UU])
    """
    function PixelCovarianceCache(nside, lmax, pixels,
                                  fields=Symbol[:QQ,:QU,:UQ,:UU])
        N = length(pixels)

        # Map the symbol description of which fields to generate into a bit array
        bitfields = [(ff in fields) for ff in FIELDMAP]

        spectra = zeros(Float64, lmax+1, 6)

        r = zeros(SVector{3,Float64}, N)
        z = zeros(Float64, N)
        cij = zeros(Float64, N)
        sij = zeros(Float64, N)
        cji = zeros(Float64, N)
        sji = zeros(Float64, N)

        F = zeros(Float64, lmax+1, 4)

        # The pixel coordinates can all be precomputed just once
        r .= pix2vec.(nside, pixels)

        return new(nside, lmax, pixels, bitfields, spectra,
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
    println(io, "    covariance blocks: $(reshape(FIELDMAP[C.fields],:))")
end

function copyspectra!(cache, spectra)
    cache.spectra .= spectra
    return cache
end

function makespectra!(cache, f, fields=[:TT,:EE])
    lmax = cache.lmax

    goodinds = map(f->in(f, SPECTRAMAP), fields)
    all(goodinds) || error("Bad field specification(s): $(fields[.!(goodinds)]...)")

    fieldinds = map(f -> findall(f .== SPECTRAMAP)[1], fields)
    for ll in 0:lmax
        Cl = f(ll)
        Cl = isfinite(Cl) ? Cl : zero(eltype(cache.spectra))
        for ff in fieldinds
            cache.spectra[ll+1,ff] = Cl
        end
    end
    return cache
end

function applybeam!(cache, fwhm)
    lmax = cache.lmax
    sigma = deg2rad(fwhm/60) / (2sqrt(2log(2)))

    fac = -0.5 * sigma*sigma
    for ll in 0:lmax
        Bl = exp(ll*(ll+1) * fac)
        cache.spectra[ll+1,:] .*= Bl
    end
    return cache
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

function pixelcovariance(nside, pixels, pixind, spec)
    lmax = size(spec,1)
    cache = PixelCovarianceCache(nside, lmax, pixels)
    updatespectra!(cache, spec)
    selectpixel!(cache, pixind)
    return cache
end

function pixelcovariance!(cache::PixelCovarianceCache, C::AbstractMatrix)
    T = eltype(cache.F)
    R = size(cache.F,1):-1:1

    @inbounds for (i,z) in enumerate(cache.z)
        Fweights!(LegendreUnitNorm(), cache.F, cache.lmax, z)

        # TT
        if cache.fields[1,1]
            tt = zero(T)
            for ll in R
                ClTT = cache.spectra[ll,1]
                tt += ClTT*cache.F[ll,1]
            end
            C[i,1] = tt # TT
        end

        # TQ and TU
        if cache.fields[1,2] | cache.fields[1,3]
            tq = zero(T)
            tu = zero(T)
            for ll in R
                ClTE = cache.spectra[ll,4]
                ClTB = cache.spectra[ll,5]
                tq -= ClTE*cache.F[ll,2]
                tu -= ClTB*cache.F[ll,2]
            end
            C[i,2] =  tq*cache.cij[i] + tu*cache.sij[i] # QT
            C[i,3] = -tq*cache.sij[i] + tu*cache.cij[i] # QU
            C[i,4] =  tq*cache.cji[i] + tu*cache.sji[i] # TQ
            C[i,7] = -tq*cache.sji[i] + tu*cache.cji[i] # TU
        end

        # QQ, QU, and UU
        if cache.fields[2,2] | cache.fields[2,3] | cache.fields[3,3]
            qq = zero(T)
            qu = zero(T)
            uu = zero(T)
            for ll in R
                ClEE = cache.spectra[ll,2]
                ClBB = cache.spectra[ll,3]
                ClEB = cache.spectra[ll,6]
                F12 = cache.F[ll,3]
                F22 = cache.F[ll,4]

                qq += ClEE*F12 - ClBB*F22
                uu += ClBB*F12 - ClEE*F22
                qu += ClEB*(F12 + F22)
            end
            cij = cache.cij[i]
            cji = cache.cji[i]
            sij = cache.sij[i]
            sji = cache.sji[i]
            C[i,5] =  qq*cij*cji + qu*(cij*sji+sij*cji) + uu*sij*sji # QQ
            C[i,6] = -qq*sij*cji + qu*(cij*cji-sij*sji) + uu*cij*sji # UQ
            C[i,8] = -qq*cij*sji + qu*(cij*cji-sij*sji) + uu*sij*cji # QU
            C[i,9] =  qq*sij*sji - qu*(cij*sji+sij*cji) + uu*cij*cji # UU
        end
    end

    return cache
end

end # module PixelCovariance

