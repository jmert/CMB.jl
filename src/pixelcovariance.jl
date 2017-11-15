"""
Collection of functions which compute the pixel-pixel covariance of the CMB
sky.

Based on equations given in Tegmark and de Oliveira-Costa (2001) *“How to
measure CMB polarization power spectra without losing information”*
arXiv:astro-ph/0012120v3
"""
module PixelCovariance

export PixelCovarianceCoeff, PixelCovarianceF!,
    PixelCovarianceCache,
    copyspectra!, makespectra!, applybeam!, selectpixel!,
    pixelcovariance, pixelcovariance!

using StaticArrays

# For computing the Legendre terms
import ..Legendre: LegendreUnitCoeff, legendre!

# For pixelcovariance() wrapper function
import ..Sphere: bearing2, cosdistance
import ..Healpix: pix2vec

import Base.@boundscheck, Base.@propagate_inbounds

"""
    struct PixelCovarianceCoeff{T<:Real}

Precomputed recursion relation coefficients for computing the pixel-pixel covariance.

# Example
```jldoctest
julia> PixelCovarianceCoeff{Float64}(2)
CMB.PixelCovariance.PixelCovarianceCoeff{Float64} for lmax = 2 with coefficients:
    λ: CMB.Legendre.LegendreNormCoeff{CMB.Legendre.LegendreUnitNorm,Float64}
    η: [0.0795775, 0.238732, 0.397887]
    α: [0.0, 0.0, 0.324874]
    β: [0.0, 0.0, 0.0331573]
```
"""
struct PixelCovarianceCoeff{T<:Real}
    λ::LegendreUnitCoeff{T}
    η::Vector{T}
    α::Vector{T}
    β::Vector{T}

    function PixelCovarianceCoeff{T}(lmax::Integer) where T
        lmax = max(lmax, 2)

        λ = LegendreUnitCoeff{T}(lmax)
        η = Vector{T}(lmax+1)
        α = Vector{T}(lmax+1)
        β = Vector{T}(lmax+1)

        η[1] = 1/(4π);  η[2] = 3/(4π)
        α[1] = zero(T); α[2] = zero(T);
        β[1] = zero(T); β[2] = zero(T);
        @inbounds for ll in 2:lmax
            lT = convert(T, ll)
            norm = convert(T, 2ll + 1) / (4π)
            invql = inv( (lT-one(T)) * lT * (lT+one(T)) * (lT+convert(T,2)) )

            η[ll+1] = norm
            α[ll+1] = 2norm * lT * sqrt(invql)
            β[ll+1] = 2norm * invql
        end

        return new(λ, η, α, β)
    end
end

# Improve printing somewhat
Base.show(io::IO, norm::PixelCovarianceCoeff{T}) where {T} =
    print(io, PixelCovarianceCoeff, "{$T}")
function Base.show(io::IO, ::MIME"text/plain", C::PixelCovarianceCoeff)
    println(io, C, " for lmax = $(length(C.η)-1) with coefficients:")
    println(io, "    λ: ", C.λ)
    println(io, "    η: ", C.η)
    println(io, "    α: ", C.α)
    println(io, "    β: ", C.β)
end

@noinline function _chkbounds_F(C, F, lmax)
    (0 ≤ lmax) || throw(DomainError())
    (lmax ≤ length(C.η)) || throw(BoundsError())
    (size(F,1)≥lmax+1 && size(F,2)≥4) || throw(DimensionMismatch())
end

function PixelCovarianceF!(C::PixelCovarianceCoeff{T}, F::AbstractMatrix{T},
        lmax::Integer, x::T) where {T<:Real}
    @boundscheck _chkbounds_F(C, F, lmax)

    # Use a local isapprox function instead of Base.isapprox. We get far fewer instructions with
    # this implementation. (Probably related to the keyword-argument penalty?)
    local ≈(x::T, y::T) where {T} = @fastmath x==y || abs(x-y) < eps(one(T))

    half = inv(convert(T, 2))

    @inbounds begin
        P = view(F, :, 1)

        # Clear out the ℓ == 0 and ℓ == 1 terms since they aren't defined, and then
        # compute the rest of the terms
        F[1,2] = zero(T)
        F[2,2] = zero(T)
        F[1,3] = zero(T)
        F[2,3] = zero(T)
        F[1,4] = zero(T)
        F[2,4] = zero(T)

        if abs(x) ≈ one(T)
        # Case where two points are antipodes

            # F10 always zero in this case
            for ll=2:lmax
                F[ll+1,2] = zero(T)
            end

            # x == +1, both F12 and F22 are constants (before applying the normalization)
            if sign(x) == one(T)
                for ll=2:lmax
                    F[ll+1,3] =  C.η[ll+1]*half
                    F[ll+1,4] = -C.η[ll+1]*half
                end

            # x == -1, the F12 and F22 terms flip signs at each ℓ (modulo normalization)
            else
                val = half # 1/2 * (-1)^ℓ for ℓ == 2
                for ll=2:lmax
                    F[ll+1,3] = C.η[ll+1]*val
                    F[ll+1,4] = C.η[ll+1]*val
                    val = -val
                end
            end

            # Fill with P^0_ℓ(x) terms
            legendre!(C.λ, P, lmax, 0, x)

        else # abs(x) ≈ one(T)
        # Case where two points are not antipodes

            y = inv(one(T) - x*x)
            xy = x * y

            # Fill with the P^2_ℓ(x) terms initially
            legendre!(C.λ, P, lmax, 2, x)

            for ll=2:lmax
                lT = convert(T, ll)
                lp2 = lT + convert(T, 2)
                lm1 = lT - one(T)
                lm4 = lT - convert(T, 4)

                F[ll+1,3] =  C.β[ll+1] * (lp2*xy*P[ll] - (lm4*y + half*lT*lm1)*P[ll+1])
                F[ll+1,4] = 2C.β[ll+1] * (lp2*y*P[ll]  - lm1*xy*P[ll+1])
            end

            # Now refill P with the P^0_ℓ(x) terms
            legendre!(C.λ, P, lmax, 0, x)

            # Compute the F10 terms with the P as is
            for ll=2:lmax
                lm1 = convert(T, ll) - one(T)
                F[ll+1,2] = C.α[ll+1] * (xy*P[ll] - (y + half*lm1)*P[ll+1])
            end

        end # abs(x) ≈ one(T)

        # Now finally apply the normalization to the P^0_ℓ(x) function
        for ll=0:lmax
            P[ll+1] = C.η[ll+1] * P[ll+1]
        end

    end # @inbounds

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

    coeff::PixelCovarianceCoeff{Float64}
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

        coeff = PixelCovarianceCoeff{Float64}(lmax)
        F = zeros(Float64, lmax+1, 4)

        # The pixel coordinates can all be precomputed just once
        r .= pix2vec.(nside, pixels)

        return new(nside, lmax, pixels, bitfields,
                   spectra,
                   1, r, z, cij, sij, cji, sji,
                   coeff, F)
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

    fieldinds = map(f -> find(f .== SPECTRAMAP)[1], fields)
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
        c,s = bearing2(cache.r[pixind], cache.r[ii])
        cache.sij[ii] = 2*c*s
        cache.cij[ii] = c*c - s*s
        c,s = bearing2(cache.r[ii], cache.r[pixind])
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
        PixelCovarianceF!(cache.coeff, cache.F, cache.lmax, z)

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

