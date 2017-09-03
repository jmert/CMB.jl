"""
Collection of functions which compute the pixel-pixel covariance of the CMB
sky.

Based on equations given in Tegmark and de Oliveira-Costa (2001) *“How to
measure CMB polarization power spectra without losing information”*
arXiv:astro-ph/0012120v3
"""
module PixelCovariance

export PixelCovarianceCoeff, PixelCovarianceF!,
    PixelCovarianceCache, updatespectra!, selectpixel!,
    pixelcovariance, pixelcovariance!

# For computing the Legendre terms
import ..Legendre: LegendreUnitCoeff, LegendreP!

# For pixelcovariance() wrapper function
import ..Sphere: bearing, cosdistance
import ..Healpix: pix2phi, pix2theta

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

    half = inv(convert(T, 2))
    y = inv(one(T) - x*x)
    xy = x * y

    @inbounds begin
        P = view(F, :, 1)

        # Fill with the P^2_ℓ(x) terms initially
        LegendreP!(C.λ, P, lmax, 2, x)

        # Clear out the ℓ == 0 and ℓ == 1 terms since they aren't defined, and then
        # compute the rest of the terms
        F[1,3] = zero(T)
        F[2,3] = zero(T)
        F[1,4] = zero(T)
        F[2,4] = zero(T)
        for ll=2:lmax
            lT = convert(T, ll)
            lp2 = lT + convert(T, 2)
            lm1 = lT - one(T)
            lm4 = lT - convert(T, 4)

            F[ll+1,3] =  C.β[ll+1] * (lp2*xy*P[ll] - (lm4*y + half*lT*lm1)*P[ll+1])
            F[ll+1,4] = 2C.β[ll+1] * (lp2*y*P[ll]  - lm1*xy*P[ll+1])
        end

        # Now refill P with the P^0_ℓ(x) terms
        LegendreP!(C.λ, P, lmax, 0, x)

        # Compute the F10 terms with the P as is
        F[1,2] = zero(T)
        F[2,2] = zero(T)
        for ll=2:lmax
            lm1 = convert(T, ll) - one(T)
            F[ll+1,2] = C.α[ll+1] * (xy*P[ll] - (y + half*lm1)*P[ll+1])
        end

        # Now finally apply the normalization to the P^0_ℓ(x) function
        for ll=0:lmax
            P[ll+1] = C.η[ll+1] * P[ll+1]
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
    θ₀::Ref{Float64}
    ϕ₀::Ref{Float64}

    θ::Vector{Float64}
    ϕ::Vector{Float64}
    z::Vector{Float64}
    αij::Vector{Float64}
    αji::Vector{Float64}
    cij::Vector{Float64}
    sij::Vector{Float64}
    cji::Vector{Float64}
    sji::Vector{Float64}

    coeff::PixelCovarianceCoeff{Float64}
    F::Matrix{Float64}

    C::Matrix{Float64}
    Cv::Dict{Symbol,typeof(view(Matrix{Float64}(1,1),:,1))}

    """
        PixelCovarianceCache(nside, lmax, pixels [, fields=Symbol[:QQ,:QU,:UU] ])
    """
    function PixelCovarianceCache(nside, lmax, pixels, fields=Symbol[:QQ,:QU,:UU])
        N = length(pixels)

        # Map the symbol description of which fields to generate into a bit array
        bitfields = [(ff in fields) for ff in FIELDMAP]
        # The off-diagonal terms must be symmetric, so guarantee symmetry.
        bitfields .|= bitfields'

        C = zeros(Float64, N, 9)
        # Name views for convienience
        Cv = Dict(
                  :TT => view(C,:,1),
                  :QT => view(C,:,2),
                  :UT => view(C,:,3),
                  :TQ => view(C,:,4),
                  :QQ => view(C,:,5),
                  :UQ => view(C,:,6),
                  :TU => view(C,:,7),
                  :QU => view(C,:,8),
                  :UU => view(C,:,9),
                 )

        spectra = zeros(Float64, lmax+1, 6)

        θ = zeros(Float64, N)
        ϕ = zeros(Float64, N)
        z = zeros(Float64, N)
        αij = zeros(Float64, N)
        αji = zeros(Float64, N)
        cij = zeros(Float64, N)
        sij = zeros(Float64, N)
        cji = zeros(Float64, N)
        sji = zeros(Float64, N)

        coeff = PixelCovarianceCoeff{Float64}(lmax)
        F = zeros(Float64, lmax+1, 4)

        # The pixel coordinates can all be precomputed just once
        θ .= pix2theta.(nside, pixels)
        ϕ .= pix2phi.(nside, pixels)

        return new(nside, lmax, pixels, bitfields,
                   spectra,
                   0, 0.0, 0.0,
                   θ, ϕ, z, αij, αji, cij, sij, cji, sji,
                   coeff, F, C, Cv)
    end
end

# Improve printing somewhat
Base.show(io::IO, norm::PixelCovarianceCache) = print(io, PixelCovarianceCache)
function Base.show(io::IO, ::MIME"text/plain", C::PixelCovarianceCache)
    println(io, C)
    println(io, "    HEALPix nside: $(C.nside)")
    println(io, "    number of pixels: $(length(C.pixels))")
    println(io, "    maximum ℓ mode: $(C.lmax)")
    println(io, "    selected pixel: ", C.pixind[], ", at (θ,ϕ) = ($(C.θ₀[]), $(C.ϕ₀[]))")
    println(io, "    covariance blocks: $(reshape(FIELDMAP[C.fields],:))")
end


function updatespectra!(cache, spectra)
    cache.spectra .= spectra
end

function selectpixel!(cache, pixind)
    cache.pixind[] = pixind
    cache.θ₀[] = pix2theta(cache.nside, cache.pixels[pixind])
    cache.ϕ₀[] = pix2phi(cache.nside, cache.pixels[pixind])

    cache.z   .= cosdistance.(cache.θ₀, cache.ϕ₀, cache.θ, cache.ϕ)
    cache.αij .= bearing.(cache.θ₀, cache.ϕ₀, cache.θ,  cache.ϕ)
    cache.αji .= bearing.(cache.θ,  cache.ϕ,  cache.θ₀, cache.ϕ₀)

    @fastmath begin
        cache.cij .= cos.(2 .* cache.αij)
        cache.sij .= sin.(2 .* cache.αij)
        cache.cji .= cos.(2 .* cache.αji)
        cache.sji .= sin.(2 .* cache.αji)
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

function pixelcovariance!(cache::PixelCovarianceCache)
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
            cache.C[i,1] = tt
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
            cache.C[i,2] =  tq*cache.cij[i] + tu*cache.sij[i]
            cache.C[i,3] = -tq*cache.sij[i] + tu*cache.cij[i]
            cache.C[i,4] =  tq*cache.cji[i] + tu*cache.sji[i]
            cache.C[i,7] = -tq*cache.sji[i] + tu*cache.cji[i]
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
            cache.C[i,5] =  qq*cij*cji + qu*(cij*sji+sij*cji) + uu*sij*sji
            cache.C[i,6] = -qq*sij*cji + qu*(cij*cji-sij*sji) + uu*cij*sji
            cache.C[i,8] = -qq*cij*sji + qu*(cij*cji-sij*sji) + uu*sij*cji
            cache.C[i,9] =  qq*sij*sji - qu*(cij*sji+sij*cji) + uu*cij*cji
        end
    end

    return cache.Cv
end

end # module PixelCovariance

