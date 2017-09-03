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

struct PixelCovarianceCache
    nside::Int
    lmax::Int
    pixels::Vector{Int}

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

    function PixelCovarianceCache(nside,lmax,pixels)
        N = length(pixels)

        spectra = Matrix{Float64}(lmax+1, 3)

        θ = Vector{Float64}(N)
        ϕ = similar(θ)

        z = similar(θ)

        αij = similar(θ)
        αji = similar(θ)
        cij = similar(θ)
        sij = similar(θ)
        cji = similar(θ)
        sji = similar(θ)

        coeff = PixelCovarianceCoeff{Float64}(lmax)
        F = Matrix{Float64}(lmax+1, 4)

        C = Matrix{Float64}(N, 9)
        Cv = Dict(
                :TT => view(C, :, 1),
                :QT => view(C, :, 2),
                :UT => view(C, :, 3),
                :TQ => view(C, :, 4),
                :QQ => view(C, :, 5),
                :UQ => view(C, :, 6),
                :TU => view(C, :, 7),
                :QU => view(C, :, 8),
                :UU => view(C, :, 9)
            )

        # The pixel coordinates can all be precomputed just once
        θ .= pix2theta.(nside, pixels)
        ϕ .= pix2phi.(nside, pixels)

        return new(nside, lmax, pixels,
                   spectra,
                   0, 0.0, 0.0,
                   θ, ϕ, z, αij, αji, cij, sij, cji, sji,
                   coeff, F, C, Cv)
    end
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
    @inbounds for (i,z) in enumerate(cache.z)
        PixelCovarianceF!(cache.coeff, cache.F, cache.lmax, z)
        tt = zero(eltype(cache.F))
        qq = zero(eltype(cache.F))
        uu = zero(eltype(cache.F))
        for ll in size(cache.F,1):-1:1
            tt += cache.spectra[ll,1]*cache.F[ll,1]
            qq += cache.spectra[ll,2]*cache.F[ll,3] - cache.spectra[ll,3]*cache.F[ll,4]
            uu += cache.spectra[ll,3]*cache.F[ll,3] - cache.spectra[ll,2]*cache.F[ll,4]
        end
        cache.C[i,1] =  tt
        cache.C[i,5] =  cache.cij[i]*qq*cache.cji[i] + cache.sij[i]*uu*cache.sji[i]
        cache.C[i,9] =  cache.sij[i]*qq*cache.sji[i] + cache.cij[i]*uu*cache.cji[i]
        cache.C[i,6] = -cache.sij[i]*qq*cache.cji[i] + cache.cij[i]*uu*cache.sji[i]
        cache.C[i,8] = -cache.cij[i]*qq*cache.sji[i] + cache.sij[i]*uu*cache.cji[i]
    end

    return cache.Cv
end

end # module PixelCovariance

