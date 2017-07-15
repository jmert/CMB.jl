"""
Collection of functions which compute the pixel-pixel covariance of the CMB
sky.

Based on equations given in Tegmark and de Oliveira-Costa (2001) *“How to
measure CMB polarization power spectra without losing information”*
arXiv:astro-ph/0012120v3
"""
module PixelCovariance
    export PixelCovarianceCoeff, PixelCovarianceF!

    # For computing the Legendre terms
    import ..Legendre: LegendreUnitCoeff, LegendreP!

    import Base.@boundscheck, Base.@propagate_inbounds

    struct PixelCovarianceCoeff{T<:Real}
        λ::LegendreUnitCoeff{T}
        η::Vector{T}
        α::Vector{T}
        β::Vector{T}

        function PixelCovarianceCoeff{T}(lmax::Integer) where T
            lmax ≥ 0 || throw(DomainError())

            λ = LegendreUnitCoeff{T}(lmax)
            η = Vector{T}(lmax+1)
            α = Vector{T}(lmax+1)
            β = Vector{T}(lmax+1)

            @inbounds for ll in 1:(lmax+1)
                lT = convert(T, ll)
                norm = convert(T, 2ll + 1) / (4π)
                invql = inv( (lT-one(T)) * lT * (lT+one(T)) * (lT+convert(T,2)) )

                η[ll] = norm
                α[ll] = 2norm * lT * sqrt(invql)
                β[ll] = 2norm * invql
            end

            return new(λ, η, α, β)
        end
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
            F10 = view(F, :, 2)
            F12 = view(F, :, 3)
            F22 = view(F, :, 4)

            # Fill with the P^2_ℓ(x) terms initially
            LegendreP!(C.λ, P, lmax, 2, x)

            # Clear out the ℓ == 0 and ℓ == 1 terms since they aren't defined, and then
            # compute the rest of the terms
            F12[1] = zero(T)
            F12[2] = zero(T)
            F22[1] = zero(T)
            F22[2] = zero(T)
            for ll=2:lmax
                lT = convert(T, ll)
                lp2 = lT + convert(T, 2)
                lm1 = lT - one(T)
                lm4 = lT - convert(T, 4)

                F12[ll+1] =  C.β[ll+1] * (lp2*xy*P[ll] - (lm4*y + half*lT*lm1)*P[ll+1])
                F22[ll+1] = 2C.β[ll+1] * (lp2*y*P[ll]  - lm1*xy*P[ll+1])
            end

            # Now refill P with the P^0_ℓ(x) terms
            LegendreP!(C.λ, P, lmax, 0, x)

            # Compute the F10 terms with the P as is
            F10[1] = zero(T)
            F10[2] = zero(T)
            for ll=2:lmax
                lm1 = convert(T, ll) - one(T)
                F10[ll+1] = C.α[ll+1] * (xy*P[ll] - (y + half*lm1)*P[ll+1])
            end

            # Now finally apply the normalization to the P^0_ℓ(x) function
            for ll=0:lmax
                P[ll+1] = C.η[ll+1] * P[ll+1]
            end
        end

        return F
    end
end

