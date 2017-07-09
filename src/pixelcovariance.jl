"""
Collection of functions which compute the pixel-pixel covariance of the CMB
sky.

Based on equations given in Tegmark and de Oliveira-Costa (2001) *“How to
measure CMB polarization power spectra without losing information”*
arXiv:astro-ph/0012120v3
"""
module PixelCovariance
    export Pl2!, F10, F10!, F12, F12!, F22, F22!

    import Base.@boundscheck, Base.@propagate_inbounds

    """
    Computes the Legendre polynomials ``P_ℓ^2(x)`` for ``0 ≤ ℓ ≤`` `lmax`.

        Pl2!(P::Vector, lmax, x) -> P
    """
    function Pl2! end

    @inline Pl2!(P::Vector{T}, lmax::Integer, x::T) where {T<:Real} =
        LegendreP!(LegendreUnitNorm(), P, lmax, 2, x)

    @inline function _chkbounds_FXX(lmax, pl, x)
        lmax >= 2 || throw(DimensionMismatch("lmax must be greater than 1"))
        (size(pl,1) ≥ lmax+1 && size(pl,2) > lmax+1) || throw(BoundsError())
    end
    @inline function _chkbounds_FXX(fxx, lmax, pl, x)
        _chkbounds_FXX(lmax, pl, x)
        size(fxx) >= (length(x),lmax+1) || throw(DimensionMismatch())
    end

    """
        F10(lmax::Integer, pl, x)
    """
    function F10 end

    """
        F10!(f10::Matrix, lmax::Integer, pl, x)
    """
    function F10! end

    @propagate_inbounds function F10(lmax::Integer, pl, x::Number)
        @boundscheck _chkbounds_FXX(lmax, pl, x)
        f10 = Array{eltype(x)}(1, lmax+1)
        @inbounds F10!(f10, pl, x)
        return reshape(f10, lmax+1)
    end

    @propagate_inbounds function F10(lmax::Integer, pl, x::Vector{T} where T<:Number)
        @boundscheck _chkbounds_FXX(lmax, pl, x)
        f10 = Array{eltype(x)}(length(x), lmax+1)
        @inbounds F10!(f10, lmax, pl, x)
        return f10
    end

    @inline function f10_kernel(ℓ, ℓℓm1by2, scale, plm1,pl,x)
        T = promote_type(typeof(plm1), typeof(pl), typeof(x))
        omxx = muladd(x, -x, one(T)) # == (1.0 - x*x)
        obyomxx = one(T) / omxx
        f10 = scale*(ℓ*x*plm1*obyomxx - (ℓ*obyomxx-ℓℓm1by2)*pl)
        return ifelse(abs(omxx)<10eps(T), zero(T), f10)
    end

    @propagate_inbounds function F10!(f10::Matrix{T1} where T1<:Number, lmax::Integer, pl, x::Union{T2,Vector{T2}} where T2<:Number)
        @boundscheck _chkbounds_FXX(f10, lmax, pl, x)

        @inbounds begin
            # unused
            f10[:,1] .= 0.0
            f10[:,2] .= 0.0

            for ℓ=2:lmax
                # Constant factors depending only on ℓ have been extracted
                # from the kernel.
                scale = 2 / sqrt((ℓ-1)*ℓ*(ℓ+1)*(ℓ+2))
                ℓℓm1by2 = ℓ*(ℓ-1) / 2
                f10[:,ℓ+1] .= f10_kernel.(ℓ, ℓℓm1by2, scale, pl[:,ℓ], pl[:,ℓ+1], x)
            end
        end

        return f10
    end

    """
        F12(lmax::Integer, pl2, x)
    """
    function F12 end

    """
        F12!(f12::Matrix, lmax::Integer, pl2, x)
    """
    function F12! end

    @propagate_inbounds function F12(lmax::Integer, pl2, x::Number)
        @boundscheck _chkbounds_FXX(lmax, pl2, x)
        f12 = Array{eltype(x)}(1, lmax+1)
        @inbounds F12!(f12, pl2, x)
        return reshape(f12, lmax+1)
    end

    @propagate_inbounds function F12(lmax::Integer, pl2, x::Vector{T} where T<:Number)
        @boundscheck _chkbounds_FXX(lmax, pl2, x)
        f12 = Array{eltype(x)}(length(x), lmax+1)
        @inbounds F12!(f12, lmax, pl2, x)
        return f12
    end

    @inline function f12_kernel(ℓ,ℓp2,ℓm4,ℓℓm1by2,scale,pl2m1,pl2,x)
        T = promote_type(typeof(pl2m1), typeof(pl2), typeof(x))
        m1toℓby2 = ifelse( ℓ>>1 == 0, 0.5, -0.5) # (-1)^ℓ / 2
        obyomxx = one(x) / muladd(x, -x, one(x))
        f12 = scale * (ℓp2*x*obyomxx*pl2m1 - (ℓm4*obyomxx + ℓℓm1by2)*pl2)
        return ifelse( abs(one(T)-x)<10eps(T), 0.5,
                      ifelse( abs(one(T)+x)<10eps(T), m1toℓby2, f12 ) )
    end

    @propagate_inbounds function F12!(f12::Matrix{T1} where T1<:Number, lmax::Integer, pl2, x::Union{T2,Vector{T2}} where T2<:Number)
        @boundscheck _chkbounds_FXX(f12, lmax, pl2, x)

        @inbounds begin
            # unused
            f12[:,1] .= 0.0
            f12[:,2] .= 0.0

            for ℓ=3:lmax
                scale = 2 / ((ℓ-1)*ℓ*(ℓ+1)*(ℓ+2))
                ℓp2 = ℓ + 2
                ℓm4 = ℓ - 4
                ℓℓm1by2 = ℓ*(ℓ-1) / 2
                f12[:,ℓ+1] .= f12_kernel.(ℓ, ℓp2, ℓm4, ℓℓm1by2, scale, pl2[:,ℓ], pl2[:,ℓ+1], x)
            end
        end

        return f12
    end

    """
        F22(lmax::Integer, pl2, x)
    """
    function F22 end

    """
        F22!(f22::Matrix, lmax::Integer, pl2, x)
    """
    function F22! end

    @propagate_inbounds function F22(lmax::Integer, pl2, x::Number)
        @boundscheck _chkbounds_FXX(lmax, pl2, x)
        f22 = Array{eltype(x)}(1, lmax+1)
        @inbounds F10!(f22, pl2, x)
        return reshape(f22, lmax+1)
    end

    @propagate_inbounds function F22(lmax::Integer, pl2, x::Vector{T} where T<:Number)
        @boundscheck _chkbounds_FXX(lmax, pl2, x)
        f22 = Array{eltype(x)}(length(x), lmax+1)
        @inbounds F22!(f22, lmax, pl2, x)
        return f22
    end

    @inline f22_kernel(ℓ,plm1,pl,x) = ifelse()

    @propagate_inbounds function F22!(f22::Matrix{T1} where T1<:Number, lmax::Integer, pl2, x::Union{T2,Vector{T2}} where T2<:Number)
        @boundscheck _chkbounds_FXX(f22, lmax, pl2, x)

        @inbounds begin
            # unused
            f22[:,1] .= 0.0
            f22[:,2] .= 0.0

            for ℓ=3:lmax
                f22[:,ℓ+1] .= f22_kernel.(ℓ, pl[:,ℓ], pl[:,ℓ+1], x)
            end
        end

        return f22
    end
end

