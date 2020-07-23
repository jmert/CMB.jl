module PixelCovariance

export Fweights!,
    pixelcovariance, pixelcovariance!

import ..Sphere: bearing2, cosdistance
import ..Healpix: pix2vec
import ..unchecked_sqrt

using BitFlags
using Legendre
using Legendre: unsafe_legendre!
using StaticArrays

import Base: @propagate_inbounds, checkindex, checkbounds_indices, OneTo, Slice

struct FweightsWork{T,N,S<:AbstractArray{T},V<:AbstractArray{T}}
    legwork::Legendre.Work{T,N,V}
    P::S
    x::V
    y::V
    xy::V
end
function FweightsWork(norm::AbstractLegendreNorm, F, z)
    T = promote_type(eltype(norm), eltype(F), eltype(z))
    legwork = Legendre.Work(norm, F, z)
    P = view(F, map(Slice, axes(z))..., :, 1)
    x = legwork.z
    y = similar(x)
    xy = similar(x)
    return FweightsWork(legwork, P, x, y, xy)
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
        work = FweightsWork(workornorm, F′, x′)
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

    P  = work.P
    x  = work.x
    y  = work.y
    xy = work.xy
    legwork = work.legwork
    T = eltype(x)
    half = one(T) / 2

    Is = map(Slice, axes(x))
    I = CartesianIndices(Is)

    @simd for ii in I
        x[ii] = x′ = convert(T, z[ii])
        y′ = inv(fma(-x′, x′, one(T)))
        y[ii], xy[ii] = y′, x′ * y′
    end

    # Calculations below implicitly assume the out-of-bounds entries have zero value, so
    # also clear out the first to ℓ entries.
    fill!(@view(F[Is..., OneTo(2), :]), zero(T))

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
@bitflag CovarianceFields TT QT UT TQ QQ UQ TU QU UU NO_FIELD=0
const TPol = TQ | TU | UT | QT
const Pol  = QQ | UU | QU | UQ

function minrow(fields::CovarianceFields)
    F = Integer(fields)
    f = (F | (F >> 0x3) | (F >> 0x6)) & 0x07
    return trailing_zeros(f) + 0x1
end

function mincol(fields::CovarianceFields)
    F = Integer(fields)
    f = (F | (F >> 0x1) | (F >> 0x2)) & 0b001001001 #= 0x49 =#
    f = (f | (f >> 0x2) | (f >> 0x4)) & 0b111 #= 0x07 =#
    return trailing_zeros(f) + 0x1
end

function pixelcovariance(pix::AbstractVector, Cl::AbstractMatrix, fields::CovarianceFields)
    T = eltype(first(pix))
    npix = length(pix)
    cov = zeros(T, npix, 9, npix)
    lmax = size(Cl, 1) - 1
    work = PixelCovWork{T}(lmax, LegendreUnitNorm())
    for ii in axes(pix, 1)
        @inbounds _pixelcovariance_impl!(work, view(cov, :, :, ii), pix, ii, Cl, fields)
    end
    cov = reshape(permutedims(reshape(cov, 3npix, 3, npix), (1, 3, 2)), 3npix, 3npix)
    return cov
end

function pixelcovariance!(cov::AbstractMatrix, pix::AbstractVector, pixind,
                          Cl::AbstractMatrix, fields::CovarianceFields)
    axes(cov, 1) == axes(pix, 1) ||
        throw(DimensionMismatch("Axes of covariance matrix and pixel vector do not match"))
    size(cov, 2) == 9 ||
        throw(DimensionMismatch("Output array expected to have 9 columns"))
    size(Cl, 2) == 6 || throw(DimensionMismatch("Expected 6 Cl spectra"))
    Base.checkbounds(pix, pixind)
    Base.require_one_based_indexing(Cl)

    return unsafe_pixelcovariance!(LegendreUnitNorm(), cov, pix, pixind, Cl, fields)
end

struct PixelCovWork{T,N,V}
    F::Matrix{T}
    Fwork::FweightsWork{T,N,V}
end
function PixelCovWork{T}(lmax::Int, norm::AbstractLegendreNorm = LegendreUnitNorm()) where {T}
    F = Matrix{T}(undef, lmax+1, 4)
    Fwork = FweightsWork(norm, F, Legendre.Scalar{T}())
    return PixelCovWork(F, Fwork)
end
Base.eltype(::PixelCovWork{T}) where {T} = T

function unsafe_pixelcovariance!(workornorm::Union{AbstractLegendreNorm,PixelCovWork},
                                 cov::AbstractMatrix, pix::AbstractVector, pixind::Int,
                                 Cl::AbstractMatrix, fields::CovarianceFields)
    if workornorm isa AbstractLegendreNorm
        T = promote_type(eltype(workornorm), eltype(cov), eltype(first(pix)))
        lmax = size(Cl, 1) - 1
        work = PixelCovWork{T}(lmax, workornorm)
        @inbounds _pixelcovariance_impl!(work, cov, pix, pixind, Cl, fields)
    else
        @inbounds _pixelcovariance_impl!(workornorm, cov, pix, pixind, Cl, fields)
    end
    return cov
end

@propagate_inbounds function _pixelcovariance_impl!(work::PixelCovWork{T},
        cov::AbstractMatrix, pix::AbstractVector, pixind::Integer,
        Cl::AbstractMatrix, fields::CovarianceFields) where {T}
    F = work.F
    Fwork = work.Fwork
    lmax = size(F, 1) - 1
    fourpi = 4 * convert(T, π)

    # N.B. In general, the spectrum causes Cl*F to decrease rapidly as ℓ → ∞; reverse
    #      the order of summation to accumulate from smallest to largest magnitude values.
    # TODO: Look into using and compare this assumption against Julia's built-in sum() which
    #       uses a divide-and-conquered summation to increase numerical precision in general
    #       cases.
    R = reverse(axes(F, 1))

    jj = pixind
    for ii in axes(pix, 1)
        z = cosdistance(pix[ii], pix[jj])
        if fields & (TPol | Pol) != NO_FIELD
            c, s = bearing2(pix[ii], pix[jj])
            sij, cij = 2*c*s, c*c - s*s
            c, s = bearing2(pix[jj], pix[ii])
            sji, cji = 2*c*s, c*c - s*s
        else
            sij, cij = zero(T), zero(T)
            sji, cji = zero(T), zero(T)
        end

        unsafe_Fweights!(Fwork, F, lmax, z)

        # TT
        if fields & TT != NO_FIELD
            tt = zero(T)
            @simd for ll in R
                ClTT = Cl[ll,1]
                F00 = F[ll, 1]
                tt = muladd(ClTT, F00, tt) # tt + ClTT*F00
            end
            tt /= fourpi
            cov[ii,1] = tt # TT
        end

        # TQ and TU
        if fields & TPol != NO_FIELD
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
            cov[ii,2] =  tq*cij + tu*sij # QT
            cov[ii,3] = -tq*sij + tu*cij # QU
            cov[ii,4] =  tq*cji + tu*sji # TQ
            cov[ii,7] = -tq*sji + tu*cji # TU
        end

        # QQ, QU, and UU
        if fields & Pol != NO_FIELD
            qq = zero(T)
            qu = zero(T)
            uu = zero(T)
            @simd for ll in R
                ClEE = Cl[ll,2]
                ClBB = Cl[ll,3]
                ClEB = Cl[ll,6]
                F12 = F[ll,3]
                F22 = F[ll,4]

                qq = muladd(ClEE, F12, muladd(-ClBB, F22, qq)) # qq + ClEE*F12 - ClBB*F22
                uu = muladd(ClBB, F12, muladd(-ClEE, F22, uu)) # uu + ClBB*F12 - ClEE*F22
                qu = muladd(ClEB, F12 + F22, qu)               # qu + ClEB*(F12 + F22)
            end
            qq /= fourpi
            uu /= fourpi
            qu /= fourpi
            cov[ii,5] =  qq*cij*cji + qu*(cij*sji+sij*cji) + uu*sij*sji # QQ
            cov[ii,6] = -qq*sij*cji + qu*(cij*cji-sij*sji) + uu*cij*sji # UQ
            cov[ii,8] = -qq*cij*sji + qu*(cij*cji-sij*sji) + uu*sij*cji # QU
            cov[ii,9] =  qq*sij*sji - qu*(cij*sji+sij*cji) + uu*cij*cji # UU
        end
    end

    return cov
end

end # module PixelCovariance
