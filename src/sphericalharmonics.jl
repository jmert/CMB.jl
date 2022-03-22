module SphericalHarmonics

# While under development, do not export these functions
#export analyze, synthesize

import ..Legendre
import ..Legendre: λlm!, unsafe_legendre!

using FFTW
using LinearAlgebra: mul!
@static if VERSION < v"1.6.0-DEV.1591"
    using Compat: cispi # Compat@v3.25
end

function centered_range(start, stop, length)
    rng = range(start, stop, length=length+1)
    Δ = step(rng)
    return range(start+Δ/2, stop-Δ/2, length=length)
end

# Brute-force synthesis of the real-space points at (θ, ϕ) given a set of alms.
# The idea is that there are no optimizations, so it is much easier to verify the
# correctness of the function at the expense of a huge hit to performance. This will
# primarily only be useful for testing development of the more specialized algorithms.
function synthesize_reference(alms::AbstractMatrix{C}, θ, ϕ) where {C<:Complex}
    axes(θ) == axes(ϕ) ||
        throw(DimensionMismatch("`θ` and `ϕ` must have same axes"))
    lmax, mmax = size(alms) .- 1

    R = real(C)
    Λ = Matrix{R}(undef, size(alms)...)
    Φ = Vector{C}(undef, mmax + 1)
    syn = zeros(R, axes(θ)...)

    @inbounds for I in eachindex(syn)
        λlm!(Λ, lmax, mmax, cos(θ[I]))       # λ_ℓ^m(cos θ) factors
        @. Φ = cispi.((0:mmax) .* (ϕ[I]/π))  # e^{imϕ} factors
                                             #   Using cispi(ϕ/π) rather than cis(ϕ) gives
                                             #   slightly more accurate results for test
                                             #   healpix rings.

        acc = zero(R)
        # Σ_{ℓ = 0}^{ℓmax}
        for ℓ in 0:lmax
            # Σ_{m = 0}^{mmax}
            acc += real(alms[ℓ+1,1] * Λ[ℓ+1,1])
            for m in 1:min(ℓ, mmax)
                # Assuming alms were sourced from real field, alms are constrained such
                # that
                #    a_{ℓ(-m)} Y_ℓ^(-m) + a_{ℓ(+m)} Y_ℓ^(+m)
                #    == 2 * Re[ a_{ℓ(+m)} Y_ℓ^(+m) ]
                acc += 2 * real(alms[ℓ+1,m+1] * Λ[ℓ+1,m+1] * Φ[m+1])
            end
        end
        syn[I] = acc
    end
    return syn
end

# Similar to above, a brute-force analysis of the given map (with corresponding coordinates
# (θ, ϕ) for each pixel into the harmonic coefficients up to order lmax and degree mmax.
#
# Note that there is no information provided on the area (in terms of the integral's
# surface area element) of pixels, so that normalization factor is not applied during
# analysis and must be applied by the caller using information it will have.
#
# For instance, an ECP grid needs an extra normalization factor of
#
#   alms .*= 2π/N_ϕ * π/N_θ
#        .*= 2π²/N_pix
#
# Without quadrature weights, an iterative approach will be necessary to achieve any
# kind of relative accuracy:
#
#   lmax = size(alms, 1) - 1
#   map = synthesize_ecp(alms, Nθ, Nϕ)
#   θ, ϕ = make_grid(Nθ, Nϕ)
#
#   alms′ = (2π^2 / prod(size(map))) * analyze_reference(map, θ, ϕ, lmax)
#   map′ = synthesize_ecp(alms′, Nθ, Nϕ)
#   δalms = (2π^2 / prod(size(map))) * analyze_reference(map - map′, θ, ϕ, lmax)
#   alms += δalms  # iterative refinement
function analyze_reference(map, θ, ϕ, lmax::Integer, mmax::Integer = lmax)
    axes(map) == axes(θ) == axes(ϕ) ||
        throw(DimensionMismatch("`map`, `θ`, and `ϕ` must all have the same axes"))

    R = eltype(map)
    C = complex(R)
    Λ = Matrix{R}(undef, lmax + 1, mmax + 1)
    Φ = Vector{C}(undef, mmax + 1)
    alms = zeros(C, lmax + 1, mmax + 1)

    @inbounds for I in eachindex(map)
        sθ, cθ = sincos(θ[I])
        λlm!(Λ, lmax, mmax, cθ)             # λ_ℓ^m(cos θ) factors
        @. Φ = cispi((0:mmax) * (-ϕ[I]/π))  # e^{-imϕ} factors
                                            #   Using cispi(ϕ/π) rather than cis(ϕ) gives
                                            #   slightly more accurate results for test
                                            #   healpix rings.

        # Σ_{ℓ = 0}^{ℓmax}
        for ℓ in 0:lmax
            # Σ_{m = 0}^{mmax}
            alms[ℓ+1,1] += map[I] * Λ[ℓ+1,1] * sθ
            for m in 1:min(ℓ, mmax)
                alms[ℓ+1,m+1] += map[I] * Λ[ℓ+1,m+1] * Φ[m+1] * sθ
            end
        end
    end
    return alms
end

"""
    alias_index(len::Int, m::Int) -> (i::Int, isconj::Bool, isnyq::Bool)

Aliases the given 0-indexed `m` for a periodic FFT frequency axis of length `len` --- with
an additional aliasing to below the Nyquist frequency approprate for use in a real-only
[`irfft()`] transform --- to the new index `i`.

To handle symmetries of the FFT, the flags `isconj` and `isnyq` indicate if the Fourier
coefficients at `i` should be conjugated or real-doubled (i.e. `c = complex(2real(c))`),
respectively.
"""
@inline function alias_index(len::Int, i::Int)
    isconj, isnyq = false, false
    nyq = max(1, len ÷ 2)
    i < nyq && return (i, isconj, isnyq)
    i = mod(i, len)
    if i > nyq
        i = len - i
        isconj = true
    elseif i == 0 || (iseven(len) && i == nyq)
        isnyq = true
    end
    return (i, isconj, isnyq)
end

"""
    alias_coeffs(coeffs, isconj::Bool, isnyq::Bool) -> coeffs

Helper for [`alias_index`](@ref) which applies the conjugation or pure-real symmetry
transformations to the Fourier coefficients `coeffs` based on the two flags `isconj` and
`isnyq`.
"""
@inline function alias_coeffs(coeffs, isconj::Bool, isnyq::Bool)
    return ifelse(isnyq, complex.(2 .* real.(coeffs)),
                  ifelse(isconj, conj.(coeffs), coeffs))
end

# Synthesize alms to an equidistant cylindrical grid covering the entire sphere.
# More advanced that `synthesize_reference` by including:
#   * Iso-latitude Legendre optimization
#   * FFT-based ring synthesis
#   * Polar-symmetry optimization
function synthesize_ecp(alms::Matrix{C}, nθ::Integer, nϕ::Integer) where {C<:Complex}
    R = real(C)
    lmax, mmax = size(alms) .- 1
    nϕr = nϕ ÷ 2 + 1 # real-symmetric FFT's Nyquist length (index)
    nθh = (nθ + 1) ÷ 2 # number of rings in northern hemisphere

    # pixel grid definition for ECP
    θr = centered_range(R(0.0), R(π), nθ)
    ϕ₀ = R(π) / nϕ

    ecp = Matrix{R}(undef, nθ, nϕ)
    Λ = zeros(R, lmax + 1, mmax + 1)
    Λw = Legendre.Work(λlm!, Λ, zeros(R))
    λ₁ = zeros(C, nϕr)  # northern ring
    λ₂ = zeros(C, nϕr)  # southern ring

    # FFTW plan built for particular alignment/length, so not all @view(ecp[:,j]) are
    # valid (especially if nϕ is odd). Instead, synthesize into a fixed array and copy
    # back to the output buffer.
    F = plan_brfft(λ₁, nϕ)
    r = Vector{R}(undef, nϕ)

    @inbounds for j in 1:nθh
        j′ = nθ - j + 1
        θ = θr[j]
        unsafe_legendre!(Λw, Λ, lmax, mmax, cos(θ))

        for m in 0:mmax
            acc₁, acc₂ = zero(C), zero(C)
            for ℓ in m:lmax
                term = alms[ℓ+1,m+1] * Λ[ℓ+1,m+1]
                acc₁ += term
                acc₂ += isodd(ℓ + m) ? -term : term
            end
            acc₁, acc₂ = (acc₁, acc₂) .* cis(m * ϕ₀)
            i, isconj, isnyq = alias_index(nϕ, m)
            acc₁, acc₂ = alias_coeffs((acc₁, acc₂), isconj, isnyq)
            λ₁[i+1] += acc₁
            λ₂[i+1] += acc₂
        end

        copyto!(@view(ecp[j,:]), mul!(r, F, λ₁))
        fill!(λ₁, zero(C))

        copyto!(@view(ecp[j′,:]), mul!(r, F, λ₂))
        fill!(λ₂, zero(C))
    end
    return ecp
end

# Analyze an equidistant cylindrical map covering the entire sphere.
# More advanced that `analyze_reference` by including:
#   * Iso-latitude Legendre optimization
#   * FFT-based ring synthesis
#   * Polar-symmetry optimization
function analyze_ecp(ecp::Matrix{R}, lmax::Integer, mmax::Integer = lmax) where {R<:Real}
    C = complex(R)
    nθ, nϕ = size(ecp)
    nϕr = nϕ ÷ 2 + 1 # real-symmetric FFT's Nyquist length (index)
    nθh = (nθ + 1) ÷ 2 # number of rings in northern hemisphere

    # pixel grid definition for ECP
    θr = centered_range(R(0.0), R(π), nθ)
    ϕ₀ = R(π) / nϕ
    ΔΩ = 2R(π)^2 / (nθ * nϕ)

    alms = fill(zero(C), lmax + 1, mmax + 1)
    Λ = zeros(R, lmax + 1, mmax + 1)
    Λw = Legendre.Work(λlm!, Λ, zero(R))
    f₁ = zeros(C, nϕr)  # northern ring
    f₂ = zeros(C, nϕr)  # southern ring

    # FFTW plan built for particular alignment/length, so not all @view(ecp[:,j]) are
    # valid (especially if nϕ is odd). Instead, copy from map to temporary vector and
    # analyze a fixed array.
    r = zeros(R, nϕ)
    F = plan_rfft(r, 1)

    @inbounds for j in 1:nθh
        j′ = nθ - j + 1
        sθ, cθ = sincos(θr[j])
        unsafe_legendre!(Λw, Λ, lmax, mmax, cθ)

        fill!(f₁, zero(C))
        mul!(f₁, F, copyto!(r, @view(ecp[j,:])))

        fill!(f₂, zero(C))
        j′ != j && mul!(f₂, F, copyto!(r, @view(ecp[j′,:])))

        sθΔΩ = sθ * ΔΩ
        for m in 0:mmax
            i, isconj, isnyq = alias_index(nϕ, m)
            a₁, a₂ = alias_coeffs((f₁[i+1], f₂[i+1]), isconj, isnyq)
            a₁, a₂ = (a₁, a₂) .* (sθΔΩ * cis(m * -ϕ₀))
            for ℓ in m:lmax
                c = isodd(ℓ+m) ? a₁ - a₂ : a₁ + a₂
                alms[ℓ+1,m+1] += c * Λ[ℓ+1,m+1]
            end
        end
    end
    return alms
end

synthesize_ring(alms, θ₀, nϕ, ϕ₀) = synthesize_ring(alms, θ₀, ϕ₀, nϕ::Integer, Val(false))
function synthesize_ring(alms::Matrix{C}, θ₀, ϕ₀, nϕ::Integer,
                         ::Val{pair}) where {C<:Complex, pair}
    R = real(C)
    lmax, mmax = size(alms) .- 1
    nϕr = nϕ ÷ 2 + 1 # real-symmetric FFT's Nyquist length (index)

    Λ = zeros(R, lmax + 1, mmax + 1)
    r₁ = Vector{R}(undef, nϕ)
    λ₁ = zeros(C, nϕr)
    if pair
        r₂ = Vector{R}(undef, nϕ)
        λ₂ = zeros(C, nϕr)
    end
    F = plan_brfft(λ₁, nϕ)

    λlm!(Λ, lmax, mmax, cos(θ₀))
    @inbounds for m in 0:mmax
        acc₁ = zero(C)
        if pair
            acc₂ = zero(C)
        end
        for ℓ in m:lmax
            term = alms[ℓ+1,m+1] * Λ[ℓ+1,m+1]
            acc₁ += term
            if pair
                acc₂ += isodd(ℓ + m) ? -term : term
            end
        end
        i, isconj, isnyq = alias_index(nϕ, m)
        rot = cis(m * ϕ₀)
        acc₁ = alias_coeffs(acc₁ * rot, isconj, isnyq)
        if pair
            acc₂ = alias_coeffs(acc₂ * rot, isconj, isnyq)
        end
        λ₁[i+1] += acc₁
        if pair
            λ₂[i+1] += acc₂
        end
    end
    mul!(r₁, F, λ₁)
    if pair
        mul!(r₂, F, λ₂)
        return (r₁, r₂)
    end
    return r₁
end

end
