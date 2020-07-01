module SphericalHarmonics

# While under development, do not export these functions
#export analyze, synthesize

using ..Legendre
using FFTW
using LinearAlgebra: mul!

@static if isdefined(Base, :sincospi)
    import Base.sincospi
else
    sincospi(x) = (sinpi(x), cospi(x))
end
@static if isdefined(Base, :cispi)
    import Base.cispi
else
    cispi(x) = complex(reverse(sincospi(x))...)
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

#= TODO: Doesn't work yet --- figure out quadrature weighting

# Similar to above, a brute-force analysis of the given map (with corresponding coordinates
# (θ, ϕ) for each pixel into the harmonic coefficients up to order lmax and degree mmax.
# Note that there is no information provided on the area (in terms of the integral's
# surface area element) of pixels, so that normalization factor is not applied during
# analysis and must be applied by the caller using information it will have.
# Again, there are essentially no optimizations made here to make the implementation easier
# to ensure it is correct, so this should only be used in testing/development.
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
        λlm!(Λ, lmax, mmax, cθ)               # λ_ℓ^m(cos θ) factors
        @. Φ = cispi.((0:mmax) .* (-ϕ[I]/π))  # e^{-imϕ} factors
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
=#

function synthesize_ring(alms::AbstractMatrix{T}, nx, ϕ₀, θ) where {T<:Complex}
    R = real(T)
    lmax, mmax = size(alms, 1) - 1, size(alms, 2) - 1
    nxr = nx ÷ 2 + 1

    Λ = zeros(real(T), lmax+1, mmax+1)
    λ = zeros(T, nxr)
    F = plan_brfft(λ, nx)
    ring = Vector{R}(undef, nx)

    λlm!(Λ, lmax, mmax, cospi(θ/π))
    fill!(λ, zero(T))

    for m in 0:mmax
        acc = zero(T)
        i, mul2, nc = alias_index(nx, m)
        for ℓ in m:lmax
            acc += alms[ℓ+1,m+1] * Λ[ℓ+1,m+1]
        end
        acc *= cis(m * ϕ₀)
        acc = mul2 ? complex(2real(acc)) : nc ? conj(acc) : acc
        λ[i+1] += acc
    end
    mul!(ring, F, λ)
end

# alias_index(len, i) -> (i′, mul2, nc)
#
# Aliases the given index `i` for a periodic FFT frequency axis of length `len`, with an
# additional aliasing to below the Nyquist frequency approprate for use in a real-only
# [irfft()] transform.
#
# - `i′` is the returned aliased index.
# - `mul2` is true if `i′` corresponds to a real-only Fourier mode (`i′ == 0` or
#   (`iseven(len) && i′ == len÷2+1`) and the original mode must be doubled to account for
#   the multiplicity of the aliased mode.
# - `nc` is true if `i′` corresponds to a negative-frequency mode in the full-length
#   (not real-only, half-length) spectrum and therefore should be conjugated before
#   summing into the aliased axis.
@inline function alias_index(len, i)
    nyq = len ÷ 2
    i′, nc = i, false
    if i′ > nyq
        i′ = mod(i′, len)
        i′, nc = i′ > nyq ? (len - i′, true) : (i′, false)
    end
    mul2 = i != 0 && i′ == 0 || (iseven(len) && i′ == nyq)
    return (i′, mul2, nc)
end

# Simple reference function which synthesizes field for an equidistant-cylindrical
# project (ECP) grid on the entire sphere.
function synthesize_ecp(alms::Matrix{T}, nx::Integer, ny::Integer) where {T<:Complex}
    lmax, mmax = size(alms) .- 1
    nxr = nx ÷ 2 + 1 # real-symmetric FFT's Nyquist length (index)

    # pixel grid definition for ECP
    θ = centered_range(0.0, 1.0π, ny)
    ϕ = centered_range(0.0, 2.0π, nx)
    ϕ₀ = first(ϕ)

    ecp = Matrix{real(T)}(undef, nx, ny) # transposed to have x-dim have stride 1
    Λ = zeros(real(T), lmax+1, mmax+1)
    λ = zeros(T, nxr)
    F = plan_brfft(λ, nx)

    for y in 1:ny
        λlm!(Λ, lmax, mmax, cos(θ[y]))
        fill!(λ, zero(T))

        for m in 0:mmax
            acc = zero(T)
            i, mul2, nc = alias_index(nx, m)
            for ℓ in m:lmax
                acc += alms[ℓ+1,m+1] * Λ[ℓ+1,m+1]
            end
            acc *= cis(m * ϕ₀)
            acc = mul2 ? complex(2real(acc)) : nc ? conj(acc) : acc
            λ[i+1] += acc
        end
        mul!(@view(ecp[:,y]), F, λ)
    end
    return permutedims(ecp)
end

end
