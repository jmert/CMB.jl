module Harmonics

# While under development, do not export these functions
#export analyze, synthesize

using ..Legendre
using FFTW
using LinearAlgebra: mul!

function centered_range(start, stop, length)
    rng = range(start, stop, length=length+1)
    Δ = step(rng)
    return range(start+Δ/2, stop-Δ/2, length=length)
end

# Brute-force synthesis of the real-space points at (θ, ϕ) given a set of alms.
# The idea is that there are no optimizations, so it is much easier to verify the
# correctness of the function at the expense of a huge hit to performance. This will
# primarily only be useful for testing development of the more specialized algorithms.
function synthesize_reference(alms::AbstractMatrix{T}, θ, ϕ) where {T<:Complex}
    R = real(T)
    axes(θ) == axes(ϕ) || throw(DimensionMismatch("θ and ϕ must have same axes"))
    lmax, mmax = size(alms) .- 1

    Λ = Matrix{R}(undef, size(alms)...)
    Φ = Vector{T}(undef, mmax + 1)
    syn = Array{R}(undef, size(θ)...)

    @inbounds for I in eachindex(syn)
        # λ_ℓ^m(cos θ) factors
        λlm!(Λ, lmax, mmax, cos(θ[I]))
        # e^{imϕ} factors
        #   Using complex(cospi(), sinpi()) rather than exp(complex()) gives slightly more
        #   accurate results for test healpix rings
        @. Φ = complex(cospi((0:mmax) * ϕ[I]/π), sinpi((0:mmax) * ϕ[I]/π))

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
            for ℓ in m:lmax
                acc += alms[ℓ+1,m+1] * Λ[ℓ+1,m+1]
            end
            # calculated aliased index
            if m >= nxr
                i = mod(m, nx)              # alias over total FFT frequency range
                i = i >= nxr ? i-nxr+1 : i  # reflect over Nyquist
                if i == 0
                    acc *= 2
                end
            else
                i = m
            end
            λ[i+1] += acc * @fastmath exp(complex(0.0, m*ϕ₀))
        end
        mul!(@view(ecp[:,y]), F, λ)
    end
    return permutedims(ecp)
end

end
