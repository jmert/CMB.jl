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

# Simple reference function which synthesizes field for an equidistant-cylindrical
# project (ECP) grid on the entire sphere.
function synthesize_ecp(alms::Matrix{T}, nx::Integer, ny::Integer) where {T<:Complex}
    lmax, mmax = last.(axes(alms)) .- 1
    nxr = nx ÷ 2 + 1 # real-symmetric FFT's Nyquist length (index)

    # pixel grid definition for ECP
    θ = centered_range(0.0, 1.0π, ny)
    ϕ = centered_range(0.0, 2.0π, nx)
    Δθ, Δϕ = step(θ), step(ϕ)
    θ₀, ϕ₀ = first(θ), first(ϕ)

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
