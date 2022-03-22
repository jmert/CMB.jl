using CMB.SphericalHarmonics
using CMB.SphericalHarmonics: centered_range,
    analyze_reference, synthesize_reference,
    analyze_ecp, synthesize_ecp,
    synthesize_ring

import Random


"""
    gen_alms(::Type{T}, lmax::Integer, mmax::Integer = lmax; seed = nothing) where T

Generates random alms up to ``(ℓ,m) = (lmax, mmax)`` of type `complex(T)`. Note that the
DC (``Y_{00}``) term is not set.

A seed value may be passed as a keyword argument to reset the random number generator before
generating alms. Values are drawn such that the return values for ``lmax < lmax′`` given
the same seed will be equal for all ``(ℓ, m)`` where ``ℓ ≤ lmax``.
"""
function gen_alms(::Type{T}, lmax::Integer, mmax::Integer = lmax; seed = nothing) where T
    if seed !== nothing
        Random.seed!(seed)
    end
    alms = zeros(complex(T), lmax + 1, mmax + 1)
    @inbounds for l in 1:lmax
        alms[l+1,1] = randn(real(T))
        alms[l+1,2:l+1] .= randn.(complex(T))
    end
    return alms
end

"""
    zerotol!(X)

Simple helper to round very small values to zero. This applies both to the absolute value
or to the real/imaginary components of a complex number.
"""
function zerotol!(X)
    ϵ = sqrt(eps(maximum(abs, X)))
    @inbounds for ii in eachindex(X)
        v = X[ii]
        if abs(v) < ϵ
            X[ii] = zero(eltype(X))
        elseif abs(real(v)) < ϵ
            X[ii] = imag(v)*im
        elseif abs(imag(v)) < ϵ
            X[ii] = real(v)
        end
    end
    return X
end


# Define the analytic expressions for the first few spherical harmonics; verify
# implementations against these.
Y00(θ, ϕ) =  sqrt(1 / 4oftype(θ, π)) * complex(true)
Y10(θ, ϕ) =  sqrt(3 / 4oftype(θ, π)) * cos(θ) * complex(true)
Y11(θ, ϕ) = -sqrt(3 / 8oftype(θ, π)) * sin(θ) * cis(ϕ)
Y20(θ, ϕ) =  sqrt(5 / 16oftype(θ, π)) * (3cos(θ)^2 - 1) * complex(true)
Y21(θ, ϕ) = -sqrt(15 / 8oftype(θ, π)) * sin(θ) * cos(θ) * cis(ϕ)
Y22(θ, ϕ) =  sqrt(15 / 32oftype(θ, π)) * sin(θ)^2 * cis(2ϕ)

n = 50
θ = repeat(centered_range(0.0, 1.0π, n), 1, 2n)
ϕ = repeat(centered_range(0.0, 2.0π, 2n)', n, 1)
θ′ = repeat(centered_range(0.0, 1.0π, n+1), 1, 2n+1)
ϕ′ = repeat(centered_range(0.0, 2.0π, 2n+1)', n+1, 1)

# Generate complex fields from running real-only synthesis on real and imaginary
# components separately.
function synthesize_reference_complex(ℓ, m, θ, ϕ, T::Type = Float64)
    alms = zeros(Complex{T}, ℓ + 1, m + 1)
    # Assign unit power to single delta (ℓ,m) mode.
    #
    # Real part
    alms[ℓ+1,m+1] = 1
    ref = synthesize_reference(alms, θ, ϕ)
    # Complex part -- Use (-im) not (+im) since the goal is to swap the internal
    #                 (a + ib) to (b + ia), which requires multiplication by -im.
    alms[ℓ+1,m+1] *= -im
    ref += synthesize_reference(alms, θ, ϕ) .* im
    # account for doubling due to assumption of real-only symmetry in synthesize_reference
    ref ./= m == 0 ? 1 : 2
    return ref
end

# Because of the disparity between analysis as integration of a field and summation over
# a discrete map, a single round of analyzing any given map suffers from non-floating
# point errors. By iterating and repeatedly analyzing the residual between input and
# analyze+(re)synthesize, the harmonic coefficient slowly converge to what we'd get if
# we had much higher resolution input maps.
#
# Generically, this is a problem that falls under the umbrella of "quadrature weights", but
# that's a complex topic that we bypass by just brute-force converging our solution.
function analyze_reference_complex_iter(map, θ, ϕ, lmax, mmax)
    R = eltype(θ)
    tol = sqrt(eps(maximum(abs, map)))
    # norm is the correct normalization for an ECP pixelization
    norm = 2R(π)^2 / length(θ)

    function iter(ecp)
        local alms
        alms = norm .* analyze_reference(ecp, θ, ϕ, lmax, mmax)
        for ii in 2:10  # allow up to 10 iterations before bailing
            ecp′ = synthesize_reference(alms, θ, ϕ)
            δecp = ecp - ecp′
            # check for convergence
            if maximum(abs, δecp) < tol
                break
            end
            alms .+= norm .* analyze_reference(δecp, θ, ϕ, lmax, mmax)
        end
        return alms
    end
    # analyze each of real and imaginary parts separately since the routines require
    # real inputs.
    alms = iter(real.(map)) .+ im .* iter(imag.(map))
    return zerotol!(alms)
end

@testset "Analytic checks — synthesis (reference)" begin
    # Validates the baseline reference implementation
    @test Y00.(θ, ϕ) ≈ synthesize_reference_complex(0, 0, θ, ϕ)
    @test Y10.(θ, ϕ) ≈ synthesize_reference_complex(1, 0, θ, ϕ)
    @test Y11.(θ, ϕ) ≈ synthesize_reference_complex(1, 1, θ, ϕ)
    @test Y20.(θ, ϕ) ≈ synthesize_reference_complex(2, 0, θ, ϕ)
    @test Y21.(θ, ϕ) ≈ synthesize_reference_complex(2, 1, θ, ϕ)
    @test Y22.(θ, ϕ) ≈ synthesize_reference_complex(2, 2, θ, ϕ)
end

@testset "Analytic checks — analysis (reference)" begin
    # Validates the baseline reference implementation
    lmax, mmax = 10, 3  # beyond Y22 just as a mild check of any obvious problems
    analyze(Y) = analyze_reference_complex_iter(Y.(θ, ϕ), θ, ϕ, lmax, mmax)
    expect(ℓ, m) = setindex!(zeros(complex(eltype(θ)), lmax+1, mmax+1), 1.0, ℓ+1, m+1)
    @test analyze(Y00) ≈ expect(0, 0)
    @test analyze(Y10) ≈ expect(1, 0)
    @test analyze(Y11) ≈ expect(1, 1)
    @test analyze(Y20) ≈ expect(2, 0)
    @test analyze(Y21) ≈ expect(2, 1)
    @test analyze(Y22) ≈ expect(2, 2)
end


function synthesize_ecp_complex(ℓ, m, nθ, nϕ, T::Type = Float64)
    alms = zeros(Complex{T}, ℓ + 1, m + 1)
    # Assign unit power to single delta (ℓ,m) mode.
    #
    # Real part
    alms[ℓ+1,m+1] = 1
    ref = synthesize_ecp(alms, nθ, nϕ)
    # Complex part -- Use (-im) not (+im) since the goal is to swap the internal
    #                 (a + ib) to (b + ia), which requires multiplication by -im.
    alms[ℓ+1,m+1] *= -im
    ref += synthesize_ecp(alms, nθ, nϕ) .* im
    # account for doubling due to assumption of real-only symmetry in synthesize_reference
    ref ./= m == 0 ? 1 : 2
    return ref
end

function analyze_ecp_complex_iter(map, lmax, mmax)
    nθ, nϕ = size(map)
    tol = sqrt(eps(maximum(abs, map)))
    function iter(ecp)
        local alms
        alms = analyze_ecp(ecp, lmax, mmax)
        for ii in 2:10  # allow up to 10 iterations before bailing
            ecp′ = synthesize_ecp(alms, nθ, nϕ)
            δecp = ecp - ecp′
            # check for convergence
            if maximum(abs, δecp) < tol
                break
            end
            alms .+= analyze_ecp(δecp, lmax, mmax)
        end
        return alms
    end
    # analyze each of real and imaginary parts separately since the routines require
    # real inputs.
    alms = iter(real.(map)) .+ im .* iter(imag.(map))
    return zerotol!(alms)
end

@testset "Analytic checks — synthesis (Fast ECP)" begin
    # Validates the optimizations:
    #   * Iso-latitude sharing Legendre polynomials
    #   * Ring-pair synthesis via polar symmetry
    #   * FFT-based ϕ synthesis
    @test Y00.(θ, ϕ) ≈ synthesize_ecp_complex(0, 0, n, 2n)
    @test Y10.(θ, ϕ) ≈ synthesize_ecp_complex(1, 0, n, 2n)
    @test Y11.(θ, ϕ) ≈ synthesize_ecp_complex(1, 1, n, 2n)
    @test Y20.(θ, ϕ) ≈ synthesize_ecp_complex(2, 0, n, 2n)
    @test Y21.(θ, ϕ) ≈ synthesize_ecp_complex(2, 1, n, 2n)
    @test Y22.(θ, ϕ) ≈ synthesize_ecp_complex(2, 2, n, 2n)
end

@testset "Analytic checks — analysis (Fast ECP)" begin
    # Validates the optimizations:
    #   * Iso-latitude sharing Legendre polynomials
    #   * Ring-pair synthesis via polar symmetry
    #   * FFT-based ϕ synthesis
    lmax, mmax = 10, 3  # beyond Y22 just as a mild check of any obvious problems
    analyze(Y) = analyze_ecp_complex_iter(Y.(θ, ϕ), lmax, mmax)
    expect(ℓ, m) = setindex!(zeros(complex(eltype(θ)), lmax+1, mmax+1), 1.0, ℓ+1, m+1)
    @test analyze(Y00) ≈ expect(0, 0)
    @test analyze(Y10) ≈ expect(1, 0)
    @test analyze(Y11) ≈ expect(1, 1)
    @test analyze(Y20) ≈ expect(2, 0)
    @test analyze(Y21) ≈ expect(2, 1)
    @test analyze(Y22) ≈ expect(2, 2)
end


# two sets of alms: one that is bandlimited below map Nyquist, and one that extends well
# above the Nyquist
lmax_lo = n ÷ 2 - 1
lmax_hi = 2n
alms_lo = gen_alms(Float64, lmax_lo, seed = 5)
alms_hi = gen_alms(Float64, lmax_hi, seed = 5)

@testset "Aliased ring synthesis" begin
    # Checks that the fast FFT algorithm correctly includes the effects of aliasing
    # high-m modes back down to low-m modes
    ref = synthesize_reference(alms_hi, θ, ϕ)
    ecp = synthesize_ecp(alms_hi, n, 2n)
    @test ref ≈ ecp
    # repeat again with one longer index to make sure even and odd cases both handled
    # correctly.
    ref = synthesize_reference(alms_hi, θ′, ϕ′)
    ecp = synthesize_ecp(alms_hi, n+1, 2n+1)
    @test ref ≈ ecp
    # nϕ == 1 is a special case in aliasing that must be handled appropriately
    alm00 = ComplexF64[1 0; 0 0]
    ref = synthesize_reference(alm00, θ[:,1:1], fill(π/1, n, 1))
    ecp = synthesize_ecp(alm00, n, 1)
    @test ref == ecp == fill(1/sqrt(4π), n, 1)
    # make sure nϕ == 2 still works
    alm11 = ComplexF64[0 0; 0 1im]
    ref = synthesize_reference(alm11, θ[:,1:2], repeat([π/2 3π/2], n, 1))
    ecp = synthesize_ecp(alm11, n, 2)
    @test ref ≈ ecp
end

@testset "Equality of ECP and per-ring synthesis" begin
    # Test that the per-ring synthesis constructs equal answers as full ECP grid
    ecp = synthesize_ecp(alms_hi, n, 2n)
    ecp_r = fill(NaN, n, 2n)
    ecp_p = fill(NaN, n, 2n)
    for j in axes(θ, 1)
        ecp_r[j,:] = synthesize_ring(alms_hi, θ[j,1], ϕ[j,1], 2n, Val(false))
        if j ≤ (n + 1) ÷ 2
            j′ = n - j + 1
            ecp_pairs  = synthesize_ring(alms_hi, θ[j,1], ϕ[j,1], 2n, Val(true))
            ecp_p[j, :] = ecp_pairs[1]
            ecp_p[j′,:] = ecp_pairs[2]
        end
    end
    @test ecp ≈ ecp_r
    @test ecp ≈ ecp_p
end
