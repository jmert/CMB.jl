using CMB.SphericalHarmonics
using CMB.SphericalHarmonics: centered_range,
    synthesize_reference, synthesize_ecp, synthesize_ring

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

# two sets of alms: one that is bandlimited below map Nyquist, and one that extends well
# above the Nyquist
lmax_lo = n ÷ 2 - 1
lmax_hi = 2n
alms_lo = zeros(ComplexF64, lmax_lo + 1, lmax_lo + 1)
alms_hi = zeros(ComplexF64, lmax_hi + 1, lmax_hi + 1)
for l in 1:lmax_hi
    alms_hi[l+1,1] = randn()
    alms_hi[l+1,2:l+1] .= randn.(ComplexF64)
    if l ≤ lmax_lo
        alms_lo[l+1,1:l+1] .= alms_hi[l+1,1:l+1]
    end
end

# Generate complex fields from running real-only analysis on real and imaginary
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
end

@testset "Round trip synthesis/analysis (reference)" begin
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
