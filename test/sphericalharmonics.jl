using Test
using CMB.SphericalHarmonics
using CMB.SphericalHarmonics: synthesize_reference, analyze_reference, centered_range
import ..NumTypes

# Define the analytic expressions for the first few spherical harmonics; verify
# implementations against these.
Y00(θ, ϕ) =  sqrt(1 / 4oftype(θ, π)) * complex(true)
Y10(θ, ϕ) =  sqrt(3 / 4oftype(θ, π)) * cos(θ) * complex(true)
Y11(θ, ϕ) = -sqrt(3 / 8oftype(θ, π)) * sin(θ) * cis(ϕ)
Y20(θ, ϕ) =  sqrt(5 / 16oftype(θ, π)) * (3cos(θ)^2 - 1) * complex(true)
Y21(θ, ϕ) = -sqrt(15 / 8oftype(θ, π)) * sin(θ) * cos(θ) * cis(ϕ)
Y22(θ, ϕ) =  sqrt(15 / 32oftype(θ, π)) * sin(θ)^2 * cis(2ϕ)

@testset "Analytic checks — synthesis (reference)" begin
    n = 50
    θ = repeat(centered_range(0.0, 1.0π, n), 1, 2n)
    ϕ = repeat(centered_range(0.0, 2.0π, 2n)', n, 1)

    # Generate complex fields from running real-only analysis on real and imaginary
    # components separately.
    function synth_ref_complex(ℓ, m, θ, ϕ, T::Type = Float64)
        alms = zeros(Complex{T}, ℓ + 1, m + 1)
        # Assign unit power to single delta (ℓ,m) mode.
        #
        # Real part
        alms[ℓ+1,m+1] = 1
        ref = complex(synthesize_reference(alms, θ, ϕ))
        # Complex part -- Use (-im) not (+im) since the goal is to swap the internal
        #                 (a + ib) to (b + ia), which requires multiplication by -im.
        alms[ℓ+1,m+1] *= -im
        ref += synthesize_reference(alms, θ, ϕ) .* im
        # account for doubling due to assumption of real-only symmetry in synthesize_reference
        ref ./= m == 0 ? 1 : 2
        return ref
    end

    @test Y00.(θ, ϕ) ≈ synth_ref_complex(0, 0, θ, ϕ)
    @test Y10.(θ, ϕ) ≈ synth_ref_complex(1, 0, θ, ϕ)
    @test Y11.(θ, ϕ) ≈ synth_ref_complex(1, 1, θ, ϕ)
    @test Y20.(θ, ϕ) ≈ synth_ref_complex(2, 0, θ, ϕ)
    @test Y21.(θ, ϕ) ≈ synth_ref_complex(2, 1, θ, ϕ)
    @test Y22.(θ, ϕ) ≈ synth_ref_complex(2, 2, θ, ϕ)
end

@testset "Analytic checks — analysis (reference)" begin
end

@testset "Round trip synthesis/analysis (reference)" begin
end
