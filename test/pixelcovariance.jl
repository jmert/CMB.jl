module PixelCovariance
    using CMB.PixelCovariance
    using CMB.Harmonics
    using Base.Test

    # The special-case `Pl2` function should give the same values as computing
    # the associated Legendre functions P_ℓ^2 in the table method.
    @testset "Equality of Pl2 and P_ℓ^2" begin
        const LMAX = 10
        const x = 2.*rand(10) .- 1
        pl  = PixelCovariance.Pl2(LMAX, x)
        pp  = Harmonics.PlmPlan(LMAX)
        plm = Harmonics.Plm(pp, x)
        # Must select m >= ℓ == 2 elements to avoid undefined behavior
        @test all(pl[:,3:end] .≈ plm[:,3,3:end])
    end
end
