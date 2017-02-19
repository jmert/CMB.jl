module PixelCovariance
    using CMB.PixelCovariance
    using CMB.Harmonics
    using Base.Test

    doc"""
    The special-case `Pl2` function should give the same values as computing
    the associated Legendre functions $P_ℓ^2$ in the table method.
    """
    function Plm_l2_equals_Pl2()
        const LMAX = 10
        const x = 2.*rand(10) .- 1
        pl  = PixelCovariance.Pl2(LMAX, x)
        pp  = Harmonics.PlmPlan(LMAX)
        plm = Harmonics.Plm(pp, x)
        # Must select m >= ℓ == 2 elements to avoid undefined behavior
        @test all(pl[:,3:end] .≈ plm[:,3,3:end])
    end

    function runtests()
        # AUX FUNCS FOR COMPUTING POLARIZATION COVARIANCE SPECTRA WEIGHTS
        Plm_l2_equals_Pl2()
    end
end
