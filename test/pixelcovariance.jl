using LinearAlgebra, Random

@testset "Polarization weights" begin
    nz = 10
    z = collect(range(-1.0, 1.0, length=nz))
    lmax = 10

    F1 = Array{Float64}(undef, lmax + 1, 4)
    FN = Array{Float64}(undef, nz, lmax + 1, 4)

    @testset "Domain Checking" begin
    end

    @testset "Scalar/Vector equality" begin
        fill!(FN, NaN) # Initialize to NaN to check that all entries are set properly
        Fweights!(LegendreUnitNorm(), FN, lmax, z)
        @test all(!isnan, FN)

        for (i, x) in enumerate(z)
            fill!(F1, NaN)
            Fweights!(LegendreUnitNorm(), F1, lmax, x)
            @test all(!isnan, F1)
            @test F1 == FN[i,:,:]
        end
    end

    @testset "Even/Odd symmetry" begin
        # Check that the input is properly symmetric
        @test z == .-(reverse(z))

        Fweights!(LegendreUnitNorm(), FN, lmax, z)
        leven = (0:2:lmax) .+ 1
        lodd  = (1:2:lmax) .+ 1

        for ii in 1:(nz÷2)
            jj = nz - ii + 1
            # F^00 is {even,odd} in z for {even,odd} ℓ
            @test all(FN[ii,leven,1] .== FN[jj,leven,1])
            @test all(FN[ii,lodd,1] .== .-FN[jj,lodd,1])
            # F^10 is {even,odd} in z for {even,odd} ℓ
            @test all(FN[ii,leven,2] .== FN[jj,leven,2])
            @test all(FN[ii,lodd,2] .== .-FN[jj,lodd,2])
            # F^12 is {even,odd} in z for {even,odd} ℓ
            @test all(FN[ii,leven,3] .== FN[jj,leven,3])
            @test all(FN[ii,lodd,3] .== .-FN[jj,lodd,3])
            # F^22 is {odd,even} in z for {even,odd} ℓ
            @test all(FN[ii,leven,4] .== .-FN[jj,leven,4])
            @test all(FN[ii,lodd,4] .== FN[jj,lodd,4])
        end
    end
end
