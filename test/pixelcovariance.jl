using LinearAlgebra, Random

# Neither isequal nor == does what exactly we want --- `isequal(NaN, NaN) == true` but
# `isequal(0.0, -0.0) == false` --- so define custom infix operator (\eqqsim).
⩳(x, y) = (isnan(x) && isnan(x)) || x == y
⩳(x::AbstractArray, y::AbstractArray) = mapreduce(⩳, &, x, y)

# Number of allocating work buffers in internal Work structures.
const nbufs_LegendreWork = 8 # z, y¹, y², pmm2, pm, plm2, plm1, pl
const nbufs_FweightsWork = 2 + nbufs_LegendreWork # y, xy; x aliased to Legendre.Work member

@testset "Polarization weights" begin
    nz = 10
    z = collect(range(-1.0, 1.0, length=nz))
    lmax = 10

    F!(F, lmax, z) = Fweights!(LegendreUnitNorm(), F, lmax, z)
    F1 = Array{Float64}(undef, lmax + 1, 4)
    FN = Array{Float64}(undef, nz, lmax + 1, 4)

    @testset "Domain Checking" begin
        # Invalid lmax
        @test_throws DomainError F!(FN, -1, z)
        # Incorrectly-sized output arrays
        @test_throws DimensionMismatch F!(zeros(nz-1, lmax+1, 4), lmax, z)
        @test_throws DimensionMismatch F!(zeros(nz,   lmax,   4), lmax, z)
        @test_throws DimensionMismatch F!(zeros(nz,   lmax+1, 3), lmax, z)
    end

    @testset "Scalar/Vector equality" begin
        fill!(FN, NaN) # Initialize to NaN to check that all entries are set properly
        Fweights!(LegendreUnitNorm(), FN, lmax, z)

        for (i, x) in enumerate(z)
            fill!(F1, NaN)
            Fweights!(LegendreUnitNorm(), F1, lmax, x)
            @test all(!isnan, F1)
            @test F1 ⩳ FN[i,:,:]
        end
    end

    @testset "Axes preservation (offset axes $axs)" for axs in ((-5:5,), (-5:5, -5:5))
        using OffsetArrays
        using Base: OneTo

        nel = prod(length.(axs))
        Xr = reshape(collect(range(-1, 1, length=nel)), length.(axs)...)
        Xo = reshape(copy(Xr), axs...)
        Fr = zeros(length.(axs)..., lmax+1, 4)
        Fo = zeros(axs..., lmax+1, 4)

        F!(Fr, lmax, Xr)
        F!(Fo, lmax, Xo)
        @test Fr == parent(Fo) # Equality of values
        @test_throws DimensionMismatch F!(Fr, lmax, Xo) # Mismatched axes
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
            @test FN[ii,leven,1] ⩳ FN[jj,leven,1]
            @test FN[ii,lodd,1] ⩳ -FN[jj,lodd,1]
            # F^10 is {even,odd} in z for {even,odd} ℓ
            @test FN[ii,leven,2] ⩳ FN[jj,leven,2]
            @test FN[ii,lodd,2] ⩳ -FN[jj,lodd,2]
            # F^12 is {even,odd} in z for {even,odd} ℓ
            @test FN[ii,leven,3] ⩳ FN[jj,leven,3]
            @test FN[ii,lodd,3] ⩳ -FN[jj,lodd,3]
            # F^22 is {odd,even} in z for {even,odd} ℓ
            @test FN[ii,leven,4] ⩳ -FN[jj,leven,4]
            @test FN[ii,lodd,4] ⩳ FN[jj,lodd,4]
        end
    end

    @testset "Preallocated work space" begin
        using CMB.PixelCovariance: unsafe_Fweights!
        x = z[2]
        norm = LegendreUnitNorm()
        work1 = PixelCovariance.FweightsWork(norm, F1, x)
        workN = PixelCovariance.FweightsWork(norm, FN, z)
        # Check equality before allocations to ensure the methods have been compiled.
        @test unsafe_Fweights!(norm, copy(F1), lmax, x) == unsafe_Fweights!(work1, F1, lmax, x)
        @test unsafe_Fweights!(norm, copy(FN), lmax, z) == unsafe_Fweights!(workN, FN, lmax, z)
        # Lower bounds on allocations required
        nb1 = sizeof(work1.x) * nbufs_FweightsWork
        nbN = sizeof(workN.x) * nbufs_FweightsWork
        @test nb1 <= @allocated unsafe_Fweights!(norm, F1, lmax, x)
        @test nbN <= @allocated unsafe_Fweights!(norm, FN, lmax, z)
        # Pre-allocated version doesn't allocate
        @test 0 == @allocated unsafe_Fweights!(work1, F1, lmax, x)
        @test 0 == @allocated unsafe_Fweights!(workN, FN, lmax, z)
    end
end

@testset "Pixel-pixel covariance" begin
    nside = 4
    lmax = 3nside - 1
    npix = nside2npix(nside)
    pix = pix2vec.(nside, 0:npix-1)
    pixind = length(pix) ÷ 2
    cov = zeros(length(pix), 9)
    Cl = ones(lmax+1, 6)
    fields = PixelCovariance.TT | PixelCovariance.TPol | PixelCovariance.Pol

    @testset "Domain Checking" begin
        # Incorrectly-sized output arrays
        @test_throws DimensionMismatch pixelcovariance!(zeros(   0, 9), pix, 1, Cl, fields)
        @test_throws DimensionMismatch pixelcovariance!(zeros(npix, 8), pix, 1, Cl, fields)
        # Incorrectly-sized spectra arrays
        @test_throws DimensionMismatch pixelcovariance!(cov, pix, 1, zeros(lmax, 5), fields)
        # Pixel indexing is inbounds
        @test_throws BoundsError pixelcovariance!(cov, pix, -1, Cl, fields)
        @test_throws BoundsError pixelcovariance!(cov, pix, npix+1, Cl, fields)
    end

    @testset "Preallocated work space" begin
        using CMB.PixelCovariance: unsafe_pixelcovariance!
        norm = LegendreUnitNorm()
        work = PixelCovariance.PixelCovWork{Float64}(lmax, norm)
        # Check equality before allocations to ensure the methods have been compiled.
        @test unsafe_pixelcovariance!(norm, copy(cov), pix, pixind, Cl, fields) ==
                unsafe_pixelcovariance!(work, cov, pix, pixind, Cl, fields)

        # Lower bound on allocations required
        nb  = sizeof(eltype(first(pix))) * nbufs_FweightsWork
        nb += sizeof(eltype(first(pix))) * 4 * (lmax + 1)
        @test nb <= @allocated unsafe_pixelcovariance!(norm, cov, pix, pixind, Cl, fields)

        # Currently broken, even on v1.5 --- haven't figured out where the allocations are
        # happening
        @test_broken nb > @allocated unsafe_pixelcovariance!(work, cov, pix, pixind, Cl, fields)
    end
end
