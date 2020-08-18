using CMB.PixelCovariance
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

    @testset "Limits at x = ±1" begin
        x₁ = [-1.0, -1.0 + 2sqrt(eps(1.0))]
        x₂ = [ 1.0,  1.0 - 2sqrt(eps(1.0))]
        Fx = zeros(2, lmax+1, 4)

        # Only a rough check of 1 part in 10000 agreement is sufficient to verify that
        # a factor-of-2 bug no longer exists.
        F!(Fx, lmax, x₁)
        @test mapreduce(<(1e-4) ∘ abs, &, diff(Fx[:,2:4,:], dims=1))
        F!(Fx, lmax, x₂)
        @test mapreduce(<(1e-4) ∘ abs, &, diff(Fx[:,2:4,:], dims=1))
    end

    @testset "Preallocated work space" begin
        using .PixelCovariance: unsafe_Fweights!
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
    using .PixelCovariance.CovarianceFields
    using .PixelCovariance.PolarizationConventions

    nside = 4
    lmax = 3nside - 1
    npix = nside2npix(nside)
    pix = pix2vec.(nside, 0:npix-1)
    pixind = length(pix) ÷ 2
    cov = zeros(length(pix), 9)
    Cl = ones(lmax+1, 6)
    ClTT = vcat(zeros(lmax, 6), [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]');
    ClEE = vcat(zeros(lmax, 6), [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]');
    ClBB = vcat(zeros(lmax, 6), [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]');
    ClTE = vcat(zeros(lmax, 6), [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]');
    ClEB = vcat(zeros(lmax, 6), [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]');
    fields = TT | TPol | Pol

    # Numerical calculation also means the output matrix will typically not be
    # perfectly symmetric, which the `isposdef()` check will implicitly require by
    # first checking `ishermitian()`. Here we construct an asymmetry test, and if
    # that is small enough, the test can proceed with an assumed symmetric matrix.
    function asymmetry(A::AbstractMatrix)
        ir, ic = axes(A)
        ir == ic || error("Not square")
        maxvalue = zero(eltype(A))
        maxerror = zero(eltype(A))
        @inbounds for jj in ic
            for ii in ir[jj:end]
                maxvalue = max(maxvalue, A[ii,jj])
                maxerror = max(maxerror, abs(A[ii,jj] - A[jj,ii]))
            end
        end
        return maxerror, maxvalue
    end
    function isnearlysymmetric(A; ulps = 0.1)
        maxerr, maxval = asymmetry(A)
        return maxerr / eps(maxval) ≤ ulps
    end

    function covblocks(cov)
        axes(cov, 1) == axes(cov, 2) || throw(DimensionMismatch())
        np, r = fldmod(size(cov, 1), 3)
        r == 0 || throw(DimensionMismatch())
        iT, iQ, iU = Iterators.partition(axes(cov, 1), np)
        return (;
                TT = cov[iT, iT], TQ = cov[iT, iQ], TU = cov[iT, iU],
                QT = cov[iQ, iT], QQ = cov[iQ, iQ], QU = cov[iQ, iU],
                UT = cov[iU, iT], UQ = cov[iU, iQ], UU = cov[iU, iU]
               )
    end

    @testset "Domain Checking" begin
        # Bad pixind shape
        @test_throws DimensionMismatch pixelcovariance!(cov, pix, ones(2, 2), Cl, fields)

        # Incorrectly-sized output arrays
        @test_throws DimensionMismatch pixelcovariance!(zeros(   0, 9), pix, 1, Cl, fields)
        @test_throws DimensionMismatch pixelcovariance!(zeros(npix, 8), pix, 1, Cl, fields)
        @test_throws DimensionMismatch pixelcovariance!(zeros(   0, 9), pix, 1:npix, Cl, fields)
        @test_throws DimensionMismatch pixelcovariance!(zeros(npix, 9), pix, 1:npix, Cl, fields)

        @test_throws DimensionMismatch pixelcovariance!(zeros(npix, npix,   9), pix, 1, Cl, fields)
        @test_throws DimensionMismatch pixelcovariance!(zeros(npix, npix-1, 9), pix, 1:npix, Cl, fields)
        @test_throws DimensionMismatch pixelcovariance!(zeros(npix, npix,   8), pix, 1:npix, Cl, fields)

        @test_throws DimensionMismatch pixelcovariance!(zeros(npix, npix, 9, 1), pix, 1:npix, Cl, fields)

        # Incorrectly-sized spectra arrays
        @test_throws DimensionMismatch pixelcovariance!(cov, pix, 1, zeros(lmax, 5), fields)

        # Pixel indexing is inbounds
        @test_throws BoundsError pixelcovariance!(cov, pix, -1, Cl, fields)
        @test_throws BoundsError pixelcovariance!(cov, pix, npix+1, Cl, fields)
    end

    @testset "Field mask to row/col indexing" begin
        using .PixelCovariance.CovarianceFields
        using .PixelCovariance.CovarianceFields:
                Field, minrow, mincol, maxrow, maxcol, field_offsets

        # There aren't that many combinations, so just exhaustively check them all.
        # "Abuse" BitMatrix to interpret the bitfield as a 3x3 matrix, and use reductions
        # to find the first selected row/column.
        b = BitMatrix(undef, 3, 3)
        alltrue = true
        for ff in 1:(2^9-1)
            b.chunks[1] = ff
            alltrue &= minrow(Field(ff)) == findfirst(==(true), dropdims(reduce(|, b, dims = 2), dims = 2))
            alltrue &= mincol(Field(ff)) == findfirst(==(true), dropdims(reduce(|, b, dims = 1), dims = 1))
            alltrue &= maxrow(Field(ff)) == findlast(==(true), dropdims(reduce(|, b, dims = 2), dims = 2))
            alltrue &= maxcol(Field(ff)) == findlast(==(true), dropdims(reduce(|, b, dims = 1), dims = 1))
        end
        @test alltrue

        # For performance, we must know that this is properly inferred and optimized by
        # the compiler. Check that return type matches expectation.
        @test @inferred(field_offsets(TT)) isa NTuple{9,Int}
        @test field_offsets(TT)     == (1, 0, 0, 0, 0, 0, 0, 0, 0)
        @test field_offsets(Pol)    == (0, 0, 0, 0, 1, 2, 0, 3, 4)
        @test field_offsets(TT|Pol) == (1, 0, 0, 0, 2, 3, 0, 4, 5)
    end

    @testset "Field calculations" begin
        # Check that specifying field types only writes to some of the output columns
        changedfields(cov) = findall(dropdims(mapreduce(!isnan, &, cov, dims = 1), dims = 1))
        for (field,changes) in zip((TT, TPol, Pol), ([1], [2, 3, 4, 7], [5, 6, 8, 9]))
            fill!(cov, NaN)
            pixelcovariance!(cov, pix, pixind, Cl, field)
            @test changedfields(cov) == changes
        end

        # Similarly, all-zeros spectra should overwrite fields with zero values
        zeroedfields(cov) = findall(dropdims(mapreduce(iszero, &, cov, dims = 1), dims = 1))
        # spectra ordering is TT, EE, BB, TE, TB, EB
        for (spec,filled) in zip((1, 4, 5, 2, 3, 6),
                                 ([1], [2, 3, 4, 7], [2, 3, 4, 7],
                                       [5, 6, 8, 9], [5, 6, 8, 9], [5, 6, 8, 9]))
            fill!(cov, NaN)
            Cl′ = zeros(lmax+1, 6)
            Cl′[:,spec] .= 1.0
            pixelcovariance!(cov, pix, pixind, Cl′, fields)
            @test zeroedfields(cov) == setdiff(1:9, filled)
        end
    end

    @testset "Preallocated work space" begin
        using .PixelCovariance: unsafe_pixelcovariance!
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

    @testset "Positive-definite covariance" begin
        # The resulting covariance matrix should always be postive-definite for any
        # combination of sub-fields (that is square and includes non-zeros on the diagonal,
        # so e.g. TT + TPol - Pol is not a valid combination to check).
        # Test this over combinations of fields and spectra.
        using LinearAlgebra: isposdef, diagind

        # Using the spectra and pixels in the outer test set doesn't actually have a
        # positive-definite matrix (only positive-semidefinite) since there are more
        # pixels (degrees of freedom in pixel space) than there are input modes
        # (degrees of freedom in harmonic space). Use spectra with higher lmax (that will
        # alias power) to fill in more degrees of freedom, the +2 accounting for the fact
        # that l == {0, 1} are unused in polarization (and l==0 uninteresting in T as well).
        Cl′ = ones(npix+2, 6)
        Cl′[1:2, :] .= 0

        # All spectra contribute

        # TT only
        C = pixelcovariance(pix, Cl′, TT)
        C = C[1:npix, 1:npix]
        @test isnearlysymmetric(C)
        @test isposdef(Symmetric(C))

        # Pol only
        C = pixelcovariance(pix, Cl′, Pol)
        C = C[npix+1:end, npix+1:end]
        @test isnearlysymmetric(C)
        @test isposdef(Symmetric(C))

        # TT + TPol + Pol
        C = pixelcovariance(pix, Cl′, TT | TPol | Pol)
        @test isnearlysymmetric(C)
        @test isposdef(Symmetric(C))

        # Drop TE and TB. This only effects the TT + TPol + Pol case.
        Cl′[:,4:5] .= 0

        C = pixelcovariance(pix, Cl′, TT | TPol | Pol)
        @test isnearlysymmetric(C)
        @test isposdef(Symmetric(C))

        # Now also drop EB. This effects both the Pol-only and (TT+TPol+Pol) cases.
        Cl′[:,6] .= 0

        # Pol only
        C = pixelcovariance(pix, Cl′, Pol)
        C = C[npix+1:end, npix+1:end]
        @test isnearlysymmetric(C)
        @test isposdef(Symmetric(C))

        # TT + TPol + Pol
        C = pixelcovariance(pix, Cl′, TT | TPol | Pol)
        @test isnearlysymmetric(C)
        @test isposdef(Symmetric(C))


        # Finally, check that both EE-only and BB-only for Pol-only covariance is still
        # positive definite.
        #
        # First EE:
        Cl′[:,2:3] .= [1.0 0.0]

        C = pixelcovariance(pix, Cl′, Pol)
        C = C[npix+1:end, npix+1:end]
        @test isnearlysymmetric(C)
        @test isposdef(Symmetric(C))

        # Then BB:
        Cl′[:,2:3] .= [0.0 1.0]

        C = pixelcovariance(pix, Cl′, Pol)
        C = C[npix+1:end, npix+1:end]
        @test isnearlysymmetric(C)
        @test isposdef(Symmetric(C))
    end

    @testset "Covariance symmetries" begin
        # The entire 3×3-block covariance matrix is symmetric, so there there are
        # symmetric relationships among each of the blocks.
        covTE = covblocks(pixelcovariance(pix, ClTT + ClTE, TT | TPol))
        covEE = covblocks(pixelcovariance(pix, ClEE, Pol))
        covBB = covblocks(pixelcovariance(pix, ClBB, Pol))
        covEB = covblocks(pixelcovariance(pix, ClEB, Pol))

        # First, the block diagonals are individually self-symmetric.
        @test covTE.TT ≈ covTE.TT'
        @test covEE.QQ ≈ covEE.QQ'
        @test covEE.UU ≈ covEE.UU'
        # The off-diagonal blocks are then transpose-symmetric with each other
        @test covTE.TQ ≈ covTE.QT'
        @test covTE.TU ≈ covTE.UT'
        @test covTE.UQ ≈ covTE.QU'

        # The patterns in QQ and UU swap with each other under a swap of the input spectra.
        @test covEE.QQ == covBB.UU
        @test covEE.UU == covBB.QQ
        # Similarly, the cross QU and UQ fields follow the same patterns in opposite
        # fields when swapping spectra, but there's an additional negation.
        @test covEE.QU == -covBB.UQ
        @test covEE.UQ == -covBB.QU
        # For a covariance with only EB power, QQ and UU are negatives of one another, and
        # the QU and UQ fields are the same
        @test covEB.QQ == -covEB.UU
        @test covEB.UQ ==  covEB.QU
    end

    @testset "Polarization conventions" begin
        covEE_hpix = covblocks(pixelcovariance(pix, ClEE, Pol, HealpixConv))
        covBB_hpix = covblocks(pixelcovariance(pix, ClBB, Pol, HealpixConv))
        covEB_hpix = covblocks(pixelcovariance(pix, ClEB, Pol, HealpixConv))
        covEE_iau  = covblocks(pixelcovariance(pix, ClEE, Pol, IAUConv))
        covBB_iau  = covblocks(pixelcovariance(pix, ClBB, Pol, IAUConv))
        covEB_iau  = covblocks(pixelcovariance(pix, ClEB, Pol, IAUConv))

        # For auto-covariances, pure E or pure B fields are invariant to polarization
        # convention...
        @test covEE_hpix.QQ == covEE_iau.QQ
        @test covEE_hpix.UU == covEE_iau.UU
        @test covBB_hpix.QQ == covBB_iau.QQ
        @test covBB_hpix.UU == covBB_iau.UU
        # ...but EB-correlated spectrum swaps roles.
        @test covEB_hpix.QQ == covEB_iau.UU

        # For cross-covariances, behavior reverses — QU/UQ are the same in EB-correlated
        # spectra...
        @test covEB_hpix.QU == covEB_iau.QU
        @test covEB_hpix.UQ == covEB_iau.UQ
        # ...and are sign-flips in the pure-E/pure-B cases.
        @test covEE_hpix.QU == -covEE_iau.QU
        @test covEE_hpix.UQ == -covEE_iau.UQ
        @test covBB_hpix.QU == -covBB_iau.QU
        @test covBB_hpix.UQ == -covBB_iau.UQ

        # Simple sign flips only work for cases where EB spectrum is zero since the fields
        # are mixed differently between the QQ/UU and QU components in the rotation.
        cov_hpix = covblocks(pixelcovariance(pix, ClEB + ClEE, Pol, HealpixConv))
        cov_iau  = covblocks(pixelcovariance(pix, ClEB + ClEE, Pol, IAUConv))
        @test !(cov_hpix.QQ ≈  cov_iau.QQ)
        @test !(cov_hpix.QQ ≈ -cov_iau.QQ)
        @test !(cov_hpix.UU ≈  cov_iau.QQ)
        @test !(cov_hpix.UU ≈ -cov_iau.QQ)
        @test !(cov_hpix.QU ≈  cov_iau.QU)
        @test !(cov_hpix.QU ≈ -cov_iau.QU)
        @test !(cov_hpix.UQ ≈  cov_iau.UQ)
        @test !(cov_hpix.UQ ≈ -cov_iau.UQ)

        # The mixed case can be reproduced as a particular linear combination of the
        # single-spectra cases, though:
        #
        # First check that same polarization convention is directly sum as expected.
        @test cov_hpix.QQ ≈ covEE_hpix.QQ + covEB_hpix.QQ
        @test cov_hpix.UU ≈ covEE_hpix.UU + covEB_hpix.UU
        @test cov_hpix.QU ≈ covEE_hpix.QU + covEB_hpix.QU
        @test cov_hpix.UQ ≈ covEE_hpix.UQ + covEB_hpix.UQ
        # Then form linear combinations of IAUConv to produce HealpixConv
        @test cov_hpix.QQ ≈  covEE_iau.QQ - covEB_iau.QQ
        @test cov_hpix.UU ≈  covEE_iau.UU - covEB_iau.UU
        @test cov_hpix.QU ≈ -covEE_iau.QU + covEB_iau.QU
        @test cov_hpix.UQ ≈ -covEE_iau.UQ + covEB_iau.UQ
    end

    @testset "Offset axes" begin
        using OffsetArrays

        ax = 0:npix-1
        pix′ = reshape(pix, ax)
        off  = first(axes(pix′, 1)) - first(axes(pix, 1))
        cov′ = reshape(copy(cov), ax, 9)

        pixelcovariance!(cov, pix, pixind, Cl, fields)
        pixelcovariance!(cov′, pix′, pixind+off, Cl, fields)
        @test cov == parent(cov′)

        # First dimension axis of cov and pix must agree
        @test_throws DimensionMismatch pixelcovariance!(cov′, pix, pixind, Cl, fields)
        @test_throws DimensionMismatch pixelcovariance!(cov, pix′, pixind, Cl, fields)

        # Trailing dimension of cov must have 1-based indexing
        @test_throws DimensionMismatch pixelcovariance!(
                reshape(cov, ax, 0:8), pix′, pixind, Cl, fields)
        # Similarly, the spectra must be 1-indexed
        @test_throws ArgumentError pixelcovariance!(
                cov′, pix′, pixind, reshape(Cl, 0:lmax, 6), fields)
        @test_throws ArgumentError pixelcovariance!(
                cov′, pix′, pixind, reshape(Cl, lmax+1, 2:7), fields)
    end
end
