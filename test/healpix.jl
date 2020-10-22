# Each of the following arrays can be initialized directly from examining Fig. 4 of
# Górski, et al (2005).
hpix4_pix     = collect(0:191)
hpix4_ring    = Int[]
hpix4_ringidx = Int[]
hpix4_isnorth        = Bool[]
hpix4_issouth        = Bool[]
hpix4_iscap          = Bool[]
hpix4_isnorthcap     = Bool[]
hpix4_issouthcap     = Bool[]
hpix4_isequbelt      = Bool[]
hpix4_isnorthequbelt = Bool[]
hpix4_issouthequbelt = Bool[]

for ii in 1:3
    append!(hpix4_ring, ii .* ones(Int, 4ii))
    append!(hpix4_ringidx, 1:4ii)
    append!(hpix4_isnorth,        fill(true,  4ii))
    append!(hpix4_issouth,        fill(false, 4ii))
    append!(hpix4_iscap,          fill(true,  4ii))
    append!(hpix4_isnorthcap,     fill(true,  4ii))
    append!(hpix4_issouthcap,     fill(false, 4ii))
    append!(hpix4_isequbelt,      fill(false, 4ii))
    append!(hpix4_isnorthequbelt, fill(false, 4ii))
    append!(hpix4_issouthequbelt, fill(false, 4ii))
end
for ii in 1:5
    append!(hpix4_ring, (3+ii) .* ones(Int, 16))
    append!(hpix4_ringidx, 1:16)
    append!(hpix4_isnorth,        fill(true,  16))
    append!(hpix4_issouth,        fill(false, 16))
    append!(hpix4_iscap,          fill(false, 16))
    append!(hpix4_isnorthcap,     fill(false, 16))
    append!(hpix4_issouthcap,     fill(false, 16))
    append!(hpix4_isequbelt,      fill(true,  16))
    append!(hpix4_isnorthequbelt, fill(true,  16))
    append!(hpix4_issouthequbelt, fill(false, 16))
end
for ii in 1:4
    append!(hpix4_ring, (8+ii) .* ones(Int, 16))
    append!(hpix4_ringidx, 1:16)
    append!(hpix4_isnorth,        fill(false, 16))
    append!(hpix4_issouth,        fill(true,  16))
    append!(hpix4_iscap,          fill(false, 16))
    append!(hpix4_isnorthcap,     fill(false, 16))
    append!(hpix4_issouthcap,     fill(false, 16))
    append!(hpix4_isequbelt,      fill(true,  16))
    append!(hpix4_isnorthequbelt, fill(false, 16))
    append!(hpix4_issouthequbelt, fill(true,  16))
end
for ii in 1:3
    jj = 3 - ii + 1
    append!(hpix4_ring, (12+ii) .* ones(Int, 4jj))
    append!(hpix4_ringidx, 1:4jj)
    append!(hpix4_isnorth,        fill(false, 4jj))
    append!(hpix4_issouth,        fill(true,  4jj))
    append!(hpix4_iscap,          fill(true,  4jj))
    append!(hpix4_isnorthcap,     fill(false, 4jj))
    append!(hpix4_issouthcap,     fill(true,  4jj))
    append!(hpix4_isequbelt,      fill(false, 4jj))
    append!(hpix4_isnorthequbelt, fill(false, 4jj))
    append!(hpix4_issouthequbelt, fill(false, 4jj))
end

@testset "Nside properties" begin
    @test nside2npix(4) == length(hpix4_pix)
    @test npix2nside(length(hpix4_pix)) == 4
    @test nside2nring(4) == hpix4_ring[end]
    @test nring2nside(hpix4_ring[end]) == 4
    @test nside2pixarea(4) == 4π / length(hpix4_pix)
end

@testset "Pixel indices" begin
    @test all(pix2ring.(4, hpix4_pix) .== hpix4_ring)
    @test all(pix2ringidx.(4, hpix4_pix) .== hpix4_ringidx)
end

@testset "Pixel classification" begin
    @test all(isnorth.(       4, hpix4_pix) .== hpix4_isnorth)
    @test all(issouth.(       4, hpix4_pix) .== hpix4_issouth)
    @test all(iscap.(         4, hpix4_pix) .== hpix4_iscap)
    @test all(isnorthcap.(    4, hpix4_pix) .== hpix4_isnorthcap)
    @test all(issouthcap.(    4, hpix4_pix) .== hpix4_issouthcap)
    @test all(isequbelt.(     4, hpix4_pix) .== hpix4_isequbelt)
    @test all(isnorthequbelt.(4, hpix4_pix) .== hpix4_isnorthequbelt)
    @test all(issouthequbelt.(4, hpix4_pix) .== hpix4_issouthequbelt)
end

@testset "Validity checks" begin
    # HEALPix pixelization is defined to be valid for pixel indices up to near limits of
    # Int64, but the corresponding Nside fits within Int32. Make sure internal type
    # promotion rules mix these facts appropriately.
    @testset "$T Nside" for T in (Int32, Int64)
        @test all(ishealpixok, T(2) .^ (0:29))
        @test all(!ishealpixok, (T(0), T(2)^30))
        @test_throws InvalidNside checkhealpix(T(0))

        @test all(p -> ishealpixok(T(4), p), hpix4_pix)
        @test all(p -> !ishealpixok(T(4), p), (-1, 192))
        @test_throws InvalidPixel checkhealpix(T(4), 192)
    end
    @test !ishealpixok(BigInt(2)^300)
    @test_throws InvalidNside checkhealpix(BigInt(2)^300)
    @test !ishealpixok(4, BigInt(2)^300)
    @test_throws InvalidPixel checkhealpix(4, BigInt(2)^300)

    for pix2fn in (pix2z, pix2theta, pix2phi, pix2ang, pix2vec)
        @test_throws InvalidNside pix2fn(5,  0)
        @test_throws InvalidPixel pix2fn(4, -1)
    end
    for ang2fn in (ang2pix,)
        @test_throws InvalidNside ang2fn(5, π/2, π/2)
        @test_throws DomainError ang2fn(4, 2π, π/2) # θ > π when required π ∈ [0,π]
    end
    for vec2fn in (vec2pix,)
        @test_throws InvalidNside vec2fn(5, Float64[0, 0, 1])
        @test_throws DimensionMismatch vec2fn(4, Float64[0, 1])
        @test_throws DimensionMismatch vec2fn(4, Float64[0, 0, 0, 1])
    end
end

@testset "Pixel identity" begin
    @test all(pix == ang2pix(4, pix2ang(4, pix)...) for pix in hpix4_pix)
    @test all(pix == vec2pix(4, pix2vec(4, pix)) for pix in hpix4_pix)
end

@testset "Wrap-around pixel ang2pix" begin
    # For pixel with centers at ϕ = 0, make sure the ϕ ≲ 0 (i.e. ϕ ≈ 2π - δϕ) are
    # handled correctly. Only occurs within the equatorial belt.
    nside = 4
    _, ϕ₀ = pix2ang(nside, nside2npixcap(nside))
    ϕ₀ = -ϕ₀ / 2
    for ring in (nside+1):2:(2nside-1)
        pix = nside2npixcap(nside) + 4nside*ring
        θ,_ = pix2ang(nside, pix)
        @test pix == ang2pix(nside, θ, ϕ₀)
        @test pix == ang2pix(nside, θ, 2π + ϕ₀)
    end
end

@testset "Internal type stability" begin
    @testset "nside query functions Nside type $T" for T in (Int32, Int64, BigInt)
        # small types promote up to Int64, otherwise type stable
        S = T <: Base.BitInteger32 ? Int64 : T
        @test npix2nside(T(length(hpix4_pix))) isa S
        @test nring2nside(T(hpix4_ring[end])) isa S
        @test nside2npix(T(4)) isa S
        @test nside2nring(T(4)) isa S
        @test nside2npixcap(T(4)) isa S
        @test nside2npixequ(T(4)) isa S
        @test nside2pixarea(T(4)) isa float(T)
    end
    @testset "pixel query functions for Nside type $T" for T in (Int32, Int64, BigInt)
        S = T <: Base.BitInteger32 ? Int64 : T
        # promoted type
        @test pix2ring(T(4), T(0)) isa S
        @test pix2ringidx(T(4), T(0)) isa S

        Tπ = convert(float(S), π)
        # π conversion done at correct precision
        # N.B. constants are all easily obtained from Fig. 4 of Górski, et al (2005).
        @test nside2pixarea(T(1)) == 4Tπ / 12
        @test pix2theta(T(4), T(88)) == Tπ / 2
        @test pix2phi(T(4), T(44)) == Tπ / 2
        @test pix2z(T(4), T(40)) == convert(float(T), 0.5)
        @test pix2theta(T(4), T(40)) ≈ Tπ/3 rtol=eps(Tπ/3)

        # also check that internal integer calculations do not overflow (by preemptively
        # promoting --- `unsafe_pix2*` variants are not similarly protected).
        # - (2^16)^2 > typemax(Int32)
        @test pix2ring(T(2)^16, T(0)) == pix2ring(S(2)^16, T(0))
        @test pix2ringidx(T(2)^16, T(0)) == pix2ringidx(T(2)^16, T(0))
        @test pix2z(T(2)^16, T(0)) == pix2z(S(2)^16, S(0))
        @test pix2theta(T(2)^16, T(0)) == pix2theta(S(2)^16, S(0))
        @test pix2phi(T(2)^16, T(0)) == pix2phi(S(2)^16, S(0))
        @test pix2ang(T(2)^16, T(0)) == pix2ang(S(2)^16, S(0))
        @test pix2vec(T(2)^16, T(0)) == pix2vec(S(2)^16, S(0))
    end
    @testset "Mixed input type promotions ($fn)" for fn in (pix2ring, pix2ringidx,
                                                pix2z, pix2theta,
                                                pix2phi, pix2ang, pix2vec)
        @test fn(4, BigInt(0)) == fn(BigInt(4), BigInt(0))
        @test fn(BigInt(4), 0) == fn(BigInt(4), BigInt(0))
    end
    @test ang2pix(4, 1f0, 1e0) == ang2pix(4, 1e0, 1e0)
    @test ang2pix(4, 1e0, 1f0) == ang2pix(4, 1e0, 1e0)
end

@testset "Accuracy of pix2vec" begin
    import LinearAlgebra: norm
    function simple_pix2vec(nside, p)
        z = Healpix.unsafe_pix2z(nside, p)
        ϕ = Healpix.unsafe_pix2phi(nside, p)
        sinθ = sqrt(fma(-z, z, one(z)))
        y, x = sinθ .* sincos(ϕ)
        return [x, y, z]
    end
    # First check that the specialized implementation is largely correct.
    @test all(pix2vec.(4, hpix4_pix) .≈ simple_pix2vec.(4, hpix4_pix))

    # Now check that we do get more accurate results than the simple implementation.
    @testset "Nside = $nside" for nside in 2 .^ (10:29) # Nside 1024 through MAX_NSIDE
        npix = nside2npix(nside)

        # use the first two and last two rings near the poles where |z| ≈ 1
        pix = [0:11; npix-12:npix-1]
        ref = [Float64.(v) for v in simple_pix2vec.(big(nside), big.(pix))]
        v0 = simple_pix2vec.(nside, pix)
        v1 = pix2vec.(nside, pix)
        # z coordinates are the same since it is only internal uses which can be optimized
        @test all(getindex.(v1, 3) .== getindex.(ref, 3))
        # the optized implementation should match more closely the vector at extended
        # precision than the simple version does
        @test all(norm.(v1 .- ref) .< norm.(v0 .- ref))
    end

    # the azimuthal angle inversion also has an optimized path, so check the equatorial
    # ring for more accurate alignment as well
    #
    # check only the first quadrant of Nside == 2^16 since that's already *many* pixels
    nside = 2^16
    pix = (0:nside-1) .+ (nside2npixequ(nside) - 4nside)
    ref = [Float64.(v) for v in simple_pix2vec.(big(nside), big.(pix))]
    v0 = simple_pix2vec.(nside, pix)
    v1 = pix2vec.(nside, pix)
    # only test that the large-scale view gets more accurate since small jitters in
    # which way rounding occurs makes testing every individual pixel hard
    @test sum(norm.(v1 .- ref)) < sum(norm.(v0 .- ref))
end
