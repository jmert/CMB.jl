@static if VERSION < v"1.6.0-DEV.1083"
    using Compat # requires v3.18 for reinterpret(reshape, ...)
end
import .Sphere: colataz, cartvec
import .Pixelizations: Vec3, VecVec3, Pixelization, pixelization

@testset "pixelization types" begin
    # constructors
    @test Pixelization("healpix") isa Pixelization{:healpix}
    @test Pixelization(:healpix) isa Pixelization{:healpix}
    @test Pixelization{:healpix}() isa Pixelization{:healpix}

    # getting the specific pixelization name
    @test pixelization(Pixelization(:healpix)) === :healpix

    # pretty-printing of generic pixelizations
    @test sprint(show, Pixelization(:healpix)) == "healpix pixelization"
end

@testset "parsing pixelization descriptions" begin
    @testset "unit vectors ($(float(I)))" for I in (Int, BigInt)
        T = float(I)
        # unit vectors on sphere, naturally created as VecVec3
        rvecs = pix2vec.(I(4), 0:191)
        @test rvecs isa VecVec3{T}
        # copy of unit vectors as 3×N matrix
        rvecs_mat = reduce(hcat, rvecs)
        @test rvecs_mat isa Matrix{T}

        # pass-through of VecVec3
        @test @inferred(parse(Pixelization, rvecs)) === rvecs

        # matrices can be reinterpreted as vector of Vec3
        rvecs_reshape = @inferred(parse(Pixelization, rvecs_mat))
        @test rvecs_reshape == rvecs
        @test eltype(rvecs_reshape) <: Vec3{T}
        wrong_shape = vcat(rvecs_mat, rvecs_mat[1:1,:])
        @test_throws DimensionMismatch("Cannot interpret 4×192 matrix as vector of pointing vectors"
                                      ) parse(Pixelization, wrong_shape)
        @test_throws DimensionMismatch("Cannot interpret 2×192 matrix as vector of pointing vectors"
                                      ) parse(Pixelization, wrong_shape[1:2, :])
    end

    @testset "HEALPix" begin
        nside = 4
        desc = Dict("type" => "healpix",
                    "nside" => nside)
        rvecs = parse(Pixelization, desc)
        @test rvecs isa VecVec3{Float64}
        @test rvecs == pix2vec.(nside, 0:nside2npix(nside)-1)

        desc["pixels"] = 0:3nside
        rvecs = parse(Pixelization, desc)
        @test rvecs isa VecVec3{Float64}
        @test rvecs == pix2vec.(nside, 0:3nside)
    end

    @testset "RA/Dec" begin
        # low-res BICEP grid
        nx, ny, res = 236, 100, 4
        lx, hx, ly, hy = -55.0, 55.0, -70.0, -45.0
        dx = (hx - lx) / (nx ÷ res)
        dy = (hy - ly) / (ny ÷ res)
        ra  = range(lx + dx/2,  hx - dx/2, length=nx÷res)
        dec = range(ly + dy/2,  hy - dy/2, length=ny÷res)
        # BICEP's canonical unraveling is as iso-latitude rings
        radec = copy(reinterpret(reshape, Float64, vec(tuple.(ra, dec'))))
        desc = Dict("type" => "radec",
                    "pixels" => radec)
        rvecs = vec(cartvec.(colataz.(dec', ra)))

        @test radec isa Matrix{Float64}
        @test parse(Pixelization, desc) == rvecs
    end
end
