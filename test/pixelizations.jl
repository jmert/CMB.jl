@static if VERSION < v"1.6.0-DEV.1083"
    import Compat # requires v3.18 for reinterpret(reshape, ...)
end
import .Sphere: colataz, cartvec
import .Pixelizations: Vec3

@testset "unit vectors" begin
    # unit vectors on sphere, naturally created as Vector{Vec3}
    rvecs = pix2vec.(4, 0:191)
    @test rvecs isa Vector{Vec3{Float64}}
    # copy of unit vectors as 3×N matrix
    rvecs_mat = reduce(hcat, rvecs)
    @test rvecs_mat isa Matrix{Float64}

    # pass-through of VecVec3
    pix = @inferred(parse_pixelization(rvecs))
    @test pix isa ArbitraryPixelization
    @test rvecs_mat == export_pixelization(pix)

    # matrices can be reinterpreted as vector of Vec3
    rvecs_reshape = @inferred(parse_pixelization(rvecs_mat))
    @test pix isa ArbitraryPixelization
    @test rvecs_mat == export_pixelization(pix)

    wrong_shape = vcat(rvecs_mat, rvecs_mat[1:1,:])
    @test_throws DimensionMismatch("Cannot interpret 4×192 matrix as vector of pointing vectors"
                                  ) parse_pixelization(wrong_shape)
    @test_throws DimensionMismatch("Cannot interpret 2×192 matrix as vector of pointing vectors"
                                  ) parse_pixelization(wrong_shape[1:2, :])

    @test @inferred(pointing(pix)) isa Vector{Vec3{Float64}}
    @test pointing(pix) == rvecs
end

@testset "HEALPix" begin
    nside = 4
    desc = Dict("type" => "healpix",
                "nside" => nside)
    pix = parse_pixelization(desc)
    @test pix isa HealpixPixelization
    @test desc == export_pixelization(pix)

    desc["pixels"] = collect(0:3nside)
    pix = parse_pixelization(desc)
    @test pix isa HealpixPixelization
    @test desc == export_pixelization(pix)

    @test @inferred(pointing(pix)) isa Vector{Vec3{Float64}}
    @test pointing(pix) == pix2vec.(nside, 0:3nside)
end

@testset "RA/Dec" begin
    # low-res BICEP grid
    nx, ny, res = 236, 100, 4
    lx, hx, ly, hy = -55.0, 55.0, -70.0, -45.0
    dx = (hx - lx) / (nx ÷ res)
    dy = (hy - ly) / (ny ÷ res)
    ra  = range(lx + dx/2,  hx - dx/2, length=nx÷res)
    dec = range(ly + dy/2,  hy - dy/2, length=ny÷res)
    desc = Dict("type" => "radec",
                "ra" => collect(ra),
                "dec" => collect(dec),
                "order" => "row")

    pix = parse_pixelization(desc)
    @test pix isa RADecPixelization
    @test desc == export_pixelization(pix)

    @test_throws ErrorException(match"RA range must span 360 degrees or less"
            ) RADecPixelization(0.0:1.0:360.0, -1:0.5:1, true)
    @test_throws ErrorException(match"Dec range must be in [-90, 90]"
            ) RADecPixelization(-1:0.5:1, 0.0:10.0:90, true)

    @test @inferred(pointing(pix)) isa Vector{Vec3{Float64}}
    @test pointing(pix) == cartvec.(Sphere.unsafe_colataz.(vec(tuple.(dec', ra))))

    pix′ = RADecPixelization(ra, dec, false)
    @test @inferred(pointing(pix′)) isa Vector{Vec3{Float64}}
    @test pointing(pix′) == cartvec.(Sphere.unsafe_colataz.(vec(tuple.(dec, ra'))))
end

@testset "custom show" begin
    show3(io::IO, x) = show(io, MIME"text/plain"(), x)

    rvecs = pix2vec.(4, 0:3)
    pix = ArbitraryPixelization(rvecs)
    msg = sprint(show3, pix)
    @test startswith(msg, "ArbitraryPixelization with $(length(rvecs)) pixels:")
    @test contains(msg, sprint(Base.print_array, rvecs))

    pix = HealpixPixelization(4, 0:3)
    msg = sprint(show3, pix)
    @test startswith(msg, "Nside = 4 HealpixPixelization with 4 pixels:")
    @test contains(msg, sprint(Base.print_array, 0:3))

    pix = RADecPixelization(-0.125:0.25:0.125, -0.125:0.25:0.125, true)
    msg = sprint(show3, pix)
    @test startswith(msg, "RADecPixelization with 2 × 2 = 4 pixels:")
    @test contains(msg, "RA ∈ [-0.25, 0.25] (step 0.25)")
    @test contains(msg, "Dec ∈ [-0.25, 0.25] (step 0.25)")
end
