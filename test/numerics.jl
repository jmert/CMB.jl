using Random
using OffsetArrays
using SparseArrays

@testset "Unchecked square root ($T)" for T in [S for S in NumTypes if float(S) <: Base.IEEEFloat]
    using CMB: unchecked_sqrt
    @test isnan(unchecked_sqrt(T(-1)))
    x = rand(T) + 1
    @test unchecked_sqrt(x) == sqrt(x)
end

@testset "Quadratic product ($T)" for T in NumTypes
    using SparseArrays
    using CMB: quadprod
    (m, n) = (10, 25)
    i, j = 7, 3
    A = sprand(T, m, n, 0.5)

    # single vector
    b = rand(T, n)
    B = sparse(collect(1:n), fill(i,n), b, n, n)
    @test A * B * A' == quadprod(A, b, i)
    @test @inferred(quadprod(A, b, i)) isa SparseMatrixCSC{T, Int}

    # block of columns
    b = rand(T, n, j)
    B[:,i:i+j-1] .= b
    @test A * B * A' == quadprod(A, b, i)
    @test @inferred(quadprod(A, b, i)) isa SparseMatrixCSC{T, Int}
end

@testset "Triangular [Un]packing" begin
    import CMB: trillen, trildim, packtril, unpacktril, packtril!, unpacktril!
    import LinearAlgebra: LowerTriangular

    test_sizes = ((i, j) for j in (0x1, 2, 5), i in (1, 2, 0x5) if j <= i)

    @testset "Sizes & Lengths, (n, m) = ($n, $m)" for (n, m) in test_sizes
        @test @inferred(trildim(trillen(n, m), m)) == (n, m)
        n == m && @test @inferred(trildim(trillen(n, m))) == (n, m)

        # only test error conditions when m != 1 since
        #      L + 1 = trilen(n, 1) + 1
        #   == L + 1 = trilen(n + 1, 1)
        if m != 1
            @test_throws ArgumentError trildim(trillen(n, m) + 1, m)
            n == m && @test_throws ArgumentError trildim(trillen(n, m) + 1)
        end
    end

    @testset "Pack ($n, $m)" for (n, m) in test_sizes
        cols = repeat((1:n)', outer = (n, 1)) |> LowerTriangular |> x -> x[:, 1:m]
        rows = repeat( 1:n,   outer = (1, n)) |> LowerTriangular |> x -> x[:, 1:m]
        @test @inferred(packtril(cols)) == cols[cols .!= 0]
        @test packtril(rows) == rows[rows .!= 0]
    end

    @testset "In-place pack ($n, $m)" for (n, m) in test_sizes
        cols = repeat((1:n)', outer = (n, 1)) |> LowerTriangular |> x -> x[:, 1:m]
        colsv = cols[cols .!= 0]
        vec = Vector{Float64}(undef, trillen(n, m))
        @test @inferred(packtril!(vec, cols)) == colsv
        @test vec == colsv

        # vector size must have correct length
        @test_throws DimensionMismatch packtril!(vec[2:end], cols)
        @test_throws DimensionMismatch packtril!([0; vec], cols)

        # vector and matrix can have different axes, just need compatible lengths
        vec′ = OffsetArray(vec, 0:length(vec)-1)
        cols′ = OffsetArray(cols, 0:n-1, 1:m)
        fill!(vec, 0.0)
        @test packtril!(vec′, cols) !== nothing  # no error
        @test vec == colsv
        fill!(vec, 0.0)
        @test packtril!(vec, cols′) !== nothing  # no error
        @test vec == colsv
    end

    @testset "Unpack ($n, $m)" for (n, m) in test_sizes
        cols = repeat((1:n)', outer = (n, 1)) |> LowerTriangular |> x -> x[:, 1:m]
        rows = repeat( 1:n,   outer = (1, n)) |> LowerTriangular |> x -> x[:, 1:m]
        colsv = cols[cols .!= 0]
        rowsv = rows[rows .!= 0]
        @test @inferred(unpacktril(colsv, m)) == cols
        @test unpacktril(rowsv, m) == rows
        if m == n
            @test @inferred(unpacktril(colsv)) == cols
            @test unpacktril(rowsv) == rows
        end
    end

    @testset "In-place unpack ($n, $m)" for (n, m) in test_sizes
        cols = repeat((1:n)', outer = (n, 1)) |> LowerTriangular |> x -> x[:, 1:m]
        colsv = cols[cols .!= 0]
        mat = zeros(n, m)  # need upper triangle to be zero initialized
        @test @inferred(unpacktril!(mat, colsv)) == cols
        @test mat == cols

        # matrix size must have correct length
        @test_throws DimensionMismatch unpacktril!(mat[2:end, 2:end], colsv)
        @test_throws DimensionMismatch unpacktril!([zeros(1, m+1); zeros(n,1) mat], colsv)

        # vector and matrix can have different axes, just need compatible lengths
        colsv′ = OffsetArray(colsv, 0:length(colsv)-1)
        mat′ = OffsetArray(mat, 0:n-1, 1:m)
        fill!(mat, 0.0)
        @test unpacktril!(mat′, colsv) !== nothing  # no error
        @test mat == cols
        fill!(mat, 0.0)
        @test unpacktril!(mat, colsv′) !== nothing  # no error
        @test mat == cols
    end
end
