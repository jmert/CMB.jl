using Random, SparseArrays

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
    i = 7
    A = sprand(T, m, n, 0.5)
    b = rand(T, n)
    Bc = sparse(collect(1:n), fill(i,n), b, n, n)
    Br = sparse(fill(i,n), collect(1:n), b, n, n)

    @test A * Bc * A' == quadprod(A, b, i, :col)
    @test @inferred(quadprod(A, b, i, :col)) isa SparseMatrixCSC
    @test A * Br * A' ≈  quadprod(A, b, i, :row)
    @test @inferred(quadprod(A, b, i, :row)) isa SparseMatrixCSC
end
