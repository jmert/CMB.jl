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
