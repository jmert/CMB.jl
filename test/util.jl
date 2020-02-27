module Util
    using Test
    using Random, SparseArrays
    using CMB.Util
    import ..NumTypes

    (m, n) = (10, 25)
    i = 7
    rng = MersenneTwister(1001)
    A = sprand(rng, m, n, 0.5)
    b = rand(rng, n)
    Bc = sparse(collect(1:n), fill(i,n), b, n, n)
    Br = sparse(fill(i,n), collect(1:n), b, n, n)

    @testset "Quadratic product ($T)" for T in NumTypes
        At = convert(SparseMatrixCSC{T}, A)
        bt = convert(Vector{T}, b)
        Bct = convert(SparseMatrixCSC{T}, Bc)
        Brt = convert(SparseMatrixCSC{T}, Br)
        @test At * Bct * At' == @inferred quadprod(At, bt, i, :col)
        @test At * Brt * At' ≈  @inferred quadprod(At, bt, i, :row)
    end
end
