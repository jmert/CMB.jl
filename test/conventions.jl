using .StokesCrossFields
using .PolarizationConventions

@testset "Stokes bit mask to row/col indexing" begin
    using .StokesCrossFields: Field, minrow, mincol, maxrow, maxcol,
            field_count, field_offsets

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
    @test @inferred(field_count(TT)) isa Int
    @test field_offsets(TT)     == (1, 0, 0, 0, 0, 0, 0, 0, 0)
    @test field_offsets(Pol)    == (0, 0, 0, 0, 1, 2, 0, 3, 4)
    @test field_offsets(TT|Pol) == (1, 0, 0, 0, 2, 3, 0, 4, 5)
    @test field_count(TT) == 1
    @test field_count(Pol) == 4
    @test field_count(TT|Pol) == 5
end

