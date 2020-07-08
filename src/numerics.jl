using SparseArrays

"""
    quadprod(A, b, n, dir=:col)

Computes the quadratic product ``ABA^\\top`` efficiently for the case where ``B`` is all zero
except for the `n`th column or row vector `b`, for `dir = :col` or `dir = :row`,
respectively.
"""
@inline function quadprod(A, b, n, dir::Symbol=:col)
    if dir == :col
        return (A * sparse(b)) * view(A, :, n)'
    elseif dir == :row
        return view(A, :, n) * (A * sparse(b))'
    else
        error("Unrecognized direction `dir = $(repr(dir))`.")
    end
end
