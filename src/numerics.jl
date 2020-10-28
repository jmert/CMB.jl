@static if VERSION < v"1.6.0-DEV.292"
    using Compat # Compat@v3.23 for sincospi()
end
using SparseArrays
using SparseArrays: AbstractSparseMatrix

# COV_EXCL_START

# TODO: When JuliaLang/julia#35792 is merged, replace with explicit version bounds
@static if isdefined(Base, :cispi)
    using Base: cispi
else
    cispi(z::Real) = complex(reverse(sincospi(z))...)
    cispi(z::Complex) = ((s, c) = sincospi(z); complex(real(c) - imag(s), imag(c) + real(s)))
end

# COV_EXCL_STOP

# Version of sqrt() which skips the domain (x < 0) check for the IEEE floating point types.
# For nonstandard number types, just fall back to a regular sqrt() since eliminating the
# domain check is probably no longer the dominant contributor to not vectorizing.
unchecked_sqrt(x::T) where {T <: Base.IEEEFloat} = Base.sqrt_llvm(x)
unchecked_sqrt(x::T) where {T <: Integer} = unchecked_sqrt(float(x))
unchecked_sqrt(x) = Base.sqrt(x)

"""
    quadprod(A::AbstractSparseMatrixCSC, b::AbstractVecOrMat, n::Integer)

Computes the quadratic product ``ABA^\\top`` efficiently for the case where ``B`` is all zero
except for a small number of columns `b` starting at the `n`th.
"""
function quadprod(A::AbstractSparseMatrix, b::AbstractVecOrMat, n::Integer)
    size(b, 1) == size(A, 2) || throw(DimensionMismatch())

    # sparse * dense naturally returns dense, but we want to dispatch to
    # a sparse-sparse matrix multiplication, so forceably sparsify.
    # - Tests with a few example matrices A show that `sparse(A * b)` is faster than
    #   `A * sparse(b)`.
    w = sparse(A * b)
    p = n + size(b, 2) - 1

    if ndims(w) == 1
        # vector outer product using column view into matrix is fast
        C = w * transpose(view(A, :, n))
    else
        # views are not fast for multiple columns; subset copies are faster
        C = w * transpose(A[:, n:p])
    end
    return C
end
