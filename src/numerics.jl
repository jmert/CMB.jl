using LinearAlgebra: mul!
using SparseArrays
using SparseArrays: AbstractSparseMatrix

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

"""
    L = trillen(n::Int, m::Int = n)

Calculates the length `L` of a vector required to store the elements of an `n` × `m`
lower-triangular matrix.

See also: [`trildim`](@ref)

# Example
```jldoctest; setup = (import CMB: trillen)
julia> trillen(3)
6

julia> trillen(3, 3)
6

julia> trillen(3, 2)
5
```
"""
function trillen(n::Integer, m::Integer = n)
    m <= n || throw(ArgumentError("require n ≤ m; got n = $n, m = $m"))
    k = n - m
    return n * (n + 1) ÷ 2 - k * (k + 1) ÷ 2
end

"""
    n, m = trildim(L::Integer[, m::Integer])

Calculates the dimensions `n` × `m` of a matrix capable of storing the `L` elements of a
lower-triangular matrix. If `m < n`, then `m` must be provided.

See also: [`trillen`](@ref)

# Example
```jldoctest; setup = (import CMB: trildim)
julia> trildim(6)
(3, 3)

julia> trildim(6, 3)
(3, 3)

julia> trildim(5)
ERROR: ArgumentError: length 5 incompatible with inferred triangle dims (2, 2)
[...]

julia> trildim(5, 2)
(3, 2)

julia> trildim(5, 1)
(5, 1)
```
"""
function trildim end

function trildim(L::Integer)
    n′ = (sqrt(8L + 1) - 1) / 2
    n = trunc(Int, n′)
    n == n′ || throw(ArgumentError("length $L incompatible with inferred triangle dims ($n, $n)"))
    return (n, n)
end
function trildim(L::Integer, m::Integer)
    n′ = (2L + m^2 - m) / 2m
    n = trunc(Int, n′)
    n == n′ || throw(ArgumentError("length $L incompatible with inferred triangle dims ($n, $m)"))
    return (n, m)
end

"""
    vec = packtril!(vec::AbstractVector, mat::AbstractMatrix)

Packs the lower triangle of the matrix `mat` into the vector `vec`.

See also: [`packtril`](@ref), [`trillen`](@ref), [`unpacktril!`](@ref)

# Example
```jldoctest; setup = (import CMB: packtril!)
julia> packtril!(zeros(6), [1 0 0;
                            2 3 0;
                            4 5 6])
6-element Vector{Float64}:
 1.0
 2.0
 4.0
 3.0
 5.0
 6.0

julia> packtril!(zeros(5), [1 0;
                            2 3;
                            4 5])
5-element Vector{Float64}:
 1.0
 2.0
 4.0
 3.0
 5.0
```
"""
function packtril!(vec::AbstractVector, mat::AbstractMatrix)
    L, (n, m) = length(vec), size(mat)
    L == trillen(n, m) || throw(DimensionMismatch(
            "cannot pack size ($n, $m) triangular matrix into length $L vector"))
    i₀, j₀, k₀ = firstindex(mat, 1), firstindex(mat, 2), firstindex(vec, 1)
    kk = 0
    @inbounds for jj in 0:m-1, ii in jj:n-1
        vec[k₀+kk] = mat[i₀+ii, j₀+jj]
        kk += 1
    end
    return vec
end

"""
    mat = unpacktril!(mat::AbstractMatrix, vec::AbstractVector)

Unpacks the vector `vec` into the lower triangle of the matrix `mat`.

See also: [`packtril!`](@ref), [`trildim`](@ref), [`unpacktril`](@ref)

# Example
```jldoctest; setup = (import CMB: unpacktril!)
julia> unpacktril!(zeros(3, 3), 1:6)
3×3 Matrix{Float64}:
 1.0  0.0  0.0
 2.0  4.0  0.0
 3.0  5.0  6.0

julia> unpacktril!(zeros(3, 2), 1:5)
3×2 Matrix{Float64}:
 1.0  0.0
 2.0  4.0
 3.0  5.0
```
"""
function unpacktril!(mat::AbstractMatrix, vec::AbstractVector)
    L, (n, m) = length(vec), size(mat)
    L == trillen(n, m) || throw(DimensionMismatch(
            "cannot unpack length $L vector into size ($n, $m) triangular matrix"))
    i₀, j₀, k₀ = firstindex(mat, 1), firstindex(mat, 2), firstindex(vec, 1)
    kk = 0
    @inbounds for jj in 0:m-1, ii in jj:n-1
        mat[i₀+ii, j₀+jj] = vec[k₀+kk]
        kk += 1
    end
    return mat
end

"""
    vec = packtril(mat::AbstractMatrix)

Packs the lower triangle of matrix `mat` into the vector `vec`.

See also: [`packtril!`](@ref), [`unpacktril`](@ref)

# Example
```jldoctest; setup = (import CMB: packtril)
julia> packtril([1 0 0;
                 2 3 0;
                 4 5 6])
6-element Vector{Int64}:
 1
 2
 4
 3
 5
 6
```
"""
function packtril(mat::AbstractMatrix)
    vec = fill!(similar(mat, trillen(size(mat)...)), zero(eltype(mat)))
    return packtril!(vec, mat)
end

"""
    mat = unpacktril(vec::AbstractVector[, m::Int])

Unpacks the vector `vec` into the lower triangle of a matrix `mat`, inferring the dimensions
based on the length `L = length(vec)` (and `m = size(mat, 2)` if necessary) with
[`trildim`](@ref).

See also: [`packtril`](@ref), [`unpacktril!`](@ref)

# Example
```jldoctest; setup = (import CMB: unpacktril)
julia> unpacktril(1:6)
3×3 Matrix{Int64}:
 1  0  0
 2  4  0
 3  5  6

julia> unpacktril(1:5, 2)
3×2 Matrix{Int64}:
 1  0
 2  4
 3  5
```
"""
function unpacktril end

function unpacktril(vec::AbstractVector)
    mat = fill!(similar(vec, trildim(length(vec))...), zero(eltype(vec)))
    return unpacktril!(mat, vec)
end
function unpacktril(vec::AbstractVector, m::Integer)
    mat = fill!(similar(vec, trildim(length(vec), m)...), zero(eltype(vec)))
    return unpacktril!(mat, vec)
end

"""
    x = cg(A, b::AbstractVector{T};
                 maxiter::Integer = length(b) * (1 + (T <: Complex)),
                 abstol::Real = zero(real(T)),
                 reltol::Real = sqrt(eps(real(T)))
          ) where {T}
    x = cg(op, b; kws...)
    (; info...) = cg!(x, A, b; kws...)
    (; info...) = cg!(x, op!, b; kws...)

Minimize the residual in the system of equations ``\\lVert Ax - b \\rVert`` for the
unknown vector ``x`` given matrix ``A`` and vector ``b`` using the conjugate gradient
method.
Optionally, the system ``A`` may be represented in "matrix-free" form as a linear operator
`op` which returns `y = op(x) == A * x` (or its in-place variant `op!(y, x)`).

For the out-of-place `cg`, the search is started with `x .= 0`.

For the in-place `cg!`, the vector `x` is not (re)initialized prior to iteration and
therefore can be used to provide an initial starting guess if an approximate solution is
known.

## Keyword Arguments

The following keywords specify the iteration stopping criteria, stopping if any condition
is met.

- `maxiter` - Maximum number of iterations to perform.
- `abstol` - The absolute error ``δ`` such that ``\\lVert r \\rVert ≤ δ``.
- `reltol` - The relative error ``ϵ`` such that
  ``\\lVert r \\rVert / \\lVert r₀ \\rVert ≤ ϵ``.

where the residual at each iteration is defined as ``r ≡ b - Ax`` and ``r₀`` is the
error prior to the first iteration (i.e. ``r₀ = b`` in the out-of-place call and
``r₀ = b - Ax₀`` for the in-place call).
Effectively, the stopping criterion is ``r ≤ \\operatorname{max}(δ, ϵ\\lVert r₀\\rVert)``.

## Returns

For the in-place call, the returned `info` is a `NamedTuple` with the following fields:

- `converged` = Boolean indicating whether convergence was reached (stopped
  based on absolute/relative error) or not (maximum iterations reached).
- `iter` = Number of iterations performed.
- `abserr` = The absolute error when iterations stopped.
- `relerr` = The relative error when iterations stopped.
- `errhist` = A vector with the history of absolute errors during iteration.

# Example
```jldoctest; setup = (import CMB: cg, cg!)
julia> A = [ 1.0 -0.2
            -0.2  2.0];

julia> b = [0.1, 0.2];

julia> cg(A, b)
2-element Vector{Float64}:
 0.12244897959183675
 0.11224489795918367

julia> Afunc!(y, x) = (y[1] =  1.0x[1] - 0.2x[2];
                       y[2] = -0.2x[1] + 2.0x[2])
Afunc! (generic function with 1 method)

julia> z = zeros(ComplexF64, 2);

julia> cg!(z, Afunc!, complex.(b))
(converged = true, iter = 2, abserr = 0.0, relerr = 0.0, errhist = [0.223606797749979, 0.07089971635974944, 0.0])

julia> z - A \\ b
2-element Vector{ComplexF64}:
 0.0 + 0.0im
 0.0 + 0.0im
```

# Extended Help

`cg` and `cg!` implement the method described by Shewchuk (1994) with two minor changes:

1. An absolute error criterion is added, as described above.

2. The every-50th iteration "restart" procedure is skipped, electing instead to only
   recalculate ``rₙ = b - Axₙ`` if convergence has potentially been reached.
   In such a case, the residual is recalculated, and only if the residual is still less
   than the desired absolute/relative error is the iteration stopped.

## References

- J.R. Shewchuk “An Introduction to the Conjugate Gradient Method Without the Agonizing
  Pain” (1994) URL: <https://www.cs.cmu.edu/~jrs/jrspapers.html>
"""
function cg end

@doc (@doc cg)
function cg! end

function cg(A::AbstractMatrix, b::AbstractVector; kws...)
    x = zeros(eltype(b), size(b)...)
    cg!(x, A, b; kws...)
    return x
end

function cg(op::Base.Callable, b::AbstractVector; kws...)
    x = zeros(eltype(b), size(b)...)
    op! = (x, y) -> (x .= op(y); x)
    cg!(x, op!, b; kws...)
    return x
end

function cg!(x::AbstractVector, A::AbstractMatrix, b::AbstractVector; kws...)
    op! = (x, y) -> mul!(x, A, y)
    return cg!(x, op!, b; kws...)
end

function cg!(x::AbstractVector, op!::Base.Callable, b::AbstractVector{T};
             maxiter::Integer = length(b) * (1 + (T <: Complex)),
             abstol::Real = zero(real(T)),
             reltol::Real = sqrt(eps(real(T)))
            ) where {T}
    r = similar(x)
    q = similar(x)

    op!(r, x); r .= b .- r

    d = copy(r)
    δ₀ = δ′ = real(r' * r)
    ε = max(abs2(reltol) * δ₀, abs2(abstol))

    h = sizehint!([δ₀], maxiter)
    iter, converged = 0, (δ′ <= ε)
    while !converged
        if iter >= maxiter
            break
        end
        op!(q, d)
        α = δ′ / (d' * q)
        x .+= α .* d
        r .-= α .* q
        δ, δ′ = δ′, real(r' * r)
        if δ′ <= ε
            # recalculate the residual from scratch to verify the accuracy
            op!(r, x); r .= b .- r
            δ′ = real(r' * r)
        end
        iter += 1
        push!(h, δ′)
        if δ′ <= ε
            converged = true
            break
        end
        β = δ′ / δ
        d .= r .+ β .* d
    end
    info = (;
            converged = converged,
            iter = iter,
            abserr = sqrt(δ′),
            relerr = sqrt(δ′ / δ₀),
            errhist = sqrt.(h),
           )
    return info
end
