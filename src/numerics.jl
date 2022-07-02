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

