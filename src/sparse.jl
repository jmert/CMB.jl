export outer, quadprod

"""
Computes the outer product between a given column of a sparse matrix and a vector.
"""
function outer end

"""
    outer(A::SparseMatrixCSC, n::Integer, w::AbstractVector)

Performs the equivalent of ``\\vec a_n \\vec w^\\dagger`` where ``\\vec a_n`` is the
column `A[:,n]`.
"""
function outer(A::SparseMatrixCSC{Tv,Ti}, n::Integer, w::AbstractVector{Tv}) where {Tv,Ti}
    colptrn = nzrange(A, n)
    rowvalA = rowvals(A)
    nzvalsA = nonzeros(A)

    nnza = length(colptrn)
    nnzw = length(w)
    numnz = nnza * nnzw

    colptr = Vector{Ti}(nnzw+1)
    rowval = Vector{Ti}(numnz)
    nzvals = Vector{Tv}(numnz)

    idx = 0
    @inbounds for jj = 1:nnzw
        colptr[jj] = idx + 1

        wv = conj(w[jj])
        iszero(wv) && continue

        for ii = colptrn
            idx += 1
            rowval[idx] = rowvalA[ii] # copy row index from A
            nzvals[idx] = wv * nzvalsA[ii] # outer product values
        end
    end
    @inbounds colptr[nnzw+1] = idx + 1
    return SparseMatrixCSC(size(A,1), nnzw, colptr, rowval, nzvals)
end

"""
    outer(w::AbstractVector, A::SparseMatrixCSC, n::Integer)

Performs the equivalent of ``\\vec w \\vec{a}_n^\\dagger`` where ``\\vec a_n`` is the
column `A[:,n]`.
"""
function outer(w::AbstractVector{Tv}, A::SparseMatrixCSC{Tv,Ti}, n::Integer) where {Tv,Ti}
    colptrn = nzrange(A, n)
    rowvalA = rowvals(A)
    nzvalsA = nonzeros(A)

    nnza = length(colptrn)
    nnzw = length(w)
    numnz = nnza * nnzw

    colptr = zeros(Ti, size(A,1)+1)
    rowval = Vector{Ti}(numnz)
    nzvals = Vector{Tv}(numnz)

    idx = 0
    @inbounds colptr[1] = 1 # col 1 always at index 1
    @inbounds for jj = colptrn
        av = conj(nzvalsA[jj])
        rv = rowvalA[jj]

        for ii = 1:nnzw
            wv = w[ii]
            iszero(wv) && continue

            idx += 1
            colptr[rv+1] += 1 # count num of entries in column
            rowval[idx] = ii
            nzvals[idx] = w[ii] * av # outer product values
        end
    end
    cumsum!(colptr, colptr) # offsets are sum of all previous

    return SparseMatrixCSC(nnzw, size(A,1), colptr, rowval, nzvals)
end

"""
    quadprod(A, b, n, dir=:col)

Computes the quadratic product ``ABA^T`` efficiently for the case where ``B`` is all zero
except for the `n`th column or row vector `b`, for `dir = :col` or `dir = :row`,
respectively.
"""
function quadprod(A, b, n, dir::Symbol=:col)
    if dir == :col
        w = A * b
        return outer(w, A, n)
    elseif dir == :row
        w = A * b
        return outer(A, n, w)
    else
        error("Unrecognized direction `dir = $(repr(dir))`.")
    end
end
