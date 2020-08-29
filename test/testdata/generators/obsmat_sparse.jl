using SparseArrays

const R = sparse(
        [
            0.1 0.2 0.0 0.0 0.0 0.0
            0.0 0.0 0.5 0.5 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 1.0
        ])

const base = abspath(joinpath(@__DIR__, ".."))

module TestJLD
    using JLD, SparseArrays
    import ..R, ..base
    jldopen(joinpath(base, "obsmat_sparse.jld"), "w") do file
        write(file, "R", R, compress = false)
    end
end

module TestJLD2
    using JLD2
    import ..R, ..base
    jldopen(joinpath(base, "obsmat_sparse.jld2"), "w") do file
        write(file, "R", R)
    end
end
