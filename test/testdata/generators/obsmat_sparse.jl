using SparseArrays

const R = sparse(
        [
            0.1 0.2 0.0 0.0 0.0 0.0
            0.0 0.0 0.5 0.5 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 0.0
            0.0 0.0 0.0 0.0 0.0 1.0
        ])

# RHS pixels are described by a complex data structure
const pixr = Dict{String,Any}(
        "type" => "dummy_pixels",
        "index" => Int64.(collect(0:5)),
        "sub" => Dict{String,Any}("extra" => Int8(1))
        )
# LHS is some very simple indexing scheme
const pixl = Int64.(collect(1:4))
const fields = "T"

const base = abspath(joinpath(@__DIR__, ".."))

#=
module TestJLD
    using JLD, SparseArrays
    import ..R, ..fields, ..pixr, ..pixl, ..base
    jldopen(joinpath(base, "obsmat_sparse.jld"), "w") do file
        write(file, "R", R)
        write(file, "pixels_right", pixr)
        write(file, "pixels_left", pixl)
        write(file, "fields", fields)
    end
end
=#

module TestJLD2
    using JLD2
    import ..R, ..fields, ..pixr, ..pixl, ..base
    jldopen(joinpath(base, "obsmat_sparse.jld2"), "w") do file
        write(file, "R", R)
        write(file, "pixels_right", pixr)
        write(file, "pixels_left", pixl)
        write(file, "fields", fields)
    end
end
