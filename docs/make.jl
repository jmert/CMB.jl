push!(LOAD_PATH, joinpath(dirname(@__FILE__), "..", "src"));
using Documenter, CMB

makedocs(
    format = :html,
    sitename = "CMB Analysis",
    authors = "Justin Willmert",
    modules = [CMB],
    pages = [
        "CMB.jl Documentation" => "index.md",
        "Manual" => [
            "HEALPix Pixelization" => "man/healpix.md",
            "Legendre Polynomials" => "man/legendre.md",
            "References" => "man/references.md"
        ],
        "API Reference" => [
            "Public" => "lib/public.md",
            "Private" => "lib/private.md"
        ]
    ]
)

