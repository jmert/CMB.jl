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
            "Coordinate Systems" => "man/coordinates.md",
            "References" => "man/references.md"
        ],
        "Library Documentation" => [
            "Public" => "lib/public.md",
            "Private" => "lib/private.md"
        ]
    ]
)

