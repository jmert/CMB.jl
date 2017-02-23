push!(LOAD_PATH, joinpath(dirname(@__FILE__), "..", "src"));
using Documenter, CMB

makedocs(
    format = :html,
    sitename = "CMB Analysis",
    authors = "Justin Willmert",
    pages = [
        "CMB.jl Documentation" => "index.md",
        "Manual" => [
            "References" => "man/references.md"
        ],
        "Library Documentation" => [
            "Public" => "lib/public.md",
            "Private" => "lib/private.md"
        ]
    ]
)

