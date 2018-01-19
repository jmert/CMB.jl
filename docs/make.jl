dodeploy = "deploy" in ARGS
dodoctest = !dodeploy

if dodeploy
    Pkg.add("Documenter")
end

push!(LOAD_PATH, joinpath(dirname(@__FILE__), "..", "src"));
using Documenter, CMB

makedocs(
    format = :html,
    sitename = "CMB Analysis",
    authors = "Justin Willmert",
    modules = [CMB],
    doctest = dodoctest,
    doctestfilters = Regex[
        r"Ptr{0x[0-9a-f]+}",
        r"[0-9\.]+ seconds \(.*\)"
    ],
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

if dodeploy
    deploydocs(
        target = "build",
        repo = "github.com/jmert/CMB.jl.git",
        julia = "0.6",
        osname = "linux",
        deps = nothing,
        make = nothing
    )
end

