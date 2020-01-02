using Documenter, CMB

doctest = "--fix" in ARGS ? :fix : true

DocMeta.setdocmeta!(CMB, :DocTestSetup, :(using CMB); recursive=true)

makedocs(
    format = Documenter.HTML(mathengine=MathJax()),
    sitename = "CMB Analysis",
    authors = "Justin Willmert",
    modules = [CMB],
    doctest = false,
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
    ],
    repo = "https://github.com/jmert/CMB.jl.git/blob/{commit}{path}#L{line}",
)

deploydocs(
    repo = "github.com/jmert/CMB.jl.git",
    target = "build",
)
