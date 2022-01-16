using Documenter, CMB
# To get the `@require`d-conditional functions
import JLD, JLD2, MAT

include("FigureExpander.jl")

doctest = "--fix"  in ARGS ? :fix :
          "--test" in ARGS ? true : false

DocMeta.setdocmeta!(CMB, :DocTestSetup, :(using CMB); recursive=true)

makedocs(
    format = Documenter.HTML(
            mathengine = Documenter.MathJax3(
                    # If and/or when Documenter is updated to do use MathJax ≥v3.2.0,
                    # this line can be deleted. (See Documenter.jl#174)
                    url = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.0/es5/tex-svg.js"
            ),
            assets = ["assets/custom.css"],
    ),
    sitename = "CMB Analysis",
    authors = "Justin Willmert",
    modules = [CMB],
    doctest = doctest,
    doctestfilters = Regex[
        r"Ptr{0x[0-9a-f]+}",
        r"[0-9\.]+ seconds \(.*\)"
    ],
    pages = [
        "CMB.jl Documentation" => "index.md",
        "Manual" => [
            "HEALPix Pixelization" => "man/healpix.md",
            "Spherical Functions" => "man/sphere.md",
            "Observing Matrices" => "man/obsmat.md",
            "Pixel-pixel Covariance" => "man/pixelcov.md",
            "References" => "man/references.md"
        ],
        "API Reference" => [
            "Public" => "lib/public.md",
            "Private" => "lib/private.md"
        ]
    ],
    repo = "https://github.com/jmert/CMB.jl/blob/{commit}{path}#L{line}",
)

deploydocs(
    repo = "github.com/jmert/CMB.jl",
    target = "build",
)
