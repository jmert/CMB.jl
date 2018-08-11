dodeploy = "deploy" in ARGS
dodoctest = "fix" in ARGS ? :fix : !dodeploy

# Only run non-nightlies on Linux
Sys.ARCH === :x86_64 && Sys.KERNEL === :Linux || exit()
get(ENV, "TRAVIS_JULIA_VERSION", "") == "nightly" && exit()

try
    import Pkg
    Pkg.add("Documenter")
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
            julia = "0.7",
            osname = "linux",
            deps = nothing,
            make = nothing
        )
    end
finally
    Pkg.rm("Documenter")
end

