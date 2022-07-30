module CMB
    import Reexport.@reexport
    import Compat.@compat # Compat@v3.21 for @compat import Mod as NewName

    include("numerics.jl")
    include("conventions.jl")

    include("sphere.jl")
    @reexport using .Sphere

    @compat import AssociatedLegendrePolynomials as Legendre
    @reexport using .Legendre
    export Legendre

    @compat import SphericalHarmonicTransforms as SphericalHarmonics
    @reexport using .SphericalHarmonics
    export SphericalHarmonics

    include("healpix.jl")
    @reexport using .Healpix

    include("pixelizations.jl")
    @reexport using .Pixelizations

    include("pixelcovariance.jl")
    @reexport using .PixelCovariance

    include("fileio.jl")
    @reexport using .Files

    include("precompile.jl")
    _precompile_()
end

