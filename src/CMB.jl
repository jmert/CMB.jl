module CMB
    import Reexport.@reexport

    include("sphere.jl")
    @reexport using .Sphere

    include("legendre.jl")
    @reexport using .Legendre

    include("harmonics.jl")
    @reexport using .Harmonics

    include("healpix.jl")
    @reexport using .Healpix

    include("mapping/mapping.jl")
    @reexport using .Mapping

    include("pixelcovariance.jl")
    @reexport using .PixelCovariance

    include("bkutils.jl")
    include("util.jl")
end

