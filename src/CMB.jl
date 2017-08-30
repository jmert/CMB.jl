module CMB
    import Reexport.@reexport

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

    include("sphere.jl")
    @reexport using .Sphere

    include("bkutils.jl")
end

