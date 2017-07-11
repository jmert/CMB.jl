module CMB
    import Reexport.@reexport

    include("healpix.jl")
    @reexport using .Healpix

    include("legendre.jl")
    @reexport using .Legendre

    include("harmonics.jl")
    @reexport using .Harmonics

    include("pixelcovariance.jl")
    @reexport using .PixelCovariance

    include("sphere.jl")
    @reexport using .Sphere
end

