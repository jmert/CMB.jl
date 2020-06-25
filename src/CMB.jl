module CMB
    import Reexport.@reexport

    include("sphere.jl")
    @reexport using .Sphere

    @reexport using Legendre
    include("harmonics.jl")
    @reexport using .Harmonics

    include("healpix.jl")
    @reexport using .Healpix

    include("pixelcovariance.jl")
    @reexport using .PixelCovariance

    include("util.jl")
    @reexport using .Util
end

