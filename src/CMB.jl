module CMB
    import Reexport.@reexport

    include("healpix.jl")
    @reexport using .Healpix

    include("harmonics.jl")
    @reexport using .Harmonics

    include("pixelcovariance.jl")
    @reexport using .PixelCovariance

    include("sphere.jl")
    @reexport using .Sphere
end

