module CMB
    import Reexport.@reexport

    include("sphere.jl")
    @reexport using .Sphere

    @reexport using Legendre
    include("sphericalharmonics.jl")
    @reexport using .SphericalHarmonics

    include("healpix.jl")
    @reexport using .Healpix

    include("pixelcovariance.jl")
    @reexport using .PixelCovariance

    include("util.jl")
    @reexport using .Util
end

