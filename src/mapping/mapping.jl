module Mapping
    # Traits
    export AbstractCoordinateSystem, AbstractPolarizationConvention,
        LatLonCoordinates, ColatAzCoordinates, UnknownCoordinates,
        IAUPolarization, HealpixPolarization, UnknownPolarization,
        coordsys, polconv

    # Abstract types
    export AbstractPixel, AbstractPixelizedMap,
        pixeltype

    # Concrete map types
    export SimplePixel, RingPixelizedMap,
        nrings, npixels, ringlength,
        eachring, eachpixel, eachringpixel,
        pixelphi, pixeltheta, pixellon, pixellat

    include("base.jl")

    # Specific definitions
    export ECPPixel, ECPMapPatchDefn, HealpixPixel, HealpixMapDefn,
        BicepMapDefn, BicepExtMapDefn

    # HealpixMapNNNN types exported in definition loop
    include("specific.jl")
end

