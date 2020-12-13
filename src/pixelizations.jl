module Pixelizations

export AbstractPixelization, ArbitraryPixelization, HealpixPixelization, RADecPixelization,
       parse_pixelization, export_pixelization, pointing

@static if VERSION < v"1.6.0-DEV.1083"
    import Compat # Compat@v3.18 for reinterpret(reshape, ...)
end
using StaticArrays
using ..Healpix
using ..Sphere
import ..Sphere: latlon, unsafe_colataz, cartvec

const Vec3{T} = SVector{3,T} where T

"""
    AbstractPixelization

The supertype of various pixelization formats.
"""
abstract type AbstractPixelization end

"""
    struct ArbitraryPixelization <: AbstractPixelization

An arbitrary pixelization format (i.e. one without any implied structure) where pixel
centers are given as a vector of unit vectors.
"""
struct ArbitraryPixelization{T,V<:AbstractVector{Vec3{T}}} <: AbstractPixelization
    rvecs::V
end
function Base.show(io::IO, mime::MIME"text/plain", pix::ArbitraryPixelization)
    println(io, ArbitraryPixelization, " with ", length(pix.rvecs), " pixels:")
    Base.print_array(io, pix.rvecs)
    return nothing
end

"""
    struct HealpixPixelization <: AbstractPixelization

A HEALPix pixelization of a particular resolution (``N_\\mathrm{side}``). May optionally
also describe a subset of the sphere, containing an explicit list of pixels to subset.
"""
struct HealpixPixelization{V} <: AbstractPixelization
    nside::Int64
    pixels::V
    function HealpixPixelization(nside::Integer, pixels::AbstractVector)
        npix = nside2npix(nside)
        minpix, maxpix = extrema(pixels)
        checkhealpix(nside, minpix)
        checkhealpix(nside, maxpix)
        return new{typeof(pixels)}(Int64(nside), pixels)
    end
    function HealpixPixelization(nside::Integer)
        n = convert(Int64, nside)
        pix = zero(n):nside2npix(n)-1
        return new{typeof(pix)}(n, pix)
    end
end
function Base.show(io::IO, mime::MIME"text/plain", pix::HealpixPixelization)
    println(io, "Nside = ", pix.nside, " ", HealpixPixelization, " with ",
            length(pix.pixels), " pixels:")
    Base.print_array(io, pix.pixels)
    return nothing
end

"""
    struct RADecPixelization <: AbstractPixelization

An equidistant cylindrical projection described in terms of a square grid of pixels by
their right ascension (RA) and declination (Dec) coordinates (in degrees).
"""
struct RADecPixelization{V<:AbstractRange} <: AbstractPixelization
    ra::V
    dec::V
    roworder::Bool
    function RADecPixelization(ra::AbstractRange, dec::AbstractRange, roworder::Bool)
        ramin, ramax = extrema(ra) .+ ((-1, 1) .* step(ra))
        decmin, decmax = extrema(dec) .+ ((-1, 1) .* step(dec))
        if ramax - ramin > 360
            error("RA range must span 360 degrees or less; got [$ramin, $ramax]")
        end
        if !(-90 ≤ decmin ≤ decmax ≤ 90)
            error("Dec range must be in [-90, 90]; got [$decmin, $decmax]")
        end
        ra′, dec′ = promote(ra, dec)
        return new{typeof(ra′)}(ra′, dec′, roworder)
    end
end
function Base.show(io::IO, mime::MIME"text/plain", pix::RADecPixelization)
    io′ = IOContext(io, :compact => true)
    ra, dec = pix.ra, pix.dec
    Δra, Δdec = step(ra), step(dec)
    nra, ndec = length(ra), length(dec)
    println(io, RADecPixelization, " with ", nra, " × ", ndec, " = ", nra*ndec, " pixels:")

    println(io′, "   RA ∈ [", first(ra)-0.5Δra, ", ", last(ra)+0.5Δra, "] (step ", Δra, ")")
    print(io′, "  Dec ∈ [", first(dec)-0.5Δdec, ", ", last(dec)+0.5Δdec, "] (step ", Δdec, ")")
    return nothing
end

"""
Parses a number of recognizable pixelization specifications, returning a subtype of
[`CMB.AbstractPixelization`](@ref AbstractPixelization).

See also [`export_pixelization`](@ref export_pixelization).
"""
function parse_pixelization end

"""
Exports a given instance of [`AbstractPixelization`](@ref AbstractPixelization)
to a "simple" data type representation suitable for export to disk.

See also [`parse_pixelization`](@ref parse_pixelization).
"""
function export_pixelization end

# Direct interpretation of raw unit vectors
# pass-through of vectors of SVectors
parse_pixelization(pixelspec::AbstractVector{Vec3{T}}) where {T} = ArbitraryPixelization(pixelspec)
# convert 3×N matrices into vector of SVectors
"""
    parse_pixelization(pixelspec::AbstractMatrix)
    parse_pixelization(pixelspec::AbstractVector{SVector{3,<:Any}}})

Interprets a collection of unit 3-vectors as a set of arbitrary pixels on the sphere,
returning an instance of [`ArbitraryPixelization`](@ref). The input may be either a vector
of `SVector{3}` unit vectors or a 3×N matrix where each column is interpreted as a single
unit vector.
"""
function parse_pixelization(pixelspec::AbstractMatrix)
    size(pixelspec, 1) == 3 ||
            throw(DimensionMismatch(string(
                    "Cannot interpret ", join(string.(size(pixelspec)), "×"),
                    " matrix as vector of pointing vectors")))
    T = eltype(pixelspec)
    return ArbitraryPixelization(reinterpret(reshape, Vec3{T}, pixelspec))
end
"""
    export_pixelization(pix::ArbitraryPixelization)

Returns a 3×N matrix that describes the given set of arbitrary pixel centers, `pix`.
See also [`parse_pixelization`](@ref parse_pixelization(::AbstractMatrix)).
"""
function export_pixelization(pix::ArbitraryPixelization)
    return reinterpret(reshape, eltype(eltype(pix.rvecs)), pix.rvecs)
end

# Parsing of dictionary specification into unit vector vectors

# Top-level dispatch mechanism
"""
    parse_pixelization(pixelspec::AbstractDict{String})

Attempts to parse an `AbstractDict` containing a pixelization description.  The dictionary
is expected to contain a field `"type"` which specifies the pixelization scheme name, which
is used to dispatch for further format-specific processing.
"""
function parse_pixelization(pixelspec::AbstractDict{String})
    haskey(pixelspec, "type") || error("pixelization specification identification requires a `\"type\"` field")
    type = Symbol(pixelspec["type"])
    pix = parse_pixelization(Val(type), pixelspec)
    return pix::AbstractPixelization
end

# HEALPix pixelization
"""
    parse_pixelization(::Val{:healpix}, pixelspec::AbstractDict{String,<:Any})

Parses a dictionary describing a HEALPix pixelization, returning an instance of
[`HealpixPixelization`](@ref).

`pixelspec` is must conform to the following format to be parsed:
```julia
Dict("type" => "healpix",
     "nside" => #= Nside value =#,
     "pixels" => #= vector of HEALPix pixel indices, or if not provided implied to be a
                    full-sky grid spanning pixels `0:nside2npix(nside)-1` =#
    )
```
"""
function parse_pixelization(::Val{:healpix}, pixelspec::AbstractDict{String})
    nside = convert(Int64, pixelspec["nside"])::Int64
    if haskey(pixelspec, "pixels")
        pixind = convert(Vector{Int64}, pixelspec["pixels"])::Vector{Int64}
        return HealpixPixelization(nside, pixind)
    else
        return HealpixPixelization(nside)
    end
end

"""
    export_pixelization(pix::HealpixPixelization)

Returns a `Dict{String}` appropriate for serializing to disk that describes the given
HEALPix pixelization, `pix`. See also
[`parse_pixelization`](@ref parse_pixelization(::Val{:healpix}, ::AbstractDict{String})).
"""
function export_pixelization(pix::HealpixPixelization)
    pixelspec = Dict{String,Any}("type" => "healpix",
                                 "nside" => pix.nside)
    if length(pix.pixels) == nside2npix(pix.nside)
        return pixelspec
    else
        pixelspec["pixels"] = convert(Vector{Int64}, pix.pixels)
        return pixelspec
    end
end

# RA/Dec pixelizations
"""
    parse_pixelization(::Val{:radec}, pixelspec::AbstractDict{String})

Parses a dictionary describing an equidistance cylindrical projection given in terms of
right ascension (RA) and declination (Dec) coordinates, returning an instance of
[`RADecPixelization`](@ref).

`pixelspec` is must conform to the following format to be parsed:
```julia
Dict("type" => "radec",
     "ra" => #= vector of uniformly-spaced pixel centers along the RA axis, in degrees =#,
     "dec" => #= vector of uniformly-spaced pixel centers along the Dec axis, in degrees =#,
     "order" => #= either "row" or "col" for row/column-major unravelling =#,
    )
```
"""
function parse_pixelization(::Val{:radec}, pixelspec::AbstractDict{String})
    ra = pixelspec["ra"]::AbstractVector{<:Real}
    dec = pixelspec["dec"]::AbstractVector{<:Real}
    order = pixelspec["order"]::String

    roworder = order == "row" ? true :
               order == "col" ? false :
               error("Unknown ordering `", order, "`")

    ra′ = range(first(ra), last(ra), length = length(ra))
    dec′ = range(first(dec), last(dec), length = length(dec))
    ra ≈ ra′ || error("RA coordinates not a uniformly spaced grid")
    dec ≈ dec′ || error("Dec coordinates not a uniformly spaced grid")
    return RADecPixelization(ra′, dec′, roworder)
end

"""
    export_pixelization(pix::RADecPixelization)

Returns a `Dict{String}` appropriate for serializing to disk that describes the given
RA/Dec pixelization, `pix`. See also
[`parse_pixelization`](@ref parse_pixelization(::Val{:radec}, ::AbstractDict{String})).
"""
function export_pixelization(pix::RADecPixelization)
    return Dict{String,Any}("type" => "radec",
                            "ra" => collect(pix.ra),
                            "dec" => collect(pix.dec),
                            "order" => pix.roworder ? "row" : "col"
                           )
end

"""
    pointing(pix::AbstractPixelization)

Returns a `Vector` of `SVector{3}` unit-vectors pointing to the centers of the pixels
described by `pix`.
"""
function pointing end

function pointing(pix::ArbitraryPixelization)
    return pix.rvecs
end

function pointing(pix::HealpixPixelization)
    nside = pix.nside
    pixels = pix.pixels
    rvec = similar(pixels, Vec3{float(eltype(pixels))})
    @inbounds for ii in eachindex(pix.pixels)
        rvec[ii] = Healpix.unsafe_pix2vec(nside, pixels[ii])
    end
    return rvec
end

function pointing(pix::RADecPixelization)
    n = length(pix.ra) * length(pix.dec)
    T = eltype(pix.ra)
    rvec = Vector{Vec3{T}}(undef, n)
    if pix.roworder
        kk = 0
        @inbounds for δ in pix.dec
            @simd ivdep for λ in pix.ra
                rvec[kk+=1] = cartvec(unsafe_colataz(δ, λ))
            end
        end
    else
        kk = 0
        @inbounds for λ in pix.ra
            @simd ivdep for δ in pix.dec
                rvec[kk+=1] = cartvec(unsafe_colataz(δ, λ))
            end
        end
    end
    return rvec
end

end
