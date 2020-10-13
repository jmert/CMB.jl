module Pixelizations

@static if VERSION < v"1.6.0-DEV.1083"
    using Compat # requires v3.18 for reinterpret(reshape, ...)
end
using StaticArrays
using ..Healpix
using ..Sphere
import ..Sphere: latlon, colataz, cartvec

const Vec3{T} = SVector{3,T} where T
const VecVec3{T} = Vector{V} where {T, V <: Vec3{T}}

abstract type AbstractPixelization end

"""
    Pixelization

A type used to signal a particular pixelization format, useful for defining dispatch over
pixelization schemes.
"""
struct Pixelization{T} <: AbstractPixelization
    function Pixelization(name)
        return new{Symbol(name)}()
    end
    function Pixelization{name}() where {name}
        return new{Symbol(name)}()
    end
end

"""
    pixelization(pix::Pixelization)

Returns the name of the given pixelization as a `Symbol`.
"""
pixelization(::Pixelization{T}) where {T} = T
Base.show(io::IO, pix::Pixelization) = print(io, pixelization(pix), " pixelization")

# Direct interpretation of raw unit vectors
# pass-through of vectors of SVectors
Base.parse(::Type{Pixelization}, pixelspec::VecVec3{T}) where {T} = pixelspec
# convert 3×N matrices into vector of SVectors
function Base.parse(::Type{Pixelization}, pixelspec::Matrix{T}) where {T}
    size(pixelspec, 1) == 3 ||
            throw(DimensionMismatch(string(
                    "Cannot interpret ", join(string.(size(pixelspec)), "×"),
                    " matrix as vector of pointing vectors")))
    if isbitstype(T)
        return reinterpret(reshape, Vec3{T}, pixelspec)
    else
        rvec = similar(Vector{Vec3{T}}, axes(pixelspec, 2))
        @inbounds for ii in axes(pixelspec, 2)
            rvec[ii] = Vec3{T}(pixelspec[1,ii], pixelspec[2,ii], pixelspec[3,ii])
        end
        return rvec
    end
end

# Parsing of dictionary specification into unit vector vectors

# Top-level dispatch mechanism
"""
    parse(CMB.Pixelization, pixelspec) -> rvec

Parses a number of recognizable pixelization specifications, returning `rvec`, a vector of
unit vectors (`Vector{SVector{3,T}} where T`) pointing to pixel centers on the unit sphere.
"""
function Base.parse(::Type{Pixelization}, pixelspec::AbstractDict)
    haskey(pixelspec, "type") || error("pixelization specification identification requires a `\"type\"` field")
    type = Symbol(pixelspec["type"])
    rvec = parse(Pixelization(type), pixelspec)
    return rvec
end

# HEALPix pixelization
function Base.parse(::Pixelization{:healpix}, pixelspec::AbstractDict)
    nside = convert(Int, pixelspec["nside"])::Int
    if haskey(pixelspec, "pixels")
        pixind = convert(Vector{Int}, pixelspec["pixels"])::Vector{Int}
        return pix2vec.(nside, pixind)
    else
        pixind = 0:nside2npix(nside)-1
        return pix2vec.(nside, pixind)
    end
end

# RA/Dec pixelizations
function Base.parse(::Pixelization{:radec}, pixelspec::AbstractDict)
    radec = pixelspec["pixels"]
    size(radec, 1) == 2 ||
            error("Cannot interpret ", join(string.(size(radec)), "×"),
            " matrix as vector of RA/Dec pairs")
    T = float(eltype(radec))
    rvec = similar(Vector{Vec3{T}}, axes(radec, 2))
    @inbounds for ii in axes(radec, 2)
        rvec[ii] = cartvec(colataz(radec[2, ii], radec[1, ii]))
    end
    return rvec
end

end
