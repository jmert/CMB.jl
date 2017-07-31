module Mapping

export
    AbstractPixelizedMap,
    RingPixelizedMap,
    AbstractRingIterator, AbstractPixelIterator, AbstractRingPixelIterator,
    RingIterator, PixelIterator,RingPixelIterator,
    pixeltype, numrings, numpixels, ringlength,
    eachring, eachpixel, eachringpixel

"""
    abstract type AbstractPixelizedMap end

Abstract supertype of all pixelized map descriptions.
"""
abstract type AbstractPixelizedMap end

abstract type AbstractRingIterator end
abstract type AbstractPixelIterator end
abstract type AbstractRingPixelIterator end

MappingIterator = Union{AbstractRingIterator,
                        AbstractPixelIterator,
                        AbstractRingPixelIterator}
Base.iteratorsize(::MappingIterator) = Base.HasLength()

"""
    pixeltype(M)

The data type of a pixel index for the map definition `M`.
"""
function pixeltype end
pixeltype(m::AbstractPixelizedMap) = pixeltype(typeof(m))

struct RingIterator{M<:AbstractPixelizedMap} <: AbstractRingIterator
    m::M
end
struct PixelIterator{M<:AbstractPixelizedMap} <: AbstractPixelIterator
    m::M
end
struct RingPixelIterator{M<:AbstractPixelizedMap} <: AbstractRingPixelIterator
    m::M
    r::Int64
end

eachring(m::M) where {M<:AbstractPixelizedMap} = RingIterator{M}(m)
eachpixel(m::M) where {M<:AbstractPixelizedMap} = PixelIterator{M}(m)
eachringpixel(m::M, r) where {M<:AbstractPixelizedMap} = RingPixelIterator{M}(m, r)

Base.start(I::RingIterator{M}) where {M} = Int64(1)
Base.start(I::PixelIterator{M}) where {M} = Int64(1)
Base.start(I::RingPixelIterator{M}) where {M} = Int64(1)

Base.next(I::RingIterator{M}, state) where {M} = (state, state+1)

Base.done(I::RingIterator{M}, state) where {M} = state > numrings(I.m)
Base.done(I::PixelIterator{M}, state) where {M} = state > numpixels(I.m)
Base.done(I::RingPixelIterator{M}, state) where {M} = state > ringlength(I.m, I.r)

Base.length(I::RingIterator{M}) where {M} = numrings(I.m)
Base.length(I::PixelIterator{M}) where {M} = numpixels(I.m)
Base.length(I::RingPixelIterator{M}) where {M} = ringlength(I.m, I.r)

Base.eltype(::Type{RingIterator{M}}) where {M} = Int64
Base.eltype(::Type{PixelIterator{M}}) where {M} = pixeltype(M)
Base.eltype(::Type{RingPixelIterator{M}}) where {M} = pixeltype(M)

struct RingPixelizedMap{T,Tp} <: AbstractPixelizedMap
    nring::Int64
    npix::Int64
    pixels::Vector{Tp}
    ringptr::Vector{Int64}
    ringtheta::Vector{T}
    pixphi::Vector{T}
end

pixeltype(m::Type{RingPixelizedMap{T,Tp}}) where {T,Tp} = Tp

numrings(m::RingPixelizedMap) = m.nring
numpixels(m::RingPixelizedMap) = m.npix
function ringlength(m::RingPixelizedMap, r::Int)
    n = numrings(m)
    r > n && throw(BoundsError())
    j = (r==n) ? numpixels(m)+1 : m.ringptr[r+1]
    return j - m.ringptr[r]
end

let RI, PI, RPI
    RI{T,Tp} = RingIterator{RingPixelizedMap{T,Tp}}
    PI{T,Tp} = PixelIterator{RingPixelizedMap{T,Tp}}
    RPI{T,Tp} = RingPixelIterator{RingPixelizedMap{T,Tp}}

    Base.next(I::PI{T,Tp},  state) where {T,Tp} = (I.m.pixels[state], state+1)
    Base.next(I::RPI{T,Tp}, state) where {T,Tp} = (I.m.pixels[I.m.ringptr[I.r]+state-1], state+1)
end

function Base.convert(::Type{RingPixelizedMap}, m::RingPixelizedMap)
    return RingPixelizedMap(m.nring,
                            m.npix,
                            copy(m.pixels),
                            copy(m.ringptr),
                            copy(m.ringtheta),
                            copy(m.pixphi))
end

function Base.convert(::Type{RingPixelizedMap{T}}, m::M) where {T,M<:AbstractPixelizedMap}
    Tp = pixeltype(M)
    nring = numrings(m)
    npix = numpixels(m)

    pixels = Vector{Tp}(npix)
    ringptr = Vector{Int64}(nring)
    ringtheta = Vector{T}(nring)
    pixphi = Vector{T}(npix)

    pixoff = 1
    @inbounds for (i,r) in enumerate(eachring(m))
        for p in eachringpixel(m, r)
            pixels[pixoff] = p
            # TODO:
            # pixphi[pixoff] = ...
            pixoff += 1
        end
        ringptr[i] = i==1 ? 1 : ringptr[i-1] + ringlength(m, r)
        # TODO:
        # ringtheta[i] = ...
    end

    return RingPixelizedMap{T,Tp}(nring, npix, pixels, ringptr, ringtheta, pixphi)
end

end
