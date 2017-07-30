module Mapping

export
    AbstractPixelizedMap,
    RingPixelizedMap,
    AbstractRingIterator, AbstractPixelIterator, AbstractRingPixelIterator,
    RingIterator, PixelIterator,RingPixelIterator,
    ringtype, pixeltype, numrings, numpixels, ringlength,
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

struct RingIterator{M<:AbstractPixelizedMap} <: AbstractRingIterator
    m::M
end
struct PixelIterator{M<:AbstractPixelizedMap} <: AbstractPixelIterator
    m::M
end
struct RingPixelIterator{M<:AbstractPixelizedMap,R} <: AbstractRingPixelIterator
    m::M
    r::R
end

eachring(m::M) where {M<:AbstractPixelizedMap} = RingIterator{M}(m)
eachpixel(m::M) where {M<:AbstractPixelizedMap} = PixelIterator{M}(m)
eachringpixel(m::M, r::R) where {M<:AbstractPixelizedMap,R} = RingPixelIterator{M,R}(m, r)

Base.length(I::RingIterator{M}) where {M} = numrings(I.m)
Base.length(I::PixelIterator{M}) where {M} = numpixels(I.m)
Base.length(I::RingPixelIterator{M}) where {M} = ringlength(I.m, I.r)

Base.eltype(::Type{RingIterator{M}}) where {M} = ringtype(M)
Base.eltype(::Type{PixelIterator{M}}) where {M} = pixeltype(M)
Base.eltype(::Type{RingPixelIterator{M}}) where {M} = pixeltype(M)

struct RingPixelizedMap{Tp,Tr} <: AbstractPixelizedMap
    pixels::AbstractVector{Tp}
    rings::AbstractVector{Tr}
    ringptr::AbstractVector{Int64}
end

ringtype(m::Type{RingPixelizedMap{Tp,Tr}}) where {Tp,Tr} = Tr
pixeltype(m::Type{RingPixelizedMap{Tp,Tr}}) where {Tp,Tr} = Tp

numrings(m::RingPixelizedMap) = length(m.rings)
numpixels(m::RingPixelizedMap) = length(m.pixels)
function ringlength(m::RingPixelizedMap, r::Int)
    n = numrings(m)
    r > n && throw(BoundsError())
    j = (r==n) ? numpixels(m) : m.ringptr[i+1]
    return j - m.ringptr[i] + 1
end

let RI, PI, RPI
    RI{Tp,Tr} = RingIterator{RingPixelizedMap{Tp,Tr}}
    PI{Tp,Tr} = PixelIterator{RingPixelizedMap{Tp,Tr}}
    RPI{Tp,Tr} = RingPixelIterator{RingPixelizedMap{Tp,Tr},Tr}

    Base.start(I::RI{Tp,Tr})  where {Tp,Tr} = Int64(1)
    Base.start(I::PI{Tp,Tr})  where {Tp,Tr} = Int64(1)
    Base.start(I::RPI{Tp,Tr}) where {Tp,Tr} = Int64(0)

    Base.next(I::RI{Tp,Tr},  state) where {Tp,Tr} = (I.m.rings[state], state+1)
    Base.next(I::PI{Tp,Tr},  state) where {Tp,Tr} = (I.m.pixels[state], state+1)
    Base.next(I::RPI{Tp,Tr}, state) where {Tp,Tr} = (I.m.pixels[I.m.ringptr[I.r]+state], state+1)

    Base.done(I::RI{Tp,Tr},  state) where {Tp,Tr} = state > numrings(I.m)
    Base.done(I::PI{Tp,Tr},  state) where {Tp,Tr} = state > numpixels(I.m)
    Base.done(I::RPI{Tp,Tr}, state) where {Tp,Tr} = state ≥ ringlength(I.m, I.r)
end

end
