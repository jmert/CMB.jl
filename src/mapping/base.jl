####
# Traits

"""
    abstract type AbstractCoordinateSystem end

Abstract supertype of all coordinate system trait types.
"""
abstract type AbstractCoordinateSystem end

"""
    struct LatLonCoordinates <: AbstractCoordinateSystem end

Type which indicates the use of latitude/longitude mapping coordinates, where latitude
measures degrees above or below the equator, and longitude measure degrees east or west of
the prime meridian. See also [`ColatAzCoordinates`](@ref).
"""
struct LatLonCoordinates <: AbstractCoordinateSystem end

"""
    struct ColatAzCoordinates <: AbstractCoordinateSystem end

Type which indicates the use of colatitude/azimuth mapping coordinates, where colatitude
measures radians south of the North Pole, and azimuth measures radians east of the prime
meridian. See also [`LatLonCoordinates`](@ref).
"""
struct ColatAzCoordinates <: AbstractCoordinateSystem end

"""
    struct UnknownCoordinates <: AbstractCoordinateSystem end

Type indicating the use of an unknown coordinate system. See also
[`LatLonCoordinates`](@ref) and [`ColatAzCoordinates`](@ref).
"""
struct UnknownCoordinates <: AbstractCoordinateSystem end

coordsys(::T) where T<:Union{C,Type{C}} where {C<:AbstractCoordinateSystem} = C()

"""
    abstract type AbstractPolarizationConvention end

Abstract supertype of all polarization convention trait types.
"""
abstract type AbstractPolarizationConvention end

"""
    struct NoPolarization <: AbstractPolarizationConvention end

Type indicating that only unpolarized data is supported. See also
[`IAUPolarization`](@ref), [`HealpixPolarization`](@ref), and [`UnknownPolarization`](@ref).
"""
struct NoPolarization <: AbstractPolarizationConvention end

"""
    struct IAUPolarization <: AbstractPolarizationConvention end

Type indicating the use of the IAU polarization conventions.
[`NoPolarization`](@ref), [`HealpixPolarization`](@ref), and [`UnknownPolarization`](@ref).
"""
struct IAUPolarization <: AbstractPolarizationConvention end

"""
    struct HealpixPolarization <: AbstractPolarizationConvention end

Type indicating the use of the Healpix polarization conventions.
[`NoPolarization`](@ref), [`IAUPolarization`](@ref), and [`UnknownPolarization`](@ref).
"""
struct HealpixPolarization <: AbstractPolarizationConvention end

"""
    struct UnknownPolarization <: AbstractPolarizationConvention end

Type indicating the use of an unknown polarization convention. See also
[`NoPolarization`](@ref), [`IAUPolarization`](@ref), and [`HealpixPolarization`](@ref).
"""
struct UnknownPolarization <: AbstractPolarizationConvention end

polconv(::T) where T<:Union{P,Type{P}} where {P<:AbstractPolarizationConvention} = P()

####
# Map interface types

abstract type AbstractPixel{Cs} end

"""
    coordsys(p::AbstractPixel)

Describes the coordinate system associated with pixel type `p`. See also
[`LatLonCoordinates`](@ref) and [`ColatAzCoordinates`](@ref).
"""
coordsys(p::AbstractPixel) = coordsys(typeof(p))
coordsys(::Type{<:AbstractPixel{Cs}}) where {Cs} = UnknownCoordinates()
coordsys(::Type{<:AbstractPixel{Cs}}) where {Cs<:AbstractCoordinateSystem} = Cs()

abstract type AbstractPixelizedMap{Tp,Pc} end

"""
    coordsys(m::AbstractPixelizedMap)

Returns the coordinate system associated with the pixel type of the given map `m`.
Equivalent to `coordsys(pixeltype(m))`.
"""
coordsys(m::AbstractPixelizedMap) = coordsys(pixeltype(m))

"""
    pixeltype(m::AbstractPixelizedMap)

Returns the pixel type of the given map `m`.
"""
pixeltype(m::AbstractPixelizedMap) = pixeltype(typeof(m))
pixeltype(::Type{<:AbstractPixelizedMap{Tp}}) where {Tp} = Tp

"""
    polconv(m::AbstractPixelizedMap)

Returns the polarization convention assumed by the given map `m`. See also
[`IAUPolarization`](@ref) and [`HealpixPolarization`](@ref).
"""
polconv(m::AbstractPixelizedMap) = polconv(typeof(m))
polconv(::Type{<:AbstractPixelizedMap{Tp,Pc} where {Tp}}) where {Pc} = UnknownPolarization()
polconv(::Type{<:AbstractPixelizedMap{Tp,Pc} where {Tp}}) where
        {Pc<:AbstractPolarizationConvention} = Pc()

function pixelphi(m::AbstractPixelizedMap{Tp}, p) where Tp<:AbstractPixel{LatLonCoordinates}
    λ = pixellon(m, p)
    if λ < 0
        λ += 360.0
    end
    return deg2rad(λ)
end

function pixeltheta(m::AbstractPixelizedMap{Tp}, p) where Tp<:AbstractPixel{LatLonCoordinates}
    δ = pixellat(m, p)
    return deg2rad(90.0 - δ)
end

function pixellon(m::AbstractPixelizedMap{Tp}, p) where Tp<:AbstractPixel{ColatAzCoordinates}
    ϕ = pixelphi(m, p)
    if ϕ > π
        ϕ -= 2π
    end
    return rad2deg(ϕ)
end

function pixellat(m::AbstractPixelizedMap{Tp}, p) where Tp<:AbstractPixel{ColatAzCoordinates}
    θ = pixeltheta(m, p)
    return rad2deg(π - θ)
end

####
# Pixel iteration

abstract type AbstractRingIterator end
abstract type AbstractPixelIterator end
abstract type AbstractRingPixelIterator end

MappingIterator = Union{AbstractRingIterator,
                        AbstractPixelIterator,
                        AbstractRingPixelIterator}
Base.iteratorsize(::I) where I<:MappingIterator = Base.HasLength()

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

Base.done(I::RingIterator{M}, state) where {M} = state > nrings(I.m)
Base.done(I::PixelIterator{M}, state) where {M} = state > npixels(I.m)
Base.done(I::RingPixelIterator{M}, state) where {M} = state > ringlength(I.m, I.r)

Base.length(I::RingIterator{M}) where {M} = nrings(I.m)
Base.length(I::PixelIterator{M}) where {M} = npixels(I.m)
Base.length(I::RingPixelIterator{M}) where {M} = ringlength(I.m, I.r)

Base.eltype(::Type{RingIterator{M}}) where {M} = Int64
Base.eltype(::Type{PixelIterator{M}}) where {M} = eltype(M)
Base.eltype(::Type{RingPixelIterator{M}}) where {M} = eltype(M)


####
# Base containers

struct SimplePixel{T,L,Cs} <: AbstractPixel{Cs}
    pixel::T
end

Base.eltype(::P) where P<:Union{Tp,Type{Tp}} where Tp<:SimplePixel{T,L,Cs} where {T,L,Cs} = T
Base.convert(::Type{T}, p::SimplePixel{T,L,Cs}) where {T,L,Cs} = p.pixel

struct RingPixelizedMap{T,Tp,Pc} <: AbstractPixelizedMap{Tp,Pc}
    mapdefn::AbstractPixelizedMap{Tp,Pc}
    pixels::Vector{T}
    ringptr::Vector{Int64}
    ringtheta::Vector{Float64}
    pixphi::Vector{Float64}
end

nrings(m::RingPixelizedMap) = length(m.ringptr) - 1
npixels(m::RingPixelizedMap) = length(m.pixels)
function ringlength(m::RingPixelizedMap, r::Int)
    n = nrings(m)
    r > nrings(m) && throw(BoundsError())
    return m.ringptr[r+1] - m.ringptr[r]
end

let RI, PI, RPI
    RI{T,Tp,Pc} = RingIterator{RingPixelizedMap{T,Tp,Pc}}
    PI{T,Tp,Pc} = PixelIterator{RingPixelizedMap{T,Tp,Pc}}
    RPI{T,Tp,Pc} = RingPixelIterator{RingPixelizedMap{T,Tp,Pc}}

    Base.next(I::PI{T,Tp,Pc},  state) where {T,Tp,Pc} = (I.m.pixels[state], state+1)
    Base.next(I::RPI{T,Tp,Pc}, state) where {T,Tp,Pc} = (I.m.pixels[I.m.ringptr[I.r]+state-1], state+1)
end

#function Base.convert(::Type{RingPixelizedMap}, m::RingPixelizedMap)
#    return RingPixelizedMap(copy(m.pixels),
#                            copy(m.ringptr),
#                            copy(m.ringtheta),
#                            copy(m.pixphi))
#end

function Base.convert(::Type{RingPixelizedMap}, m::AbstractPixelizedMap{Tp,Pc}) where {Tp,Pc}
    nring = nrings(m)
    npix = npixels(m)
    ptype = eltype(Tp)

    pixels = Vector{ptype}(npix)
    ringptr = Vector{Int64}(nring+1)
    ringtheta = Vector{Float64}(nring)
    pixphi = Vector{Float64}(npix)

    pixoff = 1
    for (i,r) in enumerate(eachring(m))
        for p in eachringpixel(m, r)
            pixels[pixoff] = convert(ptype, p)
            pixphi[pixoff] = pixelphi(m, p)
            pixoff += 1
        end
        ringptr[i] = i==1 ? 1 : ringptr[i-1] + ringlength(m, r-1)
        ringtheta[i] = pixeltheta(m, first(eachringpixel(m, r)))
    end
    ringptr[nring+1] = npix + 1

    return RingPixelizedMap{ptype,Tp,Pc}(m, pixels, ringptr, ringtheta, pixphi)
end

