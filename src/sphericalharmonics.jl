module SphericalHarmonics

export analyze, analyze!, synthesize, synthesize!

import ..Legendre
import ..Legendre: λlm!, unsafe_legendre!

using FFTW
using LinearAlgebra: mul!

@static if VERSION < v"1.6.0-DEV.1591"
    using Compat: sincospi # Compat@v3.23
    using Compat: cispi # Compat@v3.25
end

"""
    alias_index(len::Int, m::Int) -> (i::Int, isconj::Bool, isnyq::Bool)

Aliases the given 0-indexed `m` for a periodic FFT frequency axis of length `len` --- with
an additional aliasing to below the Nyquist frequency approprate for use in a real-only
[`irfft()`] transform --- to the new index `i`.

To handle symmetries of the FFT, the flags `isconj` and `isnyq` indicate if the Fourier
coefficients at `i` should be conjugated or real-doubled (i.e. `c = complex(2real(c))`),
respectively.
"""
@inline function alias_index(len::Int, i::Int)
    isconj, isnyq = false, false
    nyq = max(1, len ÷ 2)
    i < nyq && return (i, isconj, isnyq)
    i = mod(i, len)
    if i > nyq
        i = len - i
        isconj = true
    elseif i == 0 || (iseven(len) && i == nyq)
        isnyq = true
    end
    return (i, isconj, isnyq)
end

"""
    alias_coeffs(coeffs, isconj::Bool, isnyq::Bool) -> coeffs

Helper for [`alias_index`](@ref) which applies the conjugation or pure-real symmetry
transformations to the Fourier coefficients `coeffs` based on the two flags `isconj` and
`isnyq`.
"""
@inline function alias_coeffs(coeffs, isconj::Bool, isnyq::Bool)
    return ifelse(isnyq, complex.(2 .* real.(coeffs)),
                  ifelse(isconj, conj.(coeffs), coeffs))
end


"""
    RingInfo{R<:Real}

Parameters required to describe an isolatitude ring on the sphere compatible with spherical
harmonic transforms.
"""
struct RingInfo{R}
    """
    Offsets (1-indexed) to the first pixel of the isolatitude ring(s) within a
    map array. The first offset in the tuple corresponds to the colatitude ``θ``
    given by `.cθ`, while the second offset is for the ring mirrored over the
    equator (colatitude ``π - θ``). Either offset may be 0 to indicate the ring
    is not present.
    """
    offset::Tuple{Int,Int}

    """
    Stride between pixels within the ring.
    """
    stride::Int

    """
    The number of equispaced pixels in the isolatitude ring.
    """
    nϕ::Int

    """
    Cosine of the colatitude of the ring, i.e. ``cos(θ)``.
    """
    cθ::R

    """
    Azimuth of the first pixel in the isolatitude ring ``ϕ₀`` divided by ``π``.
    """
    ϕ_π::R

    """
    Solid angle area of pixels in the ring, i.e. ``sin(θ) × Δθ × Δϕ``.
    """
    ΔΩ::R
end

abstract type AbstractMapInfo{R} end

"""
    rings(mapinfo::AbstractMapInfo)
"""
function rings end

function mapinfo_counts(mapinfo::AbstractMapInfo)
    nr, np, nϕ = 0, 0, typemin(Int)
    for rinfo in rings(mapinfo)
        o₁, o₂ = rinfo.offset
        mult = (o₁ != 0) + (o₂ != 0)
        nr += mult
        np += mult * rinfo.nϕ
        nϕ = mult != 0 ? max(nϕ, rinfo.nϕ) : nϕ
    end
    return (; nring = nr, npix = np, nϕmax = nϕ)
end

"""
    nring(mapinfo::AbstractMapInfo)

Number of isolatitude rings described by `mapinfo`.
"""
nring(mapinfo::AbstractMapInfo) = mapinfo_counts(mapinfo)[:nring]

"""
    npix(mapinfo::AbstractMapInfo)

Number of pixels covered by all rings described by `mapinfo`.
"""
npix(mapinfo::AbstractMapInfo) = mapinfo_counts(mapinfo)[:npix]

"""
    nϕmax(mapinfo::AbstractMapInfo)

Number of pixels in the longest isolatitude ring described by `mapinfo`.
"""
nϕmax(mapinfo::AbstractMapInfo) = mapinfo_counts(mapinfo)[:nϕmax]

function verify_mapinfo(mapinfo::AbstractMapInfo{R}; fullsky::Bool=true) where R
    Ω = zero(R)
    pix = BitSet()
    for rinfo in rings(mapinfo)
        mult = 0
        for off in rinfo.offset
            iszero(off) && continue
            mult += 1
            ringpix = (off - 1) .+ range(1, step = rinfo.stride, length = rinfo.nϕ)
            if !isempty(intersect(ringpix, pix))
                error("verification failed: overlapping range of pixel rings encountered")
            end
            union!(pix, ringpix)
        end
        Ω += rinfo.ΔΩ * rinfo.nϕ * mult
    end
    if fullsky
        npix = length(pix)
        if pix != BitSet(1:npix)
            error("verification failed: a sparse pixel layout (non-consecutive indices) has been detected")
        end
    #    isapprox(Ω, 4R(π), atol = 0.5 * 4R(π) / npix) ||
    #        error("verification failed: got $(round(Ω/π, sigdigits=4))π surface area, expected 4π")
    end
    return nothing
end


"""
    MapInfo{R} <: AbstractMapInfo{R}

An arbitrary vector of [`RingInfo`](@ref) which describe the pixelization over the unit
sphere.
"""
struct MapInfo{R} <: AbstractMapInfo{R}
    rings::Vector{RingInfo{R}}
end

rings(mapinfo::MapInfo) = mapinfo.rings


Base.@propagate_inbounds function _unsafe_accum_phase!(f₁, f₂, alms, Λ, rinfo::RingInfo)
    C = eltype(alms)
    lmax, mmax = size(alms) .- 1
    for m in 0:mmax
        a₁, a₂ = zero(C), zero(C)
        for ℓ in m:lmax
            term = alms[ℓ+1,m+1] * Λ[ℓ+1,m+1]
            a₁ += term
            a₂ += isodd(ℓ + m) ? -term : term
        end
        i, isconj, isnyq = alias_index(rinfo.nϕ, m)
        a₁, a₂ = (a₁, a₂) .* cispi(m * rinfo.ϕ_π)
        a₁, a₂ = alias_coeffs((a₁, a₂), isconj, isnyq)
        if f₁ !== nothing
            f₁[i+1] += a₁
        end
        if f₂ !== nothing
            f₂[i+1] += a₂
        end
    end
    return nothing
end

Base.@propagate_inbounds function _unsafe_decomp_phase!(alms, f₁, f₂, Λ, rinfo::RingInfo)
    C = eltype(alms)
    lmax, mmax = size(alms) .- 1
    for m in 0:mmax
        i, isconj, isnyq = alias_index(rinfo.nϕ, m)
        a₁ = f₁ !== nothing ? f₁[i+1] : zero(C)
        a₂ = f₂ !== nothing ? f₂[i+1] : zero(C)
        a₁, a₂ = alias_coeffs((a₁, a₂), isconj, isnyq)
        a₁, a₂ = (a₁, a₂) .* (rinfo.ΔΩ * cispi(m * -rinfo.ϕ_π))
        for ℓ in m:lmax
            c = isodd(ℓ+m) ? a₁ - a₂ : a₁ + a₂
            alms[ℓ+1,m+1] += c * Λ[ℓ+1,m+1]
        end
    end
    return nothing
end


"""
    mapbuf = synthesize(info::AbstractMapInfo{R}, alms) where {R <: Real}
    synthesize!(mapbuf::AbstractArray{R}, info::AbstractMapInfo{R}, alms) where {R <: Real}

Perform the spherical harmonic synthesis transform from harmonic coefficients `alms` to
the map `mapbuf` (pixelized as described by `info`).
"""
function synthesize(info::AbstractMapInfo{R}, alms) where {R <: Real}
    mapbuf = Vector{R}(undef, npix(info))
    return synthesize!(mapbuf, info, alms)
end

@doc (@doc synthesize)
function synthesize!(mapbuf::AbstractArray{R}, info::AbstractMapInfo{R}, alms) where {R <: Real}
    C = eltype(alms)
    lmax, mmax = size(alms) .- 1

    nϕ_max = nϕmax(info)
    rv = Vector{R}(undef, nϕ_max)
    f₁ = zeros(C, (nϕ_max ÷ 2) + 1)
    f₂ = zeros(C, (nϕ_max ÷ 2) + 1)

    Λ = zeros(R, lmax + 1, mmax + 1)
    Λw = Legendre.Work(λlm!, Λ, Legendre.Scalar(zero(R)))

    for rinfo in rings(info)
        unsafe_legendre!(Λw, Λ, lmax, mmax, rinfo.cθ)

        o₁, o₂ = rinfo.offset
        nϕ = rinfo.nϕ
        nϕh = nϕ ÷ 2 + 1
        f₁′ = @view f₁[1:nϕh]
        f₂′ = @view f₂[1:nϕh]
        F = plan_brfft(f₁′, nϕ)

        _unsafe_accum_phase!(o₁ != 0 ? f₁′ : nothing,
                             o₂ != 0 ? f₂′ : nothing,
                             alms, Λ, rinfo)

        to_nϕ = Base.OneTo(nϕ)
        @inline _strided(o) = range(o, step = rinfo.stride, length = nϕ)

        rv′ = @view rv[1:nϕ]
        if o₁ != 0
            mul!(rv′, F, f₁′)
            copyto!(mapbuf, _strided(o₁), rv, to_nϕ)
            fill!(f₁′, zero(C))
        end
        if o₂ != 0
            mul!(rv′, F, f₂′)
            copyto!(mapbuf, _strided(o₂), rv, to_nϕ)
            fill!(f₂′, zero(C))
        end
    end
    return mapbuf
end

"""
    alms = analyze(info::AbstractMapInfo{R}, mapbuf::AbstractArray{R}) where {R <: Real}
    analyze!(alms, info::AbstractMapInfo{R}, mapbuf::AbstractArray{R}) where {R <: Real}

Perform the spherical harmonic analysis transform from the map `mapbuf` (pixelized as
described by `info`) to its harmonic coefficients `alms`.
"""
function analyze(info::AbstractMapInfo{R}, mapbuf::AbstractArray{R}, lmax::Int, mmax::Int = lmax) where {R <: Real}
    alms = zeros(complex(R), lmax + 1, mmax + 1)
    return analyze!(alms, info, mapbuf)
end

@doc (@doc analyze)
function analyze!(alms, info::AbstractMapInfo{R}, mapbuf::AbstractArray{R}) where {R <: Real}
    C = eltype(alms)
    lmax, mmax = size(alms) .- 1

    nϕ_max = nϕmax(info)
    rv = Vector{R}(undef, nϕ_max)
    f₁ = Vector{C}(undef, (nϕ_max ÷ 2) + 1)
    f₂ = Vector{C}(undef, (nϕ_max ÷ 2) + 1)

    Λ = zeros(R, lmax + 1, mmax + 1)
    Λw = Legendre.Work(λlm!, Λ, Legendre.Scalar(zero(R)))

    for rinfo in rings(info)
        unsafe_legendre!(Λw, Λ, lmax, mmax, rinfo.cθ)

        o₁, o₂ = rinfo.offset
        nϕ = rinfo.nϕ
        nϕh = nϕ ÷ 2 + 1
        rv′ = @view rv[1:nϕ]
        F = plan_rfft(rv′, 1)

        to_nϕ = Base.OneTo(nϕ)
        @inline _strided(o) = range(o, step = rinfo.stride, length = nϕ)

        f₁′ = @view f₁[1:nϕh]
        f₂′ = @view f₂[1:nϕh]
        if o₁ != 0
            copyto!(rv, to_nϕ, mapbuf, _strided(o₁))
            mul!(f₁′, F, rv′)
        end
        if o₂ != 0
            copyto!(rv, to_nϕ, mapbuf, _strided(o₂))
            mul!(f₂′, F, rv′)
        end
        _unsafe_decomp_phase!(alms,
                              o₁ != 0 ? f₁′ : nothing,
                              o₂ != 0 ? f₂′ : nothing,
                              Λ, rinfo)
    end
    return alms
end

function ecp_mapinfo(::Type{R}, nθ::Int, nϕ::Int) where {R <: Real}
    # rings are symmetric over the equator, so only require half
    nθh = (nθ + 1) ÷ 2
    rings = Vector{RingInfo{R}}(undef, nθh)
    # ϕ offset and ΔθΔϕ are same for all rings
    ϕ_π = one(R) / nϕ
    ΔθΔϕ = 2R(π)^2 / (nθ * nϕ)
    for ii in 1:nθh
        o₁ = ii
        o₂ = nθ - ii + 1
        offs = (o₁, isodd(nθ) && ii == nθh ? 0 : o₂)
        sθ, cθ = sincospi(R(2ii - 1) / 2nθ)
        rings[ii] = RingInfo{R}(offs, nθ, nϕ, cθ, ϕ_π, sθ * ΔθΔϕ)
    end
    return MapInfo{R}(rings)
end

end
