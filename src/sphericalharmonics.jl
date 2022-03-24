module SphericalHarmonics

export analyze, synthesize

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
struct RingInfo{R<:Real}
    """
    If `true`, the described ring has a partner ring reflected over the equator, otherwise
    `false`.
    """
    θpair::Bool

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

"""
    MapInfo{R<:Real}

Collection of [`RingInfo`](@ref) which describe the pixelization over the unit sphere.
"""
struct MapInfo{R<:Real}
    nθ::Int
    rings::Vector{RingInfo{R}}
end

function verify_mapinfo(mapinfo::MapInfo{R}; fullsky::Bool=true) where R
    nθ, npix, Ω = 0, 0, zero(R)
    allpairs = true
    for rinfo in mapinfo.rings
        allpairs &= rinfo.θpair
        if !allpairs && rinfo.θpair
            error("""
                  verification failed: encountered disjoint set of paired and unpaired rings;
                      expected all paired rings to be listed first, followed by trailing unpaired rings
                  """)
        end
        mult = rinfo.θpair ? 2 : 1
        nθ += mult
        npix += rinfo.nϕ
        Ω += rinfo.ΔΩ * rinfo.nϕ * mult
    end
    nθ == mapinfo.nθ || error("verification failed: counted $nθ rings, expected $(mapinfo.nθ)")
    if fullsky
        isapprox(Ω, 4R(π), atol = 0.5 * 4π / npix) ||
            error("verification failed: got $(round(Ω/π, sigdigits=4))π surface area, expected 4π")
    end
    return nothing
end


function _unsafe_accum_phase!(f₁, f₂, alms, Λ, rinfo::RingInfo)
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
        f₁[i+1] += a₁
        if f₂ !== nothing
            f₂[i+1] += a₂
        end
    end
    return nothing
end

function _unsafe_decomp_phase!(alms, f₁, f₂, Λ, rinfo::RingInfo)
    C = eltype(alms)
    lmax, mmax = size(alms) .- 1
    for m in 0:mmax
        i, isconj, isnyq = alias_index(rinfo.nϕ, m)
        a₁, a₂ = f₁[i+1], (f₂ !== nothing ? f₂[i+1] : zero(C))
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
    maprings = synthesize(info::MapInfo{R}, alms) where {R <: Real}
    synthesize!(maprings::Vector{Vector{R}}, info::MapInfo{R}, alms) where {R <: Real}

Perform the spherical harmonic synthesis transform from harmonic coefficients `alms` to
the map `mapbuf` (pixelized as described by `info`).
"""
function synthesize(info::MapInfo{R}, alms) where {R <: Real}
    nθ = info.nθ
    maprings = Vector{Vector{R}}(undef, nθ)
    for ii in 1:length(info.rings)
        nϕ = info.rings[ii].nϕ
        maprings[ii] = Vector{R}(undef, nϕ)
        if info.rings[ii].θpair
            maprings[nθ-ii+1] = Vector{R}(undef, nϕ)
        end
    end
    return synthesize!(maprings, info, alms)
end

@doc (@doc synthesize)
function synthesize!(maprings::Vector{Vector{R}}, info::MapInfo{R}, alms) where {R <: Real}
    C = eltype(alms)
    lmax, mmax = size(alms) .- 1

    nθ = info.nθ
    nϕ_max = mapreduce(r -> r.nϕ, max, info.rings)
    f₁ = zeros(C, nϕ_max)
    f₂ = zeros(C, nϕ_max)

    Λ = zeros(R, lmax + 1, mmax + 1)
    Λw = Legendre.Work(λlm!, Λ, Legendre.Scalar(zero(R)))

    for ii in 1:length(info.rings)
        jj = nθ - ii + 1
        ringinfo = info.rings[ii]
        unsafe_legendre!(Λw, Λ, lmax, mmax, ringinfo.cθ)

        nϕ = Int(ringinfo.nϕ)
        nϕh = nϕ ÷ 2 + 1
        f₁′ = @view f₁[1:nϕh]
        f₂′ = @view f₂[1:nϕh]
        F = plan_brfft(f₁′, nϕ)

        if ringinfo.θpair
            _unsafe_accum_phase!(f₁, f₂, alms, Λ, ringinfo)
        else
            _unsafe_accum_phase!(f₁, nothing, alms, Λ, ringinfo)
        end

        mul!(maprings[ii], F, f₁′)
        fill!(f₁′, zero(C))
        if ringinfo.θpair
            mul!(maprings[jj], F, f₂′)
            fill!(f₂′, zero(C))
        end
    end
    return maprings
end

"""
    alms = analyze(info::MapInfo{R}, maprings::Vector{Vector{R}}) where {R <: Real}
    analyze!(alms, info::MapInfo{R}, maprings::Vector{Vector{R}}) where {R <: Real}

Perform the spherical harmonic analysis transform from the map `maprings` (pixelized as
described by `info`) to its harmonic coefficients `alms`.
"""
function analyze(info::MapInfo{R}, maprings::Vector{Vector{R}}, lmax::Int, mmax::Int = lmax) where {R <: Real}
    alms = zeros(complex(R), lmax + 1, mmax + 1)
    return analyze!(alms, info, maprings)
end

@doc (@doc analyze)
function analyze!(alms, info::MapInfo{R}, maprings::Vector{Vector{R}}) where {R <: Real}
    C = eltype(alms)
    lmax, mmax = size(alms) .- 1

    nθ = info.nθ
    nϕ_max = mapreduce(r -> r.nϕ, max, info.rings)
    f₁ = Vector{C}(undef, nϕ_max)
    f₂ = Vector{C}(undef, nϕ_max)

    Λ = zeros(R, lmax + 1, mmax + 1)
    Λw = Legendre.Work(λlm!, Λ, Legendre.Scalar(zero(R)))

    for ii in 1:length(info.rings)
        jj = nθ - ii + 1
        ringinfo = info.rings[ii]
        unsafe_legendre!(Λw, Λ, lmax, mmax, ringinfo.cθ)

        nϕ = Int(ringinfo.nϕ)
        nϕh = nϕ ÷ 2 + 1
        f₁′ = @view f₁[1:nϕh]
        f₂′ = @view f₂[1:nϕh]
        F = plan_rfft(maprings[ii], 1)

        fill!(f₁′, zero(C))
        mul!(f₁′, F, maprings[ii])
        if ringinfo.θpair
            fill!(f₂′, zero(C))
            mul!(f₂′, F, maprings[jj])
            _unsafe_decomp_phase!(alms, f₁, f₂, Λ, ringinfo)
        else
            _unsafe_decomp_phase!(alms, f₁, nothing, Λ, ringinfo)
        end
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
        θpair = !(isodd(nθ) && ii == nθh)
        sθ, cθ = sincospi(R(2ii - 1) / 2nθ)
        rings[ii] = RingInfo{R}(θpair, nϕ, cθ, ϕ_π, sθ * ΔθΔϕ)
    end
    return MapInfo{R}(nθ, rings)
end

end
