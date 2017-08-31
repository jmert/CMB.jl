"""
A module of functions implementing function which interact with the HEALPix pixel
definitions. In most cases, only the RING ordering functions are being provided.

See "HEALPix: A Framework for High-Resolution Discretization and Fast Analysis of Data
Distributed on the Sphere" Górski, Hivon, & Banday et al (2005) ApJ 622:759–771
arXiv: astro-ph/0409513
"""
module Healpix

export
    npix2nside, nside2npix, nside2nring, nside2npixcap, nside2npixequ, nside2npixarea,
    isnorth, issouth, isnorthcap, issouthcap, iscap,
    isequbelt, isnorthequbelt, issouthequbelt,
    pix2ring, pix2ringidx, pix2z, pix2theta, pix2phi,
    UNSEEN

using MacroTools
using MacroTools: flatten, match, ismatch, postwalk, striplines

import Base: broadcast, broadcast!

"""
    const UNSEEN = -1.6375e+30

Special value recognized by the libhealpix/healpy routines as an unobserved/masked
pixel.
"""
const UNSEEN = -1.6375e+30

"""
    macro healpixkernel(funcexpr)

Aids in developing optimized kernel functions for the basic HEALPix routines. The
argument is expected to be a function expression of two arguments: the `nside` and a
pixel to operate upon `pix`.

```julia
pixelfunc(nside, pix) = ...
```

From the definition, a kernel function is automatically extracted with some
[single-argument] function calls involving `nside` lifted out of the kernel function and
replaced with variables; the variables are arguments to the kernel function, and the
outer function is responsible for computing the cached values.

This is all done so that optimized broadcast functions can also be automatically
generated where the cache variables are computed outside of the broadcast loop, and only
the kernel function is invoked on each iteration.
"""
macro healpixkernel(funcexpr)
    funcexpr = flatten(shortdef(striplines(funcexpr)))
    @capture(funcexpr, fn_(nside_, pix_) = body_)

    kfn = Symbol("_$(fn)_kernel")

    # List of special functions we know we can trivially transform
    kernmap = Symbol[
            :nside2npix, :nside2nring, :nside2npixcap,
            :nside2npixequ, :nside2pixarea
        ]

    # Mapping from the replacement expressions to the generated symbol
    cachemap = Dict{Any,Symbol}()

    # Identifies callable expressions where the function is in kernmap and the only
    # argument is exactly nside.
    isnside2fn(x) = begin
            isexpr(x, :call) || return false
            y = namify(x)
            y in kernmap || return false
            ismatch(:($y($nside)), x) || return false
            return true
        end

    # Permit arbitrary caching by parsing a "macro" of the form
    #   @cache(expr)
    # which is then lifted out of the kernel function.
    iscachemacro(x) = begin
            isexpr(x, :macrocall) || return false
            namify(x) == Symbol("@cache") || return false
            return true
        end

    makecache(x) = begin
            if isnside2fn(x)
                ex = x
            elseif iscachemacro(x)
                m = MacroTools.match(:(@cache(args_)), x)
                ex = m[:args]
            else
                # Break out early if we didn't encounter a cacheable state.
                return x
            end
            # At this point, ex is the expression that we actually want to
            # cache and use for later.
            if ex in keys(cachemap)
                sym = cachemap[x]
            else
                sym = gensym(namify(ex))
                cachemap[ex] = sym
            end
            return sym
        end

    cacheargs = Symbol[]
    body = MacroTools.postwalk(makecache, body)

    setup = Expr[:($v = $ex) for (ex,v) in cachemap]
    cacheargs = values(cachemap)

    quote
        @inline function $kfn($pix, $nside, $(cacheargs...))
            @fastmath begin
                $body
            end
        end

        function $fn($nside, $pix)
            $(setup...)
            return $kfn($pix, $nside, $(cacheargs...))
        end

        function Base.broadcast(::typeof($fn), $nside::Integer, $pix)
            $(setup...)
            return $kfn.($pix, $nside, $(cacheargs...))
        end

        function Base.broadcast!(::typeof($fn), A::AbstractArray, $nside::Integer, $pix)
            $(setup...)
            return A .= $kfn.($pix, $nside, $(cacheargs...))
        end

        Base.@__doc__ $fn
    end |> esc
end

Base.@pure npix2nside(npix)     = trunc(Int, sqrt(npix/12))
Base.@pure nside2npix(nside)    = 12*nside*nside
Base.@pure nside2nring(nside)   =  4*nside - 1
Base.@pure nside2npixcap(nside) =  2*nside*(nside - 1)
Base.@pure nside2npixequ(nside) =  2*nside*(3*nside + 1)
Base.@pure nside2pixarea(nside) = 4π / nside2npix(nside)

Base.@pure isnorthcap(nside, p) = p < nside2npixcap(nside)
Base.@pure issouthcap(nside, p) = p > nside2npix(nside) - nside2npixcap(nside) - 1
Base.@pure iscap(nside, p)      = isnorthcap(n, p) | issouthcap(n, p)

Base.@pure isnorthequbelt(nside, p) = (nside2npixcap(nside) ≤ p) & (p < nside2npixequ(nside))
Base.@pure issouthequbelt(nside, p) = (nside2npixequ(nside) ≤ p) & (p ≤ nside - nside2npixcap(nside) - 1)
Base.@pure isequbelt(nside, p)      = ~iscap(n, p)

Base.@pure isnorth(r, p) = p < nside2npixequ(r)
Base.@pure issouth(r, p) = p ≥ nside2npixequ(r)

# The function definitions are short enough that I wanted to keep them visually together
# without having the documentation break them up. Add documentation now

"""
    npix2nside(npix) -> Nside

Returns the equivalent `Nside` corresponding to the number of pixels `npix`. Note that no
validation is performed, so non-conformant values of `npix` will give non-conformant
Nside values.
""" npix2nside(npix)

"""
    nside2npix(nside) -> N

Returns the total number of pixels `N` in an `nside` HEALPix map.

N.B.: HEALPix pixel indexing is 0-based, so valid pixel values are in the range
0:(N-1).
""" nside2npix

"""
    nside2nring(nside) -> N

Returns the number of iso-latitude rings `N` in the `nside` HEALPix map.
""" nside2nring

"""
    nside2npixcap(nside) -> N

Returns the number of pixels `N` in the polar caps for the given `nside` HEALPix map.
""" nside2npixcap

"""
    nside2npixequ(nside) -> N

Returns the number of pixels `N` in the northern hemisphere, including the equatorial
ring, for the given `nside` HEALPix map.
""" nside2npixequ

"""
    nside2pixarea(nside) -> σ

Returns the surface area `σ` (in steradiands) of each pixel in the given `nside` HEALPix
map.
""" nside2pixarea

"""
    isnorthcap(nside, pix) -> TF

Test for whether the given pixel `pix` is in the northern polar cap for an `nside`
HEALPix map.
""" isnorthcap

"""
    issouthcap(nside, pix) -> TF

Test for whether the given pixel `pix` is in the southern polar cap for an `nside`
HEALPix map.
""" issouthcap

"""
    iscap(nside, pix) -> TF

Test for whether the given pixel `pix` is in either polar cap for an `nside` HEALPix
map.
""" iscap

"""
    isnorthequbelt(nside, pix) -> TF

Test for whether the given pixel `pix` is in the northern equatorial belt (including the
equatorial ring) for an `nside` HEALPix map.
""" isnorthequbelt

"""
    issouthequbelt(nside, pix) -> TF

Test for whether the given pixel `pix` is in the southern equatorial belt (excluding the
equatorial ring) for an `nside` HEALPix map.
""" issouthequbelt

"""
    isequbelt(nside, pix) -> TF

Test for whether the given pixel `pix` is in the equatorial belt for an `nside` HEALPix
map.
""" isequbelt

"""
    isnorth(nside, pix) -> TF

Test for whether the given pixel `pix` is in the northern hemisphere (including the
equatorial ring) for an `nside` HEALPix map.
""" isnorth

"""
    issouth(nside, pix) -> TF

Test for whether the given pixel `pix` is in the southern hemisphere (excluding the
equatorial ring) for an `nside` HEALPix map.
""" issouth

####

"""
    pix2ring(nside, p) -> i

Computes the ring index `i` for the given pixel `p`. `nside` is the Nside resolution
factor.
"""
@healpixkernel function pix2ring(nside, p)
    p′ = p < nside2npixequ(nside) ? p : (nside2npix(nside)-1) - p
    if p′ < nside2npixcap(nside)
        i′ = trunc(Int, sqrt((p′+1)/2 - sqrt(convert(Float64,(p′+1)>>1)))) + 1
    else
        i′ = div(p′-nside2npixcap(nside), 4nside) + nside
    end
    i = p < nside2npixequ(nside) ? i′ : (nside2nring(nside)+1) - i′
    return i
end

"""
    pix2ringidx(nside, p) -> j

Computes the index `j` within the ring for the given pixel `p`. `nside` is the Nside
resolution factor.
"""
@healpixkernel function pix2ringidx(nside, p)
    p′ = p < nside2npixequ(nside) ? p : (nside2npix(nside)-1) - p
    if p′ < nside2npixcap(nside)
        i′ = trunc(Int, sqrt((p′+1)/2 - sqrt(convert(Float64,(p′+1)>>1)))) + 1
        j′ = p′ + 1 - nside2npixcap(i′)
        l′ = 4i′
    else
        j′ = rem(p′-nside2npixcap(nside), 4nside) + 1
        l′ = 4nside
    end
    j = p < nside2npixequ(nside) ? j′ : (l′ - j′ + 1)
    return j
end

"""
    pix2z(nside, p) -> z

Computes the cosine of the colatitude `z` for the given pixel `p`. `nside` is the Nside
resolution factor.
"""
@healpixkernel function pix2z(nside, p)
    p′ = p < nside2npixequ(nside) ? p : (nside2npix(nside)-1) - p
    if p′ < nside2npixcap(nside)
        i′ = trunc(Int, sqrt((p′+1)/2 - sqrt(convert(Float64,(p′+1)>>1)))) + 1
        z′ = 1.0 - i′^2 / @cache(3nside^2)
    else
        i′ = div(p′-nside2npixcap(nside), 4nside) + nside
        z′ = 4/3 - 2i′ / @cache(3nside)
    end
    z = p < nside2npixequ(nside) ? z′ : -z′
    return z
end

"""
    pix2theta(nside, p) -> θ

Computes the colatitude `θ` for the given pixel `p`. `nside` is the Nside resolution
factor.
"""
@healpixkernel function pix2theta(nside, p)
    p′ = p < nside2npixequ(nside) ? p : (nside2npix(nside)-1) - p
    if p′ < nside2npixcap(nside)
        i′ = trunc(Int, sqrt((p′+1)/2 - sqrt(convert(Float64,(p′+1)>>1)))) + 1
        z′ = 1.0 - i′^2 / @cache(3nside^2)
    else
        i′ = div(p′-nside2npixcap(nside), 4nside) + nside
        z′ = 4/3 - 2i′ / @cache(3nside)
    end
    z = p < nside2npixequ(nside) ? z′ : -z′
    return acos(z)
end

"""
    pix2phi(nside, p) -> ϕ

Computes the azimuth `ϕ` for the given pixel `p`. `nside` is the Nside resolution
factor.
"""
@healpixkernel function pix2phi(nside, p)
    p′ = p < nside2npix(nside)-nside2npixcap(nside) ? p : (nside2npix(nside)-1) - p
    if p′ < nside2npixcap(nside)
        i′ = trunc(Int, sqrt((p′+1)/2 - sqrt(convert(Float64,(p′+1)>>1)))) + 1
        j′ = p′ + 1 - nside2npixcap(i′)
        ϕ′ = (π/2)/i′ * (j′ - 0.5)
    else
        i′,j′ = divrem(p′-nside2npixcap(nside), 4nside)
        i′ += nside
        j′ += 1
        ϕ′ = @cache(π/2/nside) * (j′ - 0.5*(1 + rem(i′-nside,2)))
    end
    ϕ = p < nside2npix(nside)-nside2npixcap(nside) ? ϕ′ : (2π - ϕ′)
    return ϕ
end

end # module Healpix

