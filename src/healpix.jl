#
# See "HEALPix: A Framework for High-Resolution Discretization and Fast Analysis of Data
# Distributed on the Sphere" Górski, Hivon, & Banday et al (2005) ApJ 622:759–771
# arXiv: astro-ph/0409513
#
module Healpix
    export
        nside2npix, nside2nring, nside2npixcap, nside2npixequ, nside2npixarea,
        isnorth, issouth, isnorthcap, issouthcap, iscap,
        isequbelt, isnorthequbelt, issouthequbelt,
        pix2ring, pix2ringidx, pix2z, pix2theta, pix2phi

    using MacroTools
    using MacroTools: flatten, ismatch, postwalk, striplines

    import Base: broadcast, broadcast!

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

        kernmap = Symbol[:nside2npix, :nside2nring, :nside2npixcap,
                :nside2npixequ, :nside2pixarea]
        iscacheable(x) = begin
                isexpr(x, :call) || return false
                y = namify(x)
                y in kernmap || return false
                ismatch(:($y($nside)), x)
            end
        getcachename(x) = Symbol(replace(string(namify(x)), "nside2", ""))

        cacheargs = Symbol[]
        body = MacroTools.postwalk(
            x -> ~iscacheable(x) ? x : begin
                y = getcachename(x)
                y in cacheargs || push!(cacheargs, y)
                y
            end, body)

        setup = Expr[:($n = $(Symbol("nside2$n"))($nside)) for n in cacheargs]

        @eval begin
            @inline function $kfn($pix, $nside, $(cacheargs...))
                @fastmath begin
                    $body
                end
            end

            function $fn($nside, $pix)
                $(setup...)
                return $kfn($pix, $nside, $(cacheargs...))
            end

            function Base.broadcast(::typeof($fn), $nside, $pix)
                $(setup...)
                return $kfn.($pix, $nside, $(cacheargs...))
            end

            function Base.broadcast!(::typeof($fn), A::AbstractArray, $nside, $pix)
                $(setup...)
                return A .= $kfn.($pix, $nside, $(cacheargs...))
            end

            Base.@__doc__ $fn
        end
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

    @healpixkernel function pix2z(nside, p)
        p′ = p < nside2npixequ(nside) ? p : (nside2npix(nside)-1) - p
        if p′ < nside2npixcap(nside)
            i′ = trunc(Int, sqrt((p′+1)/2 - sqrt(convert(Float64,(p′+1)>>1)))) + 1
            z′ = 1.0 - i′^2 / (3nside^2)
        else
            i′ = div(p′-nside2npixcap(nside), 4nside) + nside
            z′ = 4/3 - 2i′ / (3nside)
        end
        z = p < nside2npixequ(nside) ? z′ : -z′
        return z
    end

    @healpixkernel function pix2theta(nside, p)
        p′ = p < nside2npixequ(nside) ? p : (nside2npix(nside)-1) - p
        if p′ < nside2npixcap(nside)
            i′ = trunc(Int, sqrt((p′+1)/2 - sqrt(convert(Float64,(p′+1)>>1)))) + 1
            z′ = 1.0 - i′^2 / (3nside^2)
        else
            i′ = div(p′-nside2npixcap(nside), 4nside) + nside
            z′ = 4/3 - 2i′ / (3nside)
        end
        z = p < nside2npixequ(nside) ? z′ : -z′
        return acos(z)
    end

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
            ϕ′ = (π/2/nside) * (j′ - 0.5*(1 + rem(i′-nside,2)))
        end
        ϕ = p < nside2npix(nside)-nside2npixcap(nside) ? ϕ′ : (2π - ϕ′)
        return ϕ
    end
end

