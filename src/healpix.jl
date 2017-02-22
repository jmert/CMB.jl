#
# See "HEALPix: A Framework for High-Resolution Discretization and Fast Analysis of Data
# Distributed on the Sphere" Górski, Hivon, & Banday et al (2005) ApJ 622:759–771
# arXiv: astro-ph/0409513
#
module Healpix
    export
        nside2npix, nside2nring,
        pix2ring_ring, pix2ringidx_ring, pix2z_ring

    @inline function nside2npix(nside::N) where N<:Integer
        @boundscheck (ispow2(nside) || throw(DomainError()))
        return 12*nside*nside
    end

    @inline function nside2nring(nside::N) where N<:Integer
        @boundscheck (ispow2(nside) || throw(DomainError()))
        return 4*nside - 1
    end

    # The HEALPix grid is divided into three main sections:
    #   1. North polar cap
    #   2. Equatorial belt
    #   3. South polar cap
    # Since checking whether a pixel exist in each of these regions is going to be a
    # common activity, it'll be convenient to have a macro which generates the correct
    # three-way conditionals.
    macro pixel_region_ring(pix, nside, northcap, equatorial, southcap)
        north_cap_done  =  2nside^2 - 2nside
        equatorial_done = 10nside^2 + 2nside
        return :(
                 ifelse($(esc(pix)) < $north_cap_done, $(esc(northcap)),
                        ifelse($(esc(pix)) < $equatorial_done, $(esc(equatorial)),
                               $(esc(southcap))
                              )
                       )
                )
    end

    # Getting the southern cap pixels referenced as rings from the South Pole is useful
    # in multiple places, so this internal function returns rings in that form.
    #
    # Keep this in sync with pix2ring_ring below (modulo the update to North Pole-based
    # numbering).
    @generated function _pix2ring_ring{nside}(::Type{Val{nside}}, p)
        npix = nside2npix(nside)
        nrings = nside2nring(nside)
        fsqrt = @fastmath sqrt
        return quote
            $(Expr(:meta, :inline))
            p′ = @pixel_region_ring(p, $nside, p, p, $(npix-1)-p)
            # We can reuse the same computation with two separate interpretations for the
            # polar caps. For the northern cap, this counts rings from the North Pole,
            # while for the southern cap, from the South Pole.
            i′ = trunc(Int, $fsqrt((p′+1)/2 - $fsqrt(convert(Float64,(p′+1)>>1)))) + 1
            i = @pixel_region_ring(p, $nside,
                    i′,
                    div(p - $(2nside*(nside-1)), $(4nside)) + $nside,
                    i′
                   )
            return i
        end
    end

    @generated function pix2ring_ring{nside}(::Type{Val{nside}}, p)
        npix = nside2npix(nside)
        nrings = nside2nring(nside)
        fsqrt = @fastmath sqrt
        return quote
            $(Expr(:meta, :inline))
            p′ = @pixel_region_ring(p, $nside, p, p, $(npix-1)-p)
            # We can reuse the same computation with two separate interpretations for the
            # polar caps. For the northern cap, this counts rings from the North Pole,
            # while for the southern cap, from the South Pole.
            i′ = trunc(Int, $fsqrt((p′+1)/2 - $fsqrt(convert(Float64,(p′+1)>>1)))) + 1
            i = @pixel_region_ring(p, $nside,
                    i′,
                    div(p - $(2nside*(nside-1)), $(4nside)) + $nside,
                    $(nrings+1) - i′
                   )
            return i
        end
    end
    @inline pix2ring_ring(nside::I, p) where I<:Integer = pix2ring_ring(Val{nside}, p)

    @generated function pix2ringidx_ring{nside}(::Type{Val{nside}}, p)
        npix = nside2npix(nside)
        fsqrt = @fastmath sqrt
        return quote
            $(Expr(:meta, :inline))
            p′ = @pixel_region_ring(p, $nside, p, p, $(npix-1)-p)
            # Only used for the northern and southern caps, but with p′ reflected, we can
            # reuse the same computation with two separate interpretations. For the
            # northern cap, this counts rings from the North Pole, while for the southern
            # cap, from the South Pole.
            i = trunc(Int, $fsqrt((p′+1)/2 - $fsqrt(convert(Float64,(p′+1)>>1)))) + 1
            j = @pixel_region_ring(p, $nside,
                    p′ + 1 - 2i*(i-1),
                    mod(p′ - $(2nside*(nside-1)), $(4nside)) + 1,
                    1 + 4i - (p′ + 1 - 2i*(i-1))
                   )
            return j
        end
    end
    @inline pix2ringidx_ring(nside::I, p) where I<:Integer =
        pix2ringidx_ring(Val{nside}, p)

    @generated function pix2z_ring{nside}(::Type{Val{nside}}, p)
        north_cap_done = 2nside^2 - 2nside
        return quote
            $(Expr(:meta, :inline))
            i = _pix2ring_ring(Val{$nside}, p)
            z = @pixel_region_ring(p, $nside,
                        1.0 - i^2/$(3nside^2),
                        4/3 - 2i/$(3nside),
                       -1.0 + i^2/$(3nside^2)
                      )
            return z
        end
    end
    @inline pix2z_ring(nside::I, p) where I<:Integer = pix2z_ring(Val{nside}, p)
end
