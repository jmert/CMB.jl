#
# See "HEALPix: A Framework for High-Resolution Discretization and Fast Analysis of Data
# Distributed on the Sphere" Górski, Hivon, & Banday et al (2005) ApJ 622:759–771
# arXiv: astro-ph/0409513
#
module Healpix
    export
        nside2npix, nside2nring,
        pix2ring_ring, pix2ringidx_ring, pix2z_ring, pix2theta_ring, pix2phi_ring

    doc"""
        nside2npix(nside)

    Returns the number of pixels in an ``N_\mathrm{side} = ```nside` HEALPix map.
    """
    @inline function nside2npix(nside::N) where N<:Integer
        @boundscheck (ispow2(nside) || throw(DomainError()))
        return 12*nside*nside
    end

    doc"""
        nside2nring(nside)

    Returns the number of isolatitude rings in an ``N_\mathrm{side} = ```nside` HEALPix map.
    """
    @inline function nside2nring(nside::N) where N<:Integer
        @boundscheck (ispow2(nside) || throw(DomainError()))
        return 4*nside - 1
    end

    doc"""
        pixel_region_ring(pix, nisde, northcap, equatorial, southcap)

    Based on the ring-ordered pixel `p` and ``N_\mathrm{side}=```nside`, branches to select
    one of the expressions `northcap`, `equatorial`, and `southcap` for each of the northern
    cap, equatorial belt, and southern cap regions of a HEALPix map.

    The comparisons expand to an `ifelse` chain with comparison values interpolated in as
    constants. The use of `ifelse` provides branchless code at the cost of unconditionally
    evaluating all three expressions before selecting the value to return.
    """
    macro pixel_region_ring(pix, nside, northcap, equatorial, southcap)
        p = esc(pix)
        n = esc(nside)
        return quote
            north_cap_done  =  2*$n^2 - 2*$n
            equatorial_done = 10*$n^2 + 2*$n
            ifelse($p < north_cap_done, $(esc(northcap)),
                ifelse($p < equatorial_done, $(esc(equatorial)),
                    $(esc(southcap))
                )
            )
        end
    end

    # Getting the southern cap pixels referenced as rings from the South Pole is useful
    # in multiple places, so this internal function returns rings in that form.
    #
    # Keep this in sync with pix2ring_ring below (modulo the update to North Pole-based
    # numbering).
    @inline function _pix2ring_ring(nside, p)
        @inbounds npix = nside2npix(nside)
        @inbounds nrings = nside2nring(nside)
        @fastmath begin
            p′ = @pixel_region_ring(p, nside, p, p, (npix-1)-p)
            # We can reuse the same computation with two separate interpretations for the
            # polar caps. For the northern cap, this counts rings from the North Pole,
            # while for the southern cap, from the South Pole.
            i′ = trunc(Int, sqrt((p′+1)/2 - sqrt(convert(Float64,(p′+1)>>1)))) + 1
            i = @pixel_region_ring(p, nside,
                    i′,
                    fld(p-2nside*(nside-1), 4nside) + nside,
                    i′
                   )
        end
        return i
    end

    doc"""
        pix2ring_ring(nside, p)

    For a ring-ordered pixel `p` in an ``N_\mathrm{side}=```nside` map, returns the index of
    the isolatitude ring containing the pixel in the range 1 to
    [`nside2nring(nside)`](@ref nside2nring).
    """
    @inline function pix2ring_ring(nside, p)
        @inbounds npix = nside2npix(nside)
        @inbounds nrings = nside2nring(nside)
        @fastmath begin
            p′ = @pixel_region_ring(p, nside, p, p, npix-1-p)
            # We can reuse the same computation with two separate interpretations for the
            # polar caps. For the northern cap, this counts rings from the North Pole,
            # while for the southern cap, from the South Pole.
            i′ = trunc(Int, sqrt((p′+1)/2 - sqrt(convert(Float64,(p′+1)>>1)))) + 1
            i = @pixel_region_ring(p, nside,
                    i′,
                    fld(p-2nside*(nside-1), 4nside) + nside,
                    (nrings+1) - i′
                   )
        end
        return i
    end

    """
        pix2ringidx_ring(nside, p)
    """
    @inline function pix2ringidx_ring(nside, p)
        @inbounds npix = nside2npix(nside)
        @fastmath begin
            p′ = @pixel_region_ring(p, nside, p, p, npix-1-p)
            # Only used for the northern and southern caps, but with p′ reflected, we can reuse
            # the same computation with two separate interpretations. For the northern cap, this
            # counts rings from the North Pole, while for the southern cap, from the South Pole.
            i = trunc(Int, sqrt((p′+1)/2 - sqrt(convert(Float64,(p′+1)>>1)))) + 1
            j = @pixel_region_ring(p, nside,
                    p′ + 1 - 2i*(i-1),
                    rem(p′ - 2nside*(nside-1), 4nside) + 1,
                    1 + 4i - (p′ + 1 - 2i*(i-1))
                   )
        end
        return j
    end

    """
        pix2z_ring(nside, p)
    """
    @inline function pix2z_ring(nside, p)
        i = _pix2ring_ring(nside, p)
        @fastmath begin
            z = @pixel_region_ring(p, nside,
                        1.0 - i^2/(3nside^2),
                        4/3 - 2i/(3nside),
                       -1.0 + i^2/(3nside^2)
                      )
        end
        return z
    end

    """
        pix2theta_ring(nside, p)
    """
    @inline pix2theta_ring(nside, p) = acos(pix2z_ring(nside, p))

    """
        pix2phi_ring(nside, p)
    """
    @inline function pix2phi_ring(nside, p)
        @inbounds nring = nside2nring(nside)
        i = _pix2ring_ring(nside, p)
        j = pix2ringidx_ring(nside, p)
        @fastmath begin
            ϕ = @pixel_region_ring(p, nside,
                    (pi/2)/i * (j - 0.5),
                    (pi/2/nside) * (j - 0.5(1 + rem(i-nside,2))),
                    (pi/2)/i * (j - 0.5)
                   )
        end
        return ϕ
    end
end
