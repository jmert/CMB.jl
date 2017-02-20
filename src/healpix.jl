module Healpix
    export
        nside2npix, nside2nring,
        pix2ring_ring

    @inline function nside2npix(nside::N) where N<:Integer
        @boundscheck (ispow2(nside) || throw(DomainError()))
        return 12*nside*nside
    end

    @inline function nside2nring(nside::N) where N<:Integer
        @boundscheck (ispow2(nside) || throw(DomainError()))
        return 4*nside - 1
    end

    @generated function _pix2ring_ring_north{nside}(::Type{Val{nside}}, p)
        north_cap_done =  2nside^2 - 2nside
        return quote
            $(Expr(:meta, :inline))
            r = ifelse(p < $north_cap_done,
                       trunc(Int,sqrt((p+1)/2 - sqrt((p+1)>>1))) + 1,
                       div(p - $north_cap_done, $(4nside)) + $nside
                      )
            return r
        end
    end
    @generated function _pix2ring_ring{nside}(::Type{Val{nside}}, p)
        npix   = 12nside^2
        nrings = 4nside - 1
        north_eq_done  =  6nside^2 + 2nside
        return quote
            $(Expr(:meta, :inline))
            isnorth = p < $north_eq_done
            p = ifelse(isnorth, p, $(npix-1) - p)
            r = _pix2ring_ring_north(Val{$nside}, p)
            return ifelse(isnorth, r, $(nrings+1) - r)
        end
    end
    @inline pix2ring_ring(nside::I, p) where I<:Integer = _pix2ring_ring(Val{nside}, p)
end
