module SphericalHarmonics
    import Compat.@compat # Compat@v3.21 for @compat import Mod as NewName
    using ..Healpix

    # reexport from SphericalHarmonicTransforms
    @compat import SphericalHarmonicTransforms as SHT
    using .SHT: analyze, analyze!, synthesize, synthesize!
    export analyze, analyze!, synthesize, synthesize!

    # imports that are part of the unexported interface
    import .SHT: AbstractRingPixelization, RingPixelization, Ring,
                 ECPPixelization, rings, nring, npix, nϕmax

    @doc raw"""
        HealpixPixelization{R} <: AbstractRingPixelization{R}

    A HEALPix pixelization of the sphere with resolution factor `nside`.

    # Examples
    ```jldoctest; setup = :(import CMB.SphericalHarmonics: HealpixPixelization)
    julia> HealpixPixelization(256)
    Nside=256 CMB.SphericalHarmonics.HealpixPixelization{Float64}
    ```
    """
    struct HealpixPixelization{R} <: AbstractRingPixelization{R}
        nside::Int
    end
    HealpixPixelization(nside) = HealpixPixelization{Float64}(Int(nside))

    function Base.show(io::IO, ::MIME"text/plain", hpix::HealpixPixelization)
        if get(io, :compact, false)
            show(io, hpix)
        else
            print(io, "Nside=", hpix.nside, " ", summary(hpix))
        end
    end

    function rings(hpix::HealpixPixelization)
        R = SHT._pixeltype(hpix)
        nside = hpix.nside
        nr, npix = 2nside, nside2npix(nside)
        ΔΩ = 4R(π) / npix

        function healpixring(ii)
            if ii < nside
                p = Healpix.nside2npixcap(ii)
                nϕ = 4ii
            else
                p = Healpix.nside2npixcap(nside) + 4 * nside * (ii - nside)
                nϕ = 4nside
            end
            cosθ = R(Healpix.pix2z(nside, p))
            ϕ₀_π = R(Healpix._pix2phibypi(nside, p))
            o₁ = p + 1
            o₂ = ii < 2nside ? npix - nϕ + 1 - p : 0
            return Ring{R}((o₁, o₂), 1, nϕ, cosθ, ϕ₀_π, ΔΩ)
        end

        return Broadcast.instantiate(Broadcast.broadcasted(healpixring, 1:nr))
    end

    nring(hpix::HealpixPixelization) = nside2nring(hpix.nside)
    npix(hpix::HealpixPixelization) = nside2npix(hpix.nside)
    nϕmax(hpix::HealpixPixelization) = 4hpix.nside
end
