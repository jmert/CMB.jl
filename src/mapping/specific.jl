import ..Healpix

HealpixPixel = SimplePixel{Int32,:Healpix}
struct HealpixMapDefn{Cs,Pc} <: AbstractPixelizedMap{HealpixPixel{Cs},Pc}
  nside::Int32
end

nrings(m::HealpixMapDefn) = Healpix.nside2nring(m.nside)
npixels(m::HealpixMapDefn) = Healpix.nside2npix(m.nside)
function ringlength(m::HealpixMapDefn, r)
    nr = nrings(m)
    (r < 1 || r > nr) && throw(BoundsError())
    if r < m.nside
        return 4*r
    elseif r ≤ 3*m.nside
        return 4*m.nside
    else
        return 4*(nr - r + 1)
    end
end

function pixelphi(m::HealpixMapDefn{Tp}, p::Int32) where Tp<:ColatAzCoordinates
    return Healpix.pix2phi_ring(m.nside, p)
end

function pixeltheta(m::HealpixMapDefn{Tp}, p::Int32) where Tp<:ColatAzCoordinates
    return Healpix.pix2theta_ring(m.nside, p)
end

let RI, PI, RPI
    RI{Cs,Pc} = RingIterator{HealpixMapDefn{Cs,Pc}}
    PI{Cs,Pc} = PixelIterator{HealpixMapDefn{Cs,Pc}}
    RPI{Cs,Pc} = RingPixelIterator{HealpixMapDefn{Cs,Pc}}

    Base.next(I::PI{Cs},  state) where {Cs} = (state%Int32, state+1)
    function Base.next(I::RPI{Cs}, state) where {Cs}
        m = I.m
        r = I.r
        nr = nrings(m)
        pixoff = 0
        for ii in 1:(r-1)
            pixoff += ringlength(m, ii)
        end
        pix = (pixoff + state) % Int32
        return (pix, state+1)
    end
end

####

ECPPixel = SimplePixel{Tuple{Int32,Int32},:ECP,LatLonCoordinates}
struct ECPMapPatchDefn <: AbstractPixelizedMap{ECPPixel,IAUPolarization}
    # Required to uniquely identify the map
    lx::Float64
    hx::Float64
    ly::Float64
    hy::Float64
    ps::Float64
    racen::Float64

    # Filled in to be convenient
    xdos::Float64
    ydos::Float64
    nx::Int32
    ny::Int32
    x_tic::typeof(0.0:0.1:1.0)
    y_tic::typeof(0.0:0.1:1.0)

    function ECPMapPatchDefn(lx,hx,ly,hy,ps,racen)
        xdos = (hx - lx) * cosd(racen)
        ydos = hy - ly
        nx = round(Int32, xdos / ps)
        ny = round(Int32, ydos / ps)
        Δx = (hx - lx) / nx
        Δy = (hy - ly) / ny
        x_tic = (lx+Δx/2):Δx:(hx-Δx/2)
        y_tic = (ly+Δy/2):Δy:(hy-Δy/2)

        return new(lx, hx, ly, hy, ps, racen, xdos, ydos, nx, ny, x_tic, y_tic)
    end
end

nrings(m::ECPMapPatchDefn) = m.ny
npixels(m::ECPMapPatchDefn) = m.nx * m.ny
ringlength(m::ECPMapPatchDefn, r) = 0<r≤m.ny ? m.nx : throw(BoundsError())

function pixellon(m::ECPMapPatchDefn, p::Tuple{Int32,Int32})
    return m.x_tic[p[2]]
end

function pixellat(m::ECPMapPatchDefn, p::Tuple{Int32,Int32})
    return m.y_tic[p[1]]
end

let RI, PI, RPI
    RI = RingIterator{ECPMapPatchDefn}
    PI = PixelIterator{ECPMapPatchDefn}
    RPI = RingPixelIterator{ECPMapPatchDefn}

    Base.next(I::PI,  state) = ( (x,y)=fldmod1(state%Int32, I.m.ny); ((y,x), state+1) )
    Base.next(I::RPI, state) = ( (I.r%Int32,state%Int32), state+1 )
end

##
# Specific definitions
##

for npow in 0:13
    nside = 2^npow
    symb = Symbol("HealpixMap", nside)
    @eval begin
        export $symb
        $symb() = HealpixMapDefn{ColatAzCoordinates,HealpixPolarization}($nside)
    end
end

# Traditional BICEP Map definition
BicepMapDefn() = ECPMapPatchDefn(-55.0, 55.0, -70.0, -45.0, 0.25, -57.5)
# Expanded BICEP map definition, required to be pixel-by-pixel compatible within the overlap
# region. The definition is defined nominally as in the following line, and then we
# numerically update the definition to match by a single iteration.
let lx=-60.0, hx=+60.0, ly=-73.0, hy=-38.0
    # Required Δx to be compatible
    m₀ = BicepMapDefn();
    Δ₀ = step(m₀.x_tic)
    # Approximate bixepext definition
    m₁ = ECPMapPatchDefn(lx, hx, ly, hy, m₀.ps, m₀.racen)
    Δ₁ = step(m₁.x_tic)
    # accumulate map size mismatch caused by pixel size diff
    Δ = 0.5 * (Δ₀ - Δ₁) * m₁.nx
    # Then declare the global shortcut
    global BicepExtMapDefn() = ECPMapPatchDefn(lx-Δ, hx+Δ, ly, hy, m₀.ps, m₀.racen)
end

