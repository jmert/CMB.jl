importall .Mapping

export BicepMap

struct BicepMap <: AbstractPixelizedMap
    xr::Tuple{Float64,Float64}
    yr::Tuple{Float64,Float64}
    sq::Float64
    ps::Float64
end

function BicepMap()
    lx = -55.0; hx = +55.0;
    ly = -70.0; hy = -45.0;
    sq = -57.5; ps = 0.25;

    return BicepMap((lx,hx), (ly,hy), sq, ps)
end

_bicep_xdos(m::BicepMap) = (m.xr[2]-m.xr[1]) * cosd(m.sq)
_bicep_ydos(m::BicepMap) = (m.yr[2]-m.yr[1])
_bicep_nx(m::BicepMap) = round(Int, _bicep_xdos(m) / m.ps)
_bicep_ny(m::BicepMap) = round(Int, _bicep_ydos(m) / m.ps)

_bicep_x_tic(m::BicepMap) = (sx = (m.xr[2]-m.xr[1]) / _bicep_nx(m);
                             (m.xr[1]+sx/2):sx:(m.xr[2]-sx/2))
_bicep_y_tic(m::BicepMap) = (sy = (m.yr[2]-m.yr[1]) / _bicep_ny(m);
                             (m.yr[1]+sy/2):sy:(m.yr[2]-sy/2))

ringtype(m::Type{BicepMap}) = Int
pixeltype(m::Type{BicepMap}) = Tuple{Int,Int}

numrings(m::BicepMap) = _bicep_ny(m)
numpixels(m::BicepMap) = _bicep_ny(m) * _bicep_nx(m)
ringlength(m::BicepMap, r) = _bicep_nx(m)

let RI, PI, RPI
    RI = RingIterator{BicepMap}
    PI = PixelIterator{BicepMap}
    RPI = RingPixelIterator{BicepMap,Int}

    Base.start(I::RI)  = 1
    Base.start(I::PI)  = 1
    Base.start(I::RPI) = 1

    Base.next(I::RI, state)  = (state, state+1)
    Base.next(I::PI, state)  = ( (x,y)=fldmod1(state, _bicep_ny(I.m)); ((y,x), state+1))
    Base.next(I::RPI, state) = ( (I.r,state), state+1 )

    Base.done(I::RI,  state) = state > numrings(I.m)
    Base.done(I::PI,  state) = state > numpixels(I.m)
    Base.done(I::RPI, state) = state > ringlength(I.m, I.r)
end

