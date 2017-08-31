```@meta
DocTestSetup = quote
    using CMB
    using PyPlot

end
```

# Numerical Validation

A design goal from the outset has been to numerically validate all the routines which
are included in this package. This chapter presents the results of those tests; the
output in this chapter is automatically generated at build time using the code as it
exists in this package.

## Testing ulps

One standard test for numerical algorithms is the “ulps” accuracy; ulps stands for
“unit in the last place”, a measure of whether, within rounding to finite precision,
the answer agrees with the purely mathematical answer. We define an ulp as follows:

Let ``f`` be a function...

## Contents
```@contents
Pages = ["numericalvalidation.md"]
```

## Sphere Functions

We'll first start by systmatically sampling the entire sphere (at relatively low
resolution) to identify if there are any particular problem areas in each of the
spherical algorithms.

We'll sample at 1° resolution across the entire sphere to begin with. Arbitrarily,
we'll choose the point ``(θ₀,ϕ₀) = (π/2, π)`` to be the second coordinate in the
sampling.

```@setup sphere_ulps_grid
using CMB
using PyPlot
macro ulps(fncall)
    (fncall isa Expr && fncall.head == :call) ||
            error("Expected a function call, found $(fncall)")
    fn = esc(fncall.args[1])
    args = esc.(fncall.args[2:end])
    bigargs = Any[:(big($a)) for a in args]

    quote
        begin
            r₁ = $fn($(args...))
            r₂ = $fn($(bigargs...))
            u = (r₁-r₂) / eps(abs(r₁))
            convert(Float64, u)
        end
    end
end
```

```@example sphere_ulps_grid
function sphere_ulps_grid()
    (θ₀,ϕ₀) = (0.5π, π)
    δ = 1.0
    ε = deg2rad(δ) / 2
    Θ = linspace(0+ε,  π-ε, (1/δ)*180)
    Φ = linspace(0+ε, 2π-ε, (1/δ)*360)
    @assert step(Θ) == 2ε # hide
    @assert step(Φ) == 2ε # hide
    Δα₁ = zeros(length(Θ), length(Φ))
    Δα₂ = zeros(length(Θ), length(Φ))
    Δσ  = zeros(length(Θ), length(Φ))
    Δz  = zeros(length(Θ), length(Φ))

    for (i,θ) in enumerate(Θ)
        for (j,ϕ) in enumerate(Φ)
            Δα₁[i,j] = @ulps bearing(θ₀, ϕ₀, θ, ϕ)
            Δα₂[i,j] = @ulps bearing(θ, ϕ, θ₀, ϕ₀)
            Δσ[i,j]  = @ulps distance(θ₀, ϕ₀, θ, ϕ)
            Δz[i,j]  = @ulps cosdistance(θ₀, ϕ₀, θ, ϕ)
        end
    end
    return (Δα₁,Δα₂,Δσ,Δz)
end
Δα₁,Δα₂,Δσ,Δz = sphere_ulps_grid()
nothing # hide
```

```@example sphere_ulps_grid
args = [(:vmin,-10), (:vmax,10), (:cmap, "gray")]
figure(figsize=(6.5,8))
subplot(2,1,1); imshow(Δα₁; args...); gca()[:set_title]("bearing(θ₀, ϕ₀, θ, ϕ)")
subplot(2,1,2); imshow(Δα₂; args...); gca()[:set_title]("bearing(θ, ϕ, θ₀, ϕ₀)")
savefig("sphere_ulps_grid_bearing.png"); nothing # hide
```
![](sphere_ulps_grid_bearing.png)

```@example sphere_ulps_grid
figure(figsize=(6.5,8))
subplot(2,1,1); imshow(Δσ; args...); gca()[:set_title]("distance(θ₀, ϕ₀, θ, ϕ)")
subplot(2,1,2); imshow(Δz; args...); gca()[:set_title]("cosdistance(θ₀, ϕ₀, θ, ϕ)")
savefig("sphere_ulps_grid_distance.png"); nothing # hide
```
![](sphere_ulps_grid_distance.png)


