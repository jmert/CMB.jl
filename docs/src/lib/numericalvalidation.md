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

!!! note

    Numerical accuracy is still very much a work-in-progress. The plots and tests being
    presented here will form a foundation for future work on systematically improving
    numerical accuracy.

## Testing ulps

One standard test for numerical algorithms is the “ulps” accuracy; ulps stands for
“unit in the last place”, a measure of whether, within rounding to finite precision,
the answer agrees with the purely mathematical answer. We define an ulp as follows:

Let ``f`` be a function...

## Contents
```@contents
Pages = ["numericalvalidation.md"]
Depth = 3
```


# Sphere Functions

## Sampling the entire sphere

We'll first start by systmatically sampling the entire sphere (at relatively low
resolution) to identify if there are any particular problem areas in each of the
spherical algorithms. We'll record both the absolute and relative error, and we'll
evaluate the functions symmetrically, giving us 4 matrices of error estimates for
each case.

We'll sample at 1° resolution across the entire sphere to begin with. Arbitrarily,
we'll choose the point ``(θ₀,ϕ₀) = (π/2, π)`` to be the second coordinate in the
sampling.

```@example sphere_grid
using CMB
using CMB.Util: @absrelerr

function sphere_grid(fn)
    (θ₀,ϕ₀) = (0.5π, π)
    δ = 1.0
    ε = deg2rad(δ) / 2
    Θ = linspace(0+ε,  π-ε, (1/δ)*180)
    Φ = linspace(0+ε, 2π-ε, (1/δ)*360)
    @assert step(Θ) == 2ε # hide
    @assert step(Φ) == 2ε # hide
    Δa₁ = zeros(length(Θ), length(Φ))
    Δa₂ = zeros(length(Θ), length(Φ))
    Δr₁ = zeros(length(Θ), length(Φ))
    Δr₂ = zeros(length(Θ), length(Φ))

    for (i,θ) in enumerate(Θ)
        for (j,ϕ) in enumerate(Φ)
            (Δa₁[i,j], Δr₁[i,j]) = @absrelerr fn(θ₀, ϕ₀, θ, ϕ)
            (Δa₂[i,j], Δr₂[i,j]) = @absrelerr fn(θ, ϕ, θ₀, ϕ₀)
        end
    end
    return (Δa₁, Δa₂, Δr₁, Δr₂)
end
nothing # hide
```

We'll be using `matplotlib` via the `PyPlot` and `PyCall` interfaces to visualize the
output, and the following function is how the plots shown below have been generated:

```@example sphere_grid
using PyPlot
using PyCall
@pyimport mpl_toolkits.axes_grid1 as ag1

function plot_sphere_grid(fn)
    (Δa₁,Δa₂,Δr₁,Δr₂) = sphere_grid(fn)
    fnname = split(string(fn), ".")[end]

    fig = figure(figsize=(10,6));

    clf();
    grid = ag1.ImageGrid(fig, 211,
            nrows_ncols=(1,2),
            axes_pad=0.1,
            add_all=true,
            share_all=true,
            cbar_mode="single",
            aspect=true)
    args = [(:vmin,-1e-15), (:vmax,1e-15), (:cmap, "gray"), (:interpolation,"none")]
    grid[1][:imshow](Δa₁; args...)
    grid[2][:imshow](Δa₂; args...)
    grid[1][:set_title]("abs Δ $fnname(θ₀, ϕ₀, θ, ϕ)")
    grid[2][:set_title]("abs Δ $fnname(θ, ϕ, θ₀, ϕ₀)")
    colorbar(grid[2][:images][1], cax=grid[:cbar_axes][1])

    grid = ag1.ImageGrid(fig, 212,
            nrows_ncols=(1,2),
            axes_pad=0.1,
            add_all=true,
            share_all=true,
            cbar_mode="single",
            aspect=true)
    args = [(:vmin,-10), (:vmax,10), (:cmap, "gray"), (:interpolation,"none")]
    grid[1][:imshow](Δr₁; args...)
    grid[2][:imshow](Δr₂; args...)
    grid[1][:set_title]("rel Δ $fnname(θ₀, ϕ₀, θ, ϕ)")
    grid[2][:set_title]("rel Δ $fnname(θ, ϕ, θ₀, ϕ₀)")
    colorbar(grid[2][:images][1], cax=grid[:cbar_axes][1])

    fig[:tight_layout]()

    return nothing
end
nothing # hide
```

```@example sphere_grid
plot_sphere_grid(bearing)
savefig("sphere_grid_bearing.png"); nothing # hide
```
![](sphere_grid_bearing.png)

```@example sphere_grid
plot_sphere_grid(distance)
savefig("sphere_grid_distance.png"); nothing # hide
```
![](sphere_grid_distance.png)

```@example sphere_grid
plot_sphere_grid(cosdistance)
savefig("sphere_grid_cosdistance.png"); nothing # hide
```
![](sphere_grid_cosdistance.png)


