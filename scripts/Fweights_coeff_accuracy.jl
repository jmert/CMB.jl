using PyPlot
using Printf: @sprintf
using Statistics: mean, std

# Stated Fweights coefficients are written in the J. Willmert PhD thesis as
#
#   \begin{align}
#       \eta_\ell &\equiv \frac{2\ell+1}{4\pi}
#           \qquad
#       &
#       \chi_\ell &\equiv \frac{2\eta_\ell \ell}
#           {\sqrt{(\ell-1)\ell(\ell+1)(\ell+2)}}
#       \\
#       {}&{}
#       &
#       \gamma_\ell &\equiv \frac{2\eta_\ell}
#           {(\ell-1)\ell(\ell+1)(\ell+2)}
#   \end{align}
#
# but the convenient mathematical form is not necessarily the best for finite computation.
# In particular, χ which involves a square root:

# N.B.
#   (l-1)l(l+1)(l+2) == -2l - l^2 + 2l^3 + l^4 ==  evalpoly(l, (0, -2, -1, 2, 1))
#                    == l(-2 - l + 2l^2 + l^3) == l * evalpoly(l, (-2, -1, 2, 1))

# Unchecked sqrt that skips the domain check for (x < 0)
usqrt(x::T) where {T <: Base.IEEEFloat} = Base.sqrt_llvm(x)
usqrt(x) = sqrt(x)

if "--nopi" in ARGS
    norm = 1
else
    norm = π
end

# Compare the calculation where the factor of π is included in the normalization and
# when it isn't. The reasoning is that an irrational number must necessarily be rounded,
# and that introduces a bias at the very first step of the recurrence. Delaying the
# normalization until after summing over all terms may improve the accuracy of the total
# result.

let π = norm
    # Direct translation of the equation \chi_\ell above
    global function coeff_χ0(::Type{T}, l::Integer) where T
        η = T(2l + 1) / 4convert(T, π)
        fac = T(2l) / usqrt(evalpoly(T(l), (0, -2, -1, 2, 1)))
        return η * fac
    end

    # Second version which makes the simplification that
    #   l / usqrt(l) == usqrt(l)
    # (also simplify 2/4 ratio between η and fac)
    global function coeff_χ1(::Type{T}, l::Integer) where T
        η = T(2l + 1) / 2convert(T, π)
        fac = usqrt(T(l) / evalpoly(T(l), (-2, -1, 2, 1)))
        return η * fac
    end

    # Break up the η and fac terms - moving /2 inside the usqrt and dividing by π at the end
    global function coeff_χ2(::Type{T}, l::Integer) where T
        lT = convert(T, l)
        fac = usqrt(lT / 4evalpoly(lT, (-2, -1, 2, 1)))
        return (2lT + one(T)) * fac / convert(T, π)
    end

    # Move all terms except the division by π within the usqrt
    global function coeff_χ3(::Type{T}, l::Integer) where T
        lT = convert(T, l)
        fac1 = (2lT + one(lT))^2 * lT
        fac2 = 4 * evalpoly(lT, (-2, -1, 2, 1))
        return usqrt(fac1 / fac2) / convert(T, π)
    end

    # Move even the factor of π into the square root, this time as a multiplication by π²
    global function coeff_χ4(::Type{T}, l::Integer) where T
        lT = convert(T, l)
        fac1 = (2lT + one(lT))^2 * lT
        fac2 = 4 * convert(T, π)^2 * evalpoly(lT, (-2, -1, 2, 1))
        return usqrt(fac1 / fac2)
    end

    # Calculate the sqrt-ratio via 1-round of Newton-Rhapson refinement of the reciprocal-sqrt
    # of the reciprocal-ratio
    global function rsqrt(x)
        y = usqrt(inv(x))
        y = y * muladd(y^2, -x/2, 3/2)
    end
    global function coeff_χ5(::Type{T}, l::Integer) where T
        lT = convert(T, l)
        fac1 = (2lT + one(lT))^2 * lT
        fac2 = 4 * convert(T, π)^2 * evalpoly(lT, (-2, -1, 2, 1))
        return rsqrt(fac2 / fac1)
    end

    # Like case 5, but calcuate the ratio and inverse ratio separately
    global function rsqrt(n, d) # 1 / sqrt(n / d)
        x = n / d
        y = usqrt(d / n)
        y = y * muladd(y^2, -x/2, 3/2)
    end
    global function coeff_χ6(::Type{T}, l::Integer) where T
        lT = convert(T, l)
        fac1 = (2lT + one(lT))^2 * lT
        fac2 = 4 * convert(T, π)^2 * evalpoly(lT, (-2, -1, 2, 1))
        return rsqrt(fac2, fac1)
    end
end

ell = 2:4000
ref = coeff_χ0.(BigFloat, ell)
v0 = coeff_χ0.(Float64, ell)
v1 = coeff_χ1.(Float64, ell)
v2 = coeff_χ2.(Float64, ell)
v3 = coeff_χ3.(Float64, ell)
v3 = coeff_χ3.(Float64, ell)
v4 = coeff_χ4.(Float64, ell)
v5 = coeff_χ5.(Float64, ell)
v6 = coeff_χ6.(Float64, ell)
@assert Float64.(coeff_χ5.(BigFloat, ell)) == Float64.(coeff_χ0.(BigFloat, ell))

bins = collect(range(-2, 2, length=101))
relerr(v) = Float64.((v .- ref) ./ eps.(Float64.(ref)))
function hist_relerr(v; label::String=nothing, kws...)
    r = relerr(v)
    μ = mean(r)
    σ = std(r, mean=μ)
    label = @sprintf("\$\\langle{%s}\\rangle = %0.3f, \\sigma({%s}) = %0.3f\$", label, μ, label, σ)
    hist(relerr(v), bins=bins, histtype="step", linewidth=0.5, label=label; kws...)
end

fig = figure(figsize = (8, 6))
ax = fig.subplots(1, 1)
hist_relerr(v0, label="v0")
hist_relerr(v1, label="v1")
hist_relerr(v2, label="v2")
hist_relerr(v3, label="v3")
hist_relerr(v4, label="v4")
hist_relerr(v5, label="v5")
hist_relerr(v6, label="v6")

title("Iterations of \$χ_\\ell\$ accuracy for algorithm designs")
xlabel("Error vs BigFloat [ulps]")
ylabel("Counts")
legend()

savefig("Fweights_coeff_accuracy_πeq$norm.svg")
