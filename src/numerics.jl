using SparseArrays

# COV_EXCL_START

@static if VERSION >= v"1.6.0-DEV.292"
    using Base: sincospi
else
    sincospi(x) = (sinpi(x), cospi(x))
end

# TODO: When JuliaLang/julia#35792 is merged, replace with explicit version bounds
@static if isdefined(Base, :cispi)
    using Base: cispi
else
    cispi(z::Real) = complex(reverse(sincospi(z))...)
    cispi(z::Complex) = ((s, c) = sincospi(z); complex(real(c) - imag(s), imag(c) + real(s)))
end

# COV_EXCL_STOP

"""
    quadprod(A, b, n, dir=:col)

Computes the quadratic product ``ABA^\\top`` efficiently for the case where ``B`` is all zero
except for the `n`th column or row vector `b`, for `dir = :col` or `dir = :row`,
respectively.
"""
@inline function quadprod(A, b, n, dir::Symbol=:col)
    if dir == :col
        return (A * sparse(b)) * view(A, :, n)'
    elseif dir == :row
        return view(A, :, n) * (A * sparse(b))'
    else
        error("Unrecognized direction `dir = $(repr(dir))`.")
    end
end
