const _SQRT_TWO_PI = sqrt(2.0 * pi)

_erfc(x::Real) = ccall((:erfc, Base.Math.libm), Cdouble, (Cdouble,), Float64(x))

function _check_nonnegative_order(order::Int)
    order >= 0 || throw(ArgumentError("derivative order must be nonnegative"))
    return order
end

function _probabilists_hermite(order::Int, z::Float64)
    order == 0 && return 1.0
    order == 1 && return z

    hm2 = 1.0
    hm1 = z
    for n in 1:(order - 1)
        h = z * hm1 - n * hm2
        hm2 = hm1
        hm1 = h
    end
    return hm1
end

function _gaussian_derivative(center_value::Float64, width::Float64, x::Real, order::Int)
    z = (Float64(x) - center_value) / width
    return ((-1.0 / width)^order) * _probabilists_hermite(order, z) * exp(-0.5 * z * z)
end

function _simpson_integral(f::Function, a::Float64, b::Float64; intervals::Int = 800)
    a == b && return 0.0
    n = iseven(intervals) ? intervals : intervals + 1
    h = (b - a) / n
    total = f(a) + f(b)
    for i in 1:(n - 1)
        x = a + i * h
        total += (isodd(i) ? 4.0 : 2.0) * f(x)
    end
    return total * h / 3.0
end

"""
    Gaussian(; center=0.0, width)

Unnormalized full-line Gaussian primitive

    exp(-0.5 * ((x - center) / width)^2)
"""
struct Gaussian <: AbstractPrimitiveFunction1D
    center_value::Float64
    width::Float64

    function Gaussian(; center::Real = 0.0, width::Real)
        width_value = Float64(width)
        width_value > 0.0 || throw(ArgumentError("Gaussian requires width > 0"))
        new(Float64(center), width_value)
    end
end

value(g::Gaussian, x::Real) = exp(-0.5 * ((Float64(x) - g.center_value) / g.width)^2)
center(g::Gaussian) = g.center_value
reference_center(g::Gaussian) = g.center_value
integral_weight(g::Gaussian) = g.width * _SQRT_TWO_PI
stencil(g::Gaussian) = FunctionStencil([1.0], AbstractPrimitiveFunction1D[g])

function derivative(g::Gaussian, x::Real; order::Int = 1)
    derivative_order = _check_nonnegative_order(order)
    return _gaussian_derivative(g.center_value, g.width, x, derivative_order)
end

"""
    HalfLineGaussian(; center=0.0, width)

Unnormalized half-line Gaussian primitive. Values are zero for `x < 0`.
"""
struct HalfLineGaussian <: AbstractPrimitiveFunction1D
    center_value::Float64
    width::Float64

    function HalfLineGaussian(; center::Real = 0.0, width::Real)
        width_value = Float64(width)
        width_value > 0.0 || throw(ArgumentError("HalfLineGaussian requires width > 0"))
        new(Float64(center), width_value)
    end
end

function value(g::HalfLineGaussian, x::Real)
    xval = Float64(x)
    xval < 0.0 && return 0.0
    return exp(-0.5 * ((xval - g.center_value) / g.width)^2)
end

center(g::HalfLineGaussian) = g.center_value
reference_center(g::HalfLineGaussian) = g.center_value
stencil(g::HalfLineGaussian) = FunctionStencil([1.0], AbstractPrimitiveFunction1D[g])

function integral_weight(g::HalfLineGaussian)
    scale = (0.0 - g.center_value) / (sqrt(2.0) * g.width)
    return g.width * sqrt(pi / 2.0) * _erfc(scale)
end

function derivative(g::HalfLineGaussian, x::Real; order::Int = 1)
    derivative_order = _check_nonnegative_order(order)
    xval = Float64(x)
    xval < 0.0 && return 0.0
    return _gaussian_derivative(g.center_value, g.width, xval, derivative_order)
end

"""
    XGaussian(; alpha)

Near-origin half-line primitive

    x * exp(-0.5 * (x / alpha)^2)

Values are zero for `x < 0`. The public `center` is the peak location, so
`center(XGaussian(alpha=a)) == a`.
"""
struct XGaussian <: AbstractPrimitiveFunction1D
    alpha::Float64

    function XGaussian(; alpha::Real)
        alpha_value = Float64(alpha)
        alpha_value > 0.0 || throw(ArgumentError("XGaussian requires alpha > 0"))
        new(alpha_value)
    end
end

function value(g::XGaussian, x::Real)
    xval = Float64(x)
    xval < 0.0 && return 0.0
    return xval * exp(-0.5 * (xval / g.alpha)^2)
end

center(g::XGaussian) = g.alpha
reference_center(g::XGaussian) = g.alpha
integral_weight(g::XGaussian) = g.alpha^2
stencil(g::XGaussian) = FunctionStencil([1.0], AbstractPrimitiveFunction1D[g])

function derivative(g::XGaussian, x::Real; order::Int = 1)
    derivative_order = _check_nonnegative_order(order)
    xval = Float64(x)
    xval < 0.0 && return 0.0
    gaussian = Gaussian(center = 0.0, width = g.alpha)
    return -g.alpha^2 * derivative(gaussian, xval; order = derivative_order + 1)
end

"""
    Distorted(f, mapping)

Explicit mapped primitive wrapper. The mapped value is

    f(uofx(mapping, x)) * sqrt(dudx(mapping, x))
"""
struct Distorted{F <: AbstractPrimitiveFunction1D, M <: AbstractCoordinateMapping} <: AbstractPrimitiveFunction1D
    primitive::F
    mapping::M
end

function value(f::Distorted, x::Real)
    xval = Float64(x)
    return value(f.primitive, uofx(f.mapping, xval)) * sqrt(dudx(f.mapping, xval))
end

center(f::Distorted) = xofu(f.mapping, center(f.primitive))
reference_center(f::Distorted) = reference_center(f.primitive)
stencil(f::Distorted) = FunctionStencil([1.0], AbstractPrimitiveFunction1D[f])

function derivative(f::Distorted, x::Real; order::Int = 1)
    derivative_order = _check_nonnegative_order(order)
    derivative_order == 0 && return value(f, x)
    derivative_order == 1 || throw(ArgumentError("Distorted currently supports derivative order 0 or 1"))

    xval = Float64(x)
    uval = uofx(f.mapping, xval)
    du = dudx(f.mapping, xval)
    du2 = du2dx2(f.mapping, xval)

    return derivative(f.primitive, uval; order = 1) * du^(3 / 2) +
           0.5 * value(f.primitive, uval) * du2 / sqrt(du)
end

function _reference_bounds(f::Gaussian)
    span = 12.0 * f.width
    return f.center_value - span, f.center_value + span
end

function _reference_bounds(f::HalfLineGaussian)
    span = 12.0 * f.width
    return max(0.0, f.center_value - span), f.center_value + span
end

function _reference_bounds(f::XGaussian)
    span = 12.0 * f.alpha
    return 0.0, span
end

_reference_bounds(f::Distorted) = _reference_bounds(f.primitive)

function integral_weight(f::Distorted)
    ulo, uhi = _reference_bounds(f.primitive)
    xlo = xofu(f.mapping, ulo)
    xhi = xofu(f.mapping, uhi)
    return _simpson_integral(x -> value(f, x), xlo, xhi)
end
