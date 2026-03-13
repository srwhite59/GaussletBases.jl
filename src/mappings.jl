"""
    IdentityMapping()

Identity coordinate mapping with `uofx(map, x) == x`.
"""
struct IdentityMapping <: AbstractCoordinateMapping
end

uofx(::IdentityMapping, x::Real) = Float64(x)
xofu(::IdentityMapping, u::Real) = Float64(u)
dudx(::IdentityMapping, x::Real) = 1.0
du2dx2(::IdentityMapping, x::Real) = 0.0

"""
    AsinhMapping(; c, s, tail_spacing=10.0)
    AsinhMapping(; a, s, tail_spacing=10.0)

Forward mapping

    u(x) = x / tail_spacing + asinh(x / a) / s

The linear tail term is built in for this mapping. The keyword `a` is the
direct parameter in the asinh term. The keyword `c` is the derived near-origin
control with

    c = a * s

so `AsinhMapping(c=c0, s=s0)` is equivalent to
`AsinhMapping(a=c0 / s0, s=s0)`.
"""
struct AsinhMapping <: AbstractCoordinateMapping
    a::Float64
    s::Float64
    tail_spacing::Float64

    function AsinhMapping(a::Real, s::Real, tail_spacing::Real)
        aval = Float64(a)
        sval = Float64(s)
        tail = Float64(tail_spacing)
        aval > 0.0 || throw(ArgumentError("AsinhMapping requires a > 0"))
        sval > 0.0 || throw(ArgumentError("AsinhMapping requires s > 0"))
        tail > 0.0 || throw(ArgumentError("AsinhMapping requires tail_spacing > 0"))
        new(aval, sval, tail)
    end
end

function AsinhMapping(;
    c::Union{Nothing, Real} = nothing,
    a::Union{Nothing, Real} = nothing,
    s::Real,
    tail_spacing::Real = 10.0,
)
    (c === nothing) == (a === nothing) &&
        throw(ArgumentError("provide exactly one of c or a"))
    return AsinhMapping(c === nothing ? a : c / s, s, tail_spacing)
end

function uofx(mapping::AsinhMapping, x::Real)
    xval = Float64(x)
    return xval / mapping.tail_spacing + asinh(xval / mapping.a) / mapping.s
end

function dudx(mapping::AsinhMapping, x::Real)
    xval = Float64(x)
    return 1.0 / mapping.tail_spacing + 1.0 / (mapping.s * sqrt(xval * xval + mapping.a * mapping.a))
end

function du2dx2(mapping::AsinhMapping, x::Real)
    xval = Float64(x)
    denom = sqrt(xval * xval + mapping.a * mapping.a)
    return -xval / (mapping.s * denom^3)
end

function xofu(mapping::AsinhMapping, u::Real)
    uval = Float64(u)
    iszero(uval) && return 0.0

    signu = sign(uval)
    target = abs(uval)
    hi = max(mapping.a * sinh(mapping.s * target), mapping.tail_spacing * target)
    hi = max(hi, mapping.a)
    while uofx(mapping, hi) < target
        hi *= 2.0
    end

    lo = 0.0
    for _ in 1:200
        mid = 0.5 * (lo + hi)
        if uofx(mapping, mid) < target
            lo = mid
        else
            hi = mid
        end
    end

    return signu * 0.5 * (lo + hi)
end
