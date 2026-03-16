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

function Base.show(io::IO, mapping::AsinhMapping)
    print(
        io,
        "AsinhMapping(a=",
        mapping.a,
        ", s=",
        mapping.s,
        ", tail_spacing=",
        mapping.tail_spacing,
        ")",
    )
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

"""
    fit_asinh_mapping_for_extent(; npoints, xmax, reference_spacing=1.0, tail_spacing=10.0)

Build the symmetric one-parameter `AsinhMapping(c=s, s=s)` used in the first
mapped ordinary Cartesian hydrogen path, choosing `s` so that the outer basis
centers of an odd-count full-line basis land at `x = ±xmax`.

This helper is intentionally narrow:

- odd `npoints`
- full-line basis centered at the origin
- reference centers spaced by `reference_spacing`
- the single-parameter family `AsinhMapping(c=s, s=s, tail_spacing=...)`

The endpoint condition is

```text
u(xmax) = ((npoints - 1) / 2) * reference_spacing
```

which gives the closed-form choice

```text
s = asinh(xmax) / ( ((npoints - 1) / 2) * reference_spacing - xmax / tail_spacing ).
```
"""
function fit_asinh_mapping_for_extent(;
    npoints::Int,
    xmax::Real,
    reference_spacing::Real = 1.0,
    tail_spacing::Real = 10.0,
)
    isodd(npoints) || throw(ArgumentError("fit_asinh_mapping_for_extent requires an odd npoints"))
    npoints >= 1 || throw(ArgumentError("fit_asinh_mapping_for_extent requires npoints >= 1"))

    xmax_value = Float64(xmax)
    spacing_value = Float64(reference_spacing)
    tail_value = Float64(tail_spacing)
    xmax_value > 0.0 || throw(ArgumentError("fit_asinh_mapping_for_extent requires xmax > 0"))
    spacing_value > 0.0 || throw(ArgumentError("fit_asinh_mapping_for_extent requires reference_spacing > 0"))
    tail_value > 0.0 || throw(ArgumentError("fit_asinh_mapping_for_extent requires tail_spacing > 0"))

    uedge = 0.5 * (npoints - 1) * spacing_value
    denominator = uedge - xmax_value / tail_value
    denominator > 0.0 ||
        throw(ArgumentError("fit_asinh_mapping_for_extent requires the endpoint u-range to exceed xmax / tail_spacing"))

    s_value = asinh(xmax_value) / denominator
    return AsinhMapping(c = s_value, s = s_value, tail_spacing = tail_value)
end
