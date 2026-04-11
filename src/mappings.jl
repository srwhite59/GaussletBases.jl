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

"""
    CombinedInvsqrtMapping(; centers, core_ranges, amplitudes, tail_spacing=10.0)

Smooth multi-center inverse-sqrt-density mapping

    dudx(x) = 1 / tail_spacing + sum(amplitudes[i] / sqrt((x - centers[i])^2 + core_ranges[i]^2))

This is the narrow modernized analogue of the legacy combined multi-center
inverse-sqrt mapping used on molecular bond axes. The local densities at the
requested centers are set by solving for `amplitudes`.
"""
struct CombinedInvsqrtMapping <: AbstractCoordinateMapping
    centers::Vector{Float64}
    core_ranges::Vector{Float64}
    amplitudes::Vector{Float64}
    tail_spacing::Float64

    function CombinedInvsqrtMapping(
        centers::AbstractVector{<:Real},
        core_ranges::AbstractVector{<:Real},
        amplitudes::AbstractVector{<:Real},
        tail_spacing::Real,
    )
        length(centers) == length(core_ranges) == length(amplitudes) || throw(
            ArgumentError("CombinedInvsqrtMapping requires centers, core_ranges, and amplitudes of equal length"),
        )
        isempty(centers) && throw(ArgumentError("CombinedInvsqrtMapping requires at least one center"))
        center_values = Float64[Float64(value) for value in centers]
        core_values = Float64[Float64(value) for value in core_ranges]
        amplitude_values = Float64[Float64(value) for value in amplitudes]
        tail_value = Float64(tail_spacing)
        all(value -> value > 0.0, core_values) || throw(ArgumentError("CombinedInvsqrtMapping requires core_ranges > 0"))
        all(value -> value > 0.0, amplitude_values) || throw(ArgumentError("CombinedInvsqrtMapping requires amplitudes > 0"))
        tail_value > 0.0 || throw(ArgumentError("CombinedInvsqrtMapping requires tail_spacing > 0"))
        new(center_values, core_values, amplitude_values, tail_value)
    end
end

function Base.show(io::IO, mapping::CombinedInvsqrtMapping)
    print(
        io,
        "CombinedInvsqrtMapping(centers=",
        mapping.centers,
        ", core_ranges=",
        mapping.core_ranges,
        ", amplitudes=",
        mapping.amplitudes,
        ", tail_spacing=",
        mapping.tail_spacing,
        ")",
    )
end

function CombinedInvsqrtMapping(;
    centers::AbstractVector{<:Real},
    core_ranges::AbstractVector{<:Real},
    amplitudes::AbstractVector{<:Real},
    tail_spacing::Real = 10.0,
)
    return CombinedInvsqrtMapping(centers, core_ranges, amplitudes, tail_spacing)
end

function dudx(mapping::CombinedInvsqrtMapping, x::Real)
    xval = Float64(x)
    value = 1.0 / mapping.tail_spacing
    for index in eachindex(mapping.centers)
        delta = xval - mapping.centers[index]
        value += mapping.amplitudes[index] / sqrt(delta * delta + mapping.core_ranges[index]^2)
    end
    return value
end

function du2dx2(mapping::CombinedInvsqrtMapping, x::Real)
    xval = Float64(x)
    value = 0.0
    for index in eachindex(mapping.centers)
        delta = xval - mapping.centers[index]
        denom = sqrt(delta * delta + mapping.core_ranges[index]^2)
        value -= mapping.amplitudes[index] * delta / denom^3
    end
    return value
end

function uofx(mapping::CombinedInvsqrtMapping, x::Real)
    xval = Float64(x)
    value = xval / mapping.tail_spacing
    for index in eachindex(mapping.centers)
        value += mapping.amplitudes[index] * asinh((xval - mapping.centers[index]) / mapping.core_ranges[index])
    end
    return value
end

function xofu(mapping::CombinedInvsqrtMapping, u::Real)
    uval = Float64(u)
    center_bound = maximum(abs, mapping.centers)
    core_bound = maximum(mapping.core_ranges)
    hi = max(center_bound + core_bound, mapping.tail_spacing * abs(uval), 1.0)
    lo = -hi

    while uofx(mapping, lo) > uval
        hi = lo
        lo *= 2.0
    end
    while uofx(mapping, hi) < uval
        lo = hi
        hi *= 2.0
    end

    for _ in 1:200
        mid = 0.5 * (lo + hi)
        if uofx(mapping, mid) < uval
            lo = mid
        else
            hi = mid
        end
    end
    return 0.5 * (lo + hi)
end

"""
    fit_combined_invsqrt_mapping(; centers, core_ranges, target_spacings, tail_spacing=10.0)

Fit a positive-amplitude `CombinedInvsqrtMapping` whose local spacing at each
requested center matches `target_spacings`.

The fit solves the same linear matching condition used by the legacy
multi-center inverse-sqrt construction:

```text
1 / d_core[j] = 1 / tail_spacing + sum_i amplitudes[i] / sqrt((x_j - x_i)^2 + a_i^2)
```
"""
function fit_combined_invsqrt_mapping(;
    centers::AbstractVector{<:Real},
    core_ranges::AbstractVector{<:Real},
    target_spacings::AbstractVector{<:Real},
    tail_spacing::Real = 10.0,
)
    length(centers) == length(core_ranges) == length(target_spacings) || throw(
        ArgumentError("fit_combined_invsqrt_mapping requires centers, core_ranges, and target_spacings of equal length"),
    )
    isempty(centers) && throw(ArgumentError("fit_combined_invsqrt_mapping requires at least one center"))

    center_values = Float64[Float64(value) for value in centers]
    core_values = Float64[Float64(value) for value in core_ranges]
    spacing_values = Float64[Float64(value) for value in target_spacings]
    tail_value = Float64(tail_spacing)

    all(value -> value > 0.0, core_values) || throw(ArgumentError("fit_combined_invsqrt_mapping requires core_ranges > 0"))
    all(value -> value > 0.0, spacing_values) || throw(ArgumentError("fit_combined_invsqrt_mapping requires target_spacings > 0"))
    tail_value > 0.0 || throw(ArgumentError("fit_combined_invsqrt_mapping requires tail_spacing > 0"))

    tail_density = 1.0 / tail_value
    target_density = 1.0 ./ spacing_values
    all(value -> value > tail_density, target_density) || throw(
        ArgumentError("fit_combined_invsqrt_mapping requires each target spacing to be smaller than tail_spacing"),
    )

    matrix = zeros(Float64, length(center_values), length(center_values))
    for row in eachindex(center_values), column in eachindex(center_values)
        delta = center_values[row] - center_values[column]
        matrix[row, column] = 1.0 / sqrt(delta * delta + core_values[column]^2)
    end
    amplitudes = matrix \ (target_density .- tail_density)
    all(isfinite, amplitudes) || throw(ArgumentError("fit_combined_invsqrt_mapping produced non-finite amplitudes"))
    all(value -> value > 0.0, amplitudes) || throw(
        ArgumentError("fit_combined_invsqrt_mapping requires a positive-amplitude solution for the requested centers and spacings"),
    )
    return CombinedInvsqrtMapping(center_values, core_values, amplitudes, tail_value)
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

"""
    white_lindsey_atomic_mapping(; Z, d, tail_spacing=10.0)

Return the one-center atomic White-Lindsey-style mapping as a repo-native
`AsinhMapping`.

This is the old one-center atomic inverse-sqrt/asinh family written in the
modern repo parameterization:

```text
u(x) = x / tail_spacing + asinh(x / a) / s
a = sqrt(d / Z)
s = sqrt(d Z)
```

Equivalently, the modern near-origin parameter is

```text
c = a s = d
```

so this helper is exactly the same map as
`AsinhMapping(c=d, s=sqrt(d * Z), tail_spacing=tail_spacing)`.

This helper is only for the one-center atomic White-Lindsey contract. The
legacy multi-center `getmapping(...)` construction is a separate combined
inverse-sqrt-density path and is not represented by this helper.
"""
function white_lindsey_atomic_mapping(; Z::Real, d::Real, tail_spacing::Real = 10.0)
    zval = Float64(Z)
    dval = Float64(d)
    tail = Float64(tail_spacing)
    zval > 0.0 || throw(ArgumentError("white_lindsey_atomic_mapping requires Z > 0"))
    dval > 0.0 || throw(ArgumentError("white_lindsey_atomic_mapping requires d > 0"))
    tail > 0.0 || throw(ArgumentError("white_lindsey_atomic_mapping requires tail_spacing > 0"))
    return AsinhMapping(c = dval, s = sqrt(dval * zval), tail_spacing = tail)
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

Build the symmetric one-parameter `AsinhMapping(c=s, s=s)` used in the current
working mapped ordinary Cartesian hydrogen path, choosing `s` so that the
outer basis centers of an odd-count full-line basis land at `x = ±xmax`.

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

This helper should be read as a narrow current working choice, not as a claim
that `AsinhMapping` or the present coupled `c,s` tuning is already final for
the ordinary mapped branch.
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

"""
    fit_asinh_mapping_for_strength(; s, npoints, xmax, reference_spacing=1.0, tail_spacing=10.0)

Build an `AsinhMapping(a=a_fit, s=s, tail_spacing=tail_spacing)` with a fixed
distortion strength `s`, choosing `a_fit` so that the outer centers of an odd
full-line basis land at `x = ±xmax`.

This is the fixed-extent companion to `fit_asinh_mapping_for_extent`. It is
intentionally narrow:

- odd `npoints`
- full-line basis centered at the origin
- fixed `s`
- `a` chosen by a one-dimensional solve
- the endpoint condition

```text
u(xmax) = ((npoints - 1) / 2) * reference_spacing
```

This is again a narrow current working helper for mapped ordinary studies, not
a claim that the present mapping family or `c,s` heuristics are settled.
"""
function fit_asinh_mapping_for_strength(;
    s::Real,
    npoints::Int,
    xmax::Real,
    reference_spacing::Real = 1.0,
    tail_spacing::Real = 10.0,
)
    isodd(npoints) || throw(ArgumentError("fit_asinh_mapping_for_strength requires an odd npoints"))
    npoints >= 1 || throw(ArgumentError("fit_asinh_mapping_for_strength requires npoints >= 1"))

    s_value = Float64(s)
    xmax_value = Float64(xmax)
    spacing_value = Float64(reference_spacing)
    tail_value = Float64(tail_spacing)
    s_value > 0.0 || throw(ArgumentError("fit_asinh_mapping_for_strength requires s > 0"))
    xmax_value > 0.0 || throw(ArgumentError("fit_asinh_mapping_for_strength requires xmax > 0"))
    spacing_value > 0.0 || throw(ArgumentError("fit_asinh_mapping_for_strength requires reference_spacing > 0"))
    tail_value > 0.0 || throw(ArgumentError("fit_asinh_mapping_for_strength requires tail_spacing > 0"))

    uedge = 0.5 * (npoints - 1) * spacing_value
    uedge > xmax_value / tail_value ||
        throw(ArgumentError("fit_asinh_mapping_for_strength requires the endpoint u-range to exceed xmax / tail_spacing"))

    target(a) = xmax_value / tail_value + asinh(xmax_value / a) / s_value - uedge

    lo = 1.0e-12
    hi = max(1.0, xmax_value)
    while target(hi) > 0.0
        hi *= 2.0
    end

    for _ in 1:200
        mid = 0.5 * (lo + hi)
        if target(mid) > 0.0
            lo = mid
        else
            hi = mid
        end
    end

    return AsinhMapping(a = 0.5 * (lo + hi), s = s_value, tail_spacing = tail_value)
end
