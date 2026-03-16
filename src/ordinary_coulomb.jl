struct CoulombGaussianExpansion
    coefficients::Vector{Float64}
    exponents::Vector{Float64}
    del::Float64
    s::Float64
    c::Float64
    maxu::Float64

    function CoulombGaussianExpansion(
        coefficients::AbstractVector{<:Real},
        exponents::AbstractVector{<:Real};
        del::Real,
        s::Real,
        c::Real,
        maxu::Real,
    )
        length(coefficients) == length(exponents) ||
            throw(ArgumentError("CoulombGaussianExpansion requires matching coefficient and exponent counts"))
        all(exponent -> exponent > 0.0, exponents) ||
            throw(ArgumentError("CoulombGaussianExpansion requires positive exponents"))
        return new(
            Float64[Float64(value) for value in coefficients],
            Float64[Float64(value) for value in exponents],
            Float64(del),
            Float64(s),
            Float64(c),
            Float64(maxu),
        )
    end
end

Base.length(expansion::CoulombGaussianExpansion) = length(expansion.coefficients)

function Base.show(io::IO, expansion::CoulombGaussianExpansion)
    print(
        io,
        "CoulombGaussianExpansion(length=",
        length(expansion),
        ", del=",
        expansion.del,
        ", s=",
        expansion.s,
        ", c=",
        expansion.c,
        ", maxu=",
        expansion.maxu,
        ")",
    )
end

function (expansion::CoulombGaussianExpansion)(r::Real)
    rval = abs(Float64(r))
    return sum(expansion.coefficients .* exp.(-expansion.exponents .* (rval^2)))
end

"""
    _short_gaucoulomb_coefficients(doacc=true; del=nothing, s=nothing, c=nothing, maxu=nothing)

Near-literal internal port of the legacy `ShortGaucoulomb.jl` construction.

Purpose
-------
A deterministic Gaussian expansion of `1/r`:

    sum_i co[i] * exp(-zeta[i] * r^2)  ~=  1/r

This uses the sinh-mapping construction of White and Lindsey. The
coefficients and exponents are produced analytically from a fixed uniform grid
in the auxiliary variable `u`.

Mapping and weights
-------------------
- Grid in `u`: `u = 0.5 * del, 1.5 * del, ...` up to `maxu`
- Map to `x` by: `x(u) = c * sinh(s * u)`
- Convert to Gaussian exponent: `zeta = x(u)^2`
- Weight for each term:

      co = (2 * del / sqrt(pi)) * dx/du

  written in the original code through

      dudx(x) = 1 / (s * sqrt(x^2 + c^2))

  and then

      co = (2 * del / sqrt(pi)) / dudx(x)

The `doacc=true` preset follows the legacy high-accuracy branch with

- `del = 1.0`
- `s = 0.16`
- `c = 0.01`
- `maxu = 135.0`

The lower-cost `doacc=false` preset follows the same file with

- `del = 0.6`
- `s = 0.5`
- `c = 0.03`
- `maxu = 27.0`

The keyword overrides are kept only as a light adaptation for the current repo
wrapper. The core construction is otherwise intentionally a direct
transcription of `ShortGaucoulomb.jl`.
"""
function _short_gaucoulomb_coefficients(
    doacc::Bool = true;
    del::Union{Nothing, Real} = nothing,
    s::Union{Nothing, Real} = nothing,
    c::Union{Nothing, Real} = nothing,
    maxu::Union{Nothing, Real} = nothing,
)
    del_value = 0.6
    maxu_value = 27.0
    s_value = 0.5
    c_value = 0.03
    if doacc
        del_value = 1.0
        s_value = 0.16
        c_value = 0.01
        maxu_value = 135.0
    end

    del !== nothing && (del_value = Float64(del))
    s !== nothing && (s_value = Float64(s))
    c !== nothing && (c_value = Float64(c))
    maxu !== nothing && (maxu_value = Float64(maxu))

    del_value > 0.0 || throw(ArgumentError("_short_gaucoulomb_coefficients requires del > 0"))
    s_value > 0.0 || throw(ArgumentError("_short_gaucoulomb_coefficients requires s > 0"))
    c_value > 0.0 || throw(ArgumentError("_short_gaucoulomb_coefficients requires c > 0"))
    maxu_value > 0.0 || throw(ArgumentError("_short_gaucoulomb_coefficients requires maxu > 0"))

    dudx(x) = 1.0 / sqrt(x^2 + c_value^2) / s_value
    xofu(u) = c_value * sinh(s_value * u)
    uset = collect((0.5 * del_value):del_value:maxu_value)
    xset = xofu.(uset)
    exponents = xset .^ 2
    coefficients = 2.0 * del_value / sqrt(pi) ./ dudx.(xset)
    return coefficients, exponents, (del = del_value, s = s_value, c = c_value, maxu = maxu_value)
end

"""
    coulomb_gaussian_expansion(; doacc = true, del = nothing, s = nothing, c = nothing, maxu = nothing)

Build the deterministic Gaussian expansion

```text
1 / r ≈ sum_k coefficients[k] * exp(-exponents[k] * r^2)
```

using the internal `ShortGaucoulomb`-derived source-of-truth implementation.

By default this uses the legacy high-accuracy branch (`doacc = true`), which in
the present repo means:

- `del = 1.0`
- `s = 0.16`
- `c = 0.01`
- `maxu = 135.0`

No fitting is performed at runtime. The same parameters always produce the same
coefficient and exponent lists.
"""
function coulomb_gaussian_expansion(;
    doacc::Bool = true,
    del::Union{Nothing, Real} = nothing,
    s::Union{Nothing, Real} = nothing,
    c::Union{Nothing, Real} = nothing,
    maxu::Union{Nothing, Real} = nothing,
)
    coefficients, exponents, parameters = _short_gaucoulomb_coefficients(
        doacc;
        del = del,
        s = s,
        c = c,
        maxu = maxu,
    )
    return CoulombGaussianExpansion(
        coefficients,
        exponents;
        del = parameters.del,
        s = parameters.s,
        c = parameters.c,
        maxu = parameters.maxu,
    )
end

function _gaussian_factor(a::Gaussian, b::Gaussian, exponent::Float64, center_value::Float64)
    exponent >= 0.0 || throw(ArgumentError("gaussian_factor_matrix requires exponent >= 0"))
    exponent == 0.0 && return _gaussian_overlap(a, b)

    alpha_a = inv(a.width^2)
    alpha_b = inv(b.width^2)
    alpha_g = 2.0 * exponent
    total_alpha = alpha_a + alpha_b + alpha_g
    weighted_center =
        (alpha_a * a.center_value + alpha_b * b.center_value + alpha_g * center_value) / total_alpha
    constant_term =
        alpha_a * a.center_value^2 +
        alpha_b * b.center_value^2 +
        alpha_g * center_value^2
    return sqrt(2.0 * pi / total_alpha) *
           exp(-0.5 * (constant_term - total_alpha * weighted_center^2))
end

function _primitive_gaussian_factor_matrix(
    set::PrimitiveSet1D,
    ::_AnalyticPrimitiveMatrixBackend;
    exponent::Float64,
    center::Float64,
)
    matrix = zeros(Float64, length(set), length(set))
    for a in 1:length(set)
        pa = primitives(set)[a]
        for b in a:length(set)
            pb = primitives(set)[b]
            value_ab = _gaussian_factor(pa, pb, exponent, center)
            matrix[a, b] = value_ab
            matrix[b, a] = value_ab
        end
    end
    return matrix
end

function _primitive_gaussian_factor_matrix(
    set::PrimitiveSet1D,
    ::_NumericalPrimitiveMatrixBackend;
    exponent::Float64,
    center::Float64,
    h = nothing,
)
    exponent >= 0.0 || throw(ArgumentError("gaussian_factor_matrix requires exponent >= 0"))
    exponent == 0.0 && return overlap_matrix(set)

    xlo, xhi = _primitive_set_bounds(set)
    h_try = h === nothing ? _primitive_matrix_start_h(set) : Float64(h)
    h_try > 0.0 || throw(ArgumentError("numerical gaussian_factor_matrix requires h > 0"))

    previous = nothing
    current = nothing
    gaussian_factor = point -> exp(-exponent * (point - center)^2)
    for _ in 1:_PRIMITIVE_MATRIX_MAXITER
        points, weights = _make_midpoint_grid(xlo, xhi, h_try)
        values = _primitive_sample_matrix(set, points)
        weighted_factor = weights .* gaussian_factor.(points)
        current = _symmetrize_primitive_matrix(transpose(values) * (weighted_factor .* values))
        if previous !== nothing && norm(current - previous, Inf) <= _PRIMITIVE_MATRIX_TOL
            return current
        end
        previous = current
        h_try /= 2.0
    end
    return current
end

"""
    gaussian_factor_matrix(layer; exponent, center = 0.0)

Build the one-dimensional Gaussianized one-body matrix

```text
<phi_mu | exp(-exponent * (x - center)^2) | phi_nu>
```

for a primitive set or any basis-like layer with a visible primitive set and a
contraction matrix.

For plain full-line Gaussian primitives, the primitive-space matrix uses the
analytic Gaussian triple-product formula. Distorted or otherwise unsupported
primitive content falls back to the same adaptive midpoint sampling strategy
used elsewhere in the primitive-layer matrix code.
"""
function gaussian_factor_matrix(
    set::PrimitiveSet1D;
    exponent::Real,
    center::Real = 0.0,
)
    exponent_value = Float64(exponent)
    center_value = Float64(center)
    return _primitive_gaussian_factor_matrix(
        set,
        _select_primitive_matrix_backend(set);
        exponent = exponent_value,
        center = center_value,
    )
end

function gaussian_factor_matrix(
    basis;
    exponent::Real,
    center::Real = 0.0,
)
    primitive_matrix = gaussian_factor_matrix(
        primitive_set(basis);
        exponent = exponent,
        center = center,
    )
    return contract_primitive_matrix(basis, primitive_matrix)
end
