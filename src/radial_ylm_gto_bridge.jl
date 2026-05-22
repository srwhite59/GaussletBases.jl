"""
    RadialYlmSolidHarmonicGTOFit

Fit result for the narrow radial/Ylm-to-GTO validation bridge.

The fitted functions use the reduced-radial convention
`u(r) = r R(r)`. For angular channel `l`, each centered Gaussian radial
factor is `r^(l + 1) * exp(-beta * r^2)` in the one-dimensional `dr`
metric, corresponding to a solid-harmonic radial convention for
`(u(r) / r) * Y_lm`.

This object is not a Cartesian-GTO shell representation. Cartesian
solid-harmonic-to-Cartesian conversion is intentionally a later bridge slice.
"""
struct RadialYlmSolidHarmonicGTOFit
    radial_convention::Symbol
    l::Int
    m::Int
    exponents::Vector{Float64}
    coefficients::Matrix{Float64}
    diagnostics::NamedTuple
end

coefficients(fit::RadialYlmSolidHarmonicGTOFit) = fit.coefficients

function _radial_ylm_validate_lm(l::Integer, m::Integer)
    l_value = Int(l)
    m_value = Int(m)
    l_value >= 0 ||
        throw(ArgumentError("radial/Ylm GTO fitting requires l >= 0"))
    abs(m_value) <= l_value ||
        throw(ArgumentError("radial/Ylm GTO fitting requires -l <= m <= l"))
    return l_value, m_value
end

function _radial_ylm_validate_exponents(exponents::AbstractVector{<:Real})
    !isempty(exponents) ||
        throw(ArgumentError("radial/Ylm GTO fitting requires at least one exponent"))
    exponent_values = Float64[Float64(exponent) for exponent in exponents]
    all(isfinite, exponent_values) ||
        throw(ArgumentError("radial/Ylm GTO exponents must be finite"))
    all(>(0.0), exponent_values) ||
        throw(ArgumentError("radial/Ylm GTO exponents must be positive"))
    return exponent_values
end

function _radial_ylm_validate_grid(points, weights)
    length(points) == length(weights) ||
        throw(ArgumentError("radial/Ylm GTO fitting requires matching point and weight counts"))
    !isempty(points) ||
        throw(ArgumentError("radial/Ylm GTO fitting requires a nonempty radial grid"))
    point_values = Float64[Float64(point) for point in points]
    weight_values = Float64[Float64(weight) for weight in weights]
    all(isfinite, point_values) ||
        throw(ArgumentError("radial/Ylm GTO quadrature points must be finite"))
    all(isfinite, weight_values) ||
        throw(ArgumentError("radial/Ylm GTO quadrature weights must be finite"))
    all(>=(0.0), point_values) ||
        throw(ArgumentError("radial/Ylm GTO quadrature points must be nonnegative"))
    all(>=(0.0), weight_values) ||
        throw(ArgumentError("radial/Ylm GTO quadrature weights must be nonnegative"))
    any(>(0.0), weight_values) ||
        throw(ArgumentError("radial/Ylm GTO fitting requires at least one positive quadrature weight"))
    return point_values, weight_values
end

function _radial_ylm_column_matrix(
    values::AbstractVector{<:Real},
    expected_rows::Int,
    label::AbstractString,
)
    length(values) == expected_rows ||
        throw(ArgumentError("$label row count must be $expected_rows"))
    return reshape(Float64[Float64(value) for value in values], expected_rows, 1)
end

function _radial_ylm_column_matrix(
    values::AbstractMatrix{<:Real},
    expected_rows::Int,
    label::AbstractString,
)
    size(values, 1) == expected_rows ||
        throw(ArgumentError("$label row count must be $expected_rows"))
    return Matrix{Float64}(values)
end

function _radial_ylm_solid_harmonic_gto_design_matrix(
    points::AbstractVector{<:Real},
    l::Integer,
    exponents::AbstractVector{<:Real},
)
    l_value, _ = _radial_ylm_validate_lm(l, 0)
    exponent_values = _radial_ylm_validate_exponents(exponents)
    point_values = Float64[Float64(point) for point in points]
    all(isfinite, point_values) ||
        throw(ArgumentError("radial/Ylm GTO evaluation points must be finite"))
    all(>=(0.0), point_values) ||
        throw(ArgumentError("radial/Ylm GTO evaluation points must be nonnegative"))

    radial_power = l_value + 1
    design = Matrix{Float64}(undef, length(point_values), length(exponent_values))
    for (mu, beta) in pairs(exponent_values)
        for i in eachindex(point_values)
            r = point_values[i]
            design[i, mu] = (r ^ radial_power) * exp(-beta * r * r)
        end
    end
    return design
end

function _radial_ylm_weighted_column_norms(
    values::AbstractMatrix{<:Real},
    weights::AbstractVector{<:Real},
)
    weighted_abs2 = abs2.(values) .* reshape(Float64.(weights), :, 1)
    return sqrt.(vec(sum(weighted_abs2; dims = 1)))
end

function _fit_radial_ylm_values_to_solid_harmonic_gto(
    points::AbstractVector{<:Real},
    weights::AbstractVector{<:Real},
    reduced_radial_values::Union{AbstractVector{<:Real},AbstractMatrix{<:Real}},
    exponents::AbstractVector{<:Real};
    l::Integer,
    m::Integer = 0,
    metric_svd_cutoff::Real = 1.0e-12,
    regularization::Real = 0.0,
    radial_basis_size = nothing,
)
    l_value, m_value = _radial_ylm_validate_lm(l, m)
    point_values, weight_values = _radial_ylm_validate_grid(points, weights)
    exponent_values = _radial_ylm_validate_exponents(exponents)
    cutoff_value = Float64(metric_svd_cutoff)
    cutoff_value >= 0.0 ||
        throw(ArgumentError("metric_svd_cutoff must be nonnegative"))
    regularization_value = Float64(regularization)
    regularization_value >= 0.0 ||
        throw(ArgumentError("regularization must be nonnegative"))

    targets = _radial_ylm_column_matrix(
        reduced_radial_values,
        length(point_values),
        "reduced_radial_values",
    )
    all(isfinite, targets) ||
        throw(ArgumentError("reduced_radial_values must be finite"))

    design = _radial_ylm_solid_harmonic_gto_design_matrix(point_values, l_value, exponent_values)
    weighted_design = design .* reshape(weight_values, :, 1)
    metric = transpose(design) * weighted_design
    metric = 0.5 .* (metric .+ transpose(metric))
    rhs = transpose(design) * (targets .* reshape(weight_values, :, 1))

    decomposition = svd(metric)
    singular_values = Vector{Float64}(decomposition.S)
    max_singular_value = maximum(singular_values)
    rank_cutoff = cutoff_value * max_singular_value
    retained = singular_values .> rank_cutoff
    effective_rank = count(retained)
    effective_rank > 0 ||
        throw(ArgumentError("Gaussian radial fit metric is numerically rank zero"))
    inverse_values = zeros(Float64, length(singular_values))
    for i in eachindex(singular_values)
        if retained[i]
            sigma = singular_values[i]
            inverse_values[i] = sigma / (sigma * sigma + regularization_value * regularization_value)
        end
    end
    fit_coefficients =
        decomposition.Vt' * (Diagonal(inverse_values) * (decomposition.U' * rhs))

    fitted_values = design * fit_coefficients
    residual = targets - fitted_values
    residual_norms = _radial_ylm_weighted_column_norms(residual, weight_values)
    target_norms = _radial_ylm_weighted_column_norms(targets, weight_values)
    relative_residual_norms = Float64[
        target_norms[i] > 0.0 ? residual_norms[i] / target_norms[i] : residual_norms[i]
        for i in eachindex(residual_norms)
    ]
    retained_singular_values = singular_values[retained]
    metric_condition = maximum(retained_singular_values) / minimum(retained_singular_values)

    diagnostics = (
        radial_convention = :reduced_u_over_r,
        cartesian_conversion_status = :deferred_solid_harmonic_to_cartesian,
        l = l_value,
        m = m_value,
        radial_basis_size = radial_basis_size,
        radial_grid_size = length(point_values),
        orbital_count = size(targets, 2),
        exponent_count = length(exponent_values),
        exponent_min = minimum(exponent_values),
        exponent_max = maximum(exponent_values),
        metric_singular_values = singular_values,
        metric_condition = metric_condition,
        metric_svd_cutoff = cutoff_value,
        metric_rank_cutoff = rank_cutoff,
        effective_rank = effective_rank,
        regularization = regularization_value,
        residual_norms = residual_norms,
        relative_residual_norms = relative_residual_norms,
        target_norms = target_norms,
    )
    return RadialYlmSolidHarmonicGTOFit(
        :reduced_u_over_r,
        l_value,
        m_value,
        exponent_values,
        Matrix{Float64}(fit_coefficients),
        diagnostics,
    )
end

"""
    fit_radial_ylm_to_solid_harmonic_gto(grid, reduced_radial_values, exponents; l, m=0, ...)

Fit one or more reduced radial orbitals `u(r)` on a `RadialQuadratureGrid` to
centered solid-harmonic Gaussian radial factors
`r^(l + 1) * exp(-beta * r^2)` in the one-dimensional `dr` metric.

`reduced_radial_values` may be a vector or a matrix whose columns are separate
orbitals evaluated at `quadrature_points(grid)`. The fit uses an explicit
Gaussian radial metric with SVD truncation and optional Tikhonov
regularization.

The returned `RadialYlmSolidHarmonicGTOFit` is intentionally not a
Cartesian-GTO shell representation; Cartesian conversion is deferred to the
next bridge slice.
"""
function fit_radial_ylm_to_solid_harmonic_gto(
    grid::RadialQuadratureGrid,
    reduced_radial_values::Union{AbstractVector{<:Real},AbstractMatrix{<:Real}},
    exponents::AbstractVector{<:Real};
    l::Integer,
    m::Integer = 0,
    metric_svd_cutoff::Real = 1.0e-12,
    regularization::Real = 0.0,
)
    return _fit_radial_ylm_values_to_solid_harmonic_gto(
        quadrature_points(grid),
        quadrature_weights(grid),
        reduced_radial_values,
        exponents;
        l = l,
        m = m,
        metric_svd_cutoff = metric_svd_cutoff,
        regularization = regularization,
        radial_basis_size = nothing,
    )
end

"""
    fit_radial_ylm_to_solid_harmonic_gto(basis, radial_coefficients, exponents; l, m=0, radial_grid=nothing, ...)

Evaluate `radial_coefficients` in a `RadialBasis` and fit the resulting
reduced radial columns to the same centered solid-harmonic Gaussian radial
factors as the grid-level method.
"""
function fit_radial_ylm_to_solid_harmonic_gto(
    basis::RadialBasis,
    radial_coefficients::Union{AbstractVector{<:Real},AbstractMatrix{<:Real}},
    exponents::AbstractVector{<:Real};
    l::Integer,
    m::Integer = 0,
    radial_grid::Union{Nothing,RadialQuadratureGrid} = nothing,
    quadrature_accuracy = :high,
    metric_svd_cutoff::Real = 1.0e-12,
    regularization::Real = 0.0,
)
    grid = radial_grid === nothing ?
        radial_quadrature(basis; accuracy = quadrature_accuracy) :
        radial_grid
    basis_values = _basis_values_matrix(basis, quadrature_points(grid))
    radial_coefficient_matrix = _radial_ylm_column_matrix(
        radial_coefficients,
        length(basis),
        "radial_coefficients",
    )
    reduced_values = basis_values * radial_coefficient_matrix
    return _fit_radial_ylm_values_to_solid_harmonic_gto(
        quadrature_points(grid),
        quadrature_weights(grid),
        reduced_values,
        exponents;
        l = l,
        m = m,
        metric_svd_cutoff = metric_svd_cutoff,
        regularization = regularization,
        radial_basis_size = length(basis),
    )
end

"""
    evaluate_radial_ylm_gto_fit(fit, grid_or_points)

Evaluate a `RadialYlmSolidHarmonicGTOFit` as reduced radial values `u(r)` on
the supplied `RadialQuadratureGrid` or vector of radial points.
"""
function evaluate_radial_ylm_gto_fit(
    fit::RadialYlmSolidHarmonicGTOFit,
    points::AbstractVector{<:Real},
)
    design = _radial_ylm_solid_harmonic_gto_design_matrix(points, fit.l, fit.exponents)
    return design * fit.coefficients
end

function evaluate_radial_ylm_gto_fit(
    fit::RadialYlmSolidHarmonicGTOFit,
    grid::RadialQuadratureGrid,
)
    return evaluate_radial_ylm_gto_fit(fit, quadrature_points(grid))
end
