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

"""
    RadialYlmCartesianGTOAdapter

Centered Cartesian-GTO adapter for a `RadialYlmSolidHarmonicGTOFit`.

The `supplement` field is a repo-compatible
`CartesianGaussianShellSupplementRepresentation3D`. The fitted radial/Ylm
orbitals are obtained as linear combinations of the supplement columns through
`coefficient_map`, so downstream projection uses

    gto_overlap_matrix(working, adapter.supplement) * adapter.coefficient_map

The coefficient map explicitly converts the unnormalized reduced-radial fit
coefficients into the repo's axiswise normalized Cartesian primitive
convention. Diagnostics record the real-Ylm and primitive-normalization
contracts used for that conversion.
"""
struct RadialYlmCartesianGTOAdapter
    supplement::CartesianGaussianShellSupplementRepresentation3D
    coefficient_map::Matrix{Float64}
    source_fit::RadialYlmSolidHarmonicGTOFit
    diagnostics::NamedTuple
end

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

function _radial_ylm_center_tuple(center)
    length(center) == 3 ||
        throw(ArgumentError("radial/Ylm Cartesian adapter center must have length 3"))
    center_tuple = (
        Float64(center[1]),
        Float64(center[2]),
        Float64(center[3]),
    )
    all(isfinite, center_tuple) ||
        throw(ArgumentError("radial/Ylm Cartesian adapter center must be finite"))
    return center_tuple
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

function _radial_ylm_real_solid_harmonic_terms(channel::YlmChannel)
    l = channel.l
    m = channel.m
    l > 2 && throw(
        ArgumentError(
            "radial/Ylm Cartesian GTO adapter currently supports l = 0, 1, 2; got l = $l",
        ),
    )

    if l == 0
        return [(powers = (0, 0, 0), coefficient = inv(sqrt(4.0 * pi)))]
    elseif l == 1
        coefficient = sqrt(3.0 / (4.0 * pi))
        m == -1 && return [(powers = (0, 1, 0), coefficient = -coefficient)]
        m == 0 && return [(powers = (0, 0, 1), coefficient = coefficient)]
        return [(powers = (1, 0, 0), coefficient = -coefficient)]
    end

    coefficient_m2 = sqrt(15.0 / (4.0 * pi))
    coefficient_m1 = -sqrt(15.0 / (4.0 * pi))
    coefficient_m0 = sqrt(5.0 / (16.0 * pi))
    coefficient_p2 = sqrt(15.0 / (16.0 * pi))
    m == -2 && return [(powers = (1, 1, 0), coefficient = coefficient_m2)]
    m == -1 && return [(powers = (0, 1, 1), coefficient = coefficient_m1)]
    if m == 0
        return [
            (powers = (2, 0, 0), coefficient = -coefficient_m0),
            (powers = (0, 2, 0), coefficient = -coefficient_m0),
            (powers = (0, 0, 2), coefficient = 2.0 * coefficient_m0),
        ]
    end
    m == 1 && return [(powers = (1, 0, 1), coefficient = coefficient_m1)]
    return [
        (powers = (2, 0, 0), coefficient = coefficient_p2),
        (powers = (0, 2, 0), coefficient = -coefficient_p2),
    ]
end

function _radial_ylm_polynomial_value(
    terms,
    point::NTuple{3,Float64},
)
    value = 0.0
    x, y, z = point
    for term in terms
        lx, ly, lz = term.powers
        value += term.coefficient * (x ^ lx) * (y ^ ly) * (z ^ lz)
    end
    return value
end

function _radial_ylm_cartesian_primitive_prefactor(
    exponent::Float64,
    powers::NTuple{3,Int},
)
    return prod(
        GaussianAnalyticIntegrals.polynomial_gaussian_shell_prefactor(exponent, power)
        for power in powers
    )
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

"""
    radial_ylm_fit_cartesian_gto_adapter(fit; center=(0,0,0), label_prefix="radial_ylm")

Convert a centered `RadialYlmSolidHarmonicGTOFit` into a
`RadialYlmCartesianGTOAdapter`.

This first adapter slice supports `l = 0, 1, 2`. It expands the repo's real
`YlmChannel(l,m)` convention into Cartesian solid-harmonic polynomials at a
single center, then emits one axiswise normalized Cartesian primitive probe per
`(exponent, polynomial term)`. The returned `coefficient_map` combines those
probe columns back into the fitted radial/Ylm orbital columns.
"""
function radial_ylm_fit_cartesian_gto_adapter(
    fit::RadialYlmSolidHarmonicGTOFit;
    center = (0.0, 0.0, 0.0),
    label_prefix::AbstractString = "radial_ylm",
)
    center_tuple = _radial_ylm_center_tuple(center)
    channel = YlmChannel(fit.l, fit.m)
    terms = _radial_ylm_real_solid_harmonic_terms(channel)
    orbital_count = length(fit.exponents) * length(terms)
    coefficient_map = zeros(Float64, orbital_count, size(fit.coefficients, 2))
    orbitals = Vector{CartesianGaussianShellOrbitalRepresentation3D}(undef, orbital_count)

    row = 0
    for (exponent_index, exponent) in pairs(fit.exponents)
        for (term_index, term) in pairs(terms)
            row += 1
            powers = term.powers
            primitive_prefactor = _radial_ylm_cartesian_primitive_prefactor(exponent, powers)
            primitive_prefactor > 0.0 ||
                throw(ArgumentError("Cartesian primitive prefactor must be positive"))
            coefficient_map[row, :] .=
                (term.coefficient / primitive_prefactor) .* view(fit.coefficients, exponent_index, :)
            label = string(
                label_prefix,
                "_l",
                fit.l,
                "_m",
                fit.m,
                "_e",
                exponent_index,
                "_t",
                term_index,
                "_x",
                powers[1],
                "y",
                powers[2],
                "z",
                powers[3],
            )
            orbitals[row] = CartesianGaussianShellOrbitalRepresentation3D(
                label,
                powers,
                center_tuple,
                [Float64(exponent)],
                [1.0],
                :axiswise_normalized_cartesian_gaussian,
            )
        end
    end

    supplement_metadata = (
        source_kind = :radial_ylm_solid_harmonic_gto_fit,
        radial_convention = fit.radial_convention,
        cartesian_bridge = :centered_real_solid_harmonic_to_cartesian_gto,
        l = fit.l,
        m = fit.m,
        center = center_tuple,
        exponent_count = length(fit.exponents),
        fitted_orbital_count = size(fit.coefficients, 2),
        cartesian_orbital_count = orbital_count,
        supported_lmax = 2,
        real_ylm_convention = :GaussletBases_real_YlmChannel,
        real_ylm_normalization = :included_in_cartesian_coefficient_map,
        solid_harmonic_polynomial_normalization = :repo_real_ylm_r_power,
        cartesian_primitive_normalization = :axiswise_normalized_cartesian_gaussian,
        fit_radial_coefficients = :unnormalized_reduced_radial,
    )
    supplement = CartesianGaussianShellSupplementRepresentation3D(
        :radial_ylm_centered_cartesian_gto_bridge,
        orbitals,
        supplement_metadata,
    )
    diagnostics = (
        radial_convention = fit.radial_convention,
        l = fit.l,
        m = fit.m,
        center = center_tuple,
        exponent_count = length(fit.exponents),
        fitted_orbital_count = size(fit.coefficients, 2),
        cartesian_orbital_count = orbital_count,
        polynomial_terms = terms,
        coefficient_map_size = size(coefficient_map),
        gto_overlap_matrix_compatible = true,
        supported_lmax = 2,
        real_ylm_convention = :GaussletBases_real_YlmChannel,
        real_ylm_normalization = :included_in_cartesian_coefficient_map,
        solid_harmonic_polynomial_normalization = :repo_real_ylm_r_power,
        cartesian_primitive_normalization = :axiswise_normalized_cartesian_gaussian,
        fit_radial_coefficients = :unnormalized_reduced_radial,
    )
    return RadialYlmCartesianGTOAdapter(
        supplement,
        coefficient_map,
        fit,
        diagnostics,
    )
end
