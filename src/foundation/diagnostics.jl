function _basis_diagnostics_from_points_weights(
    basis,
    points::AbstractVector{Float64},
    weights::AbstractVector{Float64},
)
    values = _basis_values_matrix(basis, points)
    overlap = transpose(values) * (weights .* values)
    overlap_error = norm(overlap - I, Inf)

    denominators = vec(transpose(values) * weights)
    numerators = vec(transpose(values) * (points .* weights))
    moment_centers = similar(numerators)
    for i in eachindex(moment_centers)
        moment_centers[i] = denominators[i] == 0.0 ? NaN : numerators[i] / denominators[i]
    end

    center_mismatches = centers(basis) .- moment_centers
    D = sum(abs2, center_mismatches)
    return (
        overlap_error = overlap_error,
        moment_centers = moment_centers,
        center_mismatches = center_mismatches,
        D = D,
    )
end

"""
    basis_diagnostics(basis)
    basis_diagnostics(basis, grid)

Return a small diagnostics summary for `basis`.

The one-argument form chooses its own integration grid. For radial bases, that
grid is intentionally conservative and uses the same default `accuracy = :high`
quadrature profile as `radial_quadrature(basis)`, so the default diagnostics
are less likely to be dominated by outer-tail truncation.

The returned named tuple currently includes:

- `overlap_error`
- `moment_centers`
- `center_mismatches`
- `D`

For a first check, `diag.overlap_error` and `diag.D` are usually the most
useful fields to inspect.
"""
function basis_diagnostics(basis, grid)
    return _basis_diagnostics_from_points_weights(basis, quadrature_points(grid), quadrature_weights(grid))
end

function basis_diagnostics(basis::UniformBasis)
    xlo, xhi = _basis_support_bounds(basis)
    points, weights = _make_midpoint_grid(xlo, xhi, basis.spec.spacing / 8.0)
    return _basis_diagnostics_from_points_weights(basis, points, weights)
end

function _halfline_diagnostics_grid(basis::HalfLineBasis; refine::Int = 8)
    _, cutoff = _basis_support_bounds(basis)
    points, weights = _make_physical_erf_grid(mapping(basis), basis.spec.reference_spacing, cutoff; refine = refine)
    return RadialQuadratureGrid(points, weights; mapping = mapping(basis))
end

function _radial_diagnostics_grid(basis::RadialBasis; accuracy = :high, refine = nothing)
    return radial_quadrature(basis; accuracy = accuracy, refine = refine)
end

function basis_diagnostics(basis::HalfLineBasis)
    return basis_diagnostics(basis, _halfline_diagnostics_grid(basis))
end

function basis_diagnostics(basis::RadialBasis; accuracy = :high)
    return basis_diagnostics(basis, _radial_diagnostics_grid(basis; accuracy = accuracy))
end
