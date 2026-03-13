function _primitive_physical_bounds(primitive::Gaussian)
    return _reference_bounds(primitive)
end

function _primitive_physical_bounds(primitive::HalfLineGaussian)
    return _reference_bounds(primitive)
end

function _primitive_physical_bounds(primitive::XGaussian)
    return _reference_bounds(primitive)
end

function _primitive_physical_bounds(primitive::Distorted)
    ulo, uhi = _reference_bounds(primitive.primitive)
    return xofu(primitive.mapping, ulo), xofu(primitive.mapping, uhi)
end

function _basis_support_bounds(basis)
    lowers = Float64[]
    uppers = Float64[]
    for primitive in primitives(basis)
        xlo, xhi = _primitive_physical_bounds(primitive)
        push!(lowers, xlo)
        push!(uppers, xhi)
    end
    return minimum(lowers), maximum(uppers)
end

function _make_midpoint_grid(xmin::Float64, xmax::Float64, h::Float64)
    xmax > xmin || throw(ArgumentError("grid upper bound must exceed lower bound"))
    h > 0.0 || throw(ArgumentError("grid spacing must be positive"))

    edges = Float64[xmin]
    x = xmin
    while x + h < xmax
        x += h
        push!(edges, x)
    end
    push!(edges, xmax)

    points = Float64[]
    weights = Float64[]
    for i in 1:(length(edges) - 1)
        push!(points, 0.5 * (edges[i] + edges[i + 1]))
        push!(weights, edges[i + 1] - edges[i])
    end
    return points, weights
end

function _make_physical_erf_grid(
    mapping_value::AbstractCoordinateMapping,
    reference_spacing::Float64,
    cutoff::Float64;
    refine::Int,
)
    refine > 0 || throw(ArgumentError("refine must be positive"))
    cutoff > 0.0 || throw(ArgumentError("quadrature cutoff must be positive"))

    h = reference_spacing / refine
    reference_cutoff = _identity_mapping(mapping_value) ? cutoff : uofx(mapping_value, cutoff)
    reference_points, reference_weights = _make_erf_grid(; h = h, rmax = reference_cutoff)

    physical_points = Float64[]
    physical_weights = Float64[]
    for i in eachindex(reference_points)
        x = xofu(mapping_value, reference_points[i])
        push!(physical_points, x)
        push!(physical_weights, reference_weights[i] / dudx(mapping_value, x))
    end
    return physical_points, physical_weights
end

"""
    RadialQuadratureGrid(points, weights; mapping=nothing)

Quadrature grid on the physical half line used for radial and half-line
diagnostics and integrals.
"""
struct RadialQuadratureGrid{M}
    point_data::Vector{Float64}
    weight_data::Vector{Float64}
    mapping_value::M
end

function RadialQuadratureGrid(
    points::AbstractVector{<:Real},
    weights::AbstractVector{<:Real};
    mapping = nothing,
)
    length(points) == length(weights) ||
        throw(ArgumentError("quadrature point and weight counts must match"))
    mapping === nothing || mapping isa AbstractCoordinateMapping ||
        throw(ArgumentError("quadrature mapping must be nothing or an AbstractCoordinateMapping"))
    point_data = Float64[Float64(point) for point in points]
    weight_data = Float64[Float64(weight) for weight in weights]
    return RadialQuadratureGrid{typeof(mapping)}(point_data, weight_data, mapping)
end

"""
    quadrature_points(grid)

Return the physical-space quadrature points of `grid`.
"""
quadrature_points(grid::RadialQuadratureGrid) = grid.point_data

"""
    quadrature_weights(grid)

Return the physical-space quadrature weights of `grid`.
"""
quadrature_weights(grid::RadialQuadratureGrid) = grid.weight_data

"""
    radial_quadrature(basis::RadialBasis; refine=8, rmax)

Build a fine physical-space quadrature grid matched to `basis`.

`rmax` is an explicit physical-space cutoff for the quadrature. Construction
grids used while building the basis are separate internal objects and are not
set by this API.
"""
function radial_quadrature(basis::RadialBasis; refine::Int = 8, rmax = nothing)
    rmax === nothing && throw(ArgumentError("radial_quadrature requires explicit rmax"))
    cutoff = Float64(rmax)
    cutoff > 0.0 || throw(ArgumentError("radial_quadrature requires rmax > 0"))
    points, weights = _make_physical_erf_grid(mapping(basis), basis.spec.reference_spacing, cutoff; refine = refine)
    return RadialQuadratureGrid(points, weights; mapping = mapping(basis))
end

function _moment_center_from_points_weights(
    f::AbstractFunction1D,
    points::AbstractVector{Float64},
    weights::AbstractVector{Float64},
)
    values = Float64[value(f, point) for point in points]
    denominator = sum(weights .* values)
    denominator == 0.0 && return NaN
    numerator = sum(weights .* points .* values)
    return numerator / denominator
end

"""
    moment_center(f, grid)

Return the first-moment center of `f` evaluated on `grid`.
"""
function moment_center(f::AbstractFunction1D, grid)
    return _moment_center_from_points_weights(f, quadrature_points(grid), quadrature_weights(grid))
end

function _primitive_sample_matrix(
    primitive_data::Vector{AbstractPrimitiveFunction1D},
    points::AbstractVector{Float64},
)
    samples = zeros(Float64, length(points), length(primitive_data))
    for mu in eachindex(primitive_data)
        samples[:, mu] = [value(primitive_data[mu], point) for point in points]
    end
    return samples
end

function _basis_values_matrix(basis, points::AbstractVector{Float64})
    primitive_values = _primitive_sample_matrix(primitives(basis), points)
    return primitive_values * stencil_matrix(basis)
end

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

function _radial_diagnostics_grid(basis::RadialBasis; refine::Int = 8)
    cutoff =
        if basis.spec.rmax === nothing
            _, xhi = _basis_support_bounds(basis)
            xhi
        else
            basis.spec.rmax
        end
    points, weights = _make_physical_erf_grid(mapping(basis), basis.spec.reference_spacing, cutoff; refine = refine)
    return RadialQuadratureGrid(points, weights; mapping = mapping(basis))
end

function basis_diagnostics(basis::HalfLineBasis)
    return basis_diagnostics(basis, _halfline_diagnostics_grid(basis))
end

function basis_diagnostics(basis::RadialBasis)
    return basis_diagnostics(basis, _radial_diagnostics_grid(basis))
end
