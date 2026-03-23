const _QUADRATURE_GRID_SIGMA = 6.5
const _QUADRATURE_GRID_S0 = 3.0
const _QUADRATURE_MEDIUM_REFINES = (8, 16, 24, 32, 48)
const _QUADRATURE_HIGH_REFINES = (24, 32, 48, 64)
const _QUADRATURE_VERYHIGH_REFINES = (48, 64, 96, 128)
const _QUADRATURE_PUBLIC_OVERLAP_TOL = 1.0e-5

function _quadrature_accuracy_profile(accuracy)
    accuracy isa Symbol ||
        throw(ArgumentError("radial_quadrature requires accuracy to be one of :medium, :high, or :veryhigh"))

    if accuracy == :medium
        return (
            name = :medium,
            refines = _QUADRATURE_MEDIUM_REFINES,
            overlap_tol = 1.0e-6,
            use_stability = false,
            overlap_change_tol = Inf,
            inverse_radius_change_tol = Inf,
            moment_center_change_tol = Inf,
        )
    elseif accuracy == :high
        return (
            name = :high,
            refines = _QUADRATURE_HIGH_REFINES,
            overlap_tol = 1.0e-6,
            use_stability = true,
            overlap_change_tol = 1.0e-10,
            inverse_radius_change_tol = 2.0e-10,
            moment_center_change_tol = 1.0e-12,
        )
    elseif accuracy == :veryhigh
        return (
            name = :veryhigh,
            refines = _QUADRATURE_VERYHIGH_REFINES,
            overlap_tol = 1.0e-6,
            use_stability = true,
            overlap_change_tol = 5.0e-11,
            inverse_radius_change_tol = 5.0e-11,
            moment_center_change_tol = 5.0e-13,
        )
    end

    throw(ArgumentError("radial_quadrature requires accuracy to be one of :medium, :high, or :veryhigh"))
end

function _quadrature_refine_schedule(profile; refine)
    if refine === nothing
        return Int[profile.refines...]
    end

    refine isa Integer ||
        throw(ArgumentError("radial_quadrature requires refine to be an integer or nothing"))

    refine_value = Int(refine)
    refine_value > 0 || throw(ArgumentError("radial_quadrature requires refine > 0"))

    refines = Int[refine_value]
    while length(refines) < length(profile.refines)
        push!(refines, 2 * refines[end])
    end
    return refines
end

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

function _radial_quadrature_tail_bound(basis::RadialBasis)
    quadrature_umax = _radial_quadrature_umax(basis)
    return xofu(mapping(basis), quadrature_umax)
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
    sigma::Real = _QUADRATURE_GRID_SIGMA,
    s0::Real = _QUADRATURE_GRID_S0,
)
    refine > 0 || throw(ArgumentError("refine must be positive"))
    cutoff > 0.0 || throw(ArgumentError("quadrature cutoff must be positive"))

    h = reference_spacing / refine
    reference_cutoff = _identity_mapping(mapping_value) ? cutoff : uofx(mapping_value, cutoff)
    reference_points, reference_weights = _make_erf_grid(; h = h, rmax = reference_cutoff, sigma = sigma, s0 = s0)

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

function Base.show(io::IO, grid::RadialQuadratureGrid)
    print(io, "RadialQuadratureGrid(length=", length(grid.point_data))
    if !isempty(grid.point_data)
        print(io, ", rmin=", first(grid.point_data), ", rmax=", last(grid.point_data))
    end
    print(io, ", mapping=")
    show(io, grid.mapping_value)
    print(io, ")")
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
    radial_quadrature(basis::RadialBasis; accuracy=:high, refine=nothing, quadrature_rmax=nothing)

Build a fine physical-space quadrature grid matched to `basis`.

With no keywords, the routine chooses a conservative quadrature cutoff from the
retained radial basis support itself and uses the `:high` accuracy profile.
Construction/setup extents used while building the basis are separate internal
objects and are not set by this API.

`accuracy` may be `:medium`, `:high`, or `:veryhigh`. `:medium` reproduces the
older cheaper overlap-focused behavior. `:high` is the default and also checks
that simple basis diagnostics have stabilized across refinement. `:veryhigh`
pushes the same checks farther. `refine` is an optional expert starting
resolution hint. `quadrature_rmax` is an optional explicit physical-space
cutoff override kept for expert compatibility. If an explicit cutoff is too
short to cover the retained basis support, the routine warns rather than
silently reporting good overlap on a truncated grid. On the automatic default
path, the routine returns quietly once it reaches the repo's public-quality
overlap regime even if the stricter internal refinement-stability target has
not been met yet.
"""
function radial_quadrature(
    basis::RadialBasis;
    accuracy = :high,
    refine = nothing,
    quadrature_rmax = nothing,
    rmax = nothing,
)
    quadrature_rmax === nothing || rmax === nothing ||
        throw(ArgumentError("provide at most one of quadrature_rmax or rmax"))
    profile = _quadrature_accuracy_profile(accuracy)
    refines = _quadrature_refine_schedule(profile; refine = refine)

    tail_bound = _radial_quadrature_tail_bound(basis)
    cutoff_input = quadrature_rmax === nothing ? rmax : quadrature_rmax
    cutoff = cutoff_input === nothing ? tail_bound : Float64(cutoff_input)
    cutoff > 0.0 || throw(ArgumentError("radial_quadrature requires quadrature_rmax > 0"))
    explicit_cutoff = cutoff_input !== nothing

    previous_metrics = nothing
    latest_grid = nothing
    latest_metrics = nothing
    latest_refine = first(refines)
    for refine_try in refines
        grid = _radial_quadrature_grid(basis, cutoff; refine = refine_try)
        metrics = _quadrature_quality_metrics(
            basis,
            quadrature_points(grid),
            quadrature_weights(grid),
        )
        latest_grid = grid
        latest_metrics = metrics
        latest_refine = refine_try

        if _quadrature_profile_satisfied(metrics, previous_metrics, profile)
            return grid
        end
        previous_metrics = metrics
    end

    if !explicit_cutoff &&
       profile.name != :veryhigh &&
       latest_metrics.overlap_error <= _QUADRATURE_PUBLIC_OVERLAP_TOL
        return latest_grid
    end

    if explicit_cutoff && cutoff < tail_bound
        @warn(
            "radial_quadrature did not reach the requested accuracy; the requested quadrature_rmax may be truncating basis tails",
            accuracy = profile.name,
            refine_start = first(refines),
            best_refine = latest_refine,
            best_overlap_error = latest_metrics.overlap_error,
            requested_quadrature_rmax = cutoff,
            conservative_tail_bound = tail_bound,
        )
    else
        @warn(
            "radial_quadrature did not reach the requested accuracy",
            accuracy = profile.name,
            refine_start = first(refines),
            best_refine = latest_refine,
            best_overlap_error = latest_metrics.overlap_error,
            used_quadrature_rmax = cutoff,
        )
    end
    return latest_grid
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

function _basis_overlap_error_from_points_weights(
    basis,
    points::AbstractVector{Float64},
    weights::AbstractVector{Float64},
)
    values = _basis_values_matrix(basis, points)
    overlap = transpose(values) * (weights .* values)
    return norm(overlap - I, Inf)
end

function _quadrature_quality_metrics(
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

    safe_points = max.(points, eps(Float64))
    inverse_radius_matrix = transpose(values) * ((weights ./ safe_points) .* values)

    return (
        overlap_error = overlap_error,
        overlap = overlap,
        moment_centers = moment_centers,
        inverse_radius_matrix = inverse_radius_matrix,
    )
end

function _quadrature_profile_satisfied(metrics, previous_metrics, profile)
    metrics.overlap_error <= profile.overlap_tol || return false
    profile.use_stability || return true
    previous_metrics === nothing && return false

    overlap_change = norm(metrics.overlap - previous_metrics.overlap, Inf)
    inverse_radius_change = norm(
        metrics.inverse_radius_matrix - previous_metrics.inverse_radius_matrix,
        Inf,
    )
    moment_center_change = maximum(abs.(metrics.moment_centers .- previous_metrics.moment_centers))

    return overlap_change <= profile.overlap_change_tol &&
           inverse_radius_change <= profile.inverse_radius_change_tol &&
           moment_center_change <= profile.moment_center_change_tol
end

function _radial_quadrature_grid(
    basis::RadialBasis,
    cutoff::Float64;
    refine::Int,
)
    points, weights = _make_physical_erf_grid(mapping(basis), basis.spec.reference_spacing, cutoff; refine = refine)
    return RadialQuadratureGrid(points, weights; mapping = mapping(basis))
end
