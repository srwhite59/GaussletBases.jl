# PQS raw-source axis transform fact builders.

"""
    pqs_source_axis_transform_facts_from_pgdg_axes(axis_sources;
        source_intervals, source_mode_dims, enforce_symmetric_odd = false)

Build materialized `CartesianRawProductSources.AxisSourceTransformFact`s from
PGDG intermediates or mapped ordinary gausslet bundle axis sources. This is a
narrow source-axis transform builder for PQS raw-source plans; it calls the
low-level doside/PGDG source transform kernel directly and does not call the
old raw-product-box wrapper.
"""
function pqs_source_axis_transform_facts_from_pgdg_axes(
    axis_sources;
    source_intervals,
    source_mode_dims,
    enforce_symmetric_odd::Bool = false,
)
    sources = _pqs_source_axis_transform_axis_tuple(axis_sources, "axis_sources")
    intervals = _pqs_source_axis_transform_interval_tuple(source_intervals)
    requested_dims = _pqs_source_axis_transform_dim_tuple(source_mode_dims)
    axis_symbols = (:x, :y, :z)
    axes = ntuple(
        axis -> _pqs_source_axis_transform_from_pgdg_axis(
            sources[axis],
            intervals[axis],
            requested_dims[axis],
            axis,
            axis_symbols[axis],
            enforce_symmetric_odd,
        ),
        3,
    )
    facts = ntuple(axis -> axes[axis].fact, 3)
    resolved_dims = ntuple(axis -> facts[axis].source_mode_dim, 3)
    overlap_errors =
        ntuple(axis -> axes[axis].coefficient_overlap_error, 3)
    max_overlap_error = maximum(overlap_errors)
    return (;
        object_kind = :pqs_source_axis_transform_fact_set,
        status = :available_pqs_source_axis_transform_facts,
        transform_source = :repo_owned_pgdg_doside_source_axis_transform,
        source_intervals = intervals,
        source_mode_dims_requested = requested_dims,
        source_mode_dims = resolved_dims,
        source_mode_count = prod(resolved_dims),
        source_mode_dims_adjusted = resolved_dims != requested_dims,
        enforce_symmetric_odd,
        axis_transform_facts = facts,
        axis_coefficient_shapes =
            ntuple(axis -> size(facts[axis].coefficient_matrix), 3),
        axis_overlap_errors = overlap_errors,
        max_axis_overlap_error = max_overlap_error,
        source_product_modes_orthogonal = max_overlap_error <= 1.0e-8,
        shell_realization_materialized = false,
        lowdin_cleanup_used = false,
        ida_data_materialized = false,
        hamiltonian_data_materialized = false,
        driver_route_materialized = false,
        artifacts_materialized = false,
    )
end

function _pqs_source_axis_transform_from_pgdg_axis(
    axis_source,
    source_interval::UnitRange{Int},
    source_mode_dim::Int,
    axis::Int,
    axis_symbol::Symbol,
    enforce_symmetric_odd::Bool,
)
    side = ParentGaussletBases._nested_doside_1d(
        axis_source,
        source_interval,
        source_mode_dim;
        enforce_symmetric_odd,
    )
    coefficient_matrix = Matrix{Float64}(side.local_coefficients)
    coefficient_overlap_error = norm(
        transpose(coefficient_matrix) *
        side.local_overlap *
        coefficient_matrix -
        Matrix{Float64}(I, side.retained_count, side.retained_count),
        Inf,
    )
    fact = CRPS.axis_source_transform_fact(
        axis,
        source_interval,
        side.retained_count,
        coefficient_matrix;
        metadata = (;
            transform_source = :repo_owned_pgdg_doside_source_axis_transform,
            axis_symbol,
            requested_source_mode_dim = source_mode_dim,
            resolved_source_mode_dim = side.retained_count,
            source_mode_dim_adjusted = side.retained_count != source_mode_dim,
            enforce_symmetric_odd,
            coefficient_overlap_error,
            shell_realization_materialized = false,
            lowdin_cleanup_used = false,
            ida_data_materialized = false,
            hamiltonian_data_materialized = false,
            driver_route_materialized = false,
            artifacts_materialized = false,
        ),
    )
    return (;
        fact,
        coefficient_overlap_error,
    )
end

function _pqs_source_axis_transform_axis_tuple(value, name::AbstractString)
    if value isa NamedTuple && all(key -> haskey(value, key), (:x, :y, :z))
        return (value.x, value.y, value.z)
    end
    value isa Tuple && length(value) == 3 ||
        throw(ArgumentError("$(name) must provide x, y, and z axis sources"))
    return value
end

function _pqs_source_axis_transform_interval_tuple(value)
    intervals = _pqs_source_axis_transform_axis_tuple(value, "source_intervals")
    return ntuple(
        axis -> _pqs_source_axis_transform_interval(intervals[axis], axis),
        3,
    )
end

function _pqs_source_axis_transform_interval(value::UnitRange{<:Integer}, axis)
    first(value) <= last(value) ||
        throw(ArgumentError("source interval on axis $(axis) must be nonempty"))
    return Int(first(value)):Int(last(value))
end

function _pqs_source_axis_transform_interval(_value, axis)
    throw(ArgumentError("source interval on axis $(axis) must be a UnitRange"))
end

function _pqs_source_axis_transform_dim_tuple(value)
    dims = _pqs_source_axis_transform_axis_tuple(value, "source_mode_dims")
    return ntuple(axis -> _pqs_source_axis_transform_dim(dims[axis], axis), 3)
end

function _pqs_source_axis_transform_dim(value::Integer, axis)
    value > 0 ||
        throw(ArgumentError("source mode dimension on axis $(axis) must be positive"))
    return Int(value)
end

function _pqs_source_axis_transform_dim(_value, axis)
    throw(ArgumentError("source mode dimension on axis $(axis) must be an integer"))
end
