# Metadata-only axis transform fact defaults.
#
# These facts document what a raw product source would need on each axis. They
# do not construct source-mode coefficient matrices. Numerical axis transforms
# remain owned by legacy adapters or future materialization modules.

function _default_axis_transform_fact(
    axis::Int,
    source_interval::UnitRange{Int},
    source_mode_dim::Int,
)
    return AxisSourceTransformFact(
        axis,
        source_interval,
        source_mode_dim,
        :not_materialized,
        nothing,
        (; object_kind = :axis_source_transform_fact, materialized = false),
    )
end

function _default_axis_transform_facts(
    source_intervals::NTuple{3,UnitRange{Int}},
    source_mode_dims::NTuple{3,Int},
)
    return ntuple(
        axis -> _default_axis_transform_fact(
            axis,
            source_intervals[axis],
            source_mode_dims[axis],
        ),
        3,
    )
end

"""
    axis_source_transform_fact(axis, source_interval, source_mode_dim, coefficient_matrix)

Build one externally materialized source-axis transform fact. The coefficient
matrix is validated as `source support rows x source modes`.
"""
function axis_source_transform_fact(
    axis::Int,
    source_interval::UnitRange{Int},
    source_mode_dim::Int,
    coefficient_matrix;
    metadata = (;),
)
    axis in 1:3 ||
        throw(ArgumentError("source-axis transform axis must be 1, 2, or 3"))
    source_mode_dim > 0 ||
        throw(ArgumentError("source-axis transform source_mode_dim must be positive"))
    matrix = _validated_axis_transform_matrix(
        coefficient_matrix,
        source_interval,
        source_mode_dim,
        axis,
    )
    return AxisSourceTransformFact(
        axis,
        source_interval,
        source_mode_dim,
        :materialized,
        matrix,
        merge(
            (;
                object_kind = :axis_source_transform_fact,
                materialized = true,
                coefficient_matrix_shape = size(matrix),
                coefficient_source = :externally_supplied_source_axis_transform,
            ),
            NamedTuple(metadata),
        ),
    )
end

function _axis_transform_facts(
    source_intervals::NTuple{3,UnitRange{Int}},
    source_mode_dims::NTuple{3,Int};
    axis_transform_matrices = nothing,
    axis_transform_facts = nothing,
)
    if !isnothing(axis_transform_matrices) && !isnothing(axis_transform_facts)
        throw(
            ArgumentError(
                "provide either axis_transform_matrices or axis_transform_facts, not both",
            ),
        )
    end
    !isnothing(axis_transform_facts) &&
        return _validated_axis_transform_facts(
            axis_transform_facts,
            source_intervals,
            source_mode_dims,
        )
    !isnothing(axis_transform_matrices) &&
        return _materialized_axis_transform_facts(
            axis_transform_matrices,
            source_intervals,
            source_mode_dims,
        )
    return _default_axis_transform_facts(source_intervals, source_mode_dims)
end

function _materialized_axis_transform_facts(
    axis_transform_matrices,
    source_intervals::NTuple{3,UnitRange{Int}},
    source_mode_dims::NTuple{3,Int},
)
    matrices = _axis_transform_tuple(axis_transform_matrices)
    return ntuple(
        axis -> axis_source_transform_fact(
            axis,
            source_intervals[axis],
            source_mode_dims[axis],
            matrices[axis],
        ),
        3,
    )
end

function _validated_axis_transform_facts(
    axis_transform_facts,
    source_intervals::NTuple{3,UnitRange{Int}},
    source_mode_dims::NTuple{3,Int},
)
    facts = _axis_transform_tuple(axis_transform_facts)
    return ntuple(
        axis -> _validated_axis_transform_fact(
            facts[axis],
            axis,
            source_intervals[axis],
            source_mode_dims[axis],
        ),
        3,
    )
end

function _validated_axis_transform_fact(
    fact,
    axis::Int,
    source_interval::UnitRange{Int},
    source_mode_dim::Int,
)
    fact isa AxisSourceTransformFact ||
        throw(ArgumentError("axis_transform_facts entries must be AxisSourceTransformFact"))
    fact.axis == axis ||
        throw(ArgumentError("axis_transform_facts axis ordering mismatch"))
    fact.source_interval == source_interval ||
        throw(ArgumentError("axis_transform_facts source interval mismatch"))
    fact.source_mode_dim == source_mode_dim ||
        throw(ArgumentError("axis_transform_facts source mode dimension mismatch"))
    if fact.coefficient_status === :materialized
        _validated_axis_transform_matrix(
            fact.coefficient_matrix,
            source_interval,
            source_mode_dim,
            axis,
        )
    elseif fact.coefficient_status === :not_materialized
        isnothing(fact.coefficient_matrix) ||
            throw(ArgumentError("not-materialized axis transform facts cannot carry a coefficient matrix"))
    else
        throw(ArgumentError("unsupported axis transform coefficient status"))
    end
    return fact
end

function _axis_transform_tuple(axis_transforms)
    if axis_transforms isa NamedTuple
        all(name -> haskey(axis_transforms, name), (:x, :y, :z)) ||
            throw(ArgumentError("axis transforms named tuple must contain x, y, and z"))
        return (
            getfield(axis_transforms, :x),
            getfield(axis_transforms, :y),
            getfield(axis_transforms, :z),
        )
    end
    axis_transforms isa Tuple && length(axis_transforms) == 3 ||
        throw(ArgumentError("axis transforms must be a 3-tuple or named tuple with x/y/z"))
    return axis_transforms
end

function _validated_axis_transform_matrix(
    coefficient_matrix,
    source_interval::UnitRange{Int},
    source_mode_dim::Int,
    axis::Int,
)
    coefficient_matrix isa AbstractMatrix{<:Real} ||
        throw(ArgumentError("source-axis transform coefficient matrix must be real"))
    expected_shape = (length(source_interval), source_mode_dim)
    size(coefficient_matrix) == expected_shape ||
        throw(
            ArgumentError(
                "source-axis transform matrix for axis $(axis) must have shape $(expected_shape)",
            ),
        )
    return Matrix{Float64}(coefficient_matrix)
end
