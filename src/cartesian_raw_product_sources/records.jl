# Metadata-only raw product source records.

"""
    AxisSourceTransformFact

Metadata record for one source-axis transform. In this module's initial
contract pass, coefficient matrices are not materialized.
"""
struct AxisSourceTransformFact
    axis::Int
    source_interval::UnitRange{Int}
    source_mode_dim::Int
    coefficient_status::Symbol
    coefficient_matrix::Any
    metadata::NamedTuple
end

"""
    RawProductBoxPlan

Metadata-only raw product-box source facts. Source-mode dimensions are total
source-mode lengths, not interior counts.
"""
struct RawProductBoxPlan
    source_key::Symbol
    source_cpb::CPB.CoordinateProductBox
    source_intervals::NTuple{3,UnitRange{Int}}
    source_shape::NTuple{3,Int}
    source_mode_dims::NTuple{3,Int}
    source_mode_count::Int
    source_mode_indices::Tuple{Vararg{NTuple{3,Int}}}
    source_mode_column_indices::Tuple{Vararg{Int}}
    source_mode_ordering::Symbol
    axis_transform_facts::NTuple{3,AxisSourceTransformFact}
    materialized::Bool
    metadata::NamedTuple
end

source_cpb(plan::RawProductBoxPlan) = plan.source_cpb
source_mode_dims(plan::RawProductBoxPlan) = plan.source_mode_dims
source_mode_count(plan::RawProductBoxPlan) = plan.source_mode_count
source_mode_indices(plan::RawProductBoxPlan) = plan.source_mode_indices
axis_transform_facts(plan::RawProductBoxPlan) = plan.axis_transform_facts

function _normalize_source_mode_dims(dims)
    length(dims) == 3 ||
        throw(ArgumentError("source_mode_dims must contain exactly three dimensions"))
    normalized = map(dims) do dim
        dim isa Integer ||
            throw(ArgumentError("source_mode_dims entries must be integers"))
        dim > 0 ||
            throw(ArgumentError("source_mode_dims entries must be positive"))
        Int(dim)
    end
    return Tuple(normalized)::NTuple{3,Int}
end

function raw_product_box_plan(
    source_box;
    source_mode_dims,
    source_key::Symbol = :raw_product_source,
    source_mode_ordering::Symbol = :x_major_y_major_z_fast,
    metadata = (;),
)
    source_box isa CPB.CoordinateProductBox ||
        throw(ArgumentError("source_cpb must be a CartesianCPB.CoordinateProductBox"))
    normalized_dims = _normalize_source_mode_dims(source_mode_dims)
    indices = source_mode_indices(normalized_dims; source_mode_ordering)
    count = length(indices)
    intervals = CPB.intervals(source_box)
    return RawProductBoxPlan(
        source_key,
        source_box,
        intervals,
        CPB.shape(source_box),
        normalized_dims,
        count,
        indices,
        Tuple(1:count),
        source_mode_ordering,
        _default_axis_transform_facts(intervals, normalized_dims),
        false,
        NamedTuple(metadata),
    )
end
