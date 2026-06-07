# Metadata-only raw product source records.

"""
    AxisSourceTransformFact

Metadata record for one source-axis transform fact.

This is a status/provenance record, not a numerical transform builder. In the
current raw-source fact layer, axis coefficient matrices are not materialized:
default facts use `coefficient_status = :not_materialized` and
`coefficient_matrix = nothing`.

Future adapters may attach externally built axis transforms, but that numerical
work belongs outside `CartesianRawProductSources`.
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

Metadata-only raw product-box source facts.

A `RawProductBoxPlan` records the filled/source CPB, explicit source-mode
dimensions, source-mode ordering, source-mode tuple order, and per-axis
transform facts for a raw product source space. Source-mode dimensions are
total source-mode lengths, not interior counts, and are not inferred from the
CPB shape.

This object owns no retained rule, no boundary-mode selection, no shell
projection/Lowdin realization, no IDA weights, and no pair/operator blocks.
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

"""
    source_cpb(plan)

Return the coordinate product box used as the raw product source support.
"""
source_cpb(plan::RawProductBoxPlan) = plan.source_cpb

"""
    source_mode_dims(plan)

Return the three total source-mode lengths attached to a raw product source
plan.
"""
source_mode_dims(plan::RawProductBoxPlan) = plan.source_mode_dims

"""
    source_mode_count(plan)

Return `prod(source_mode_dims(plan))`.
"""
source_mode_count(plan::RawProductBoxPlan) = plan.source_mode_count

"""
    source_mode_indices(plan)

Return the deterministic source-mode tuples for the plan's ordering.
"""
source_mode_indices(plan::RawProductBoxPlan) = plan.source_mode_indices

"""
    axis_transform_facts(plan)

Return the per-axis metadata-only transform facts. These are not numerical
axis transforms unless a later adapter explicitly marks them as such.
"""
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
