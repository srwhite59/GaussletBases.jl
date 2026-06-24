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
dimensions, source-mode ordering, deterministic source-mode order, and per-axis
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
    source_mode_indices::Vector{NTuple{3,Int}}
    source_mode_column_indices::Vector{Int}
    source_mode_ordering::Symbol
    axis_transform_facts::NTuple{3,AxisSourceTransformFact}
    materialized::Bool
    metadata::NamedTuple
end

"""
    PQSBoundaryProductModeRetainedRule

Metadata-only retained source-mode rule for the first PQS source-box route.

The rule selects boundary product modes from a raw product source space: a mode
is retained when at least one local axis mode is the first or last mode on that
axis. It records deterministic source-mode and column indices only. It does not
materialize shell rows, shell projection, Lowdin cleanup, coefficient matrices,
pair blocks, IDA data, Hamiltonians, or artifacts.
"""
struct PQSBoundaryProductModeRetainedRule
    source_key::Symbol
    source_mode_dims::NTuple{3,Int}
    source_mode_ordering::Symbol
    retained_rule_kind::Symbol
    retained_mode_indices::Vector{NTuple{3,Int}}
    retained_column_indices::Vector{Int}
    retained_count::Int
    transform_kind::Symbol
    shell_realization_materialized::Bool
    lowdin_cleanup_used::Bool
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
    retained_mode_indices(rule)

Return deterministic retained source-mode tuples for a retained rule.
"""
retained_mode_indices(rule::PQSBoundaryProductModeRetainedRule) =
    rule.retained_mode_indices

"""
    retained_column_indices(rule)

Return source-mode column indices retained by a retained rule.
"""
retained_column_indices(rule::PQSBoundaryProductModeRetainedRule) =
    rule.retained_column_indices

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
    axis_transform_matrices = nothing,
    axis_transform_facts = nothing,
    metadata = (;),
)
    source_box isa CPB.CoordinateProductBox ||
        throw(ArgumentError("source_cpb must be a CartesianCPB.CoordinateProductBox"))
    normalized_dims = _normalize_source_mode_dims(source_mode_dims)
    indices = source_mode_indices(normalized_dims; source_mode_ordering)
    count = length(indices)
    intervals = CPB.intervals(source_box)
    transform_facts = _axis_transform_facts(
        intervals,
        normalized_dims;
        axis_transform_matrices,
        axis_transform_facts,
    )
    return RawProductBoxPlan(
        source_key,
        source_box,
        intervals,
        CPB.shape(source_box),
        normalized_dims,
        count,
        indices,
        collect(1:count),
        source_mode_ordering,
        transform_facts,
        any(fact -> fact.coefficient_status === :materialized, transform_facts),
        NamedTuple(metadata),
    )
end

function _is_boundary_product_mode(
    mode::NTuple{3,Int},
    dims::NTuple{3,Int},
)
    return any(axis -> mode[axis] == 1 || mode[axis] == dims[axis], 1:3)
end

"""
    pqs_boundary_product_mode_retained_rule(raw_plan)
    pqs_boundary_product_mode_retained_rule(source_mode_dims; ...)

Build the metadata-only PQS retained source-mode boundary rule.

The retained columns follow the same deterministic source-mode ordering used by
`source_mode_indices`. For `(5, 5, 5)`, the rule retains `98` boundary product
modes.
"""
function pqs_boundary_product_mode_retained_rule(
    source_mode_dims;
    source_mode_ordering::Symbol = _SUPPORTED_SOURCE_MODE_ORDERING,
    source_key::Symbol = :pqs_raw_product_source,
    metadata = (;),
)
    normalized_dims = _normalize_source_mode_dims(source_mode_dims)
    indices = source_mode_indices(
        normalized_dims;
        source_mode_ordering,
    )
    retained_modes = NTuple{3,Int}[]
    retained_columns = Int[]
    for (column_index, mode) in pairs(indices)
        if _is_boundary_product_mode(mode, normalized_dims)
            push!(retained_modes, mode)
            push!(retained_columns, column_index)
        end
    end
    return PQSBoundaryProductModeRetainedRule(
        source_key,
        normalized_dims,
        source_mode_ordering,
        :boundary_comx_product_mode_selection,
        retained_modes,
        retained_columns,
        length(retained_modes),
        :source_mode_column_selector,
        false,
        false,
        NamedTuple(metadata),
    )
end

function pqs_boundary_product_mode_retained_rule(
    raw_plan::RawProductBoxPlan;
    metadata = (;),
)
    return pqs_boundary_product_mode_retained_rule(
        raw_plan.source_mode_dims;
        source_mode_ordering = raw_plan.source_mode_ordering,
        source_key = raw_plan.source_key,
        metadata = _merge_rule_metadata(
            raw_plan.metadata,
            NamedTuple(metadata),
        ),
    )
end

function _merge_rule_metadata(left::NamedTuple, right::NamedTuple)
    return merge(left, right)
end
