# Intermediate/final retained spaces and shell-realization records.

function _positive_or_nothing(value, label::AbstractString)
    isnothing(value) && return nothing
    value isa Integer && value > 0 ||
        throw(ArgumentError("$label must be a positive integer or nothing"))
    return Int(value)
end

function _source_mode_dims_or_nothing(source_mode_dims)
    isnothing(source_mode_dims) && return nothing
    length(source_mode_dims) == 3 ||
        throw(ArgumentError("source_mode_dims must have length 3"))
    dims = (Int(source_mode_dims[1]), Int(source_mode_dims[2]), Int(source_mode_dims[3]))
    all(>(0), dims) || throw(ArgumentError("source_mode_dims must be positive"))
    return dims
end

"""
    boundary_product_mode_count(source_mode_dims)

Count product modes with at least one boundary index along a 3D product source.

For source dimensions `(nx, ny, nz)`, this returns the full product count minus
the strictly interior product count. This is the PQS counting rule for boundary
COMX-product mode selection.
"""
function boundary_product_mode_count(source_mode_dims)
    dims = _source_mode_dims_or_nothing(source_mode_dims)
    total = prod(dims)
    interior = prod(max(dim - 2, 0) for dim in dims)
    return total - interior
end

"""
    IntermediateRetainedSpace

Planned or materialized retained source-space object downstream of a
`LoweringSource` and upstream of shell realization.
"""
struct IntermediateRetainedSpace
    space_kind::Symbol
    lowering_source::LoweringSource
    retained_rule::Symbol
    dimension::Union{Int,Nothing}
    source_mode_dims::Union{NTuple{3,Int},Nothing}
    materialized::Bool
    metadata::NamedTuple
end

"""
    intermediate_retained_space(lowering_source;
                                retained_rule,
                                dimension = nothing,
                                source_mode_dims = nothing,
                                materialized = false,
                                metadata = (;))

Represent a retained source-space object before final shell realization.

For PQS, this is where boundary COMX-product modes live. For direct pieces, this
may be an identity-like source-mode space. This object is metadata/planning
unless `materialized = true`.
"""
function intermediate_retained_space(
    lowering_source::LoweringSource;
    space_kind::Symbol = :intermediate_retained_space,
    retained_rule::Symbol,
    dimension = nothing,
    source_mode_dims = nothing,
    materialized::Bool = false,
    metadata = (;),
)
    dims = _source_mode_dims_or_nothing(source_mode_dims)
    computed_dimension =
        !isnothing(dimension) ? _positive_or_nothing(dimension, "dimension") :
        retained_rule === :pqs_boundary_comx_product_modes && !isnothing(dims) ?
        boundary_product_mode_count(dims) :
        retained_rule === :identity_source_modes && !isnothing(dims) ?
        prod(dims) :
        nothing

    return IntermediateRetainedSpace(
        space_kind,
        lowering_source,
        retained_rule,
        computed_dimension,
        dims,
        materialized,
        NamedTuple(metadata),
    )
end

"""
    ShellRealization

Plan or record for realizing an intermediate retained space on
shellification-owned support.
"""
struct ShellRealization
    realization_kind::Symbol
    intermediate::IntermediateRetainedSpace
    owned_region::ShellificationRegion
    status::Symbol
    final_dimension::Union{Int,Nothing}
    metadata::NamedTuple
end

"""
    shell_realization(intermediate, owned_region; realization_kind,
                      status = :planned,
                      final_dimension = intermediate.dimension,
                      metadata = (;))

Represent how an intermediate retained space is realized on shellification-owned
support.

This records the realization contract only. It does not build projection,
Lowdin, or embedding matrices.
"""
function shell_realization(
    intermediate::IntermediateRetainedSpace,
    owned_region::ShellificationRegion;
    realization_kind::Symbol,
    status::Symbol = :planned,
    final_dimension = intermediate.dimension,
    metadata = (;),
)
    dim = _positive_or_nothing(final_dimension, "final_dimension")
    return ShellRealization(
        realization_kind,
        intermediate,
        owned_region,
        status,
        dim,
        NamedTuple(metadata),
    )
end

"""
    trivial_shell_realization(intermediate, owned_region;
                              status = :direct_or_trivial,
                              final_dimension = intermediate.dimension,
                              metadata = (;))

Construct a direct/trivial shell realization.

Use this for direct or identity-like retained spaces where no PQS shell
projection/Lowdin cleanup is planned.
"""
function trivial_shell_realization(
    intermediate::IntermediateRetainedSpace,
    owned_region::ShellificationRegion;
    status::Symbol = :direct_or_trivial,
    final_dimension = intermediate.dimension,
    metadata = (;),
)
    return shell_realization(
        intermediate,
        owned_region;
        realization_kind = :direct_or_trivial_embedding,
        status,
        final_dimension,
        metadata,
    )
end

"""
    pqs_shell_realization(intermediate, owned_region;
                          status = :planned_shell_projection_lowdin,
                          final_dimension = intermediate.dimension,
                          metadata = (;))

Construct a PQS shell-realization plan.

The intermediate space must use `:pqs_boundary_comx_product_modes`. This records
that shell projection and Lowdin cleanup are planned; it does not materialize
those transforms.
"""
function pqs_shell_realization(
    intermediate::IntermediateRetainedSpace,
    owned_region::ShellificationRegion;
    status::Symbol = :planned_shell_projection_lowdin,
    final_dimension = intermediate.dimension,
    metadata = (;),
)
    intermediate.retained_rule === :pqs_boundary_comx_product_modes ||
        throw(ArgumentError("PQS shell realization expects a PQS boundary-mode intermediate space"))
    return shell_realization(
        intermediate,
        owned_region;
        realization_kind = :shell_projection_lowdin,
        status,
        final_dimension,
        metadata,
    )
end

"""
    FinalRetainedUnit

Column-owning retained unit used by pair planning. It preserves links back to
the lowering source, intermediate space, and shell realization.
"""
struct FinalRetainedUnit
    unit_key::Symbol
    role::Symbol
    lowering_source::LoweringSource
    intermediate::IntermediateRetainedSpace
    shell_realization::ShellRealization
    column_range::Union{UnitRange{Int},Nothing}
    dimension::Union{Int,Nothing}
    metadata::NamedTuple
end

function _column_range_or_nothing(column_range)
    isnothing(column_range) && return nothing
    column_range isa AbstractUnitRange{<:Integer} ||
        throw(ArgumentError("column_range must be an integer unit range or nothing"))
    isempty(column_range) && throw(ArgumentError("column_range must be nonempty"))
    return Int(first(column_range)):Int(last(column_range))
end

"""
    final_retained_unit(unit_key, role, lowering_source, intermediate, realization;
                        column_range = nothing, dimension = nothing,
                        metadata = (;))

Construct the column-owning unit used by pair planning.

A final retained unit must be downstream of a lowering source, an intermediate
retained space, and a shell realization. Pair inventories should be built from
these objects, not directly from shellification regions or CPBs.
"""
function final_retained_unit(
    unit_key::Symbol,
    role::Symbol,
    lowering_source::LoweringSource,
    intermediate::IntermediateRetainedSpace,
    shell_realization::ShellRealization;
    column_range = nothing,
    dimension = nothing,
    metadata = (;),
)
    intermediate.lowering_source === lowering_source ||
        throw(ArgumentError("FinalRetainedUnit intermediate does not belong to lowering_source"))
    shell_realization.intermediate === intermediate ||
        throw(ArgumentError("FinalRetainedUnit realization does not belong to intermediate space"))
    shell_realization.owned_region === lowering_source.owned_region ||
        throw(ArgumentError("FinalRetainedUnit realization region does not match lowering source"))

    range = _column_range_or_nothing(column_range)
    dim =
        !isnothing(dimension) ? _positive_or_nothing(dimension, "dimension") :
        !isnothing(range) ? length(range) :
        !isnothing(shell_realization.final_dimension) ? shell_realization.final_dimension :
        intermediate.dimension

    if !isnothing(range) && !isnothing(dim) && length(range) != dim
        throw(ArgumentError("FinalRetainedUnit column_range length does not match dimension"))
    end

    return FinalRetainedUnit(
        unit_key,
        role,
        lowering_source,
        intermediate,
        shell_realization,
        range,
        dim,
        NamedTuple(metadata),
    )
end

role(unit::FinalRetainedUnit) = unit.role
lowering_recipe(unit::FinalRetainedUnit) = lowering_recipe(unit.lowering_source)
source_cpbs(unit::FinalRetainedUnit) = source_cpbs(unit.lowering_source)
owned_support(unit::FinalRetainedUnit) = owned_support(unit.lowering_source)
