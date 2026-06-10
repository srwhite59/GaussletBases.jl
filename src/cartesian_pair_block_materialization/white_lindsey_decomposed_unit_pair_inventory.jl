# Route-owned decomposed White-Lindsey unit-pair inventory.
#
# This validates retained-unit pair metadata that later acceptance assembly can
# consume. It does not build local blocks, assign ranges, assemble matrices, or
# call the legacy direct Cartesian operator path.

const _WHITE_LINDSEY_ACCEPTANCE_ONE_BODY_TERMS = (
    :overlap,
    :kinetic,
    :electron_nuclear_by_center,
)

function white_lindsey_decomposed_unit_pair_inventory(
    plan::CUP.UnitPairPlan;
    supported_terms = _WHITE_LINDSEY_ACCEPTANCE_ONE_BODY_TERMS,
    metadata = (;),
)
    return white_lindsey_decomposed_unit_pair_inventory(
        CUP.unit_pairs(plan);
        source_kind = :cartesian_unit_pair_plan,
        supported_terms,
        metadata,
    )
end

function white_lindsey_decomposed_unit_pair_inventory(
    pairs::AbstractVector{<:CUP.UnitPairRecord};
    kwargs...,
)
    return white_lindsey_decomposed_unit_pair_inventory(Tuple(pairs); kwargs...)
end

function white_lindsey_decomposed_unit_pair_inventory(
    pairs::Tuple{Vararg{CUP.UnitPairRecord}};
    source_kind::Symbol = :unit_pair_records,
    supported_terms = _WHITE_LINDSEY_ACCEPTANCE_ONE_BODY_TERMS,
    metadata = (;),
)
    term_tuple = Tuple(supported_terms)
    if isempty(pairs)
        return _white_lindsey_decomposed_unit_pair_inventory_result(
            :blocked_white_lindsey_decomposed_unit_pair_inventory,
            :empty_decomposed_wl_unit_pair_inventory,
            ();
            source_kind,
            supported_terms = term_tuple,
            metadata,
        )
    end

    units = _white_lindsey_decomposed_inventory_units(pairs)
    blocker = _white_lindsey_decomposed_inventory_blocker(units)
    status =
        isnothing(blocker) ?
        :available_white_lindsey_decomposed_unit_pair_inventory :
        :blocked_white_lindsey_decomposed_unit_pair_inventory

    return _white_lindsey_decomposed_unit_pair_inventory_result(
        status,
        blocker,
        pairs;
        source_kind,
        supported_terms = term_tuple,
        metadata,
        units,
    )
end

function white_lindsey_decomposed_unit_pair_inventory(
    source;
    source_kind::Symbol = :missing_unit_pair_inventory,
    supported_terms = _WHITE_LINDSEY_ACCEPTANCE_ONE_BODY_TERMS,
    metadata = (;),
)
    isnothing(source) ||
        return _white_lindsey_decomposed_unit_pair_inventory_result(
            :blocked_white_lindsey_decomposed_unit_pair_inventory,
            :unsupported_decomposed_wl_unit_pair_inventory_source,
            ();
            source_kind,
            supported_terms = Tuple(supported_terms),
            metadata,
        )
    return _white_lindsey_decomposed_unit_pair_inventory_result(
        :blocked_white_lindsey_decomposed_unit_pair_inventory,
        :missing_decomposed_wl_unit_pair_inventory_source,
        ();
        source_kind,
        supported_terms = Tuple(supported_terms),
        metadata,
    )
end

function _white_lindsey_decomposed_inventory_units(pairs)
    by_key = Dict{Symbol,Any}()
    ordered = Any[]
    for pair in pairs
        for unit in (pair.left_unit, pair.right_unit)
            if !haskey(by_key, unit.unit_key)
                by_key[unit.unit_key] = unit
                push!(ordered, unit)
            end
        end
    end
    return Tuple(ordered)
end

function _white_lindsey_decomposed_inventory_blocker(units)
    all(
        unit -> unit.unit_kind === :white_lindsey_boundary_stratum_retained_unit,
        units,
    ) || return :unsupported_decomposed_wl_unit_kind
    all(unit -> unit.dimension_status === :available && !isnothing(unit.dimension), units) ||
        return :missing_retained_unit_dimension
    all(unit -> unit.column_range_status === :available && !isnothing(unit.column_range), units) ||
        return :missing_retained_unit_column_ranges
    all(
        unit -> length(unit.column_range) == unit.dimension,
        units,
    ) || return :retained_unit_range_dimension_mismatch
    return nothing
end

function _white_lindsey_decomposed_inventory_global_dimension(units, blocker)
    !isnothing(blocker) && return nothing
    isempty(units) && return nothing
    return maximum(last(unit.column_range) for unit in units)
end

function _white_lindsey_decomposed_inventory_unit_summaries(units)
    return Tuple(
        (;
            unit_key = unit.unit_key,
            unit_index = unit.unit_index,
            unit_kind = unit.unit_kind,
            stratum_kind =
                _white_lindsey_unit_metadata_value(unit, :stratum_kind),
            source_cpb_index = unit.source_cpb_index,
            source_cpb_count = length(unit.source_cpbs),
            dimension_status = unit.dimension_status,
            dimension = unit.dimension,
            column_range_status = unit.column_range_status,
            column_range = unit.column_range,
        ) for unit in units
    )
end

function _white_lindsey_decomposed_inventory_pair_summaries(pairs, blocker)
    return Tuple(
        (;
            pair_key = pair.pair_key,
            pair_index = pair.pair_index,
            pair_family = pair.pair_family,
            left_unit_key = pair.left_unit_key,
            right_unit_key = pair.right_unit_key,
            left_column_range =
                isnothing(blocker) ? pair.left_unit.column_range : nothing,
            right_column_range =
                isnothing(blocker) ? pair.right_unit.column_range : nothing,
            left_dimension =
                isnothing(blocker) ? pair.left_unit.dimension : nothing,
            right_dimension =
                isnothing(blocker) ? pair.right_unit.dimension : nothing,
        ) for pair in pairs
    )
end

function _white_lindsey_decomposed_unit_pair_inventory_result(
    status::Symbol,
    blocker,
    pairs;
    source_kind::Symbol,
    supported_terms,
    metadata,
    units = (),
)
    available = status === :available_white_lindsey_decomposed_unit_pair_inventory
    unit_summaries = _white_lindsey_decomposed_inventory_unit_summaries(units)
    pair_summaries =
        _white_lindsey_decomposed_inventory_pair_summaries(pairs, blocker)
    retained_dimension =
        _white_lindsey_decomposed_inventory_global_dimension(units, blocker)

    return (;
        object_kind = :white_lindsey_decomposed_unit_pair_inventory,
        status,
        blocker,
        source_kind,
        supported_terms,
        term_compatibility = (;
            overlap = :overlap in supported_terms,
            kinetic = :kinetic in supported_terms,
            electron_nuclear_by_center =
                :electron_nuclear_by_center in supported_terms,
        ),
        unit_count = length(units),
        pair_count = length(pairs),
        unit_keys = Tuple(unit.unit_key for unit in units),
        pair_keys = Tuple(pair.pair_key for pair in pairs),
        unit_summaries,
        pair_summaries,
        retained_dimension,
        retained_dimension_status =
            available ?
            :available_from_decomposed_wl_unit_column_ranges :
            :not_available,
        retained_unit_column_ranges_materialized = available,
        decomposed_unit_pair_column_ranges_available = available,
        decomposed_wl_unit_pair_inventory_available = available,
        route_owned_inventory = true,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        local_operator_blocks_materialized = false,
        global_matrices_materialized = false,
        hamiltonian_data_materialized = false,
        route_driver_wiring = false,
        metadata = NamedTuple(metadata),
    )
end
