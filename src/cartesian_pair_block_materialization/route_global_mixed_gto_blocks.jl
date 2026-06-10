# Route-global mixed gausslet/GTO block placement over decomposed WL retained
# units. This contracts provider-level CPB/GTO local rows with existing
# White--Lindsey unit coefficient maps and places the retained rows into the
# route-owned gausslet retained range. It does not assemble a Hamiltonian, use a
# direct Cartesian fallback, construct a full-parent CPB, or apply PQS
# transforms.

function route_global_mixed_gto_blocks_from_decomposed_units(
    inventory,
    parent,
    supplement;
    expansion = nothing,
    center_records = (),
    metadata = (;),
)
    blocker = _route_global_mixed_gto_blocks_input_blocker(inventory, supplement)
    if !isnothing(blocker)
        return _route_global_mixed_gto_blocks_result(
            :blocked_route_global_mixed_gto_blocks,
            blocker,
            inventory,
            supplement,
            nothing,
            nothing,
            (),
            nothing;
            metadata,
        )
    end

    provider_module = _route_global_mixed_gto_provider_module()
    isnothing(provider_module) &&
        return _route_global_mixed_gto_blocks_result(
            :blocked_route_global_mixed_gto_blocks,
            :missing_cpb_gto_provider_module,
            inventory,
            supplement,
            nothing,
            nothing,
            (),
            nothing;
            metadata,
        )

    retained_dimension = Int(inventory.retained_dimension)
    gto_count = length(supplement.orbitals)
    center_records_tuple = Tuple(center_records)
    mixed_overlap = zeros(Float64, retained_dimension, gto_count)
    mixed_kinetic = zeros(Float64, retained_dimension, gto_count)
    mixed_nuclear = Tuple(
        zeros(Float64, retained_dimension, gto_count) for _ in center_records_tuple
    )
    covered = falses(retained_dimension)
    placed_unit_keys = Symbol[]
    direct_core_unit_count = 0
    boundary_unit_count = 0
    gto_blocks = nothing
    gto_nuclear_blocks = ()

    for unit in inventory.retained_units
        unit_result = _route_global_mixed_gto_unit_blocks(
            unit,
            parent,
            supplement,
            provider_module;
            expansion,
            center_records = center_records_tuple,
        )
        unit_result.status === :materialized_route_global_mixed_gto_unit_blocks ||
            return _route_global_mixed_gto_blocks_result(
                :blocked_route_global_mixed_gto_blocks,
                unit_result.blocker,
                inventory,
                supplement,
                nothing,
                nothing,
                (),
                nothing;
                metadata = merge(
                    NamedTuple(metadata),
                    (;
                        blocked_unit_key = unit_result.unit_key,
                        blocked_unit_stratum_kind = unit_result.stratum_kind,
                    ),
                ),
            )

        row_range = unit_result.column_range
        any(covered[row_range]) &&
            return _route_global_mixed_gto_blocks_result(
                :blocked_route_global_mixed_gto_blocks,
                :duplicate_mixed_gto_route_global_row_coverage,
                inventory,
                supplement,
                nothing,
                nothing,
                (),
                nothing;
                metadata,
            )
        covered[row_range] .= true
        mixed_overlap[row_range, :] .= unit_result.mixed_overlap
        mixed_kinetic[row_range, :] .= unit_result.mixed_kinetic
        for index in eachindex(mixed_nuclear)
            mixed_nuclear[index][row_range, :] .= unit_result.mixed_nuclear[index]
        end
        push!(placed_unit_keys, unit_result.unit_key)
        unit_result.stratum_kind === :direct_core ?
            (direct_core_unit_count += 1) :
            (boundary_unit_count += 1)
        if isnothing(gto_blocks)
            gto_blocks = unit_result.gto_blocks
            gto_nuclear_blocks = unit_result.gto_nuclear_blocks
        end
    end

    all(covered) ||
        return _route_global_mixed_gto_blocks_result(
            :blocked_route_global_mixed_gto_blocks,
            :missing_mixed_gto_route_global_row_coverage,
            inventory,
            supplement,
            (overlap = mixed_overlap, kinetic = mixed_kinetic),
            gto_blocks,
            mixed_nuclear,
            gto_nuclear_blocks;
            covered,
            placed_unit_keys = Tuple(placed_unit_keys),
            direct_core_unit_count,
            boundary_unit_count,
            metadata,
        )

    return _route_global_mixed_gto_blocks_result(
        :materialized_route_global_mixed_gto_blocks,
        nothing,
        inventory,
        supplement,
        (overlap = mixed_overlap, kinetic = mixed_kinetic),
        gto_blocks,
        mixed_nuclear,
        gto_nuclear_blocks;
        covered,
        placed_unit_keys = Tuple(placed_unit_keys),
        direct_core_unit_count,
        boundary_unit_count,
        metadata,
    )
end

function _route_global_mixed_gto_blocks_input_blocker(inventory, supplement)
    _route_global_combined_gto_property(inventory, :status, nothing) ===
        :available_white_lindsey_decomposed_unit_pair_inventory ||
        return :missing_decomposed_wl_unit_pair_inventory
    retained_dimension =
        _route_global_combined_gto_property(inventory, :retained_dimension, nothing)
    retained_dimension isa Integer && retained_dimension > 0 ||
        return :missing_decomposed_wl_retained_dimension
    units = _route_global_combined_gto_property(inventory, :retained_units, ())
    units isa Tuple && !isempty(units) ||
        return :missing_decomposed_wl_retained_units
    _route_global_combined_gto_has_orbitals(supplement) ||
        return :missing_gto_supplement_orbitals
    isempty(supplement.orbitals) && return :empty_gto_supplement_orbitals
    return nothing
end

function _route_global_mixed_gto_provider_module()
    parent = Base.parentmodule(@__MODULE__)
    isdefined(parent, :CartesianCPBBlockProviders) || return nothing
    return getproperty(parent, :CartesianCPBBlockProviders)
end

function _route_global_mixed_gto_unit_blocks(
    unit,
    parent,
    supplement,
    provider_module;
    expansion,
    center_records,
)
    column_range = _route_global_combined_gto_property(unit, :column_range, nothing)
    column_range isa AbstractUnitRange{<:Integer} ||
        return _route_global_mixed_gto_unit_result(
            :blocked_route_global_mixed_gto_unit_blocks,
            :missing_retained_unit_column_range,
            unit,
            column_range,
            nothing,
            nothing,
            (),
            nothing,
        )
    source_cpbs = _route_global_combined_gto_property(unit, :source_cpbs, ())
    length(source_cpbs) == 1 ||
        return _route_global_mixed_gto_unit_result(
            :blocked_route_global_mixed_gto_unit_blocks,
            :missing_decomposed_wl_mixed_gto_unit_cpb_source,
            unit,
            column_range,
            nothing,
            nothing,
            (),
            nothing,
        )
    source_cpb = only(source_cpbs)

    coefficients = white_lindsey_boundary_stratum_unit_coefficients(unit)
    coefficients.status ===
        :blocked_white_lindsey_boundary_stratum_unit_coefficients &&
        return _route_global_mixed_gto_unit_result(
            :blocked_route_global_mixed_gto_unit_blocks,
            coefficients.blocker,
            unit,
            column_range,
            nothing,
            nothing,
            (),
            nothing,
        )
    coefficients.coefficient_maps_materialized ||
        return _route_global_mixed_gto_unit_result(
            :blocked_route_global_mixed_gto_unit_blocks,
            :missing_white_lindsey_unit_coefficient_map,
            unit,
            column_range,
            nothing,
            nothing,
            (),
            nothing,
        )

    bundle = provider_module.cpb_gto_supplement_local_operator_bundle(
        parent,
        source_cpb,
        supplement;
        expansion,
        center_records,
    )
    bundle_summary = provider_module.summary(bundle)
    bundle_summary.status ===
        :materialized_cpb_gto_supplement_local_operator_bundle ||
        return _route_global_mixed_gto_unit_result(
            :blocked_route_global_mixed_gto_unit_blocks,
            :mixed_gto_unit_bundle_blocked,
            unit,
            column_range,
            nothing,
            nothing,
            (),
            nothing,
        )

    support_coefficients = _white_lindsey_support_coefficient_matrix(
        coefficients.coefficient_matrix,
        coefficients.support_indices,
        coefficients.coefficient_space,
        :mixed_gto,
    )
    mixed_overlap = _route_global_mixed_gto_contract_unit_block(
        support_coefficients,
        bundle.mixed_blocks.overlap.dense_block,
        column_range,
        :mixed_gto_overlap_dimension_mismatch,
        unit,
    )
    mixed_kinetic = _route_global_mixed_gto_contract_unit_block(
        support_coefficients,
        bundle.mixed_blocks.kinetic.dense_block,
        column_range,
        :mixed_gto_kinetic_dimension_mismatch,
        unit,
    )
    mixed_nuclear = Tuple(
        _route_global_mixed_gto_contract_unit_block(
            support_coefficients,
            block.dense_block,
            column_range,
            :mixed_gto_electron_nuclear_dimension_mismatch,
            unit,
        ) for block in bundle.mixed_nuclear_by_center_blocks
    )
    any(block -> block isa Symbol, (mixed_overlap, mixed_kinetic)) &&
        return _route_global_mixed_gto_unit_result(
            :blocked_route_global_mixed_gto_unit_blocks,
            mixed_overlap isa Symbol ? mixed_overlap : mixed_kinetic,
            unit,
            column_range,
            nothing,
            nothing,
            (),
            nothing,
        )
    nuclear_blocker = findfirst(block -> block isa Symbol, mixed_nuclear)
    isnothing(nuclear_blocker) ||
        return _route_global_mixed_gto_unit_result(
            :blocked_route_global_mixed_gto_unit_blocks,
            mixed_nuclear[nuclear_blocker],
            unit,
            column_range,
            nothing,
            nothing,
            (),
            nothing,
        )

    return _route_global_mixed_gto_unit_result(
        :materialized_route_global_mixed_gto_unit_blocks,
        nothing,
        unit,
        column_range,
        mixed_overlap,
        mixed_kinetic,
        mixed_nuclear,
        bundle.gto_blocks;
        gto_nuclear_blocks = bundle.gto_nuclear_by_center_blocks,
    )
end

function _route_global_mixed_gto_contract_unit_block(
    support_coefficients,
    local_block,
    column_range,
    blocker::Symbol,
    unit,
)
    local_block isa AbstractMatrix || return blocker
    size(local_block, 1) == size(support_coefficients, 1) || return blocker
    size(support_coefficients, 2) == length(column_range) || return blocker
    return Matrix{Float64}(transpose(support_coefficients) * local_block)
end

function _route_global_mixed_gto_unit_result(
    status::Symbol,
    blocker,
    unit,
    column_range,
    mixed_overlap,
    mixed_kinetic,
    mixed_nuclear,
    gto_blocks;
    gto_nuclear_blocks = (),
)
    return (;
        object_kind = :route_global_mixed_gto_unit_blocks,
        status,
        blocker,
        unit_key = _route_global_combined_gto_property(unit, :unit_key, :unavailable),
        stratum_kind =
            _route_global_combined_gto_nested_property(
                unit,
                (:metadata, :stratum_kind),
                :unavailable,
            ),
        column_range,
        mixed_overlap,
        mixed_kinetic,
        mixed_nuclear,
        gto_blocks,
        gto_nuclear_blocks,
    )
end

function _route_global_mixed_gto_blocks_result(
    status::Symbol,
    blocker,
    inventory,
    supplement,
    mixed_blocks,
    gto_blocks,
    mixed_nuclear_by_center_matrices,
    gto_nuclear_by_center_blocks;
    covered = falses(
        _route_global_combined_gto_property(inventory, :retained_dimension, 0),
    ),
    placed_unit_keys = (),
    direct_core_unit_count = 0,
    boundary_unit_count = 0,
    metadata = (;),
)
    materialized = status === :materialized_route_global_mixed_gto_blocks
    row_range = _route_global_mixed_gto_row_range(covered)
    mixed_nuclear_blocks = Tuple(
        (;
            dense_block = matrix,
            by_center = true,
            nuclear_charge_applied = false,
            centers_summed = false,
        ) for matrix in mixed_nuclear_by_center_matrices
    )
    return (;
        object_kind = :route_global_mixed_gto_blocks_from_decomposed_units,
        status,
        blocker,
        mixed_blocks =
            isnothing(mixed_blocks) ? (;) :
            (;
                overlap = (; dense_block = mixed_blocks.overlap),
                kinetic = (; dense_block = mixed_blocks.kinetic),
            ),
        gto_blocks = isnothing(gto_blocks) ? (;) : gto_blocks,
        mixed_nuclear_by_center_blocks = mixed_nuclear_blocks,
        gto_nuclear_by_center_blocks =
            isnothing(gto_nuclear_by_center_blocks) ? () :
            gto_nuclear_by_center_blocks,
        decomposed_inventory_status =
            _route_global_combined_gto_property(inventory, :status, :unavailable),
        decomposed_inventory_source_kind =
            _route_global_combined_gto_property(inventory, :source_kind, :unavailable),
        retained_dimension =
            _route_global_combined_gto_property(inventory, :retained_dimension, :unavailable),
        gto_supplement_orbital_count =
            _route_global_combined_gto_orbital_count(supplement),
        mixed_gausslet_row_range = row_range,
        mixed_gausslet_row_count =
            row_range isa AbstractUnitRange ? length(row_range) : 0,
        mixed_gausslet_row_coverage_status =
            materialized ?
            :full_mixed_gto_route_global_row_coverage :
            :partial_mixed_gto_route_global_row_coverage,
        missing_mixed_gausslet_row_ranges =
            _route_global_mixed_gto_missing_row_ranges(covered),
        placed_unit_count = length(placed_unit_keys),
        placed_unit_keys,
        direct_core_unit_count,
        boundary_unit_count,
        direct_core_units_represented = direct_core_unit_count > 0,
        boundary_units_represented = boundary_unit_count > 0,
        route_global_mixed_gto_blocks_materialized = materialized,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        pqs_transforms_materialized = false,
        exports_or_artifacts = false,
        metadata = NamedTuple(metadata),
    )
end

function _route_global_mixed_gto_row_range(covered::AbstractVector{Bool})
    any(covered) || return :unavailable
    first_row = findfirst(identity, covered)
    last_row = findlast(identity, covered)
    return first_row:last_row
end

function _route_global_mixed_gto_missing_row_ranges(covered::AbstractVector{Bool})
    missing = UnitRange{Int}[]
    index = firstindex(covered)
    while index <= lastindex(covered)
        if covered[index]
            index += 1
            continue
        end
        start = index
        while index <= lastindex(covered) && !covered[index]
            index += 1
        end
        push!(missing, start:(index - 1))
    end
    return Tuple(missing)
end
