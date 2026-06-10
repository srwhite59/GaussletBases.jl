# Metadata-only route-global layout for a decomposed WL gausslet sector plus a
# small Cartesian GTO supplement sector. This does not assemble overlap,
# Hamiltonian, route-driver data, exports, or artifacts.

function route_global_combined_gto_basis_layout(inventory, supplement)
    inventory_status = _route_global_combined_gto_property(
        inventory,
        :status,
        :missing_decomposed_wl_inventory,
    )
    if inventory_status !== :available_white_lindsey_decomposed_unit_pair_inventory
        return _route_global_combined_gto_basis_layout_result(
            :blocked_route_global_combined_gto_basis_layout,
            :missing_decomposed_wl_unit_pair_inventory,
            inventory,
            supplement,
        )
    end

    retained_dimension = _route_global_combined_gto_property(
        inventory,
        :retained_dimension,
        nothing,
    )
    retained_dimension isa Integer && retained_dimension > 0 ||
        return _route_global_combined_gto_basis_layout_result(
            :blocked_route_global_combined_gto_basis_layout,
            :missing_decomposed_wl_retained_dimension,
            inventory,
            supplement,
        )

    _route_global_combined_gto_has_orbitals(supplement) ||
        return _route_global_combined_gto_basis_layout_result(
            :blocked_route_global_combined_gto_basis_layout,
            :missing_gto_supplement_orbitals,
            inventory,
            supplement,
        )
    isempty(supplement.orbitals) &&
        return _route_global_combined_gto_basis_layout_result(
            :blocked_route_global_combined_gto_basis_layout,
            :empty_gto_supplement_orbitals,
            inventory,
            supplement,
        )

    return _route_global_combined_gto_basis_layout_result(
        :available_route_global_combined_gto_basis_layout,
        nothing,
        inventory,
        supplement,
    )
end

function _route_global_combined_gto_basis_layout_result(
    status::Symbol,
    blocker,
    inventory,
    supplement,
)
    available = status === :available_route_global_combined_gto_basis_layout
    retained_dimension = _route_global_combined_gto_property(
        inventory,
        :retained_dimension,
        :unavailable,
    )
    gausslet_dimension =
        available ? Int(retained_dimension) :
        retained_dimension isa Integer && retained_dimension > 0 ?
        Int(retained_dimension) :
        :unavailable
    gto_count = _route_global_combined_gto_orbital_count(supplement)
    gto_available = gto_count isa Integer && gto_count > 0

    gausslet_range =
        gausslet_dimension isa Integer ? (1:gausslet_dimension) : :unavailable
    gto_range =
        gausslet_dimension isa Integer && gto_available ?
        ((gausslet_dimension + 1):(gausslet_dimension + gto_count)) :
        :unavailable
    total_dimension =
        gausslet_dimension isa Integer && gto_available ?
        gausslet_dimension + gto_count :
        :unavailable
    block_layout =
        available ?
        _route_global_combined_gto_block_layout(gausslet_range, gto_range) :
        :unavailable

    return (;
        object_kind = :route_global_combined_gto_basis_layout,
        status,
        blocker,
        layout_kind = :decomposed_wl_gausslet_plus_gto_supplement,
        gausslet_basis_kind = :decomposed_white_lindsey_retained_gausslet,
        gto_basis_kind = :cartesian_gaussian_shell_supplement,
        gausslet_retained_range = gausslet_range,
        gausslet_retained_dimension = gausslet_dimension,
        gausslet_retained_dimension_source =
            gausslet_dimension isa Integer ?
            :decomposed_wl_unit_pair_inventory_retained_dimension :
            :unavailable,
        decomposed_unit_count =
            _route_global_combined_gto_property(inventory, :unit_count, :unavailable),
        decomposed_unit_pair_count =
            _route_global_combined_gto_property(inventory, :pair_count, :unavailable),
        decomposed_inventory_status =
            _route_global_combined_gto_property(inventory, :status, :unavailable),
        decomposed_inventory_source_kind =
            _route_global_combined_gto_property(inventory, :source_kind, :unavailable),
        gto_supplement_range = gto_range,
        gto_supplement_orbital_count = gto_count,
        gto_supplement_orbital_labels =
            _route_global_combined_gto_orbital_labels(supplement),
        gto_supplement_kind =
            _route_global_combined_gto_property(supplement, :supplement_kind, :unavailable),
        total_combined_dimension = total_dimension,
        combined_basis_dimension = total_dimension,
        block_layout_keys =
            available ?
            (:gausslet_gausslet, :gausslet_gto, :gto_gausslet, :gto_gto) :
            (),
        block_layout,
        gausslet_gausslet_block_layout =
            available ? block_layout.gausslet_gausslet : :unavailable,
        gausslet_gto_block_layout =
            available ? block_layout.gausslet_gto : :unavailable,
        gto_gausslet_block_layout =
            available ? block_layout.gto_gausslet : :unavailable,
        gto_gto_block_layout =
            available ? block_layout.gto_gto : :unavailable,
        mixed_cpb_gto_blocks_orientation = :gausslet_rows_by_gto_columns,
        gto_gto_blocks_orientation = :gto_rows_by_gto_columns,
        by_center_nuclear_blocks_separated = true,
        nuclear_charge_application_stage = :future_combined_hamiltonian_assembly,
        combined_overlap_layout_available = available,
        combined_hamiltonian_layout_available = available,
        combined_overlap_matrix_materialized = false,
        combined_hamiltonian_matrix_materialized = false,
        route_global_combined_matrix_materialized = false,
        gto_route_global_blocks_materialized = false,
        hamiltonian_assembly = false,
        hamiltonian_data_materialized = false,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        pqs_transforms_materialized = false,
        exports_or_artifacts = false,
    )
end

function _route_global_combined_gto_block_layout(gausslet_range, gto_range)
    return (;
        gausslet_gausslet = (;
            row_sector = :gausslet,
            column_sector = :gausslet,
            row_range = gausslet_range,
            column_range = gausslet_range,
            source = :decomposed_wl_route_global_gausslet_blocks,
        ),
        gausslet_gto = (;
            row_sector = :gausslet,
            column_sector = :gto,
            row_range = gausslet_range,
            column_range = gto_range,
            source = :provider_level_mixed_cpb_gto_blocks,
        ),
        gto_gausslet = (;
            row_sector = :gto,
            column_sector = :gausslet,
            row_range = gto_range,
            column_range = gausslet_range,
            source = :transpose_of_provider_level_mixed_cpb_gto_blocks,
        ),
        gto_gto = (;
            row_sector = :gto,
            column_sector = :gto,
            row_range = gto_range,
            column_range = gto_range,
            source = :provider_level_gto_gto_blocks,
        ),
    )
end

function _route_global_combined_gto_property(object, field::Symbol, default)
    return hasproperty(object, field) ? getproperty(object, field) : default
end

function _route_global_combined_gto_has_orbitals(supplement)
    return hasproperty(supplement, :orbitals)
end

function _route_global_combined_gto_orbital_count(supplement)
    _route_global_combined_gto_has_orbitals(supplement) || return :unavailable
    return length(supplement.orbitals)
end

function _route_global_combined_gto_orbital_labels(supplement)
    _route_global_combined_gto_has_orbitals(supplement) || return :unavailable
    return Tuple(
        hasproperty(orbital, :label) ? orbital.label : :unavailable
        for orbital in supplement.orbitals
    )
end
