# Private one-center atom + GTO supplement construction seam for the decomposed
# WL route. This wires existing route-global/final-basis pieces together from a
# mapped parent basis, shellification-backed retained-unit inventory, center
# records, and a GTO supplement representation. It does not solve RHF, add a
# public API, use the low-order seed path, or accept raw GTO density-density as
# final-basis electron-electron data.

function _white_lindsey_decomposed_atom_gto_final_basis_route(
    basis,
    supplement;
    center_records,
    nuclear_positions = Tuple(record.location for record in center_records),
    q::Integer = 5,
    ns::Integer = q,
    expansion = ParentGaussletBases.coulomb_gaussian_expansion(doacc = false),
    shellification_policy = CSH.OneCenterShellification(core_side = ns, q = q),
    lowering_policy = CTL.WhiteLindseyLowering(),
    gausslet_backend::Symbol = :numerical_reference,
    refinement_levels::Integer = 0,
    build_density_density::Bool = true,
    metadata = (;),
)
    axis_inputs = _white_lindsey_atom_gto_parent_axis_inputs(
        basis,
        expansion;
        gausslet_backend,
        refinement_levels,
    )
    parent = ParentGaussletBases.CartesianParentGaussletBases.CartesianParentGaussletBasis3D(
        basis,
    )
    supplement_representation = ParentGaussletBases.basis_representation(supplement)
    parent_axes = ntuple(_ -> ParentGaussletBases.centers(basis), 3)
    inventory_source = white_lindsey_shellification_decomposed_unit_pair_inventory(
        parent_axes,
        nuclear_positions;
        shellification_policy,
        lowering_policy,
        metadata = merge(
            NamedTuple(metadata),
            (;
                q,
                ns,
                route_seam =
                    :white_lindsey_decomposed_atom_gto_final_basis_route,
            ),
        ),
        parent_axis_counts = axis_inputs.parent_axis_counts,
        parent_axis_bundle_object = axis_inputs.parent_axis_bundle_object,
    )
    inventory = inventory_source.inventory
    inventory.status === :available_white_lindsey_decomposed_unit_pair_inventory ||
        return _white_lindsey_atom_gto_route_result(
            :blocked_decomposed_atom_gto_final_basis_route,
            inventory.blocker,
            axis_inputs,
            parent,
            supplement_representation,
            inventory_source,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing;
            metadata,
        )

    layout = route_global_combined_gto_basis_layout(inventory, supplement_representation)
    layout.status === :available_route_global_combined_gto_basis_layout ||
        return _white_lindsey_atom_gto_route_result(
            :blocked_decomposed_atom_gto_final_basis_route,
            layout.blocker,
            axis_inputs,
            parent,
            supplement_representation,
            inventory_source,
            layout,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing;
            metadata,
        )

    overlap = route_global_decomposed_wl_overlap_matrix(
        inventory;
        parent_axis_counts = axis_inputs.parent_axis_counts,
        parent_axis_bundle_object = axis_inputs.parent_axis_bundle_object,
        overlap_1d = axis_inputs.overlap_1d,
    )
    kinetic = route_global_decomposed_wl_kinetic_matrix(
        inventory;
        parent_axis_counts = axis_inputs.parent_axis_counts,
        parent_axis_bundle_object = axis_inputs.parent_axis_bundle_object,
        overlap_1d = axis_inputs.overlap_1d,
        kinetic_1d = axis_inputs.kinetic_1d,
    )
    nuclear = route_global_electron_nuclear_by_center_matrices(
        inventory;
        parent_axis_counts = axis_inputs.parent_axis_counts,
        parent_axis_bundle_object = axis_inputs.parent_axis_bundle_object,
        coulomb_expansion = expansion,
        center_records,
    )
    mixed = route_global_mixed_gto_blocks_from_decomposed_units(
        inventory,
        parent,
        supplement_representation;
        expansion,
        center_records,
        include_moment_blocks = build_density_density,
        metadata,
    )
    mixed.status === :materialized_route_global_mixed_gto_blocks ||
        return _white_lindsey_atom_gto_route_result(
            :blocked_decomposed_atom_gto_final_basis_route,
            mixed.blocker,
            axis_inputs,
            parent,
            supplement_representation,
            inventory_source,
            layout,
            overlap,
            kinetic,
            nuclear,
            mixed,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing;
            metadata,
        )

    combined = route_global_combined_gto_one_electron_matrices(
        layout;
        overlap_result = overlap,
        kinetic_result = kinetic,
        electron_nuclear_by_center_results = nuclear,
        gto_bundle = mixed,
        mixed_gausslet_row_range = mixed.mixed_gausslet_row_range,
        metadata,
    )
    combined.status === :materialized_route_global_combined_gto_one_electron_matrices ||
        return _white_lindsey_atom_gto_route_result(
            :blocked_decomposed_atom_gto_final_basis_route,
            combined.blocker,
            axis_inputs,
            parent,
            supplement_representation,
            inventory_source,
            layout,
            overlap,
            kinetic,
            nuclear,
            mixed,
            combined,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing;
            metadata,
        )

    projection = route_global_combined_gto_final_basis_projection(combined; metadata)
    projection.status === :materialized_route_global_combined_gto_final_basis_projection ||
        return _white_lindsey_atom_gto_route_result(
            :blocked_decomposed_atom_gto_final_basis_route,
            projection.blocker,
            axis_inputs,
            parent,
            supplement_representation,
            inventory_source,
            layout,
            overlap,
            kinetic,
            nuclear,
            mixed,
            combined,
            projection,
            nothing,
            nothing,
            nothing,
            nothing,
            nothing;
            metadata,
        )

    build_density_density || return _white_lindsey_atom_gto_route_result(
        :materialized_decomposed_atom_gto_final_basis_one_electron_route,
        nothing,
        axis_inputs,
        parent,
        supplement_representation,
        inventory_source,
        layout,
        overlap,
        kinetic,
        nuclear,
        mixed,
        combined,
        projection,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing;
        metadata,
    )

    position_x = route_global_decomposed_wl_position_x_matrix(
        inventory;
        parent_axis_counts = axis_inputs.parent_axis_counts,
        overlap_1d = axis_inputs.overlap_1d,
        position_1d = axis_inputs.position_1d,
    )
    position_y = route_global_decomposed_wl_position_y_matrix(
        inventory;
        parent_axis_counts = axis_inputs.parent_axis_counts,
        overlap_1d = axis_inputs.overlap_1d,
        position_1d = axis_inputs.position_1d,
    )
    position_z = route_global_decomposed_wl_position_z_matrix(
        inventory;
        parent_axis_counts = axis_inputs.parent_axis_counts,
        overlap_1d = axis_inputs.overlap_1d,
        position_1d = axis_inputs.position_1d,
    )
    x2_x = route_global_decomposed_wl_x2_x_matrix(
        inventory;
        parent_axis_counts = axis_inputs.parent_axis_counts,
        overlap_1d = axis_inputs.overlap_1d,
        x2_1d = axis_inputs.x2_1d,
    )
    x2_y = route_global_decomposed_wl_x2_y_matrix(
        inventory;
        parent_axis_counts = axis_inputs.parent_axis_counts,
        overlap_1d = axis_inputs.overlap_1d,
        x2_1d = axis_inputs.x2_1d,
    )
    x2_z = route_global_decomposed_wl_x2_z_matrix(
        inventory;
        parent_axis_counts = axis_inputs.parent_axis_counts,
        overlap_1d = axis_inputs.overlap_1d,
        x2_1d = axis_inputs.x2_1d,
    )
    moments = route_global_combined_gto_residual_moment_matrices(
        layout;
        position_x_result = position_x,
        position_y_result = position_y,
        position_z_result = position_z,
        x2_x_result = x2_x,
        x2_y_result = x2_y,
        x2_z_result = x2_z,
        gto_bundle = mixed,
        mixed_gausslet_row_range = mixed.mixed_gausslet_row_range,
        metadata,
    )
    moments.status === :materialized_route_global_combined_gto_residual_moment_matrices ||
        return _white_lindsey_atom_gto_route_result(
            :materialized_decomposed_atom_gto_final_basis_one_electron_route,
            moments.blocker,
            axis_inputs,
            parent,
            supplement_representation,
            inventory_source,
            layout,
            overlap,
            kinetic,
            nuclear,
            mixed,
            combined,
            projection,
            moments,
            nothing,
            nothing,
            nothing,
            nothing;
            metadata,
        )

    residual_mwg = route_global_residual_gto_mwg_representation(
        projection,
        combined,
        supplement_representation;
        raw_moment_matrices = moments,
        metadata,
    )
    residual_mwg.status === :materialized_route_global_residual_gto_mwg_representation ||
        return _white_lindsey_atom_gto_route_result(
            :materialized_decomposed_atom_gto_final_basis_one_electron_route,
            residual_mwg.blocker,
            axis_inputs,
            parent,
            supplement_representation,
            inventory_source,
            layout,
            overlap,
            kinetic,
            nuclear,
            mixed,
            combined,
            projection,
            moments,
            residual_mwg,
            nothing,
            nothing,
            nothing;
            metadata,
        )

    gausslet_density = route_global_decomposed_wl_density_density_matrix(
        inventory;
        parent_axis_counts = axis_inputs.parent_axis_counts,
        parent_axis_bundle_object = axis_inputs.parent_axis_bundle_object,
        coulomb_expansion = expansion,
    )
    gausslet_density.status === :materialized_route_global_density_density_interaction_matrix ||
        return _white_lindsey_atom_gto_route_result(
            :materialized_decomposed_atom_gto_final_basis_one_electron_route,
            gausslet_density.blocker,
            axis_inputs,
            parent,
            supplement_representation,
            inventory_source,
            layout,
            overlap,
            kinetic,
            nuclear,
            mixed,
            combined,
            projection,
            moments,
            residual_mwg,
            gausslet_density,
            nothing,
            nothing;
            metadata,
        )

    final_density = route_global_combined_gto_final_basis_density_density_matrix(
        combined,
        projection;
        gausslet_density_density_result = gausslet_density,
        residual_mwg_representation = residual_mwg,
        decomposed_inventory = inventory,
        parent_axis_counts = axis_inputs.parent_axis_counts,
        parent_axis_bundle_object = axis_inputs.parent_axis_bundle_object,
        coulomb_expansion = expansion,
        metadata,
    )
    status =
        final_density.status ===
        :materialized_route_global_combined_gto_final_basis_density_density_matrix ?
        :materialized_decomposed_atom_gto_final_basis_density_density_route :
        :materialized_decomposed_atom_gto_final_basis_one_electron_route
    blocker = status ===
        :materialized_decomposed_atom_gto_final_basis_density_density_route ?
        nothing : final_density.blocker
    return _white_lindsey_atom_gto_route_result(
        status,
        blocker,
        axis_inputs,
        parent,
        supplement_representation,
        inventory_source,
        layout,
        overlap,
        kinetic,
        nuclear,
        mixed,
        combined,
        projection,
        moments,
        residual_mwg,
        gausslet_density,
        final_density,
        nothing;
        metadata,
    )
end

function _white_lindsey_atom_gto_parent_axis_inputs(
    basis,
    expansion;
    gausslet_backend::Symbol,
    refinement_levels::Integer,
)
    source_1d = ParentGaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        center = 0.0,
        backend = gausslet_backend,
        refinement_levels,
    )
    pgdg = source_1d.pgdg_intermediate
    return (;
        parent_axis_counts = ntuple(_ -> length(ParentGaussletBases.centers(basis)), 3),
        parent_axis_bundle_object = (; x = source_1d, y = source_1d, z = source_1d),
        overlap_1d = (; x = pgdg.overlap, y = pgdg.overlap, z = pgdg.overlap),
        kinetic_1d = (; x = pgdg.kinetic, y = pgdg.kinetic, z = pgdg.kinetic),
        position_1d = (; x = pgdg.position, y = pgdg.position, z = pgdg.position),
        x2_1d = (; x = pgdg.x2, y = pgdg.x2, z = pgdg.x2),
    )
end

function _white_lindsey_atom_gto_route_result(
    status::Symbol,
    blocker,
    axis_inputs,
    parent,
    supplement,
    inventory_source,
    layout,
    overlap,
    kinetic,
    nuclear,
    mixed,
    combined,
    projection,
    moments,
    residual_mwg,
    gausslet_density,
    final_density,
    unused;
    metadata,
)
    inventory = isnothing(inventory_source) ? nothing : inventory_source.inventory
    final_dimension =
        _route_global_combined_gto_property(projection, :final_dimension, :unavailable)
    retained_dimension =
        _route_global_combined_gto_property(inventory, :retained_dimension, :unavailable)
    return (;
        object_kind = :white_lindsey_decomposed_atom_gto_final_basis_route,
        status,
        blocker,
        construction_seam =
            :white_lindsey_decomposed_atom_gto_final_basis_route,
        shellification_backed_decomposed_wl_inventory =
            _route_global_combined_gto_property(
                inventory_source,
                :shellification_backed_decomposed_wl_inventory,
                false,
            ),
        low_order_materialized_seed_inventory_used = false,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        raw_gto_density_density_accepted_as_final_operator = false,
        generalized_overlap_final_solve_used = false,
        parent_axis_counts = axis_inputs.parent_axis_counts,
        retained_gausslet_dimension = retained_dimension,
        unit_count = _route_global_combined_gto_property(inventory, :unit_count, :unavailable),
        pair_count = _route_global_combined_gto_property(inventory, :pair_count, :unavailable),
        raw_supplement_count =
            _route_global_combined_gto_orbital_count(supplement),
        retained_supplement_count =
            _route_global_combined_gto_property(
                projection,
                :retained_supplement_count,
                :unavailable,
            ),
        dropped_supplement_count =
            _route_global_combined_gto_property(
                projection,
                :dropped_supplement_count,
                :unavailable,
            ),
        final_dimension,
        final_overlap_identity_error =
            _route_global_combined_gto_property(
                projection,
                :final_overlap_identity_error,
                :unavailable,
            ),
        inventory_status =
            _route_global_combined_gto_property(inventory, :status, :unavailable),
        layout_status =
            _route_global_combined_gto_property(layout, :status, :unavailable),
        mixed_status =
            _route_global_combined_gto_property(mixed, :status, :unavailable),
        combined_one_electron_status =
            _route_global_combined_gto_property(combined, :status, :unavailable),
        final_basis_projection_status =
            _route_global_combined_gto_property(projection, :status, :unavailable),
        residual_moment_status =
            _route_global_combined_gto_property(moments, :status, :unavailable),
        residual_mwg_status =
            _route_global_combined_gto_property(residual_mwg, :status, :unavailable),
        gausslet_density_density_status =
            _route_global_combined_gto_property(
                gausslet_density,
                :status,
                :unavailable,
            ),
        final_density_density_status =
            _route_global_combined_gto_property(
                final_density,
                :status,
                :unavailable,
            ),
        parent,
        supplement,
        axis_inputs,
        inventory_source,
        inventory,
        layout,
        overlap,
        kinetic,
        electron_nuclear_by_center = nuclear,
        mixed_gto_blocks = mixed,
        combined_one_electron = combined,
        final_basis_projection = projection,
        residual_moment_matrices = moments,
        residual_mwg_representation = residual_mwg,
        gausslet_density_density = gausslet_density,
        final_density_density = final_density,
        metadata = NamedTuple(metadata),
    )
end
