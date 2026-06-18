# Optional package precompile workloads.
#
# The Cartesian/PQS driver workload is intentionally opt-in: it executes a real
# staged H2 PQS source-box route, which is useful for repeated ladder work but
# too expensive to impose on every package user during normal precompilation.

function _gaussletbases_precompile_flag(name::AbstractString)
    value = lowercase(get(ENV, name, "0"))
    return value in ("1", "true", "yes", "on")
end

function _cartesian_driver_precompile_inputs()
    system_inputs = (;
        atom_symbols = ("H", "H"),
        nuclear_charges = (1, 1),
        atom_locations = ((0.0, 0.0, -2.0), (0.0, 0.0, 2.0)),
        bond_axis = :z,
        bond_length = 4.0,
        radius = 4.0,
        parent_axis_counts = nothing,
        map_backend = :pgdg_localized_experimental,
    )
    spacing_inputs = (;
        q = 5,
        n_s = 5,
        reference_spacing = 1.0,
        tail_spacing = 10.0,
        q_to_core_spacing_rule = :standard_pqs_ns_equals_q,
        core_spacing = 0.5,
        xmax_parallel = 6.0,
        xmax_transverse = 4.0,
    )
    parent_inputs = (;
        parent_axis_bundle_backend = :pgdg_localized_experimental,
        parent_axis_family = :G10,
        parent_mapping_rule = :identity_mapping,
        parent_mapping_Z = nothing,
        parent_mapping_d = nothing,
        parent_mapping_tail_spacing = 10.0,
    )
    private_rhf_inputs = (;
        run_private_rhf = false,
        private_rhf_electron_count = nothing,
        private_rhf_fixture_role = :route_smoke,
        private_rhf_mixing_kind = :fock_diis,
        private_rhf_max_iterations = 25,
        private_rhf_density_atol = 1.0e-8,
        private_rhf_energy_atol = 1.0e-10,
        private_rhf_residual_atol = 1.0e-8,
        private_rhf_trace_atol = 1.0e-8,
        private_rhf_idempotency_atol = 1.0e-8,
        private_rhf_max_history = nothing,
        private_rhf_diis_start_iteration = 2,
        private_rhf_diis_regularization = 1.0e-12,
        private_rhf_diis_coefficient_max_abs = 25.0,
    )
    route_inputs = (;
        route_family = :pqs_source_box,
        route_kind = :bond_aligned_diatomic_independent_pqs_source_box_core_shell,
        route_shape = (:atom_contact_core, :shared_shell_1, :shared_shell_2),
        product_body_rule = :centered_single_z_slab,
        pqs_retained_rule = :boundary_comx_product_mode_selection,
        product_retained_rule = :product_doside_retained_unit,
        terms = (
            :overlap,
            :position_x,
            :position_y,
            :position_z,
            :x2_x,
            :x2_y,
            :x2_z,
            :kinetic,
        ),
        pair_factor_normalization = :density_normalized,
        white_lindsey_route_shape =
            (:standard_cartesian_units, :low_order_comx_coarsening),
        white_lindsey_mapping_rule = :standard_unit_backbone_mapping_family,
        white_lindsey_nesting_rule = :unit_box_low_order_comx_coarsening,
        white_lindsey_retained_rule = :low_order_unit_comx_retained_basis,
        white_lindsey_operator_rule = :low_order_unit_operator_blocks,
        supplement_policy = :none,
        run_final_basis = false,
        run_h1 = false,
        run_h1_j = false,
        private_rhf_inputs,
    )
    materialization_inputs = (;
        materialize_route = false,
        save_basis_artifact = false,
        save_ham_artifact = false,
        basisfile = nothing,
        hamfile = nothing,
        materializer_backend = nothing,
        materializer_nside = nothing,
        white_lindsey_expansion = nothing,
        white_lindsey_Z = nothing,
        residual_gto_provider_blocks = :none,
    )
    return (; system_inputs, spacing_inputs, parent_inputs, route_inputs, materialization_inputs)
end

function _cartesian_driver_precompile_workload(inputs)
    old_timing_enabled = timing_enabled()
    old_timing_live = timing_live_enabled()
    set_timing!(false)
    set_timing_live!(false)
    try
        system = cartesian_system(inputs.system_inputs)
        recipe = cartesian_recipe(inputs.route_inputs)
        parent =
            cartesian_parent(system, inputs.spacing_inputs, inputs.parent_inputs, recipe)
        shells = cartesian_shells(parent, inputs.spacing_inputs, recipe)
        units = cartesian_units(parent, shells, recipe)
        transforms = cartesian_transforms(units, recipe)
        pairs = cartesian_pair_terms(units, transforms, recipe)
        assembly = cartesian_assembly(parent, shells, units, transforms, pairs, recipe)
        report = cartesian_report(system, parent, assembly, recipe)
        return cartesian_materialization(report, inputs.materialization_inputs)
    finally
        set_timing!(old_timing_enabled)
        set_timing_live!(old_timing_live)
    end
end

if _gaussletbases_precompile_flag("GAUSSLETBASES_PRECOMPILE_CARTESIAN_DRIVER")
    @setup_workload begin
        cartesian_driver_inputs = _cartesian_driver_precompile_inputs()
        @compile_workload begin
            _cartesian_driver_precompile_workload(cartesian_driver_inputs)
        end
    end
end
