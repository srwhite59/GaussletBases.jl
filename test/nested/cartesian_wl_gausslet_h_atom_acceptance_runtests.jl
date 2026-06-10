# Runtime role: decomposed White-Lindsey H acceptance solve.
#
# The former post-CPB acceptance test used one CPB covering the full parent
# window. That is not a decomposed White-Lindsey acceptance path. This file keeps
# the active acceptance path on real q = 5, ns = 5 decomposed retained-unit and
# unit-pair inventories with retained column ranges.

using Test
using LinearAlgebra
using GaussletBases

include("cartesian_white_lindsey_adapter_fixture_helpers.jl")

const WLAcceptanceReadinessCPBM = GaussletBases.CartesianPairBlockMaterialization
const WLAcceptanceReadinessCPB = GaussletBases.CartesianCPB
const WLAcceptanceReadinessCPGB = GaussletBases.CartesianParentGaussletBases
const WLAcceptanceReadinessCBP = GaussletBases.CartesianCPBBlockProviders

function _wl_acceptance_parent_axis_inputs()
    doside_source_1d = _lw_adapter_doside_source_1d()
    parent_axis_bundle_object = (;
        x = doside_source_1d,
        y = doside_source_1d,
        z = doside_source_1d,
    )
    pgdg = doside_source_1d.pgdg_intermediate
    return (;
        parent_axis_counts = (7, 7, 7),
        parent_axis_bundle_object,
        overlap_1d = (; x = pgdg.overlap, y = pgdg.overlap, z = pgdg.overlap),
        kinetic_1d = (; x = pgdg.kinetic, y = pgdg.kinetic, z = pgdg.kinetic),
        spacing_parameter = 0.25,
        distortion_strength = 1.0,
        tail_spacing = 10.0,
    )
end

function _wl_acceptance_center_records()
    return (
        (;
            center_key = :proton_a,
            center_index = 1,
            nuclear_charge = 1.0,
            location = (0.0, 0.0, 0.0),
        ),
    )
end

function _wl_h2plus_acceptance_center_records()
    return (
        (;
            center_key = :proton_a,
            center_index = 1,
            nuclear_charge = 1.0,
            location = (0.0, 0.0, -1.0),
        ),
        (;
            center_key = :proton_b,
            center_index = 2,
            nuclear_charge = 1.0,
            location = (0.0, 0.0, 1.0),
        ),
    )
end

function _wl_h_acceptance_gto_supplement()
    return basis_representation(
        legacy_atomic_gaussian_supplement("H", "cc-pVTZ"; lmax = 0),
    )
end

function _wl_acceptance_direct_core_cpb()
    return WLAcceptanceReadinessCPB.cpb(
        2:6,
        2:6,
        2:6;
        role = :decomposed_wl_gto_readiness_direct_core_cpb,
    )
end

function _wl_acceptance_parent_basis(axis_inputs)
    axis = axis_inputs.parent_axis_bundle_object.x.pgdg_intermediate.basis
    return WLAcceptanceReadinessCPGB.CartesianParentGaussletBasis3D(axis)
end

function _wl_overlap_metric_diagnostics(overlap_matrix, decomposed_inventory)
    sym_s = Symmetric((overlap_matrix + transpose(overlap_matrix)) ./ 2)
    symmetry_error = maximum(abs.(overlap_matrix - transpose(overlap_matrix)))
    diagonal_values = diag(overlap_matrix)
    overlap_eigenvalues = eigvals(sym_s)
    overlap_minimum_eigenvalue = minimum(overlap_eigenvalues)
    overlap_maximum_eigenvalue = maximum(overlap_eigenvalues)
    overlap_condition_estimate =
        overlap_maximum_eigenvalue / overlap_minimum_eigenvalue
    near_zero_eigenvalue_count = count(abs.(overlap_eigenvalues) .<= 1.0e-10)
    negative_eigenvalue_count = count(overlap_eigenvalues .< -1.0e-10)
    rank_estimate = count(overlap_eigenvalues .> 1.0e-10)
    zero_diagonal_count = count(abs.(diagonal_values) .<= 1.0e-10)
    unit_range_start =
        minimum(first(summary.column_range) for summary in decomposed_inventory.unit_summaries)
    unit_range_stop =
        maximum(last(summary.column_range) for summary in decomposed_inventory.unit_summaries)
    missing_prefix_column_count = unit_range_start - 1
    return (;
        matrix_size = size(overlap_matrix),
        symmetry_error,
        diagonal_minimum = minimum(diagonal_values),
        diagonal_maximum = maximum(diagonal_values),
        overlap_minimum_eigenvalue,
        overlap_maximum_eigenvalue,
        overlap_condition_estimate,
        near_zero_eigenvalue_count,
        negative_eigenvalue_count,
        rank_estimate,
        zero_diagonal_count,
        decomposed_unit_column_range_span = unit_range_start:unit_range_stop,
        missing_prefix_column_count,
        boundary_inventory_only =
            missing_prefix_column_count > 0 &&
            zero_diagonal_count >= missing_prefix_column_count,
        blocker =
            missing_prefix_column_count > 0 ?
            :missing_decomposed_wl_interior_retained_operator_inventory :
            overlap_minimum_eigenvalue <= 0.0 ?
            :decomposed_wl_overlap_metric_not_positive_definite :
            nothing,
    )
end

function _wl_lowest_one_electron_energy(
    hamiltonian_matrix,
    overlap_matrix,
    overlap_diagnostics,
)
    sym_h = Symmetric((hamiltonian_matrix + transpose(hamiltonian_matrix)) ./ 2)
    overlap_identity =
        isapprox(overlap_matrix, Matrix{Float64}(I, size(overlap_matrix)...);
            atol = 1.0e-9,
            rtol = 1.0e-9,
        )
    if !isnothing(overlap_diagnostics.blocker)
        return (;
            status = :blocked_decomposed_wl_one_electron_solve,
            blocker = overlap_diagnostics.blocker,
            energy = NaN,
            solve_kind = :blocked_overlap_metric,
            overlap_minimum_eigenvalue =
                overlap_diagnostics.overlap_minimum_eigenvalue,
            overlap_maximum_eigenvalue =
                overlap_diagnostics.overlap_maximum_eigenvalue,
            overlap_condition_estimate =
                overlap_diagnostics.overlap_condition_estimate,
            overlap_identity,
        )
    end
    sym_s = Symmetric((overlap_matrix + transpose(overlap_matrix)) ./ 2)
    values = overlap_identity ? eigvals(sym_h) : eigen(sym_h, sym_s).values
    return (;
        status = :materialized_decomposed_wl_one_electron_solve,
        blocker = nothing,
        energy = minimum(values),
        solve_kind =
            overlap_identity ? :ordinary_symmetric : :generalized_symmetric,
        overlap_minimum_eigenvalue =
            overlap_diagnostics.overlap_minimum_eigenvalue,
        overlap_maximum_eigenvalue =
            overlap_diagnostics.overlap_maximum_eigenvalue,
        overlap_condition_estimate =
            overlap_diagnostics.overlap_condition_estimate,
        overlap_identity,
    )
end

function _wl_decomposed_h_atom_acceptance_report()
    adapter = WLAcceptanceReadinessCPBM.white_lindsey_boundary_stratum_one_body_adapter_summary()
    local_terms = adapter.supported_one_body_terms
    route_global_terms = WLAcceptanceReadinessCPBM.route_global_safe_one_body_terms()
    route_global_supported_terms =
        WLAcceptanceReadinessCPBM._route_global_one_body_supported_terms()
    seed_report = GaussletBases._white_lindsey_low_order_materialized_seed_report()
    axis_inputs = _wl_acceptance_parent_axis_inputs()
    decomposed_inventory =
        WLAcceptanceReadinessCPBM.white_lindsey_decomposed_unit_pair_inventory(
            seed_report,
        )
    overlap_global = nothing
    kinetic_global = nothing
    by_center_global = nothing
    hamiltonian = nothing
    overlap_diagnostics = nothing
    solve = nothing
    elapsed_seconds = @elapsed begin
        overlap_global =
            WLAcceptanceReadinessCPBM.route_global_decomposed_wl_overlap_matrix(
                seed_report;
                parent_axis_counts = axis_inputs.parent_axis_counts,
                parent_axis_bundle_object =
                    axis_inputs.parent_axis_bundle_object,
                overlap_1d = axis_inputs.overlap_1d,
            )
        kinetic_global =
            WLAcceptanceReadinessCPBM.route_global_decomposed_wl_kinetic_matrix(
                seed_report;
                parent_axis_counts = axis_inputs.parent_axis_counts,
                parent_axis_bundle_object =
                    axis_inputs.parent_axis_bundle_object,
                overlap_1d = axis_inputs.overlap_1d,
                kinetic_1d = axis_inputs.kinetic_1d,
            )
        by_center_global =
        WLAcceptanceReadinessCPBM.route_global_electron_nuclear_by_center_matrices(
            seed_report;
                parent_axis_counts = axis_inputs.parent_axis_counts,
            parent_axis_bundle_object =
                    axis_inputs.parent_axis_bundle_object,
            coulomb_expansion = coulomb_gaussian_expansion(doacc = false),
            center_records = _wl_acceptance_center_records(),
        )
        hamiltonian =
            WLAcceptanceReadinessCPBM.white_lindsey_decomposed_one_electron_hamiltonian(
                kinetic_global,
                by_center_global,
            )
        overlap_diagnostics =
            _wl_overlap_metric_diagnostics(overlap_global.matrix, decomposed_inventory)
        solve = _wl_lowest_one_electron_energy(
            hamiltonian.matrix,
            overlap_global.matrix,
            overlap_diagnostics,
        )
    end
    h_reference = -0.5
    old_full_window_reference = -0.4832079279118124

    return (;
        object_kind = :decomposed_wl_h_atom_acceptance_report,
        acceptance_suite =
            :decomposed_wl_gausslet_one_electron_acceptance,
        acceptance_fixture = :h_atom,
        status =
            solve.status === :materialized_decomposed_wl_one_electron_solve ?
            :materialized_decomposed_wl_h_atom_acceptance :
            :blocked_decomposed_wl_h_atom_acceptance,
        blocker = solve.blocker,
        q = 5,
        ns = 5,
        n_s = 5,
        core_spacing = axis_inputs.spacing_parameter,
        distortion_strength = axis_inputs.distortion_strength,
        tail_spacing = axis_inputs.tail_spacing,
        parent_axis_counts = axis_inputs.parent_axis_counts,
        decomposed_wl_units_consumed = true,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        supported_decomposed_one_body_terms = local_terms,
        route_global_safe_one_body_terms = route_global_terms,
        route_global_supported_one_body_terms = route_global_supported_terms,
        decomposed_overlap_available = :overlap in local_terms,
        decomposed_kinetic_available = :kinetic in local_terms,
        decomposed_electron_nuclear_by_center_available =
            :electron_nuclear_by_center in local_terms,
        decomposed_electron_nuclear_by_center_selector_available =
            :electron_nuclear_by_center in local_terms,
        route_global_overlap_adapter_available = :overlap in route_global_terms,
        route_global_kinetic_adapter_available = :kinetic in route_global_terms,
        route_global_electron_nuclear_by_center_adapter_available =
            :electron_nuclear_by_center in route_global_supported_terms,
        route_global_overlap_status = overlap_global.status,
        route_global_kinetic_status = kinetic_global.status,
        route_global_overlap_matrix_materialized =
            overlap_global.global_overlap_matrix_materialized,
        route_global_kinetic_matrix_materialized =
            kinetic_global.global_kinetic_matrix_materialized,
        unit_inventory_audit_source =
            decomposed_inventory.source_kind,
        terminal_shellification_unit_inventory_exposed = true,
        terminal_shellification_unit_inventory_granularity =
            :terminal_region_units,
        terminal_shellification_pair_inventory_exposed = false,
        terminal_shellification_pair_inventory_status =
            :deferred_terminal_shellification_pair_inventory,
        terminal_shellification_pair_materialization_status =
            :deferred_terminal_shellification_pair_materialization,
        decomposed_unit_pair_inventory_status = decomposed_inventory.status,
        decomposed_unit_pair_inventory_blocker = decomposed_inventory.blocker,
        decomposed_unit_pair_inventory_source_kind =
            decomposed_inventory.source_kind,
        route_owned_decomposed_unit_pair_inventory_available =
            decomposed_inventory.decomposed_wl_unit_pair_inventory_available,
        decomposed_unit_count = decomposed_inventory.unit_count,
        decomposed_unit_pair_count = decomposed_inventory.pair_count,
        decomposed_unit_key_sample =
            Tuple(Iterators.take(decomposed_inventory.unit_keys, 6)),
        decomposed_unit_pair_key_sample =
            Tuple(Iterators.take(decomposed_inventory.pair_keys, 6)),
        retained_unit_column_ranges_materialized =
            decomposed_inventory.retained_unit_column_ranges_materialized,
        retained_dimension_from_decomposed_unit_inventory_available =
            !isnothing(decomposed_inventory.retained_dimension),
        decomposed_unit_pair_column_ranges_available =
            decomposed_inventory.decomposed_unit_pair_column_ranges_available,
        decomposed_unit_pair_inventory_retained_dimension =
            decomposed_inventory.retained_dimension,
        retained_global_dimension_source =
            decomposed_inventory.retained_dimension_status,
        route_global_by_center_acceptance_matrix_available =
            by_center_global.route_global_by_center_matrices_materialized,
        route_global_by_center_acceptance_matrix_status =
            by_center_global.status,
        route_global_by_center_acceptance_matrix_blocker =
            by_center_global.blocker,
        route_global_by_center_matrix_count =
            by_center_global.by_center_matrix_count,
        route_global_by_center_center_count =
            by_center_global.center_count,
        route_global_by_center_nuclear_charge_applied =
            by_center_global.nuclear_charge_applied,
        route_global_by_center_centers_summed =
            by_center_global.centers_summed,
        hamiltonian_status = hamiltonian.status,
        hamiltonian_matrix_materialized =
            hamiltonian.hamiltonian_data_materialized,
        nuclear_charge_application_stage =
            hamiltonian.charge_application_stage,
        nuclear_charge_applied_at_hamiltonian_assembly =
            hamiltonian.nuclear_charge_applied_at_hamiltonian_assembly,
        center_summation_stage = hamiltonian.center_summation_stage,
        centers_summed_at_hamiltonian_assembly =
            hamiltonian.centers_summed_at_hamiltonian_assembly,
        retained_dimension = hamiltonian.retained_dimension,
        overlap_matrix_dimension = size(overlap_global.matrix, 1),
        hamiltonian_matrix_dimension = size(hamiltonian.matrix, 1),
        overlap_minimum_eigenvalue = solve.overlap_minimum_eigenvalue,
        overlap_maximum_eigenvalue = solve.overlap_maximum_eigenvalue,
        overlap_condition_estimate = solve.overlap_condition_estimate,
        overlap_symmetry_error = overlap_diagnostics.symmetry_error,
        overlap_diagonal_minimum = overlap_diagnostics.diagonal_minimum,
        overlap_diagonal_maximum = overlap_diagnostics.diagonal_maximum,
        overlap_near_zero_eigenvalue_count =
            overlap_diagnostics.near_zero_eigenvalue_count,
        overlap_negative_eigenvalue_count =
            overlap_diagnostics.negative_eigenvalue_count,
        overlap_rank_estimate = overlap_diagnostics.rank_estimate,
        overlap_zero_diagonal_count = overlap_diagnostics.zero_diagonal_count,
        decomposed_unit_column_range_span =
            overlap_diagnostics.decomposed_unit_column_range_span,
        missing_interior_retained_column_count =
            overlap_diagnostics.missing_prefix_column_count,
        boundary_inventory_only = overlap_diagnostics.boundary_inventory_only,
        solve_kind = solve.solve_kind,
        solve_status = solve.status,
        solve_blocker = solve.blocker,
        lowest_h_atom_energy = solve.energy,
        h_atom_exact_energy = h_reference,
        h_atom_distance_from_exact = solve.energy - h_reference,
        h_atom_distance_from_full_window_transition_baseline =
            solve.energy - old_full_window_reference,
        elapsed_seconds,
        fixed_block_operator_matrices_available =
            :overlap in route_global_terms && :kinetic in route_global_terms,
        fixed_block_operator_matrices_used = false,
        fixed_block_operator_matrix_source_rejected =
            :nested_fixed_block_is_not_decomposed_wl_acceptance_path,
        acceptance_energy_materialized =
            solve.status === :materialized_decomposed_wl_one_electron_solve,
        h_atom_acceptance_active =
            solve.status === :materialized_decomposed_wl_one_electron_solve,
        acceptance_fixture_active =
            solve.status === :materialized_decomposed_wl_one_electron_solve,
    )
end

function _wl_decomposed_h2plus_acceptance_report()
    adapter = WLAcceptanceReadinessCPBM.white_lindsey_boundary_stratum_one_body_adapter_summary()
    local_terms = adapter.supported_one_body_terms
    route_global_terms = WLAcceptanceReadinessCPBM.route_global_safe_one_body_terms()
    route_global_supported_terms =
        WLAcceptanceReadinessCPBM._route_global_one_body_supported_terms()
    seed_report = GaussletBases._white_lindsey_low_order_materialized_seed_report()
    axis_inputs = _wl_acceptance_parent_axis_inputs()
    center_records = _wl_h2plus_acceptance_center_records()
    decomposed_inventory =
        WLAcceptanceReadinessCPBM.white_lindsey_decomposed_unit_pair_inventory(
            seed_report,
        )
    overlap_global = nothing
    kinetic_global = nothing
    by_center_global = nothing
    hamiltonian = nothing
    overlap_diagnostics = nothing
    solve = nothing
    elapsed_seconds = @elapsed begin
        overlap_global =
            WLAcceptanceReadinessCPBM.route_global_decomposed_wl_overlap_matrix(
                seed_report;
                parent_axis_counts = axis_inputs.parent_axis_counts,
                parent_axis_bundle_object =
                    axis_inputs.parent_axis_bundle_object,
                overlap_1d = axis_inputs.overlap_1d,
            )
        kinetic_global =
            WLAcceptanceReadinessCPBM.route_global_decomposed_wl_kinetic_matrix(
                seed_report;
                parent_axis_counts = axis_inputs.parent_axis_counts,
                parent_axis_bundle_object =
                    axis_inputs.parent_axis_bundle_object,
                overlap_1d = axis_inputs.overlap_1d,
                kinetic_1d = axis_inputs.kinetic_1d,
            )
        by_center_global =
        WLAcceptanceReadinessCPBM.route_global_electron_nuclear_by_center_matrices(
            seed_report;
                parent_axis_counts = axis_inputs.parent_axis_counts,
            parent_axis_bundle_object =
                    axis_inputs.parent_axis_bundle_object,
            coulomb_expansion = coulomb_gaussian_expansion(doacc = false),
            center_records,
        )
        hamiltonian =
            WLAcceptanceReadinessCPBM.white_lindsey_decomposed_one_electron_hamiltonian(
                kinetic_global,
                by_center_global,
            )
        overlap_diagnostics =
            _wl_overlap_metric_diagnostics(overlap_global.matrix, decomposed_inventory)
        solve = _wl_lowest_one_electron_energy(
            hamiltonian.matrix,
            overlap_global.matrix,
            overlap_diagnostics,
        )
    end
    bond_length = 2.0
    nuclear_repulsion = 1.0 / bond_length
    exact_electronic = -1.1026342144949465
    exact_total = -0.6026342144949465
    total_energy = solve.energy + nuclear_repulsion
    old_full_window_total_reference = -0.5971828374927926

    return (;
        object_kind = :decomposed_wl_h2plus_acceptance_report,
        acceptance_suite =
            :decomposed_wl_gausslet_one_electron_acceptance,
        acceptance_fixture = :h2plus,
        status =
            solve.status === :materialized_decomposed_wl_one_electron_solve ?
            :materialized_decomposed_wl_h2plus_acceptance :
            :blocked_decomposed_wl_h2plus_acceptance,
        blocker = solve.blocker,
        q = 5,
        ns = 5,
        n_s = 5,
        bond_length,
        center_count = length(center_records),
        center_keys = Tuple(record.center_key for record in center_records),
        center_locations = Tuple(record.location for record in center_records),
        nuclear_charges = Tuple(record.nuclear_charge for record in center_records),
        nuclear_repulsion,
        core_spacing = axis_inputs.spacing_parameter,
        distortion_strength = axis_inputs.distortion_strength,
        tail_spacing = axis_inputs.tail_spacing,
        parent_axis_counts = axis_inputs.parent_axis_counts,
        decomposed_wl_units_consumed = true,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        supported_decomposed_one_body_terms = local_terms,
        route_global_safe_one_body_terms = route_global_terms,
        route_global_supported_one_body_terms = route_global_supported_terms,
        decomposed_overlap_available = :overlap in local_terms,
        decomposed_kinetic_available = :kinetic in local_terms,
        decomposed_electron_nuclear_by_center_available =
            :electron_nuclear_by_center in local_terms,
        route_global_overlap_adapter_available = :overlap in route_global_terms,
        route_global_kinetic_adapter_available = :kinetic in route_global_terms,
        route_global_electron_nuclear_by_center_adapter_available =
            :electron_nuclear_by_center in route_global_supported_terms,
        route_global_overlap_status = overlap_global.status,
        route_global_kinetic_status = kinetic_global.status,
        route_global_overlap_matrix_materialized =
            overlap_global.global_overlap_matrix_materialized,
        route_global_kinetic_matrix_materialized =
            kinetic_global.global_kinetic_matrix_materialized,
        decomposed_unit_pair_inventory_status = decomposed_inventory.status,
        decomposed_unit_pair_inventory_blocker = decomposed_inventory.blocker,
        decomposed_unit_pair_inventory_source_kind =
            decomposed_inventory.source_kind,
        route_owned_decomposed_unit_pair_inventory_available =
            decomposed_inventory.decomposed_wl_unit_pair_inventory_available,
        decomposed_unit_count = decomposed_inventory.unit_count,
        decomposed_unit_pair_count = decomposed_inventory.pair_count,
        decomposed_unit_key_sample =
            Tuple(Iterators.take(decomposed_inventory.unit_keys, 6)),
        retained_unit_column_ranges_materialized =
            decomposed_inventory.retained_unit_column_ranges_materialized,
        decomposed_unit_pair_column_ranges_available =
            decomposed_inventory.decomposed_unit_pair_column_ranges_available,
        decomposed_unit_pair_inventory_retained_dimension =
            decomposed_inventory.retained_dimension,
        route_global_by_center_acceptance_matrix_available =
            by_center_global.route_global_by_center_matrices_materialized,
        route_global_by_center_acceptance_matrix_status =
            by_center_global.status,
        route_global_by_center_acceptance_matrix_blocker =
            by_center_global.blocker,
        route_global_by_center_matrix_count =
            by_center_global.by_center_matrix_count,
        route_global_by_center_center_count =
            by_center_global.center_count,
        route_global_by_center_center_keys = by_center_global.center_keys,
        route_global_by_center_nuclear_charge_applied =
            by_center_global.nuclear_charge_applied,
        route_global_by_center_centers_summed =
            by_center_global.centers_summed,
        hamiltonian_status = hamiltonian.status,
        hamiltonian_matrix_materialized =
            hamiltonian.hamiltonian_data_materialized,
        nuclear_charge_application_stage =
            hamiltonian.charge_application_stage,
        nuclear_charge_applied_at_hamiltonian_assembly =
            hamiltonian.nuclear_charge_applied_at_hamiltonian_assembly,
        center_summation_stage = hamiltonian.center_summation_stage,
        centers_summed_at_hamiltonian_assembly =
            hamiltonian.centers_summed_at_hamiltonian_assembly,
        retained_dimension = hamiltonian.retained_dimension,
        overlap_matrix_dimension = size(overlap_global.matrix, 1),
        hamiltonian_matrix_dimension = size(hamiltonian.matrix, 1),
        overlap_minimum_eigenvalue = solve.overlap_minimum_eigenvalue,
        overlap_maximum_eigenvalue = solve.overlap_maximum_eigenvalue,
        overlap_condition_estimate = solve.overlap_condition_estimate,
        overlap_symmetry_error = overlap_diagnostics.symmetry_error,
        overlap_near_zero_eigenvalue_count =
            overlap_diagnostics.near_zero_eigenvalue_count,
        overlap_negative_eigenvalue_count =
            overlap_diagnostics.negative_eigenvalue_count,
        overlap_rank_estimate = overlap_diagnostics.rank_estimate,
        decomposed_unit_column_range_span =
            overlap_diagnostics.decomposed_unit_column_range_span,
        missing_interior_retained_column_count =
            overlap_diagnostics.missing_prefix_column_count,
        boundary_inventory_only = overlap_diagnostics.boundary_inventory_only,
        solve_kind = solve.solve_kind,
        solve_status = solve.status,
        solve_blocker = solve.blocker,
        electronic_energy = solve.energy,
        exact_electronic_energy = exact_electronic,
        electronic_energy_distance_from_exact = solve.energy - exact_electronic,
        total_bo_energy = total_energy,
        exact_total_bo_energy = exact_total,
        total_bo_energy_distance_from_exact = total_energy - exact_total,
        total_bo_energy_distance_from_full_window_transition_baseline =
            total_energy - old_full_window_total_reference,
        elapsed_seconds,
        fixed_block_operator_matrices_used = false,
        fixed_block_operator_matrix_source_rejected =
            :nested_fixed_block_is_not_decomposed_wl_acceptance_path,
        acceptance_energy_materialized =
            solve.status === :materialized_decomposed_wl_one_electron_solve,
        h2plus_acceptance_active =
            solve.status === :materialized_decomposed_wl_one_electron_solve,
        acceptance_fixture_active =
            solve.status === :materialized_decomposed_wl_one_electron_solve,
    )
end

function _wl_decomposed_h_gto_supplement_readiness_report()
    seed_report = GaussletBases._white_lindsey_low_order_materialized_seed_report()
    axis_inputs = _wl_acceptance_parent_axis_inputs()
    decomposed_inventory =
        WLAcceptanceReadinessCPBM.white_lindsey_decomposed_unit_pair_inventory(
            seed_report,
        )
    route_global_terms = WLAcceptanceReadinessCPBM.route_global_safe_one_body_terms()
    route_global_supported_terms =
        WLAcceptanceReadinessCPBM._route_global_one_body_supported_terms()
    supplement = _wl_h_acceptance_gto_supplement()
    parent = _wl_acceptance_parent_basis(axis_inputs)
    direct_core_cpb = _wl_acceptance_direct_core_cpb()
    expansion = coulomb_gaussian_expansion(doacc = false)
    center_records = _wl_acceptance_center_records()
    pgdg = axis_inputs.parent_axis_bundle_object.x.pgdg_intermediate
    required_local_terms = (
        :mixed_gto_overlap,
        :mixed_gto_kinetic,
        :mixed_gto_electron_nuclear_by_center,
        :gto_overlap,
        :gto_kinetic,
        :gto_electron_nuclear_by_center,
    )
    bundle = WLAcceptanceReadinessCBP.cpb_gto_supplement_local_operator_bundle(
        parent,
        direct_core_cpb,
        supplement;
        expansion,
        center_records,
    )
    bundle_summary = WLAcceptanceReadinessCBP.summary(bundle)
    mixed_overlap_summary =
        WLAcceptanceReadinessCBP.summary(bundle.mixed_blocks.overlap)
    mixed_kinetic_summary =
        WLAcceptanceReadinessCBP.summary(bundle.mixed_blocks.kinetic)
    mixed_nuclear_summaries = Tuple(
        WLAcceptanceReadinessCBP.summary(block) for
        block in bundle.mixed_nuclear_by_center_blocks
    )
    gto_overlap_summary =
        WLAcceptanceReadinessCBP.summary(bundle.gto_blocks.overlap)
    gto_kinetic_summary =
        WLAcceptanceReadinessCBP.summary(bundle.gto_blocks.kinetic)
    gto_nuclear_summaries = Tuple(
        WLAcceptanceReadinessCBP.summary(block) for
        block in bundle.gto_nuclear_by_center_blocks
    )
    gto_shapes = (;
        overlap = gto_overlap_summary.dense_block_shape,
        kinetic = gto_kinetic_summary.dense_block_shape,
        nuclear_by_center =
            Tuple(summary.dense_block_shape for summary in gto_nuclear_summaries),
    )
    mixed_shapes = (;
        overlap = mixed_overlap_summary.dense_block_shape,
        kinetic = mixed_kinetic_summary.dense_block_shape,
        nuclear_by_center =
            Tuple(summary.dense_block_shape for summary in mixed_nuclear_summaries),
    )
    mixed_gto_blocker =
        nothing
    mixed_gto_blocker_source =
        :_cpb_mixed_gto_parent_axis_representations
    mixed_gto_rejected_helper = :unavailable
    mixed_gto_blocker_detail =
        :mapped_working_gaussian_proxy_axis_representation_available
    mixed_gto_required_source =
        :mapped_ordinary_working_gaussian_proxy_axis_representation
    combined_layout =
        WLAcceptanceReadinessCPBM.route_global_combined_gto_basis_layout(
            decomposed_inventory,
            supplement,
        )
    route_global_combined_basis_layout_status = combined_layout.status
    route_global_combined_basis_layout_blocker = combined_layout.blocker
    readiness_blocker =
        route_global_combined_basis_layout_status ===
        :available_route_global_combined_gto_basis_layout ?
        :missing_route_global_combined_gto_matrix_assembly :
        route_global_combined_basis_layout_blocker
    decomposed_unit_range_start = minimum(
        first(summary.column_range) for summary in decomposed_inventory.unit_summaries
    )
    decomposed_unit_range_stop = maximum(
        last(summary.column_range) for summary in decomposed_inventory.unit_summaries
    )
    decomposed_unit_column_range_span =
        decomposed_unit_range_start:decomposed_unit_range_stop
    return (;
        object_kind = :decomposed_wl_h_gto_supplement_acceptance_readiness_report,
        status = :blocked_decomposed_wl_gto_supplement_acceptance_readiness,
        blocker = readiness_blocker,
        acceptance_suite =
            :decomposed_wl_gausslet_plus_gto_one_electron_acceptance,
        acceptance_fixture = :h_atom_gto_supplement_readiness,
        base_decomposed_wl_fixture = :h_atom,
        q = 5,
        ns = 5,
        n_s = 5,
        parent_axis_counts = axis_inputs.parent_axis_counts,
        decomposed_wl_gausslet_base_path_available =
            decomposed_inventory.status ===
            :available_white_lindsey_decomposed_unit_pair_inventory &&
            :overlap in route_global_terms &&
            :kinetic in route_global_terms &&
            :electron_nuclear_by_center in route_global_supported_terms,
        decomposed_unit_pair_inventory_status = decomposed_inventory.status,
        decomposed_unit_pair_inventory_blocker = decomposed_inventory.blocker,
        decomposed_unit_count = decomposed_inventory.unit_count,
        decomposed_unit_pair_count = decomposed_inventory.pair_count,
        decomposed_retained_dimension = decomposed_inventory.retained_dimension,
        decomposed_unit_column_range_span,
        route_global_overlap_available = :overlap in route_global_terms,
        route_global_kinetic_available = :kinetic in route_global_terms,
        route_global_electron_nuclear_by_center_available =
            :electron_nuclear_by_center in route_global_supported_terms,
        parent_pgdg_intermediate_available = !isnothing(pgdg),
        parent_pgdg_overlap_available = hasproperty(pgdg, :overlap),
        parent_pgdg_kinetic_available = hasproperty(pgdg, :kinetic),
        parent_pgdg_position_available = hasproperty(pgdg, :position),
        parent_pgdg_x2_available = hasproperty(pgdg, :x2),
        parent_pgdg_gaussian_factor_terms_available =
            hasproperty(pgdg, :gaussian_factor_terms),
        parent_pgdg_pair_factor_terms_available =
            hasproperty(pgdg, :pair_factor_terms),
        parent_pgdg_weights_available = hasproperty(pgdg, :weights),
        mixed_gto_axis_contract_status =
            :available_mixed_gto_axis_representation_adapter,
        mixed_gto_axis_contract_blocker = mixed_gto_blocker,
        mixed_gto_rejected_helper,
        mixed_gto_axis_representation_source = mixed_gto_required_source,
        mixed_gto_required_source,
        supplement_source_available = true,
        supplement_kind = supplement.supplement_kind,
        supplement_atom = supplement.metadata.atom,
        supplement_basis_name = supplement.metadata.basis_name,
        supplement_lmax = supplement.metadata.lmax,
        supplement_orbital_count = length(supplement.orbitals),
        supplement_orbital_labels =
            Tuple(orbital.label for orbital in supplement.orbitals),
        cpb_local_gto_bundle_status = bundle_summary.status,
        cpb_local_gto_bundle_blocker = bundle_summary.blocker,
        cpb_local_gto_bundle_terms_available =
            bundle_summary.provider_level_local_blocks_materialized,
        cpb_local_gto_required_terms = required_local_terms,
        cpb_local_gto_included_terms = bundle_summary.included_terms,
        cpb_local_gto_missing_terms = bundle_summary.missing_terms,
        cpb_local_direct_core_support_count =
            WLAcceptanceReadinessCPB.support_count(direct_core_cpb),
        mixed_gausslet_gto_overlap_available =
            mixed_overlap_summary.dense_block_available,
        mixed_gausslet_gto_kinetic_available =
            mixed_kinetic_summary.dense_block_available,
        mixed_gausslet_gto_nuclear_by_center_available =
            all(summary.dense_block_available for summary in mixed_nuclear_summaries),
        mixed_gausslet_gto_blocker = mixed_gto_blocker,
        mixed_gausslet_gto_blocker_source = mixed_gto_blocker_source,
        mixed_gausslet_gto_rejected_helper = mixed_gto_rejected_helper,
        mixed_gausslet_gto_blocker_detail = mixed_gto_blocker_detail,
        mixed_gausslet_gto_required_source = mixed_gto_required_source,
        gto_gto_overlap_available =
            gto_overlap_summary.dense_block_available,
        gto_gto_kinetic_available =
            gto_kinetic_summary.dense_block_available,
        gto_gto_nuclear_by_center_available = all(
            summary -> summary.dense_block_available,
            gto_nuclear_summaries,
        ),
        mixed_gausslet_gto_shapes = mixed_shapes,
        gto_gto_shapes = gto_shapes,
        route_global_combined_basis_layout_status,
        route_global_combined_basis_layout_blocker,
        route_global_combined_basis_layout_kind = combined_layout.layout_kind,
        gausslet_retained_range = combined_layout.gausslet_retained_range,
        gausslet_retained_dimension =
            combined_layout.gausslet_retained_dimension,
        gto_supplement_range = combined_layout.gto_supplement_range,
        gto_supplement_orbital_count =
            combined_layout.gto_supplement_orbital_count,
        total_combined_dimension = combined_layout.total_combined_dimension,
        combined_basis_dimension = combined_layout.combined_basis_dimension,
        combined_block_layout_keys = combined_layout.block_layout_keys,
        gausslet_gausslet_block_layout =
            combined_layout.gausslet_gausslet_block_layout,
        gausslet_gto_block_layout = combined_layout.gausslet_gto_block_layout,
        gto_gausslet_block_layout = combined_layout.gto_gausslet_block_layout,
        gto_gto_block_layout = combined_layout.gto_gto_block_layout,
        mixed_cpb_gto_blocks_orientation =
            combined_layout.mixed_cpb_gto_blocks_orientation,
        gto_gto_blocks_orientation = combined_layout.gto_gto_blocks_orientation,
        final_combined_overlap_layout_available =
            combined_layout.combined_overlap_layout_available,
        final_combined_hamiltonian_layout_available =
            combined_layout.combined_hamiltonian_layout_available,
        route_global_combined_overlap_matrix_materialized =
            combined_layout.combined_overlap_matrix_materialized,
        route_global_combined_hamiltonian_matrix_materialized =
            combined_layout.combined_hamiltonian_matrix_materialized,
        gto_route_global_blocks_materialized = false,
        gto_hamiltonian_assembly_materialized = false,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        exports_or_artifacts = false,
    )
end

@testset "decomposed WL gausslet-only H atom acceptance" begin
    report = _wl_decomposed_h_atom_acceptance_report()
    println("decomposed WL gausslet-only H atom acceptance: ", report)

    @test report.status == :materialized_decomposed_wl_h_atom_acceptance
    @test isnothing(report.blocker)
    @test report.acceptance_suite ==
          :decomposed_wl_gausslet_one_electron_acceptance
    @test report.acceptance_fixture == :h_atom
    @test report.q == 5
    @test report.ns == 5
    @test report.n_s == 5
    @test report.decomposed_overlap_available
    @test report.decomposed_kinetic_available
    @test report.decomposed_electron_nuclear_by_center_available
    @test report.decomposed_electron_nuclear_by_center_selector_available
    @test report.route_global_overlap_adapter_available
    @test report.route_global_kinetic_adapter_available
    @test report.route_global_electron_nuclear_by_center_adapter_available
    @test report.route_global_overlap_status ==
          :materialized_route_global_overlap_matrix
    @test report.route_global_kinetic_status ==
          :materialized_route_global_kinetic_matrix
    @test report.route_global_overlap_matrix_materialized
    @test report.route_global_kinetic_matrix_materialized
    @test report.terminal_shellification_unit_inventory_exposed
    @test report.terminal_shellification_unit_inventory_granularity ==
          :terminal_region_units
    @test !report.terminal_shellification_pair_inventory_exposed
    @test report.terminal_shellification_pair_inventory_status ==
          :deferred_terminal_shellification_pair_inventory
    @test report.decomposed_unit_pair_inventory_status ==
          :available_white_lindsey_decomposed_unit_pair_inventory
    @test isnothing(report.decomposed_unit_pair_inventory_blocker)
    @test report.decomposed_unit_pair_inventory_source_kind ==
          :white_lindsey_low_order_materialized_seed_ranges
    @test report.route_owned_decomposed_unit_pair_inventory_available
    @test report.decomposed_unit_count == 27
    @test report.decomposed_unit_pair_count == 378
    @test first(report.decomposed_unit_key_sample) == :white_lindsey_seed_direct_core
    @test report.retained_unit_column_ranges_materialized
    @test report.retained_dimension_from_decomposed_unit_inventory_available
    @test report.decomposed_unit_pair_column_ranges_available
    @test report.decomposed_unit_pair_inventory_retained_dimension == 223
    @test report.retained_global_dimension_source ==
          :available_from_decomposed_wl_unit_column_ranges
    @test report.route_global_by_center_acceptance_matrix_available
    @test report.route_global_by_center_acceptance_matrix_status ==
          :materialized_route_global_electron_nuclear_by_center_matrix_set
    @test isnothing(report.route_global_by_center_acceptance_matrix_blocker)
    @test report.route_global_by_center_matrix_count == 1
    @test report.route_global_by_center_center_count == 1
    @test !report.route_global_by_center_nuclear_charge_applied
    @test !report.route_global_by_center_centers_summed
    @test report.hamiltonian_status ==
          :materialized_white_lindsey_decomposed_one_electron_hamiltonian
    @test report.hamiltonian_matrix_materialized
    @test report.nuclear_charge_application_stage ==
          :white_lindsey_hamiltonian_assembly
    @test report.nuclear_charge_applied_at_hamiltonian_assembly
    @test report.center_summation_stage == :white_lindsey_hamiltonian_assembly
    @test report.centers_summed_at_hamiltonian_assembly
    @test report.retained_dimension == 223
    @test report.overlap_matrix_dimension == report.retained_dimension
    @test report.hamiltonian_matrix_dimension == report.retained_dimension
    @test report.overlap_minimum_eigenvalue > 0.0
    @test report.overlap_condition_estimate < 1.0000001
    @test report.overlap_symmetry_error < 1.0e-12
    @test isapprox(report.overlap_diagonal_minimum, 1.0; atol = 1.0e-12)
    @test isapprox(report.overlap_diagonal_maximum, 1.0; atol = 1.0e-12)
    @test report.overlap_near_zero_eigenvalue_count == 0
    @test report.overlap_negative_eigenvalue_count == 0
    @test report.overlap_rank_estimate == report.retained_dimension
    @test report.overlap_zero_diagonal_count == 0
    @test report.decomposed_unit_column_range_span == 1:223
    @test report.missing_interior_retained_column_count == 0
    @test !report.boundary_inventory_only
    @test report.solve_status == :materialized_decomposed_wl_one_electron_solve
    @test isnothing(report.solve_blocker)
    @test report.solve_kind == :ordinary_symmetric
    @test report.lowest_h_atom_energy > report.h_atom_exact_energy
    @test report.lowest_h_atom_energy < -0.45
    @test isapprox(
        report.lowest_h_atom_energy,
        -0.4788666674548281;
        atol = 1.0e-10,
    )
    @test report.fixed_block_operator_matrices_available
    @test !report.fixed_block_operator_matrices_used
    @test report.decomposed_wl_units_consumed
    @test !report.full_parent_window_cpb_used
    @test !report.direct_cartesian_product_assembly_used
    @test !report.ordinary_cartesian_ida_operators_used
    @test report.acceptance_energy_materialized
    @test report.h_atom_acceptance_active
    @test report.acceptance_fixture_active
    @test report.elapsed_seconds >= 0.0
end

@testset "decomposed WL H GTO supplement acceptance readiness" begin
    report = _wl_decomposed_h_gto_supplement_readiness_report()
    println("decomposed WL H GTO supplement readiness: ", report)

    @test report.status ==
          :blocked_decomposed_wl_gto_supplement_acceptance_readiness
    @test report.blocker ==
          :missing_route_global_combined_gto_matrix_assembly
    @test report.acceptance_suite ==
          :decomposed_wl_gausslet_plus_gto_one_electron_acceptance
    @test report.acceptance_fixture == :h_atom_gto_supplement_readiness
    @test report.base_decomposed_wl_fixture == :h_atom
    @test report.q == 5
    @test report.ns == 5
    @test report.n_s == 5
    @test report.parent_axis_counts == (7, 7, 7)
    @test report.decomposed_wl_gausslet_base_path_available
    @test report.decomposed_unit_pair_inventory_status ==
          :available_white_lindsey_decomposed_unit_pair_inventory
    @test isnothing(report.decomposed_unit_pair_inventory_blocker)
    @test report.decomposed_unit_count == 27
    @test report.decomposed_unit_pair_count == 378
    @test report.decomposed_retained_dimension == 223
    @test report.decomposed_unit_column_range_span == 1:223
    @test report.route_global_overlap_available
    @test report.route_global_kinetic_available
    @test report.route_global_electron_nuclear_by_center_available
    @test report.parent_pgdg_intermediate_available
    @test report.parent_pgdg_overlap_available
    @test report.parent_pgdg_kinetic_available
    @test report.parent_pgdg_position_available
    @test report.parent_pgdg_x2_available
    @test report.parent_pgdg_gaussian_factor_terms_available
    @test report.parent_pgdg_pair_factor_terms_available
    @test report.parent_pgdg_weights_available
    @test report.mixed_gto_axis_contract_status ==
          :available_mixed_gto_axis_representation_adapter
    @test isnothing(report.mixed_gto_axis_contract_blocker)
    @test report.mixed_gto_rejected_helper == :unavailable
    @test report.mixed_gto_axis_representation_source ==
          :mapped_ordinary_working_gaussian_proxy_axis_representation
    @test report.mixed_gto_required_source ==
          :mapped_ordinary_working_gaussian_proxy_axis_representation
    @test report.supplement_source_available
    @test report.supplement_kind == :atomic_cartesian_shell
    @test report.supplement_atom == "H"
    @test report.supplement_basis_name == "cc-pVTZ"
    @test report.supplement_lmax == 0
    @test report.supplement_orbital_count == 3
    @test report.supplement_orbital_labels == ("s1", "s2", "s3")
    @test report.cpb_local_gto_bundle_status ==
          :materialized_cpb_gto_supplement_local_operator_bundle
    @test isnothing(report.cpb_local_gto_bundle_blocker)
    @test report.cpb_local_gto_bundle_terms_available
    @test report.cpb_local_gto_missing_terms == ()
    for term in (
        :mixed_gto_overlap,
        :mixed_gto_kinetic,
        :mixed_gto_electron_nuclear_by_center,
        :gto_overlap,
        :gto_kinetic,
        :gto_electron_nuclear_by_center,
    )
        @test term in report.cpb_local_gto_included_terms
    end
    @test report.cpb_local_direct_core_support_count == 125
    @test report.mixed_gausslet_gto_overlap_available
    @test report.mixed_gausslet_gto_kinetic_available
    @test report.mixed_gausslet_gto_nuclear_by_center_available
    @test isnothing(report.mixed_gausslet_gto_blocker)
    @test report.mixed_gausslet_gto_blocker_source ==
          :_cpb_mixed_gto_parent_axis_representations
    @test report.mixed_gausslet_gto_rejected_helper == :unavailable
    @test report.mixed_gausslet_gto_blocker_detail ==
          :mapped_working_gaussian_proxy_axis_representation_available
    @test report.mixed_gausslet_gto_required_source ==
          :mapped_ordinary_working_gaussian_proxy_axis_representation
    @test report.gto_gto_overlap_available
    @test report.gto_gto_kinetic_available
    @test report.gto_gto_nuclear_by_center_available
    @test report.mixed_gausslet_gto_shapes.overlap == (125, 3)
    @test report.mixed_gausslet_gto_shapes.kinetic == (125, 3)
    @test report.mixed_gausslet_gto_shapes.nuclear_by_center == ((125, 3),)
    @test report.gto_gto_shapes.overlap == (3, 3)
    @test report.gto_gto_shapes.kinetic == (3, 3)
    @test report.gto_gto_shapes.nuclear_by_center == ((3, 3),)
    @test report.route_global_combined_basis_layout_status ==
          :available_route_global_combined_gto_basis_layout
    @test isnothing(report.route_global_combined_basis_layout_blocker)
    @test report.route_global_combined_basis_layout_kind ==
          :decomposed_wl_gausslet_plus_gto_supplement
    @test report.gausslet_retained_range == 1:223
    @test report.gausslet_retained_dimension == 223
    @test report.gto_supplement_range == 224:226
    @test report.gto_supplement_orbital_count == 3
    @test report.total_combined_dimension == 226
    @test report.combined_basis_dimension == 226
    @test report.combined_block_layout_keys ==
          (:gausslet_gausslet, :gausslet_gto, :gto_gausslet, :gto_gto)
    @test report.gausslet_gausslet_block_layout.row_range == 1:223
    @test report.gausslet_gausslet_block_layout.column_range == 1:223
    @test report.gausslet_gto_block_layout.row_range == 1:223
    @test report.gausslet_gto_block_layout.column_range == 224:226
    @test report.gto_gausslet_block_layout.row_range == 224:226
    @test report.gto_gausslet_block_layout.column_range == 1:223
    @test report.gto_gto_block_layout.row_range == 224:226
    @test report.gto_gto_block_layout.column_range == 224:226
    @test report.mixed_cpb_gto_blocks_orientation ==
          :gausslet_rows_by_gto_columns
    @test report.gto_gto_blocks_orientation == :gto_rows_by_gto_columns
    @test report.final_combined_overlap_layout_available
    @test report.final_combined_hamiltonian_layout_available
    @test !report.route_global_combined_overlap_matrix_materialized
    @test !report.route_global_combined_hamiltonian_matrix_materialized
    @test !report.gto_route_global_blocks_materialized
    @test !report.gto_hamiltonian_assembly_materialized
    @test !report.full_parent_window_cpb_used
    @test !report.direct_cartesian_product_assembly_used
    @test !report.ordinary_cartesian_ida_operators_used
    @test !report.exports_or_artifacts
end

@testset "decomposed WL gausslet-only H2+ acceptance" begin
    report = _wl_decomposed_h2plus_acceptance_report()
    println("decomposed WL gausslet-only H2+ acceptance: ", report)

    @test report.status == :materialized_decomposed_wl_h2plus_acceptance
    @test isnothing(report.blocker)
    @test report.acceptance_suite ==
          :decomposed_wl_gausslet_one_electron_acceptance
    @test report.acceptance_fixture == :h2plus
    @test report.q == 5
    @test report.ns == 5
    @test report.n_s == 5
    @test report.bond_length == 2.0
    @test report.center_count == 2
    @test report.center_keys == (:proton_a, :proton_b)
    @test report.center_locations == ((0.0, 0.0, -1.0), (0.0, 0.0, 1.0))
    @test report.nuclear_charges == (1.0, 1.0)
    @test report.nuclear_repulsion == 0.5
    @test report.decomposed_overlap_available
    @test report.decomposed_kinetic_available
    @test report.decomposed_electron_nuclear_by_center_available
    @test report.route_global_overlap_adapter_available
    @test report.route_global_kinetic_adapter_available
    @test report.route_global_electron_nuclear_by_center_adapter_available
    @test report.route_global_overlap_status ==
          :materialized_route_global_overlap_matrix
    @test report.route_global_kinetic_status ==
          :materialized_route_global_kinetic_matrix
    @test report.route_global_overlap_matrix_materialized
    @test report.route_global_kinetic_matrix_materialized
    @test report.decomposed_unit_pair_inventory_status ==
          :available_white_lindsey_decomposed_unit_pair_inventory
    @test isnothing(report.decomposed_unit_pair_inventory_blocker)
    @test report.decomposed_unit_pair_inventory_source_kind ==
          :white_lindsey_low_order_materialized_seed_ranges
    @test report.route_owned_decomposed_unit_pair_inventory_available
    @test report.decomposed_unit_count == 27
    @test report.decomposed_unit_pair_count == 378
    @test first(report.decomposed_unit_key_sample) == :white_lindsey_seed_direct_core
    @test report.retained_unit_column_ranges_materialized
    @test report.decomposed_unit_pair_column_ranges_available
    @test report.decomposed_unit_pair_inventory_retained_dimension == 223
    @test report.route_global_by_center_acceptance_matrix_available
    @test report.route_global_by_center_acceptance_matrix_status ==
          :materialized_route_global_electron_nuclear_by_center_matrix_set
    @test isnothing(report.route_global_by_center_acceptance_matrix_blocker)
    @test report.route_global_by_center_matrix_count == 2
    @test report.route_global_by_center_center_count == 2
    @test report.route_global_by_center_center_keys == (:proton_a, :proton_b)
    @test !report.route_global_by_center_nuclear_charge_applied
    @test !report.route_global_by_center_centers_summed
    @test report.hamiltonian_status ==
          :materialized_white_lindsey_decomposed_one_electron_hamiltonian
    @test report.hamiltonian_matrix_materialized
    @test report.nuclear_charge_application_stage ==
          :white_lindsey_hamiltonian_assembly
    @test report.nuclear_charge_applied_at_hamiltonian_assembly
    @test report.center_summation_stage == :white_lindsey_hamiltonian_assembly
    @test report.centers_summed_at_hamiltonian_assembly
    @test report.retained_dimension == 223
    @test report.overlap_matrix_dimension == report.retained_dimension
    @test report.hamiltonian_matrix_dimension == report.retained_dimension
    @test report.overlap_minimum_eigenvalue > 0.0
    @test report.overlap_condition_estimate < 1.0000001
    @test report.overlap_symmetry_error < 1.0e-12
    @test report.overlap_near_zero_eigenvalue_count == 0
    @test report.overlap_negative_eigenvalue_count == 0
    @test report.overlap_rank_estimate == report.retained_dimension
    @test report.decomposed_unit_column_range_span == 1:223
    @test report.missing_interior_retained_column_count == 0
    @test !report.boundary_inventory_only
    @test report.solve_status == :materialized_decomposed_wl_one_electron_solve
    @test isnothing(report.solve_blocker)
    @test report.solve_kind == :ordinary_symmetric
    @test report.electronic_energy > report.exact_electronic_energy
    @test report.total_bo_energy > report.exact_total_bo_energy
    @test report.total_bo_energy < -0.52
    @test isapprox(
        report.total_bo_energy,
        -0.533841728044377;
        atol = 1.0e-10,
    )
    @test !report.fixed_block_operator_matrices_used
    @test report.decomposed_wl_units_consumed
    @test !report.full_parent_window_cpb_used
    @test !report.direct_cartesian_product_assembly_used
    @test !report.ordinary_cartesian_ida_operators_used
    @test report.acceptance_energy_materialized
    @test report.h2plus_acceptance_active
    @test report.acceptance_fixture_active
    @test report.elapsed_seconds >= 0.0
end
