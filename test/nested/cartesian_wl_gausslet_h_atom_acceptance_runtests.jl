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
        h2plus_acceptance_active = false,
        h2plus_acceptance_blocker =
            :deferred_until_single_center_decomposed_wl_acceptance_review,
    )
end

@testset "decomposed WL gausslet-only H atom acceptance" begin
    report = _wl_decomposed_h_atom_acceptance_report()
    println("decomposed WL gausslet-only H atom acceptance: ", report)

    @test report.status == :materialized_decomposed_wl_h_atom_acceptance
    @test isnothing(report.blocker)
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
    @test !report.h2plus_acceptance_active
    @test report.h2plus_acceptance_blocker ==
          :deferred_until_single_center_decomposed_wl_acceptance_review
end
