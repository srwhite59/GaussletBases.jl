# Runtime role: decomposed White-Lindsey He acceptance readiness audit.
#
# This file exercises the gausslet-only one-center Z = 2 decomposed WL
# one-electron path and then stops at the current missing decomposed
# two-electron density-density/IDA route. It must not fall back to full-parent
# CPBs, direct Cartesian product operators, or ordinary Cartesian IDA operators.

using Test
using LinearAlgebra
using GaussletBases

include("cartesian_white_lindsey_adapter_fixture_helpers.jl")

const WLHeAcceptanceCPBM = GaussletBases.CartesianPairBlockMaterialization

function _wl_he_parent_axis_inputs()
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

function _wl_he_center_records()
    return (
        (;
            center_key = :helium_nucleus,
            center_index = 1,
            nuclear_charge = 2.0,
            location = (0.0, 0.0, 0.0),
        ),
    )
end

function _wl_he_overlap_metric_diagnostics(overlap_matrix, decomposed_inventory)
    sym_s = Symmetric((overlap_matrix + transpose(overlap_matrix)) ./ 2)
    eigenvalues = eigvals(sym_s)
    minimum_eigenvalue = minimum(eigenvalues)
    maximum_eigenvalue = maximum(eigenvalues)
    return (;
        blocker = minimum_eigenvalue <= 0.0 ?
            :decomposed_wl_overlap_metric_not_positive_definite :
            nothing,
        minimum_eigenvalue,
        maximum_eigenvalue,
        condition_estimate = maximum_eigenvalue / minimum_eigenvalue,
        symmetry_error = maximum(abs.(overlap_matrix - transpose(overlap_matrix))),
        rank_estimate = count(eigenvalues .> 1.0e-10),
        retained_dimension = decomposed_inventory.retained_dimension,
    )
end

function _wl_he_lowest_one_electron_energy(
    hamiltonian_matrix,
    overlap_matrix,
    overlap_diagnostics,
)
    if !isnothing(overlap_diagnostics.blocker)
        return (;
            status = :blocked_decomposed_wl_he_one_electron_solve,
            blocker = overlap_diagnostics.blocker,
            energy = NaN,
            solve_kind = :blocked_overlap_metric,
        )
    end
    overlap_identity =
        isapprox(overlap_matrix, Matrix{Float64}(I, size(overlap_matrix)...);
            atol = 1.0e-9,
            rtol = 1.0e-9,
        )
    sym_h = Symmetric((hamiltonian_matrix + transpose(hamiltonian_matrix)) ./ 2)
    values = overlap_identity ?
        eigvals(sym_h) :
        eigen(
            sym_h,
            Symmetric((overlap_matrix + transpose(overlap_matrix)) ./ 2),
        ).values
    return (;
        status = :materialized_decomposed_wl_he_one_electron_solve,
        blocker = nothing,
        energy = minimum(values),
        solve_kind =
            overlap_identity ? :ordinary_symmetric : :generalized_symmetric,
        overlap_identity,
    )
end

function _wl_he_density_density_route_audit(seed_report)
    route_global_supported_terms =
        WLHeAcceptanceCPBM._route_global_one_body_supported_terms()
    route_function_available =
        isdefined(WLHeAcceptanceCPBM, :route_global_decomposed_wl_density_density_matrix)
    seed_electron_electron_materialized =
        hasproperty(seed_report, :electron_electron_materialized) ?
        seed_report.electron_electron_materialized :
        false
    available = route_function_available && seed_electron_electron_materialized
    return (;
        status =
            available ?
            :available_decomposed_wl_density_density_interaction_route :
            :blocked_decomposed_wl_density_density_interaction_route,
        blocker =
            available ? nothing : :missing_decomposed_wl_density_density_interaction_route,
        route_global_density_density_function_available = route_function_available,
        seed_electron_electron_materialized,
        electron_electron_term_advertised =
            :electron_electron_density_density in route_global_supported_terms,
        ordinary_cartesian_ida_operators_used = false,
        direct_cartesian_product_assembly_used = false,
        full_parent_window_cpb_used = false,
    )
end

function _wl_decomposed_he_atom_acceptance_audit()
    seed_report = GaussletBases._white_lindsey_low_order_materialized_seed_report()
    axis_inputs = _wl_he_parent_axis_inputs()
    center_records = _wl_he_center_records()
    decomposed_inventory =
        WLHeAcceptanceCPBM.white_lindsey_decomposed_unit_pair_inventory(
            seed_report;
            metadata = (;
                parent_axis_counts = axis_inputs.parent_axis_counts,
                parent_axis_bundle_object = axis_inputs.parent_axis_bundle_object,
            ),
        )
    overlap_global = nothing
    kinetic_global = nothing
    by_center_global = nothing
    hamiltonian = nothing
    overlap_diagnostics = nothing
    solve = nothing
    density_density_audit = nothing
    elapsed_seconds = @elapsed begin
        overlap_global =
            WLHeAcceptanceCPBM.route_global_decomposed_wl_overlap_matrix(
                seed_report;
                parent_axis_counts = axis_inputs.parent_axis_counts,
                parent_axis_bundle_object =
                    axis_inputs.parent_axis_bundle_object,
                overlap_1d = axis_inputs.overlap_1d,
            )
        kinetic_global =
            WLHeAcceptanceCPBM.route_global_decomposed_wl_kinetic_matrix(
                seed_report;
                parent_axis_counts = axis_inputs.parent_axis_counts,
                parent_axis_bundle_object =
                    axis_inputs.parent_axis_bundle_object,
                overlap_1d = axis_inputs.overlap_1d,
                kinetic_1d = axis_inputs.kinetic_1d,
            )
        by_center_global =
            WLHeAcceptanceCPBM.route_global_electron_nuclear_by_center_matrices(
                seed_report;
                parent_axis_counts = axis_inputs.parent_axis_counts,
                parent_axis_bundle_object =
                    axis_inputs.parent_axis_bundle_object,
                coulomb_expansion = coulomb_gaussian_expansion(doacc = false),
                center_records,
            )
        hamiltonian =
            WLHeAcceptanceCPBM.white_lindsey_decomposed_one_electron_hamiltonian(
                kinetic_global,
                by_center_global,
            )
        overlap_diagnostics =
            _wl_he_overlap_metric_diagnostics(
                overlap_global.matrix,
                decomposed_inventory,
            )
        solve = _wl_he_lowest_one_electron_energy(
            hamiltonian.matrix,
            overlap_global.matrix,
            overlap_diagnostics,
        )
        density_density_audit = _wl_he_density_density_route_audit(seed_report)
    end
    two_electron_route_available =
        density_density_audit.status ===
        :available_decomposed_wl_density_density_interaction_route
    one_electron_contribution = 2.0 * solve.energy
    return (;
        object_kind = :decomposed_wl_he_atom_acceptance_audit,
        acceptance_suite = :decomposed_wl_gausslet_two_electron_acceptance,
        acceptance_fixture = :he_atom,
        status =
            two_electron_route_available ?
            :ready_decomposed_wl_he_atom_two_electron_acceptance :
            :blocked_decomposed_wl_he_atom_two_electron_acceptance,
        blocker = density_density_audit.blocker,
        q = 5,
        ns = 5,
        n_s = 5,
        center_count = length(center_records),
        center_keys = Tuple(record.center_key for record in center_records),
        center_locations = Tuple(record.location for record in center_records),
        nuclear_charges = Tuple(record.nuclear_charge for record in center_records),
        core_spacing = axis_inputs.spacing_parameter,
        distortion_strength = axis_inputs.distortion_strength,
        tail_spacing = axis_inputs.tail_spacing,
        parent_axis_counts = axis_inputs.parent_axis_counts,
        decomposed_unit_count = decomposed_inventory.unit_count,
        decomposed_unit_pair_count = decomposed_inventory.pair_count,
        retained_dimension = decomposed_inventory.retained_dimension,
        route_global_overlap_status = overlap_global.status,
        route_global_kinetic_status = kinetic_global.status,
        route_global_by_center_status = by_center_global.status,
        route_global_by_center_center_count = by_center_global.center_count,
        route_global_by_center_nuclear_charge_applied =
            by_center_global.nuclear_charge_applied,
        route_global_by_center_centers_summed = by_center_global.centers_summed,
        hamiltonian_status = hamiltonian.status,
        nuclear_charge_application_stage =
            hamiltonian.charge_application_stage,
        centers_summed_at_hamiltonian_assembly =
            hamiltonian.centers_summed_at_hamiltonian_assembly,
        overlap_minimum_eigenvalue = overlap_diagnostics.minimum_eigenvalue,
        overlap_maximum_eigenvalue = overlap_diagnostics.maximum_eigenvalue,
        overlap_condition_estimate = overlap_diagnostics.condition_estimate,
        overlap_symmetry_error = overlap_diagnostics.symmetry_error,
        overlap_rank_estimate = overlap_diagnostics.rank_estimate,
        one_electron_solve_status = solve.status,
        one_electron_solve_blocker = solve.blocker,
        one_electron_solve_kind = solve.solve_kind,
        lowest_one_electron_orbital_energy = solve.energy,
        closed_shell_one_electron_energy_contribution = one_electron_contribution,
        electron_electron_status = density_density_audit.status,
        electron_electron_blocker = density_density_audit.blocker,
        electron_electron_contribution = :unavailable,
        total_he_energy = :unavailable,
        exact_he_reference_energy = -2.9037243770341196,
        total_he_energy_distance_from_exact = :unavailable,
        ida_scf_energy_interpretation_available = two_electron_route_available,
        elapsed_seconds,
        decomposed_wl_units_consumed = true,
        full_parent_window_cpb_used = false,
        direct_cartesian_product_assembly_used = false,
        ordinary_cartesian_ida_operators_used = false,
        generalized_overlap_final_solve = false,
        gto_supplement_used = false,
        pqs_transforms_materialized = false,
        exports_or_artifacts = false,
        density_density_route_global_function_available =
            density_density_audit.route_global_density_density_function_available,
        seed_electron_electron_materialized =
            density_density_audit.seed_electron_electron_materialized,
        electron_electron_term_advertised =
            density_density_audit.electron_electron_term_advertised,
    )
end

@testset "decomposed WL gausslet-only He atom acceptance readiness" begin
    report = _wl_decomposed_he_atom_acceptance_audit()
    println("decomposed WL gausslet-only He atom acceptance readiness: ", report)

    @test report.status == :blocked_decomposed_wl_he_atom_two_electron_acceptance
    @test report.blocker ==
          :missing_decomposed_wl_density_density_interaction_route
    @test report.acceptance_suite ==
          :decomposed_wl_gausslet_two_electron_acceptance
    @test report.acceptance_fixture == :he_atom
    @test report.q == 5
    @test report.ns == 5
    @test report.center_count == 1
    @test report.center_keys == (:helium_nucleus,)
    @test report.nuclear_charges == (2.0,)
    @test report.retained_dimension == 223
    @test report.decomposed_unit_count == 27
    @test report.decomposed_unit_pair_count == 378
    @test report.route_global_overlap_status ==
          :materialized_route_global_overlap_matrix
    @test report.route_global_kinetic_status ==
          :materialized_route_global_kinetic_matrix
    @test report.route_global_by_center_status ==
          :materialized_route_global_electron_nuclear_by_center_matrix_set
    @test !report.route_global_by_center_nuclear_charge_applied
    @test !report.route_global_by_center_centers_summed
    @test report.hamiltonian_status ==
          :materialized_white_lindsey_decomposed_one_electron_hamiltonian
    @test report.nuclear_charge_application_stage ==
          :white_lindsey_hamiltonian_assembly
    @test report.centers_summed_at_hamiltonian_assembly
    @test report.overlap_minimum_eigenvalue > 0.0
    @test report.overlap_condition_estimate < 1.0000001
    @test report.overlap_rank_estimate == report.retained_dimension
    @test report.one_electron_solve_status ==
          :materialized_decomposed_wl_he_one_electron_solve
    @test isnothing(report.one_electron_solve_blocker)
    @test report.one_electron_solve_kind == :ordinary_symmetric
    @test isfinite(report.lowest_one_electron_orbital_energy)
    @test report.lowest_one_electron_orbital_energy < -1.0
    @test isfinite(report.closed_shell_one_electron_energy_contribution)
    @test report.electron_electron_status ==
          :blocked_decomposed_wl_density_density_interaction_route
    @test report.electron_electron_blocker ==
          :missing_decomposed_wl_density_density_interaction_route
    @test report.electron_electron_contribution == :unavailable
    @test report.total_he_energy == :unavailable
    @test !report.ida_scf_energy_interpretation_available
    @test !report.density_density_route_global_function_available
    @test !report.seed_electron_electron_materialized
    @test !report.electron_electron_term_advertised
    @test report.decomposed_wl_units_consumed
    @test !report.full_parent_window_cpb_used
    @test !report.direct_cartesian_product_assembly_used
    @test !report.ordinary_cartesian_ida_operators_used
    @test !report.generalized_overlap_final_solve
    @test !report.gto_supplement_used
    @test !report.pqs_transforms_materialized
    @test !report.exports_or_artifacts
    @test report.elapsed_seconds >= 0.0
end
