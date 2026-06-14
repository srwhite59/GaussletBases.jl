using Test
using JLD2

const _H2_PHYSICAL_PQS_DRIVER =
    normpath(joinpath(@__DIR__, "..", "..", "bin", "cartesian_ham_builder.jl"))
const _H2_PHYSICAL_PQS_INPUT =
    normpath(joinpath(@__DIR__, "..", "driver_inputs", "h2_pqs_q5_physical_gausslet_r4.jl"))

@testset "cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4 target artifact" begin
    mktempdir() do dir
        outfile = joinpath(dir, "h2_pqs_q5_physical_gausslet_r4.jld2")
        tsvfile = joinpath(dir, "h2_pqs_q5_physical_gausslet_r4.tsv")
        saved_args = copy(ARGS)
        empty!(ARGS)
        append!(ARGS, [_H2_PHYSICAL_PQS_INPUT, "outfile=$(repr(outfile))", "tsvfile=$(repr(tsvfile))"])
        try
            include(_H2_PHYSICAL_PQS_DRIVER)
        finally
            empty!(ARGS)
            append!(ARGS, saved_args)
        end

        @test isfile(outfile)
        jldopen(outfile, "r") do file
            @test file["config/route_kind"] ===
                  :bond_aligned_diatomic_physical_gausslet_core_shell_pqs
            @test file["config/supplement_policy"] === :none
            @test file["config/comparison_ready"] == true
            @test file["config/run_final_basis"] == true
            @test file["supplement_preflight/status"] === :not_requested
            @test file["supplement_preflight/blocker"] === nothing
            @test file["supplement_preflight/retained_transform_kind"] === :pqs
            @test file["supplement_preflight/gausslet_final_dimension"] == 463
            @test file["supplement_preflight/supplement_policy"] === :none
            @test Tuple(file["supplement_preflight/required_fact_labels"]) == ()
            @test Tuple(file["supplement_preflight/missing_fact_labels"]) == ()
            @test file["supplement_preflight/matrices_materialized"] == false
            @test file["supplement_preflight/supplemented_values_materialized"] == false
            @test file["comparison/ready"] == true
            @test file["comparison/blocker"] === nothing
            @test file["comparison/reference_label"] ==
                  "WL/QW H2 R=4 gausslet-only 463"
            @test file["comparison/wl_h1_lowest"] ≈ -0.7946609179724673
            @test abs(file["comparison/delta_h1"]) < 1.0e-10
            @test file["comparison/wl_h1_self_coulomb"] ≈ 0.45696639804337047
            @test abs(file["comparison/delta_h1_j"]) < 1.0e-10
            @test file["comparison/wl_rhf_electronic_energy"] ≈
                  -1.1589518556651142
            @test abs(file["comparison/delta_rhf_electronic_energy"]) < 1.0e-9
            @test file["comparison/wl_rhf_nuclear_repulsion"] ≈ 0.25
            @test file["comparison/wl_rhf_total_with_nuclear_repulsion"] ≈
                  -0.9089518556651142
            @test abs(file["comparison/delta_rhf_total_with_nuclear_repulsion"]) <
                  1.0e-9
            @test file["comparison/wl_reference_candidate_status"] ===
                  :available_wl_h2_gausslet_only_reference_candidate
            @test file["comparison/wl_reference_candidate_blocker"] === nothing
            @test file["comparison/wl_reference_final_dimension"] == 463
            @test file["comparison/wl_reference_retained_transform_kind"] ===
                  :white_lindsey_old_qw_gausslet_retained_transform
            @test file["comparison/wl_reference_supplement_policy"] === :none
            @test file["comparison/wl_reference_label"] ==
                  "WL/QW H2 R=4 gausslet-only 463"
            @test file["comparison/old_supplemented_wl_qw_scalar_references_blocked"] ==
                  true

            @test file["target/status"] ===
                  :available_physical_gausslet_core_shell_target_inventory
            @test file["target/blocker"] === nothing
            @test Tuple(file["target/support_units"]) ==
                  (:atom_contact_core, :shared_shell_1, :shared_shell_2)
            @test Tuple(file["target/support_counts"]) == (275, 578, 362)
            @test Tuple(file["target/retained_units"]) ==
                  (:atom_contact_core, :shared_shell_1, :shared_shell_2)
            @test Tuple(file["target/retained_counts"]) == (251, 98, 114)
            @test Tuple(file["target/retained_order"]) ==
                  (:atom_contact_core, :shared_shell_1, :shared_shell_2)
            @test file["target/expected_final_dimension"] == 463
            @test file["target/retained_atom_core_interiors"] == true
            @test file["target/source_plan_role"] ===
                  :atom_contact_core_plus_pqs_shared_shells
            @test file["target/source_plan_status"] ===
                  :available_pqs_diatomic_physical_gausslet_core_shell_source_plan

            @test file["route/artifact_role"] ===
                  :physical_gausslet_endpoint_target
            @test file["route/source_plan_status"] ===
                  :available_pqs_diatomic_physical_gausslet_core_shell_source_plan
            @test file["route/final_basis_status"] ===
                  :available_pqs_physical_gausslet_final_basis
            @test file["route/h1_status"] ===
                  :materialized_pqs_physical_gausslet_h1_solve
            @test file["route/h1_materialized"] == true
            @test file["route/h1_j_status"] ===
                  :materialized_pqs_physical_gausslet_h1_j_payload
            @test file["route/h1_j_materialized"] == true
            @test file["route/private_rhf_input_contract_status"] ===
                  :available_pqs_physical_gausslet_rhf_input_contract

            @test file["basis/retained_atom_core_interiors"] == true
            @test file["basis/source_plan_role"] ===
                  :atom_contact_core_plus_pqs_shared_shells
            @test file["basis/final_dimension"] == 463
            @test file["basis/final_overlap_identity_error"] < 1e-10

            @test file["physics/endpoint_ready"] == true
            @test isfinite(file["physics/h1_lowest"])
            @test file["physics/h1_lowest"] < 0
            @test file["physics/h1_hamiltonian_matrix_finite"] == true
            @test file["physics/h1_hamiltonian_symmetry_error"] < 1e-8
            @test file["density_interaction/status"] ===
                  :materialized_pqs_physical_gausslet_pre_final_density_interaction
            @test file["density_interaction/density_gauge"] ===
                  :pre_final_localized_positive_weight
            @test file["density_interaction/raw_pair_factor_convention"] ===
                  :raw_numerator
            @test file["density_interaction/support_weight_count"] == 1215
            @test file["density_interaction/support_weights_all_positive"] == true
            @test file["density_interaction/support_raw_pair_shape"] == (1215, 1215)
            @test file["density_interaction/support_raw_pair_finite"] == true
            @test file["density_interaction/pre_final_pair_matrix_shape"] == (463, 463)
            @test file["density_interaction/pre_final_pair_matrix_finite"] == true
            @test file["density_interaction/pre_final_pair_matrix_symmetry_error"] < 1e-8
            @test isfinite(file["density_interaction/h1_j_self_coulomb"])
            @test file["private_rhf/input_contract_status"] ===
                  :available_pqs_physical_gausslet_rhf_input_contract
            @test file["private_rhf/input_contract_blocker"] === nothing
            @test file["private_rhf/input_contract_available"] == true
            @test file["private_rhf/electron_count"] == 2
            @test file["private_rhf/occupation_policy"] === :closed_shell_rhf
            @test file["private_rhf/occupation_nocc"] == 1
            @test file["private_rhf/h1_matrix_available"] == true
            @test file["private_rhf/h1_matrix_finite"] == true
            @test file["private_rhf/h1_matrix_symmetry_error"] < 1e-8
            @test file["private_rhf/density_interaction_available"] == true
            @test file["private_rhf/final_to_pre_final_transform_available"] == true
            @test file["private_rhf/pre_final_pair_matrix_available"] == true
            @test file["private_rhf/pre_final_pair_matrix_finite"] == true
            @test file["private_rhf/pre_final_pair_matrix_symmetry_error"] < 1e-8
            @test file["private_rhf/executed"] == true
            @test file["private_rhf/iteration_count"] > 0
            if file["private_rhf/materialized"]
                @test file["route/private_rhf_execution_status"] ===
                      :materialized_pqs_physical_gausslet_private_rhf_execution
                @test file["route/private_rhf_materialized"] == true
                @test file["physics/endpoint_blocker"] === nothing
                @test file["private_rhf/execution_status"] ===
                      :materialized_pqs_physical_gausslet_private_rhf_execution
                @test file["private_rhf/execution_blocker"] === nothing
                @test file["private_rhf/converged"] == true
                @test isfinite(file["private_rhf/total_energy"])
                @test isfinite(file["private_rhf/one_body_energy"])
                @test isfinite(file["private_rhf/two_body_energy"])
                @test file["private_rhf/final_density_one_step_consistency_status"] ===
                      :reviewed_recomputed
                @test file["comparison/pqs_rhf_total_with_nuclear_repulsion"] ≈
                      file["private_rhf/total_energy"] +
                      file["comparison/wl_rhf_nuclear_repulsion"]
            else
                @test file["route/private_rhf_execution_status"] ===
                      :blocked_pqs_physical_gausslet_private_rhf_execution
                @test file["route/private_rhf_materialized"] == false
                @test file["physics/endpoint_blocker"] ===
                      :physical_gausslet_private_rhf_not_converged
                @test file["private_rhf/execution_status"] ===
                      :blocked_pqs_physical_gausslet_private_rhf_execution
                @test file["private_rhf/execution_blocker"] ===
                      :physical_gausslet_private_rhf_not_converged
                @test file["private_rhf/converged"] == false
                @test file["private_rhf/total_energy"] === nothing ||
                      isfinite(file["private_rhf/total_energy"])
            end
            @test file["private_rhf/requested"] == true
        end
    end

    mktempdir() do dir
        outfile = joinpath(dir, "h2_pqs_q5_physical_gausslet_r4_mwg_preflight.jld2")
        saved_args = copy(ARGS)
        empty!(ARGS)
        append!(
            ARGS,
            [
                _H2_PHYSICAL_PQS_INPUT,
                "outfile=$(repr(outfile))",
                "save_tsv=false",
                "supplement_policy=:mwg_residual_gto",
                "comparison_ready=false",
                "comparison_blocker=:mwg_residual_gto_preflight_only",
                "run_final_basis=false",
                "run_h1=false",
                "run_h1_j=false",
                "run_private_rhf=false",
            ],
        )
        try
            include(_H2_PHYSICAL_PQS_DRIVER)
        finally
            empty!(ARGS)
            append!(ARGS, saved_args)
        end

        @test isfile(outfile)
        jldopen(outfile, "r") do file
            @test file["config/supplement_policy"] === :mwg_residual_gto
            @test file["route/supplement_preflight_status"] ===
                  :blocked_pqs_physical_gausslet_mwg_residual_gto_preflight
            @test file["route/supplement_preflight_blocker"] ===
                  :missing_provider_gto_supplement_blocks
            @test file["supplement_preflight/status"] ===
                  :blocked_pqs_physical_gausslet_mwg_residual_gto_preflight
            @test file["supplement_preflight/blocker"] ===
                  :missing_provider_gto_supplement_blocks
            @test file["supplement_preflight/fixture_label"] ===
                  :h2_r4_physical_gausslet_q5
            @test Tuple(file["supplement_preflight/support_counts"]) ==
                  (275, 578, 362)
            @test Tuple(file["supplement_preflight/retained_counts"]) ==
                  (251, 98, 114)
            @test Tuple(file["supplement_preflight/retained_order"]) ==
                  (:atom_contact_core, :shared_shell_1, :shared_shell_2)
            @test file["supplement_preflight/retained_transform_kind"] === :pqs
            @test file["supplement_preflight/gausslet_final_dimension"] == 463
            @test file["supplement_preflight/supplement_policy"] ===
                  :mwg_residual_gto
            @test Tuple(file["supplement_preflight/available_fact_labels"]) ==
                  (
                      :physical_gausslet_core_shell_target_inventory,
                      :pqs_retained_transform_kind,
                      :gausslet_only_final_dimension,
                  )
            @test Tuple(file["supplement_preflight/required_fact_labels"]) ==
                  (
                      :provider_gto_supplement_blocks,
                      :mixed_gausslet_gto_blocks,
                      :gto_gto_blocks,
                      :combined_raw_moment_matrices,
                      :residual_mwg_representation,
                      :combined_density_density_readiness,
                  )
            @test Tuple(file["supplement_preflight/missing_fact_labels"]) ==
                  (
                      :missing_provider_gto_supplement_blocks,
                      :missing_mixed_gausslet_gto_blocks,
                      :missing_gto_gto_blocks,
                      :missing_combined_raw_moment_matrices,
                      :missing_residual_mwg_representation,
                      :missing_combined_density_density_readiness,
                  )
            @test file["supplement_preflight/matrices_materialized"] == false
            @test file["supplement_preflight/supplemented_values_materialized"] == false
            @test file["route/source_plan_status"] ===
                  :blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan
        end
    end
end
