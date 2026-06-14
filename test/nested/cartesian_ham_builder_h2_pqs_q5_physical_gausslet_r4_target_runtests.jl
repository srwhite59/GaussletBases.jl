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
            @test file["config/comparison_ready"] == false
            @test file["config/run_final_basis"] == true

            @test file["parent/parent_axis_counts"] == (x = 9, y = 9, z = 15)

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
            @test file["target/source_plan_blocker"] === nothing
            @test file["target/source_plan_candidate_status"] ===
                  :available_physical_gausslet_source_plan_candidate
            @test file["target/source_plan_candidate_source"] ===
                  :source_backed_fixed_source_oracle
            @test file["target/source_plan_candidate_counts_match"] == true
            @test file["target/source_plan_authority_status"] ===
                  :private_source_backed_adapter_authority
            @test file["target/supplement_policy"] === :none

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

            @test file["physics/endpoint_ready"] == false
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
                @test file["physics/endpoint_blocker"] ===
                      :missing_h2_gausslet_only_reference_comparison
                @test file["private_rhf/execution_status"] ===
                      :materialized_pqs_physical_gausslet_private_rhf_execution
                @test file["private_rhf/execution_blocker"] === nothing
                @test file["private_rhf/converged"] == true
                @test isfinite(file["private_rhf/total_energy"])
                @test isfinite(file["private_rhf/one_body_energy"])
                @test isfinite(file["private_rhf/two_body_energy"])
                @test file["private_rhf/final_density_one_step_consistency_status"] ===
                      :reviewed_recomputed
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
            @test file["comparison/ready"] == false
            @test file["private_rhf/requested"] == true
        end
    end
end
