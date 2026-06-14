using Test
using JLD2

const _H2_FAKE_PQS_DRIVER =
    normpath(joinpath(@__DIR__, "..", "..", "bin", "cartesian_ham_builder.jl"))
const _H2_FAKE_PQS_INPUT =
    normpath(joinpath(@__DIR__, "..", "driver_inputs", "h2_fake_pqs_q5_wl_source_backed_r4.jl"))

@testset "cartesian_ham_builder_h2_fake_pqs_wl_source_backed_r4 reproduction artifact" begin
    mktempdir() do dir
        outfile = joinpath(dir, "h2_fake_pqs_q5_wl_source_backed_r4.jld2")
        tsvfile = joinpath(dir, "h2_fake_pqs_q5_wl_source_backed_r4.tsv")
        saved_args = copy(ARGS)
        empty!(ARGS)
        append!(ARGS, [_H2_FAKE_PQS_INPUT, "outfile=$(repr(outfile))", "tsvfile=$(repr(tsvfile))"])
        try
            include(_H2_FAKE_PQS_DRIVER)
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
            @test file["comparison/ready"] == true
            @test file["comparison/role"] === :fake_pqs_wl_reproduction
            @test file["comparison/blocker"] === nothing
            @test file["comparison/reference_label"] ==
                  "WL/QW H2 R=4 gausslet-only 463 reproduction"
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
            @test file["comparison/wl_reference_final_dimension"] == 463
            @test file["comparison/wl_reference_retained_transform_kind"] ===
                  :white_lindsey_old_qw_gausslet_retained_transform
            @test file["comparison/old_supplemented_wl_qw_scalar_references_blocked"] ==
                  true
            @test file["fake_pqs/enabled"] == true
            @test file["fake_pqs/source"] === :source_backed_fixed_source_oracle
            @test file["fake_pqs/reason"] ===
                  :wl_qw_fixed_source_retained_transform_imported
            @test file["fake_pqs/independent_pqs_transform"] == false
            @test file["fake_pqs/temporary_object"] == true
            @test file["fake_pqs/delete_after_independent_pqs"] == true
            @test file["fake_pqs/warning"] ===
                  :retained_transform_imported_from_wl_qw_fixed_source_oracle

            @test file["target/status"] ===
                  :available_physical_gausslet_core_shell_target_inventory
            @test file["target/blocker"] === nothing
            @test Tuple(file["target/support_counts"]) == (275, 578, 362)
            @test Tuple(file["target/retained_counts"]) == (251, 98, 114)
            @test Tuple(file["target/retained_order"]) ==
                  (:atom_contact_core, :shared_shell_1, :shared_shell_2)
            @test file["target/expected_final_dimension"] == 463
            @test file["target/retained_atom_core_interiors"] == true
            @test file["target/source_plan_status"] ===
                  :available_pqs_diatomic_physical_gausslet_core_shell_source_plan
            @test file["target/source_plan_authority_status"] ===
                  :fake_pqs_private_source_backed_adapter_authority

            route_fingerprint = (;
                role = file["route/artifact_role"],
                source = file["route/source_plan_status"],
                final_basis = file["route/final_basis_status"],
                h1 = file["route/h1_status"],
                h1j = file["route/h1_j_status"],
                rhf_input = file["route/private_rhf_input_contract_status"],
                fake = file["route/fake_pqs_enabled"],
                independent = file["route/independent_pqs_transform"],
            )
            @test route_fingerprint == (;
                role = :fake_pqs_source_backed_wl_reproduction,
                source = :available_pqs_diatomic_physical_gausslet_core_shell_source_plan,
                final_basis = :available_pqs_physical_gausslet_final_basis,
                h1 = :materialized_pqs_physical_gausslet_h1_solve,
                h1j = :materialized_pqs_physical_gausslet_h1_j_payload,
                rhf_input = :available_pqs_physical_gausslet_rhf_input_contract,
                fake = true,
                independent = false,
            )
            @test file["route/comparison_role"] === :fake_pqs_wl_reproduction
            @test file["route/source"] === :source_backed_fixed_source_oracle
            @test file["route/warning"] ===
                  :retained_transform_imported_from_wl_qw_fixed_source_oracle

            @test file["basis/retained_atom_core_interiors"] == true
            @test file["basis/source_plan_role"] ===
                  :fake_pqs_source_backed_fixed_source_adapter
            @test file["basis/final_dimension"] == 463
            @test file["basis/final_overlap_identity_error"] < 1e-10

            @test file["physics/endpoint_role"] === :fake_pqs_wl_reproduction
            @test file["physics/independent_pqs_transform"] == false
            @test file["physics/endpoint_ready"] == false
            @test file["physics/endpoint_blocker"] ===
                  :fake_pqs_source_backed_wl_reproduction_not_independent_pqs
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
            @test file["density_interaction/support_raw_pair_shape"] == (1215, 1215)
            @test file["density_interaction/pre_final_pair_matrix_shape"] == (463, 463)
            @test file["density_interaction/pre_final_pair_matrix_finite"] == true
            @test file["density_interaction/pre_final_pair_matrix_symmetry_error"] < 1e-8
            @test isfinite(file["density_interaction/h1_j_self_coulomb"])
            @test (
                file["private_rhf/input_contract_status"],
                file["private_rhf/electron_count"],
                file["private_rhf/occupation_policy"],
            ) == (
                :available_pqs_physical_gausslet_rhf_input_contract,
                2,
                :closed_shell_rhf,
            )
            @test file["private_rhf/executed"] == true
            @test file["private_rhf/iteration_count"] > 0
            @test file["private_rhf/materialized"] == true
            @test file["private_rhf/converged"] == true
            @test isfinite(file["private_rhf/total_energy"])
            @test file["comparison/pqs_rhf_total_with_nuclear_repulsion"] ≈
                  file["private_rhf/total_energy"] +
                  file["comparison/wl_rhf_nuclear_repulsion"]
            @test file["private_rhf/requested"] == true
        end
    end

    mktempdir() do dir
        outfile = joinpath(dir, "h2_fake_pqs_q5_wl_source_backed_r4_mwg_preflight.jld2")
        saved_args = copy(ARGS)
        empty!(ARGS)
        append!(
            ARGS,
            [
                _H2_FAKE_PQS_INPUT,
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
            include(_H2_FAKE_PQS_DRIVER)
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
            @test (
                file["supplement_request/status"],
                file["supplement_request/supplement_policy"],
                Tuple(file["supplement_request/missing_fact_labels"]),
                file["supplement_request/fixture_label"],
                file["supplement_request/basis_name"],
                file["supplement_request/lmax"],
                Tuple(file["supplement_request/atom_symbols"]),
                Tuple(file["supplement_request/nuclear_charges"]),
                file["supplement_request/bond_axis"],
            ) == (
                :available_pqs_physical_gausslet_supplement_request,
                :mwg_residual_gto,
                (:missing_provider_gto_supplement_blocks,),
                :h2_r4_physical_gausslet_q5,
                "H/cc-pVTZ",
                1,
                ("H", "H"),
                (1, 1),
                :z,
            )
            @test file["supplement_request/bond_length"] ≈ 4.0
            @test file["supplement_representation/status"] ===
                  :available_pqs_physical_gausslet_gto_supplement_representation
            @test file["supplement_representation/orbital_count"] == 18
            @test file["supplement_preflight/status"] ===
                  :blocked_pqs_physical_gausslet_mwg_residual_gto_preflight
            @test file["supplement_preflight/blocker"] ===
                  :missing_provider_gto_supplement_blocks
            @test Tuple(file["supplement_preflight/support_counts"]) ==
                  (275, 578, 362)
            @test Tuple(file["supplement_preflight/retained_counts"]) ==
                  (251, 98, 114)
            @test file["supplement_preflight/gausslet_final_dimension"] == 463
            @test Tuple(file["supplement_preflight/required_fact_labels"]) ==
                  (
                      :gto_supplement_representation,
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
            @test file["route/source_plan_status"] ===
                  :blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan
        end
    end
end
