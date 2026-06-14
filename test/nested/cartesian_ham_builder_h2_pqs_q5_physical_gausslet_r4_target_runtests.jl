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
            @test file["config/run_final_basis"] == false

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
                  :blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan
            @test file["target/source_plan_blocker"] ===
                  :source_plan_candidate_not_route_authority
            @test file["target/source_plan_candidate_status"] ===
                  :available_physical_gausslet_source_plan_candidate
            @test file["target/source_plan_candidate_source"] ===
                  :source_backed_fixed_source_oracle
            @test file["target/source_plan_candidate_counts_match"] == true
            @test file["target/source_plan_authority_status"] ===
                  :candidate_not_route_authority
            @test file["target/supplement_policy"] === :none

            @test file["route/artifact_role"] ===
                  :physical_gausslet_endpoint_target
            @test file["route/source_plan_status"] ===
                  :blocked_pqs_diatomic_physical_gausslet_core_shell_source_plan
            @test file["route/final_basis_status"] !==
                  :available_pqs_complete_core_shell_final_basis
            @test file["route/h1_status"] !==
                  :materialized_pqs_complete_core_shell_final_h1_solve
            @test file["route/h1_materialized"] == false
            @test file["route/h1_j_materialized"] == false
            @test file["route/private_rhf_materialized"] == false

            @test file["basis/retained_atom_core_interiors"] == true
            @test file["basis/source_plan_role"] ===
                  :atom_contact_core_plus_pqs_shared_shells
            @test !haskey(file, "basis/final_dimension")
            @test !haskey(file, "basis/final_overlap_identity_error")

            @test file["physics/endpoint_ready"] == false
            @test file["physics/endpoint_blocker"] ===
                  :source_plan_candidate_not_route_authority
            @test !haskey(file, "physics/h1_lowest")
            @test file["comparison/ready"] == false
            @test file["private_rhf/requested"] == false
            @test file["private_rhf/materialized"] == false
        end
    end
end
