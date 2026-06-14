using Test
using JLD2

const _H2_PQS_DRIVER =
    normpath(joinpath(@__DIR__, "..", "..", "bin", "cartesian_ham_builder.jl"))
const _H2_PQS_INPUT =
    normpath(joinpath(@__DIR__, "..", "driver_inputs", "h2_pqs_q5_gausslet_only_r4.jl"))

function _h2_tuple_locations(value)
    return Tuple(Tuple(Float64(component) for component in location) for location in value)
end

@testset "cartesian_ham_builder_h2_pqs_q5_gausslet_only_r4 readiness artifact" begin
    mktempdir() do dir
        outfile = joinpath(dir, "h2_pqs_q5_gausslet_only_r4.jld2")
        tsvfile = joinpath(dir, "h2_pqs_q5_gausslet_only_r4.tsv")
        println("h2_readiness_artifact_path=", outfile)
        saved_args = copy(ARGS)
        empty!(ARGS)
        append!(ARGS, [_H2_PQS_INPUT, "outfile=$(repr(outfile))", "tsvfile=$(repr(tsvfile))"])
        try
            include(_H2_PQS_DRIVER)
        finally
            empty!(ARGS)
            append!(ARGS, saved_args)
        end

        @test isfile(outfile)
        jldopen(outfile, "r") do file
            @test Tuple(file["system/atom_symbols"]) == ("H", "H")
            @test Tuple(file["system/nuclear_charges"]) == (1, 1)
            @test _h2_tuple_locations(file["system/atom_locations"]) ==
                  ((0.0, 0.0, -2.0), (0.0, 0.0, 2.0))
            @test file["system/bond_axis"] === :z
            @test file["system/bond_length"] ≈ 4.0 atol = 0.0 rtol = 0.0

            @test file["config/route_family"] === :pqs_source_box
            @test file["config/route_kind"] ===
                  :bond_aligned_diatomic_fixed_q_complete_core_shell
            @test file["config/q"] == 5
            @test file["config/n_s"] == 5
            @test file["config/supplement_policy"] === :none
            @test file["config/comparison_ready"] == false
            @test file["config/run_final_basis"] == true

            @test file["comparison/ready"] == false
            @test file["comparison/blocker"] ===
                  :supplemented_reference_not_comparable_to_gausslet_only
            @test file["comparison/reference_label"] ==
                  "WL/QW H2 R=4 supplemented reference not used"
            @test !haskey(file, "comparison/wl_rhf_total")
            @test !haskey(file, "comparison/delta_rhf")
            @test file["route/artifact_role"] === :source_box_diagnostic

            @test file["parent/parent_axis_counts"] == (x = 9, y = 9, z = 15)
            @test file["parent/parent_axis_counts_source"] ===
                  :constructed_parent_axis_probe
            @test file["parent/parent_materialization_blocker"] === nothing
            @test file["parent/parent_basis_object_available"] == true
            @test file["parent/parent_qw_basis_object_available"] == true
            @test file["parent/parent_axis_bundle_object_available"] == true
            @test file["parent/parent_basis_object_type_label"] ==
                  "CartesianParentGaussletBasis3D"
            @test file["parent/parent_qw_basis_object_type_label"] ==
                  "BondAlignedDiatomicQWBasis3D"
            @test file["parent/parent_axis_bundle_object_type_label"] ==
                  "_CartesianNestedAxisBundles3D"

            @test file["private_rhf/requested"] == false
            @test file["private_rhf/materialized"] == false
            @test file["route/h1_j_materialized"] == false
            @test file["route/h1_materialized"] == true
            @test file["route/private_rhf_materialized"] == false
            @test file["route/public_api"] == false
            @test file["route/exports_materialized"] == false
            @test file["route/artifacts_materialized"] == false
            @test file["route/readiness_status"] ===
                  :blocked_diatomic_complete_core_shell_ham_readiness
            @test file["route/readiness_blocker"] ===
                  :missing_diatomic_complete_core_shell_hamiltonian_handoff_payload
            @test file["route/source_plan_status"] ===
                  :available_pqs_diatomic_complete_core_shell_source_plan
            @test file["route/final_basis_status"] ===
                  :available_pqs_complete_core_shell_final_basis
            @test file["route/h1_status"] ===
                  :materialized_pqs_complete_core_shell_final_h1_solve
            @test file["basis/final_dimension"] == 221
            @test file["basis/retained_atom_core_interiors"] == false
            @test file["basis/source_plan_role"] ===
                  :boundary_source_box_diagnostic
            @test isfinite(file["basis/final_overlap_identity_error"])
            @test file["basis/final_overlap_identity_error"] < 1e-10
            @test file["physics/endpoint_ready"] == false
            @test file["physics/endpoint_blocker"] ===
                  :retained_atom_core_interiors_missing
            @test isfinite(file["physics/h1_lowest"])
            @test file["physics/h1_hamiltonian_matrix_finite"] == true
            @test isfinite(file["physics/h1_hamiltonian_symmetry_error"])
            @test file["physics/h1_hamiltonian_symmetry_error"] < 1e-10
        end
    end
end
