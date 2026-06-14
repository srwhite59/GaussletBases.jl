# Default nested runner: fast/default nested contract tests only.
# Slow/manual integration and audit tests live in integration_runtests.jl.
#
# For mixed one-body consumer edits, prefer the tiny smoke test directly:
#   julia --project=. test/nested/cartesian_pair_block_one_body_consumer_smoke_runtests.jl
# Heavier mixed one-body contract/boundary files are for semantic changes or
# baton closeout, not routine small-edit validation. White-Lindsey boundary and
# oracle comparisons are gate tests only when LW adapter behavior changes or at
# baton closeout.

using Test
using LinearAlgebra
using SparseArrays
using GaussletBases

include("pqs_standard_source_box_route_setup_runtests.jl")
include("pqs_standard_parent_axis_readiness_runtests.jl")
include("pqs_source_box_route_driver_crc_print_line_runtests.jl")
include("cartesian_route_core_examples_runtests.jl")
include("cartesian_selected_terminal_lowering_contract_inventory_runtests.jl")
include("cartesian_retained_units_contract_runtests.jl")
include("cartesian_retained_unit_transform_contracts_runtests.jl")
include("cartesian_terminal_route_retained_units_fingerprint_runtests.jl")
include("cartesian_unit_pairs_contract_runtests.jl")
include("cartesian_pair_operator_plans_contract_runtests.jl")
include("cartesian_final_basis_realization_contract_runtests.jl")
include("cartesian_pair_block_one_body_consumer_smoke_runtests.jl")
include("cartesian_pair_block_materialization_contract_runtests.jl")
include("cartesian_terminal_route_flat_glue_cleanup_runtests.jl")
include("cartesian_pair_stage_fingerprint_helpers_runtests.jl")

@testset "Cartesian nested owned-unit coverage audit" begin
    dense_unit = GaussletBases._CartesianNestedOwnedUnit3D(
        :endcap_a,
        [1, 2],
        [1.0 0.0; 0.0 1.0];
        metadata = (side = :left,),
    )
    sparse_unit = GaussletBases._CartesianNestedOwnedUnit3D(
        :panel_b,
        [3, 4],
        sparse([1, 2], [1, 1], [0.5, 0.5], 2, 1);
        metadata = (side = :right,),
    )
    exact = GaussletBases._nested_owned_unit_coverage_audit(
        [dense_unit, sparse_unit],
        [1, 2, 3, 4],
    )

    @test dense_unit.coefficient_matrix isa Matrix{Float64}
    @test sparse_unit.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test dense_unit.role == :endcap_a
    @test dense_unit.metadata.side == :left
    @test exact.expected_support_count == 4
    @test exact.owned_support_count == 4
    @test exact.duplicate_count == 0
    @test exact.missing_count == 0
    @test exact.outside_count == 0
    @test exact.retained_count == 3
    @test exact.coverage_ok

    duplicate_unit = GaussletBases._CartesianNestedOwnedUnit3D(:duplicate_panel, [2, 3], ones(2, 1))
    duplicate = GaussletBases._nested_owned_unit_coverage_audit(
        [dense_unit, duplicate_unit],
        [1, 2, 3],
    )
    @test duplicate.duplicate_count == 1
    @test duplicate.missing_count == 0
    @test duplicate.outside_count == 0
    @test !duplicate.coverage_ok

    missing = GaussletBases._nested_owned_unit_coverage_audit([dense_unit], [1, 2, 3])
    @test missing.missing_count == 1
    @test !missing.coverage_ok

    outside = GaussletBases._nested_owned_unit_coverage_audit([dense_unit], [1])
    @test outside.outside_count == 1
    @test !outside.coverage_ok

    zero_retained = GaussletBases._CartesianNestedOwnedUnit3D(:empty_panel, [1, 2], zeros(2, 0))
    nonfinite = GaussletBases._CartesianNestedOwnedUnit3D(:bad_panel, [1], [Inf;;])
    @test_throws DimensionMismatch GaussletBases._CartesianNestedOwnedUnit3D(:bad_rows, [1, 2], ones(1, 1))
    @test_throws ArgumentError GaussletBases._nested_owned_unit_coverage_audit([zero_retained], [1, 2])
    @test_throws ArgumentError GaussletBases._nested_owned_unit_coverage_audit([nonfinite], [1])
    @test_throws ArgumentError GaussletBases._nested_owned_unit_coverage_audit([dense_unit], [1, 1, 2])
end

include("bond_aligned_diatomic_high_order_recipe_policy_metadata_runtests.jl")
include("global_timing_macro_surface_runtests.jl")

include("cartesian_parent_gausslet_basis_identity_runtests.jl")

include("cartesian_contracted_parent_scaffold_runtests.jl")

include("qw_residual_space_keep_policy_runtests.jl")
