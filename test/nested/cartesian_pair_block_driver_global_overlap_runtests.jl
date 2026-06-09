# Runtime role: tiny smoke / private driver-facing global overlap hook test.
#
# This verifies the overlap-only driver-facing hook over structured
# pair-block materialization state. It does not cover Hamiltonians, Coulomb,
# IDA/MWG, PQS Lowdin/projection, exports, route-driver wiring, or
# White-Lindsey oracle fixtures.

using Test
using GaussletBases

const CPBMDriverOverlap = GaussletBases.CartesianPairBlockMaterialization
const CPBDriverOverlap = GaussletBases.CartesianCPB
const CTLDriverOverlap = GaussletBases.CartesianTerminalLowering
const CRUDriverOverlap = GaussletBases.CartesianRetainedUnits
const CRTCDriverOverlap =
    GaussletBases.CartesianRetainedUnitTransformContracts
const CUPDriverOverlap = GaussletBases.CartesianUnitPairs
const CPOPDriverOverlap = GaussletBases.CartesianPairOperatorPlans

function _driver_overlap_pair_operator_plan()
    lowering_plan = CTLDriverOverlap.TerminalLoweringPlan(
        CTLDriverOverlap.PQSLowering(q = 2),
        (),
        (),
        (; status = :available_terminal_lowering_plan, materialized = false),
        (; fixture = :driver_global_overlap_hook),
    )
    retained_plan = CRUDriverOverlap.RetainedUnitPlan(
        CRUDriverOverlap.MetadataOnlyRetainedUnits(),
        lowering_plan,
        (),
        (; status = :available_retained_unit_plan, materialized = false),
        (; fixture = :driver_global_overlap_hook),
    )
    unit_pair_plan = CUPDriverOverlap.UnitPairPlan(
        CUPDriverOverlap.MetadataOnlyUnitPairs(),
        retained_plan,
        (),
        nothing,
        (; status = :available_unit_pair_plan, materialized = false),
        (; fixture = :driver_global_overlap_hook),
    )
    transform_plan = CRTCDriverOverlap.RetainedUnitTransformContractPlan(
        CRTCDriverOverlap.MetadataOnlyRetainedUnitTransformContracts(),
        retained_plan,
        (),
        (; status = :available_transform_contract_plan, materialized = false),
        (; fixture = :driver_global_overlap_hook),
    )
    return CPOPDriverOverlap.PairOperatorPlan(
        CPOPDriverOverlap.MetadataOnlyPairOperatorPlans(),
        unit_pair_plan,
        transform_plan,
        (),
        nothing,
        (; status = :available_pair_operator_plan, materialized = false),
        (; fixture = :driver_global_overlap_hook),
    )
end

function _driver_overlap_plan()
    direct_source = CPBDriverOverlap.cpb(
        1:1,
        1:2,
        1:1;
        role = :driver_global_overlap_direct_source,
    )
    record = CPBMDriverOverlap.PairBlockMaterializationRecord(
        (:direct_diag, :direct_diag),
        1,
        :direct_direct,
        :synthetic_source_operator_path,
        (; left = :synthetic_left_transform, right = :synthetic_right_transform),
        (; left = :synthetic_left_realization, right = :synthetic_right_realization),
        :synthetic_final_block_path,
        (:overlap,),
        :direct_direct_pair_block_materialization_pilot,
        :ready_metadata_only_not_materialized,
        nothing,
        false,
        (;
            left_source_cpbs = (direct_source,),
            right_source_cpbs = (direct_source,),
            pair_index = 1,
            selector_family = :direct_direct,
            left_column_range = 1:2,
            right_column_range = 1:2,
        ),
    )
    return CPBMDriverOverlap.PairBlockMaterializationPlan(
        CPBMDriverOverlap.MetadataOnlyPairBlockMaterialization(),
        _driver_overlap_pair_operator_plan(),
        (record,),
        (; status = :available_pair_block_materialization_plan),
        (; fixture = :driver_global_overlap_hook),
    )
end

function _driver_overlap_inputs()
    return (;
        parent_axis_counts = (2, 2, 2),
        overlap_1d = (;
            x = [1.0 0.2; 0.2 1.1],
            y = [1.2 0.3; 0.3 1.3],
            z = [1.4 0.4; 0.4 1.5],
        ),
    )
end

function _driver_overlap_expected_matrix()
    overlap = _driver_overlap_inputs().overlap_1d
    scale = overlap.x[1, 1] * overlap.z[1, 1]
    return scale .* overlap.y
end

function _test_driver_overlap_nonclaim_flags(result)
    @test !result.route_driver_wiring
    @test !result.hamiltonian_data_materialized
    @test !result.global_hamiltonian_data_materialized
    @test !result.coulomb_materialized
    @test !result.ida_mwg_data_materialized
    @test !result.pqs_lowdin_materialized
    @test !result.pqs_shell_projection_materialized
    @test !result.artifacts_materialized
    @test !result.exports_materialized
    @test !result.full_white_lindsey_route_assembled
end

@testset "CartesianPairBlockMaterialization driver global overlap hook" begin
    plan = _driver_overlap_plan()
    expected = CPBMDriverOverlap.route_state_global_overlap_matrix(
        plan;
        global_dimension = 2,
        inputs = _driver_overlap_inputs(),
    )
    direct = CPBMDriverOverlap.driver_global_overlap_result(
        plan;
        global_dimension = 2,
        inputs = _driver_overlap_inputs(),
    )
    wrapped = CPBMDriverOverlap.driver_global_overlap_result(
        (; pair_block_materialization_plan = plan);
        global_dimension = 2,
        factors = _driver_overlap_inputs(),
    )
    terminal_route_state_source =
        (; terminal_route_state = (; pair_block_materialization_plan = plan))
    terminal_route_state = CPBMDriverOverlap.driver_global_overlap_result(
        terminal_route_state_source;
        global_dimension = 2,
        inputs = _driver_overlap_inputs(),
    )
    expected_matrix = _driver_overlap_expected_matrix()

    @test direct.status === :materialized_route_global_overlap_matrix
    @test direct.global_one_body_term_matrix_materialized
    @test direct.global_overlap_matrix_materialized
    @test size(terminal_route_state.global_matrix_result.matrix) == (2, 2)
    @test terminal_route_state.status === :materialized_route_global_overlap_matrix
    @test terminal_route_state.global_overlap_matrix_materialized
    @test terminal_route_state.global_one_body_term_matrix_materialized
    @test terminal_route_state.global_matrix_result.matrix ≈ expected_matrix
    @test terminal_route_state.global_matrix_result.matrix[1, 1] ≈ 1.68
    @test terminal_route_state.global_matrix_result.matrix[1, 2] ≈ 0.42
    @test terminal_route_state.global_matrix_result.matrix[2, 1] ≈ 0.42
    @test terminal_route_state.global_matrix_result.matrix[2, 2] ≈ 1.82
    @test direct.global_matrix_result.matrix ≈ expected.global_matrix_result.matrix
    @test direct.global_matrix_result.matrix ≈ expected_matrix
    @test wrapped.global_matrix_result.matrix ≈ expected.global_matrix_result.matrix
    @test terminal_route_state.global_matrix_result.matrix ≈
          expected.global_matrix_result.matrix
    @test direct.pair_block_materialization_plan === plan

    missing_plan = CPBMDriverOverlap.driver_global_overlap_result(
        (; route_stage = :missing_pair_block_plan);
        global_dimension = 2,
        inputs = _driver_overlap_inputs(),
    )
    @test missing_plan.status === :blocked_route_global_overlap_matrix
    @test missing_plan.blocker === :missing_pair_block_materialization_plan
    @test !missing_plan.global_one_body_term_matrix_materialized
    @test !missing_plan.global_overlap_matrix_materialized

    missing_dimension = CPBMDriverOverlap.driver_global_overlap_result(
        (; pair_block_materialization_plan = plan);
        inputs = _driver_overlap_inputs(),
    )
    @test missing_dimension.status === :blocked_route_global_overlap_matrix
    @test missing_dimension.blocker === :missing_global_dimension
    @test !missing_dimension.global_one_body_term_matrix_materialized

    _test_driver_overlap_nonclaim_flags(direct)
    _test_driver_overlap_nonclaim_flags(missing_plan)
end
