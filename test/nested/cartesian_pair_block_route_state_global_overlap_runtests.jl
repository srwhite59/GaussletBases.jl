# Runtime role: tiny smoke / private route-state global overlap adapter test.
#
# This verifies the overlap-only route-state adapter over structured
# pair-block materialization state. It does not cover Hamiltonians, Coulomb,
# IDA/MWG, PQS Lowdin/projection, exports, route-driver wiring, or
# White-Lindsey oracle fixtures.

using Test
using GaussletBases

const CPBMRouteStateOverlap = GaussletBases.CartesianPairBlockMaterialization
const CPBRouteStateOverlap = GaussletBases.CartesianCPB
const CTLRouteStateOverlap = GaussletBases.CartesianTerminalLowering
const CRURouteStateOverlap = GaussletBases.CartesianRetainedUnits
const CRTCRouteStateOverlap =
    GaussletBases.CartesianRetainedUnitTransformContracts
const CUPRouteStateOverlap = GaussletBases.CartesianUnitPairs
const CPOPRouteStateOverlap = GaussletBases.CartesianPairOperatorPlans

function _route_state_overlap_pair_operator_plan()
    lowering_plan = CTLRouteStateOverlap.TerminalLoweringPlan(
        CTLRouteStateOverlap.PQSLowering(q = 2),
        (),
        (),
        (; status = :available_terminal_lowering_plan, materialized = false),
        (; fixture = :route_state_global_overlap_adapter),
    )
    retained_plan = CRURouteStateOverlap.RetainedUnitPlan(
        CRURouteStateOverlap.MetadataOnlyRetainedUnits(),
        lowering_plan,
        (),
        (; status = :available_retained_unit_plan, materialized = false),
        (; fixture = :route_state_global_overlap_adapter),
    )
    unit_pair_plan = CUPRouteStateOverlap.UnitPairPlan(
        CUPRouteStateOverlap.MetadataOnlyUnitPairs(),
        retained_plan,
        (),
        nothing,
        (; status = :available_unit_pair_plan, materialized = false),
        (; fixture = :route_state_global_overlap_adapter),
    )
    transform_plan = CRTCRouteStateOverlap.RetainedUnitTransformContractPlan(
        CRTCRouteStateOverlap.MetadataOnlyRetainedUnitTransformContracts(),
        retained_plan,
        (),
        (; status = :available_transform_contract_plan, materialized = false),
        (; fixture = :route_state_global_overlap_adapter),
    )
    return CPOPRouteStateOverlap.PairOperatorPlan(
        CPOPRouteStateOverlap.MetadataOnlyPairOperatorPlans(),
        unit_pair_plan,
        transform_plan,
        (),
        nothing,
        (; status = :available_pair_operator_plan, materialized = false),
        (; fixture = :route_state_global_overlap_adapter),
    )
end

function _route_state_overlap_record(; metadata = (;))
    return CPBMRouteStateOverlap.PairBlockMaterializationRecord(
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
        NamedTuple(metadata),
    )
end

function _route_state_overlap_plan()
    direct_source = CPBRouteStateOverlap.cpb(
        1:1,
        1:2,
        1:1;
        role = :route_state_global_overlap_direct_source,
    )
    record = _route_state_overlap_record(
        metadata = (;
            left_source_cpbs = (direct_source,),
            right_source_cpbs = (direct_source,),
            pair_index = 1,
            selector_family = :direct_direct,
            left_column_range = 1:2,
            right_column_range = 1:2,
        ),
    )
    return CPBMRouteStateOverlap.PairBlockMaterializationPlan(
        CPBMRouteStateOverlap.MetadataOnlyPairBlockMaterialization(),
        _route_state_overlap_pair_operator_plan(),
        (record,),
        (; status = :available_pair_block_materialization_plan),
        (; fixture = :route_state_global_overlap_adapter),
    )
end

function _route_state_overlap_inputs()
    return (;
        parent_axis_counts = (2, 2, 2),
        overlap_1d = (;
            x = [1.0 0.2; 0.2 1.1],
            y = [1.2 0.3; 0.3 1.3],
            z = [1.4 0.4; 0.4 1.5],
        ),
    )
end

@testset "CartesianPairBlockMaterialization route-state global overlap adapter" begin
    plan = _route_state_overlap_plan()
    expected = CPBMRouteStateOverlap.route_global_overlap_matrix(
        plan;
        global_dimension = 2,
        inputs = _route_state_overlap_inputs(),
    )
    direct = CPBMRouteStateOverlap.route_state_global_overlap_matrix(
        plan;
        global_dimension = 2,
        inputs = _route_state_overlap_inputs(),
    )
    wrapped = CPBMRouteStateOverlap.route_state_global_overlap_matrix(
        (; pair_block_materialization_plan = plan);
        global_dimension = 2,
        factors = _route_state_overlap_inputs(),
    )
    route_state = CPBMRouteStateOverlap.route_state_global_overlap_matrix(
        (; terminal_route_state = (; pair_block_materialization_plan = plan));
        global_dimension = 2,
        inputs = _route_state_overlap_inputs(),
    )

    @test direct.status === :materialized_route_global_overlap_matrix
    @test direct.global_overlap_matrix_materialized
    @test direct.global_matrix_result.matrix ≈ expected.global_matrix_result.matrix
    @test wrapped.global_matrix_result.matrix ≈ expected.global_matrix_result.matrix
    @test route_state.global_matrix_result.matrix ≈
          expected.global_matrix_result.matrix
    @test direct.pair_block_materialization_plan === plan

    missing_plan = CPBMRouteStateOverlap.route_state_global_overlap_matrix(
        (; route_state = :missing_pair_block_plan);
        global_dimension = 2,
        inputs = _route_state_overlap_inputs(),
    )
    @test missing_plan.status === :blocked_route_global_overlap_matrix
    @test missing_plan.blocker === :missing_pair_block_materialization_plan
    @test !missing_plan.global_one_body_term_matrix_materialized
    @test !missing_plan.global_overlap_matrix_materialized

    missing_dimension =
        CPBMRouteStateOverlap.route_state_global_overlap_matrix(
            (; pair_block_materialization_plan = plan);
            inputs = _route_state_overlap_inputs(),
        )
    @test missing_dimension.status === :blocked_route_global_overlap_matrix
    @test missing_dimension.blocker === :missing_global_dimension
    @test !missing_dimension.global_one_body_term_matrix_materialized

    @test !direct.route_driver_wiring
    @test !direct.hamiltonian_data_materialized
    @test !direct.global_hamiltonian_data_materialized
    @test !direct.coulomb_materialized
    @test !direct.ida_mwg_data_materialized
    @test !direct.pqs_lowdin_materialized
    @test !direct.pqs_shell_projection_materialized
    @test !direct.artifacts_materialized
    @test !direct.exports_materialized
end
