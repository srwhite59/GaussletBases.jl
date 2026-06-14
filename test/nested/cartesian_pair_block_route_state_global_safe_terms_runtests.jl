# Runtime role: private route-state global safe-term adapter smoke.
# This verifies overlap, kinetic, position, and x2 route-state adapters over structured
# pair-block materialization state. It does not cover Hamiltonians, Coulomb,
# IDA/MWG, PQS Lowdin/projection, exports, route-driver wiring, or
# White-Lindsey oracle fixtures.

using Test
using GaussletBases

const CPBMRouteStateSafeTerms = GaussletBases.CartesianPairBlockMaterialization
const CPBRouteStateSafeTerms = GaussletBases.CartesianCPB
const CTLRouteStateSafeTerms = GaussletBases.CartesianTerminalLowering
const CRURouteStateSafeTerms = GaussletBases.CartesianRetainedUnits
const CRTCRouteStateSafeTerms =
    GaussletBases.CartesianRetainedUnitTransformContracts
const CUPRouteStateSafeTerms = GaussletBases.CartesianUnitPairs
const CPOPRouteStateSafeTerms = GaussletBases.CartesianPairOperatorPlans

function _route_state_safe_terms_pair_operator_plan()
    lowering_plan = CTLRouteStateSafeTerms.TerminalLoweringPlan(
        CTLRouteStateSafeTerms.PQSLowering(q = 2),
        (),
        (),
        (; status = :available_terminal_lowering_plan, materialized = false),
        (; fixture = :route_state_global_safe_terms_adapter),
    )
    retained_plan = CRURouteStateSafeTerms.RetainedUnitPlan(
        CRURouteStateSafeTerms.MetadataOnlyRetainedUnits(),
        lowering_plan,
        (),
        (; status = :available_retained_unit_plan, materialized = false),
        (; fixture = :route_state_global_safe_terms_adapter),
    )
    unit_pair_plan = CUPRouteStateSafeTerms.UnitPairPlan(
        CUPRouteStateSafeTerms.MetadataOnlyUnitPairs(),
        retained_plan,
        (),
        nothing,
        (; status = :available_unit_pair_plan, materialized = false),
        (; fixture = :route_state_global_safe_terms_adapter),
    )
    transform_plan = CRTCRouteStateSafeTerms.RetainedUnitTransformContractPlan(
        CRTCRouteStateSafeTerms.MetadataOnlyRetainedUnitTransformContracts(),
        retained_plan,
        (),
        (; status = :available_transform_contract_plan, materialized = false),
        (; fixture = :route_state_global_safe_terms_adapter),
    )
    return CPOPRouteStateSafeTerms.PairOperatorPlan(
        CPOPRouteStateSafeTerms.MetadataOnlyPairOperatorPlans(),
        unit_pair_plan,
        transform_plan,
        (),
        nothing,
        (; status = :available_pair_operator_plan, materialized = false),
        (; fixture = :route_state_global_safe_terms_adapter),
    )
end

function _route_state_safe_terms_record(; metadata = (;))
    return CPBMRouteStateSafeTerms.PairBlockMaterializationRecord(
        (:direct_diag, :direct_diag),
        1,
        :direct_direct,
        :synthetic_source_operator_path,
        (; left = :synthetic_left_transform, right = :synthetic_right_transform),
        (; left = :synthetic_left_realization, right = :synthetic_right_realization),
        :synthetic_final_block_path,
        (:overlap, :kinetic, :position_x, :position_y, :position_z,
            :x2_x, :x2_y, :x2_z),
        :direct_direct_pair_block_materialization_pilot,
        :ready_metadata_only_not_materialized,
        nothing,
        false,
        NamedTuple(metadata),
    )
end

function _route_state_safe_terms_plan()
    direct_source = CPBRouteStateSafeTerms.cpb(
        1:1,
        1:2,
        1:1;
        role = :route_state_global_safe_terms_direct_source,
    )
    record = _route_state_safe_terms_record(
        metadata = (;
            left_source_cpbs = (direct_source,),
            right_source_cpbs = (direct_source,),
            pair_index = 1,
            selector_family = :direct_direct,
            left_column_range = 1:2,
            right_column_range = 1:2,
        ),
    )
    return CPBMRouteStateSafeTerms.PairBlockMaterializationPlan(
        CPBMRouteStateSafeTerms.MetadataOnlyPairBlockMaterialization(),
        _route_state_safe_terms_pair_operator_plan(),
        (record,),
        (; status = :available_pair_block_materialization_plan),
        (; fixture = :route_state_global_safe_terms_adapter),
    )
end

function _route_state_safe_terms_inputs()
    return (;
        parent_axis_counts = (2, 2, 2),
        overlap_1d = (;
            x = [1.0 0.2; 0.2 1.1],
            y = [1.2 0.3; 0.3 1.3],
            z = [1.4 0.4; 0.4 1.5],
        ),
        kinetic_1d = (;
            x = [2.0 0.5; 0.5 2.2],
            y = [3.0 0.6; 0.6 3.3],
            z = [4.0 0.7; 0.7 4.4],
        ),
        position_1d = (;
            x = [5.0 0.7; 0.8 6.0],
            y = [7.0 0.9; 1.0 8.0],
            z = [9.0 1.1; 1.2 10.0],
        ),
        x2_1d = (;
            x = [11.0 1.3; 1.4 12.0],
            y = [13.0 1.5; 1.6 14.0],
            z = [15.0 1.7; 1.8 16.0],
        ),
    )
end

function _route_state_matrix_function(prefix::Symbol, term::Symbol)
    return getproperty(
        CPBMRouteStateSafeTerms,
        Symbol(String(prefix), "_", String(term), "_matrix"),
    )
end

function _test_route_state_nonclaim_flags(result)
    for field in (
        :route_driver_wiring,
        :hamiltonian_data_materialized,
        :global_hamiltonian_data_materialized,
        :coulomb_materialized,
        :ida_mwg_data_materialized,
        :pqs_lowdin_materialized,
        :pqs_shell_projection_materialized,
        :artifacts_materialized,
        :exports_materialized,
    )
        @test !getproperty(result, field)
    end
end

function _test_route_state_scalar_adapter(
    expected_builder,
    adapter_builder;
    materialized_status,
    blocked_status,
    materialized_flag,
)
    plan = _route_state_safe_terms_plan()
    inputs = _route_state_safe_terms_inputs()
    expected = expected_builder(plan; global_dimension = 2, inputs)
    direct = adapter_builder(plan; global_dimension = 2, inputs)

    @test direct.status === materialized_status
    @test getproperty(direct, materialized_flag)
    @test direct.global_matrix_result.matrix ≈ expected.global_matrix_result.matrix
    @test direct.pair_block_materialization_plan === plan

    missing_plan = adapter_builder(
        (; route_state = :missing_pair_block_plan);
        global_dimension = 2,
        inputs,
    )
    @test missing_plan.status === blocked_status
    @test missing_plan.blocker === :missing_pair_block_materialization_plan
    @test !missing_plan.global_one_body_term_matrix_materialized
    @test !getproperty(missing_plan, materialized_flag)

    missing_dimension = adapter_builder(
        (; pair_block_materialization_plan = plan);
        inputs,
    )
    @test missing_dimension.status === blocked_status
    @test missing_dimension.blocker === :missing_global_dimension
    @test !missing_dimension.global_one_body_term_matrix_materialized
    _test_route_state_nonclaim_flags(direct)
end

function _test_route_state_axis_family(
    axis_adapter,
    axis_expected,
    terms,
)
    plan = _route_state_safe_terms_plan()
    inputs = _route_state_safe_terms_inputs()
    axis_result = axis_adapter(plan; axis = :x, global_dimension = 2, inputs)
    expected_axis = axis_expected(plan; global_dimension = 2, inputs)
    @test axis_result.global_matrix_result.matrix ≈
          expected_axis.global_matrix_result.matrix

    for term in terms
        expected =
            _route_state_matrix_function(:route_global, term)(
                plan;
                global_dimension = 2,
                inputs,
            )
        direct =
            _route_state_matrix_function(:route_state_global, term)(
                plan;
                global_dimension = 2,
                inputs,
            )
        @test direct.status ===
              Symbol("materialized_route_global_", String(term), "_matrix")
        @test getproperty(
            direct,
            Symbol("global_", String(term), "_matrix_materialized"),
        )
        @test direct.global_matrix_result.matrix ≈
              expected.global_matrix_result.matrix
        @test direct.pair_block_materialization_plan === plan
        _test_route_state_nonclaim_flags(direct)
    end

    wrapped = _route_state_matrix_function(:route_state_global, terms[2])(
        (; pair_block_materialization_plan = plan);
        global_dimension = 2,
        factors = inputs,
    )
    wrapped_expected =
        _route_state_matrix_function(:route_global, terms[2])(
            plan;
            global_dimension = 2,
            inputs,
        )
    @test wrapped.global_matrix_result.matrix ≈
          wrapped_expected.global_matrix_result.matrix
end

@testset "CartesianPairBlockMaterialization route-state global overlap adapter" begin
    _test_route_state_scalar_adapter(
        CPBMRouteStateSafeTerms.route_global_overlap_matrix,
        CPBMRouteStateSafeTerms.route_state_global_overlap_matrix;
        materialized_status = :materialized_route_global_overlap_matrix,
        blocked_status = :blocked_route_global_overlap_matrix,
        materialized_flag = :global_overlap_matrix_materialized,
    )
end

@testset "CartesianPairBlockMaterialization route-state global kinetic adapter" begin
    _test_route_state_scalar_adapter(
        CPBMRouteStateSafeTerms.route_global_kinetic_matrix,
        CPBMRouteStateSafeTerms.route_state_global_kinetic_matrix;
        materialized_status = :materialized_route_global_kinetic_matrix,
        blocked_status = :blocked_route_global_kinetic_matrix,
        materialized_flag = :global_kinetic_matrix_materialized,
    )
end

@testset "CartesianPairBlockMaterialization route-state global position adapter" begin
    _test_route_state_axis_family(
        CPBMRouteStateSafeTerms.route_state_global_position_matrix,
        CPBMRouteStateSafeTerms.route_global_position_x_matrix,
        (:position_x, :position_y, :position_z),
    )
end

@testset "CartesianPairBlockMaterialization route-state global x2 adapter" begin
    _test_route_state_axis_family(
        CPBMRouteStateSafeTerms.route_state_global_x2_matrix,
        CPBMRouteStateSafeTerms.route_global_x2_x_matrix,
        (:x2_x, :x2_y, :x2_z),
    )
end
