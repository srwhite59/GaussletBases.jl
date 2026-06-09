# Runtime role: tiny smoke / private route-state global safe-term adapter test.
#
# This verifies overlap, kinetic, position, and x2 route-state adapters over structured
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
        (
            :overlap,
            :kinetic,
            :position_x,
            :position_y,
            :position_z,
            :x2_x,
            :x2_y,
            :x2_z,
        ),
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

function _route_state_position_global_matrix(term::Symbol, source; kwargs...)
    term === :position_x &&
        return CPBMRouteStateOverlap.route_global_position_x_matrix(
            source;
            kwargs...,
        )
    term === :position_y &&
        return CPBMRouteStateOverlap.route_global_position_y_matrix(
            source;
            kwargs...,
        )
    term === :position_z &&
        return CPBMRouteStateOverlap.route_global_position_z_matrix(
            source;
            kwargs...,
        )
    throw(ArgumentError("expected route-state position term"))
end

function _route_state_position_adapter(term::Symbol, source; kwargs...)
    term === :position_x &&
        return CPBMRouteStateOverlap.route_state_global_position_x_matrix(
            source;
            kwargs...,
        )
    term === :position_y &&
        return CPBMRouteStateOverlap.route_state_global_position_y_matrix(
            source;
            kwargs...,
        )
    term === :position_z &&
        return CPBMRouteStateOverlap.route_state_global_position_z_matrix(
            source;
            kwargs...,
        )
    throw(ArgumentError("expected route-state position term"))
end

function _route_state_x2_global_matrix(term::Symbol, source; kwargs...)
    term === :x2_x &&
        return CPBMRouteStateOverlap.route_global_x2_x_matrix(source; kwargs...)
    term === :x2_y &&
        return CPBMRouteStateOverlap.route_global_x2_y_matrix(source; kwargs...)
    term === :x2_z &&
        return CPBMRouteStateOverlap.route_global_x2_z_matrix(source; kwargs...)
    throw(ArgumentError("expected route-state x2 term"))
end

function _route_state_x2_adapter(term::Symbol, source; kwargs...)
    term === :x2_x &&
        return CPBMRouteStateOverlap.route_state_global_x2_x_matrix(
            source;
            kwargs...,
        )
    term === :x2_y &&
        return CPBMRouteStateOverlap.route_state_global_x2_y_matrix(
            source;
            kwargs...,
        )
    term === :x2_z &&
        return CPBMRouteStateOverlap.route_state_global_x2_z_matrix(
            source;
            kwargs...,
        )
    throw(ArgumentError("expected route-state x2 term"))
end

function _test_route_state_nonclaim_flags(result)
    @test !result.route_driver_wiring
    @test !result.hamiltonian_data_materialized
    @test !result.global_hamiltonian_data_materialized
    @test !result.coulomb_materialized
    @test !result.ida_mwg_data_materialized
    @test !result.pqs_lowdin_materialized
    @test !result.pqs_shell_projection_materialized
    @test !result.artifacts_materialized
    @test !result.exports_materialized
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

    _test_route_state_nonclaim_flags(direct)
end

@testset "CartesianPairBlockMaterialization route-state global kinetic adapter" begin
    plan = _route_state_overlap_plan()
    expected = CPBMRouteStateOverlap.route_global_kinetic_matrix(
        plan;
        global_dimension = 2,
        inputs = _route_state_overlap_inputs(),
    )
    direct = CPBMRouteStateOverlap.route_state_global_kinetic_matrix(
        plan;
        global_dimension = 2,
        inputs = _route_state_overlap_inputs(),
    )
    wrapped = CPBMRouteStateOverlap.route_state_global_kinetic_matrix(
        (; pair_block_materialization_plan = plan);
        global_dimension = 2,
        factors = _route_state_overlap_inputs(),
    )
    route_state = CPBMRouteStateOverlap.route_state_global_kinetic_matrix(
        (; terminal_route_state = (; pair_block_materialization_plan = plan));
        global_dimension = 2,
        inputs = _route_state_overlap_inputs(),
    )

    @test direct.status === :materialized_route_global_kinetic_matrix
    @test direct.global_kinetic_matrix_materialized
    @test direct.global_matrix_result.matrix ≈ expected.global_matrix_result.matrix
    @test wrapped.global_matrix_result.matrix ≈ expected.global_matrix_result.matrix
    @test route_state.global_matrix_result.matrix ≈
          expected.global_matrix_result.matrix
    @test direct.pair_block_materialization_plan === plan

    missing_plan = CPBMRouteStateOverlap.route_state_global_kinetic_matrix(
        (; route_state = :missing_pair_block_plan);
        global_dimension = 2,
        inputs = _route_state_overlap_inputs(),
    )
    @test missing_plan.status === :blocked_route_global_kinetic_matrix
    @test missing_plan.blocker === :missing_pair_block_materialization_plan
    @test !missing_plan.global_one_body_term_matrix_materialized
    @test !missing_plan.global_kinetic_matrix_materialized

    missing_dimension =
        CPBMRouteStateOverlap.route_state_global_kinetic_matrix(
            (; pair_block_materialization_plan = plan);
            inputs = _route_state_overlap_inputs(),
        )
    @test missing_dimension.status === :blocked_route_global_kinetic_matrix
    @test missing_dimension.blocker === :missing_global_dimension
    @test !missing_dimension.global_one_body_term_matrix_materialized

    _test_route_state_nonclaim_flags(direct)
end

@testset "CartesianPairBlockMaterialization route-state global position adapter" begin
    plan = _route_state_overlap_plan()

    axis_adapter = CPBMRouteStateOverlap.route_state_global_position_matrix(
        plan;
        axis = :x,
        global_dimension = 2,
        inputs = _route_state_overlap_inputs(),
    )
    axis_expected = CPBMRouteStateOverlap.route_global_position_x_matrix(
        plan;
        global_dimension = 2,
        inputs = _route_state_overlap_inputs(),
    )
    @test axis_adapter.global_matrix_result.matrix ≈
          axis_expected.global_matrix_result.matrix

    for term in (:position_x, :position_y, :position_z)
        expected = _route_state_position_global_matrix(
            term,
            plan;
            global_dimension = 2,
            inputs = _route_state_overlap_inputs(),
        )
        direct = _route_state_position_adapter(
            term,
            plan;
            global_dimension = 2,
            inputs = _route_state_overlap_inputs(),
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

    wrapped = CPBMRouteStateOverlap.route_state_global_position_y_matrix(
        (; pair_block_materialization_plan = plan);
        global_dimension = 2,
        factors = _route_state_overlap_inputs(),
    )
    wrapped_expected = CPBMRouteStateOverlap.route_global_position_y_matrix(
        plan;
        global_dimension = 2,
        inputs = _route_state_overlap_inputs(),
    )
    @test wrapped.global_matrix_result.matrix ≈
          wrapped_expected.global_matrix_result.matrix

    route_state = CPBMRouteStateOverlap.route_state_global_position_z_matrix(
        (; terminal_route_state = (; pair_block_materialization_plan = plan));
        global_dimension = 2,
        inputs = _route_state_overlap_inputs(),
    )
    route_state_expected =
        CPBMRouteStateOverlap.route_global_position_z_matrix(
            plan;
            global_dimension = 2,
            inputs = _route_state_overlap_inputs(),
        )
    @test route_state.global_matrix_result.matrix ≈
          route_state_expected.global_matrix_result.matrix

    missing_plan = CPBMRouteStateOverlap.route_state_global_position_x_matrix(
        (; route_state = :missing_pair_block_plan);
        global_dimension = 2,
        inputs = _route_state_overlap_inputs(),
    )
    @test missing_plan.status === :blocked_route_global_position_x_matrix
    @test missing_plan.blocker === :missing_pair_block_materialization_plan

    missing_dimension =
        CPBMRouteStateOverlap.route_state_global_position_x_matrix(
            (; pair_block_materialization_plan = plan);
            inputs = _route_state_overlap_inputs(),
        )
    @test missing_dimension.status === :blocked_route_global_position_x_matrix
    @test missing_dimension.blocker === :missing_global_dimension
end

@testset "CartesianPairBlockMaterialization route-state global x2 adapter" begin
    plan = _route_state_overlap_plan()

    axis_adapter = CPBMRouteStateOverlap.route_state_global_x2_matrix(
        plan;
        axis = :x,
        global_dimension = 2,
        inputs = _route_state_overlap_inputs(),
    )
    axis_expected = CPBMRouteStateOverlap.route_global_x2_x_matrix(
        plan;
        global_dimension = 2,
        inputs = _route_state_overlap_inputs(),
    )
    @test axis_adapter.global_matrix_result.matrix ≈
          axis_expected.global_matrix_result.matrix

    for term in (:x2_x, :x2_y, :x2_z)
        expected = _route_state_x2_global_matrix(
            term,
            plan;
            global_dimension = 2,
            inputs = _route_state_overlap_inputs(),
        )
        direct = _route_state_x2_adapter(
            term,
            plan;
            global_dimension = 2,
            inputs = _route_state_overlap_inputs(),
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

    wrapped = CPBMRouteStateOverlap.route_state_global_x2_y_matrix(
        (; pair_block_materialization_plan = plan);
        global_dimension = 2,
        factors = _route_state_overlap_inputs(),
    )
    wrapped_expected = CPBMRouteStateOverlap.route_global_x2_y_matrix(
        plan;
        global_dimension = 2,
        inputs = _route_state_overlap_inputs(),
    )
    @test wrapped.global_matrix_result.matrix ≈
          wrapped_expected.global_matrix_result.matrix

    route_state = CPBMRouteStateOverlap.route_state_global_x2_z_matrix(
        (; terminal_route_state = (; pair_block_materialization_plan = plan));
        global_dimension = 2,
        inputs = _route_state_overlap_inputs(),
    )
    route_state_expected = CPBMRouteStateOverlap.route_global_x2_z_matrix(
        plan;
        global_dimension = 2,
        inputs = _route_state_overlap_inputs(),
    )
    @test route_state.global_matrix_result.matrix ≈
          route_state_expected.global_matrix_result.matrix

    missing_plan = CPBMRouteStateOverlap.route_state_global_x2_x_matrix(
        (; route_state = :missing_pair_block_plan);
        global_dimension = 2,
        inputs = _route_state_overlap_inputs(),
    )
    @test missing_plan.status === :blocked_route_global_x2_x_matrix
    @test missing_plan.blocker === :missing_pair_block_materialization_plan

    missing_dimension = CPBMRouteStateOverlap.route_state_global_x2_x_matrix(
        (; pair_block_materialization_plan = plan);
        inputs = _route_state_overlap_inputs(),
    )
    @test missing_dimension.status === :blocked_route_global_x2_x_matrix
    @test missing_dimension.blocker === :missing_global_dimension
end
