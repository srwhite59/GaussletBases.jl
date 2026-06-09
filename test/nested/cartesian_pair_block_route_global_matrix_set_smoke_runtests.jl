# Runtime role: tiny smoke / routine route-global matrix-set adapter test.
#
# This keeps route-global safe one-body matrix-set coverage separate from the
# broader individual-term adapter contract. It does not cover Hamiltonians,
# Coulomb materialization, IDA/MWG, PQS Lowdin/projection, exports,
# route-driver wiring, or White-Lindsey oracle fixtures.

using Test
using GaussletBases

const CPBMRouteGlobalSet = GaussletBases.CartesianPairBlockMaterialization
const CPBRouteGlobalSet = GaussletBases.CartesianCPB
const CTLRouteGlobalSet = GaussletBases.CartesianTerminalLowering
const CRURouteGlobalSet = GaussletBases.CartesianRetainedUnits
const CRTCRouteGlobalSet =
    GaussletBases.CartesianRetainedUnitTransformContracts
const CUPRouteGlobalSet = GaussletBases.CartesianUnitPairs
const CPOPRouteGlobalSet = GaussletBases.CartesianPairOperatorPlans

function _route_global_set_pair_operator_plan()
    lowering_plan = CTLRouteGlobalSet.TerminalLoweringPlan(
        CTLRouteGlobalSet.PQSLowering(q = 2),
        (),
        (),
        (; status = :available_terminal_lowering_plan, materialized = false),
        (; fixture = :route_global_matrix_set_smoke),
    )
    retained_plan = CRURouteGlobalSet.RetainedUnitPlan(
        CRURouteGlobalSet.MetadataOnlyRetainedUnits(),
        lowering_plan,
        (),
        (; status = :available_retained_unit_plan, materialized = false),
        (; fixture = :route_global_matrix_set_smoke),
    )
    unit_pair_plan = CUPRouteGlobalSet.UnitPairPlan(
        CUPRouteGlobalSet.MetadataOnlyUnitPairs(),
        retained_plan,
        (),
        nothing,
        (; status = :available_unit_pair_plan, materialized = false),
        (; fixture = :route_global_matrix_set_smoke),
    )
    transform_plan = CRTCRouteGlobalSet.RetainedUnitTransformContractPlan(
        CRTCRouteGlobalSet.MetadataOnlyRetainedUnitTransformContracts(),
        retained_plan,
        (),
        (; status = :available_transform_contract_plan, materialized = false),
        (; fixture = :route_global_matrix_set_smoke),
    )
    return CPOPRouteGlobalSet.PairOperatorPlan(
        CPOPRouteGlobalSet.MetadataOnlyPairOperatorPlans(),
        unit_pair_plan,
        transform_plan,
        (),
        nothing,
        (; status = :available_pair_operator_plan, materialized = false),
        (; fixture = :route_global_matrix_set_smoke),
    )
end

function _route_global_set_record(
    pair_key::Tuple{Symbol,Symbol},
    pair_index::Int,
    pair_family::Symbol,
    materialization_path::Symbol;
    supported_terms = (
        :overlap,
        :position_x,
    ),
    metadata = (;),
)
    return CPBMRouteGlobalSet.PairBlockMaterializationRecord(
        pair_key,
        pair_index,
        pair_family,
        :synthetic_source_operator_path,
        (; left = :synthetic_left_transform, right = :synthetic_right_transform),
        (; left = :synthetic_left_realization, right = :synthetic_right_realization),
        :synthetic_final_block_path,
        supported_terms,
        materialization_path,
        :ready_metadata_only_not_materialized,
        nothing,
        false,
        NamedTuple(metadata),
    )
end

function _route_global_set_plan()
    direct_left = CPBRouteGlobalSet.cpb(
        1:1,
        1:2,
        1:1;
        role = :route_global_matrix_set_direct_left_source,
    )
    records = (
        _route_global_set_record(
            (:direct_diag, :direct_diag),
            1,
            :direct_direct,
            :direct_direct_pair_block_materialization_pilot;
            metadata = (;
                left_source_cpbs = (direct_left,),
                right_source_cpbs = (direct_left,),
                pair_index = 1,
                selector_family = :direct_direct,
                left_column_range = 1:2,
                right_column_range = 1:2,
            ),
        ),
    )
    return CPBMRouteGlobalSet.PairBlockMaterializationPlan(
        CPBMRouteGlobalSet.MetadataOnlyPairBlockMaterialization(),
        _route_global_set_pair_operator_plan(),
        records,
        (; status = :available_pair_block_materialization_plan),
        (; fixture = :route_global_matrix_set_smoke),
    )
end

function _route_global_set_inputs()
    return (;
        parent_axis_counts = (2, 2, 2),
        overlap_1d = (;
            x = [1.0 0.2; 0.2 1.1],
            y = [1.2 0.3; 0.3 1.3],
            z = [1.4 0.4; 0.4 1.5],
        ),
        position_1d = (;
            x = [5.0 0.7; 0.8 6.0],
            y = [7.0 0.9; 1.0 8.0],
            z = [9.0 1.1; 1.2 10.0],
        ),
    )
end

@testset "CartesianPairBlockMaterialization route global matrix-set smoke" begin
    plan = _route_global_set_plan()
    matrix_set = CPBMRouteGlobalSet.route_global_safe_one_body_matrices(
        (; pair_block_materialization_plan = plan);
        terms = (:overlap, :position_x),
        global_dimension = 4,
        factors = _route_global_set_inputs(),
        metadata = (; fixture = :route_global_matrix_set_smoke),
    )

    @test matrix_set.object_kind ===
          :cartesian_pair_block_route_global_safe_one_body_matrix_set
    @test matrix_set.status ===
          :materialized_route_global_safe_one_body_matrix_set
    @test isnothing(matrix_set.blocker)
    @test matrix_set.terms == (:overlap, :position_x)
    @test length(matrix_set.term_results) == 2
    @test matrix_set.materialized_term_count == 2
    @test matrix_set.blocked_term_count == 0
    @test matrix_set.all_terms_materialized
    @test matrix_set.global_one_body_term_matrices_materialized
    @test matrix_set.summary.materialized_terms == (:overlap, :position_x)
    @test matrix_set.summary.blocked_terms == ()

    overlap_result =
        CPBMRouteGlobalSet.route_global_one_body_matrix_set_result(
            matrix_set,
            :overlap,
        )
    position_x_result =
        CPBMRouteGlobalSet.route_global_one_body_matrix_set_result(
            matrix_set,
            :position_x,
        )
    @test overlap_result.status === :materialized_route_global_overlap_matrix
    @test position_x_result.status ===
          :materialized_route_global_position_x_matrix
    @test isnothing(
        CPBMRouteGlobalSet.route_global_one_body_matrix_set_result(
            matrix_set,
            :coulomb,
        ),
    )

    partial_set = CPBMRouteGlobalSet.route_global_one_body_matrix_set(
        plan;
        terms = (:overlap, :coulomb),
        global_dimension = 4,
        inputs = _route_global_set_inputs(),
    )
    @test partial_set.status === :partial_route_global_safe_one_body_matrix_set
    @test partial_set.blocker === :some_route_global_safe_one_body_terms_blocked
    @test partial_set.materialized_term_count == 1
    @test partial_set.blocked_term_count == 1
    @test partial_set.summary.materialized_terms == (:overlap,)
    @test partial_set.summary.blocked_terms == (:coulomb,)
    @test partial_set.summary.blocked_term_blocker_counts ==
          ((; blocker = :unsupported_route_global_one_body_term, count = 1),)

    missing_dimension_set =
        CPBMRouteGlobalSet.route_global_one_body_matrix_set(
            plan;
            terms = (:overlap, :position_x),
            inputs = _route_global_set_inputs(),
        )
    @test missing_dimension_set.status ===
          :blocked_route_global_safe_one_body_matrix_set
    @test missing_dimension_set.blocker === :missing_global_dimension
    @test missing_dimension_set.materialized_term_count == 0
    @test missing_dimension_set.blocked_term_count == 2
    @test missing_dimension_set.summary.blocked_term_blocker_counts ==
          ((; blocker = :missing_global_dimension, count = 2),)

    @test !matrix_set.route_driver_wiring
    @test !matrix_set.operator_blocks_materialized
    @test !matrix_set.hamiltonian_data_materialized
    @test !matrix_set.artifacts_materialized
    @test !matrix_set.exports_materialized
    @test !matrix_set.global_operator_blocks_materialized
    @test !matrix_set.global_hamiltonian_data_materialized
    @test !matrix_set.global_artifacts_materialized
    @test !matrix_set.coulomb_materialized
    @test !matrix_set.density_density_materialized
    @test !matrix_set.ida_mwg_data_materialized
    @test !matrix_set.pqs_lowdin_materialized
    @test !matrix_set.pqs_shell_projection_materialized
    @test !matrix_set.full_white_lindsey_route_assembled
end
