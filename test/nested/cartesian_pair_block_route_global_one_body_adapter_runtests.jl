# Runtime role: tiny smoke / route-shaped global overlap adapter contract.
#
# This verifies the private route-global one-body adapter for overlap and
# kinetic only. It does not cover Hamiltonians, Coulomb, IDA/MWG, PQS
# Lowdin/projection, exports, route-driver wiring, or White-Lindsey oracle
# fixtures.

using Test
using GaussletBases

const CPBMRouteGlobal = GaussletBases.CartesianPairBlockMaterialization
const CPBRouteGlobal = GaussletBases.CartesianCPB
const CTLRouteGlobal = GaussletBases.CartesianTerminalLowering
const CRURouteGlobal = GaussletBases.CartesianRetainedUnits
const CRTCRouteGlobal =
    GaussletBases.CartesianRetainedUnitTransformContracts
const CUPRouteGlobal = GaussletBases.CartesianUnitPairs
const CPOPRouteGlobal = GaussletBases.CartesianPairOperatorPlans

function _route_global_pair_operator_plan()
    lowering_plan = CTLRouteGlobal.TerminalLoweringPlan(
        CTLRouteGlobal.PQSLowering(q = 2),
        (),
        (),
        (; status = :available_terminal_lowering_plan, materialized = false),
        (; fixture = :route_global_one_body_adapter),
    )
    retained_plan = CRURouteGlobal.RetainedUnitPlan(
        CRURouteGlobal.MetadataOnlyRetainedUnits(),
        lowering_plan,
        (),
        (; status = :available_retained_unit_plan, materialized = false),
        (; fixture = :route_global_one_body_adapter),
    )
    unit_pair_plan = CUPRouteGlobal.UnitPairPlan(
        CUPRouteGlobal.MetadataOnlyUnitPairs(),
        retained_plan,
        (),
        nothing,
        (; status = :available_unit_pair_plan, materialized = false),
        (; fixture = :route_global_one_body_adapter),
    )
    transform_plan = CRTCRouteGlobal.RetainedUnitTransformContractPlan(
        CRTCRouteGlobal.MetadataOnlyRetainedUnitTransformContracts(),
        retained_plan,
        (),
        (; status = :available_transform_contract_plan, materialized = false),
        (; fixture = :route_global_one_body_adapter),
    )
    return CPOPRouteGlobal.PairOperatorPlan(
        CPOPRouteGlobal.MetadataOnlyPairOperatorPlans(),
        unit_pair_plan,
        transform_plan,
        (),
        nothing,
        (; status = :available_pair_operator_plan, materialized = false),
        (; fixture = :route_global_one_body_adapter),
    )
end

function _route_global_record(
    pair_key::Tuple{Symbol,Symbol},
    pair_index::Int,
    pair_family::Symbol,
    materialization_path::Symbol;
    supported_terms = (:overlap, :kinetic),
    metadata = (;),
)
    return CPBMRouteGlobal.PairBlockMaterializationRecord(
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

function _route_global_plan()
    direct_left = CPBRouteGlobal.cpb(
        1:1,
        1:2,
        1:1;
        role = :route_global_direct_left_source,
    )
    direct_right = CPBRouteGlobal.cpb(
        2:2,
        1:2,
        1:1;
        role = :route_global_direct_right_source,
    )
    records = (
        _route_global_record(
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
        _route_global_record(
            (:direct_left, :direct_right),
            2,
            :direct_direct,
            :direct_direct_pair_block_materialization_pilot;
            metadata = (;
                left_source_cpbs = (direct_left,),
                right_source_cpbs = (direct_right,),
                pair_index = 2,
                selector_family = :direct_direct,
                left_column_range = 1:2,
                right_column_range = 3:4,
            ),
        ),
        _route_global_record(
            (:pqs_left, :pqs_right),
            3,
            :pqs_pqs,
            :pqs_source_pair_preflight;
            metadata = (;
                left_source_mode_dims = (2, 2, 2),
                right_source_mode_dims = (2, 2, 2),
                left_source_mode_count = 8,
                right_source_mode_count = 8,
                source_mode_ordering = :x_major_y_major_z_fast,
                left_source_mode_ordering = :x_major_y_major_z_fast,
                right_source_mode_ordering = :x_major_y_major_z_fast,
            ),
        ),
        _route_global_record(
            (:lw_left, :lw_right),
            4,
            :white_lindsey_boundary_stratum,
            :white_lindsey_boundary_stratum_adapter_preflight,
        ),
        _route_global_record(
            (:unsupported_left, :unsupported_right),
            5,
            :unsupported_pair_family,
            :synthetic_unsupported_pair_block_path,
        ),
    )
    return CPBMRouteGlobal.PairBlockMaterializationPlan(
        CPBMRouteGlobal.MetadataOnlyPairBlockMaterialization(),
        _route_global_pair_operator_plan(),
        records,
        (; status = :available_pair_block_materialization_plan),
        (; fixture = :route_global_one_body_adapter),
    )
end

function _route_global_overlap_1d()
    return (;
        x = [1.0 0.2; 0.2 1.1],
        y = [1.2 0.3; 0.3 1.3],
        z = [1.4 0.4; 0.4 1.5],
    )
end

function _route_global_inputs()
    return (;
        parent_axis_counts = (2, 2, 2),
        overlap_1d = _route_global_overlap_1d(),
        kinetic_1d = (;
            x = [2.0 0.5; 0.5 2.2],
            y = [3.0 0.6; 0.6 3.3],
            z = [4.0 0.7; 0.7 4.4],
        ),
    )
end

@testset "CartesianPairBlockMaterialization route global overlap adapter" begin
    plan = _route_global_plan()
    adapter = CPBMRouteGlobal.route_global_overlap_matrix(
        plan;
        global_dimension = 4,
        inputs = _route_global_inputs(),
        metadata = (; fixture = :route_global_one_body_adapter),
    )
    alias_adapter = CPBMRouteGlobal.route_global_one_body_matrix(
        (; pair_block_materialization_plan = plan);
        term = :overlap,
        global_dimension = 4,
        factors = _route_global_inputs(),
    )

    expected = [
        1.68 0.42 0.336 0.084
        0.42 1.82 0.084 0.364
        0.336 0.084 0.0 0.0
        0.084 0.364 0.0 0.0
    ]
    @test adapter.object_kind ===
          :cartesian_pair_block_route_global_one_body_matrix_adapter
    @test adapter.status === :materialized_route_global_overlap_matrix
    @test isnothing(adapter.blocker)
    @test adapter.term === :overlap
    @test adapter.global_one_body_term_matrix_materialized
    @test adapter.global_overlap_matrix_materialized
    @test !adapter.global_kinetic_matrix_materialized
    @test adapter.operator_matrix_materialized
    @test adapter.global_operator_assembled
    @test adapter.global_matrix_result.matrix ≈ expected
    @test alias_adapter.global_matrix_result.matrix ≈ expected

    @test adapter.materialized_local_block_count == 3
    @test adapter.skipped_local_block_count == 2
    @test adapter.placeable_record_count == 2
    @test adapter.blocked_placement_count == 3
    @test adapter.placed_block_count == 3
    @test adapter.skipped_block_count == 3
    @test adapter.global_matrix_result.global_dimension == 4

    pqs_record = only(
        record for record in adapter.placement_plan.blocked_records
        if record.selector_family === :pqs_source_pair
    )
    @test pqs_record.blocker === :source_space_block_requires_shell_realization
    @test !any(value -> value > 8.0, adapter.global_matrix_result.matrix)

    missing_dimension = CPBMRouteGlobal.route_global_overlap_matrix(
        plan;
        inputs = _route_global_inputs(),
    )
    @test missing_dimension.status === :blocked_route_global_overlap_matrix
    @test missing_dimension.blocker === :missing_global_dimension
    @test !missing_dimension.global_one_body_term_matrix_materialized

    unsupported = CPBMRouteGlobal.route_global_one_body_matrix(
        plan;
        term = :position_x,
        global_dimension = 4,
        inputs = _route_global_inputs(),
    )
    @test unsupported.status === :blocked_route_global_one_body_matrix
    @test unsupported.blocker === :unsupported_route_global_one_body_term
    @test !unsupported.global_one_body_term_matrix_materialized

    no_ranges_plan = CPBMRouteGlobal.PairBlockMaterializationPlan(
        CPBMRouteGlobal.MetadataOnlyPairBlockMaterialization(),
        _route_global_pair_operator_plan(),
        (
            _route_global_record(
                (:direct_left, :direct_right),
                1,
                :direct_direct,
                :direct_direct_pair_block_materialization_pilot;
                metadata = (;
                    left_source_cpbs = (
                        CPBRouteGlobal.cpb(1:1, 1:2, 1:1),
                    ),
                    right_source_cpbs = (
                        CPBRouteGlobal.cpb(2:2, 1:2, 1:1),
                    ),
                ),
            ),
        ),
        (; status = :available_pair_block_materialization_plan),
        (; fixture = :route_global_one_body_missing_ranges),
    )
    no_ranges = CPBMRouteGlobal.route_global_overlap_matrix(
        no_ranges_plan;
        global_dimension = 4,
        inputs = _route_global_inputs(),
    )
    @test no_ranges.status === :blocked_route_global_overlap_matrix
    @test no_ranges.blocker === :no_placeable_overlap_blocks
    @test no_ranges.placement_plan.blocker === :missing_column_ranges

    @test !adapter.route_driver_wiring
    @test !adapter.operator_blocks_materialized
    @test !adapter.hamiltonian_data_materialized
    @test !adapter.artifacts_materialized
    @test !adapter.exports_materialized
    @test !adapter.global_operator_blocks_materialized
    @test !adapter.global_hamiltonian_data_materialized
    @test !adapter.global_artifacts_materialized
    @test !adapter.coulomb_materialized
    @test !adapter.density_density_materialized
    @test !adapter.ida_mwg_data_materialized
    @test !adapter.pqs_lowdin_materialized
    @test !adapter.pqs_shell_projection_materialized
    @test !adapter.full_white_lindsey_route_assembled
end

@testset "CartesianPairBlockMaterialization route global kinetic adapter" begin
    plan = _route_global_plan()
    adapter = CPBMRouteGlobal.route_global_kinetic_matrix(
        plan;
        global_dimension = 4,
        inputs = _route_global_inputs(),
    )
    alias_adapter = CPBMRouteGlobal.route_global_one_body_matrix(
        (; pair_block_materialization_plan = plan);
        term = :kinetic,
        global_dimension = 4,
        factors = _route_global_inputs(),
    )

    expected = [
        12.36 2.88 2.64 0.618
        2.88 13.46 0.618 2.874
        2.64 0.618 0.0 0.0
        0.618 2.874 0.0 0.0
    ]
    @test adapter.object_kind ===
          :cartesian_pair_block_route_global_one_body_matrix_adapter
    @test adapter.status === :materialized_route_global_kinetic_matrix
    @test isnothing(adapter.blocker)
    @test adapter.term === :kinetic
    @test adapter.global_one_body_term_matrix_materialized
    @test !adapter.global_overlap_matrix_materialized
    @test adapter.global_kinetic_matrix_materialized
    @test adapter.operator_matrix_materialized
    @test adapter.global_operator_assembled
    @test adapter.global_matrix_result.matrix ≈ expected
    @test alias_adapter.global_matrix_result.matrix ≈ expected

    @test adapter.materialized_local_block_count == 3
    @test adapter.skipped_local_block_count == 2
    @test adapter.placeable_record_count == 2
    @test adapter.blocked_placement_count == 3
    @test adapter.placed_block_count == 3
    @test adapter.skipped_block_count == 3
    @test adapter.global_matrix_result.global_dimension == 4

    pqs_record = only(
        record for record in adapter.placement_plan.blocked_records
        if record.selector_family === :pqs_source_pair
    )
    @test pqs_record.blocker === :source_space_block_requires_shell_realization
    @test !any(value -> value > 13.5, adapter.global_matrix_result.matrix)

    missing_dimension = CPBMRouteGlobal.route_global_kinetic_matrix(
        plan;
        inputs = _route_global_inputs(),
    )
    @test missing_dimension.status === :blocked_route_global_kinetic_matrix
    @test missing_dimension.blocker === :missing_global_dimension
    @test !missing_dimension.global_one_body_term_matrix_materialized

    @test !adapter.route_driver_wiring
    @test !adapter.operator_blocks_materialized
    @test !adapter.hamiltonian_data_materialized
    @test !adapter.artifacts_materialized
    @test !adapter.exports_materialized
    @test !adapter.global_operator_blocks_materialized
    @test !adapter.global_hamiltonian_data_materialized
    @test !adapter.global_artifacts_materialized
    @test !adapter.coulomb_materialized
    @test !adapter.density_density_materialized
    @test !adapter.ida_mwg_data_materialized
    @test !adapter.pqs_lowdin_materialized
    @test !adapter.pqs_shell_projection_materialized
    @test !adapter.full_white_lindsey_route_assembled
end
