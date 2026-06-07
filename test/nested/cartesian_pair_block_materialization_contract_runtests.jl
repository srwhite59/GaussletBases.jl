using Test
using GaussletBases

const CPBM = GaussletBases.CartesianPairBlockMaterialization
const CPOPForPairBlocks = GaussletBases.CartesianPairOperatorPlans
const CUPForPairBlocks = GaussletBases.CartesianUnitPairs
const CRTCForPairBlocks = GaussletBases.CartesianRetainedUnitTransformContracts
const CRUForPairBlocks = GaussletBases.CartesianRetainedUnits
const CTLForPairBlocks = GaussletBases.CartesianTerminalLowering
const CRCForPairBlocks = GaussletBases.CartesianRouteCore
const CPBForPairBlocks = GaussletBases.CartesianCPB

function _pair_block_count(counts, field::Symbol, value)
    matches = Tuple(entry for entry in counts if getproperty(entry, field) == value)
    isempty(matches) && return 0
    return only(matches).count
end

function _pair_block_minimal_lowering_plan()
    return CTLForPairBlocks.TerminalLoweringPlan(
        CTLForPairBlocks.PQSLowering(q = 3),
        (),
        (),
        (;
            object_kind = :synthetic_terminal_lowering_plan_summary,
            status = :available_terminal_lowering_plan,
            policy_kind = :pair_block_materialization_synthetic_lowering,
            terminal_region_count = 0,
            available_contract_count = 0,
            selected_contract_count = 0,
            selected_contract_kinds = (),
            all_terminal_regions_have_selected_contract = true,
            materialized = false,
            retained_spaces_materialized = false,
            final_retained_units_materialized = false,
            pair_inventory_materialized = false,
            operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
        (; fixture = :cartesian_pair_block_materialization_contract),
    )
end

function _pair_block_retained_unit(
    unit_key::Symbol,
    unit_index::Int,
    unit_kind::Symbol,
    lowering_kind::Symbol,
    retained_rule::Symbol,
    realization_rule;
    owned_support = nothing,
    source_cpbs = (),
    source_cpb_index = nothing,
)
    return CRUForPairBlocks.RetainedUnitRecord(
        unit_key,
        unit_index,
        unit_kind,
        Symbol(unit_key, "_contract"),
        Symbol(unit_key, "_terminal_region"),
        :synthetic_terminal_region,
        :synthetic_terminal_region,
        lowering_kind,
        retained_rule,
        realization_rule,
        owned_support,
        Tuple(source_cpbs),
        source_cpb_index,
        :not_materialized,
        nothing,
        :not_materialized,
        nothing,
        nothing,
        false,
        (; route_core_sidecar_status = :not_materialized),
    )
end

function _pair_block_retained_plan()
    direct_source = CPBForPairBlocks.filled_cpb(
        1:2,
        2:3,
        1:2;
        role = :pair_block_direct_source_cpb,
    )
    direct = _pair_block_retained_unit(
        :pair_block_direct_unit,
        1,
        :direct_cpb_retained_unit,
        :direct_core_identity_cpb,
        :direct_source_modes,
        :direct_or_trivial_embedding;
        owned_support = CRCForPairBlocks.owned_cpb(direct_source),
        source_cpbs = (direct_source,),
    )
    pqs = _pair_block_retained_unit(
        :pair_block_pqs_unit,
        2,
        :pqs_shell_retained_unit,
        :pqs_filled_source_cpb,
        :pqs_boundary_comx_product_modes,
        :shell_projection_lowdin,
    )
    units = (direct, pqs)
    return CRUForPairBlocks.RetainedUnitPlan(
        CRUForPairBlocks.MetadataOnlyRetainedUnits(),
        _pair_block_minimal_lowering_plan(),
        units,
        (;
            object_kind = :synthetic_retained_unit_plan_summary,
            status = :available_retained_unit_plan,
            retained_unit_count = length(units),
            materialized = false,
            transforms_materialized = false,
            coefficient_maps_materialized = false,
            pair_inventory_materialized = false,
            operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
        (; fixture = :cartesian_pair_block_materialization_contract),
    )
end

function _pair_block_rectangular_retained_plan()
    left_source = CPBForPairBlocks.cpb(
        1:2,
        1:2,
        1:1;
        role = :pair_block_direct_left_source_cpb,
    )
    right_source = CPBForPairBlocks.cpb(
        3:3,
        2:3,
        1:3;
        role = :pair_block_direct_right_source_cpb,
    )
    left_direct = _pair_block_retained_unit(
        :pair_block_direct_left_unit,
        1,
        :direct_cpb_retained_unit,
        :direct_core_identity_cpb,
        :direct_source_modes,
        :direct_or_trivial_embedding;
        owned_support = CRCForPairBlocks.owned_cpb(left_source),
        source_cpbs = (left_source,),
    )
    right_direct = _pair_block_retained_unit(
        :pair_block_direct_right_unit,
        2,
        :direct_cpb_retained_unit,
        :direct_core_identity_cpb,
        :direct_source_modes,
        :direct_or_trivial_embedding;
        owned_support = CRCForPairBlocks.owned_cpb(right_source),
        source_cpbs = (right_source,),
    )
    pqs = _pair_block_retained_unit(
        :pair_block_rectangular_pqs_unit,
        3,
        :pqs_shell_retained_unit,
        :pqs_filled_source_cpb,
        :pqs_boundary_comx_product_modes,
        :shell_projection_lowdin,
    )
    units = (left_direct, right_direct, pqs)
    return CRUForPairBlocks.RetainedUnitPlan(
        CRUForPairBlocks.MetadataOnlyRetainedUnits(),
        _pair_block_minimal_lowering_plan(),
        units,
        (;
            object_kind = :synthetic_retained_unit_plan_summary,
            status = :available_retained_unit_plan,
            retained_unit_count = length(units),
            materialized = false,
            transforms_materialized = false,
            coefficient_maps_materialized = false,
            pair_inventory_materialized = false,
            operator_blocks_materialized = false,
            hamiltonian_data_materialized = false,
        ),
        (; fixture = :cartesian_pair_block_rectangular_direct_direct),
    )
end

function _pair_block_record_for(plan, left_key::Symbol, right_key::Symbol)
    return only(
        record for record in CPBM.pair_block_materialization_records(plan)
        if record.pair_key == (left_key, right_key)
    )
end

function _pair_block_overlap_axes()
    overlap_x = [
        1.00 0.11 0.13 0.17
        0.11 1.20 0.19 0.23
        0.13 0.19 1.40 0.29
        0.17 0.23 0.29 1.60
    ]
    overlap_y = [
        1.70 0.31 0.37 0.41
        0.31 1.90 0.43 0.47
        0.37 0.43 2.10 0.53
        0.41 0.47 0.53 2.30
    ]
    overlap_z = [
        2.50 0.59 0.61 0.67
        0.59 2.70 0.71 0.73
        0.61 0.71 2.90 0.79
        0.67 0.73 0.79 3.10
    ]
    return overlap_x, overlap_y, overlap_z
end

function _pair_block_position_axes()
    position_x = [
        0.20 1.10 1.30 1.70
        1.10 0.40 1.90 2.30
        1.30 1.90 0.60 2.90
        1.70 2.30 2.90 0.80
    ]
    position_y = [
        3.10 0.21 0.23 0.27
        0.21 3.30 0.29 0.31
        0.23 0.29 3.50 0.37
        0.27 0.31 0.37 3.70
    ]
    position_z = [
        4.10 0.41 0.43 0.47
        0.41 4.30 0.49 0.51
        0.43 0.49 4.50 0.57
        0.47 0.51 0.57 4.70
    ]
    return position_x, position_y, position_z
end

function _pair_block_x2_axes()
    x2_x = [
        5.10 0.61 0.63 0.67
        0.61 5.30 0.69 0.71
        0.63 0.69 5.50 0.77
        0.67 0.71 0.77 5.70
    ]
    x2_y = [
        6.10 0.81 0.83 0.87
        0.81 6.30 0.89 0.91
        0.83 0.89 6.50 0.97
        0.87 0.91 0.97 6.70
    ]
    x2_z = [
        7.10 1.01 1.03 1.07
        1.01 7.30 1.09 1.11
        1.03 1.09 7.50 1.17
        1.07 1.11 1.17 7.70
    ]
    return x2_x, x2_y, x2_z
end

function _pair_block_kinetic_axes()
    kinetic_x = [
        -1.10 0.12 0.14 0.18
        0.12 -1.30 0.20 0.24
        0.14 0.20 -1.50 0.30
        0.18 0.24 0.30 -1.70
    ]
    kinetic_y = [
        -2.10 0.32 0.34 0.38
        0.32 -2.30 0.40 0.44
        0.34 0.40 -2.50 0.50
        0.38 0.44 0.50 -2.70
    ]
    kinetic_z = [
        -3.10 0.52 0.54 0.58
        0.52 -3.30 0.60 0.64
        0.54 0.60 -3.50 0.70
        0.58 0.64 0.70 -3.70
    ]
    return kinetic_x, kinetic_y, kinetic_z
end

function _pair_block_cpb_states(source_cpb)
    ix, iy, iz = CPBForPairBlocks.intervals(source_cpb)
    return Tuple((x, y, z) for x in ix for y in iy for z in iz)
end

function _pair_block_expected_product(left_cpb, right_cpb, operator_x, operator_y, operator_z)
    left_states = _pair_block_cpb_states(left_cpb)
    right_states = _pair_block_cpb_states(right_cpb)
    return [
        operator_x[left[1], right[1]] *
        operator_y[left[2], right[2]] *
        operator_z[left[3], right[3]]
        for left in left_states, right in right_states
    ]
end

function _pair_block_expected_overlap(left_cpb, right_cpb, overlap_x, overlap_y, overlap_z)
    return _pair_block_expected_product(left_cpb, right_cpb, overlap_x, overlap_y, overlap_z)
end

function _pair_block_expected_position(
    left_cpb,
    right_cpb,
    axis,
    overlap_x,
    overlap_y,
    overlap_z,
    position_x,
    position_y,
    position_z,
)
    operator_x, operator_y, operator_z =
        axis === :x ? (position_x, overlap_y, overlap_z) :
        axis === :y ? (overlap_x, position_y, overlap_z) :
        (overlap_x, overlap_y, position_z)
    return _pair_block_expected_product(left_cpb, right_cpb, operator_x, operator_y, operator_z)
end

function _pair_block_expected_x2(
    left_cpb,
    right_cpb,
    axis,
    overlap_x,
    overlap_y,
    overlap_z,
    x2_x,
    x2_y,
    x2_z,
)
    operator_x, operator_y, operator_z =
        axis === :x ? (x2_x, overlap_y, overlap_z) :
        axis === :y ? (overlap_x, x2_y, overlap_z) :
        (overlap_x, overlap_y, x2_z)
    return _pair_block_expected_product(left_cpb, right_cpb, operator_x, operator_y, operator_z)
end

function _pair_block_expected_kinetic(
    left_cpb,
    right_cpb,
    overlap_x,
    overlap_y,
    overlap_z,
    kinetic_x,
    kinetic_y,
    kinetic_z,
)
    return _pair_block_expected_product(
        left_cpb,
        right_cpb,
        kinetic_x,
        overlap_y,
        overlap_z,
    ) +
           _pair_block_expected_product(
        left_cpb,
        right_cpb,
        overlap_x,
        kinetic_y,
        overlap_z,
    ) +
           _pair_block_expected_product(
        left_cpb,
        right_cpb,
        overlap_x,
        overlap_y,
        kinetic_z,
    )
end

@testset "CartesianPairBlockMaterialization unavailable summary" begin
    summary = CPBM.unavailable_summary(:not_selected, :not_selected_route)

    @test summary.object_kind == :cartesian_pair_block_materialization_plan_summary
    @test summary.status == :not_selected
    @test summary.blocker == :not_selected_route
    @test summary.pair_operator_plan_count == 0
    @test summary.pair_block_record_count == 0
    @test summary.ready_record_count == 0
    @test summary.blocked_record_count == 0
    @test summary.materialization_path_counts == ()
    @test summary.readiness_status_counts == ()
    @test summary.blocker_counts == ()
    @test !summary.materialized
    @test !summary.source_operator_blocks_materialized
    @test !summary.final_pair_blocks_materialized
    @test !summary.operator_blocks_materialized
    @test !summary.hamiltonian_data_materialized
    @test !summary.artifacts_materialized
end

@testset "CartesianPairBlockMaterialization direct/direct preflight" begin
    retained_plan = _pair_block_retained_plan()
    unit_pair_plan = CUPForPairBlocks.unit_pair_plan(retained_plan)
    transform_plan =
        CRTCForPairBlocks.retained_unit_transform_contract_plan(retained_plan)
    pair_operator_plan =
        CPOPForPairBlocks.pair_operator_plan(
            unit_pair_plan,
            transform_plan;
            route_core_sidecars = false,
        )
    materialization_plan =
        CPBM.pair_block_materialization_plan(pair_operator_plan)
    materialization_summary = CPBM.summary(materialization_plan)

    @test materialization_plan isa CPBM.PairBlockMaterializationPlan
    @test length(CPBM.pair_block_materialization_records(materialization_plan)) == 3
    @test materialization_summary.status == :blocked_pair_block_materialization_plan
    @test materialization_summary.pair_operator_plan_count == 3
    @test materialization_summary.pair_block_record_count == 3
    @test materialization_summary.ready_record_count == 1
    @test materialization_summary.blocked_record_count == 2
    @test _pair_block_count(
        materialization_summary.materialization_path_counts,
        :materialization_path,
        :direct_direct_pair_block_materialization_pilot,
    ) == 1
    @test _pair_block_count(
        materialization_summary.materialization_path_counts,
        :materialization_path,
        :deferred_pair_block_materialization_path,
    ) == 2
    @test _pair_block_count(
        materialization_summary.readiness_status_counts,
        :readiness_status,
        :ready_metadata_only_not_materialized,
    ) == 1
    @test _pair_block_count(
        materialization_summary.readiness_status_counts,
        :readiness_status,
        :blocked_pair_block_materialization_not_implemented,
    ) == 2
    @test _pair_block_count(
        materialization_summary.blocker_counts,
        :blocker,
        :non_direct_direct_pair_block_materialization_not_implemented,
    ) == 2

    direct_record = _pair_block_record_for(
        materialization_plan,
        :pair_block_direct_unit,
        :pair_block_direct_unit,
    )
    direct_pqs_record = _pair_block_record_for(
        materialization_plan,
        :pair_block_direct_unit,
        :pair_block_pqs_unit,
    )
    pqs_pqs_record = _pair_block_record_for(
        materialization_plan,
        :pair_block_pqs_unit,
        :pair_block_pqs_unit,
    )

    @test direct_record.materialization_path ==
          :direct_direct_pair_block_materialization_pilot
    @test direct_record.readiness_status == :ready_metadata_only_not_materialized
    @test direct_record.blocker === nothing
    @test !direct_record.materialized

    overlap_x, overlap_y, overlap_z = _pair_block_overlap_axes()
    overlap_result = CPBM.direct_direct_overlap_block(
        direct_record;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
    )
    direct_source = only(direct_record.metadata.left_source_cpbs)
    expected_overlap =
        _pair_block_expected_overlap(
            direct_source,
            direct_source,
            overlap_x,
            overlap_y,
            overlap_z,
        )

    @test overlap_result isa CPBM.PairBlockMaterializationResult
    @test overlap_result.term == :overlap
    @test overlap_result.pair_key == (:pair_block_direct_unit, :pair_block_direct_unit)
    @test size(overlap_result.block) == (8, 8)
    @test overlap_result.block ≈ expected_overlap
    @test overlap_result.materialized
    @test overlap_result.source_operator_blocks_materialized
    @test overlap_result.final_pair_blocks_materialized
    @test !overlap_result.operator_blocks_materialized
    @test !overlap_result.hamiltonian_data_materialized
    @test !overlap_result.artifacts_materialized

    @test direct_pqs_record.readiness_status ==
          :blocked_pair_block_materialization_not_implemented
    @test direct_pqs_record.blocker ==
          :non_direct_direct_pair_block_materialization_not_implemented
    @test_throws ArgumentError CPBM.direct_direct_overlap_block(
        direct_pqs_record;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
    )
    position_x, position_y, position_z = _pair_block_position_axes()
    @test_throws ArgumentError CPBM.direct_direct_position_block(
        direct_pqs_record;
        axis = :x,
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        position_1d = (; x = position_x, y = position_y, z = position_z),
    )
    x2_x, x2_y, x2_z = _pair_block_x2_axes()
    @test_throws ArgumentError CPBM.direct_direct_x2_block(
        direct_pqs_record;
        axis = :x,
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        x2_1d = (; x = x2_x, y = x2_y, z = x2_z),
    )
    kinetic_x, kinetic_y, kinetic_z = _pair_block_kinetic_axes()
    @test_throws ArgumentError CPBM.direct_direct_kinetic_block(
        direct_pqs_record;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        kinetic_1d = (; x = kinetic_x, y = kinetic_y, z = kinetic_z),
    )
    @test pqs_pqs_record.readiness_status ==
          :blocked_pair_block_materialization_not_implemented
    @test pqs_pqs_record.blocker ==
          :non_direct_direct_pair_block_materialization_not_implemented
    @test !materialization_summary.materialized
    @test !materialization_summary.source_operator_blocks_materialized
    @test !materialization_summary.final_pair_blocks_materialized
    @test !materialization_summary.operator_blocks_materialized
    @test !materialization_summary.hamiltonian_data_materialized
    @test !materialization_summary.artifacts_materialized
end

@testset "CartesianPairBlockMaterialization direct/direct overlap selector" begin
    retained_plan = _pair_block_rectangular_retained_plan()
    unit_pair_plan = CUPForPairBlocks.unit_pair_plan(retained_plan)
    transform_plan =
        CRTCForPairBlocks.retained_unit_transform_contract_plan(retained_plan)
    pair_operator_plan =
        CPOPForPairBlocks.pair_operator_plan(
            unit_pair_plan,
            transform_plan;
            route_core_sidecars = false,
        )
    materialization_plan =
        CPBM.pair_block_materialization_plan(pair_operator_plan)
    materialization_summary = CPBM.summary(materialization_plan)

    @test length(CPBM.pair_block_materialization_records(materialization_plan)) == 6
    @test materialization_summary.ready_record_count == 3
    @test materialization_summary.blocked_record_count == 3

    cross_record = _pair_block_record_for(
        materialization_plan,
        :pair_block_direct_left_unit,
        :pair_block_direct_right_unit,
    )
    overlap_x, overlap_y, overlap_z = _pair_block_overlap_axes()
    cross_result = CPBM.direct_direct_overlap_block(
        cross_record;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
    )
    left_source = only(cross_record.metadata.left_source_cpbs)
    right_source = only(cross_record.metadata.right_source_cpbs)
    expected_cross =
        _pair_block_expected_overlap(
            left_source,
            right_source,
            overlap_x,
            overlap_y,
            overlap_z,
        )

    @test size(cross_result.block) == (4, 6)
    @test cross_result.block ≈ expected_cross

    position_x, position_y, position_z = _pair_block_position_axes()
    expected_position_terms = (;
        x = :position_x,
        y = :position_y,
        z = :position_z,
    )
    for axis in (:x, :y, :z)
        position_result = CPBM.direct_direct_position_block(
            cross_record;
            axis,
            parent_axis_counts = (4, 4, 4),
            overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
            position_1d = (; x = position_x, y = position_y, z = position_z),
        )
        expected_position =
            _pair_block_expected_position(
                left_source,
                right_source,
                axis,
                overlap_x,
                overlap_y,
                overlap_z,
                position_x,
                position_y,
                position_z,
            )

        @test position_result.term == getproperty(expected_position_terms, axis)
        @test size(position_result.block) == (4, 6)
        @test position_result.block ≈ expected_position
        @test position_result.materialized
        @test position_result.source_operator_blocks_materialized
        @test position_result.final_pair_blocks_materialized
        @test !position_result.operator_blocks_materialized
        @test !position_result.hamiltonian_data_materialized
        @test !position_result.artifacts_materialized
    end

    x2_x, x2_y, x2_z = _pair_block_x2_axes()
    expected_x2_terms = (;
        x = :x2_x,
        y = :x2_y,
        z = :x2_z,
    )
    for axis in (:x, :y, :z)
        x2_result = CPBM.direct_direct_x2_block(
            cross_record;
            axis,
            parent_axis_counts = (4, 4, 4),
            overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
            x2_1d = (; x = x2_x, y = x2_y, z = x2_z),
        )
        expected_x2 =
            _pair_block_expected_x2(
                left_source,
                right_source,
                axis,
                overlap_x,
                overlap_y,
                overlap_z,
                x2_x,
                x2_y,
                x2_z,
            )

        @test x2_result.term == getproperty(expected_x2_terms, axis)
        @test size(x2_result.block) == (4, 6)
        @test x2_result.block ≈ expected_x2
        @test x2_result.materialized
        @test x2_result.source_operator_blocks_materialized
        @test x2_result.final_pair_blocks_materialized
        @test !x2_result.operator_blocks_materialized
        @test !x2_result.hamiltonian_data_materialized
        @test !x2_result.artifacts_materialized
    end

    kinetic_x, kinetic_y, kinetic_z = _pair_block_kinetic_axes()
    kinetic_result = CPBM.direct_direct_kinetic_block(
        cross_record;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        kinetic_1d = (; x = kinetic_x, y = kinetic_y, z = kinetic_z),
    )
    expected_kinetic =
        _pair_block_expected_kinetic(
            left_source,
            right_source,
            overlap_x,
            overlap_y,
            overlap_z,
            kinetic_x,
            kinetic_y,
            kinetic_z,
        )

    @test kinetic_result.term == :kinetic
    @test size(kinetic_result.block) == (4, 6)
    @test kinetic_result.block ≈ expected_kinetic
    @test kinetic_result.materialized
    @test kinetic_result.source_operator_blocks_materialized
    @test kinetic_result.final_pair_blocks_materialized
    @test !kinetic_result.operator_blocks_materialized
    @test !kinetic_result.hamiltonian_data_materialized
    @test !kinetic_result.artifacts_materialized

    selector_overlap = CPBM.direct_direct_one_body_block(
        cross_record,
        :overlap;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
    )
    position_y_result = CPBM.direct_direct_position_block(
        cross_record;
        axis = :y,
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        position_1d = (; x = position_x, y = position_y, z = position_z),
    )
    selector_position_y = CPBM.direct_direct_one_body_block(
        cross_record,
        :position_y;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        position_1d = (; x = position_x, y = position_y, z = position_z),
    )
    x2_z_result = CPBM.direct_direct_x2_block(
        cross_record;
        axis = :z,
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        x2_1d = (; x = x2_x, y = x2_y, z = x2_z),
    )
    selector_x2_z = CPBM.direct_direct_one_body_block(
        cross_record,
        :x2_z;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        x2_1d = (; x = x2_x, y = x2_y, z = x2_z),
    )
    selector_kinetic = CPBM.direct_direct_one_body_block(
        cross_record,
        :kinetic;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        kinetic_1d = (; x = kinetic_x, y = kinetic_y, z = kinetic_z),
    )

    @test selector_overlap.block ≈ cross_result.block
    @test selector_position_y.block ≈ position_y_result.block
    @test selector_x2_z.block ≈ x2_z_result.block
    @test selector_kinetic.block ≈ kinetic_result.block
    @test_throws ArgumentError CPBM.direct_direct_one_body_block(
        cross_record,
        :coulomb;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
    )
    @test_throws ArgumentError CPBM.direct_direct_one_body_block(
        cross_record,
        :position_x;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
    )
    @test_throws ArgumentError CPBM.direct_direct_one_body_block(
        cross_record,
        :x2_x;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
    )
    @test_throws ArgumentError CPBM.direct_direct_one_body_block(
        cross_record,
        :kinetic;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
    )

    batch_result = CPBM.direct_direct_overlap_blocks(
        materialization_plan;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
    )
    position_batch_result = CPBM.direct_direct_position_blocks(
        materialization_plan;
        axis = :z,
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        position_1d = (; x = position_x, y = position_y, z = position_z),
    )
    x2_batch_result = CPBM.direct_direct_x2_blocks(
        materialization_plan;
        axis = :z,
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        x2_1d = (; x = x2_x, y = x2_y, z = x2_z),
    )
    kinetic_batch_result = CPBM.direct_direct_kinetic_blocks(
        materialization_plan;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        kinetic_1d = (; x = kinetic_x, y = kinetic_y, z = kinetic_z),
    )
    selector_kinetic_batch = CPBM.direct_direct_one_body_blocks(
        materialization_plan,
        :kinetic;
        parent_axis_counts = (4, 4, 4),
        overlap_1d = (; x = overlap_x, y = overlap_y, z = overlap_z),
        kinetic_1d = (; x = kinetic_x, y = kinetic_y, z = kinetic_z),
    )
    batch_cross = only(
        result for result in batch_result.materialized_results
        if result.pair_key == (:pair_block_direct_left_unit, :pair_block_direct_right_unit)
    )
    materialized_keys = Set(result.pair_key for result in batch_result.materialized_results)

    @test batch_result isa CPBM.PairBlockMaterializationBatchResult
    @test batch_result.term == :overlap
    @test batch_result.materialized_count == 3
    @test batch_result.skipped_count == 3
    @test materialized_keys == Set((
        (:pair_block_direct_left_unit, :pair_block_direct_left_unit),
        (:pair_block_direct_left_unit, :pair_block_direct_right_unit),
        (:pair_block_direct_right_unit, :pair_block_direct_right_unit),
    ))
    @test all(
        skipped -> skipped.blocker ==
                   :non_direct_direct_pair_block_materialization_not_implemented,
        batch_result.skipped_records,
    )
    @test batch_cross.block ≈ expected_cross
    @test batch_result.materialized
    @test batch_result.source_operator_blocks_materialized
    @test batch_result.final_pair_blocks_materialized
    @test !batch_result.operator_blocks_materialized
    @test !batch_result.hamiltonian_data_materialized
    @test !batch_result.artifacts_materialized

    @test position_batch_result.term == :position_z
    @test position_batch_result.materialized_count == 3
    @test position_batch_result.skipped_count == 3
    @test position_batch_result.materialized
    @test position_batch_result.source_operator_blocks_materialized
    @test position_batch_result.final_pair_blocks_materialized
    @test !position_batch_result.operator_blocks_materialized
    @test !position_batch_result.hamiltonian_data_materialized
    @test !position_batch_result.artifacts_materialized

    @test x2_batch_result.term == :x2_z
    @test x2_batch_result.materialized_count == 3
    @test x2_batch_result.skipped_count == 3
    @test x2_batch_result.materialized
    @test x2_batch_result.source_operator_blocks_materialized
    @test x2_batch_result.final_pair_blocks_materialized
    @test !x2_batch_result.operator_blocks_materialized
    @test !x2_batch_result.hamiltonian_data_materialized
    @test !x2_batch_result.artifacts_materialized

    @test kinetic_batch_result.term == :kinetic
    @test kinetic_batch_result.materialized_count == 3
    @test kinetic_batch_result.skipped_count == 3
    @test kinetic_batch_result.materialized
    @test kinetic_batch_result.source_operator_blocks_materialized
    @test kinetic_batch_result.final_pair_blocks_materialized
    @test !kinetic_batch_result.operator_blocks_materialized
    @test !kinetic_batch_result.hamiltonian_data_materialized
    @test !kinetic_batch_result.artifacts_materialized

    @test selector_kinetic_batch.materialized_count ==
          kinetic_batch_result.materialized_count
    @test selector_kinetic_batch.skipped_count == kinetic_batch_result.skipped_count
    @test selector_kinetic_batch.term == kinetic_batch_result.term
end
