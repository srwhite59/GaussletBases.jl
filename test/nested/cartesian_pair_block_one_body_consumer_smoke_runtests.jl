# Runtime role: tiny smoke / routine per-pass.
#
# This is the preferred mixed one-body consumer validation for small edits.
# It uses synthetic metadata only, checks compact counts/flags, and avoids the
# real White-Lindsey adapter fixture. Use the contract/boundary tests for
# semantic changes or closeout validation.

using Test
using GaussletBases

const CPBMSmoke = GaussletBases.CartesianPairBlockMaterialization
const CPBSmoke = GaussletBases.CartesianCPB
const CTLSmoke = GaussletBases.CartesianTerminalLowering
const CRUSmoke = GaussletBases.CartesianRetainedUnits
const CRTCSmoke = GaussletBases.CartesianRetainedUnitTransformContracts
const CUPSmoke = GaussletBases.CartesianUnitPairs
const CPOPSmoke = GaussletBases.CartesianPairOperatorPlans

function _mixed_consumer_smoke_count(counts, field::Symbol, value)
    matches = Tuple(entry for entry in counts if getproperty(entry, field) == value)
    isempty(matches) && return 0
    return only(matches).count
end

function _mixed_consumer_smoke_pair_operator_plan()
    lowering_plan = CTLSmoke.TerminalLoweringPlan(
        CTLSmoke.PQSLowering(q = 2),
        (),
        (),
        (; status = :available_terminal_lowering_plan, materialized = false),
        (; fixture = :mixed_one_body_consumer_smoke),
    )
    retained_plan = CRUSmoke.RetainedUnitPlan(
        CRUSmoke.MetadataOnlyRetainedUnits(),
        lowering_plan,
        (),
        (; status = :available_retained_unit_plan, materialized = false),
        (; fixture = :mixed_one_body_consumer_smoke),
    )
    unit_pair_plan = CUPSmoke.UnitPairPlan(
        CUPSmoke.MetadataOnlyUnitPairs(),
        retained_plan,
        (),
        nothing,
        (; status = :available_unit_pair_plan, materialized = false),
        (; fixture = :mixed_one_body_consumer_smoke),
    )
    transform_plan = CRTCSmoke.RetainedUnitTransformContractPlan(
        CRTCSmoke.MetadataOnlyRetainedUnitTransformContracts(),
        retained_plan,
        (),
        (; status = :available_transform_contract_plan, materialized = false),
        (; fixture = :mixed_one_body_consumer_smoke),
    )
    return CPOPSmoke.PairOperatorPlan(
        CPOPSmoke.MetadataOnlyPairOperatorPlans(),
        unit_pair_plan,
        transform_plan,
        (),
        nothing,
        (; status = :available_pair_operator_plan, materialized = false),
        (; fixture = :mixed_one_body_consumer_smoke),
    )
end

function _mixed_consumer_smoke_record(
    pair_key::Tuple{Symbol,Symbol},
    pair_index::Int,
    pair_family::Symbol,
    materialization_path::Symbol;
    metadata = (;),
)
    return CPBMSmoke.PairBlockMaterializationRecord(
        pair_key,
        pair_index,
        pair_family,
        :synthetic_source_operator_path,
        (; left = :synthetic_left_transform, right = :synthetic_right_transform),
        (; left = :synthetic_left_realization, right = :synthetic_right_realization),
        :synthetic_final_block_path,
        (:overlap,),
        materialization_path,
        :ready_metadata_only_not_materialized,
        nothing,
        false,
        NamedTuple(metadata),
    )
end

function _mixed_consumer_smoke_plan()
    direct_left = CPBSmoke.cpb(
        1:1,
        1:2,
        1:1;
        role = :mixed_consumer_smoke_direct_left_source,
    )
    direct_right = CPBSmoke.cpb(
        2:2,
        1:2,
        1:1;
        role = :mixed_consumer_smoke_direct_right_source,
    )
    records = (
        _mixed_consumer_smoke_record(
            (:direct_left, :direct_right),
            1,
            :direct_direct,
            :direct_direct_pair_block_materialization_pilot;
            metadata = (;
                left_source_cpbs = (direct_left,),
                right_source_cpbs = (direct_right,),
            ),
        ),
        _mixed_consumer_smoke_record(
            (:pqs_left, :pqs_right),
            2,
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
        _mixed_consumer_smoke_record(
            (:lw_left, :lw_right),
            3,
            :white_lindsey_boundary_stratum,
            :white_lindsey_boundary_stratum_adapter_preflight,
        ),
        _mixed_consumer_smoke_record(
            (:unsupported_left, :unsupported_right),
            4,
            :unsupported_pair_family,
            :synthetic_unsupported_pair_block_path,
        ),
    )
    return CPBMSmoke.PairBlockMaterializationPlan(
        CPBMSmoke.MetadataOnlyPairBlockMaterialization(),
        _mixed_consumer_smoke_pair_operator_plan(),
        records,
        (; status = :available_pair_block_materialization_plan),
        (; fixture = :mixed_one_body_consumer_smoke),
    )
end

function _mixed_consumer_smoke_overlap_1d()
    return (;
        x = [1.0 0.2; 0.2 1.1],
        y = [1.2 0.3; 0.3 1.3],
        z = [1.4 0.4; 0.4 1.5],
    )
end

@testset "CartesianPairBlockMaterialization mixed one-body consumer smoke" begin
    consumption = CPBMSmoke._one_body_pair_block_consumption(
        _mixed_consumer_smoke_plan(),
        :overlap;
        inputs = (;
            parent_axis_counts = (2, 2, 2),
            overlap_1d = _mixed_consumer_smoke_overlap_1d(),
        ),
    )

    @test consumption.object_kind ==
          :cartesian_pair_block_mixed_one_body_consumption
    @test consumption.status ==
          :partially_materialized_mixed_one_body_pair_block_batch
    @test consumption.term == :overlap
    @test consumption.materialized_count == 2
    @test consumption.skipped_count == 2
    @test consumption.direct_direct_materialized
    @test consumption.pqs_source_pair_materialized
    @test !consumption.white_lindsey_materialized
    @test consumption.source_space_only_result_count == 1
    @test consumption.final_local_block_result_count == 1
    @test consumption.source_operator_blocks_materialized
    @test consumption.final_pair_blocks_materialized
    @test !consumption.operator_blocks_materialized
    @test !consumption.hamiltonian_data_materialized
    @test !consumption.artifacts_materialized
    @test !consumption.route_driver_wiring
    @test !consumption.factors_constructed

    summary = consumption.summary
    @test summary.pair_block_record_count == 4
    @test _mixed_consumer_smoke_count(
        summary.materialized_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 1
    @test _mixed_consumer_smoke_count(
        summary.materialized_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 1
    @test _mixed_consumer_smoke_count(
        summary.skipped_selector_family_counts,
        :selector_family,
        :white_lindsey_boundary_stratum,
    ) == 1
    @test _mixed_consumer_smoke_count(
        summary.skipped_selector_family_counts,
        :selector_family,
        :unsupported,
    ) == 1
    @test _mixed_consumer_smoke_count(
        summary.skipped_blocker_counts,
        :blocker,
        :missing_white_lindsey_unit_pair,
    ) == 1
    @test _mixed_consumer_smoke_count(
        summary.skipped_blocker_counts,
        :blocker,
        :unsupported_pair_block_materialization_path,
    ) == 1

    direct_result = only(
        result for result in consumption.batch_result.materialized_results
        if result.metadata.selector_family == :direct_direct
    )
    pqs_result = only(
        result for result in consumption.batch_result.materialized_results
        if result.metadata.selector_family == :pqs_source_pair
    )
    @test direct_result.final_pair_blocks_materialized
    @test pqs_result.source_operator_blocks_materialized
    @test !pqs_result.final_pair_blocks_materialized
end
