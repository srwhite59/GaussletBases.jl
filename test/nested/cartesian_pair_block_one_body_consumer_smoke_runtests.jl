# Runtime role: tiny smoke / routine per-pass.
#
# This is the preferred mixed one-body consumer gate for small status/count
# edits. It uses synthetic metadata only, checks compact summaries/flags, and
# avoids the real White-Lindsey adapter fixture. Use the block-set
# contract/boundary tests for semantic changes or baton closeout validation.

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
    return sum(entry -> hasproperty(entry, :count) ? entry.count : 1, matches)
end

function _mixed_consumer_smoke_compact_summary(summary)
    forbidden_fields = (
        :term_batch_results,
        :batch_result,
        :block,
        :blocks,
        :materialized_results,
        :skipped_records,
        :preflight_summary,
        :block_set_summary,
        :term_summaries,
        :consumption_summary,
    )
    return !any(field -> hasproperty(summary, field), forbidden_fields)
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

    block_set_consumption = CPBMSmoke._one_body_pair_block_set_consumption(
        _mixed_consumer_smoke_plan();
        terms = (:overlap, :kinetic),
        inputs = (;
            parent_axis_counts = (2, 2, 2),
            overlap_1d = _mixed_consumer_smoke_overlap_1d(),
        ),
        materialize_terms = (:overlap,),
    )
    block_set_summary =
        CPBMSmoke._one_body_pair_block_set_consumption_summary(
            block_set_consumption,
        )
    @test block_set_summary.object_kind ==
          :cartesian_pair_block_mixed_one_body_block_set_consumption_summary
    @test block_set_summary.status ==
          :partially_materialized_mixed_one_body_block_set_consumption
    @test isnothing(block_set_summary.blocker)
    @test block_set_summary.requested_terms == (:overlap, :kinetic)
    @test block_set_summary.requested_materialize_terms == (:overlap,)
    @test block_set_summary.materialized_terms == (:overlap,)
    @test block_set_summary.deferred_terms == (:kinetic,)
    @test block_set_summary.preflight_status ==
          :partially_ready_mixed_one_body_block_set_preflight
    @test block_set_summary.block_set_summary_status ==
          :partially_deferred_mixed_one_body_block_set
    @test block_set_summary.supplied_term_count == 1
    @test block_set_summary.deferred_term_count == 1
    @test block_set_summary.total_materialized_count == 2
    @test block_set_summary.total_skipped_count == 2
    @test _mixed_consumer_smoke_count(
        block_set_summary.materialized_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 1
    @test _mixed_consumer_smoke_count(
        block_set_summary.materialized_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 1
    @test _mixed_consumer_smoke_count(
        block_set_summary.skipped_selector_family_counts,
        :selector_family,
        :white_lindsey_boundary_stratum,
    ) == 1
    @test _mixed_consumer_smoke_count(
        block_set_summary.skipped_selector_family_counts,
        :selector_family,
        :unsupported,
    ) == 1
    @test _mixed_consumer_smoke_count(
        block_set_summary.skipped_blocker_counts,
        :blocker,
        :missing_white_lindsey_unit_pair,
    ) == 1
    @test _mixed_consumer_smoke_count(
        block_set_summary.skipped_blocker_counts,
        :blocker,
        :unsupported_pair_block_materialization_path,
    ) == 1
    @test block_set_summary.term_batch_results_available_on_consumption
    @test !block_set_summary.term_batch_results_stored_in_summary
    @test !block_set_summary.factors_constructed
    @test block_set_summary.source_operator_blocks_materialized
    @test block_set_summary.final_pair_blocks_materialized
    @test !block_set_summary.operator_blocks_materialized
    @test !block_set_summary.hamiltonian_data_materialized
    @test !block_set_summary.artifacts_materialized
    @test !block_set_summary.route_driver_wiring
    @test !block_set_summary.coulomb_materialized
    @test !block_set_summary.ida_mwg_data_materialized
    @test !block_set_summary.pqs_lowdin_materialized
    @test !block_set_summary.full_white_lindsey_route_assembled
    @test _mixed_consumer_smoke_compact_summary(block_set_summary)
    @test all(
        term_status -> !hasproperty(term_status, :batch_result) &&
                       !hasproperty(term_status, :block),
        block_set_summary.term_statuses,
    )

    block_set_view =
        CPBMSmoke._one_body_pair_block_set_view(block_set_consumption)
    @test block_set_view.object_kind ==
          :cartesian_pair_block_mixed_one_body_block_set_view
    @test block_set_view.summary_object_kind ==
          :cartesian_pair_block_mixed_one_body_block_set_consumption_summary
    @test block_set_view.view_status == :available_mixed_one_body_block_set_view
    @test block_set_view.status ==
          :partially_materialized_mixed_one_body_block_set_consumption
    @test isnothing(block_set_view.blocker)
    @test block_set_view.requested_terms == (:overlap, :kinetic)
    @test block_set_view.requested_materialize_terms == (:overlap,)
    @test block_set_view.materialized_terms == (:overlap,)
    @test block_set_view.deferred_terms == (:kinetic,)
    @test block_set_view.term_statuses == block_set_summary.term_statuses
    @test block_set_view.total_materialized_count == 2
    @test block_set_view.total_skipped_count == 2
    @test _mixed_consumer_smoke_count(
        block_set_view.materialized_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 1
    @test _mixed_consumer_smoke_count(
        block_set_view.materialized_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 1
    @test _mixed_consumer_smoke_count(
        block_set_view.skipped_selector_family_counts,
        :selector_family,
        :white_lindsey_boundary_stratum,
    ) == 1
    @test _mixed_consumer_smoke_count(
        block_set_view.skipped_selector_family_counts,
        :selector_family,
        :unsupported,
    ) == 1
    @test _mixed_consumer_smoke_count(
        block_set_view.skipped_blocker_counts,
        :blocker,
        :missing_white_lindsey_unit_pair,
    ) == 1
    @test _mixed_consumer_smoke_count(
        block_set_view.skipped_blocker_counts,
        :blocker,
        :unsupported_pair_block_materialization_path,
    ) == 1
    @test !block_set_view.term_batch_results_stored_in_view
    @test !block_set_view.batch_result_objects_stored_in_view
    @test !block_set_view.matrix_fields_stored_in_view
    @test !block_set_view.nested_preflight_or_block_set_summaries_stored_in_view
    @test !block_set_view.term_batch_results_stored_in_summary
    @test !block_set_view.factors_constructed
    @test block_set_view.source_operator_blocks_materialized
    @test block_set_view.final_pair_blocks_materialized
    @test !block_set_view.operator_blocks_materialized
    @test !block_set_view.hamiltonian_data_materialized
    @test !block_set_view.artifacts_materialized
    @test !block_set_view.global_operator_blocks_materialized
    @test !block_set_view.global_hamiltonian_data_materialized
    @test !block_set_view.global_artifacts_materialized
    @test !block_set_view.route_driver_wiring
    @test !block_set_view.coulomb_materialized
    @test !block_set_view.ida_mwg_data_materialized
    @test !block_set_view.pqs_lowdin_materialized
    @test !block_set_view.full_white_lindsey_route_assembled
    @test _mixed_consumer_smoke_compact_summary(block_set_view)

    local_collection =
        CPBMSmoke._one_body_local_block_collection(block_set_consumption)
    @test local_collection.object_kind ==
          :cartesian_pair_block_local_one_body_block_collection
    @test local_collection.status ==
          :partially_materialized_local_one_body_block_collection
    @test isnothing(local_collection.blocker)
    @test local_collection.terms == (:overlap, :kinetic)
    @test local_collection.requested_materialize_terms == (:overlap,)
    @test local_collection.materialized_terms == (:overlap,)
    @test local_collection.deferred_terms == (:kinetic,)
    @test local_collection.entry_count == 4
    @test local_collection.materialized_entry_count == 2
    @test local_collection.skipped_entry_count == 2
    @test local_collection.deferred_term_count == 1
    @test local_collection.source_space_entry_count == 1
    @test local_collection.final_local_entry_count == 1
    @test local_collection.term_separated_entries
    @test local_collection.pair_separated_entries
    @test !local_collection.block_set_results_summed
    @test !local_collection.block_matrices_copied_into_collection
    @test local_collection.source_operator_blocks_materialized
    @test local_collection.final_pair_blocks_materialized
    @test !local_collection.operator_blocks_materialized
    @test !local_collection.hamiltonian_data_materialized
    @test !local_collection.artifacts_materialized
    @test !local_collection.global_operator_blocks_materialized
    @test !local_collection.global_hamiltonian_data_materialized
    @test !local_collection.global_artifacts_materialized
    @test !local_collection.local_operator_assembled
    @test !local_collection.global_operator_assembled
    @test !local_collection.route_driver_wiring
    @test !local_collection.coulomb_materialized
    @test !local_collection.density_density_materialized
    @test !local_collection.ida_mwg_data_materialized
    @test !local_collection.pqs_lowdin_materialized
    @test !local_collection.pqs_shell_projection_materialized
    @test !local_collection.full_white_lindsey_route_assembled

    @test count(
        entry -> entry.selector_family === :direct_direct,
        local_collection.materialized_entries,
    ) == 1
    @test count(
        entry -> entry.selector_family === :pqs_source_pair,
        local_collection.materialized_entries,
    ) == 1
    @test count(
        entry -> entry.selector_family === :white_lindsey_boundary_stratum,
        local_collection.skipped_entries,
    ) == 1
    @test count(
        entry -> entry.selector_family === :unsupported,
        local_collection.skipped_entries,
    ) == 1

    direct_collection_entry = only(
        entry for entry in local_collection.materialized_entries
        if entry.selector_family === :direct_direct
    )
    pqs_collection_entry = only(
        entry for entry in local_collection.materialized_entries
        if entry.selector_family === :pqs_source_pair
    )
    lw_collection_entry = only(
        entry for entry in local_collection.skipped_entries
        if entry.selector_family === :white_lindsey_boundary_stratum
    )
    unsupported_collection_entry = only(
        entry for entry in local_collection.skipped_entries
        if entry.selector_family === :unsupported
    )
    @test direct_collection_entry.block_set_term === :overlap
    @test direct_collection_entry.result_term === :overlap
    @test isnothing(direct_collection_entry.source_space_term)
    @test direct_collection_entry.block_space === :final_local_space
    @test direct_collection_entry.result_available
    @test direct_collection_entry.result.term === :overlap
    @test pqs_collection_entry.block_set_term === :overlap
    @test pqs_collection_entry.result_term === :source_overlap
    @test pqs_collection_entry.source_space_term === :source_overlap
    @test pqs_collection_entry.block_space === :source_space
    @test pqs_collection_entry.result_available
    @test pqs_collection_entry.result.term === :source_overlap
    @test lw_collection_entry.block_set_term === :overlap
    @test isnothing(lw_collection_entry.result_term)
    @test lw_collection_entry.blocker === :missing_white_lindsey_unit_pair
    @test lw_collection_entry.skipped_record_available
    @test lw_collection_entry.skipped_record.blocker ===
          :missing_white_lindsey_unit_pair
    @test unsupported_collection_entry.block_set_term === :overlap
    @test isnothing(unsupported_collection_entry.result_term)
    @test unsupported_collection_entry.blocker ===
          :unsupported_pair_block_materialization_path
    @test unsupported_collection_entry.skipped_record_available
    @test unsupported_collection_entry.skipped_record.blocker ===
          :unsupported_pair_block_materialization_path

    overlap_collection_entries =
        CPBMSmoke._one_body_local_block_collection_entries_for_term(
            local_collection,
            :overlap,
        )
    overlap_collection_materialized =
        CPBMSmoke._one_body_local_block_collection_materialized_entries_for_term(
            local_collection,
            :overlap,
        )
    overlap_collection_skipped =
        CPBMSmoke._one_body_local_block_collection_skipped_entries_for_term(
            local_collection,
            :overlap,
        )
    overlap_collection_status =
        CPBMSmoke._one_body_local_block_collection_term_status(
            local_collection,
            :overlap,
        )
    @test length(overlap_collection_entries) == 4
    @test length(overlap_collection_materialized) == 2
    @test length(overlap_collection_skipped) == 2
    @test count(
        entry -> entry.selector_family === :pqs_source_pair,
        overlap_collection_materialized,
    ) == 1
    @test count(
        entry -> entry.result_term === :source_overlap,
        overlap_collection_materialized,
    ) == 1
    @test CPBMSmoke._one_body_local_block_collection_entries_for_term(
        local_collection,
        :source_overlap,
    ) == ()
    @test overlap_collection_status.status ==
          :partially_materialized_local_one_body_collection_term
    @test overlap_collection_status.blocker === :missing_white_lindsey_unit_pair
    @test overlap_collection_status.requested
    @test overlap_collection_status.requested_materialization
    @test overlap_collection_status.materialized_term
    @test !overlap_collection_status.deferred_term
    @test overlap_collection_status.entry_count == 4
    @test overlap_collection_status.materialized_entry_count == 2
    @test overlap_collection_status.skipped_entry_count == 2
    @test _mixed_consumer_smoke_count(
        overlap_collection_status.materialized_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 1
    @test _mixed_consumer_smoke_count(
        overlap_collection_status.materialized_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 1
    @test _mixed_consumer_smoke_count(
        overlap_collection_status.skipped_selector_family_counts,
        :selector_family,
        :white_lindsey_boundary_stratum,
    ) == 1
    @test _mixed_consumer_smoke_count(
        overlap_collection_status.skipped_blocker_counts,
        :blocker,
        :missing_white_lindsey_unit_pair,
    ) == 1
    @test _mixed_consumer_smoke_count(
        overlap_collection_status.block_space_counts,
        :block_space,
        :final_local_space,
    ) == 1
    @test _mixed_consumer_smoke_count(
        overlap_collection_status.block_space_counts,
        :block_space,
        :source_space,
    ) == 1
    @test _mixed_consumer_smoke_count(
        overlap_collection_status.block_space_counts,
        :block_space,
        :not_materialized,
    ) == 2
    @test !overlap_collection_status.entries_stored_in_status
    @test !overlap_collection_status.matrix_fields_stored_in_status
    @test !overlap_collection_status.block_set_results_summed
    @test !overlap_collection_status.block_matrices_copied_into_status
    @test !overlap_collection_status.local_operator_assembled
    @test !overlap_collection_status.global_operator_assembled
    @test !overlap_collection_status.route_driver_wiring
    @test !overlap_collection_status.operator_blocks_materialized
    @test !overlap_collection_status.hamiltonian_data_materialized
    @test !overlap_collection_status.artifacts_materialized
    @test !overlap_collection_status.coulomb_materialized
    @test !overlap_collection_status.density_density_materialized
    @test !overlap_collection_status.ida_mwg_data_materialized
    @test !overlap_collection_status.pqs_lowdin_materialized
    @test !overlap_collection_status.pqs_shell_projection_materialized
    @test !overlap_collection_status.full_white_lindsey_route_assembled

    kinetic_collection_status =
        CPBMSmoke._one_body_local_block_collection_term_status(
            local_collection,
            :kinetic,
        )
    @test CPBMSmoke._one_body_local_block_collection_entries_for_term(
        local_collection,
        :kinetic,
    ) == ()
    @test CPBMSmoke._one_body_local_block_collection_materialized_entries_for_term(
        local_collection,
        :kinetic,
    ) == ()
    @test CPBMSmoke._one_body_local_block_collection_skipped_entries_for_term(
        local_collection,
        :kinetic,
    ) == ()
    @test kinetic_collection_status.status ==
          :deferred_metadata_only_local_one_body_collection_term
    @test kinetic_collection_status.requested
    @test !kinetic_collection_status.requested_materialization
    @test !kinetic_collection_status.materialized_term
    @test kinetic_collection_status.deferred_term
    @test kinetic_collection_status.entry_count == 0
    @test !kinetic_collection_status.route_driver_wiring
    @test !kinetic_collection_status.hamiltonian_data_materialized
    @test !kinetic_collection_status.coulomb_materialized
    @test !kinetic_collection_status.ida_mwg_data_materialized
    @test !kinetic_collection_status.pqs_lowdin_materialized
    @test !kinetic_collection_status.pqs_shell_projection_materialized
    @test !kinetic_collection_status.full_white_lindsey_route_assembled

    position_collection_status =
        CPBMSmoke._one_body_local_block_collection_term_status(
            local_collection,
            :position_x,
        )
    @test CPBMSmoke._one_body_local_block_collection_entries_for_term(
        local_collection,
        :position_x,
    ) == ()
    @test position_collection_status.status === :term_not_requested
    @test position_collection_status.blocker === :term_not_requested
    @test !position_collection_status.requested
    @test !position_collection_status.requested_materialization
    @test position_collection_status.entry_count == 0
    @test !position_collection_status.route_driver_wiring
    @test !position_collection_status.hamiltonian_data_materialized
    @test !position_collection_status.coulomb_materialized
    @test !position_collection_status.ida_mwg_data_materialized
    @test !position_collection_status.pqs_lowdin_materialized
    @test !position_collection_status.pqs_shell_projection_materialized
    @test !position_collection_status.full_white_lindsey_route_assembled

    direct_collection_pair_entries =
        CPBMSmoke._one_body_local_block_collection_entries_for_pair(
            local_collection,
            (:direct_left, :direct_right),
        )
    direct_collection_pair_materialized =
        CPBMSmoke._one_body_local_block_collection_materialized_entries_for_pair(
            local_collection,
            (:direct_left, :direct_right),
        )
    direct_collection_pair_skipped =
        CPBMSmoke._one_body_local_block_collection_skipped_entries_for_pair(
            local_collection,
            (:direct_left, :direct_right),
        )
    direct_collection_pair_status =
        CPBMSmoke._one_body_local_block_collection_pair_status(
            local_collection,
            (:direct_left, :direct_right),
        )
    @test length(direct_collection_pair_entries) == 1
    @test length(direct_collection_pair_materialized) == 1
    @test direct_collection_pair_skipped == ()
    @test only(direct_collection_pair_materialized).block_space ===
          :final_local_space
    @test direct_collection_pair_status.status ===
          :materialized_local_one_body_collection_pair
    @test isnothing(direct_collection_pair_status.blocker)
    @test direct_collection_pair_status.block_set_terms == (:overlap,)
    @test direct_collection_pair_status.result_terms == (:overlap,)
    @test direct_collection_pair_status.source_space_terms == ()
    @test direct_collection_pair_status.materialized_entry_count == 1
    @test direct_collection_pair_status.skipped_entry_count == 0
    @test direct_collection_pair_status.final_pair_blocks_materialized
    @test direct_collection_pair_status.source_operator_blocks_materialized
    @test !direct_collection_pair_status.operator_blocks_materialized
    @test !direct_collection_pair_status.hamiltonian_data_materialized
    @test !direct_collection_pair_status.artifacts_materialized
    @test !direct_collection_pair_status.route_driver_wiring
    @test !direct_collection_pair_status.coulomb_materialized
    @test !direct_collection_pair_status.ida_mwg_data_materialized
    @test !direct_collection_pair_status.pqs_lowdin_materialized
    @test !direct_collection_pair_status.pqs_shell_projection_materialized
    @test !direct_collection_pair_status.full_white_lindsey_route_assembled

    pqs_collection_pair_entries =
        CPBMSmoke._one_body_local_block_collection_entries_for_pair(
            local_collection,
            (:pqs_left, :pqs_right),
        )
    pqs_collection_pair_status =
        CPBMSmoke._one_body_local_block_collection_pair_status(
            local_collection,
            (:pqs_left, :pqs_right),
        )
    @test length(pqs_collection_pair_entries) == 1
    @test CPBMSmoke._one_body_local_block_collection_skipped_entries_for_pair(
        local_collection,
        (:pqs_left, :pqs_right),
    ) == ()
    @test only(pqs_collection_pair_entries).block_set_term === :overlap
    @test only(pqs_collection_pair_entries).result_term === :source_overlap
    @test only(pqs_collection_pair_entries).block_space === :source_space
    @test pqs_collection_pair_status.status ===
          :materialized_local_one_body_collection_pair
    @test pqs_collection_pair_status.block_set_terms == (:overlap,)
    @test pqs_collection_pair_status.result_terms == (:source_overlap,)
    @test pqs_collection_pair_status.source_space_terms == (:source_overlap,)
    @test pqs_collection_pair_status.source_operator_blocks_materialized
    @test !pqs_collection_pair_status.final_pair_blocks_materialized
    @test !pqs_collection_pair_status.operator_blocks_materialized
    @test !pqs_collection_pair_status.hamiltonian_data_materialized
    @test !pqs_collection_pair_status.artifacts_materialized
    @test !pqs_collection_pair_status.route_driver_wiring
    @test !pqs_collection_pair_status.coulomb_materialized
    @test !pqs_collection_pair_status.ida_mwg_data_materialized
    @test !pqs_collection_pair_status.pqs_lowdin_materialized
    @test !pqs_collection_pair_status.pqs_shell_projection_materialized
    @test !pqs_collection_pair_status.full_white_lindsey_route_assembled

    lw_collection_pair_entries =
        CPBMSmoke._one_body_local_block_collection_entries_for_pair(
            local_collection,
            (:lw_left, :lw_right),
        )
    lw_collection_pair_status =
        CPBMSmoke._one_body_local_block_collection_pair_status(
            local_collection,
            (:lw_left, :lw_right),
        )
    @test length(lw_collection_pair_entries) == 1
    @test CPBMSmoke._one_body_local_block_collection_materialized_entries_for_pair(
        local_collection,
        (:lw_left, :lw_right),
    ) == ()
    @test only(lw_collection_pair_entries).blocker ===
          :missing_white_lindsey_unit_pair
    @test lw_collection_pair_status.status ===
          :skipped_local_one_body_collection_pair
    @test lw_collection_pair_status.blocker === :missing_white_lindsey_unit_pair
    @test lw_collection_pair_status.block_set_terms == (:overlap,)
    @test lw_collection_pair_status.result_terms == ()
    @test lw_collection_pair_status.source_space_terms == ()
    @test lw_collection_pair_status.skipped_entry_count == 1
    @test !lw_collection_pair_status.hamiltonian_data_materialized
    @test !lw_collection_pair_status.coulomb_materialized
    @test !lw_collection_pair_status.ida_mwg_data_materialized
    @test !lw_collection_pair_status.pqs_lowdin_materialized
    @test !lw_collection_pair_status.full_white_lindsey_route_assembled

    unsupported_collection_pair_status =
        CPBMSmoke._one_body_local_block_collection_pair_status(
            local_collection,
            (:unsupported_left, :unsupported_right),
        )
    @test unsupported_collection_pair_status.status ===
          :skipped_local_one_body_collection_pair
    @test unsupported_collection_pair_status.blocker ===
          :unsupported_pair_block_materialization_path
    @test unsupported_collection_pair_status.block_set_terms == (:overlap,)
    @test unsupported_collection_pair_status.skipped_entry_count == 1
    @test !unsupported_collection_pair_status.hamiltonian_data_materialized
    @test !unsupported_collection_pair_status.coulomb_materialized
    @test !unsupported_collection_pair_status.ida_mwg_data_materialized
    @test !unsupported_collection_pair_status.pqs_lowdin_materialized
    @test !unsupported_collection_pair_status.full_white_lindsey_route_assembled

    absent_collection_pair_status =
        CPBMSmoke._one_body_local_block_collection_pair_status(
            local_collection,
            (:missing_left, :missing_right),
        )
    @test CPBMSmoke._one_body_local_block_collection_entries_for_pair(
        local_collection,
        (:missing_left, :missing_right),
    ) == ()
    @test absent_collection_pair_status.status === :pair_key_not_found
    @test absent_collection_pair_status.blocker === :pair_key_not_found
    @test absent_collection_pair_status.block_set_terms == ()
    @test absent_collection_pair_status.result_terms == ()
    @test absent_collection_pair_status.source_space_terms == ()
    @test absent_collection_pair_status.entry_count == 0
    @test !absent_collection_pair_status.hamiltonian_data_materialized
    @test !absent_collection_pair_status.coulomb_materialized
    @test !absent_collection_pair_status.ida_mwg_data_materialized
    @test !absent_collection_pair_status.pqs_lowdin_materialized
    @test !absent_collection_pair_status.full_white_lindsey_route_assembled

    overlap_results =
        CPBMSmoke._one_body_pair_block_results_for_term(
            block_set_consumption,
            :overlap,
        )
    overlap_skips =
        CPBMSmoke._one_body_pair_block_skips_for_term(
            block_set_consumption,
            :overlap,
        )
    overlap_status =
        CPBMSmoke._one_body_pair_block_term_status(
            block_set_consumption,
            :overlap,
        )
    @test length(overlap_results) == 2
    @test length(overlap_skips) == 2
    @test count(
        result -> result.metadata.selector_family === :direct_direct,
        overlap_results,
    ) == 1
    @test count(
        result -> result.metadata.selector_family === :pqs_source_pair,
        overlap_results,
    ) == 1
    @test _mixed_consumer_smoke_count(
        overlap_skips,
        :selector_family,
        :white_lindsey_boundary_stratum,
    ) == 1
    @test _mixed_consumer_smoke_count(
        overlap_skips,
        :selector_family,
        :unsupported,
    ) == 1
    @test _mixed_consumer_smoke_count(
        overlap_skips,
        :blocker,
        :missing_white_lindsey_unit_pair,
    ) == 1
    @test _mixed_consumer_smoke_count(
        overlap_skips,
        :blocker,
        :unsupported_pair_block_materialization_path,
    ) == 1
    @test overlap_status.status ==
          :partially_materialized_mixed_one_body_pair_block_batch
    @test overlap_status.requested
    @test overlap_status.requested_materialization
    @test overlap_status.batch_result_supplied
    @test overlap_status.materialized_count == 2
    @test overlap_status.skipped_count == 2
    @test !overlap_status.term_batch_results_stored_in_status
    @test !overlap_status.batch_result_objects_stored_in_status
    @test !overlap_status.matrix_fields_stored_in_status
    @test _mixed_consumer_smoke_compact_summary(overlap_status)

    kinetic_status =
        CPBMSmoke._one_body_pair_block_term_status(
            block_set_consumption,
            :kinetic,
        )
    @test CPBMSmoke._one_body_pair_block_results_for_term(
        block_set_consumption,
        :kinetic,
    ) == ()
    @test CPBMSmoke._one_body_pair_block_skips_for_term(
        block_set_consumption,
        :kinetic,
    ) == ()
    @test kinetic_status.status ==
          :deferred_metadata_only_mixed_one_body_pair_block_batch
    @test kinetic_status.requested
    @test !kinetic_status.requested_materialization
    @test !kinetic_status.batch_result_supplied
    @test kinetic_status.materialized_count == 0
    @test kinetic_status.skipped_count == 0
    @test _mixed_consumer_smoke_compact_summary(kinetic_status)

    position_status =
        CPBMSmoke._one_body_pair_block_term_status(
            block_set_consumption,
            :position_x,
        )
    @test CPBMSmoke._one_body_pair_block_results_for_term(
        block_set_consumption,
        :position_x,
    ) == ()
    @test CPBMSmoke._one_body_pair_block_skips_for_term(
        block_set_consumption,
        :position_x,
    ) == ()
    @test position_status.status == :term_not_requested
    @test position_status.blocker == :term_not_requested
    @test !position_status.requested
    @test !position_status.requested_materialization
    @test position_status.materialized_count == 0
    @test position_status.skipped_count == 0
    @test _mixed_consumer_smoke_compact_summary(position_status)

    direct_pair_results =
        CPBMSmoke._one_body_pair_block_results_for_pair(
            block_set_consumption,
            (:direct_left, :direct_right),
        )
    direct_pair_skips =
        CPBMSmoke._one_body_pair_block_skips_for_pair(
            block_set_consumption,
            (:direct_left, :direct_right),
        )
    direct_pair_status =
        CPBMSmoke._one_body_pair_block_pair_status(
            block_set_consumption,
            (:direct_left, :direct_right),
        )
    @test length(direct_pair_results) == 1
    @test direct_pair_skips == ()
    @test only(direct_pair_results).term === :overlap
    @test only(direct_pair_results).metadata.selector_family === :direct_direct
    @test direct_pair_status.status == :materialized_mixed_one_body_pair_key
    @test isnothing(direct_pair_status.blocker)
    @test direct_pair_status.materialized_count == 1
    @test direct_pair_status.skipped_count == 0
    @test direct_pair_status.materialized_terms_for_pair == (:overlap,)
    @test direct_pair_status.skipped_terms_for_pair == ()
    @test direct_pair_status.final_pair_blocks_materialized
    @test !direct_pair_status.operator_blocks_materialized
    @test _mixed_consumer_smoke_compact_summary(direct_pair_status)

    pqs_pair_results =
        CPBMSmoke._one_body_pair_block_results_for_pair(
            block_set_consumption,
            (:pqs_left, :pqs_right),
        )
    pqs_pair_status =
        CPBMSmoke._one_body_pair_block_pair_status(
            block_set_consumption,
            (:pqs_left, :pqs_right),
        )
    @test length(pqs_pair_results) == 1
    @test CPBMSmoke._one_body_pair_block_skips_for_pair(
        block_set_consumption,
        (:pqs_left, :pqs_right),
    ) == ()
    @test only(pqs_pair_results).term === :source_overlap
    @test only(pqs_pair_results).metadata.selector_family === :pqs_source_pair
    @test pqs_pair_status.status == :materialized_mixed_one_body_pair_key
    @test pqs_pair_status.source_operator_blocks_materialized
    @test !pqs_pair_status.final_pair_blocks_materialized
    @test pqs_pair_status.materialized_terms_for_pair == (:overlap,)
    @test _mixed_consumer_smoke_compact_summary(pqs_pair_status)

    lw_pair_skips =
        CPBMSmoke._one_body_pair_block_skips_for_pair(
            block_set_consumption,
            (:lw_left, :lw_right),
        )
    lw_pair_status =
        CPBMSmoke._one_body_pair_block_pair_status(
            block_set_consumption,
            (:lw_left, :lw_right),
        )
    @test CPBMSmoke._one_body_pair_block_results_for_pair(
        block_set_consumption,
        (:lw_left, :lw_right),
    ) == ()
    @test length(lw_pair_skips) == 1
    @test only(lw_pair_skips).blocker === :missing_white_lindsey_unit_pair
    @test lw_pair_status.status == :skipped_mixed_one_body_pair_key
    @test lw_pair_status.blocker == :missing_white_lindsey_unit_pair
    @test lw_pair_status.materialized_count == 0
    @test lw_pair_status.skipped_count == 1
    @test lw_pair_status.skipped_terms_for_pair == (:overlap,)
    @test _mixed_consumer_smoke_count(
        lw_pair_status.skipped_blocker_counts,
        :blocker,
        :missing_white_lindsey_unit_pair,
    ) == 1
    @test _mixed_consumer_smoke_compact_summary(lw_pair_status)

    unsupported_pair_skips =
        CPBMSmoke._one_body_pair_block_skips_for_pair(
            block_set_consumption,
            (:unsupported_left, :unsupported_right),
        )
    unsupported_pair_status =
        CPBMSmoke._one_body_pair_block_pair_status(
            block_set_consumption,
            (:unsupported_left, :unsupported_right),
        )
    @test CPBMSmoke._one_body_pair_block_results_for_pair(
        block_set_consumption,
        (:unsupported_left, :unsupported_right),
    ) == ()
    @test length(unsupported_pair_skips) == 1
    @test only(unsupported_pair_skips).blocker ===
          :unsupported_pair_block_materialization_path
    @test unsupported_pair_status.status == :skipped_mixed_one_body_pair_key
    @test unsupported_pair_status.blocker ==
          :unsupported_pair_block_materialization_path
    @test unsupported_pair_status.skipped_count == 1
    @test _mixed_consumer_smoke_compact_summary(unsupported_pair_status)

    absent_pair_status =
        CPBMSmoke._one_body_pair_block_pair_status(
            block_set_consumption,
            (:missing_left, :missing_right),
        )
    @test CPBMSmoke._one_body_pair_block_results_for_pair(
        block_set_consumption,
        (:missing_left, :missing_right),
    ) == ()
    @test CPBMSmoke._one_body_pair_block_skips_for_pair(
        block_set_consumption,
        (:missing_left, :missing_right),
    ) == ()
    @test absent_pair_status.status == :pair_key_not_found
    @test absent_pair_status.blocker == :pair_key_not_found
    @test absent_pair_status.materialized_count == 0
    @test absent_pair_status.skipped_count == 0
    @test !absent_pair_status.matrix_fields_stored_in_status
    @test _mixed_consumer_smoke_compact_summary(absent_pair_status)

    direct_lookup =
        CPBMSmoke._one_body_pair_block_lookup(
            block_set_consumption,
            :overlap,
            (:direct_left, :direct_right),
        )
    @test direct_lookup.status ==
          :materialized_mixed_one_body_pair_block_lookup
    @test direct_lookup.result_available
    @test !direct_lookup.skipped_record_available
    @test direct_lookup.selector_family === :direct_direct
    @test direct_lookup.result.term === :overlap
    @test direct_lookup.result.metadata.selector_family === :direct_direct
    @test direct_lookup.final_pair_blocks_materialized
    @test direct_lookup.explicit_result_field_contains_matrix
    @test !direct_lookup.hidden_matrix_fields_stored_in_lookup
    @test _mixed_consumer_smoke_compact_summary(direct_lookup)

    pqs_lookup =
        CPBMSmoke._one_body_pair_block_lookup(
            block_set_consumption,
            :overlap,
            (:pqs_left, :pqs_right),
        )
    @test pqs_lookup.status ==
          :materialized_mixed_one_body_pair_block_lookup
    @test pqs_lookup.result_available
    @test pqs_lookup.selector_family === :pqs_source_pair
    @test pqs_lookup.result.term === :source_overlap
    @test pqs_lookup.source_operator_blocks_materialized
    @test !pqs_lookup.final_pair_blocks_materialized
    @test pqs_lookup.explicit_result_field_contains_matrix
    @test !pqs_lookup.hidden_matrix_fields_stored_in_lookup
    @test _mixed_consumer_smoke_compact_summary(pqs_lookup)

    lw_lookup =
        CPBMSmoke._one_body_pair_block_lookup(
            block_set_consumption,
            :overlap,
            (:lw_left, :lw_right),
        )
    @test lw_lookup.status == :skipped_mixed_one_body_pair_block_lookup
    @test !lw_lookup.result_available
    @test lw_lookup.skipped_record_available
    @test lw_lookup.selector_family === :white_lindsey_boundary_stratum
    @test lw_lookup.blocker == :missing_white_lindsey_unit_pair
    @test lw_lookup.skipped_record.blocker === :missing_white_lindsey_unit_pair
    @test !lw_lookup.explicit_result_field_contains_matrix
    @test !lw_lookup.hidden_matrix_fields_stored_in_lookup
    @test _mixed_consumer_smoke_compact_summary(lw_lookup)

    deferred_lookup =
        CPBMSmoke._one_body_pair_block_lookup(
            block_set_consumption,
            :kinetic,
            (:direct_left, :direct_right),
        )
    @test deferred_lookup.status ==
          :deferred_metadata_only_mixed_one_body_pair_block_lookup
    @test deferred_lookup.requested
    @test !deferred_lookup.requested_materialization
    @test !deferred_lookup.result_available
    @test !deferred_lookup.skipped_record_available
    @test isnothing(deferred_lookup.blocker)
    @test _mixed_consumer_smoke_compact_summary(deferred_lookup)

    unrequested_lookup =
        CPBMSmoke._one_body_pair_block_lookup(
            block_set_consumption,
            :position_x,
            (:direct_left, :direct_right),
        )
    @test unrequested_lookup.status == :term_not_requested
    @test unrequested_lookup.blocker == :term_not_requested
    @test !unrequested_lookup.requested
    @test !unrequested_lookup.result_available
    @test !unrequested_lookup.skipped_record_available
    @test _mixed_consumer_smoke_compact_summary(unrequested_lookup)

    absent_lookup =
        CPBMSmoke._one_body_pair_block_lookup(
            block_set_consumption,
            :overlap,
            (:missing_left, :missing_right),
        )
    @test absent_lookup.status == :pair_key_not_found
    @test absent_lookup.blocker == :pair_key_not_found
    @test absent_lookup.requested
    @test absent_lookup.requested_materialization
    @test !absent_lookup.result_available
    @test !absent_lookup.skipped_record_available
    @test !absent_lookup.hidden_matrix_fields_stored_in_lookup
    @test _mixed_consumer_smoke_compact_summary(absent_lookup)

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

    direct_entry =
        CPBMSmoke._one_body_local_block_collection_entry(direct_result)
    @test direct_entry.object_kind ==
          :cartesian_pair_block_local_one_body_block_collection_entry
    @test direct_entry.entry_kind === :materialized_result
    @test direct_entry.status ===
          :materialized_local_one_body_block_collection_entry
    @test direct_entry.term === :overlap
    @test direct_entry.pair_key == (:direct_left, :direct_right)
    @test direct_entry.selector_family === :direct_direct
    @test direct_entry.block_space === :final_local_space
    @test direct_entry.block_shape == (2, 2)
    @test direct_entry.block_size == 4
    @test direct_entry.result_available
    @test !direct_entry.skipped_record_available
    @test direct_entry.result.term === :overlap
    @test isnothing(direct_entry.skipped_record)
    @test direct_entry.final_pair_blocks_materialized
    @test !direct_entry.operator_blocks_materialized
    @test !direct_entry.hamiltonian_data_materialized
    @test !direct_entry.artifacts_materialized
    @test !direct_entry.block_copied_into_entry
    @test !direct_entry.local_operator_assembled
    @test !direct_entry.global_operator_assembled
    @test !direct_entry.route_driver_wiring
    @test !direct_entry.coulomb_materialized
    @test !direct_entry.ida_mwg_data_materialized
    @test !direct_entry.pqs_lowdin_materialized
    @test !direct_entry.pqs_shell_projection_materialized
    @test !direct_entry.full_white_lindsey_route_assembled

    pqs_entry = CPBMSmoke._one_body_local_block_collection_entry(pqs_result)
    @test pqs_entry.entry_kind === :materialized_result
    @test pqs_entry.term === :source_overlap
    @test pqs_entry.pair_key == (:pqs_left, :pqs_right)
    @test pqs_entry.selector_family === :pqs_source_pair
    @test pqs_entry.block_space === :source_space
    @test pqs_entry.block_shape == (8, 8)
    @test pqs_entry.block_size == 64
    @test pqs_entry.result_available
    @test !pqs_entry.skipped_record_available
    @test pqs_entry.result.term === :source_overlap
    @test pqs_entry.source_operator_blocks_materialized
    @test !pqs_entry.final_pair_blocks_materialized
    @test !pqs_entry.operator_blocks_materialized
    @test !pqs_entry.hamiltonian_data_materialized
    @test !pqs_entry.artifacts_materialized
    @test !pqs_entry.block_copied_into_entry
    @test !pqs_entry.route_driver_wiring
    @test !pqs_entry.coulomb_materialized
    @test !pqs_entry.ida_mwg_data_materialized
    @test !pqs_entry.pqs_lowdin_materialized
    @test !pqs_entry.pqs_shell_projection_materialized
    @test !pqs_entry.full_white_lindsey_route_assembled

    lw_skip = only(
        skip for skip in consumption.batch_result.skipped_records
        if skip.selector_family === :white_lindsey_boundary_stratum
    )
    unsupported_skip = only(
        skip for skip in consumption.batch_result.skipped_records
        if skip.selector_family === :unsupported
    )
    lw_entry =
        CPBMSmoke._one_body_local_block_collection_skipped_entry(lw_skip)
    @test lw_entry.entry_kind === :skipped_record
    @test lw_entry.status === :skipped_local_one_body_block_collection_entry
    @test lw_entry.term === :overlap
    @test lw_entry.pair_key == (:lw_left, :lw_right)
    @test lw_entry.selector_family === :white_lindsey_boundary_stratum
    @test lw_entry.blocker === :missing_white_lindsey_unit_pair
    @test lw_entry.block_space === :not_materialized
    @test isnothing(lw_entry.block_shape)
    @test lw_entry.block_size == 0
    @test !lw_entry.result_available
    @test lw_entry.skipped_record_available
    @test isnothing(lw_entry.result)
    @test lw_entry.skipped_record.blocker === :missing_white_lindsey_unit_pair
    @test !lw_entry.source_operator_blocks_materialized
    @test !lw_entry.final_pair_blocks_materialized
    @test !lw_entry.operator_blocks_materialized
    @test !lw_entry.hamiltonian_data_materialized
    @test !lw_entry.artifacts_materialized
    @test !lw_entry.route_driver_wiring
    @test !lw_entry.coulomb_materialized
    @test !lw_entry.ida_mwg_data_materialized
    @test !lw_entry.pqs_lowdin_materialized
    @test !lw_entry.pqs_shell_projection_materialized
    @test !lw_entry.full_white_lindsey_route_assembled

    unsupported_entry =
        CPBMSmoke._one_body_local_block_collection_skipped_entry(
            unsupported_skip,
        )
    @test unsupported_entry.entry_kind === :skipped_record
    @test unsupported_entry.term === :overlap
    @test unsupported_entry.pair_key == (:unsupported_left, :unsupported_right)
    @test unsupported_entry.selector_family === :unsupported
    @test unsupported_entry.blocker ===
          :unsupported_pair_block_materialization_path
    @test unsupported_entry.block_space === :not_materialized
    @test !unsupported_entry.result_available
    @test unsupported_entry.skipped_record_available
    @test isnothing(unsupported_entry.result)
    @test unsupported_entry.skipped_record.blocker ===
          :unsupported_pair_block_materialization_path
    @test !unsupported_entry.hamiltonian_data_materialized
    @test !unsupported_entry.artifacts_materialized
    @test !unsupported_entry.route_driver_wiring
    @test !unsupported_entry.coulomb_materialized
    @test !unsupported_entry.ida_mwg_data_materialized
    @test !unsupported_entry.pqs_lowdin_materialized
    @test !unsupported_entry.pqs_shell_projection_materialized
    @test !unsupported_entry.full_white_lindsey_route_assembled
end
