# Runtime role: tiny smoke / route-shaped adapter contract.
#
# This verifies the private route-local one-body adapter without driver wiring,
# global assembly, Hamiltonians, Coulomb, IDA/MWG, PQS Lowdin/projection, or
# White-Lindsey oracle fixtures.

using Test
using GaussletBases

const CPBMRouteAdapter = GaussletBases.CartesianPairBlockMaterialization
const CPBRouteAdapter = GaussletBases.CartesianCPB
const CTLRouteAdapter = GaussletBases.CartesianTerminalLowering
const CRURouteAdapter = GaussletBases.CartesianRetainedUnits
const CRTCRouteAdapter = GaussletBases.CartesianRetainedUnitTransformContracts
const CUPRouteAdapter = GaussletBases.CartesianUnitPairs
const CPOPRouteAdapter = GaussletBases.CartesianPairOperatorPlans

function _route_adapter_count(counts, field::Symbol, value)
    matches = Tuple(entry for entry in counts if getproperty(entry, field) == value)
    isempty(matches) && return 0
    return sum(entry -> hasproperty(entry, :count) ? entry.count : 1, matches)
end

function _route_adapter_pair_operator_plan()
    lowering_plan = CTLRouteAdapter.TerminalLoweringPlan(
        CTLRouteAdapter.PQSLowering(q = 2),
        (),
        (),
        (; status = :available_terminal_lowering_plan, materialized = false),
        (; fixture = :route_one_body_adapter),
    )
    retained_plan = CRURouteAdapter.RetainedUnitPlan(
        CRURouteAdapter.MetadataOnlyRetainedUnits(),
        lowering_plan,
        (),
        (; status = :available_retained_unit_plan, materialized = false),
        (; fixture = :route_one_body_adapter),
    )
    unit_pair_plan = CUPRouteAdapter.UnitPairPlan(
        CUPRouteAdapter.MetadataOnlyUnitPairs(),
        retained_plan,
        (),
        nothing,
        (; status = :available_unit_pair_plan, materialized = false),
        (; fixture = :route_one_body_adapter),
    )
    transform_plan =
        CRTCRouteAdapter.RetainedUnitTransformContractPlan(
            CRTCRouteAdapter.MetadataOnlyRetainedUnitTransformContracts(),
            retained_plan,
            (),
            (; status = :available_transform_contract_plan, materialized = false),
            (; fixture = :route_one_body_adapter),
        )
    return CPOPRouteAdapter.PairOperatorPlan(
        CPOPRouteAdapter.MetadataOnlyPairOperatorPlans(),
        unit_pair_plan,
        transform_plan,
        (),
        nothing,
        (; status = :available_pair_operator_plan, materialized = false),
        (; fixture = :route_one_body_adapter),
    )
end

function _route_adapter_record(
    pair_key::Tuple{Symbol,Symbol},
    pair_index::Int,
    pair_family::Symbol,
    materialization_path::Symbol;
    metadata = (;),
)
    return CPBMRouteAdapter.PairBlockMaterializationRecord(
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

function _route_adapter_plan()
    direct_left = CPBRouteAdapter.cpb(
        1:1,
        1:2,
        1:1;
        role = :route_adapter_direct_left_source,
    )
    direct_right = CPBRouteAdapter.cpb(
        2:2,
        1:2,
        1:1;
        role = :route_adapter_direct_right_source,
    )
    records = (
        _route_adapter_record(
            (:direct_left, :direct_right),
            1,
            :direct_direct,
            :direct_direct_pair_block_materialization_pilot;
            metadata = (;
                left_source_cpbs = (direct_left,),
                right_source_cpbs = (direct_right,),
            ),
        ),
        _route_adapter_record(
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
        _route_adapter_record(
            (:lw_left, :lw_right),
            3,
            :white_lindsey_boundary_stratum,
            :white_lindsey_boundary_stratum_adapter_preflight,
        ),
        _route_adapter_record(
            (:unsupported_left, :unsupported_right),
            4,
            :unsupported_pair_family,
            :synthetic_unsupported_pair_block_path,
        ),
    )
    return CPBMRouteAdapter.PairBlockMaterializationPlan(
        CPBMRouteAdapter.MetadataOnlyPairBlockMaterialization(),
        _route_adapter_pair_operator_plan(),
        records,
        (; status = :available_pair_block_materialization_plan),
        (; fixture = :route_one_body_adapter),
    )
end

function _route_adapter_overlap_1d()
    return (;
        x = [1.0 0.2; 0.2 1.1],
        y = [1.2 0.3; 0.3 1.3],
        z = [1.4 0.4; 0.4 1.5],
    )
end

@testset "CartesianPairBlockMaterialization route local one-body adapter" begin
    plan = _route_adapter_plan()
    inputs = (;
        parent_axis_counts = (2, 2, 2),
        overlap_1d = _route_adapter_overlap_1d(),
    )

    adapter = CPBMRouteAdapter.route_local_one_body_block_collection(
        plan;
        terms = (:overlap,),
        inputs,
        materialize_terms = (:overlap,),
    )

    @test adapter.object_kind ==
          :cartesian_pair_block_route_local_one_body_block_adapter
    @test adapter.status ==
          :partially_materialized_local_one_body_block_collection
    @test isnothing(adapter.blocker)
    @test adapter.requested_terms == (:overlap,)
    @test adapter.materialized_terms == (:overlap,)
    @test adapter.deferred_terms == ()
    @test adapter.total_materialized_count == 2
    @test adapter.total_skipped_count == 2
    @test adapter.entry_count == adapter.local_block_collection_summary.entry_count
    @test adapter.materialized_entry_count ==
          adapter.local_block_collection_summary.materialized_entry_count
    @test adapter.skipped_entry_count ==
          adapter.local_block_collection_summary.skipped_entry_count
    @test adapter.materialized_entry_count == 2
    @test adapter.skipped_entry_count == 2
    @test adapter.source_space_entry_count == 1
    @test adapter.final_local_entry_count == 1
    @test adapter.term_separated_entries
    @test adapter.pair_separated_entries
    @test adapter.result_terms_remain_separated
    @test !adapter.block_set_results_summed
    @test !adapter.factors_constructed
    @test adapter.source_operator_blocks_materialized
    @test adapter.final_pair_blocks_materialized

    summary = adapter.local_block_collection_summary
    @test _route_adapter_count(
        summary.materialized_selector_family_counts,
        :selector_family,
        :direct_direct,
    ) == 1
    @test _route_adapter_count(
        summary.materialized_selector_family_counts,
        :selector_family,
        :pqs_source_pair,
    ) == 1
    @test _route_adapter_count(
        summary.skipped_selector_family_counts,
        :selector_family,
        :white_lindsey_boundary_stratum,
    ) == 1
    @test _route_adapter_count(
        summary.skipped_selector_family_counts,
        :selector_family,
        :unsupported,
    ) == 1
    @test _route_adapter_count(
        summary.skipped_blocker_counts,
        :blocker,
        :missing_white_lindsey_unit_pair,
    ) == 1
    @test _route_adapter_count(
        summary.skipped_blocker_counts,
        :blocker,
        :unsupported_pair_block_materialization_path,
    ) == 1

    pqs_entry = only(
        entry for entry in adapter.local_block_collection.materialized_entries
        if entry.selector_family === :pqs_source_pair
    )
    direct_entry = only(
        entry for entry in adapter.local_block_collection.materialized_entries
        if entry.selector_family === :direct_direct
    )
    @test pqs_entry.result_term === :source_overlap
    @test pqs_entry.block_space === :source_space
    @test direct_entry.result_term === :overlap
    @test direct_entry.block_space === :final_local_space

    wrapped = (; pair_block_materialization_plan = plan)
    wrapped_adapter = CPBMRouteAdapter.route_local_one_body_block_collection(
        wrapped;
        terms = (:overlap,),
        inputs,
        materialize_terms = (:overlap,),
    )
    @test wrapped_adapter.status == adapter.status
    @test wrapped_adapter.total_materialized_count == adapter.total_materialized_count
    @test wrapped_adapter.total_skipped_count == adapter.total_skipped_count

    pair_operator_adapter = CPBMRouteAdapter.route_local_one_body_block_collection(
        _route_adapter_pair_operator_plan();
        terms = (:overlap,),
        inputs,
        materialize_terms = (:overlap,),
    )
    @test pair_operator_adapter.status == :empty_local_one_body_block_collection
    @test pair_operator_adapter.entry_count == 0
    @test pair_operator_adapter.total_materialized_count == 0
    @test pair_operator_adapter.total_skipped_count == 0

    deferred = CPBMRouteAdapter.route_local_one_body_block_collection(
        plan;
        terms = (:overlap,),
        inputs,
    )
    @test deferred.status ==
          :deferred_metadata_only_local_one_body_block_collection
    @test deferred.materialized_entry_count == 0
    @test deferred.skipped_entry_count == 0
    @test deferred.deferred_terms == (:overlap,)

    @test !adapter.local_operator_assembled
    @test !adapter.global_operator_assembled
    @test !adapter.route_driver_wiring
    @test !adapter.operator_blocks_materialized
    @test !adapter.hamiltonian_data_materialized
    @test !adapter.artifacts_materialized
    @test !adapter.global_operator_blocks_materialized
    @test !adapter.global_hamiltonian_data_materialized
    @test !adapter.global_artifacts_materialized
    @test !adapter.exports_materialized
    @test !adapter.coulomb_materialized
    @test !adapter.density_density_materialized
    @test !adapter.ida_mwg_data_materialized
    @test !adapter.pqs_lowdin_materialized
    @test !adapter.pqs_shell_projection_materialized
    @test !adapter.full_white_lindsey_route_assembled
end
