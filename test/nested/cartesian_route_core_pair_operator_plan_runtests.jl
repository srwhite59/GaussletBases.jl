using Test
using GaussletBases

const CRC = GaussletBases.CartesianRouteCore

function _crc_pair_plan_complete_shell_region()
    outer = CRC.filled_cpb(1:5, 1:5, 1:5; role = :outer_box)
    inner = CRC.filled_cpb(2:4, 2:4, 2:4; role = :inner_exclusion_box)
    support = CRC.complete_shell_support(outer, inner)
    return CRC.shellification_region(:synthetic_shell, support)
end

function _crc_pair_plan_pqs_unit(unit_key::Symbol)
    region = _crc_pair_plan_complete_shell_region()
    source = CRC.filled_cpb(1:5, 1:5, 1:5; role = Symbol(unit_key, "_source_cpb"))
    lowering = CRC.pqs_filled_source_lowering(region, source)
    intermediate = CRC.intermediate_retained_space(
        lowering;
        retained_rule = :pqs_boundary_comx_product_modes,
        source_mode_dims = CRC.shape(source),
    )
    realization = CRC.pqs_shell_realization(intermediate, region)
    return CRC.final_retained_unit(
        unit_key,
        :synthetic_pqs_shell,
        lowering,
        intermediate,
        realization;
        dimension = intermediate.dimension,
    )
end

function _crc_pair_plan_direct_unit(unit_key::Symbol)
    source = CRC.slab_cpb(1:5, 1:5, 3:3; role = Symbol(unit_key, "_source_cpb"))
    support = CRC.owned_cpb(source)
    region = CRC.shellification_region(:synthetic_direct_slab, support)
    lowering = CRC.lowering_source(:direct_identity_cpb, region, (source,))
    intermediate = CRC.intermediate_retained_space(
        lowering;
        retained_rule = :identity_source_modes,
        source_mode_dims = CRC.shape(source),
    )
    realization = CRC.trivial_shell_realization(intermediate, region)
    return CRC.final_retained_unit(
        unit_key,
        :synthetic_direct_slab,
        lowering,
        intermediate,
        realization;
        dimension = intermediate.dimension,
    )
end

function _crc_pair_plan_aggregate_unit(unit_key::Symbol)
    source = CRC.filled_cpb(1:3, 1:3, 1:3; role = Symbol(unit_key, "_source_cpb"))
    support = CRC.owned_cpb(source)
    region = CRC.shellification_region(:synthetic_atom_local_child, support)
    lowering = CRC.lowering_source(
        :white_lindsey_atom_local_child_shellification,
        region,
        (source,);
        metadata = (; child_source_cpbs_enumerated = false),
    )
    intermediate = CRC.intermediate_retained_space(
        lowering;
        retained_rule = :atom_local_child_shellification_sequence,
        source_mode_dims = CRC.shape(source),
        materialized = false,
        metadata = (; child_source_cpbs_enumerated = false),
    )
    realization = CRC.trivial_shell_realization(
        intermediate,
        region;
        final_dimension = nothing,
    )
    return CRC.final_retained_unit(
        unit_key,
        :synthetic_atom_local_child,
        lowering,
        intermediate,
        realization,
    )
end

function _crc_pair_plan_by_key(inventory, pair_key)
    return only(plan for plan in CRC.pair_operator_plans(inventory)
                if plan.pair.pair_key == pair_key)
end

function _crc_pair_count(counts, field::Symbol, value)
    matches = Tuple(entry for entry in counts if getproperty(entry, field) === value)
    isempty(matches) && return 0
    return only(matches).pair_count
end

@testset "CartesianRouteCore pair operator plan contract" begin
    pqs_unit = _crc_pair_plan_pqs_unit(:pqs_shell)
    direct_unit = _crc_pair_plan_direct_unit(:direct_slab)
    aggregate_unit = _crc_pair_plan_aggregate_unit(:atom_local_child)

    pair_inventory = CRC.unit_pair_inventory((pqs_unit, direct_unit, aggregate_unit))
    plan_inventory = CRC.pair_operator_plan_inventory(
        pair_inventory;
        supported_terms = (:overlap, :kinetic),
    )

    @test plan_inventory isa CRC.PairOperatorPlanInventory
    @test CRC.pair_operator_plan_count(plan_inventory) == 6
    @test length(CRC.pair_operator_plans(plan_inventory)) ==
          length(CRC.pair_entries(pair_inventory))
    @test CRC.supported_terms(plan_inventory) == (:overlap, :kinetic)
    @test CRC.materialization_status(plan_inventory) ==
          :blocked_metadata_only_not_materialized
    @test CRC.blocker(plan_inventory) ==
          :aggregate_subtree_operator_plan_required
    @test !plan_inventory.materialized
    @test sum(entry.pair_count for entry in
              CRC.pair_operator_source_path_counts(plan_inventory)) == 6
    @test _crc_pair_count(
        CRC.pair_operator_source_path_counts(plan_inventory),
        :source_operator_path,
        :pqs_source_cpb_1d_factor_path,
    ) == 2
    @test _crc_pair_count(
        CRC.pair_operator_source_path_counts(plan_inventory),
        :source_operator_path,
        :direct_identity_cpb_path,
    ) == 1
    @test _crc_pair_count(
        CRC.pair_operator_source_path_counts(plan_inventory),
        :source_operator_path,
        :aggregate_subtree_adapter_required,
    ) == 3
    @test _crc_pair_count(
        CRC.pair_operator_final_block_path_counts(plan_inventory),
        :final_block_path,
        :source_block_realization_then_final_block,
    ) == 2
    @test _crc_pair_count(
        CRC.pair_operator_final_block_path_counts(plan_inventory),
        :final_block_path,
        :source_block_direct_to_final_block,
    ) == 1
    @test _crc_pair_count(
        CRC.pair_operator_final_block_path_counts(plan_inventory),
        :final_block_path,
        :blocked_final_pair_block_path,
    ) == 3
    @test _crc_pair_count(
        CRC.pair_operator_materialization_status_counts(plan_inventory),
        :materialization_status,
        :metadata_only_not_materialized,
    ) == 3
    @test _crc_pair_count(
        CRC.pair_operator_materialization_status_counts(plan_inventory),
        :materialization_status,
        :blocked_metadata_only_not_materialized,
    ) == 3
    @test _crc_pair_count(
        CRC.pair_operator_blocker_counts(plan_inventory),
        :blocker,
        nothing,
    ) == 3
    @test _crc_pair_count(
        CRC.pair_operator_blocker_counts(plan_inventory),
        :blocker,
        :aggregate_subtree_operator_plan_required,
    ) == 3
    synthetic_family_counts = CRC.pair_operator_plan_family_counts(plan_inventory)
    @test sum(entry.pair_count for entry in synthetic_family_counts) == 6
    @test all(entry -> !entry.materialized, synthetic_family_counts)
    blocked_readiness =
        CRC.pair_operator_materialization_readiness(plan_inventory)
    @test !blocked_readiness.ready
    @test blocked_readiness.status == :blocked_pair_operator_materialization
    @test blocked_readiness.blocker ==
          :aggregate_subtree_operator_plan_required
    @test blocked_readiness.requirements ==
          (
              :nonempty_pair_operator_plan_inventory,
              :no_blocked_typed_pair_operator_plans,
              :no_pending_pair_operator_source_paths,
              :no_already_materialized_typed_pair_operator_plans,
          )
    @test blocked_readiness.plan_count == 6
    @test blocked_readiness.blocked_count == 3
    @test blocked_readiness.materialized_count == 0

    direct_ready_inventory = CRC.pair_operator_plan_inventory(
        CRC.unit_pair_inventory((direct_unit,)),
    )
    direct_ready_readiness =
        CRC.pair_operator_materialization_readiness(direct_ready_inventory)
    @test direct_ready_readiness.ready
    @test direct_ready_readiness.status == :ready_pair_operator_materialization
    @test isnothing(direct_ready_readiness.blocker)
    @test direct_ready_readiness.plan_count == 1
    @test direct_ready_readiness.blocked_count == 0
    @test direct_ready_readiness.materialized_count == 0

    direct_materialized_inventory = CRC.pair_operator_plan_inventory(
        CRC.unit_pair_inventory((direct_unit,));
        materialized = true,
    )
    direct_materialized_readiness =
        CRC.pair_operator_materialization_readiness(
            direct_materialized_inventory,
        )
    @test !direct_materialized_readiness.ready
    @test direct_materialized_readiness.status ==
          :blocked_already_materialized_pair_operator_plans
    @test direct_materialized_readiness.blocker ==
          :typed_pair_operator_plans_already_materialized
    @test direct_materialized_readiness.plan_count == 1
    @test direct_materialized_readiness.blocked_count == 0
    @test direct_materialized_readiness.materialized_count == 1

    pqs_plan = _crc_pair_plan_by_key(plan_inventory, (:pqs_shell, :pqs_shell))
    @test pqs_plan isa CRC.PairOperatorPlan
    @test CRC.source_operator_path(pqs_plan) ==
          :pqs_source_cpb_1d_factor_path
    @test CRC.realization_path(pqs_plan) == (;
        left = :shell_projection_lowdin_planned,
        right = :shell_projection_lowdin_planned,
    )
    @test CRC.final_block_path(pqs_plan) ==
          :source_block_realization_then_final_block
    @test CRC.materialization_status(pqs_plan) ==
          :metadata_only_not_materialized
    @test isnothing(CRC.blocker(pqs_plan))
    @test !pqs_plan.source_operator_plan.materialized
    @test !pqs_plan.final_pair_block_plan.materialized

    direct_plan = _crc_pair_plan_by_key(plan_inventory, (:direct_slab, :direct_slab))
    @test CRC.source_operator_path(direct_plan) == :direct_identity_cpb_path
    @test CRC.left_realization_path(direct_plan) ==
          :identity_or_trivial_embedding
    @test CRC.right_realization_path(direct_plan) ==
          :identity_or_trivial_embedding
    @test CRC.final_block_path(direct_plan) ==
          :source_block_direct_to_final_block
    @test CRC.materialization_status(direct_plan) ==
          :metadata_only_not_materialized
    @test isnothing(CRC.blocker(direct_plan))
    @test !direct_plan.source_operator_plan.materialized
    @test !direct_plan.final_pair_block_plan.materialized

    aggregate_plan =
        _crc_pair_plan_by_key(plan_inventory, (:atom_local_child, :atom_local_child))
    @test CRC.source_operator_path(aggregate_plan) ==
          :aggregate_subtree_adapter_required
    @test CRC.source_operator_path(aggregate_plan) !=
          :pqs_source_cpb_1d_factor_path
    @test CRC.source_operator_path(aggregate_plan) !=
          :direct_source_cpb_1d_factor_path
    @test CRC.source_operator_path(aggregate_plan) != :direct_identity_cpb_path
    @test CRC.final_block_path(aggregate_plan) ==
          :blocked_final_pair_block_path
    @test CRC.materialization_status(aggregate_plan) ==
          :blocked_metadata_only_not_materialized
    @test CRC.blocker(aggregate_plan) ==
          :aggregate_subtree_operator_plan_required
    @test !aggregate_plan.source_operator_plan.materialized
    @test !aggregate_plan.final_pair_block_plan.materialized

    mixed_aggregate_plan =
        _crc_pair_plan_by_key(plan_inventory, (:direct_slab, :atom_local_child))
    @test CRC.source_operator_path(mixed_aggregate_plan) ==
          :aggregate_subtree_adapter_required
    @test CRC.blocker(mixed_aggregate_plan) ==
          :aggregate_subtree_operator_plan_required
end
