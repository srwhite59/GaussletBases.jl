# Integration/slow test. Do not include in default nested runner.

using Test
using GaussletBases

function _terminal_geometry_box_count(box)
    return prod(length(box[axis]) for axis in 1:3)
end

function _terminal_geometry_linear_indices_for_box(box, parent_dims)
    inds = Set{Int}()
    nx, ny, _ = parent_dims
    for k in box[3], j in box[2], i in box[1]
        push!(inds, i + (j - 1) * nx + (k - 1) * nx * ny)
    end
    return inds
end

function _terminal_geometry_region_indices(region, parent_dims)
    inds = _terminal_geometry_linear_indices_for_box(region.outer_box, parent_dims)
    if !isnothing(region.inner_exclusion_box)
        setdiff!(
            inds,
            _terminal_geometry_linear_indices_for_box(
                region.inner_exclusion_box,
                parent_dims,
            ),
        )
    end
    return inds
end

function _terminal_geometry_independent_coverage(plan)
    seen = Set{Int}()
    duplicates = 0
    for region in plan.regions
        inds = _terminal_geometry_region_indices(region, plan.parent_dims)
        @test length(inds) == region.support_count
        for idx in inds
            if idx in seen
                duplicates += 1
            else
                push!(seen, idx)
            end
        end
    end
    return (;
        covered = length(seen),
        duplicates,
        expected = prod(plan.parent_dims),
        region_support = sum(region.support_count for region in plan.regions; init = 0),
    )
end

function _terminal_geometry_assert_full_partition(plan)
    coverage = _terminal_geometry_independent_coverage(plan)
    @test coverage.covered == coverage.expected
    @test coverage.duplicates == 0
    @test coverage.region_support == coverage.expected
    @test plan.coverage.coverage_complete
    @test plan.coverage.covered_site_count == coverage.expected
    @test plan.coverage.duplicate_site_count == 0
    @test all(region.terminal for region in plan.regions)
    @test !plan.aggregate_atom_boxes_emitted
    @test plan.diagnostics.private_development_only
    @test plan.diagnostics.shellification_geometry_only
    @test !plan.diagnostics.coordinate_product_box_lowering_materialized
    @test !plan.diagnostics.retained_spaces_materialized
    @test !plan.diagnostics.coefficient_maps_materialized
    @test !plan.diagnostics.operator_blocks_materialized
    @test !plan.diagnostics.pair_operator_blocks_materialized
    @test !plan.diagnostics.hamiltonian_data_materialized
    @test !plan.diagnostics.public_default_behavior_changed
    @test !any(
        endswith(String(region.region_kind), "_cpb")
        for region in plan.regions
    )
    return coverage
end

function _terminal_geometry_assert_private_summary_contract(plan)
    summary =
        GaussletBases._cartesian_terminal_shellification_geometry_private_summary(
            plan,
        )
    @test summary.object_kind ==
          :cartesian_terminal_shellification_geometry_private_summary
    @test summary.source_kind == :terminal_cartesian_shellification_geometry
    @test summary.private_development_only
    @test summary.system_kind == plan.system_kind
    @test summary.parent_dims == plan.parent_dims
    @test summary.parent_box == plan.parent_box
    @test summary.nuclear_indices == plan.nuclear_indices
    @test summary.bond_axis == plan.bond_axis
    @test summary.core_side == plan.core_side
    @test summary.q == plan.q
    @test summary.ordered_terminal_region_roles == plan.region_roles
    @test summary.ordered_terminal_region_kinds == plan.region_kinds
    @test summary.region_count == plan.region_count
    @test summary.terminal_region_count == plan.terminal_region_count
    @test summary.region_support_counts ==
          Tuple(region.support_count for region in plan.regions)
    @test summary.total_support_count == prod(plan.parent_dims)
    @test summary.coverage_status == :coverage_complete
    @test summary.coverage_complete
    @test summary.coverage.expected_parent_site_count == prod(plan.parent_dims)
    @test summary.coverage.region_support_count == prod(plan.parent_dims)
    @test summary.coverage.duplicate_site_count == 0
    @test summary.coverage.missing_site_count == 0
    @test !summary.shellification_regions_are_cpbs
    @test !summary.shellification_regions_are_lowering_sources
    @test !summary.lowering_applied_by_summary
    @test !summary.retained_spaces_materialized
    @test !summary.coefficient_maps_materialized
    @test !summary.operator_blocks_materialized
    @test !summary.pair_operator_blocks_materialized
    @test !summary.hamiltonian_data_materialized
    @test !summary.public_default_behavior_changed
    @test all(!region.shellification_region_is_cpb for region in summary.regions)
    @test all(
        !region.shellification_region_is_lowering_source
        for region in summary.regions
    )
    @test all(
        region.lowering_status == :planned_not_lowered
        for region in summary.regions
    )
    return summary
end

function _terminal_geometry_assert_scaffold_contract(
    plan;
    route_family::Symbol = :white_lindsey_low_order,
)
    scaffold = GaussletBases._cartesian_terminal_shellification_geometry_scaffold(
        plan;
        route_family,
    )
    @test scaffold.object_kind == :cartesian_terminal_shellification_scaffold3d
    @test scaffold.status == :planned_metadata_only
    @test scaffold.private_development_only
    @test scaffold.source_kind == :terminal_cartesian_shellification_geometry
    @test scaffold.route_family == route_family
    @test scaffold.system_classification == plan.system_kind
    @test scaffold.shellification_role == :terminal_cartesian_shellification_geometry
    @test scaffold.shellification_stage == :route_neutral_spatial_planning
    @test scaffold.lowering_stage == :not_lowered_by_shellification_plan
    @test scaffold.parent_box == plan.parent_box
    @test scaffold.working_box == plan.parent_box
    @test scaffold.full_parent_working_box
    @test scaffold.parent_dims == plan.parent_dims
    @test scaffold.nuclear_indices == plan.nuclear_indices
    @test scaffold.bond_axis == plan.bond_axis
    @test scaffold.core_side == plan.core_side
    @test scaffold.q == plan.q
    @test scaffold.region_count == plan.region_count
    @test scaffold.ordered_region_roles == plan.region_roles
    @test scaffold.ordered_region_kinds == plan.region_kinds
    @test scaffold.ordered_region_boxes ==
          Tuple(region.outer_box for region in plan.regions)
    @test scaffold.coverage.object_kind ==
          :cartesian_terminal_shellification_scaffold_coverage3d
    @test scaffold.coverage.coverage_complete
    @test scaffold.coverage.expected_support_count == prod(plan.parent_dims)
    @test scaffold.coverage.region_support_count == prod(plan.parent_dims)
    @test scaffold.coverage.covered_support_count == prod(plan.parent_dims)
    @test scaffold.coverage.duplicate_count == 0
    @test scaffold.coverage.missing_count == 0
    @test scaffold.coverage.outside_count == 0
    @test scaffold.terminal_geometry_summary.object_kind ==
          :cartesian_terminal_shellification_geometry_private_summary
    @test scaffold.diagnostics.terminal_geometry_authority
    @test scaffold.diagnostics.route_neutral_spatial_planning
    @test !scaffold.diagnostics.lowering_applied_by_plan
    @test !scaffold.diagnostics.materialization_behavior_changed
    @test !scaffold.diagnostics.public_default_behavior_changed
    @test !scaffold.diagnostics.shellification_regions_are_cpbs
    @test !scaffold.diagnostics.shellification_regions_are_lowering_sources
    @test !scaffold.diagnostics.retained_spaces_materialized
    @test !scaffold.diagnostics.coefficient_maps_materialized
    @test !scaffold.diagnostics.operator_blocks_materialized
    @test !scaffold.diagnostics.pair_operator_blocks_materialized
    @test !scaffold.diagnostics.hamiltonian_data_materialized
    @test !scaffold.diagnostics.cartesian_shells_behavior_changed
    @test all(region.source_backed == false for region in scaffold.regions)
    @test all(!region.shellification_region_is_cpb for region in scaffold.regions)
    @test all(
        !region.shellification_region_is_lowering_source
        for region in scaffold.regions
    )
    @test all(
        !region.coordinate_product_box_lowering_materialized
        for region in scaffold.regions
    )
    @test all(
        region.lowering_status == :planned_not_lowered
        for region in scaffold.regions
    )
    return scaffold
end

function _terminal_geometry_assert_region_unit_inventory_contract(scaffold)
    inventory =
        GaussletBases._cartesian_terminal_shellification_region_unit_inventory(
            scaffold,
        )
    @test inventory.object_kind == :cartesian_terminal_region_unit_inventory
    @test inventory.status == :available_terminal_region_unit_inventory
    @test inventory.inventory_source == :terminal_shellification_scaffold
    @test inventory.private_development_only
    @test inventory.terminal_region_count == scaffold.region_count
    @test inventory.unit_count == scaffold.region_count
    @test length(inventory.terminal_region_units) == scaffold.region_count
    @test inventory.unit_records == inventory.terminal_region_units
    @test inventory.terminal_region_roles == scaffold.ordered_region_roles
    @test inventory.terminal_region_kinds == scaffold.ordered_region_kinds
    @test inventory.support_counts ==
          Tuple(region.support_count for region in scaffold.regions)
    @test inventory.all_units_from_terminal_regions
    @test !inventory.final_retained_unit_inventory_available
    @test !inventory.pair_inventory_available
    @test inventory.pair_inventory_status ==
          :not_available_terminal_region_metadata_only
    @test !inventory.coefficient_maps_materialized
    @test !inventory.transform_contracts_materialized
    @test !inventory.retained_spaces_materialized
    @test !inventory.operator_blocks_materialized
    @test !inventory.pair_operator_blocks_materialized
    @test !inventory.hamiltonian_data_materialized
    @test !inventory.artifacts_materialized
    @test inventory.diagnostics.terminal_region_metadata_only
    @test inventory.diagnostics.all_units_from_terminal_regions
    @test inventory.diagnostics.no_aggregate_atom_box_units
    @test inventory.diagnostics.no_atom_growth_route_units
    @test !inventory.diagnostics.shellification_regions_are_cpbs
    @test !inventory.diagnostics.final_retained_unit_inventory_available
    @test !inventory.diagnostics.pair_inventory_available
    @test !inventory.diagnostics.materialization_behavior_changed
    @test !inventory.diagnostics.public_default_behavior_changed

    @test all(
        record.unit_key ∉ (:left_atom_box, :right_atom_box)
        for record in inventory.terminal_region_units
    )
    @test all(
        record.unit_role ∉ (:left_atom_box, :right_atom_box)
        for record in inventory.terminal_region_units
    )
    @test all(
        record.lowering_family_planned !=
        :white_lindsey_atom_local_child_shellification
        for record in inventory.terminal_region_units
    )

    for (record, region) in zip(inventory.terminal_region_units, scaffold.regions)
        @test record.object_kind == :cartesian_terminal_region_unit_record
        @test record.unit_index == region.order_index
        @test record.terminal_region_order_index == region.order_index
        @test record.terminal_region_role == region.role
        @test record.terminal_region_kind == region.region_kind
        @test record.outer_box == region.outer_box
        @test record.inner_exclusion_box == region.inner_exclusion_box
        @test record.support_count == region.support_count
        @test !record.owned_support_is_cpb
        @test !record.shellification_region_is_cpb
        @test !record.shellification_region_is_lowering_source
        @test record.lowering_recipe_status == :planned_not_materialized
        @test record.source_cpb_plan_status == :planned_not_materialized
        @test record.retained_space_status == :not_materialized
        @test !record.coefficient_maps_materialized
        @test !record.transform_contracts_materialized
        @test !record.operator_blocks_materialized
        @test !record.pair_operator_blocks_materialized
        @test !record.hamiltonian_data_materialized
        @test record.final_retained_unit_status ==
              :not_available_terminal_region_metadata_only
    end

    complete_shell_units =
        filter(
            record -> record.terminal_region_kind == :complete_shell,
            inventory.terminal_region_units,
        )
    @test inventory.complete_shell_unit_count == length(complete_shell_units)
    @test inventory.lw_complete_shell_region_count == length(complete_shell_units)
    @test inventory.lw_complete_shell_cpb_count == 26 * length(complete_shell_units)
    @test inventory.lw_complete_shell_cpb_family_counts == (
        facet_cpb = 6 * length(complete_shell_units),
        edge_cpb = 12 * length(complete_shell_units),
        corner_cpb = 8 * length(complete_shell_units),
    )
    for record in complete_shell_units
        @test record.unit_kind in (
            :atom_local_complete_shell_unit,
            :shared_molecular_complete_shell_unit,
        )
        @test !record.owned_support_is_cpb
        @test record.owned_support_status ==
              :owned_shell_support_outer_minus_inner_exclusion
        @test !record.identity_lowering_planned
        @test record.lw_complete_shell_lowering.facet_count == 6
        @test record.lw_complete_shell_lowering.edge_count == 12
        @test record.lw_complete_shell_lowering.corner_count == 8
        @test record.lw_complete_shell_lowering.total_source_cpb_count == 26
        @test record.lw_complete_shell_lowering.status ==
              :planned_not_materialized
        @test record.pqs_complete_shell_lowering.source_cpb_count == 1
        @test record.pqs_complete_shell_lowering.source_cpb_plan_kind ==
              :filled_source_cpb
        @test record.pqs_complete_shell_lowering.retained_rule ==
              :boundary_comx_product_mode_selection
        @test !record.pqs_complete_shell_lowering.face_edge_corner_decomposition_required
        @test record.pqs_filled_source_cpb_plan ===
              record.pqs_complete_shell_lowering
    end

    direct_units =
        filter(
            record ->
                record.terminal_region_kind in (
                    :direct_core,
                    :direct_midpoint_slab,
                    :outer_mismatch_slab,
                ),
            inventory.terminal_region_units,
    )
    for record in direct_units
        @test record.owned_support_status == :owned_support_equals_region_box
        @test record.identity_lowering_planned
        @test record.source_cpb_plan.source_cpb_plan_box == record.outer_box
        @test record.source_cpb_plan.source_cpb_plan_equals_owned_support
        @test record.source_cpb_plan.source_cpb_count == 1
        @test record.source_cpb_plan.status == :planned_not_materialized
        @test !record.source_cpb_plan.coefficient_maps_materialized
        @test !record.source_cpb_plan.operator_blocks_materialized
    end

    central_units =
        filter(
            record ->
                record.terminal_region_kind == :central_distorted_product_box,
            inventory.terminal_region_units,
        )
    for record in central_units
        @test record.unit_kind == :central_distorted_product_box_unit
        @test !record.identity_lowering_planned
        @test record.lowering_family_planned == :distorted_product_box_comx
        @test !isnothing(record.source_mode_shape)
        @test !isnothing(record.q)
        @test !isnothing(record.L)
        @test !isnothing(record.aspect_ratio)
        @test record.source_cpb_plan.source_mode_shape == record.source_mode_shape
        @test record.source_cpb_plan.retained_rule ==
              :distorted_product_comx_all_axes
    end

    return inventory
end

function _terminal_geometry_assert_lowering_contract_inventory_contract(inventory)
    lowering_inventory =
        GaussletBases._cartesian_terminal_region_lowering_contract_inventory(
            inventory,
        )
    @test lowering_inventory.object_kind ==
          :cartesian_terminal_region_lowering_contract_inventory
    @test lowering_inventory.status ==
          :available_terminal_region_lowering_contract_inventory
    @test lowering_inventory.inventory_source == :terminal_region_unit_inventory
    @test lowering_inventory.source_object_kind == inventory.object_kind
    @test lowering_inventory.private_development_only
    @test lowering_inventory.terminal_region_unit_count == inventory.unit_count
    @test lowering_inventory.lowering_contract_count ==
          length(lowering_inventory.lowering_contracts)
    @test lowering_inventory.contract_records ==
          lowering_inventory.lowering_contracts
    @test lowering_inventory.unit_keys == inventory.unit_keys
    @test lowering_inventory.unit_roles == inventory.unit_roles
    @test lowering_inventory.unit_kinds == inventory.unit_kinds
    @test lowering_inventory.terminal_region_roles ==
          inventory.terminal_region_roles
    @test lowering_inventory.terminal_region_kinds ==
          inventory.terminal_region_kinds
    @test lowering_inventory.support_counts == inventory.support_counts
    @test all(
        entry.lowering_contract_count >= 1
        for entry in lowering_inventory.contract_counts_by_unit
    )
    @test !lowering_inventory.final_retained_unit_inventory_available
    @test !lowering_inventory.pair_inventory_available
    @test lowering_inventory.pair_inventory_status ==
          :not_available_lowering_contract_metadata_only
    @test !lowering_inventory.coefficient_maps_materialized
    @test !lowering_inventory.transform_contracts_materialized
    @test !lowering_inventory.retained_spaces_materialized
    @test !lowering_inventory.operator_blocks_materialized
    @test !lowering_inventory.pair_operator_blocks_materialized
    @test !lowering_inventory.hamiltonian_data_materialized
    @test !lowering_inventory.artifacts_materialized
    @test lowering_inventory.diagnostics.lowering_contracts_metadata_only
    @test lowering_inventory.diagnostics.all_units_have_lowering_contracts
    @test !lowering_inventory.diagnostics.final_retained_unit_inventory_available
    @test !lowering_inventory.diagnostics.pair_inventory_available
    @test !lowering_inventory.diagnostics.coefficient_maps_materialized
    @test !lowering_inventory.diagnostics.transform_contracts_materialized
    @test !lowering_inventory.diagnostics.retained_spaces_materialized
    @test !lowering_inventory.diagnostics.operator_blocks_materialized
    @test !lowering_inventory.diagnostics.pair_operator_blocks_materialized
    @test !lowering_inventory.diagnostics.hamiltonian_data_materialized
    @test !lowering_inventory.diagnostics.artifacts_materialized
    @test !lowering_inventory.diagnostics.materialization_behavior_changed
    @test !lowering_inventory.diagnostics.public_default_behavior_changed

    for contract in lowering_inventory.lowering_contracts
        @test contract.object_kind == :cartesian_terminal_region_lowering_contract
        @test contract.status == :planned_not_materialized
        @test contract.unit_key in inventory.unit_keys
        @test contract.outer_box == contract.box
        @test !contract.owned_support_is_cpb
        @test !contract.shellification_region_is_cpb
        @test !contract.shellification_region_is_lowering_source
        @test !contract.coefficient_maps_materialized
        @test !contract.transform_contracts_materialized
        @test !contract.retained_spaces_materialized
        @test !contract.operator_blocks_materialized
        @test !contract.pair_operator_blocks_materialized
        @test !contract.hamiltonian_data_materialized
        @test !contract.artifacts_materialized
        @test contract.final_retained_unit_status == :not_materialized
        @test !contract.final_retained_unit_records_materialized
    end

    direct_contracts =
        filter(
            contract ->
                contract.lowering_contract_kind in (
                    :direct_core_identity_cpb,
                    :direct_slab_identity_cpb,
                    :direct_boundary_slab_identity_cpb,
                ),
            lowering_inventory.lowering_contracts,
        )
    for contract in direct_contracts
        @test contract.identity_like_source_contract
        @test contract.source_cpb_plan_box == contract.outer_box
        @test contract.source_cpb_plan_equals_owned_support
        @test contract.source_cpb_count == 1
        @test contract.source_cpb_plan_kind == :direct_coordinate_product_source
        @test contract.retained_rule == :direct_source_modes
        @test !contract.face_edge_corner_decomposition_required
        @test contract.final_unit_count_planned == 1
    end

    complete_shell_units =
        filter(
            record -> record.terminal_region_kind == :complete_shell,
            inventory.terminal_region_units,
        )
    lw_contracts =
        filter(
            contract ->
                contract.lowering_contract_kind == :white_lindsey_boundary_strata,
            lowering_inventory.lowering_contracts,
        )
    @test length(lw_contracts) == length(complete_shell_units)
    for contract in lw_contracts
        @test contract.owned_support_status ==
              :owned_shell_support_outer_minus_inner_exclusion
        @test !contract.owned_support_is_cpb
        @test contract.source_cpb_plan_kind == :complete_shell_boundary_strata
        @test contract.source_cpb_family_counts ==
              (facet_cpb = 6, edge_cpb = 12, corner_cpb = 8)
        @test contract.source_cpb_count == 26
        @test !contract.coefficient_maps_materialized
        @test contract.final_unit_granularity ==
              :white_lindsey_boundary_strata_children_planned
        @test contract.final_unit_count_planned == 26
    end

    pqs_contracts =
        filter(
            contract -> contract.lowering_contract_kind == :pqs_filled_source_cpb,
            lowering_inventory.lowering_contracts,
        )
    @test length(pqs_contracts) == length(complete_shell_units)
    for contract in pqs_contracts
        @test contract.source_cpb_count == 1
        @test contract.source_cpb_plan_box == contract.outer_box
        @test contract.source_cpb_plan_kind == :filled_source_cpb
        @test contract.retained_rule == :boundary_comx_product_mode_selection
        @test contract.intermediate_retained_space_status ==
              :planned_not_materialized
        @test contract.shell_realization_status ==
              :projection_lowdin_planned_not_materialized
        @test !contract.face_edge_corner_decomposition_required
        @test !contract.identity_like_source_contract
        @test contract.final_unit_count_planned == 1
    end

    central_contracts =
        filter(
            contract ->
                contract.lowering_contract_kind == :distorted_product_box_comx,
            lowering_inventory.lowering_contracts,
        )
    central_units =
        filter(
            record ->
                record.terminal_region_kind == :central_distorted_product_box,
            inventory.terminal_region_units,
        )
    @test length(central_contracts) == length(central_units)
    for contract in central_contracts
        matching_unit = only(
            record for record in central_units
            if record.unit_key == contract.unit_key
        )
        @test !contract.identity_like_source_contract
        @test contract.source_cpb_plan_box == contract.outer_box
        @test contract.source_mode_shape == matching_unit.source_mode_shape
        @test contract.q == matching_unit.q
        @test contract.L == matching_unit.L
        @test contract.aspect_ratio == matching_unit.aspect_ratio
        @test contract.retained_rule == :distorted_product_comx_all_axes
        @test contract.final_unit_count_planned == 1
    end

    return lowering_inventory
end

@testset "terminal Cartesian shellification one-center centered complete shells" begin
    axes = (collect(1:9), collect(1:9), collect(1:9))
    plan = GaussletBases._cartesian_terminal_shellification_geometry(
        axes,
        (5, 5, 5);
        core_side = 5,
    )

    @test plan.object_kind == :cartesian_terminal_shellification_geometry_plan
    @test plan.system_kind == :one_center
    @test plan.nuclear_indices == ((5, 5, 5),)
    @test plan.region_count == 3
    @test plan.region_roles == (:atom_local_core, :atom_local_shell, :atom_local_shell)
    @test plan.region_kinds == (:direct_core, :complete_shell, :complete_shell)
    @test plan.regions[1].support_count == 5^3
    @test plan.regions[2].support_count == 7^3 - 5^3
    @test plan.regions[3].support_count == 9^3 - 7^3
    _terminal_geometry_assert_full_partition(plan)
    summary = _terminal_geometry_assert_private_summary_contract(plan)
    expected_dependencies = (
        :plan_lowerable_direct_core,
        :plan_lowerable_complete_shell,
        :plan_lowerable_complete_shell,
    )
    dependency_counts = summary.materialization_dependency_counts
    @test summary.ordered_materialization_dependencies == expected_dependencies
    @test dependency_counts.plan_lowerable_direct_core_count == 1
    @test dependency_counts.plan_lowerable_complete_shell_count == 2
    @test dependency_counts.plan_lowerable_shared_complete_shell_count == 0
    @test dependency_counts.plan_lowerable_direct_slab_count == 0
    @test dependency_counts.planned_distorted_product_box_lowering_count == 0
    @test dependency_counts.unsupported_terminal_shellification_region_count == 0
    @test summary.central_gap_region_count == 0
    scaffold = _terminal_geometry_assert_scaffold_contract(plan)
    @test scaffold.spatial_policy_order == :single_center_outward
    @test scaffold.materialization_status == :ready_supported_terminal_subset
    @test scaffold.independently_lowerable_region_count == plan.region_count
    @test scaffold.deferred_region_count == 0
    @test scaffold.unsupported_region_count == 0
    @test all(region.independently_lowerable for region in scaffold.regions)
    @test all(
        region.retirement_target == :already_plan_lowered_region
        for region in scaffold.regions
    )
    inventory = _terminal_geometry_assert_region_unit_inventory_contract(scaffold)
    _terminal_geometry_assert_lowering_contract_inventory_contract(inventory)
    @test inventory.unit_kinds == (
        :direct_core_unit,
        :atom_local_complete_shell_unit,
        :atom_local_complete_shell_unit,
    )
end

@testset "terminal Cartesian shellification one-center outer mismatch slabs" begin
    axes = (collect(1:11), collect(1:9), collect(1:7))
    plan = GaussletBases._cartesian_terminal_shellification_geometry(
        axes,
        (4, 4, 4);
        core_side = 3,
    )

    @test plan.system_kind == :one_center
    @test :outer_mismatch_slab in plan.region_kinds
    @test any(region.role == :x_high_outer_mismatch_slab for region in plan.regions)
    @test any(region.role == :y_high_outer_mismatch_slab for region in plan.regions)
    _terminal_geometry_assert_full_partition(plan)
    summary = _terminal_geometry_assert_private_summary_contract(plan)
    dependency_counts = summary.materialization_dependency_counts
    @test dependency_counts.plan_lowerable_direct_boundary_slab_count > 0
    @test all(
        region.materialization_dependency == :plan_lowerable_direct_boundary_slab
        for region in summary.regions
        if region.region_kind == :outer_mismatch_slab
    )
    scaffold = _terminal_geometry_assert_scaffold_contract(plan)
    inventory = _terminal_geometry_assert_region_unit_inventory_contract(scaffold)
    _terminal_geometry_assert_lowering_contract_inventory_contract(inventory)
    mismatch_units =
        filter(
            record -> record.terminal_region_kind == :outer_mismatch_slab,
            inventory.terminal_region_units,
        )
    @test !isempty(mismatch_units)
    @test all(
        record.unit_kind == :outer_mismatch_boundary_slab_unit
        for record in mismatch_units
    )
    @test all(
        record.lowering_family_planned == :direct_boundary_slab_identity_cpb
        for record in mismatch_units
    )
end

@testset "terminal Cartesian shellification coordinate snapping" begin
    axes = (
        collect(range(-1.0, 1.0; length = 9)),
        collect(range(-2.0, 2.0; length = 9)),
        collect(range(-3.0, 3.0; length = 9)),
    )
    plan = GaussletBases._cartesian_terminal_shellification_geometry(
        axes,
        (0.1, -0.2, 0.4);
        core_side = 3,
    )

    @test plan.system_kind == :one_center
    @test plan.nuclear_indices == ((5, 5, 6),)
    _terminal_geometry_assert_full_partition(plan)
end

@testset "terminal Cartesian shellification diatomic with midpoint and shared shells" begin
    axes = (collect(1:15), collect(1:11), collect(1:11))
    plan = GaussletBases._cartesian_terminal_shellification_geometry(
        axes,
        ((5, 6, 6), (11, 6, 6));
        core_side = 3,
    )

    @test plan.system_kind == :bond_aligned_diatomic
    @test plan.bond_axis == :x
    @test plan.region_roles[1:2] == (:atom_local_core, :atom_local_core)
    @test any(region.role == :midpoint_slab for region in plan.regions)
    @test any(region.role == :shared_molecular_shell for region in plan.regions)
    @test all(region.role ∉ (:left_atom_box, :right_atom_box) for region in plan.regions)
    _terminal_geometry_assert_full_partition(plan)
    summary = _terminal_geometry_assert_private_summary_contract(plan)
    dependency_counts = summary.materialization_dependency_counts
    @test dependency_counts.plan_lowerable_shared_complete_shell_count > 0
    @test dependency_counts.plan_lowerable_direct_slab_count == 1
    scaffold = _terminal_geometry_assert_scaffold_contract(plan)
    @test scaffold.spatial_policy_order == :atom_outward
    @test scaffold.materialization_status == :ready_supported_terminal_subset
    @test scaffold.independently_lowerable_region_count == plan.region_count
    @test scaffold.deferred_region_count == 0
    midpoint_scaffold_regions =
        filter(region -> region.region_kind == :direct_midpoint_slab, scaffold.regions)
    @test length(midpoint_scaffold_regions) == 1
    @test only(midpoint_scaffold_regions).lowering_family == :direct_midpoint_slab
    @test only(midpoint_scaffold_regions).independently_lowerable
    inventory = _terminal_geometry_assert_region_unit_inventory_contract(scaffold)
    _terminal_geometry_assert_lowering_contract_inventory_contract(inventory)
    @test any(
        record.unit_kind == :shared_molecular_complete_shell_unit
        for record in inventory.terminal_region_units
    )
    midpoint_units =
        filter(
            record -> record.terminal_region_kind == :direct_midpoint_slab,
            inventory.terminal_region_units,
        )
    @test length(midpoint_units) == 1
    @test only(midpoint_units).unit_kind == :direct_midpoint_slab_unit
    @test only(midpoint_units).lowering_family_planned ==
          :direct_slab_identity_cpb
end

@testset "terminal Cartesian shellification diatomic without midpoint" begin
    axes = (collect(1:15), collect(1:11), collect(1:11))
    plan = GaussletBases._cartesian_terminal_shellification_geometry(
        axes,
        ((5, 6, 6), (10, 6, 6));
        core_side = 3,
    )

    @test plan.system_kind == :bond_aligned_diatomic
    @test !any(region.role == :midpoint_slab for region in plan.regions)
    @test any(region.role == :shared_molecular_shell for region in plan.regions)
    _terminal_geometry_assert_full_partition(plan)
end

@testset "terminal Cartesian shellification central gap below q becomes midpoint slabs" begin
    axes = (collect(1:15), collect(1:7), collect(1:7))
    plan = GaussletBases._cartesian_terminal_shellification_geometry(
        axes,
        ((3, 4, 4), (13, 4, 4));
        core_side = 3,
        q = 7,
    )

    midpoint_slabs = filter(region -> region.role == :midpoint_slab, plan.regions)
    @test length(midpoint_slabs) == 5
    @test all(region.region_kind == :direct_midpoint_slab for region in midpoint_slabs)
    @test all(region.metadata.central_gap_width == 5 for region in midpoint_slabs)
    @test all(region.support_count == 5^2 for region in midpoint_slabs)
    @test !any(region.role == :central_distorted_product_box for region in plan.regions)
    _terminal_geometry_assert_full_partition(plan)
    summary = _terminal_geometry_assert_private_summary_contract(plan)
    dependency_counts = summary.materialization_dependency_counts
    @test summary.central_gap_region_count == 5
    @test summary.central_midpoint_slab_count == 5
    @test summary.central_distorted_product_box_count == 0
    @test dependency_counts.plan_lowerable_direct_slab_count == 5
    @test dependency_counts.planned_distorted_product_box_lowering_count == 0
    @test all(
        region.materialization_dependency == :plan_lowerable_direct_slab
        for region in summary.regions
        if region.region_kind == :direct_midpoint_slab
    )
    scaffold = _terminal_geometry_assert_scaffold_contract(plan)
    @test scaffold.materialization_status == :ready_supported_terminal_subset
    @test scaffold.central_gap_region_count == 5
    @test scaffold.central_midpoint_slab_count == 5
    @test scaffold.central_distorted_product_box_count == 0
    midpoint_scaffold_regions =
        filter(region -> region.region_kind == :direct_midpoint_slab, scaffold.regions)
    @test length(midpoint_scaffold_regions) == 5
    @test all(region.independently_lowerable for region in midpoint_scaffold_regions)
    @test all(
        region.materialization_dependency == :plan_lowerable_direct_slab
        for region in midpoint_scaffold_regions
    )
    inventory = _terminal_geometry_assert_region_unit_inventory_contract(scaffold)
    _terminal_geometry_assert_lowering_contract_inventory_contract(inventory)
    midpoint_units =
        filter(
            record -> record.terminal_region_kind == :direct_midpoint_slab,
            inventory.terminal_region_units,
        )
    @test length(midpoint_units) == 5
    @test all(
        record.unit_kind == :direct_midpoint_slab_unit
        for record in midpoint_units
    )
    @test all(
        record.lowering_family_planned == :direct_slab_identity_cpb
        for record in midpoint_units
    )
end

@testset "terminal Cartesian shellification central gap at least q uses one region" begin
    axes = (collect(1:21), collect(1:7), collect(1:7))
    plan = GaussletBases._cartesian_terminal_shellification_geometry(
        axes,
        ((3, 4, 4), (19, 4, 4));
        core_side = 3,
        q = 3,
    )

    central_boxes =
        filter(region -> region.role == :central_distorted_product_box, plan.regions)
    @test length(central_boxes) == 1
    central_box = only(central_boxes)
    @test central_box.region_kind == :central_distorted_product_box
    @test central_box.support_count == 11 * 5 * 5
    @test central_box.metadata.bond_axis == :x
    @test central_box.metadata.central_gap_width == 11
    @test central_box.metadata.q == 3
    @test central_box.metadata.L == 7
    @test central_box.metadata.source_mode_shape == (7, 3, 3)
    @test central_box.metadata.lowering_hint == :distorted_comx_all_axes
    @test !any(region.role == :midpoint_slab for region in plan.regions)
    _terminal_geometry_assert_full_partition(plan)
    summary = _terminal_geometry_assert_private_summary_contract(plan)
    dependency_counts = summary.materialization_dependency_counts
    distorted_metadata = only(summary.central_distorted_product_box_metadata)
    @test summary.central_gap_region_count == 1
    @test summary.central_midpoint_slab_count == 0
    @test summary.central_distorted_product_box_count == 1
    @test dependency_counts.plan_lowerable_direct_slab_count == 0
    @test dependency_counts.planned_distorted_product_box_lowering_count == 1
    @test distorted_metadata.q == 3
    @test distorted_metadata.L == 7
    @test distorted_metadata.source_mode_shape == (7, 3, 3)
    @test only(
        region for region in summary.regions
        if region.region_kind == :central_distorted_product_box
    ).materialization_dependency == :planned_distorted_product_box_lowering
    scaffold = _terminal_geometry_assert_scaffold_contract(plan)
    @test scaffold.materialization_status ==
          :deferred_pending_distorted_product_box_lowering
    @test scaffold.central_gap_region_count == 1
    @test scaffold.central_midpoint_slab_count == 0
    @test scaffold.central_distorted_product_box_count == 1
    @test scaffold.deferred_region_count == 1
    central_scaffold_region = only(
        region for region in scaffold.regions
        if region.region_kind == :central_distorted_product_box
    )
    @test !central_scaffold_region.independently_lowerable
    @test central_scaffold_region.materialization_dependency ==
          :planned_distorted_product_box_lowering
    @test central_scaffold_region.lowering_family ==
          :distorted_comx_product_box_deferred
    @test central_scaffold_region.lowering_hint == :distorted_comx_all_axes
    @test central_scaffold_region.missing_independent_lowering_reason ==
          :distorted_product_box_lowering_pending
    @test central_scaffold_region.retirement_target ==
          :pending_distorted_product_box_lowering_support
    @test central_scaffold_region.metadata.q == 3
    @test central_scaffold_region.metadata.L == 7
    @test central_scaffold_region.metadata.source_mode_shape == (7, 3, 3)
    inventory = _terminal_geometry_assert_region_unit_inventory_contract(scaffold)
    _terminal_geometry_assert_lowering_contract_inventory_contract(inventory)
    central_unit = only(
        record for record in inventory.terminal_region_units
        if record.terminal_region_kind == :central_distorted_product_box
    )
    @test central_unit.source_mode_shape == (7, 3, 3)
    @test central_unit.q == 3
    @test central_unit.L == 7
    @test central_unit.source_cpb_plan.source_cpb_plan_kind ==
          :central_distorted_product_box_source
    @test !central_unit.identity_lowering_planned
end

@testset "terminal Cartesian shellification invalid inputs and unsupported geometries" begin
    axes = (collect(1:9), collect(1:9), collect(1:9))
    @test_throws ArgumentError GaussletBases._cartesian_terminal_shellification_geometry(
        axes,
        (5, 5, 5);
        core_side = 4,
    )
    @test_throws ArgumentError GaussletBases._cartesian_terminal_shellification_geometry(
        axes,
        (5, 5, 5);
        q = 0,
    )
    @test_throws ArgumentError GaussletBases._cartesian_terminal_shellification_geometry(
        axes,
        ((4, 4, 4), (6, 5, 4));
        core_side = 3,
    )
    @test_throws ArgumentError GaussletBases._cartesian_terminal_shellification_geometry(
        axes,
        ((4, 5, 5), (6, 5, 5));
        core_side = 3,
    )
end
