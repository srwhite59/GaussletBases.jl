using Test
using GaussletBases

const CSH = GaussletBases.CartesianShellification
const CRCForShellification = GaussletBases.CartesianRouteCore

function _shellification_module_compact_fingerprint(plan)
    raw = CSH.raw_plan(plan)
    return (;
        plan_type = typeof(plan),
        raw_region_count = raw.region_count,
        typed_region_count = length(CSH.terminal_regions(plan)),
        coverage_complete = CSH.coverage(plan).coverage_complete,
        summary_region_count = CSH.summary(plan).region_count,
        region_kinds = Tuple(region.region_kind for region in CSH.terminal_regions(plan)),
        support_counts =
            Tuple(CRCForShellification.support_count(region.owned_support)
                  for region in CSH.terminal_regions(plan)),
        route_region_roles = Tuple(
            CRCForShellification.role(region.route_region)
            for region in CSH.terminal_regions(plan)
        ),
        raw_region_roles = raw.region_roles,
    )
end

function _shellification_raw_compact_fingerprint(raw)
    return (;
        object_kind = raw.object_kind,
        status = raw.status,
        system_kind = raw.system_kind,
        parent_dims = raw.parent_dims,
        nuclear_indices = raw.nuclear_indices,
        bond_axis = raw.bond_axis,
        core_side = raw.core_side,
        q = raw.q,
        region_count = raw.region_count,
        terminal_region_count = raw.terminal_region_count,
        region_roles = raw.region_roles,
        region_kinds = raw.region_kinds,
        support_counts = Tuple(region.support_count for region in raw.regions),
        coverage_complete = raw.coverage.coverage_complete,
        shellification_geometry_only =
            raw.diagnostics.shellification_geometry_only,
        coordinate_product_box_lowering_materialized =
            raw.diagnostics.coordinate_product_box_lowering_materialized,
        retained_spaces_materialized =
            raw.diagnostics.retained_spaces_materialized,
        coefficient_maps_materialized =
            raw.diagnostics.coefficient_maps_materialized,
        operator_blocks_materialized =
            raw.diagnostics.operator_blocks_materialized,
        pair_operator_blocks_materialized =
            raw.diagnostics.pair_operator_blocks_materialized,
        hamiltonian_data_materialized =
            raw.diagnostics.hamiltonian_data_materialized,
    )
end

function _shellification_summary_compact_fingerprint(summary)
    return (;
        object_kind = summary.object_kind,
        status = summary.status,
        source_kind = summary.source_kind,
        source_object_kind = summary.source_object_kind,
        system_kind = summary.system_kind,
        parent_dims = summary.parent_dims,
        nuclear_indices = summary.nuclear_indices,
        bond_axis = summary.bond_axis,
        core_side = summary.core_side,
        q = summary.q,
        ordered_terminal_region_roles =
            summary.ordered_terminal_region_roles,
        ordered_terminal_region_kinds =
            summary.ordered_terminal_region_kinds,
        ordered_materialization_dependencies =
            summary.ordered_materialization_dependencies,
        region_count = summary.region_count,
        terminal_region_count = summary.terminal_region_count,
        region_support_counts = summary.region_support_counts,
        total_support_count = summary.total_support_count,
        coverage_status = summary.coverage_status,
        coverage_complete = summary.coverage_complete,
        central_gap_region_count = summary.central_gap_region_count,
        central_midpoint_slab_count = summary.central_midpoint_slab_count,
        central_distorted_product_box_count =
            summary.central_distorted_product_box_count,
        shellification_regions_are_cpbs =
            summary.shellification_regions_are_cpbs,
        shellification_regions_are_lowering_sources =
            summary.shellification_regions_are_lowering_sources,
        lowering_applied_by_summary = summary.lowering_applied_by_summary,
        retained_spaces_materialized = summary.retained_spaces_materialized,
        coefficient_maps_materialized = summary.coefficient_maps_materialized,
        operator_blocks_materialized = summary.operator_blocks_materialized,
        pair_operator_blocks_materialized =
            summary.pair_operator_blocks_materialized,
        hamiltonian_data_materialized = summary.hamiltonian_data_materialized,
        public_default_behavior_changed =
            summary.public_default_behavior_changed,
    )
end

function _shellification_scaffold_compact_fingerprint(scaffold)
    return (;
        object_kind = scaffold.object_kind,
        status = scaffold.status,
        materialization_status = scaffold.materialization_status,
        source_kind = scaffold.source_kind,
        route_family = scaffold.route_family,
        system_classification = scaffold.system_classification,
        shellification_role = scaffold.shellification_role,
        shellification_stage = scaffold.shellification_stage,
        lowering_stage = scaffold.lowering_stage,
        spatial_policy_order = scaffold.spatial_policy_order,
        parent_dims = scaffold.parent_dims,
        nuclear_indices = scaffold.nuclear_indices,
        bond_axis = scaffold.bond_axis,
        core_side = scaffold.core_side,
        q = scaffold.q,
        region_count = scaffold.region_count,
        ordered_region_roles = scaffold.ordered_region_roles,
        ordered_region_kinds = scaffold.ordered_region_kinds,
        ordered_materialization_dependencies =
            scaffold.ordered_materialization_dependencies,
        support_counts =
            Tuple(region.support_count for region in scaffold.regions),
        independently_lowerable_region_count =
            scaffold.independently_lowerable_region_count,
        deferred_region_count = scaffold.deferred_region_count,
        unsupported_region_count = scaffold.unsupported_region_count,
        central_gap_region_count = scaffold.central_gap_region_count,
        central_midpoint_slab_count = scaffold.central_midpoint_slab_count,
        central_distorted_product_box_count =
            scaffold.central_distorted_product_box_count,
        coverage_complete = scaffold.coverage.coverage_complete,
        terminal_geometry_summary_status =
            scaffold.terminal_geometry_summary.status,
        lowering_applied_by_plan =
            scaffold.diagnostics.lowering_applied_by_plan,
        materialization_behavior_changed =
            scaffold.diagnostics.materialization_behavior_changed,
        public_default_behavior_changed =
            scaffold.diagnostics.public_default_behavior_changed,
        shellification_regions_are_cpbs =
            scaffold.diagnostics.shellification_regions_are_cpbs,
        shellification_regions_are_lowering_sources =
            scaffold.diagnostics.shellification_regions_are_lowering_sources,
        retained_spaces_materialized =
            scaffold.diagnostics.retained_spaces_materialized,
        coefficient_maps_materialized =
            scaffold.diagnostics.coefficient_maps_materialized,
        operator_blocks_materialized =
            scaffold.diagnostics.operator_blocks_materialized,
        pair_operator_blocks_materialized =
            scaffold.diagnostics.pair_operator_blocks_materialized,
        hamiltonian_data_materialized =
            scaffold.diagnostics.hamiltonian_data_materialized,
    )
end

@testset "CartesianShellification module facade" begin
    parent_axes = (collect(1:9), collect(1:9), collect(1:13))
    policy = CSH.AtomOutwardShellification(;
        core_side = 3,
        q = 3,
        bond_axis = :z,
    )
    plan = CSH.shellify(parent_axes, ((5, 5, 4), (5, 5, 10)), policy)
    raw = CSH.raw_plan(plan)
    fingerprint = _shellification_module_compact_fingerprint(plan)
    module_raw = CSH.raw_terminal_geometry(
        parent_axes,
        ((5, 5, 4), (5, 5, 10));
        core_side = 3,
        q = 3,
        bond_axis = :z,
    )
    wrapper_raw = GaussletBases._cartesian_terminal_shellification_geometry(
        parent_axes,
        ((5, 5, 4), (5, 5, 10));
        core_side = 3,
        q = 3,
        bond_axis = :z,
    )
    module_summary = CSH.private_summary(raw)
    wrapper_summary =
        GaussletBases._cartesian_terminal_shellification_geometry_private_summary(raw)
    module_scaffold = CSH.scaffold(raw; route_family = :white_lindsey_low_order)
    wrapper_scaffold =
        GaussletBases._cartesian_terminal_shellification_geometry_scaffold(
            raw;
            route_family = :white_lindsey_low_order,
        )
    plan_scaffold = CSH.scaffold(plan; route_family = :white_lindsey_low_order)

    @test plan isa CSH.ShellificationPlan
    @test plan.policy === policy
    @test fingerprint.raw_region_count == raw.region_count
    @test fingerprint.typed_region_count == raw.region_count
    @test fingerprint.summary_region_count == raw.region_count
    @test fingerprint.coverage_complete
    @test fingerprint.raw_region_roles == fingerprint.route_region_roles
    @test fingerprint.support_counts ==
          Tuple(region.support_count for region in raw.regions)
    @test all(region -> region isa CSH.TerminalRegion, CSH.terminal_regions(plan))
    @test all(
        region -> region.owned_support isa CRCForShellification.OwnedSupport,
        CSH.terminal_regions(plan),
    )
    @test all(
        region -> region.route_region isa CRCForShellification.ShellificationRegion,
        CSH.terminal_regions(plan),
    )
    @test !raw.aggregate_atom_boxes_emitted
    @test raw.diagnostics.shellification_geometry_only
    @test !raw.diagnostics.coordinate_product_box_lowering_materialized
    @test _shellification_raw_compact_fingerprint(module_raw) ==
          _shellification_raw_compact_fingerprint(wrapper_raw)
    @test _shellification_raw_compact_fingerprint(raw) ==
          _shellification_raw_compact_fingerprint(wrapper_raw)
    @test _shellification_summary_compact_fingerprint(module_summary) ==
          _shellification_summary_compact_fingerprint(wrapper_summary)
    @test _shellification_summary_compact_fingerprint(CSH.summary(plan)) ==
          _shellification_summary_compact_fingerprint(wrapper_summary)
    @test _shellification_scaffold_compact_fingerprint(module_scaffold) ==
          _shellification_scaffold_compact_fingerprint(wrapper_scaffold)
    @test _shellification_scaffold_compact_fingerprint(plan_scaffold) ==
          _shellification_scaffold_compact_fingerprint(wrapper_scaffold)
    @test !module_summary.shellification_regions_are_cpbs
    @test !module_summary.shellification_regions_are_lowering_sources
    @test !module_summary.lowering_applied_by_summary
    @test !module_scaffold.diagnostics.lowering_applied_by_plan
    @test !module_scaffold.diagnostics.retained_spaces_materialized
    @test !module_scaffold.diagnostics.operator_blocks_materialized

    one_center = CSH.shellify(
        parent_axes,
        (5, 5, 7);
        policy = CSH.OneCenterShellification(; core_side = 3, q = 3),
    )
    @test one_center.policy isa CSH.OneCenterShellification
    @test CSH.raw_plan(one_center).system_kind == :one_center
    @test CSH.coverage(one_center).coverage_complete
    @test length(CSH.terminal_regions(one_center)) ==
          CSH.raw_plan(one_center).region_count
end
