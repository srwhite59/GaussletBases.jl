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
