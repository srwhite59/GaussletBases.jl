using Test
using GaussletBases

function _one_center_shellification_plan_fixture(; count::Int = 9, nside::Int = 5)
    expansion = coulomb_gaussian_expansion(doacc = false)
    basis = build_basis(MappedUniformBasisSpec(:G10;
        count,
        mapping = white_lindsey_atomic_mapping(; Z = 4.0, d = 0.15),
        reference_spacing = 1.0,
    ))
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    sequence = build_one_center_atomic_full_parent_shell_sequence(
        bundle;
        expansion,
        nside,
    )
    specific_plan = GaussletBases._cartesian_shellification_plan_one_center_low_order(
        count;
        nside,
    )
    plan = GaussletBases._cartesian_shellification_plan(
        :one_center_full_parent_low_order,
        count;
        nside,
    )
    audit = GaussletBases._nested_shell_sequence_contract_audit(
        sequence,
        (count, count, count),
    )
    diagnostics = one_center_atomic_nested_structure_diagnostics(
        sequence;
        parent_side_count = count,
        nside,
    )
    return (;
        count,
        nside,
        basis,
        bundle,
        expansion,
        sequence,
        specific_plan,
        plan,
        audit,
        diagnostics,
    )
end

function _bond_aligned_diatomic_shellification_plan_fixture(;
    nside::Int = 5,
    shared_shell_layer_policy::Symbol = :endcap_panel_owned,
    xmax_parallel::Float64 = 6.0,
    xmax_transverse::Float64 = 4.0,
    min_unsplit_parallel_to_transverse_ratio_for_split::Float64 = 3.0,
)
    basis = bond_aligned_homonuclear_qw_basis(
        family = :G10,
        bond_length = 1.4,
        core_spacing = 0.7,
        xmax_parallel = xmax_parallel,
        xmax_transverse = xmax_transverse,
        bond_axis = :z,
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(basis, expansion)
    source = GaussletBases._nested_bond_aligned_diatomic_source(
        basis,
        bundles;
        bond_axis = :z,
        nside,
        term_coefficients = Float64.(expansion.coefficients),
        packet_kernel = :factorized_direct,
        min_unsplit_parallel_to_transverse_ratio_for_split =
            min_unsplit_parallel_to_transverse_ratio_for_split,
        shared_shell_layer_policy = shared_shell_layer_policy,
        shared_shell_endcap_panel_q = 4,
        shared_shell_endcap_panel_L = 4,
    )
    specific_plan =
        GaussletBases._cartesian_shellification_plan_bond_aligned_diatomic_low_order(
            source,
        )
    plan = GaussletBases._cartesian_shellification_plan(
        :bond_aligned_diatomic_active_source_low_order,
        source,
    )
    audit = GaussletBases._nested_source_contract_audit(source)
    return (; nside, basis, bundles, expansion, source, specific_plan, plan, audit)
end

@testset "one-center low-order shellification plan matches legacy sequence" begin
    fixture = _one_center_shellification_plan_fixture()
    count = fixture.count
    nside = fixture.nside
    sequence = fixture.sequence
    specific_plan = fixture.specific_plan
    plan = fixture.plan
    audit = fixture.audit
    diagnostics = fixture.diagnostics

    @test plan == specific_plan
    @test plan.object_kind == :cartesian_shellification_plan3d
    @test plan.status == :planned_metadata_only
    @test plan.private_development_only
    @test plan.source_kind == :one_center_atomic_full_parent_shellification_plan
    @test plan.route_family == :white_lindsey_low_order
    @test plan.system_classification == :one_center
    @test plan.shellification_stage == :route_neutral_spatial_planning
    @test plan.lowering_stage == :not_lowered_by_shellification_plan
    @test plan.parent_box == sequence.working_box
    @test plan.working_box == sequence.working_box
    @test plan.full_parent_working_box
    @test plan.nside == nside
    @test plan.shell_layer_count == length(sequence.shell_layers)
    @test plan.shell_layer_count == diagnostics.shell_layer_count
    @test plan.core_side_count == diagnostics.core_side_count
    @test plan.direct_core_point_count == length(sequence.core_indices)
    @test plan.retained_dimension == size(sequence.coefficient_matrix, 2)
    @test plan.retained_dimension == diagnostics.total_actual_gausslet_count
    @test plan.region_count == plan.shell_layer_count + 1
    @test plan.shell_region_count == plan.shell_layer_count
    @test plan.direct_core_region_count == 1
    @test plan.ordered_region_roles ==
          Tuple(vcat(fill(:low_order_complete_shell, plan.shell_layer_count), :direct_core))
    @test plan.ordered_materialization_dependencies == Tuple(
        vcat(
            fill(:plan_lowerable_complete_shell, plan.shell_layer_count),
            :plan_lowerable_direct_core,
        ),
    )
    @test plan.materialization_dependency_counts.plan_lowerable_region_count ==
          plan.region_count
    @test plan.materialization_dependency_counts.source_backed_region_count == 0
    @test plan.materialization_dependency_counts.plan_lowerable_complete_shell_count ==
          plan.shell_region_count
    @test plan.materialization_dependency_counts.plan_lowerable_direct_core_count == 1
    @test plan.materialization_dependency_counts.source_box_direct_adapter_region_count == 0

    shell_regions =
        Tuple(region for region in plan.regions if region.role == :low_order_complete_shell)
    @test length(shell_regions) == length(sequence.shell_layers)
    for (region, shell, column_range) in
        zip(shell_regions, sequence.shell_layers, sequence.layer_column_ranges)
        @test region.object_kind == :cartesian_shellification_region3d
        @test region.lowering_family == :white_lindsey_complete_shell
        @test region.lowering_status == :planned_not_lowered
        @test region.materialization_dependency == :plan_lowerable_complete_shell
        @test region.box == shell.provenance.source_box
        @test region.next_inner_box == shell.provenance.next_inner_box
        @test region.box_shape == Tuple(length.(shell.provenance.source_box))
        @test region.source_point_count == shell.provenance.source_point_count
        @test region.retained_count == shell.provenance.retained_fixed_count
        @test region.retained_count == length(column_range)
    end

    direct_core_region = only(region for region in plan.regions if region.role == :direct_core)
    @test direct_core_region.lowering_family == :direct_product_core
    @test direct_core_region.lowering_status == :planned_not_lowered
    @test direct_core_region.materialization_dependency == :plan_lowerable_direct_core
    @test direct_core_region.box == plan.direct_core_box
    @test direct_core_region.box_shape == (nside, nside, nside)
    @test direct_core_region.source_point_count == nside^3
    @test direct_core_region.retained_count == length(sequence.core_column_range)
    @test direct_core_region.retained_count == length(sequence.core_indices)

    @test plan.coverage.status == :coverage_complete
    @test plan.coverage.coverage_complete
    @test plan.coverage.parent_point_count == count^3
    @test plan.coverage.parent_point_count == audit.expected_support_count
    @test plan.coverage.total_source_point_count == audit.support_count
    @test plan.coverage.missing_point_count == 0
    @test audit.full_parent_working_box
    @test audit.missing_row_count == 0
    @test audit.ownership_group_count_min == 1
    @test audit.ownership_group_count_max == 1
    @test audit.ownership_unowned_row_count == 0
    @test audit.ownership_multi_owned_row_count == 0

    @test plan.diagnostics.route_neutral_spatial_planning
    @test !plan.diagnostics.lowering_applied_by_plan
    @test !plan.diagnostics.materialization_behavior_changed
    @test !plan.diagnostics.public_default_behavior_changed
    @test !plan.diagnostics.pqs_production_source_box_materialization_claimed
    @test !plan.diagnostics.mwg_ida_semantics_changed

    plan_summary = GaussletBases._cartesian_shellification_plan_private_summary(plan)
    @test plan_summary.ordered_materialization_dependencies ==
          plan.ordered_materialization_dependencies
    @test plan_summary.materialization_dependency_counts ==
          plan.materialization_dependency_counts
    @test plan_summary.plan_lowerable_region_count == plan.region_count
    @test plan_summary.source_backed_region_count == 0
    @test plan_summary.source_box_direct_adapter_region_count == 0
end

@testset "one-center low-order shellification plan materializer matches legacy sequence" begin
    fixture = _one_center_shellification_plan_fixture()
    legacy = fixture.sequence
    plan = fixture.plan
    materialized = GaussletBases._cartesian_materialize_shellification_low_order(
        plan,
        fixture.bundle;
        expansion = fixture.expansion,
    )
    audit = GaussletBases._nested_shell_sequence_contract_audit(
        materialized,
        (fixture.count, fixture.count, fixture.count),
    )

    @test materialized.working_box == legacy.working_box
    @test materialized.working_box == plan.working_box
    @test length(materialized.shell_layers) == length(legacy.shell_layers)
    @test length(materialized.shell_layers) == plan.shell_layer_count
    @test materialized.core_indices == legacy.core_indices
    @test materialized.core_states == legacy.core_states
    @test materialized.core_column_range == legacy.core_column_range
    @test materialized.layer_column_ranges == legacy.layer_column_ranges
    @test materialized.support_indices == legacy.support_indices
    @test materialized.support_states == legacy.support_states
    @test size(materialized.coefficient_matrix) == size(legacy.coefficient_matrix)
    @test materialized.coefficient_matrix ≈ legacy.coefficient_matrix atol = 1.0e-10 rtol = 1.0e-10

    shell_regions =
        Tuple(region for region in plan.regions if region.role == :low_order_complete_shell)
    @test length(shell_regions) == length(materialized.shell_layers)
    for (region, shell, legacy_shell, column_range) in zip(
        shell_regions,
        materialized.shell_layers,
        legacy.shell_layers,
        materialized.layer_column_ranges,
    )
        @test shell isa GaussletBases._CartesianNestedCompleteShell3D
        @test shell.provenance == region.provenance
        @test shell.provenance == legacy_shell.provenance
        @test shell.support_indices == legacy_shell.support_indices
        @test shell.support_states == legacy_shell.support_states
        @test shell.face_column_ranges == legacy_shell.face_column_ranges
        @test shell.edge_column_ranges == legacy_shell.edge_column_ranges
        @test shell.corner_column_ranges == legacy_shell.corner_column_ranges
        @test size(shell.coefficient_matrix) == size(legacy_shell.coefficient_matrix)
        @test shell.coefficient_matrix ≈ legacy_shell.coefficient_matrix atol = 1.0e-10 rtol = 1.0e-10
        @test length(column_range) == region.retained_count
    end

    direct_core_region = only(region for region in plan.regions if region.role == :direct_core)
    @test direct_core_region.box == plan.direct_core_box
    @test length(materialized.core_indices) == direct_core_region.retained_count
    @test length(materialized.core_column_range) == direct_core_region.retained_count

    @test materialized.packet.overlap ≈ legacy.packet.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.kinetic ≈ legacy.packet.kinetic atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.position_x ≈ legacy.packet.position_x atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.position_y ≈ legacy.packet.position_y atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.position_z ≈ legacy.packet.position_z atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.x2_x ≈ legacy.packet.x2_x atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.x2_y ≈ legacy.packet.x2_y atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.x2_z ≈ legacy.packet.x2_z atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.weights ≈ legacy.packet.weights atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.gaussian_sum ≈ legacy.packet.gaussian_sum atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.pair_sum ≈ legacy.packet.pair_sum atol = 1.0e-10 rtol = 1.0e-10

    @test audit.full_parent_working_box
    @test audit.missing_row_count == 0
    @test audit.ownership_group_count_min == 1
    @test audit.ownership_group_count_max == 1
    @test audit.ownership_unowned_row_count == 0
    @test audit.ownership_multi_owned_row_count == 0
end

@testset "bond-aligned diatomic source-backed shellification materializer matches active sequence" begin
    fixture = _bond_aligned_diatomic_shellification_plan_fixture()
    source = fixture.source
    plan = fixture.plan
    legacy = source.sequence
    materialized =
        GaussletBases._cartesian_materialize_source_backed_shellification_low_order(
            plan,
            source;
            packet_kernel = :factorized_direct,
            term_coefficients = Float64.(fixture.expansion.coefficients),
        )
    audit = GaussletBases._nested_shell_sequence_contract_audit(
        materialized,
        Tuple(length.(source.geometry.parent_box)),
    )

    @test materialized.working_box == legacy.working_box
    @test materialized.core_indices == legacy.core_indices
    @test materialized.core_states == legacy.core_states
    @test materialized.core_column_range == legacy.core_column_range
    @test materialized.layer_column_ranges == legacy.layer_column_ranges
    @test materialized.support_indices == legacy.support_indices
    @test materialized.support_states == legacy.support_states
    @test length(materialized.shell_layers) == length(source.shared_shell_layers)
    @test all(
        materialized.shell_layers[index] === source.shared_shell_layers[index] for
        index in eachindex(source.shared_shell_layers)
    )
    @test size(materialized.coefficient_matrix) == size(legacy.coefficient_matrix)
    @test materialized.coefficient_matrix ≈ legacy.coefficient_matrix atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.overlap ≈ legacy.packet.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.kinetic ≈ legacy.packet.kinetic atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.position_x ≈ legacy.packet.position_x atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.position_y ≈ legacy.packet.position_y atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.position_z ≈ legacy.packet.position_z atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.x2_x ≈ legacy.packet.x2_x atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.x2_y ≈ legacy.packet.x2_y atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.x2_z ≈ legacy.packet.x2_z atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.weights ≈ legacy.packet.weights atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.gaussian_sum ≈ legacy.packet.gaussian_sum atol = 1.0e-10 rtol = 1.0e-10
    @test materialized.packet.pair_sum ≈ legacy.packet.pair_sum atol = 1.0e-10 rtol = 1.0e-10

    @test audit.full_parent_working_box
    @test audit.missing_row_count == 0
    @test audit.ownership_group_count_min == 1
    @test audit.ownership_group_count_max == 1
    @test audit.ownership_unowned_row_count == 0
    @test audit.ownership_multi_owned_row_count == 0
end

@testset "atom-outward atom-local child subplans lower to oracle child sequences" begin
    fixture = _bond_aligned_diatomic_shellification_plan_fixture(
        xmax_parallel = 12.0,
        xmax_transverse = 8.0,
        min_unsplit_parallel_to_transverse_ratio_for_split = 1.0,
    )
    source = fixture.source
    @test source.geometry.did_split
    @test length(source.geometry.child_boxes) == 2
    @test length(source.child_sequences) == 2

    protect_rows =
        GaussletBases._nested_diatomic_resolve_core_near_nucleus_protect_rows(
            :auto,
            source.nside,
        )
    for (child_index, atom_side) in enumerate((:left, :right))
        oracle = source.child_sequences[child_index]
        child_box = source.geometry.child_boxes[child_index]
        plan = GaussletBases._cartesian_shellification_plan_atom_local_child_low_order(
            fixture.bundles,
            child_box;
            order_index = child_index,
            atom_side = atom_side,
            bond_axis = source.geometry.bond_axis,
            nside = source.nside,
            retention_policy = source.child_shell_retention_contract,
            reference_fudge_factor = 1.2,
            core_near_nucleus_protect_rows = protect_rows,
        )

        @test plan.object_kind == :cartesian_atom_local_child_shellification_plan3d
        @test plan.source_kind == :atom_outward_atom_local_child_shellification_plan
        @test plan.role == :atom_local_subtree
        @test plan.atom_side == atom_side
        @test plan.order_index == child_index
        @test plan.outer_box == child_box
        @test plan.working_box == child_box
        @test plan.bond_axis == source.geometry.bond_axis
        @test plan.nside == source.nside
        @test plan.retention_policy == source.child_shell_retention_contract
        @test plan.core_policy == :diatomic_nonuniform_bond_axis_core
        @test !plan.source_backed
        @test plan.missing_independent_lowering_reason === nothing
        @test plan.retirement_target == :already_plan_lowered_region
        @test plan.coverage_metadata.coverage_complete
        @test plan.coverage_metadata.expected_support_count ==
              GaussletBases._cartesian_shellification_box_point_count(child_box)
        @test plan.coverage_metadata.support_count == length(plan.support_indices)
        @test all(region -> !region.source_backed, plan.shell_regions)
        @test plan.core_region.source_backed == false
        @test plan.core_region.core_policy == :diatomic_nonuniform_bond_axis_core
        @test plan.diagnostics.atom_outward_policy
        @test plan.diagnostics.active_source_oracle_only
        @test !plan.diagnostics.active_source_authority
        @test !plan.diagnostics.route_behavior_changed

        materialized =
            GaussletBases._cartesian_materialize_atom_local_child_shellification_low_order(
                plan,
                fixture.basis,
                fixture.bundles;
                term_coefficients = Float64.(fixture.expansion.coefficients),
                packet_kernel = :factorized_direct,
                build_packet = false,
            )
        materialized_audit = GaussletBases._nested_shell_sequence_contract_audit(
            materialized,
            Tuple(length.(source.geometry.parent_box)),
        )
        oracle_audit = GaussletBases._nested_shell_sequence_contract_audit(
            oracle,
            Tuple(length.(source.geometry.parent_box)),
        )

        @test materialized.working_box == oracle.working_box
        @test materialized.working_box == child_box
        @test materialized.core_indices == oracle.core_indices
        @test materialized.core_states == oracle.core_states
        @test materialized.core_column_range == oracle.core_column_range
        @test materialized.layer_column_ranges == oracle.layer_column_ranges
        @test materialized.support_indices == oracle.support_indices
        @test materialized.support_states === oracle.support_states
        @test size(materialized.coefficient_matrix) == size(oracle.coefficient_matrix)
        @test isapprox(
            materialized.coefficient_matrix,
            oracle.coefficient_matrix;
            atol = 1.0e-10,
            rtol = 1.0e-10,
        )
        @test materialized.packet === oracle.packet
        @test isnothing(materialized.packet)
        @test materialized_audit == oracle_audit

        @test length(materialized.shell_layers) == length(oracle.shell_layers)
        @test length(materialized.shell_layers) == plan.shell_layer_count
        for (shell, oracle_shell, region, column_range) in zip(
            materialized.shell_layers,
            oracle.shell_layers,
            plan.shell_regions,
            materialized.layer_column_ranges,
        )
            @test region.lowering_family == :white_lindsey_adaptive_complete_shell
            @test region.retirement_target == :already_plan_lowered_region
            @test shell.provenance == oracle_shell.provenance
            @test shell.support_indices == oracle_shell.support_indices
            @test shell.support_states == oracle_shell.support_states
            @test shell.face_column_ranges == oracle_shell.face_column_ranges
            @test shell.edge_column_ranges == oracle_shell.edge_column_ranges
            @test shell.corner_column_ranges == oracle_shell.corner_column_ranges
            @test size(shell.coefficient_matrix) == size(oracle_shell.coefficient_matrix)
            @test isapprox(
                shell.coefficient_matrix,
                oracle_shell.coefficient_matrix;
                atol = 1.0e-10,
                rtol = 1.0e-10,
            )
            @test length(column_range) == size(shell.coefficient_matrix, 2)
        end
    end
end

@testset "atom-outward direct midpoint slab region lowers to oracle slab block" begin
    fixture = _bond_aligned_diatomic_shellification_plan_fixture(
        xmax_parallel = 12.0,
        xmax_transverse = 8.0,
        min_unsplit_parallel_to_transverse_ratio_for_split = 1.0,
    )
    source = fixture.source
    @test source.geometry.did_split
    @test !isnothing(source.geometry.shared_midpoint_box)
    @test !isnothing(source.midpoint_slab_column_range)

    midpoint_box = source.geometry.shared_midpoint_box
    region = GaussletBases._cartesian_shellification_plan_direct_midpoint_slab_region3d(
        fixture.bundles,
        midpoint_box;
        order_index = 2,
        bond_axis = source.geometry.bond_axis,
        split_index = source.geometry.split_index,
        column_range = source.midpoint_slab_column_range,
    )

    @test region.object_kind == :cartesian_atom_outward_direct_midpoint_slab_region3d
    @test region.role == :midpoint_slab
    @test region.order == 2
    @test region.order_index == 2
    @test region.box == midpoint_box
    @test region.outer_box == midpoint_box
    @test region.box_shape == Tuple(length.(midpoint_box))
    @test region.region_kind == :direct_contact_or_midpoint_slab
    @test region.lowering_family == :direct_midpoint_slab
    @test region.materialization_dependency == :plan_lowerable_direct_slab
    @test !region.source_backed
    @test region.missing_independent_lowering_reason === nothing
    @test region.retirement_target == :already_plan_lowered_region
    @test region.provenance.source == :atom_outward_direct_midpoint_slab_region
    @test region.provenance.box_source == :geometry_input
    @test region.provenance.column_range_source == :oracle_only_column_range
    @test region.provenance.split_index == source.geometry.split_index
    @test region.provenance.bond_axis == source.geometry.bond_axis
    @test region.column_range == source.midpoint_slab_column_range
    @test region.source_point_count ==
          GaussletBases._cartesian_shellification_box_point_count(midpoint_box)
    @test region.support_count == length(region.support_indices)
    @test region.retained_count == region.source_point_count
    @test region.retained_count == length(source.midpoint_slab_column_range)

    slab_data = GaussletBases._cartesian_materialize_direct_box_region(
        region,
        fixture.bundles,
    )
    expected = GaussletBases._nested_direct_box_coefficients(
        fixture.bundles,
        midpoint_box,
    )
    legacy_block =
        source.sequence.coefficient_matrix[:, source.midpoint_slab_column_range]

    @test slab_data.support_indices == expected.support_indices
    @test slab_data.support_indices == region.support_indices
    @test size(slab_data.coefficient_matrix) == size(expected.coefficient_matrix)
    @test slab_data.coefficient_matrix == expected.coefficient_matrix
    @test size(slab_data.coefficient_matrix) == size(legacy_block)
    @test isapprox(
        slab_data.coefficient_matrix,
        legacy_block;
        atol = 1.0e-10,
        rtol = 1.0e-10,
    )
end

@testset "atom-outward shared complete-shell regions lower to oracle layers" begin
    fixture = _bond_aligned_diatomic_shellification_plan_fixture(
        shared_shell_layer_policy = :complete_rectangular,
    )
    source = fixture.source
    @test !isempty(source.shared_shell_layers)
    @test all(
        layer -> layer isa GaussletBases._CartesianNestedCompleteShell3D,
        source.shared_shell_layers,
    )

    for (layer_index, oracle) in enumerate(source.shared_shell_layers)
        outer_box = oracle.provenance.source_box
        inner_box = oracle.provenance.next_inner_box
        column_range = source.sequence.layer_column_ranges[layer_index]
        region =
            GaussletBases._cartesian_shellification_plan_shared_complete_shell_region3d(
                fixture.bundles,
                outer_box,
                inner_box;
                order_index = layer_index,
                bond_axis = source.geometry.bond_axis,
                nside = source.nside,
                retention_policy = source.shared_shell_retention_contract,
                shared_shell_angular_resolution_scale = 1.4,
                retained_count = size(oracle.coefficient_matrix, 2),
                column_range = column_range,
            )

        @test region.object_kind ==
              :cartesian_atom_outward_shared_complete_shell_region3d
        @test region.role == :regular_shared_molecular_shell
        @test region.order == layer_index
        @test region.order_index == layer_index
        @test region.outer_box == outer_box
        @test region.inner_exclusion_box == inner_box
        @test region.next_inner_box == inner_box
        @test region.region_kind == :complete_rectangular_shell
        @test region.lowering_family == :white_lindsey_adaptive_complete_shell
        @test region.materialization_dependency ==
              :plan_lowerable_shared_complete_shell
        @test !region.source_backed
        @test region.missing_independent_lowering_reason === nothing
        @test region.retirement_target == :already_plan_lowered_region
        @test region.source_point_count == oracle.provenance.source_point_count
        @test region.support_count == length(oracle.support_indices)
        @test region.retained_count == size(oracle.coefficient_matrix, 2)
        @test region.retained_count == length(column_range)
        @test region.column_range == column_range
        @test region.provenance.source == :atom_outward_shared_complete_shell_region
        @test region.provenance.box_source == :geometry_input
        @test region.provenance.retained_count_source == :oracle_only
        @test region.provenance.column_range_source == :oracle_only

        materialized =
            GaussletBases._cartesian_materialize_shared_complete_shell_region(
                region,
                fixture.basis,
                fixture.bundles;
                term_coefficients = Float64.(fixture.expansion.coefficients),
                packet_kernel = :factorized_direct,
            )

        @test materialized.provenance == oracle.provenance
        @test materialized.support_indices == oracle.support_indices
        @test materialized.support_states == oracle.support_states
        @test materialized.face_column_ranges == oracle.face_column_ranges
        @test materialized.edge_column_ranges == oracle.edge_column_ranges
        @test materialized.corner_column_ranges == oracle.corner_column_ranges
        @test size(materialized.coefficient_matrix) == size(oracle.coefficient_matrix)
        @test isapprox(
            materialized.coefficient_matrix,
            oracle.coefficient_matrix;
            atol = 1.0e-10,
            rtol = 1.0e-10,
        )
    end
end

@testset "atom-outward complete-rectangular split assembly matches oracle sequence" begin
    fixture = _bond_aligned_diatomic_shellification_plan_fixture(
        shared_shell_layer_policy = :complete_rectangular,
        xmax_parallel = 12.0,
        xmax_transverse = 8.0,
        min_unsplit_parallel_to_transverse_ratio_for_split = 1.0,
    )
    source = fixture.source
    oracle = source.sequence
    @test source.geometry.did_split
    @test length(source.geometry.child_boxes) == 2
    @test length(source.child_sequences) == 2
    @test !isempty(source.shared_shell_layers)
    @test all(
        layer -> layer isa GaussletBases._CartesianNestedCompleteShell3D,
        source.shared_shell_layers,
    )

    protect_rows =
        GaussletBases._nested_diatomic_resolve_core_near_nucleus_protect_rows(
            :auto,
            source.nside,
        )
    left_child_plan =
        GaussletBases._cartesian_shellification_plan_atom_local_child_low_order(
            fixture.bundles,
            source.geometry.child_boxes[1];
            order_index = 1,
            atom_side = :left,
            bond_axis = source.geometry.bond_axis,
            nside = source.nside,
            retention_policy = source.child_shell_retention_contract,
            reference_fudge_factor = 1.2,
            core_near_nucleus_protect_rows = protect_rows,
        )
    right_child_plan =
        GaussletBases._cartesian_shellification_plan_atom_local_child_low_order(
            fixture.bundles,
            source.geometry.child_boxes[2];
            order_index = 3,
            atom_side = :right,
            bond_axis = source.geometry.bond_axis,
            nside = source.nside,
            retention_policy = source.child_shell_retention_contract,
            reference_fudge_factor = 1.2,
            core_near_nucleus_protect_rows = protect_rows,
        )
    midpoint_region =
        isnothing(source.geometry.shared_midpoint_box) ?
        nothing :
        GaussletBases._cartesian_shellification_plan_direct_midpoint_slab_region3d(
            fixture.bundles,
            source.geometry.shared_midpoint_box;
            order_index = 2,
            bond_axis = source.geometry.bond_axis,
            split_index = source.geometry.split_index,
            column_range = source.midpoint_slab_column_range,
        )
    shared_regions = Tuple(
        GaussletBases._cartesian_shellification_plan_shared_complete_shell_region3d(
            fixture.bundles,
            layer.provenance.source_box,
            layer.provenance.next_inner_box;
            order_index = layer_index,
            bond_axis = source.geometry.bond_axis,
            nside = source.nside,
            retention_policy = source.shared_shell_retention_contract,
            shared_shell_angular_resolution_scale = 1.4,
            retained_count = size(layer.coefficient_matrix, 2),
            column_range = source.sequence.layer_column_ranges[layer_index],
        ) for (layer_index, layer) in enumerate(source.shared_shell_layers)
    )

    materialization =
        GaussletBases._cartesian_materialize_split_complete_rectangular_shellification_low_order(
            left_child_plan,
            right_child_plan,
            shared_regions,
            fixture.basis,
            fixture.bundles;
            midpoint_region,
            term_coefficients = Float64.(fixture.expansion.coefficients),
            packet_kernel = :factorized_direct,
        )
    materialized = materialization.sequence

    @test materialization.object_kind ==
          :cartesian_atom_outward_split_complete_rectangular_materialization
    @test materialization.private_development_only
    @test !materialization.active_source_authority
    @test !materialization.route_behavior_changed
    @test !(:source in propertynames(materialization))
    @test materialization.child_column_ranges == Tuple(source.child_column_ranges)
    @test materialization.midpoint_slab_column_range ==
          source.midpoint_slab_column_range

    @test materialized.working_box == oracle.working_box
    @test materialized.core_indices == oracle.core_indices
    @test materialized.core_states == oracle.core_states
    @test materialized.core_column_range == oracle.core_column_range
    @test materialized.layer_column_ranges == oracle.layer_column_ranges
    @test materialized.support_indices == oracle.support_indices
    @test materialized.support_states == oracle.support_states
    @test length(materialized.shell_layers) == length(oracle.shell_layers)
    @test length(materialized.shell_layers) == length(source.shared_shell_layers)
    @test size(materialized.coefficient_matrix) == size(oracle.coefficient_matrix)
    @test isapprox(
        materialized.coefficient_matrix,
        oracle.coefficient_matrix;
        atol = 1.0e-10,
        rtol = 1.0e-10,
    )

    for (layer, oracle_layer, region, column_range) in zip(
        materialization.shared_shell_layers,
        source.shared_shell_layers,
        shared_regions,
        materialized.layer_column_ranges,
    )
        @test region.column_range == column_range
        @test layer.provenance == oracle_layer.provenance
        @test layer.support_indices == oracle_layer.support_indices
        @test layer.support_states == oracle_layer.support_states
        @test layer.face_column_ranges == oracle_layer.face_column_ranges
        @test layer.edge_column_ranges == oracle_layer.edge_column_ranges
        @test layer.corner_column_ranges == oracle_layer.corner_column_ranges
        @test size(layer.coefficient_matrix) == size(oracle_layer.coefficient_matrix)
        @test isapprox(
            layer.coefficient_matrix,
            oracle_layer.coefficient_matrix;
            atol = 1.0e-10,
            rtol = 1.0e-10,
        )
    end
end

@testset "atom-growth-backed complete-rectangular scaffold maps construction plan" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    basis = bond_aligned_homonuclear_qw_basis(
        family = :G10,
        bond_length = 5.0,
        core_spacing = 0.7,
        xmax_parallel = 8.0,
        xmax_transverse = 4.0,
        bond_axis = :z,
    )
    bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(basis, expansion)
    anatomy = GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(
        basis,
        bundles;
        protected_atom_side_count = 5,
    )
    construction_plan =
        GaussletBases._nested_bond_aligned_diatomic_atom_growth_construction_plan(
            anatomy,
        )
    nside = 5
    retention = GaussletBases._nested_resolve_complete_shell_retention(nside)
    protect_rows =
        GaussletBases._nested_diatomic_resolve_core_near_nucleus_protect_rows(
            :auto,
            nside,
        )

    scaffold =
        GaussletBases._cartesian_shellification_plan_atom_growth_complete_rectangular_low_order(
            construction_plan,
            bundles;
            nside,
            child_retention_policy = retention,
            shared_retention_policy = retention,
            reference_fudge_factor = 1.2,
            core_near_nucleus_protect_rows = protect_rows,
            shared_shell_angular_resolution_scale = 1.4,
        )

    @test construction_plan.region_order == [
        :outer_mismatch_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :left_atom_box,
        :right_atom_box,
        :contact_cap,
    ]
    @test scaffold.object_kind == :cartesian_atom_growth_shellification_plan3d
    @test scaffold.source_kind ==
          :bond_aligned_diatomic_atom_growth_construction_plan
    @test scaffold.spatial_policy_order == :atom_outward
    @test scaffold.construction_region_order == Tuple(construction_plan.region_order)
    @test scaffold.assembly_shell_order == :outside_in
    @test scaffold.assembly_core_order ==
          (
              :outer_mismatch_shared_molecular_shell,
              :left_atom_box,
              :contact_cap,
              :right_atom_box,
          )
    @test scaffold.ordered_region_roles == Tuple(construction_plan.region_order)
    @test [region.support_count for region in scaffold.regions] ==
          [length(region.support_indices) for region in construction_plan.regions]
    @test all(region -> !region.source_backed, scaffold.regions)
    @test !(:source in propertynames(scaffold))
    @test scaffold.diagnostics.atom_growth_construction_plan_authority
    @test !scaffold.diagnostics.active_source_authority
    @test !scaffold.diagnostics.active_source_oracle_comparison_run
    @test !scaffold.diagnostics.materialization_behavior_changed
    @test !scaffold.diagnostics.public_default_behavior_changed

    counts = scaffold.materialization_dependency_counts
    @test counts.source_backed_region_count == 0
    @test counts.atom_local_child_plan_count == 2
    @test counts.direct_contact_slab_count == 1
    @test counts.shared_complete_shell_count == 1
    @test counts.outer_mismatch_boundary_slab_set_count == 1
    @test counts.unsupported_outer_mismatch_count == 0
    @test scaffold.unsupported_region_count == 0
    @test scaffold.materialization_status == :ready_supported_complete_rectangular_subset

    left_region = only(
        region for region in construction_plan.regions if region.role == :left_atom_box
    )
    right_region = only(
        region for region in construction_plan.regions if region.role == :right_atom_box
    )
    contact_region = only(
        region for region in construction_plan.regions if region.role == :contact_cap
    )
    shared_region = only(
        region for region in construction_plan.regions
        if region.role == :regular_shared_molecular_shell
    )
    mismatch_construction_region = only(
        region for region in construction_plan.regions
        if region.role == :outer_mismatch_shared_molecular_shell
    )
    mismatch_region = only(
        region for region in scaffold.regions
        if region.role == :outer_mismatch_shared_molecular_shell
    )

    @test scaffold.left_child_plan.outer_box == left_region.box
    @test scaffold.left_child_plan.support_count == length(left_region.support_indices)
    @test scaffold.left_child_plan.source_backed == false
    @test scaffold.right_child_plan.outer_box == right_region.box
    @test scaffold.right_child_plan.support_count == length(right_region.support_indices)
    @test scaffold.right_child_plan.source_backed == false
    @test scaffold.contact_cap_region.box == contact_region.box
    @test scaffold.contact_cap_region.support_count ==
          length(contact_region.support_indices)
    @test scaffold.contact_cap_region.source_backed == false
    @test length(scaffold.shared_complete_shell_regions) == 1
    @test scaffold.shared_complete_shell_regions[1].outer_box == shared_region.box
    @test scaffold.shared_complete_shell_regions[1].inner_exclusion_box ==
          shared_region.inner_exclusion_box
    @test scaffold.shared_complete_shell_regions[1].support_count ==
          length(shared_region.support_indices)
    @test scaffold.shared_complete_shell_regions[1].source_backed == false

    for region in scaffold.regions
        @test region.independently_lowerable
        @test region.missing_independent_lowering_reason === nothing
        @test region.retirement_target == :already_plan_lowered_region
        @test region.support_count_matches_lowering_piece
    end
    @test mismatch_region.support_count == 98
    @test mismatch_region.materialization_dependency ==
          :plan_lowerable_outer_mismatch_boundary_slab_set
    @test mismatch_region.lowering_piece_object_kind ==
          :cartesian_outer_mismatch_boundary_slab_set3d
    slab_set = only(scaffold.outer_mismatch_boundary_slab_sets)
    @test slab_set.role == :outer_mismatch_shared_molecular_shell
    @test slab_set.source_authority ==
          :_BondAlignedDiatomicAtomGrowthConstructionPlan3D
    @test slab_set.source_backed == false
    @test slab_set.independently_lowerable
    @test slab_set.support_count == length(mismatch_construction_region.support_indices)
    @test slab_set.retained_count == slab_set.support_count
    @test slab_set.support_count_matches_retained_count
    @test slab_set.slab_piece_roles ==
          (:outer_mismatch_z_low_slab, :outer_mismatch_z_high_slab)
    @test map(piece -> piece.support_count, slab_set.slab_pieces) == (49, 49)
    @test slab_set.slab_column_ranges == (1:49, 50:98)
    @test slab_set.aggregate_support_coverage.coverage_ok
    @test slab_set.aggregate_support_coverage.duplicate_count == 0
    @test slab_set.aggregate_support_coverage.missing_count == 0
    @test slab_set.aggregate_support_coverage.outside_count == 0
    piece_support_indices = reduce(
        vcat,
        (piece.support_indices for piece in slab_set.slab_pieces);
        init = Int[],
    )
    @test length(unique(piece_support_indices)) == length(piece_support_indices)
    @test sort(piece_support_indices) ==
          sort(mismatch_construction_region.support_indices)

    slab_materialization =
        GaussletBases._cartesian_materialize_outer_mismatch_boundary_slab_set(
            slab_set,
            bundles,
        )
    parent_dimension = prod(length.(construction_plan.anatomy.recipe.parent_box))
    @test slab_materialization.status ==
          :materialized_outer_mismatch_boundary_slab_set
    @test slab_materialization.support_indices == piece_support_indices
    @test size(slab_materialization.coefficient_matrix) ==
          (parent_dimension, slab_set.retained_count)
    local_direct_block =
        Matrix(slab_materialization.coefficient_matrix[piece_support_indices, :])
    @test local_direct_block == [
        row == column ? 1.0 : 0.0 for row in 1:slab_set.retained_count,
        column in 1:slab_set.retained_count
    ]

    materialization =
        GaussletBases._cartesian_materialize_atom_growth_complete_rectangular_shellification_low_order(
            scaffold,
            basis,
            bundles,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :factorized_direct,
        )
    @test materialization.status == :materialized_supported_complete_rectangular_low_order
    @test materialization.blocked_reason === nothing
    @test materialization.unsupported_region_count == 0
    @test length(materialization.outer_mismatch_boundary_slab_sets) == 1
    @test !materialization.active_source_authority
    @test !materialization.route_behavior_changed
    assembly = materialization.assembly
    sequence = materialization.sequence
    audit = GaussletBases._nested_shell_sequence_contract_audit(
        sequence,
        Tuple(length.(construction_plan.anatomy.recipe.parent_box)),
    )
    @test assembly.status == :materialized_atom_growth_complete_rectangular_low_order
    @test assembly.core_block_order == scaffold.assembly_core_order
    expected_core_support = vcat(
        piece_support_indices,
        scaffold.left_child_plan.support_indices,
        scaffold.contact_cap_region.support_indices,
        scaffold.right_child_plan.support_indices,
    )
    @test assembly.core_support_indices == expected_core_support
    @test sequence.core_indices == expected_core_support
    @test sequence.working_box == construction_plan.anatomy.recipe.parent_box
    @test audit.full_parent_working_box
    @test audit.missing_row_count == 0
    @test audit.ownership_unowned_row_count == 0
    @test audit.ownership_multi_owned_row_count == 0
    @test assembly.outer_mismatch_column_ranges == (1:98,)
    @test assembly.outer_mismatch_slab_column_ranges == ((1:49, 50:98),)
    contact_start = last(assembly.child_column_ranges[1]) + 1
    @test assembly.contact_cap_column_range ==
          contact_start:(contact_start + length(contact_region.support_indices) - 1)
    @test assembly.shared_shell_column_ranges == Tuple(sequence.layer_column_ranges)
    @test length(assembly.shared_shell_layers) ==
          length(scaffold.shared_complete_shell_regions)
    direct_outer_columns = only(assembly.outer_mismatch_column_ranges)
    local_outer_block =
        Matrix(sequence.coefficient_matrix[piece_support_indices, direct_outer_columns])
    @test local_outer_block == local_direct_block
end

@testset "bond-aligned diatomic low-order shellification plan matches active source" begin
    fixture = _bond_aligned_diatomic_shellification_plan_fixture()
    nside = fixture.nside
    source = fixture.source
    sequence = source.sequence
    specific_plan = fixture.specific_plan
    plan = fixture.plan
    audit = fixture.audit
    geometry = source.geometry

    @test plan == specific_plan
    @test plan.object_kind == :cartesian_shellification_plan3d
    @test plan.status == :planned_metadata_only
    @test plan.private_development_only
    @test plan.source_kind == :bond_aligned_diatomic_active_source_shellification_plan
    @test plan.route_family == :white_lindsey_low_order
    @test plan.system_classification == :bond_aligned_diatomic
    @test plan.shellification_stage == :route_neutral_spatial_planning
    @test plan.lowering_stage == :not_lowered_by_shellification_plan
    @test plan.parent_box == geometry.parent_box
    @test plan.working_box == sequence.working_box
    @test plan.split_working_box == geometry.working_box
    @test plan.bond_axis == geometry.bond_axis == :z
    @test plan.nside == nside
    @test plan.split_status == (geometry.did_split ? :split : :no_split)
    @test plan.did_split == geometry.did_split
    @test plan.shared_shell_layer_count == length(source.shared_shell_layers)
    @test plan.child_sequence_count == length(source.child_sequences)
    @test plan.child_column_ranges == Tuple(source.child_column_ranges)
    @test plan.midpoint_slab_present == !isnothing(source.midpoint_slab_column_range)
    @test plan.midpoint_slab_column_range == source.midpoint_slab_column_range
    @test plan.merged_sequence_core_column_range == sequence.core_column_range
    @test plan.merged_sequence_shell_layer_column_ranges ==
          Tuple(sequence.layer_column_ranges)
    @test plan.merged_sequence_shell_layer_count == length(sequence.shell_layers)
    @test plan.merged_sequence_retained_dimension == size(sequence.coefficient_matrix, 2)
    @test plan.merged_sequence_support_count == length(sequence.support_indices)
    @test plan.retained_dimension == size(sequence.coefficient_matrix, 2)
    @test plan.support_count == length(sequence.support_indices)
    @test plan.region_count == length(plan.regions)
    @test plan.ordered_region_roles == Tuple(region.role for region in plan.regions)
    @test plan.ordered_region_boxes == Tuple(region.box for region in plan.regions)
    @test plan.ordered_materialization_dependencies ==
          Tuple(region.materialization_dependency for region in plan.regions)
    @test plan.materialization_dependency_counts.plan_lowerable_region_count == 0
    @test plan.materialization_dependency_counts.source_backed_region_count ==
          plan.region_count
    @test plan.materialization_dependency_counts.source_backed_shared_shell_layer_count ==
          plan.shared_shell_region_count
    @test plan.materialization_dependency_counts.source_backed_child_sequence_count ==
          plan.child_subtree_region_count
    @test plan.materialization_dependency_counts.source_box_direct_adapter_region_count ==
          plan.midpoint_slab_region_count

    shared_regions =
        Tuple(region for region in plan.regions if region.role == :shared_outer_shell)
    @test length(shared_regions) == length(source.shared_shell_layers)
    @test length(shared_regions) == length(sequence.layer_column_ranges)
    @test all(
        sequence.shell_layers[index] === source.shared_shell_layers[index] for
        index in eachindex(source.shared_shell_layers)
    )
    for (index, (region, layer, column_range)) in
        enumerate(zip(shared_regions, source.shared_shell_layers, sequence.layer_column_ranges))
        layer_provenance = layer.provenance
        source_box =
            hasproperty(layer_provenance, :source_box) ?
            layer_provenance.source_box :
            layer_provenance.current_box
        next_inner_box =
            hasproperty(layer_provenance, :next_inner_box) ?
            layer_provenance.next_inner_box :
            layer_provenance.inner_box
        @test region.object_kind == :cartesian_shellification_region3d
        @test region.order == index
        @test region.lowering_status == :planned_not_lowered
        @test region.materialization_dependency == :source_backed_shared_shell_layer
        @test region.box == source_box
        @test region.next_inner_box == next_inner_box
        @test region.box_shape == Tuple(length.(source_box))
        @test region.support_count == length(layer.support_indices)
        @test region.retained_count == size(layer.coefficient_matrix, 2)
        @test region.retained_count == length(column_range)
        @test region.column_range == column_range
        @test region.provenance.shared_layer_index == index
        @test region.provenance.layer_kind == Symbol(nameof(typeof(layer)))
        @test region.provenance.sequence_column_range == column_range
        if layer isa GaussletBases._CartesianNestedEndcapPanelShellLayer3D
            @test region.lowering_family == :white_lindsey_endcap_panel_owned_shell
            @test region.source_point_count == layer.owned_units.audit.expected_support_count
            @test region.provenance.owned_unit_summary.coverage_ok
            @test region.provenance.owned_unit_summary.retained_count ==
                  size(layer.coefficient_matrix, 2)
        else
            @test region.lowering_family == :white_lindsey_complete_shell
            @test region.source_point_count == layer.provenance.source_point_count
        end
    end

    child_regions = Tuple(
        region for region in plan.regions
        if region.role in (:atom_local_subtree, :unsplit_child_subtree)
    )
    @test length(child_regions) == length(source.child_sequences)
    @test plan.child_subtree_region_count == length(source.child_sequences)
    for (child_index, (region, child_sequence, column_range)) in enumerate(
        zip(child_regions, source.child_sequences, source.child_column_ranges),
    )
        @test region.object_kind == :cartesian_shellification_region3d
        @test region.lowering_status == :planned_not_lowered
        @test region.lowering_family == :white_lindsey_child_shell_sequence
        @test region.materialization_dependency == :source_backed_child_sequence
        @test region.role == (geometry.did_split ? :atom_local_subtree : :unsplit_child_subtree)
        @test region.box == child_sequence.working_box
        @test region.box_shape == Tuple(length.(child_sequence.working_box))
        @test region.source_point_count == length(child_sequence.support_indices)
        @test region.support_count == length(child_sequence.support_indices)
        @test region.retained_count == size(child_sequence.coefficient_matrix, 2)
        @test region.retained_count == length(column_range)
        @test region.column_range == column_range
        @test region.provenance.child_index == child_index
        @test region.provenance.child_sequence_working_box == child_sequence.working_box
        @test region.provenance.sequence_core_column_range ==
              child_sequence.core_column_range
        @test region.provenance.merged_sequence_column_range == column_range
        if geometry.did_split
            @test region.provenance.geometry_child_box == geometry.child_boxes[child_index]
            @test region.provenance.geometry_child_box_matches_sequence
        else
            @test isnothing(region.provenance.geometry_child_box)
            @test region.provenance.geometry_child_box_matches_sequence
        end
    end

    midpoint_regions =
        Tuple(region for region in plan.regions if region.role == :midpoint_slab)
    if isnothing(source.midpoint_slab_column_range)
        @test isempty(midpoint_regions)
        @test plan.midpoint_slab_region_count == 0
        @test plan.shared_midpoint_box === nothing
    else
        region = only(midpoint_regions)
        midpoint_box = source.geometry.shared_midpoint_box
        @test plan.midpoint_slab_region_count == 1
        @test region.box == midpoint_box
        @test region.lowering_family == :direct_midpoint_slab
        @test region.materialization_dependency ==
              :source_box_direct_in_source_backed_adapter
        @test region.column_range == source.midpoint_slab_column_range
        @test region.source_point_count ==
              GaussletBases._cartesian_shellification_box_point_count(midpoint_box)
        @test region.retained_count == length(source.midpoint_slab_column_range)
    end

    @test plan.coverage.object_kind == :cartesian_shellification_plan_coverage3d
    @test plan.coverage.status == :coverage_complete
    @test plan.coverage.coverage_complete
    @test plan.coverage.parent_point_count ==
          GaussletBases._cartesian_shellification_box_point_count(geometry.parent_box)
    @test plan.coverage.expected_support_count == audit.expected_support_count
    @test plan.coverage.support_count == audit.support_count
    @test plan.coverage.region_support_count == audit.support_count
    @test plan.coverage.region_support_count_matches_sequence
    @test plan.coverage.missing_row_count == 0
    @test plan.coverage.ownership_group_count_min == 1
    @test plan.coverage.ownership_group_count_max == 1
    @test plan.coverage.ownership_unowned_row_count == 0
    @test plan.coverage.ownership_multi_owned_row_count == 0
    @test audit.missing_row_count == 0
    @test audit.ownership_group_count_min == 1
    @test audit.ownership_group_count_max == 1
    @test audit.ownership_unowned_row_count == 0
    @test audit.ownership_multi_owned_row_count == 0

    @test plan.diagnostics.route_neutral_spatial_planning
    @test plan.diagnostics.describes_existing_active_source
    @test !plan.diagnostics.lowering_applied_by_plan
    @test !plan.diagnostics.materialization_behavior_changed
    @test !plan.diagnostics.public_default_behavior_changed
    @test !plan.diagnostics.shellification_rewrite
    @test !plan.diagnostics.pqs_production_source_box_materialization_claimed
    @test !plan.diagnostics.mwg_ida_semantics_changed

    plan_summary = GaussletBases._cartesian_shellification_plan_private_summary(plan)
    @test plan_summary.ordered_materialization_dependencies ==
          plan.ordered_materialization_dependencies
    @test plan_summary.materialization_dependency_counts ==
          plan.materialization_dependency_counts
    @test plan_summary.plan_lowerable_region_count == 0
    @test plan_summary.source_backed_region_count == plan.region_count
    @test plan_summary.source_box_direct_adapter_region_count ==
          plan.midpoint_slab_region_count
end
