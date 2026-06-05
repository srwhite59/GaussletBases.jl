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
    plan = GaussletBases._cartesian_shellification_plan_one_center_low_order(
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
    return (; count, nside, basis, bundle, expansion, sequence, plan, audit, diagnostics)
end

function _bond_aligned_diatomic_shellification_plan_fixture(;
    nside::Int = 5,
    shared_shell_layer_policy::Symbol = :endcap_panel_owned,
)
    basis = bond_aligned_homonuclear_qw_basis(
        family = :G10,
        bond_length = 1.4,
        core_spacing = 0.7,
        xmax_parallel = 6.0,
        xmax_transverse = 4.0,
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
        shared_shell_layer_policy,
        shared_shell_endcap_panel_q = 4,
        shared_shell_endcap_panel_L = 4,
    )
    plan = GaussletBases._cartesian_shellification_plan_bond_aligned_diatomic_low_order(
        source,
    )
    audit = GaussletBases._nested_source_contract_audit(source)
    return (; nside, basis, bundles, expansion, source, plan, audit)
end

@testset "one-center low-order shellification plan matches legacy sequence" begin
    fixture = _one_center_shellification_plan_fixture()
    count = fixture.count
    nside = fixture.nside
    sequence = fixture.sequence
    plan = fixture.plan
    audit = fixture.audit
    diagnostics = fixture.diagnostics

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

    shell_regions =
        Tuple(region for region in plan.regions if region.role == :low_order_complete_shell)
    @test length(shell_regions) == length(sequence.shell_layers)
    for (region, shell, column_range) in
        zip(shell_regions, sequence.shell_layers, sequence.layer_column_ranges)
        @test region.object_kind == :cartesian_shellification_region3d
        @test region.lowering_family == :white_lindsey_complete_shell
        @test region.lowering_status == :planned_not_lowered
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

@testset "bond-aligned diatomic low-order shellification plan matches active source" begin
    fixture = _bond_aligned_diatomic_shellification_plan_fixture()
    nside = fixture.nside
    source = fixture.source
    sequence = source.sequence
    plan = fixture.plan
    audit = fixture.audit
    geometry = source.geometry

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
end
