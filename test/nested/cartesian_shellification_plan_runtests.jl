using Test
using GaussletBases

function _one_center_shellification_plan_fixture(; count::Int = 9, nside::Int = 5)
    expansion = coulomb_gaussian_expansion(doacc = false)
    basis = build_basis(MappedUniformBasisSpec(:G10;
        count,
        mapping = white_lindsey_atomic_mapping(; Z = 4.0, d = 0.15),
        reference_spacing = 1.0,
    ))
    sequence = build_one_center_atomic_full_parent_shell_sequence(
        basis;
        exponents = expansion.exponents,
        gausslet_backend = :numerical_reference,
        refinement_levels = 0,
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
    return (; count, nside, sequence, plan, audit, diagnostics)
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
