@testset "Bond-aligned diatomic endcap-panel shared shell source policy" begin
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
    term_coefficients = Float64.(expansion.coefficients)

    default_source = GaussletBases._nested_bond_aligned_diatomic_source(
        basis,
        bundles;
        bond_axis = :z,
        nside = 5,
        term_coefficients,
        packet_kernel = :factorized_direct,
    )
    endcap_source = GaussletBases._nested_bond_aligned_diatomic_source(
        basis,
        bundles;
        bond_axis = :z,
        nside = 5,
        term_coefficients,
        packet_kernel = :factorized_direct,
        shared_shell_layer_policy = :endcap_panel_owned,
        shared_shell_endcap_panel_q = 4,
        shared_shell_endcap_panel_L = 4,
    )

    @test length(default_source.shared_shell_layers) == 1
    @test all(layer isa GaussletBases._CartesianNestedCompleteShell3D for layer in default_source.shared_shell_layers)
    @test length(endcap_source.shared_shell_layers) == length(default_source.shared_shell_layers)
    @test only(endcap_source.shared_shell_layers) isa GaussletBases._CartesianNestedEndcapPanelShellLayer3D
    @test length(endcap_source.child_sequences) == length(default_source.child_sequences) == 1
    @test size(only(endcap_source.child_sequences).coefficient_matrix, 2) ==
        size(only(default_source.child_sequences).coefficient_matrix, 2)

    default_shared_columns = [size(layer.coefficient_matrix, 2) for layer in default_source.shared_shell_layers]
    endcap_shared_columns = [size(layer.coefficient_matrix, 2) for layer in endcap_source.shared_shell_layers]
    @test default_shared_columns == [130]
    @test endcap_shared_columns == [96]
    @test size(endcap_source.sequence.coefficient_matrix, 2) ==
        size(default_source.sequence.coefficient_matrix, 2) - sum(default_shared_columns) +
        sum(endcap_shared_columns)

    layer = only(endcap_source.shared_shell_layers)
    @test layer.provenance.support_contract == :thin_endcap_box_perimeter
    @test layer.provenance.coefficient_contract == :product_doside
    @test layer.provenance.q == 4
    @test layer.provenance.L == 4
    @test layer.provenance.packet_kernel == :factorized_direct
    @test layer.owned_units.audit.coverage_ok
    @test layer.owned_units.audit.expected_support_count == 314
    @test layer.owned_units.audit.owned_support_count == 314
    @test layer.owned_units.audit.duplicate_count == 0
    @test layer.owned_units.audit.missing_count == 0
    @test layer.owned_units.audit.outside_count == 0
    @test layer.owned_units.audit.retained_count == 96
    @test length(layer.support_indices) == 314
    @test size(layer.coefficient_matrix) == (prod(GaussletBases._nested_axis_lengths(bundles)), 96)
    @test all(isfinite, layer.packet.overlap)
    @test norm(layer.packet.overlap - I, Inf) < 1.0e-8
    CCP = GaussletBases.CartesianContractedParents
    CP = GaussletBases.CartesianParentGaussletBases
    CCPM = GaussletBases.CartesianContractedParentMetrics
    pgdg_x = GaussletBases._nested_axis_pgdg(bundles, :x)
    pgdg_y = GaussletBases._nested_axis_pgdg(bundles, :y)
    pgdg_z = GaussletBases._nested_axis_pgdg(bundles, :z)
    axis_metrics = (
        x = (
            overlap = pgdg_x.overlap,
            position = pgdg_x.position,
            x2 = pgdg_x.x2,
            weights = pgdg_x.weights,
            centers = pgdg_x.centers,
            source = :nested_pgdg_axis,
        ),
        y = (
            overlap = pgdg_y.overlap,
            position = pgdg_y.position,
            x2 = pgdg_y.x2,
            weights = pgdg_y.weights,
            centers = pgdg_y.centers,
            source = :nested_pgdg_axis,
        ),
        z = (
            overlap = pgdg_z.overlap,
            position = pgdg_z.position,
            x2 = pgdg_z.x2,
            weights = pgdg_z.weights,
            centers = pgdg_z.centers,
            source = :nested_pgdg_axis,
        ),
    )
    safe_axis_data = (
        x = merge(axis_metrics.x, (kinetic = pgdg_x.kinetic,)),
        y = merge(axis_metrics.y, (kinetic = pgdg_y.kinetic,)),
        z = merge(axis_metrics.z, (kinetic = pgdg_z.kinetic,)),
    )
    dims = GaussletBases._nested_axis_lengths(bundles)
    pre_packet_source = CCP._cartesian_endcap_panel_pre_packet_build_source(
        layer.owned_units,
        layer.coefficient_matrix,
        layer.unit_column_ranges,
        layer.support_indices,
        dims,
    )
    post_layer_units = [
        GaussletBases._nested_product_staged_unit_from_owned_unit(
            owned_unit;
            column_range,
            dims,
        ) for (owned_unit, column_range) in
            zip(layer.owned_units.units, layer.unit_column_ranges)
    ]
    post_layer_sidecar = GaussletBases._CartesianNestedProductStagedByCenterSidecar3D(
        dims,
        post_layer_units,
        (; source = :post_layer_test_sidecar, support_contract = :product_owned_units),
        (
            parent_dimension = prod(dims),
            final_dimension = size(layer.coefficient_matrix, 2),
            unit_count = length(post_layer_units),
            product_unit_count = length(post_layer_units),
            generic_unit_count = 0,
            support_counts = Int[unit.diagnostics.support_count for unit in post_layer_units],
            max_support_count = maximum(unit.diagnostics.support_count for unit in post_layer_units),
        ),
    )
    post_layer_source = CCP._cartesian_packet_build_source(post_layer_sidecar)
    @test pre_packet_source.parent_dimension == post_layer_source.parent_dimension == prod(dims)
    @test pre_packet_source.contracted_dimension == post_layer_source.contracted_dimension == 96
    @test pre_packet_source.payload_kind_counts == post_layer_source.payload_kind_counts
    @test Dict(pre_packet_source.payload_kind_counts)[:product_doside] == 6
    @test [payload.payload_kind for payload in pre_packet_source.resolved_payloads] ==
          [payload.payload_kind for payload in post_layer_source.resolved_payloads]
    @test [payload.column_range for payload in pre_packet_source.resolved_payloads] ==
          [payload.column_range for payload in post_layer_source.resolved_payloads]
    @test [payload.support_indices for payload in pre_packet_source.resolved_payloads] ==
          [payload.support_indices for payload in post_layer_source.resolved_payloads]
    @test pre_packet_source.candidate_packet_fields == post_layer_source.candidate_packet_fields
    @test pre_packet_source.missing_packet_fields == post_layer_source.missing_packet_fields
    @test pre_packet_source.diagnostics.packet_construction_consumes_source == false
    @test pre_packet_source.diagnostics.source_object_builds_packet_matrices == false
    @test pre_packet_source.diagnostics.nested_shell_packet_remains_authoritative
    @test !pre_packet_source.diagnostics.numerical_packet_matrices_built
    @test !pre_packet_source.diagnostics.operator_data_available
    @test !pre_packet_source.diagnostics.packet_operator_data_checked
    pre_packet_shadow = CCPM._cartesian_packet_build_source_safe_field_shadow(
        pre_packet_source,
        safe_axis_data,
    )
    layer_parent = CP.cartesian_parent_gausslet_basis(basis)
    layer_contracted_parent = CCP.cartesian_contracted_parent(
        layer_parent,
        layer.coefficient_matrix;
        units = CCP._contracted_parent_units_from_staged_sidecar(post_layer_sidecar),
        metadata = (; source = :endcap_panel_layer_pre_packet_shadow_test),
    )
    layer_metric_packet = CCPM.cartesian_contracted_parent_metric_packet(
        layer_contracted_parent;
        axis_metrics,
        construction_path = :product_staged_metric_contraction,
    )
    @test pre_packet_shadow.diagnostics.source_driven_shadow_only
    @test !pre_packet_shadow.diagnostics.construction_adoption
    @test pre_packet_shadow.diagnostics.current_builder_authority ==
          :nested_shell_packet
    @test pre_packet_shadow.diagnostics.product_unit_count == 6
    @test pre_packet_shadow.diagnostics.support_dense_unit_count == 0
    @test pre_packet_shadow.overlap ≈ layer.packet.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test pre_packet_shadow.position_x ≈ layer.packet.position_x atol = 1.0e-10 rtol = 1.0e-10
    @test pre_packet_shadow.position_y ≈ layer.packet.position_y atol = 1.0e-10 rtol = 1.0e-10
    @test pre_packet_shadow.position_z ≈ layer.packet.position_z atol = 1.0e-10 rtol = 1.0e-10
    @test pre_packet_shadow.x2_x ≈ layer.packet.x2_x atol = 1.0e-10 rtol = 1.0e-10
    @test pre_packet_shadow.x2_y ≈ layer.packet.x2_y atol = 1.0e-10 rtol = 1.0e-10
    @test pre_packet_shadow.x2_z ≈ layer.packet.x2_z atol = 1.0e-10 rtol = 1.0e-10
    @test pre_packet_shadow.weights ≈ layer.packet.weights atol = 1.0e-10 rtol = 1.0e-10
    @test pre_packet_shadow.kinetic ≈ layer.packet.kinetic atol = 1.0e-10 rtol = 1.0e-10
    @test pre_packet_shadow.first_moments ≈ layer_metric_packet.first_moments atol = 1.0e-10 rtol = 1.0e-10
    sequence_core_coefficients =
        endcap_source.sequence.coefficient_matrix[:, endcap_source.sequence.core_column_range]
    pre_sequence_source = CCP._cartesian_nested_sequence_pre_packet_build_source(
        dims,
        endcap_source.sequence.core_indices,
        sequence_core_coefficients,
        endcap_source.sequence.core_column_range,
        endcap_source.sequence.shell_layers,
        endcap_source.sequence.layer_column_ranges,
        endcap_source.sequence.coefficient_matrix,
        endcap_source.sequence.support_indices,
    )
    @test pre_sequence_source.parent_dimension == prod(dims)
    @test pre_sequence_source.contracted_dimension == size(endcap_source.sequence.coefficient_matrix, 2)
    @test Dict(pre_sequence_source.payload_kind_counts)[:product_doside] == 6
    @test Dict(pre_sequence_source.payload_kind_counts)[:support_dense] == 1
    @test pre_sequence_source.diagnostics.packet_construction_consumes_source == false
    @test pre_sequence_source.diagnostics.source_object_builds_packet_matrices == false
    @test pre_sequence_source.diagnostics.nested_shell_packet_remains_authoritative
    @test !pre_sequence_source.diagnostics.numerical_packet_matrices_built
    @test !pre_sequence_source.diagnostics.operator_data_available
    @test !pre_sequence_source.diagnostics.packet_operator_data_checked
    layer_region = CCP.cartesian_shell_region(
        layer;
        parent_dimension = prod(GaussletBases._nested_axis_lengths(bundles)),
    )
    @test layer_region isa CCP.CartesianShellRegion3D
    @test layer_region.region_family == :endcap_panel_shared_exterior
    @test layer_region.role == :shared_endcap_panel_shell_layer
    @test layer_region.status == :transitional
    @test layer_region.box == layer.provenance.current_box
    @test layer_region.inner_exclusion_box == layer.provenance.inner_box
    @test layer_region.support_summary.entry_count == length(layer.support_indices)
    @test layer_region.support_summary.outside_count == 0
    @test layer_region.ownership_coverage_contract == :boundary_only
    @test layer_region.retention.retention_rule == :old_endcap_panel_product_split
    @test layer_region.retention.cleanup_rule == :locally_orthonormal_product_doside
    @test layer_region.retention.preferred_contraction_rule ==
          :old_endcap_panel_product_split
    @test layer_region.retention.expected_unit_family == :product_owned_unit
    @test layer_region.retention.metric_capability == :product_staged_metric_contraction
    @test isempty(layer_region.retention.missing_payload_fields)
    @test layer_region.current_route_consumes
    @test !layer_region.descriptor_drives_builder
    @test !layer_region.descriptor_only
    @test layer_region.geometry.q == 4
    @test layer_region.geometry.L == 4
    @test layer_region.geometry.unit_count == 6
    @test layer_region.diagnostics.transitional_current_active_implementation

    @test size(default_source.sequence.coefficient_matrix) == (539, 347)
    @test size(endcap_source.sequence.coefficient_matrix) == (539, 313)
    @test all(isfinite, endcap_source.sequence.packet.overlap)
    @test norm(endcap_source.sequence.packet.overlap - I, Inf) < 1.0e-8

    fixed_block = GaussletBases._nested_fixed_block(endcap_source)
    @test fixed_block isa GaussletBases._NestedFixedBlock3D
    @test size(fixed_block.coefficient_matrix) == (539, 313)
    @test fixed_block.staged_by_center_sidecar[] isa
          GaussletBases._CartesianNestedProductStagedByCenterSidecar3D
    @test fixed_block.staged_by_center_sidecar[].diagnostics.product_unit_count == 6
    @test fixed_block.staged_by_center_sidecar[].diagnostics.final_dimension == 313
    @test fixed_block.staged_by_center_sidecar[].diagnostics.max_support_count <= 225
    @test all(isfinite, fixed_block.overlap)
    @test norm(fixed_block.overlap - I, Inf) < 1.0e-8
    @test all(isfinite, fixed_block.weights)

    CCS = GaussletBases.CartesianCarriedSpaces
    contracted_parent = CCP.cartesian_contracted_parent(fixed_block)
    sidecar_units = fixed_block.staged_by_center_sidecar[].units
    parent_dim = CP.parent_dimension(CCP.contracted_parent_basis(contracted_parent))
    rule_built_units = [
        CCP.cartesian_contraction_unit_from_rule(
            CCP.cartesian_contraction_rule(sidecar_unit; parent_dimension = parent_dim),
            sidecar_unit,
        ) for sidecar_unit in sidecar_units
    ]
    @test CCP.contracted_parent_metadata(contracted_parent).staged_by_center_path ==
        :product_staged_factorized
    @test length(CCP.contracted_parent_units(contracted_parent)) == length(sidecar_units)
    @test length(rule_built_units) == length(sidecar_units)
    @test first(CCP.contracted_parent_units(contracted_parent)).metadata.staged_by_center_unit ===
        first(sidecar_units)
    product_sidecar_unit = first(unit for unit in sidecar_units if unit.kind == :product_doside)
    support_sidecar_unit = first(unit for unit in sidecar_units if unit.kind == :support_dense)
    product_resolved_payload = CCP._cartesian_resolved_contraction_payload(
        product_sidecar_unit;
        parent_dimension = parent_dim,
    )
    support_resolved_payload = CCP._cartesian_resolved_contraction_payload(
        support_sidecar_unit;
        parent_dimension = parent_dim,
    )
    @test product_resolved_payload isa CCP._CartesianResolvedContractionPayload3D
    @test support_resolved_payload isa CCP._CartesianResolvedContractionPayload3D
    @test product_resolved_payload.metric_path == :product_staged_metric_contraction
    @test support_resolved_payload.metric_path == :support_local_product
    @test product_resolved_payload.ready_for_metric_execution
    @test support_resolved_payload.ready_for_metric_execution
    @test product_resolved_payload.payload_kind == :product_doside
    @test support_resolved_payload.payload_kind == :support_dense
    @test product_resolved_payload.column_range == product_sidecar_unit.column_range
    @test support_resolved_payload.column_range == support_sidecar_unit.column_range
    @test product_resolved_payload.support_indices == product_sidecar_unit.support_indices
    @test support_resolved_payload.support_indices == support_sidecar_unit.support_indices
    @test product_resolved_payload.support_states == product_sidecar_unit.support_states
    @test support_resolved_payload.support_states == support_sidecar_unit.support_states
    @test product_resolved_payload.payload === product_sidecar_unit
    @test support_resolved_payload.payload === support_sidecar_unit
    @test isempty(product_resolved_payload.missing_fields)
    @test isempty(support_resolved_payload.missing_fields)
    product_sidecar_entries = CCPM._staged_unit_entries(product_sidecar_unit)
    support_sidecar_entries = CCPM._staged_unit_entries(support_sidecar_unit)
    @test length(product_sidecar_entries) == length(product_sidecar_unit.column_range)
    @test length(support_sidecar_entries) == length(support_sidecar_unit.column_range)
    @test sum(length, product_sidecar_entries) ==
          count(!iszero, Matrix{Float64}(product_sidecar_unit.coefficient_matrix))
    @test sum(length, support_sidecar_entries) ==
          count(!iszero, Matrix{Float64}(support_sidecar_unit.coefficient_matrix))
    @test product_resolved_payload.diagnostics.linear_vector_path ==
          :product_staged_axis_projection
    @test support_resolved_payload.diagnostics.linear_vector_path ==
          :support_local_fallback
    @test product_resolved_payload.diagnostics.block_role == :product
    @test support_resolved_payload.diagnostics.block_role == :fallback
    for (sidecar_unit, rule_unit, adapted_unit) in zip(
        sidecar_units,
        rule_built_units,
        CCP.contracted_parent_units(contracted_parent),
    )
        @test rule_unit.role == sidecar_unit.role
        @test rule_unit.support_indices == sidecar_unit.support_indices
        @test rule_unit.column_range == sidecar_unit.column_range
        @test adapted_unit.role == rule_unit.role
        @test adapted_unit.support_indices == rule_unit.support_indices
        @test adapted_unit.column_range == rule_unit.column_range
        @test adapted_unit.metadata.staged_by_center_unit === sidecar_unit
        @test adapted_unit.metadata.contraction_rule.kind == sidecar_unit.kind
        @test adapted_unit.metadata.rule_driven_unit_creation
        @test adapted_unit.metadata.rule_family ==
              adapted_unit.metadata.contraction_rule.rule_family
        @test adapted_unit.metadata.rule_kind ==
              adapted_unit.metadata.contraction_rule.kind
        @test adapted_unit.metadata.rule_metric_capability ==
              adapted_unit.metadata.contraction_rule.metric_capability
        adapted_resolved_payload = CCP._cartesian_resolved_contraction_payload(
            adapted_unit;
            parent_dimension = parent_dim,
        )
        @test adapted_resolved_payload.payload === sidecar_unit
        @test adapted_resolved_payload.column_range == sidecar_unit.column_range
        @test adapted_resolved_payload.ready_for_metric_execution
    end
    contraction_rules = [
        CCP.contraction_unit_rule(
            unit;
            parent_dimension = parent_dim,
        ) for unit in CCP.contracted_parent_units(contracted_parent)
    ]
    product_rules = filter(
        rule -> CCP.contraction_rule_family(rule) == :product_owned_unit,
        contraction_rules,
    )
    support_rules = filter(
        rule -> CCP.contraction_rule_family(rule) == :support_dense_fallback,
        contraction_rules,
    )
    @test length(product_rules) == 6
    @test length(support_rules) >= 1
    first_product_rule = first(product_rules)
    @test first_product_rule.kind == :product_doside
    @test first_product_rule.support_summary.parent_dimension == 539
    @test first_product_rule.support_summary.entry_count ==
          length(first_product_rule.support_indices)
    @test first_product_rule.support_summary.duplicate_count == 0
    @test first_product_rule.support_summary.outside_count == 0
    @test CCP.contraction_rule_retained_dimension(first_product_rule) ==
          length(first_product_rule.column_range)
    @test CCP.contraction_rule_transform_rule(first_product_rule) ==
          :two_active_axis_product_doside
    @test CCP.contraction_rule_cleanup_rule(first_product_rule) ==
          :locally_orthonormal_product_doside
    @test CCP.contraction_rule_metric_capability(first_product_rule) ==
          :product_staged_metric_contraction
    @test first_product_rule.local_geometry.axis_function_index_count ==
          CCP.contraction_rule_retained_dimension(first_product_rule)
    @test first_product_rule.diagnostics.coefficient_contract == :product_doside
    first_support_rule = first(support_rules)
    @test first_support_rule.kind == :support_dense
    @test first_support_rule.support_summary.outside_count == 0
    @test CCP.contraction_rule_transform_rule(first_support_rule) ==
          :explicit_support_dense_coefficients
    @test CCP.contraction_rule_cleanup_rule(first_support_rule) ==
          :external_or_already_cleaned
    @test CCP.contraction_rule_metric_capability(first_support_rule) ==
          :support_local_product
    rule_inventory = CCP.contracted_parent_rule_inventory(contracted_parent)
    rule_family_counts = Dict(rule_inventory.rule_family_counts)
    @test rule_inventory.rule_count == length(CCP.contracted_parent_units(contracted_parent))
    @test rule_inventory.unit_count == length(CCP.contracted_parent_units(contracted_parent))
    @test rule_family_counts[:product_owned_unit] == 6
    @test rule_family_counts[:support_dense_fallback] == length(support_rules)
    @test rule_inventory.parent_dimension == 539
    @test rule_inventory.contracted_dimension == 313
    @test rule_inventory.total_retained_dimension == 313
    @test rule_inventory.support_summary.parent_dimension == 539
    @test rule_inventory.support_summary.outside_count == 0
    @test rule_inventory.support_summary.missing_count == 0
    @test rule_inventory.support_summary.support_complete
    @test Set(rule_inventory.metric_capabilities) ==
          Set([:product_staged_metric_contraction, :support_local_product])
    @test rule_inventory.every_unit_has_rule_metadata
    @test rule_inventory.every_unit_rule_derivable
    @test !rule_inventory.any_metadata_only_rule
    @test !rule_inventory.any_prototype_rule
    @test rule_inventory.diagnostics.parent_level_unit_inventory
    @test rule_inventory.diagnostics.all_rules_have_column_ranges
    dispatch_shadow = CCPM._contracted_parent_metric_dispatch_shadow_plan(contracted_parent)
    resolved_payloads = [
        CCP._cartesian_resolved_contraction_payload(unit; parent_dimension = parent_dim)
        for unit in CCP.contracted_parent_units(contracted_parent)
    ]
    @test dispatch_shadow.comparison.agree
    @test isempty(dispatch_shadow.comparison.mismatch_fields)
    @test dispatch_shadow.payload_plan.plan_supported
    @test dispatch_shadow.rule_plan.plan_supported
    @test dispatch_shadow.payload_plan.unit_count == rule_inventory.rule_count
    @test dispatch_shadow.rule_plan.unit_count == rule_inventory.rule_count
    @test dispatch_shadow.resolved_plan.unit_count == rule_inventory.rule_count
    @test dispatch_shadow.payload_plan.product_unit_count == 6
    @test dispatch_shadow.rule_plan.product_unit_count == 6
    @test dispatch_shadow.resolved_plan.product_unit_count == 6
    @test dispatch_shadow.payload_plan.support_fallback_unit_count == length(support_rules)
    @test dispatch_shadow.rule_plan.support_fallback_unit_count == length(support_rules)
    @test dispatch_shadow.resolved_plan.support_fallback_unit_count == length(support_rules)
    expected_product_blocks = 6 * (6 + 1) ÷ 2
    expected_total_blocks = rule_inventory.rule_count * (rule_inventory.rule_count + 1) ÷ 2
    @test dispatch_shadow.payload_plan.product_product_block_count == expected_product_blocks
    @test dispatch_shadow.rule_plan.product_product_block_count == expected_product_blocks
    @test dispatch_shadow.payload_plan.fallback_block_count ==
          expected_total_blocks - expected_product_blocks
    @test dispatch_shadow.rule_plan.fallback_block_count ==
          expected_total_blocks - expected_product_blocks
    @test dispatch_shadow.resolved_plan.fallback_block_count ==
          expected_total_blocks - expected_product_blocks
    @test dispatch_shadow.payload_plan.unsupported_unit_count == 0
    @test dispatch_shadow.rule_plan.unsupported_unit_count == 0
    @test dispatch_shadow.resolved_plan.unsupported_unit_count == 0
    @test dispatch_shadow.rule_plan.prototype_rule_count == 0
    @test dispatch_shadow.resolved_plan.prototype_rule_count == 0
    @test dispatch_shadow.resolved_comparison.agree
    @test isempty(dispatch_shadow.resolved_comparison.mismatch_fields)
    @test all(payload -> payload.ready_for_metric_execution, resolved_payloads)
    @test [payload.diagnostics.block_role for payload in resolved_payloads] ==
          [path.block_role for path in dispatch_shadow.payload_plan.unit_paths]
    @test [payload.diagnostics.linear_vector_path for payload in resolved_payloads] ==
          [path.linear_vector_path for path in dispatch_shadow.payload_plan.unit_paths]
    @test [payload.diagnostics.metric_capability for payload in resolved_payloads] ==
          [path.metric_capability for path in dispatch_shadow.payload_plan.unit_paths]
    @test [path.path for path in dispatch_shadow.payload_plan.block_paths] ==
          [path.path for path in dispatch_shadow.rule_plan.block_paths]
    packet_build_plan = CCP._cartesian_packet_build_plan(
        fixed_block.staged_by_center_sidecar[],
    )
    packet_build_source = packet_build_plan.source
    packet_payload_counts = Dict(packet_build_source.payload_kind_counts)
    pre_sequence_payload_counts = Dict(pre_sequence_source.payload_kind_counts)
    @test packet_build_plan isa CCP._CartesianPacketBuildPlan3D
    @test packet_build_source isa CCP._CartesianPacketBuildSource3D
    @test packet_build_source.parent_dimension == parent_dim
    @test packet_build_source.contracted_dimension == 313
    @test pre_sequence_source isa CCP._CartesianPacketBuildSource3D
    @test pre_sequence_source.parent_dimension == packet_build_source.parent_dimension
    @test pre_sequence_source.contracted_dimension == packet_build_source.contracted_dimension
    @test pre_sequence_source.payload_kind_counts == packet_build_source.payload_kind_counts
    @test pre_sequence_payload_counts[:product_doside] == 6
    @test pre_sequence_payload_counts[:support_dense] == length(support_rules)
    @test length(packet_build_source.resolved_payloads) == length(resolved_payloads)
    @test length(pre_sequence_source.resolved_payloads) ==
          length(packet_build_source.resolved_payloads)
    @test [payload.payload_kind for payload in pre_sequence_source.resolved_payloads] ==
          [payload.payload_kind for payload in packet_build_source.resolved_payloads]
    @test [payload.column_range for payload in pre_sequence_source.resolved_payloads] ==
          [payload.column_range for payload in packet_build_source.resolved_payloads]
    @test [payload.support_indices for payload in pre_sequence_source.resolved_payloads] ==
          [payload.support_indices for payload in packet_build_source.resolved_payloads]
    @test [payload.payload for payload in packet_build_source.resolved_payloads] ==
          [payload.payload for payload in resolved_payloads]
    @test [payload.payload_kind for payload in packet_build_source.resolved_payloads] ==
          [payload.payload_kind for payload in resolved_payloads]
    @test [payload.column_range for payload in packet_build_source.resolved_payloads] ==
          [payload.column_range for payload in resolved_payloads]
    @test all(payload -> payload.ready_for_metric_execution, packet_build_source.resolved_payloads)
    @test packet_build_source.column_ranges ==
          [unit.column_range for unit in sidecar_units]
    @test packet_build_source.column_coverage.entry_count == 313
    @test packet_build_source.column_coverage.unique_count == 313
    @test packet_build_source.column_coverage.duplicate_count == 0
    @test packet_build_source.column_coverage.missing_count == 0
    @test packet_build_source.column_coverage.outside_count == 0
    @test packet_build_source.support_union_summary.parent_dimension == parent_dim
    @test packet_build_source.support_union_summary.outside_count == 0
    @test pre_sequence_source.column_coverage.entry_count ==
          packet_build_source.column_coverage.entry_count
    @test pre_sequence_source.column_coverage.unique_count ==
          packet_build_source.column_coverage.unique_count
    @test pre_sequence_source.column_coverage.duplicate_count ==
          packet_build_source.column_coverage.duplicate_count
    @test pre_sequence_source.column_coverage.missing_count ==
          packet_build_source.column_coverage.missing_count
    @test pre_sequence_source.column_coverage.outside_count ==
          packet_build_source.column_coverage.outside_count
    @test pre_sequence_source.support_union_summary.entry_count ==
          packet_build_source.support_union_summary.entry_count
    @test pre_sequence_source.support_union_summary.unique_count ==
          packet_build_source.support_union_summary.unique_count
    @test pre_sequence_source.support_union_summary.outside_count ==
          packet_build_source.support_union_summary.outside_count
    @test packet_payload_counts[:product_doside] == 6
    @test packet_payload_counts[:support_dense] == length(support_rules)
    @test packet_build_source.candidate_packet_fields == (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :weights,
        :first_moments,
        :kinetic,
    )
    @test !(:x2_x in packet_build_source.missing_packet_fields)
    @test !(:x2_y in packet_build_source.missing_packet_fields)
    @test !(:x2_z in packet_build_source.missing_packet_fields)
    @test :nuclear_one_body in packet_build_source.missing_packet_fields
    @test :local_coulomb_one_body in packet_build_source.missing_packet_fields
    @test :local_ecp_one_body in packet_build_source.missing_packet_fields
    @test :gaussian_local_terms in packet_build_source.missing_packet_fields
    @test :mwg_interaction in packet_build_source.missing_packet_fields
    @test packet_build_source.axis_operator_requirements.kinetic ==
          (:overlap, :kinetic)
    @test packet_build_source.diagnostics.metadata_only
    @test packet_build_source.diagnostics.descriptor_does_not_drive_builder
    @test packet_build_source.diagnostics.current_builder_authority ==
          :nested_shell_packet
    @test !packet_build_source.diagnostics.numerical_packet_matrices_built
    @test !packet_build_source.diagnostics.operator_data_available
    @test !packet_build_source.diagnostics.packet_operator_data_checked
    @test !packet_build_source.diagnostics.overlap_matrix_built
    @test !packet_build_source.diagnostics.kinetic_matrix_built
    @test packet_build_source.diagnostics.column_ranges_cover_contract
    @test packet_build_source.diagnostics.column_layout_ready_for_candidate_fields
    @test packet_build_source.diagnostics.support_union_summary_informational
    @test !packet_build_source.diagnostics.parent_support_complete_required
    @test packet_build_source.diagnostics.overlapping_payload_support_allowed
    @test pre_sequence_source.candidate_packet_fields ==
          packet_build_source.candidate_packet_fields
    @test pre_sequence_source.missing_packet_fields ==
          packet_build_source.missing_packet_fields
    @test pre_sequence_source.diagnostics.packet_construction_consumes_source == false
    @test pre_sequence_source.diagnostics.source_object_builds_packet_matrices == false
    @test pre_sequence_source.diagnostics.nested_shell_packet_remains_authoritative
    @test !pre_sequence_source.diagnostics.numerical_packet_matrices_built
    @test !pre_sequence_source.diagnostics.operator_data_available
    @test !pre_sequence_source.diagnostics.packet_operator_data_checked
    @test !packet_build_source.diagnostics.full_packet_builder_ready
    @test packet_build_plan.current_builder_authority == :nested_shell_packet
    @test !packet_build_plan.descriptor_drives_builder
    @test !packet_build_plan.numerical_packet_matrices_built
    @test packet_build_plan.diagnostics.metadata_only
    @test packet_build_plan.diagnostics.current_nested_shell_packet_authoritative
    @test !packet_build_plan.diagnostics.fixed_block_construction_changed
    @test !packet_build_plan.diagnostics.metric_packet_execution_changed
    @test !packet_build_plan.diagnostics.qwhamiltonian_changed
    @test !packet_build_plan.diagnostics.backend_policy_changed
    @test !packet_build_plan.diagnostics.quadrature_policy_changed
    @test !packet_build_plan.diagnostics.ida_positive_weight_semantics_changed
    @test !packet_build_plan.diagnostics.cr2_science_status_changed
    @test !packet_build_plan.diagnostics.operator_data_available
    @test !packet_build_plan.diagnostics.packet_operator_data_checked
    support_metric_packet = CCPM.cartesian_contracted_parent_metric_packet(
        contracted_parent;
        axis_metrics,
    )
    product_metric_packet = CCPM.cartesian_contracted_parent_metric_packet(
        contracted_parent;
        axis_metrics,
        construction_path = :product_staged_metric_contraction,
    )
    resolved_metric_packet = CCPM._resolved_payload_product_staged_metric_packet(
        contracted_parent;
        axis_metrics,
    )
    @test CP.parent_dimension(CCP.contracted_parent_basis(contracted_parent)) == 539
    @test support_metric_packet.diagnostics.construction_path == :support_local_product
    @test product_metric_packet.diagnostics.construction_path == :product_staged_metric_contraction
    @test !(:resolved_payload_count in propertynames(product_metric_packet.diagnostics))
    @test !(:default_metric_execution_changed in propertynames(product_metric_packet.diagnostics))
    @test resolved_metric_packet.diagnostics.construction_path ==
          :resolved_payload_product_staged_metric_contraction
    @test resolved_metric_packet.diagnostics.source == :resolved_payload_metric_shadow
    @test !resolved_metric_packet.diagnostics.default_metric_execution_changed
    @test resolved_metric_packet.diagnostics.resolved_payload_count ==
          rule_inventory.rule_count
    @test product_metric_packet.diagnostics.dense_parent_matrix_used == false
    @test resolved_metric_packet.diagnostics.dense_parent_matrix_used == false
    @test product_metric_packet.diagnostics.product_unit_count == 6
    @test resolved_metric_packet.diagnostics.product_unit_count ==
          product_metric_packet.diagnostics.product_unit_count
    @test product_metric_packet.diagnostics.generic_unit_count >= 1
    @test resolved_metric_packet.diagnostics.generic_unit_count ==
          product_metric_packet.diagnostics.generic_unit_count
    @test product_metric_packet.diagnostics.product_block_count > 0
    @test resolved_metric_packet.diagnostics.product_block_count ==
          product_metric_packet.diagnostics.product_block_count
    @test product_metric_packet.diagnostics.fallback_block_count > 0
    @test resolved_metric_packet.diagnostics.fallback_block_count ==
          product_metric_packet.diagnostics.fallback_block_count
    @test resolved_metric_packet.overlap ≈ product_metric_packet.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test resolved_metric_packet.weights ≈ product_metric_packet.weights atol = 1.0e-12 rtol = 1.0e-12
    @test resolved_metric_packet.centers ≈ product_metric_packet.centers atol = 1.0e-12 rtol = 1.0e-12
    @test resolved_metric_packet.first_moments ≈ product_metric_packet.first_moments atol = 1.0e-12 rtol = 1.0e-12
    @test product_metric_packet.overlap ≈ support_metric_packet.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test product_metric_packet.weights ≈ support_metric_packet.weights atol = 1.0e-10 rtol = 1.0e-10
    @test product_metric_packet.centers ≈ support_metric_packet.centers atol = 1.0e-10 rtol = 1.0e-10
    @test product_metric_packet.first_moments ≈ support_metric_packet.first_moments atol = 1.0e-10 rtol = 1.0e-10
    @test product_metric_packet.overlap ≈ fixed_block.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test product_metric_packet.weights ≈ fixed_block.weights atol = 1.0e-10 rtol = 1.0e-10
    @test product_metric_packet.centers ≈ fixed_block.fixed_centers atol = 1.0e-10 rtol = 1.0e-10
    pre_sequence_shadow = CCPM._cartesian_packet_build_source_safe_field_shadow(
        pre_sequence_source,
        safe_axis_data,
    )
    @test pre_sequence_shadow.diagnostics.source_driven_shadow_only
    @test !pre_sequence_shadow.diagnostics.construction_adoption
    @test pre_sequence_shadow.diagnostics.current_builder_authority ==
          :nested_shell_packet
    @test pre_sequence_shadow.diagnostics.product_unit_count == 6
    @test pre_sequence_shadow.diagnostics.support_dense_unit_count == length(support_rules)
    @test pre_sequence_shadow.overlap ≈ fixed_block.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test pre_sequence_shadow.position_x ≈ fixed_block.position_x atol = 1.0e-10 rtol = 1.0e-10
    @test pre_sequence_shadow.position_y ≈ fixed_block.position_y atol = 1.0e-10 rtol = 1.0e-10
    @test pre_sequence_shadow.position_z ≈ fixed_block.position_z atol = 1.0e-10 rtol = 1.0e-10
    @test pre_sequence_shadow.x2_x ≈ fixed_block.x2_x atol = 1.0e-10 rtol = 1.0e-10
    @test pre_sequence_shadow.x2_y ≈ fixed_block.x2_y atol = 1.0e-10 rtol = 1.0e-10
    @test pre_sequence_shadow.x2_z ≈ fixed_block.x2_z atol = 1.0e-10 rtol = 1.0e-10
    @test pre_sequence_shadow.weights ≈ fixed_block.weights atol = 1.0e-10 rtol = 1.0e-10
    @test pre_sequence_shadow.kinetic ≈ fixed_block.kinetic atol = 1.0e-10 rtol = 1.0e-10
    @test pre_sequence_shadow.first_moments ≈ product_metric_packet.first_moments atol = 1.0e-10 rtol = 1.0e-10
    safe_field_shadow = CCPM._cartesian_packet_build_source_safe_field_shadow(
        packet_build_source,
        safe_axis_data,
    )
    @test safe_field_shadow.diagnostics.source ==
          :cartesian_packet_build_source_safe_field_shadow
    @test safe_field_shadow.diagnostics.source_driven_shadow_only
    @test !safe_field_shadow.diagnostics.construction_adoption
    @test safe_field_shadow.diagnostics.current_builder_authority ==
          :nested_shell_packet
    @test safe_field_shadow.diagnostics.nested_shell_packet_remains_authoritative
    @test !safe_field_shadow.diagnostics.descriptor_drives_builder
    @test !safe_field_shadow.diagnostics.numerical_packet_authority_changed
    @test !safe_field_shadow.diagnostics.fixed_block_construction_changed
    @test !safe_field_shadow.diagnostics.metric_packet_execution_changed
    @test !safe_field_shadow.diagnostics.qwhamiltonian_changed
    @test !safe_field_shadow.diagnostics.backend_policy_changed
    @test !safe_field_shadow.diagnostics.quadrature_policy_changed
    @test !safe_field_shadow.diagnostics.ida_positive_weight_semantics_changed
    @test !safe_field_shadow.diagnostics.cr2_science_status_changed
    @test safe_field_shadow.diagnostics.requested_fields ==
          packet_build_source.candidate_packet_fields
    @test safe_field_shadow.diagnostics.missing_packet_fields ==
          packet_build_source.missing_packet_fields
    @test safe_field_shadow.diagnostics.x2_built
    @test !safe_field_shadow.diagnostics.gaussian_terms_built
    @test !safe_field_shadow.diagnostics.nuclear_one_body_built
    @test !safe_field_shadow.diagnostics.local_coulomb_one_body_built
    @test !safe_field_shadow.diagnostics.local_ecp_one_body_built
    @test !safe_field_shadow.diagnostics.pair_mwg_interaction_built
    @test safe_field_shadow.diagnostics.product_unit_count == 6
    @test safe_field_shadow.diagnostics.support_dense_unit_count == length(support_rules)
    @test safe_field_shadow.diagnostics.low_order_product_block_count == expected_product_blocks
    @test safe_field_shadow.diagnostics.kinetic_product_block_count == expected_product_blocks
    @test safe_field_shadow.diagnostics.low_order_fallback_block_count ==
          expected_total_blocks - expected_product_blocks
    @test safe_field_shadow.diagnostics.x2_product_block_count == expected_product_blocks
    @test safe_field_shadow.diagnostics.x2_fallback_block_count ==
          expected_total_blocks - expected_product_blocks
    @test safe_field_shadow.diagnostics.kinetic_fallback_block_count ==
          expected_total_blocks - expected_product_blocks
    @test safe_field_shadow.overlap ≈ fixed_block.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test safe_field_shadow.position_x ≈ fixed_block.position_x atol = 1.0e-10 rtol = 1.0e-10
    @test safe_field_shadow.position_y ≈ fixed_block.position_y atol = 1.0e-10 rtol = 1.0e-10
    @test safe_field_shadow.position_z ≈ fixed_block.position_z atol = 1.0e-10 rtol = 1.0e-10
    @test safe_field_shadow.x2_x ≈ fixed_block.x2_x atol = 1.0e-10 rtol = 1.0e-10
    @test safe_field_shadow.x2_y ≈ fixed_block.x2_y atol = 1.0e-10 rtol = 1.0e-10
    @test safe_field_shadow.x2_z ≈ fixed_block.x2_z atol = 1.0e-10 rtol = 1.0e-10
    @test safe_field_shadow.weights ≈ fixed_block.weights atol = 1.0e-10 rtol = 1.0e-10
    @test safe_field_shadow.kinetic ≈ fixed_block.kinetic atol = 1.0e-10 rtol = 1.0e-10
    @test safe_field_shadow.first_moments ≈ product_metric_packet.first_moments atol = 1.0e-10 rtol = 1.0e-10
    @test_throws ArgumentError CCPM._cartesian_packet_build_source_safe_field_shadow(
        packet_build_source,
        safe_axis_data;
        fields = (:gaussian_sum,),
    )
    kinetic_axis_ops = (
        x = (overlap = pgdg_x.overlap, kinetic = pgdg_x.kinetic),
        y = (overlap = pgdg_y.overlap, kinetic = pgdg_y.kinetic),
        z = (overlap = pgdg_z.overlap, kinetic = pgdg_z.kinetic),
    )
    kinetic_product_units = [
        unit for unit in sidecar_units if unit.kind == :product_doside
    ]
    @test length(kinetic_product_units) ==
          fixed_block.staged_by_center_sidecar[].diagnostics.product_unit_count
    kinetic_shadow_packet = CCPM._product_doside_retained_kinetic_shadow_matrix(
        sidecar_units,
        kinetic_axis_ops;
        final_dimension = size(fixed_block.kinetic, 1),
    )
    @test kinetic_shadow_packet.product_units == kinetic_product_units
    @test size(kinetic_shadow_packet.kinetic) == size(fixed_block.kinetic)
    @test kinetic_shadow_packet.diagnostics.product_only_shadow
    @test !kinetic_shadow_packet.diagnostics.full_kinetic_packet
    @test kinetic_shadow_packet.diagnostics.support_dense_blocks_absent
    @test kinetic_shadow_packet.diagnostics.mixed_blocks_absent
    @test !kinetic_shadow_packet.diagnostics.default_execution_changed
    @test !kinetic_shadow_packet.diagnostics.qwhamiltonian_consumes
    @test !kinetic_shadow_packet.diagnostics.backend_policy_changed
    @test !kinetic_shadow_packet.diagnostics.quadrature_policy_changed
    @test !kinetic_shadow_packet.diagnostics.cr2_science_status_changed
    @test kinetic_shadow_packet.diagnostics.product_unit_count ==
          length(kinetic_product_units)
    @test kinetic_shadow_packet.diagnostics.final_dimension ==
          size(fixed_block.kinetic, 1)
    kinetic_shadow_block_count = 0
    kinetic_shadow_errors = Float64[]
    for right_index in eachindex(kinetic_product_units)
        right = kinetic_product_units[right_index]
        for left_index in 1:right_index
            left = kinetic_product_units[left_index]
            kinetic_shadow_block =
                kinetic_shadow_packet.kinetic[left.column_range, right.column_range]
            kinetic_reference_block =
                fixed_block.kinetic[left.column_range, right.column_range]
            push!(
                kinetic_shadow_errors,
                maximum(abs.(kinetic_shadow_block .- kinetic_reference_block)),
            )
            @test kinetic_shadow_block ≈ kinetic_reference_block atol = 1.0e-10 rtol = 1.0e-10
            if left_index != right_index
                kinetic_mirror_reference =
                    fixed_block.kinetic[right.column_range, left.column_range]
                @test transpose(kinetic_shadow_block) ≈ kinetic_mirror_reference atol = 1.0e-10 rtol = 1.0e-10
            end
            kinetic_shadow_block_count += 1
        end
    end
    @test kinetic_shadow_block_count ==
          length(kinetic_product_units) * (length(kinetic_product_units) + 1) ÷ 2
    @test kinetic_shadow_packet.diagnostics.product_block_count ==
          kinetic_shadow_block_count
    @test maximum(kinetic_shadow_errors) <= 1.0e-10
    product_columns = Int[]
    for unit in kinetic_product_units
        append!(product_columns, unit.column_range)
    end
    product_columns = sort(unique(product_columns))
    nonproduct_columns = setdiff(1:size(fixed_block.kinetic, 1), product_columns)
    @test !isempty(nonproduct_columns)
    @test maximum(abs.(kinetic_shadow_packet.kinetic[nonproduct_columns, :])) == 0.0
    @test maximum(abs.(kinetic_shadow_packet.kinetic[:, nonproduct_columns])) == 0.0
    full_kinetic_shadow_packet = CCPM._staged_retained_kinetic_shadow_matrix(
        sidecar_units,
        kinetic_axis_ops;
        final_dimension = size(fixed_block.kinetic, 1),
    )
    expected_kinetic_total_blocks = length(sidecar_units) * (length(sidecar_units) + 1) ÷ 2
    @test size(full_kinetic_shadow_packet.kinetic) == size(fixed_block.kinetic)
    @test full_kinetic_shadow_packet.diagnostics.full_private_shadow_matrix
    @test !full_kinetic_shadow_packet.diagnostics.product_only_shadow
    @test !full_kinetic_shadow_packet.diagnostics.production_adoption
    @test !full_kinetic_shadow_packet.diagnostics.default_execution_changed
    @test !full_kinetic_shadow_packet.diagnostics.metric_packet_execution_changed
    @test !full_kinetic_shadow_packet.diagnostics.fixed_block_construction_changed
    @test !full_kinetic_shadow_packet.diagnostics.qwhamiltonian_consumes
    @test !full_kinetic_shadow_packet.diagnostics.backend_policy_changed
    @test !full_kinetic_shadow_packet.diagnostics.quadrature_policy_changed
    @test !full_kinetic_shadow_packet.diagnostics.ida_positive_weight_semantics_changed
    @test !full_kinetic_shadow_packet.diagnostics.cr2_science_status_changed
    @test full_kinetic_shadow_packet.diagnostics.product_unit_count ==
          length(kinetic_product_units)
    @test full_kinetic_shadow_packet.diagnostics.generic_unit_count ==
          length(sidecar_units) - length(kinetic_product_units)
    @test full_kinetic_shadow_packet.diagnostics.product_block_count ==
          kinetic_shadow_block_count
    @test full_kinetic_shadow_packet.diagnostics.fallback_block_count ==
          expected_kinetic_total_blocks - kinetic_shadow_block_count
    @test full_kinetic_shadow_packet.diagnostics.total_block_count ==
          expected_kinetic_total_blocks
    @test full_kinetic_shadow_packet.diagnostics.final_dimension ==
          size(fixed_block.kinetic, 1)
    @test maximum(abs.(full_kinetic_shadow_packet.kinetic .- fixed_block.kinetic)) <= 1.0e-10

    carried = CCS.cartesian_carried_space(fixed_block)
    carried_parent = CCS.carried_space_parent(carried)
    carried_contracted_parent = CCS.carried_space_contracted_parent(carried)
    carried_representation = CCS.carried_space_representation(carried)
    carried_diagnostics = CCS.carried_space_diagnostics(carried)
    @test carried_parent isa CP.CartesianParentGaussletBasis3D
    @test carried_contracted_parent isa CCP.CartesianContractedParent3D
    @test carried_representation isa CartesianBasisRepresentation3D
    @test CP.parent_dimension(carried_parent) == 539
    @test carried_diagnostics.parent_axis_counts == CP.parent_axis_counts(carried_parent)
    @test carried_diagnostics.has_contracted_parent
    @test carried_diagnostics.contracted_dimension == 313
    @test carried_diagnostics.contracted_dimension_matches_representation
    @test carried_diagnostics.contracted_parent_dimension_matches_parent
    @test carried_diagnostics.has_staged_sidecar
    @test carried_diagnostics.staged_by_center_path == :product_staged_factorized
    @test carried_diagnostics.dense_parent_matrix_used == false
    @test carried_diagnostics.heavy_metric_packet_built == false
    @test CCS.carried_space_provenance(carried).input_kind == :nested_fixed_block
    @test CCP.contracted_parent_metadata(carried_contracted_parent).staged_by_center_sidecar ===
        fixed_block.staged_by_center_sidecar[]
    @test first(CCP.contracted_parent_units(carried_contracted_parent)).metadata.staged_by_center_unit ===
        first(sidecar_units)
    carried_metric_packet = CCPM.cartesian_contracted_parent_metric_packet(
        carried_contracted_parent;
        axis_metrics,
        construction_path = :product_staged_metric_contraction,
    )
    @test carried_metric_packet.diagnostics.construction_path ==
        :product_staged_metric_contraction
    @test carried_metric_packet.diagnostics.dense_parent_matrix_used == false
    @test carried_metric_packet.overlap ≈ product_metric_packet.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test carried_metric_packet.weights ≈ product_metric_packet.weights atol = 1.0e-10 rtol = 1.0e-10
    @test carried_metric_packet.centers ≈ product_metric_packet.centers atol = 1.0e-10 rtol = 1.0e-10

    @test_throws ArgumentError GaussletBases._nested_bond_aligned_diatomic_source(
        basis,
        bundles;
        bond_axis = :z,
        nside = 5,
        term_coefficients,
        packet_kernel = :factorized_direct,
        shared_shell_layer_policy = :bad_policy,
    )
end
