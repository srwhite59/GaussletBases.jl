include("pqs_source_metadata_real_artifact_acceptance_runtests.jl")
include("pqs_component_route_report_adapter_runtests.jl")
include("pqs_standard_source_box_route_setup_runtests.jl")
include("pqs_standard_parent_axis_readiness_runtests.jl")
include("pqs_explicit_core_spacing_parent_axis_probe_runtests.jl")
include("pqs_route_axis_count_selection_runtests.jl")
include("pqs_raw_product_box_plan_probe_runtests.jl")
include("pqs_source_box_route_skeleton_runtests.jl")
include("pqs_source_box_route_driver_report_runtests.jl")
include("pqs_source_box_route_driver_crc_print_line_runtests.jl")
include("cartesian_route_core_examples_runtests.jl")
include("cartesian_shellification_module_runtests.jl")
include("cartesian_driver_module_boundary_runtests.jl")
include("cartesian_terminal_shellification_geometry_runtests.jl")
include("cartesian_selected_terminal_lowering_contract_inventory_runtests.jl")
include("cartesian_route_core_selected_terminal_lowering_sidecar_runtests.jl")
include("cartesian_pair_stage_fingerprint_helpers_runtests.jl")
include("cartesian_shellification_plan_runtests.jl")
include("cartesian_ham_builder_one_center_config_smoke_runtests.jl")
include("cartesian_ham_builder_diatomic_config_smoke_runtests.jl")
include("cartesian_route_diatomic_materializer_probe_runtests.jl")
include("white_lindsey_materialized_seed_runtests.jl")

@testset "Cartesian nested face first primitive" begin
    function _fixed_a_nested_test_basis(count::Int; a::Float64 = 0.25, xmax::Float64 = 10.0, tail_spacing::Float64 = 10.0)
        endpoint = (count - 1) / 2
        s = asinh(xmax / a) / (endpoint - xmax / tail_spacing)
        basis = build_basis(MappedUniformBasisSpec(:G10;
            count,
            mapping = AsinhMapping(; a, s, tail_spacing),
            reference_spacing = 1.0,
        ))
        return basis, s
    end

    expansion = coulomb_gaussian_expansion(doacc = false)
    basis, s = _fixed_a_nested_test_basis(13)
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    pgdg = bundle.pgdg_intermediate
    interval = 2:(length(basis) - 1)
    side = GaussletBases._nested_doside_1d(bundle, interval, 4)
    pgdg_exact_bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :pgdg_localized_experimental,
        refinement_levels = 0,
    )
    source_axis = GaussletBases._cartesian_source_box_axis_transform(
        pgdg_exact_bundle,
        interval,
        4;
        axis = :x,
        enforce_symmetric_odd = false,
    )
    source_axis_plan = GaussletBases._cartesian_source_box_axis_transform_plan(
        GaussletBases._CartesianNestedAxisBundles3D(
            pgdg_exact_bundle,
            pgdg_exact_bundle,
            pgdg_exact_bundle,
        ),
        (interval, interval, interval),
        (4, 4, 5);
        enforce_symmetric_odd = false,
    )
    pgdg_exact_bundles = GaussletBases._CartesianNestedAxisBundles3D(
        pgdg_exact_bundle,
        pgdg_exact_bundle,
        pgdg_exact_bundle,
    )
    cubic_raw_box_plan = GaussletBases._cartesian_raw_product_box_plan(
        pgdg_exact_bundles,
        (interval, interval, interval),
        (5, 5, 5);
        enforce_symmetric_odd = false,
    )
    rectangular_direct_axis_plan =
        GaussletBases._cartesian_source_box_axis_transform_plan(
            pgdg_exact_bundles,
            (interval, interval, interval),
            (5, 5, 7);
            enforce_symmetric_odd = false,
        )
    rectangular_raw_box_plan = GaussletBases._cartesian_raw_product_box_plan(
        pgdg_exact_bundles,
        (interval, interval, interval),
        (5, 5, 7);
        enforce_symmetric_odd = false,
    )

    @test s > 0.0
    @test side isa GaussletBases._CartesianNestedDoSide1D
    @test side.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test side.interval == interval
    @test side.retained_count == 3
    @test size(side.local_coefficients) == (length(interval), 3)
    @test size(side.coefficient_matrix) == (length(basis), 3)
    @test maximum(abs.(side.coefficient_matrix[1:(first(interval) - 1), :])) == 0.0
    @test maximum(abs.(side.coefficient_matrix[(last(interval) + 1):end, :])) == 0.0
    @test norm(transpose(side.local_coefficients) * side.local_overlap * side.local_coefficients - I, Inf) < 1.0e-10
    @test norm(transpose(side.coefficient_matrix) * pgdg.overlap * side.coefficient_matrix - I, Inf) < 1.0e-10
    @test issorted(side.localized_centers)
    @test length(side.localized_weights) == 3
    @test any(abs.(side.localized_centers) .< 1.0e-10)
    @test source_axis.object_kind == :cartesian_source_box_axis_transform_1d
    @test source_axis.axis == :x
    @test source_axis.interval == interval
    @test source_axis.source_mode_dim_requested == 4
    @test source_axis.source_mode_dim == 4
    @test !source_axis.source_mode_dim_adjusted
    @test source_axis.integration_contract == :pgdg_exact
    @test source_axis.integration_contract_label == "pgdg-exact"
    @test source_axis.diagnostics.pgdg_backend == :pgdg_localized_experimental
    @test source_axis.diagnostics.exact_with_respect_to_pgdg_proxy_basis
    @test !source_axis.diagnostics.numerical_reference_fallback
    @test source_axis.diagnostics.coefficient_overlap_error < 1.0e-10
    @test source_axis_plan.object_kind ==
          :cartesian_source_box_axis_transform_plan_3d
    @test source_axis_plan.source_box == (interval, interval, interval)
    @test source_axis_plan.source_mode_dims == (4, 4, 5)
    @test source_axis_plan.source_mode_count == 80
    @test !source_axis_plan.diagnostics.source_mode_dims_adjusted
    @test source_axis_plan.diagnostics.integration_contract == :pgdg_exact
    @test source_axis_plan.diagnostics.max_axis_overlap_error < 1.0e-10
    @test cubic_raw_box_plan.object_kind == :cartesian_raw_product_box_plan_3d
    @test cubic_raw_box_plan.source_box == (interval, interval, interval)
    @test cubic_raw_box_plan.axis_intervals == (interval, interval, interval)
    @test cubic_raw_box_plan.source_mode_dims == (5, 5, 5)
    @test cubic_raw_box_plan.source_mode_count == 125
    @test length(cubic_raw_box_plan.source_mode_indices) == 125
    @test first(cubic_raw_box_plan.source_mode_indices) == (1, 1, 1)
    @test cubic_raw_box_plan.source_mode_indices[2] == (1, 1, 2)
    @test cubic_raw_box_plan.source_mode_indices[6] == (1, 2, 1)
    @test last(cubic_raw_box_plan.source_mode_indices) == (5, 5, 5)
    @test cubic_raw_box_plan.source_mode_column_indices == collect(1:125)
    @test cubic_raw_box_plan.source_mode_ordering == :x_major_y_major_z_fast
    @test all(
        axis -> size(cubic_raw_box_plan.axis_local_coefficients[axis]) ==
                (length(interval), 5),
        1:3,
    )
    @test cubic_raw_box_plan.axis_transform_plan.source_mode_dims == (5, 5, 5)
    @test cubic_raw_box_plan.diagnostics.source_mode_dims_are_total_lengths
    @test cubic_raw_box_plan.diagnostics.deterministic_given_box_and_dims
    @test cubic_raw_box_plan.diagnostics.integration_contract == :pgdg_exact
    @test !cubic_raw_box_plan.diagnostics.numerical_reference_fallback
    @test !cubic_raw_box_plan.diagnostics.retained_rule_attached
    @test !cubic_raw_box_plan.diagnostics.packet_adoption
    @test cubic_raw_box_plan.diagnostics.max_axis_overlap_error < 1.0e-10
    @test cubic_raw_box_plan.diagnostics.source_product_modes_orthogonal
    @test rectangular_raw_box_plan.source_mode_dims == (5, 5, 7)
    @test rectangular_raw_box_plan.source_mode_count == 175
    @test length(rectangular_raw_box_plan.source_mode_indices) == 175
    @test rectangular_raw_box_plan.source_mode_indices[7] == (1, 1, 7)
    @test rectangular_raw_box_plan.source_mode_indices[8] == (1, 2, 1)
    @test last(rectangular_raw_box_plan.source_mode_indices) == (5, 5, 7)
    @test all(
        axis -> rectangular_raw_box_plan.axis_transform_plan.axes[axis].source_mode_dim ==
                rectangular_direct_axis_plan.axes[axis].source_mode_dim,
        1:3,
    )
    @test all(
        axis -> rectangular_raw_box_plan.axis_local_coefficients[axis] ≈
                rectangular_direct_axis_plan.axes[axis].local_coefficients,
        1:3,
    )
    @test rectangular_raw_box_plan.diagnostics.max_axis_overlap_error < 1.0e-10
    @test rectangular_raw_box_plan.diagnostics.source_product_modes_orthogonal

    face_lo = GaussletBases._nested_xy_face_product(
        pgdg,
        interval,
        interval,
        1;
        retain_x = 4,
        retain_y = 3,
    )
    face_hi = GaussletBases._nested_xy_face_product(
        pgdg,
        interval,
        interval,
        length(basis);
        retain_x = 4,
        retain_y = 3,
    )
    face_overlap = GaussletBases._nested_xy_face_overlap(face_lo, pgdg.overlap)
    face_cross = GaussletBases._nested_xy_face_cross_overlap(face_lo, face_hi, pgdg.overlap)

    @test face_lo isa GaussletBases._CartesianNestedXYFace3D
    @test face_lo.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test size(face_lo.coefficient_matrix) == (length(basis)^3, 9)
    @test length(face_lo.support_indices) == length(interval)^2
    @test isempty(intersect(face_lo.support_indices, face_hi.support_indices))
    @test norm(face_overlap - I, Inf) < 1.0e-10
    @test norm(face_cross, Inf) < 1.0e-10
end

@testset "Cartesian nested owned-unit coverage audit" begin
    dense_unit = GaussletBases._CartesianNestedOwnedUnit3D(
        :endcap_a,
        [1, 2],
        [1.0 0.0; 0.0 1.0];
        metadata = (side = :left,),
    )
    sparse_unit = GaussletBases._CartesianNestedOwnedUnit3D(
        :panel_b,
        [3, 4],
        sparse([1, 2], [1, 1], [0.5, 0.5], 2, 1);
        metadata = (side = :right,),
    )
    exact = GaussletBases._nested_owned_unit_coverage_audit(
        [dense_unit, sparse_unit],
        [1, 2, 3, 4],
    )

    @test dense_unit.coefficient_matrix isa Matrix{Float64}
    @test sparse_unit.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test dense_unit.role == :endcap_a
    @test dense_unit.metadata.side == :left
    @test exact.expected_support_count == 4
    @test exact.owned_support_count == 4
    @test exact.duplicate_count == 0
    @test exact.missing_count == 0
    @test exact.outside_count == 0
    @test exact.retained_count == 3
    @test exact.coverage_ok

    duplicate_unit = GaussletBases._CartesianNestedOwnedUnit3D(:duplicate_panel, [2, 3], ones(2, 1))
    duplicate = GaussletBases._nested_owned_unit_coverage_audit(
        [dense_unit, duplicate_unit],
        [1, 2, 3],
    )
    @test duplicate.duplicate_count == 1
    @test duplicate.missing_count == 0
    @test duplicate.outside_count == 0
    @test !duplicate.coverage_ok

    missing = GaussletBases._nested_owned_unit_coverage_audit([dense_unit], [1, 2, 3])
    @test missing.missing_count == 1
    @test !missing.coverage_ok

    outside = GaussletBases._nested_owned_unit_coverage_audit([dense_unit], [1])
    @test outside.outside_count == 1
    @test !outside.coverage_ok

    zero_retained = GaussletBases._CartesianNestedOwnedUnit3D(:empty_panel, [1, 2], zeros(2, 0))
    nonfinite = GaussletBases._CartesianNestedOwnedUnit3D(:bad_panel, [1], [Inf;;])
    @test_throws DimensionMismatch GaussletBases._CartesianNestedOwnedUnit3D(:bad_rows, [1, 2], ones(1, 1))
    @test_throws ArgumentError GaussletBases._nested_owned_unit_coverage_audit([zero_retained], [1, 2])
    @test_throws ArgumentError GaussletBases._nested_owned_unit_coverage_audit([nonfinite], [1])
    @test_throws ArgumentError GaussletBases._nested_owned_unit_coverage_audit([dense_unit], [1, 1, 2])
end

include("pqs_projected_q_shell_local_layer_integration_runtests.jl")

include("cartesian_endcap_panel_owned_shell_producer_runtests.jl")

@testset "Bond-aligned diatomic atom-growth anatomy policy" begin
    parent29 = (1:29, 1:29, 1:29)
    recipe29 = GaussletBases._nested_bond_aligned_diatomic_atom_growth_recipe(
        parent29;
        bond_axis = :z,
        atom_axis_indices = (10, 20),
        protected_atom_side_count = 5,
    )
    anatomy29 = GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(recipe29)

    @test recipe29 isa GaussletBases._BondAlignedDiatomicAtomGrowthRecipe3D
    @test recipe29.contact_cap_policy == :single_shared_contact_cap
    @test recipe29.mismatch_absorption_policy == :outermost_shared_molecular_shell
    @test anatomy29.atom_side_count_ladder == [5, 7, 9]
    @test anatomy29.final_atom_side_count == 9
    @test anatomy29.contact_gap_count == 1
    @test anatomy29.contact_policy == :single_shared_contact_cap
    @test anatomy29.left_atom_box == (11:19, 11:19, 6:14)
    @test anatomy29.right_atom_box == (11:19, 11:19, 16:24)
    @test anatomy29.contact_box == (11:19, 11:19, 15:15)
    @test anatomy29.inner_atom_contact_box == (11:19, 11:19, 6:24)
    @test anatomy29.outer_regular_start_box == (6:24, 6:24, 1:29)
    @test anatomy29.regular_shared_shell_count == 5
    @test anatomy29.outer_mismatch_low_counts == (5, 5, 0)
    @test anatomy29.outer_mismatch_high_counts == (5, 5, 0)
    @test anatomy29.support_coverage.status == :full_parent_covered
    @test anatomy29.support_coverage.expected_support_count == 29^3
    @test anatomy29.support_coverage.atom_contact_support_count == 9 * 9 * 19
    @test anatomy29.support_coverage.shared_molecular_support_count == 29^3 - 9 * 9 * 19
    @test anatomy29.support_coverage.covered_support_count == 29^3
    @test anatomy29.support_coverage.duplicate_count == 0
    @test anatomy29.support_coverage.missing_count == 0
    @test anatomy29.support_coverage.outside_count == 0
    @test anatomy29.support_coverage.coverage_ok

    even_parent = (1:20, 1:20, 1:20)
    even_anatomy = GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(
        even_parent;
        bond_axis = :z,
        atom_axis_indices = (7, 14),
        protected_atom_side_count = 4,
    )

    @test even_anatomy.atom_side_count_ladder == [4, 6]
    @test even_anatomy.final_atom_side_count == 6
    @test even_anatomy.contact_gap_count == 0
    @test even_anatomy.contact_policy == :touching_atom_boxes
    @test isnothing(even_anatomy.contact_box)
    @test even_anatomy.left_atom_box == (8:13, 8:13, 5:10)
    @test even_anatomy.right_atom_box == (8:13, 8:13, 11:16)
    @test even_anatomy.left_atom_box[1:2] == even_anatomy.right_atom_box[1:2]
    @test 7 - first(even_anatomy.left_atom_box[3]) == 2
    @test last(even_anatomy.left_atom_box[3]) - 7 == 3
    @test 14 - first(even_anatomy.right_atom_box[3]) == 3
    @test last(even_anatomy.right_atom_box[3]) - 14 == 2
    @test even_anatomy.support_coverage.status == :full_parent_covered
    @test even_anatomy.support_coverage.coverage_ok

    asymmetric_parent = (1:11, 1:13, 1:25)
    asymmetric = GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(
        asymmetric_parent;
        bond_axis = :z,
        atom_axis_indices = (8, 17),
        protected_atom_side_count = 5,
    )

    @test asymmetric.contact_policy == :touching_atom_boxes
    @test asymmetric.inner_atom_contact_box == (2:10, 3:11, 4:21)
    @test asymmetric.outer_regular_start_box == (1:11, 2:12, 3:22)
    @test asymmetric.regular_shared_shell_count == 1
    @test asymmetric.outer_mismatch_low_counts == (0, 1, 2)
    @test asymmetric.outer_mismatch_high_counts == (0, 1, 3)
    @test asymmetric.recipe.mismatch_absorption_policy == :outermost_shared_molecular_shell
    @test asymmetric.support_coverage.expected_support_count == 11 * 13 * 25
    @test asymmetric.support_coverage.covered_support_count == 11 * 13 * 25
    @test asymmetric.support_coverage.coverage_ok

    expansion = coulomb_gaussian_expansion(doacc = false)
    real_basis = bond_aligned_homonuclear_qw_basis(
        family = :G10,
        bond_length = 3.0,
        core_spacing = 0.7,
        xmax_parallel = 8.0,
        xmax_transverse = 4.0,
        bond_axis = :z,
    )
    real_bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(real_basis, expansion)
    real_recipe = GaussletBases._nested_bond_aligned_diatomic_atom_growth_recipe(
        real_basis,
        real_bundles;
        protected_atom_side_count = 3,
    )
    real_anatomy = GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(real_recipe)
    @test real_recipe.parent_box == (1:7, 1:7, 1:13)
    @test real_recipe.atom_axis_indices == (5, 9)
    @test real_anatomy.contact_policy == :single_shared_contact_cap
    @test real_anatomy.support_coverage.coverage_ok

    @test_throws ArgumentError GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(
        (1:7, 1:7, 1:13);
        bond_axis = :z,
        atom_axis_indices = (3, 11),
        protected_atom_side_count = 3,
    )
end

@testset "Bond-aligned diatomic atom-growth construction plan" begin
    parent29 = (1:29, 1:29, 1:29)
    anatomy29 = GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(
        parent29;
        bond_axis = :z,
        atom_axis_indices = (10, 20),
        protected_atom_side_count = 5,
    )
    plan29 = GaussletBases._nested_bond_aligned_diatomic_atom_growth_construction_plan(
        anatomy29,
    )

    @test plan29 isa GaussletBases._BondAlignedDiatomicAtomGrowthConstructionPlan3D
    @test plan29.region_order == [
        :outer_mismatch_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :left_atom_box,
        :right_atom_box,
        :contact_cap,
    ]
    @test [region.order_index for region in plan29.regions] == collect(1:9)
    @test plan29.outer_mismatch_is_outermost
    @test plan29.middle_contact_clean
    @test plan29.support_coverage.expected_support_count == 29^3
    @test plan29.support_coverage.region_support_count == 29^3
    @test plan29.support_coverage.covered_support_count == 29^3
    @test plan29.support_coverage.duplicate_count == 0
    @test plan29.support_coverage.missing_count == 0
    @test plan29.support_coverage.outside_count == 0
    @test plan29.support_coverage.coverage_ok
    @test length(plan29.regions[1].support_indices) == 29^3 - 19 * 19 * 29
    @test [
        length(region.support_indices)
        for region in plan29.regions
        if region.role == :regular_shared_molecular_shell
    ] == [2666, 2178, 1738, 1346, 1002]
    @test length(plan29.regions[end - 2].support_indices) == 9^3
    @test length(plan29.regions[end - 1].support_indices) == 9^3
    @test length(plan29.regions[end].support_indices) == 9 * 9
    @test plan29.regions[1].metadata.mismatch_absorption_policy ==
        :outermost_shared_molecular_shell
    @test plan29.regions[2].metadata.shell_offset == 5
    @test plan29.regions[6].metadata.shell_offset == 1
    @test plan29.regions[end].metadata.contact_policy == :single_shared_contact_cap
    CCP = GaussletBases.CartesianContractedParents
    atom_region = CCP.cartesian_shell_region(
        plan29.regions[end - 2];
        parent_dimension = 29^3,
    )
    shared_region = CCP.cartesian_shell_region(
        plan29.regions[2];
        parent_dimension = 29^3,
    )
    contact_region = CCP.cartesian_shell_region(
        plan29.regions[end];
        parent_dimension = 29^3,
    )
    mismatch_region = CCP.cartesian_shell_region(
        plan29.regions[1];
        parent_dimension = 29^3,
    )
    @test atom_region isa CCP.CartesianShellRegion3D
    @test atom_region.region_family == :atom_core_cube
    @test atom_region.status == :clean
    @test atom_region.role == :left_atom_box
    @test atom_region.support_summary.entry_count == 9^3
    @test atom_region.support_summary.outside_count == 0
    @test atom_region.ownership_coverage_contract == :disjoint_partition_piece
    @test atom_region.retention.retention_rule == :protected_atom_cubic_shell
    @test atom_region.retention.preferred_contraction_rule == :complete_shell_sequence
    @test atom_region.retention.metric_capability == :support_local_product
    @test :coefficient_matrix in atom_region.retention.missing_payload_fields
    @test !atom_region.current_route_consumes
    @test !atom_region.descriptor_drives_builder
    @test atom_region.descriptor_only
    @test shared_region.region_family == :rectangular_molecular_shell
    @test shared_region.retention.retention_rule == :policy_selected_shared_exterior
    @test shared_region.retention.metric_capability == :metadata_only_policy_dependent
    @test contact_region.region_family == :shared_midpoint_slab_cap
    @test contact_region.status == :clean
    @test contact_region.retention.preferred_contraction_rule == :contact_cap_owned_slab
    @test !contact_region.current_route_consumes
    @test mismatch_region.region_family == :outer_boundary_shell
    @test mismatch_region.retention.retention_rule ==
          :outermost_mismatch_shared_molecular_shell
    @test mismatch_region.retention.preferred_contraction_rule ==
          :outer_mismatch_boundary_slab_set
    @test mismatch_region.geometry.mismatch_low_counts == (5, 5, 0)
    @test !mismatch_region.current_route_consumes

    even_plan = GaussletBases._nested_bond_aligned_diatomic_atom_growth_construction_plan(
        (1:20, 1:20, 1:20);
        bond_axis = :z,
        atom_axis_indices = (7, 14),
        protected_atom_side_count = 4,
    )
    @test even_plan.support_coverage.coverage_ok
    @test even_plan.middle_contact_clean
    @test !(:contact_cap in even_plan.region_order)
    @test even_plan.region_order[end - 1:end] == [:left_atom_box, :right_atom_box]

    expansion = coulomb_gaussian_expansion(doacc = false)
    real_basis = bond_aligned_homonuclear_qw_basis(
        family = :G10,
        bond_length = 5.0,
        core_spacing = 0.7,
        xmax_parallel = 8.0,
        xmax_transverse = 4.0,
        bond_axis = :z,
    )
    real_bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(real_basis, expansion)
    real_anatomy = GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(
        real_basis,
        real_bundles;
        protected_atom_side_count = 5,
    )
    real_plan = GaussletBases._nested_bond_aligned_diatomic_atom_growth_construction_plan(
        real_anatomy,
    )
    diagnostic_anatomy = bond_aligned_diatomic_nested_geometry_diagnostics(
        real_basis;
        nside = 5,
    ).atom_growth_anatomy.anatomy
    diagnostic_plan = GaussletBases._nested_bond_aligned_diatomic_atom_growth_construction_plan(
        diagnostic_anatomy,
    )
    @test real_plan.region_order == [
        :outer_mismatch_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :left_atom_box,
        :right_atom_box,
        :contact_cap,
    ]
    @test [length(region.support_indices) for region in real_plan.regions] ==
        [98, 362, 125, 125, 25]
    @test real_plan.support_coverage.expected_support_count == 7 * 7 * 15
    @test real_plan.support_coverage.region_support_count == 7 * 7 * 15
    @test real_plan.support_coverage.covered_support_count == 7 * 7 * 15
    @test real_plan.support_coverage.duplicate_count == 0
    @test real_plan.support_coverage.missing_count == 0
    @test real_plan.support_coverage.outside_count == 0
    @test real_plan.support_coverage.coverage_ok
    @test real_plan.outer_mismatch_is_outermost
    @test real_plan.middle_contact_clean
    @test real_plan.regions[1].metadata.low_counts == (0, 0, 1)
    @test real_plan.regions[1].metadata.high_counts == (0, 0, 1)
    @test diagnostic_plan.region_order == real_plan.region_order
    @test diagnostic_plan.support_coverage.coverage_ok
end

@testset "Bond-aligned diatomic high-order recipe policy metadata" begin
    plan29 = GaussletBases._nested_bond_aligned_diatomic_atom_growth_construction_plan(
        (1:29, 1:29, 1:29);
        bond_axis = :z,
        atom_axis_indices = (10, 20),
        protected_atom_side_count = 5,
    )
    policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        plan29;
        q_min = 4,
        atom_q = 4,
        shared_q = 4,
    )
    diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy_diagnostics(
            policy,
        )

    @test policy isa GaussletBases._BondAlignedDiatomicHighOrderRecipePolicy3D
    @test policy.recipe_label == :mixed_atom_cubic_shared_endcap_panel
    @test policy.q_min == 4
    @test policy.q_region_counts == Dict(4 => 9)
    @test policy.buildable_region_count == 7
    @test policy.planned_region_count == 2
    @test policy.experimental_region_count == 0
    @test policy.metadata.active_builder_uses_policy == false
    @test policy.metadata.source_builder_changed == false
    @test diagnostics.region_count == length(plan29.regions)
    @test diagnostics.active_builder_uses_policy == false
    @test diagnostics.support_coverage.coverage_ok
    @test diagnostics.outer_mismatch_is_outermost
    @test diagnostics.middle_contact_clean

    roles = [choice.region_role for choice in policy.region_choices]
    families = [choice.recipe_family for choice in policy.region_choices]
    categories = [choice.region_category for choice in policy.region_choices]
    statuses = [choice.buildability_status for choice in policy.region_choices]
    support_counts = [choice.support_count for choice in policy.region_choices]
    @test roles == plan29.region_order
    @test families[1] == :outermost_mismatch_shared_molecular_shell
    @test all(==(:shared_endcap_panel_exterior), families[2:6])
    @test all(==(:protected_atom_cubic_shell), families[7:8])
    @test families[9] == :shared_contact_cap
    @test categories[1] == :outer_mismatch
    @test all(==(:shared_exterior), categories[2:6])
    @test all(==(:atom_local), categories[7:8])
    @test categories[9] == :contact_cap
    @test statuses[1] == :planned_only
    @test all(==(:buildable_now), statuses[2:8])
    @test statuses[9] == :planned_only
    @test support_counts == [length(region.support_indices) for region in plan29.regions]
    @test all(choice.q == 4 && choice.order == 4 for choice in policy.region_choices)

    annulus_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        plan29;
        q_min = 4,
        atom_q = 4,
        shared_q = 5,
        contact_q = 4,
        outer_mismatch_q = 4,
        shared_exterior_family = :transverse_annulus_exterior,
    )
    annulus_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy_diagnostics(
            annulus_policy,
        )
    @test annulus_policy.recipe_label == :mixed_atom_cubic_shared_transverse_annulus
    @test annulus_policy.q_region_counts == Dict(4 => 4, 5 => 5)
    @test annulus_policy.buildable_region_count == 2
    @test annulus_policy.planned_region_count == 2
    @test annulus_policy.experimental_region_count == 5
    @test all(
        choice.recipe_family == :transverse_annulus_exterior &&
        choice.q == 5 &&
        choice.implementation_status == :promising_experimental &&
        choice.buildability_status == :planned_experimental
        for choice in annulus_policy.region_choices
        if choice.region_category == :shared_exterior
    )
    @test annulus_diagnostics.region_choices[2].recipe_family ==
        :transverse_annulus_exterior
    @test annulus_diagnostics.region_choices[2].q == 5

    @test_throws ArgumentError GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        plan29;
        q_min = 5,
        atom_q = 4,
    )
end

@testset "Bond-aligned diatomic high-order recipe realization audit" begin
    plan29 = GaussletBases._nested_bond_aligned_diatomic_atom_growth_construction_plan(
        (1:29, 1:29, 1:29);
        bond_axis = :z,
        atom_axis_indices = (10, 20),
        protected_atom_side_count = 5,
    )
    policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        plan29;
        q_min = 4,
    )
    audit =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_realization_audit(
            policy,
        )
    diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_realization_diagnostics(
            audit,
        )

    @test audit isa GaussletBases._BondAlignedDiatomicHighOrderRecipeRealizationAudit3D
    @test diagnostics.recipe_label == :mixed_atom_cubic_shared_endcap_panel
    @test diagnostics.region_count == 9
    @test diagnostics.mapped_region_count == 9
    @test diagnostics.missing_region_count == 0
    @test diagnostics.buildable_without_mapped_primitive_count == 0
    @test diagnostics.active_builder_consumed_region_count == 0
    @test diagnostics.descriptor_region_count == 2
    @test diagnostics.exact_descriptor_region_count == 2
    @test diagnostics.ready_for_opt_in_builder == true
    @test diagnostics.active_builder_uses_policy == false
    @test diagnostics.support_coverage.coverage_ok

    realizations = diagnostics.region_realizations
    @test realizations[1].role == :outer_mismatch_shared_molecular_shell
    @test realizations[1].mapped_primitive ==
        :_nested_diatomic_high_order_outer_mismatch_descriptor
    @test realizations[1].mapped_primitive_status == :construction_piece_descriptor
    @test isnothing(realizations[1].missing_implementation)
    @test realizations[1].realization_descriptor.parent_support_count ==
        length(plan29.regions[1].support_indices)
    @test realizations[1].realization_descriptor.owned_unit_count == 4
    @test realizations[1].realization_descriptor.primitive_family ==
        :outer_mismatch_boundary_slab_set
    @test realizations[1].realization_descriptor.support_coverage.coverage_ok
    @test realizations[1].realization_descriptor.exact_full_coverage
    @test all(
        piece.primitive_family == :outer_mismatch_boundary_slab
        for piece in realizations[1].realization_descriptor.pieces
    )
    @test all(
        realization.mapped_primitive == :_nested_endcap_panel_shell_layer &&
        realization.mapped_primitive_status == :existing_internal_primitive &&
        isnothing(realization.missing_implementation) &&
        realization.existing_opt_in_route == :shared_shell_layer_policy_endcap_panel_owned
        for realization in realizations
        if realization.region_category == :shared_exterior
    )
    shared_realizations = [
        realization for realization in realizations
        if realization.region_category == :shared_exterior
    ]
    @test all(
        !isnothing(realization.metadata.realization_notes.projected_q_shell_candidate)
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.mapped_primitive ==
        :_nested_projected_q_shell_layer
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.primitive_family ==
        :projected_q_shell
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.realization_primitive ==
        :projected_q_shell_boundary_comx_product_modes
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.support_contract ==
        :projected_q_shell_raw_boundary
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.coefficient_contract ==
        :full_block_boundary_comx_product_mode_projection
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.cleanup_contract ==
        :full_rank_symmetric_lowdin
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.mode_selection_rule ==
        :any_axis_mode_index_first_or_last
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.sidecar_status ==
        :not_yet_optimized_product_staged_for_pqs
        for realization in shared_realizations
    )
    @test all(
        !realization.metadata.realization_notes.projected_q_shell_candidate.active_builder_consumes &&
        !realization.metadata.realization_notes.projected_q_shell_candidate.source_builder_consumes &&
        !realization.metadata.realization_notes.projected_q_shell_candidate.fixed_block_consumes &&
        !realization.metadata.realization_notes.projected_q_shell_candidate.qw_consumes &&
        !realization.metadata.realization_notes.projected_q_shell_candidate.hamiltonian_consumes
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.current_transitional_implementation ==
        :_nested_endcap_panel_shell_layer
        for realization in shared_realizations
    )
    @test all(
        realization.mapped_primitive == :_nested_bond_aligned_diatomic_sequence_for_box &&
        realization.mapped_primitive_status == :existing_internal_primitive &&
        isnothing(realization.missing_implementation)
        for realization in realizations
        if realization.region_category == :atom_local
    )
    @test realizations[end].role == :contact_cap
    @test realizations[end].mapped_primitive ==
        :_nested_diatomic_high_order_contact_cap_descriptor
    @test realizations[end].mapped_primitive_status == :construction_piece_descriptor
    @test isnothing(realizations[end].missing_implementation)
    @test realizations[end].realization_descriptor.parent_support_count ==
        length(plan29.regions[end].support_indices)
    @test realizations[end].realization_descriptor.owned_unit_count == 1
    @test realizations[end].realization_descriptor.primitive_family == :contact_cap_owned_slab
    @test realizations[end].realization_descriptor.support_coverage.coverage_ok
    @test realizations[end].realization_descriptor.exact_full_coverage
    @test all(
        !(
            realization.buildability_status == :buildable_now &&
            isnothing(realization.mapped_primitive)
        )
        for realization in realizations
    )
    @test all(!realization.active_builder_consumes for realization in realizations)

    annulus_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        plan29;
        q_min = 4,
        atom_q = 4,
        shared_q = 5,
        contact_q = 4,
        outer_mismatch_q = 4,
        shared_exterior_family = :transverse_annulus_exterior,
    )
    annulus_audit =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_realization_audit(
            annulus_policy,
        )
    annulus_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_realization_diagnostics(
            annulus_audit,
        )
    @test annulus_diagnostics.recipe_label == :mixed_atom_cubic_shared_transverse_annulus
    @test annulus_diagnostics.mapped_region_count == 4
    @test annulus_diagnostics.missing_region_count == 5
    @test annulus_diagnostics.buildable_without_mapped_primitive_count == 0
    @test annulus_diagnostics.descriptor_region_count == 2
    @test annulus_diagnostics.exact_descriptor_region_count == 2
    @test annulus_diagnostics.ready_for_opt_in_builder == false
    @test all(
        realization.recipe_family == :transverse_annulus_exterior &&
        realization.q == 5 &&
        realization.mapped_primitive_status == :missing_experimental_primitive &&
        realization.missing_implementation == :transverse_annulus_owned_unit_producer &&
        !(:projected_q_shell_candidate in propertynames(realization.metadata.realization_notes))
        for realization in annulus_diagnostics.region_realizations
        if realization.region_category == :shared_exterior
    )
end

@testset "Bond-aligned diatomic high-order recipe opt-in source construction" begin
    basis = bond_aligned_homonuclear_qw_basis(
        family = :G10,
        bond_length = 5.0,
        core_spacing = 0.7,
        xmax_parallel = 8.0,
        xmax_transverse = 4.0,
        bond_axis = :z,
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(basis, expansion)
    policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        basis,
        bundles;
        protected_atom_side_count = 5,
        q_min = 4,
    )
    realization_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_realization_diagnostics(
            GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_realization_audit(
                policy,
            ),
        )
    @test realization_diagnostics.ready_for_opt_in_builder
    @test !realization_diagnostics.active_builder_uses_policy
    @test realization_diagnostics.active_builder_consumed_region_count == 0

    construction =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :factorized_direct,
            build_sequence_packet = false,
        )
    diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction_diagnostics(
            construction,
        )

    @test construction isa
          GaussletBases._BondAlignedDiatomicHighOrderRecipeSourceConstruction3D
    @test diagnostics.recipe_label == :mixed_atom_cubic_shared_endcap_panel
    @test diagnostics.active_builder_consumes
    @test diagnostics.active_builder_uses_policy
    @test diagnostics.metadata.default_source_builder_changed == false
    @test diagnostics.consumed_region_count == diagnostics.region_count == 5
    @test diagnostics.unsupported_region_count == 0
    @test diagnostics.parent_dimension == 7 * 7 * 15
    @test diagnostics.fixed_dimension == 469
    @test diagnostics.support_coverage.expected_support_count == 7 * 7 * 15
    @test diagnostics.support_coverage.coverage_ok
    @test isnothing(construction.sequence.packet)

    region_builds = diagnostics.region_builds
    @test [build.role for build in region_builds] == [
        :outer_mismatch_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :left_atom_box,
        :right_atom_box,
        :contact_cap,
    ]
    @test [build.primitive_family for build in region_builds] == [
        :outer_mismatch_boundary_slab_set,
        :shared_endcap_panel_shell_layer,
        :atom_local_complete_shell_sequence,
        :atom_local_complete_shell_sequence,
        :contact_cap_owned_slab,
    ]
    @test [build.built_support_count for build in region_builds] == [98, 362, 125, 125, 25]
    @test [build.retained_count for build in region_builds] == [98, 96, 125, 125, 25]
    @test [build.column_range for build in region_builds] ==
          [1:98, 374:469, 99:223, 224:348, 349:373]
    @test all(build.built && build.active_builder_consumes for build in region_builds)
    @test all(build.support_coverage.coverage_ok for build in region_builds)
    @test region_builds[2].metadata.support_contract == :thin_endcap_box_perimeter
    @test region_builds[2].metadata.coefficient_contract == :product_doside
    @test region_builds[5].metadata.descriptor_scope == :middle_contact_cap
    CCP = GaussletBases.CartesianContractedParents
    inventory = CCP.cartesian_shell_region_inventory(
        construction;
        parent_dimension = diagnostics.parent_dimension,
    )
    @test inventory isa CCP.CartesianShellRegionInventory3D
    @test inventory.region_count == diagnostics.region_count == length(region_builds)
    @test inventory.region_order == [build.role for build in region_builds]
    @test [region.role for region in inventory.regions] == inventory.region_order
    @test [region.status for region in inventory.regions] ==
          [:clean, :transitional, :clean, :clean, :clean]
    @test [region.retention.preferred_contraction_rule for region in inventory.regions] == [
        :outer_mismatch_boundary_slab_set,
        :old_endcap_panel_product_split,
        :complete_shell_sequence,
        :complete_shell_sequence,
        :contact_cap_owned_slab,
    ]
    @test inventory.status_counts.clean == 4
    @test inventory.status_counts.transitional == 1
    @test inventory.current_route_consumes_count == diagnostics.region_count
    @test inventory.descriptor_only_count == 0
    @test inventory.support_summary.region_support_entry_count ==
          diagnostics.support_coverage.covered_support_count
    @test inventory.support_summary.support_complete_by_region_counts
    @test inventory.support_summary.count_only_summaries_for_all_regions
    @test all(isnothing(region.support_indices) for region in inventory.regions)
    @test all(region.current_route_consumes for region in inventory.regions)
    @test all(!region.descriptor_drives_builder for region in inventory.regions)
    @test all(!region.descriptor_only for region in inventory.regions)
    @test inventory.regions[2].ownership_coverage_contract == :boundary_only
    @test inventory.regions[2].retention.metric_capability ==
          :product_staged_metric_contraction
    @test isempty(inventory.regions[2].retention.missing_payload_fields)
    @test inventory.regions[5].ownership_coverage_contract == :disjoint_partition_piece
    @test isempty(inventory.regions[5].retention.missing_payload_fields)
    @test diagnostics.metadata.q_policy == :atom_growth_endcap_panel_shared_q_variable
    @test diagnostics.metadata.q4_acceptance_fixture
    @test diagnostics.metadata.non_shared_q_policy == :fixed_q4_order4
    @test diagnostics.metadata.shared_q_values == (4,)
    @test diagnostics.metadata.shared_order_values == (4,)
    @test diagnostics.metadata.shared_shell_realization == :endcap_panel_owned
    @test !diagnostics.metadata.projected_q_shell_opt_in

    be2_basis = bond_aligned_homonuclear_qw_basis(
        family = :G10,
        bond_length = 5.0,
        core_spacing = 0.15,
        xmax_parallel = 10.5,
        xmax_transverse = 8.0,
        bond_axis = :z,
        nuclear_charge = 4.0,
    )
    be2_bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(be2_basis, expansion)
    @test GaussletBases._nested_axis_lengths(be2_bundles) == (15, 15, 27)
    be2_policy0 = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        be2_basis,
        be2_bundles;
        protected_atom_side_count = 5,
        q_min = 4,
    )
    be2_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        be2_policy0.construction_plan;
        q_min = 4,
        atom_q = 4,
        atom_order = 4,
        shared_q = 5,
        shared_order = 5,
        contact_q = 4,
        contact_order = 4,
        outer_mismatch_q = 4,
        outer_mismatch_order = 4,
    )
    be2_retention = GaussletBases._nested_resolve_complete_shell_retention(5)
    be2_shared_dimensions = [
        GaussletBases._nested_diatomic_projected_q_shell_adaptive_source_dimensions(
            be2_basis,
            be2_bundles,
            region,
            be2_retention;
            bond_axis = :z,
            nside = 5,
            selected_q = choice.q,
            shared_shell_angular_resolution_scale = 1.4,
        )
        for (choice, region) in zip(
            be2_policy.region_choices,
            be2_policy.construction_plan.regions,
        )
        if choice.recipe_family == :shared_endcap_panel_exterior
    ]
    be2_source_box_dimension_plans = [
        GaussletBases._nested_diatomic_source_box_dimension_plan(
            be2_basis,
            be2_bundles,
            region.box,
            something(region.inner_exclusion_box),
            be2_retention;
            bond_axis = :z,
            nside = 5,
            selected_q = choice.q,
            shared_shell_angular_resolution_scale = 1.4,
            support_count = length(region.support_indices),
        )
        for (choice, region) in zip(
            be2_policy.region_choices,
            be2_policy.construction_plan.regions,
        )
        if choice.recipe_family == :shared_endcap_panel_exterior
    ]
    be2_shared_regions = [
        region
        for (choice, region) in zip(
            be2_policy.region_choices,
            be2_policy.construction_plan.regions,
        )
        if choice.recipe_family == :shared_endcap_panel_exterior
    ]
    pqs_plan = GaussletBases._nested_projected_q_shell_source_mode_plan(
        (5, 5, 6);
        bond_axis = :z,
        selected_q = 5,
        physical_box_lengths = (11, 11, 21),
        support_count = 1002,
    )
    @test pqs_plan.raw_source_dims == (5, 5, 6)
    @test pqs_plan.source_mode_dims == (5, 5, 6)
    @test pqs_plan.axis_selector_retained_counts == (3, 3, 4)
    @test pqs_plan.raw_q == 5
    @test pqs_plan.raw_L == 6
    @test pqs_plan.raw_q_matches_selected_q
    @test pqs_plan.physical_box_lengths == (11, 11, 21)
    @test pqs_plan.support_count == 1002
    @test pqs_plan.pqs_retained_count == 114
    @test pqs_plan.decomposition_status == :adaptive_broad_support_q_local_modes
    @test !pqs_plan.broad_parent_boundary_reference
    @test !pqs_plan.excluded_from_mvp_gate
    pqs_mismatch_plan = GaussletBases._nested_projected_q_shell_source_mode_plan(
        (4, 4, 5);
        bond_axis = :z,
        selected_q = 5,
        physical_box_lengths = (9, 9, 9),
        support_count = 488,
    )
    @test !pqs_mismatch_plan.raw_q_matches_selected_q
    @test pqs_mismatch_plan.decomposition_status == :adaptive_raw_q_mismatch
    @test pqs_mismatch_plan.excluded_from_mvp_gate
    @test !(:adaptive_retain in propertynames(pqs_plan))
    @test [length.(region.box) for region in be2_shared_regions] ==
          [(15, 15, 25), (13, 13, 23), (11, 11, 21)]
    @test [dims.raw_source_dims for dims in be2_shared_dimensions] ==
          [(5, 5, 5), (5, 5, 5), (5, 5, 6)]
    @test [plan.dimension_policy for plan in be2_source_box_dimension_plans] ==
          fill(:diatomic_adaptive_angular_source_box, 3)
    @test [plan.source_mode_dims for plan in be2_source_box_dimension_plans] ==
          [(5, 5, 5), (5, 5, 5), (5, 5, 6)]
    @test [plan.side_dimensions for plan in be2_source_box_dimension_plans] ==
          [(5, 5, 5), (5, 5, 5), (5, 5, 6)]
    @test [plan.raw_L for plan in be2_source_box_dimension_plans] == [5, 5, 6]
    @test all(plan -> plan.total_source_dimensions_primary, be2_source_box_dimension_plans)
    @test all(
        plan -> plan.diagnostics.source_mode_dims_are_total_lengths,
        be2_source_box_dimension_plans,
    )
    @test all(
        plan -> plan.diagnostics.axis_selector_retained_counts_are_diagnostic,
        be2_source_box_dimension_plans,
    )
    @test [dims.pqs_retained_count for dims in be2_shared_dimensions] ==
          [98, 98, 114]
    @test [
        dims.source_box_dimension_plan.source_mode_dims
        for dims in be2_shared_dimensions
    ] == [(5, 5, 5), (5, 5, 5), (5, 5, 6)]
    @test all(
        dims -> !(:adaptive_retain in propertynames(dims)),
        be2_shared_dimensions,
    )
    @test all(dims -> dims.raw_q == 5, be2_shared_dimensions)
    @test all(dims -> dims.raw_q_matches_selected_q, be2_shared_dimensions)
    @test all(
        dims -> dims.decomposition_status == :adaptive_broad_support_q_local_modes,
        be2_shared_dimensions,
    )
    @test all(dims -> !dims.broad_parent_boundary_reference, be2_shared_dimensions)
    @test all(dims -> !dims.excluded_from_mvp_gate, be2_shared_dimensions)
    @test [
        GaussletBases._nested_diatomic_projected_q_shell_retained_count(length.(region.box))
        for region in be2_shared_regions
    ] == [1738, 1346, 1002]
    @test [
        dims.pqs_retained_count !=
        GaussletBases._nested_diatomic_projected_q_shell_retained_count(length.(region.box))
        for (dims, region) in zip(be2_shared_dimensions, be2_shared_regions)
    ] == [true, true, true]

    pqs_construction =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :support_reference,
            shared_shell_realization = :projected_q_shell,
        )
    pqs_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction_diagnostics(
            pqs_construction,
        )
    @test pqs_diagnostics.metadata.shared_shell_realization == :projected_q_shell
    @test pqs_diagnostics.metadata.projected_q_shell_opt_in
    @test pqs_diagnostics.metadata.default_source_builder_changed == false
    @test pqs_diagnostics.metadata.packet_kernel == :support_reference
    @test pqs_diagnostics.metadata.build_sequence_packet
    @test pqs_diagnostics.support_coverage.coverage_ok
    @test pqs_diagnostics.support_coverage.expected_support_count == 7 * 7 * 15
    @test pqs_diagnostics.support_coverage.covered_support_count == 7 * 7 * 15
    @test pqs_diagnostics.support_coverage.duplicate_count == 0
    @test pqs_diagnostics.support_coverage.missing_count == 0
    @test pqs_diagnostics.support_coverage.outside_count == 0
    @test [build.role for build in pqs_diagnostics.region_builds] ==
          [build.role for build in region_builds]
    @test [build.primitive_family for build in pqs_diagnostics.region_builds] == [
        :outer_mismatch_boundary_slab_set,
        :projected_q_shell,
        :atom_local_complete_shell_sequence,
        :atom_local_complete_shell_sequence,
        :contact_cap_owned_slab,
    ]
    @test [
        build.built_support_count for build in pqs_diagnostics.region_builds
    ] == [build.built_support_count for build in region_builds]
    @test pqs_diagnostics.region_builds[2].mapped_primitive ==
          :_nested_projected_q_shell_layer
    @test pqs_diagnostics.region_builds[2].metadata.support_contract ==
          :projected_q_shell_raw_boundary
    @test pqs_diagnostics.region_builds[2].metadata.coefficient_contract ==
          :full_block_boundary_comx_product_mode_projection
    @test pqs_diagnostics.region_builds[2].metadata.seed_contract ==
          :raw_boundary_projection_of_boundary_comx_product_modes_from_full_local_block_transform
    @test pqs_diagnostics.region_builds[2].metadata.cleanup_contract ==
          :full_rank_symmetric_lowdin
    @test pqs_diagnostics.region_builds[2].metadata.cleanup_method ==
          :projected_boundary_symmetric_lowdin
    @test !pqs_diagnostics.region_builds[2].metadata.pqs_product_staged_sidecar_available
    @test !pqs_diagnostics.region_builds[2].metadata.factorized_direct_allowed
    @test !pqs_diagnostics.region_builds[2].metadata.active_default_builder_changed
    @test pqs_diagnostics.region_builds[2].metadata.selected_q == 4
    @test pqs_diagnostics.region_builds[2].metadata.raw_source_dims == (5, 5, 6)
    @test pqs_diagnostics.region_builds[2].metadata.source_mode_dims == (5, 5, 6)
    @test !(
        :adaptive_retain in
        propertynames(pqs_diagnostics.region_builds[2].metadata)
    )
    @test pqs_diagnostics.region_builds[2].metadata.raw_q == 5
    @test pqs_diagnostics.region_builds[2].metadata.raw_L == 6
    @test !pqs_diagnostics.region_builds[2].metadata.raw_q_matches_selected_q
    @test pqs_diagnostics.region_builds[2].metadata.physical_box_lengths == (7, 7, 13)
    @test pqs_diagnostics.region_builds[2].metadata.support_count == 362
    @test pqs_diagnostics.region_builds[2].metadata.pqs_retained_count == 114
    @test pqs_diagnostics.region_builds[2].metadata.decomposition_status ==
          :adaptive_raw_q_mismatch
    @test !pqs_diagnostics.region_builds[2].metadata.broad_parent_boundary_reference
    @test pqs_diagnostics.region_builds[2].metadata.excluded_from_mvp_gate
    @test pqs_diagnostics.region_builds[2].metadata.policy_q == 4
    @test pqs_diagnostics.region_builds[2].metadata.policy_order == 4
    pqs_source_descriptor =
        pqs_diagnostics.region_builds[2].metadata.pqs_staged_unit_descriptor
    @test pqs_source_descriptor.kind == :projected_q_shell
    @test length.(pqs_source_descriptor.current_box) == (7, 7, 13)
    @test all(
        axis ->
            first(pqs_source_descriptor.inner_box[axis]) ==
            first(pqs_source_descriptor.current_box[axis]) + 1 &&
            last(pqs_source_descriptor.inner_box[axis]) ==
            last(pqs_source_descriptor.current_box[axis]) - 1,
        1:3,
    )
    @test pqs_source_descriptor.bond_axis == :z
    @test pqs_source_descriptor.q == 5
    @test pqs_source_descriptor.L == 6
    @test pqs_source_descriptor.support_count == 362
    @test pqs_source_descriptor.mode_count == 114
    @test pqs_source_descriptor.retained_count == 114
    @test pqs_source_descriptor.cleanup_method == :projected_boundary_symmetric_lowdin
    @test pqs_source_descriptor.cleanup_matrix_size == (114, 114)
    @test pqs_source_descriptor.cleanup_rank_count == 114
    @test pqs_source_descriptor.cleanup_rank_drop_count == 0
    @test pqs_source_descriptor.selection_rule == :any_axis_mode_index_first_or_last
    @test all(
        mode -> any(axis -> mode[axis] == 1 || mode[axis] == (5, 5, 6)[axis], 1:3),
        pqs_source_descriptor.boundary_mode_indices,
    )
    pqs_route_fact_diagnostic =
        CCPM._pqs_pqs_product_route_descriptor_diagnostic(
            pqs_construction,
            _pqs_axis_metrics(bundles),
        )
    @test pqs_route_fact_diagnostic.status == :descriptor_unavailable
    @test pqs_route_fact_diagnostic.descriptor === nothing
    @test :second_pqs_raw_plan in pqs_route_fact_diagnostic.missing
    @test :middle_product_doside_unit in pqs_route_fact_diagnostic.missing
    @test pqs_route_fact_diagnostic.diagnostics.pqs_descriptor_count == 1
    @test pqs_route_fact_diagnostic.diagnostics.pqs_raw_plan_convertible_count == 1
    @test pqs_route_fact_diagnostic.diagnostics.product_doside_unit_count == 0
    @test pqs_route_fact_diagnostic.diagnostics.direct_or_support_body_piece_count == 4
    @test !pqs_route_fact_diagnostic.diagnostics.descriptor_emitted
    @test pqs_route_fact_diagnostic.diagnostics.raw_plan_convertibility_checked
    @test pqs_route_fact_diagnostic.diagnostics.raw_plan_conversion_failure_count == 0
    @test pqs_route_fact_diagnostic.diagnostics.current_route_contains_pqs_descriptors
    @test !pqs_route_fact_diagnostic.diagnostics.current_route_contains_explicit_product_doside_body_unit
    @test !pqs_route_fact_diagnostic.diagnostics.packet_adoption
    @test !pqs_route_fact_diagnostic.diagnostics.fixed_block_construction_changed
    @test !pqs_route_fact_diagnostic.diagnostics.qwhamiltonian_changed
    @test !pqs_route_fact_diagnostic.diagnostics.shell_projection_used
    @test !pqs_route_fact_diagnostic.diagnostics.lowdin_cleanup_used
    @test !pqs_route_fact_diagnostic.diagnostics.support_local_pqs_oracle_used
    @test pqs_route_fact_diagnostic.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !pqs_route_fact_diagnostic.diagnostics.ida_weight_division_allowed
    @test !pqs_route_fact_diagnostic.diagnostics.direct_support_reinterpreted_as_product_doside
    @test any(
        mismatch ->
            hasproperty(mismatch, :reason) &&
            mismatch.reason == :direct_or_support_piece_not_product_doside &&
            hasproperty(mismatch, :primitive_family) &&
            mismatch.primitive_family == :contact_cap_owned_slab,
        pqs_route_fact_diagnostic.mismatches,
    )
    @test any(
        mismatch ->
            hasproperty(mismatch, :reason) &&
            mismatch.reason == :shared_pqs_descriptors_are_not_route_left_right_group,
        pqs_route_fact_diagnostic.mismatches,
    )
    pqs_retained_unit_audit =
        CCPM._pqs_route_retained_unit_fact_audit(pqs_construction)
    @test pqs_retained_unit_audit.object_kind ==
          :pqs_route_retained_unit_fact_audit
    @test pqs_retained_unit_audit.status == :audit_only
    @test length(pqs_retained_unit_audit.unit_facts) == 5
    @test pqs_retained_unit_audit.summary.classification_counts.product_box_constructible == 2
    @test pqs_retained_unit_audit.summary.classification_counts.needs_direct_support_retained_unit_kind == 2
    @test pqs_retained_unit_audit.summary.classification_counts.out_of_scope == 1
    @test pqs_retained_unit_audit.summary.product_box_constructible_slab_rule_count == 3
    pqs_retained_facts_by_role =
        Dict(fact.role => fact for fact in pqs_retained_unit_audit.unit_facts)
    contact_fact = pqs_retained_facts_by_role[:contact_cap]
    @test contact_fact.classification == :product_box_constructible
    @test contact_fact.primitive_family == :contact_cap_owned_slab
    @test !contact_fact.raw_product_box_operator_contract
    @test contact_fact.product_box_construction_rule_available
    @test !contact_fact.product_doside_unit
    @test contact_fact.safe_term_capability ==
          :product_box_rule_available_not_instantiated
    @test contact_fact.coefficient_scope == :support_local_direct_rows
    @test contact_fact.construction_rule.rule_kind ==
          :identity_selector_product_doside_slab
    @test contact_fact.construction_rule.fixed_axis == :z
    @test contact_fact.construction_rule.fixed_index == 8
    @test contact_fact.construction_rule.active_axes == (:x, :y)
    @test contact_fact.construction_rule.active_intervals == (2:6, 2:6)
    @test contact_fact.construction_rule.retained_count == 25
    outer_fact =
        pqs_retained_facts_by_role[:outer_mismatch_shared_molecular_shell]
    @test outer_fact.classification == :product_box_constructible
    @test outer_fact.primitive_family == :outer_mismatch_boundary_slab_set
    @test !outer_fact.raw_product_box_operator_contract
    @test outer_fact.product_box_construction_rule_available
    @test !outer_fact.product_doside_unit
    @test outer_fact.construction_rule.rule_kind ==
          :identity_selector_product_doside_boundary_slab_set
    @test outer_fact.construction_rule.boundary_slab_set
    @test outer_fact.construction_rule.slab_piece_count == 2
    @test all(
        piece_rule -> piece_rule.fixed_axis == :z,
        outer_fact.construction_rule.slab_piece_rules,
    )
    @test sum(
        piece_rule -> piece_rule.support_count,
        outer_fact.construction_rule.slab_piece_rules,
    ) == outer_fact.support_count
    outer_mismatch_product_units =
        CCPM._pqs_outer_mismatch_product_doside_units(pqs_construction)
    @test outer_mismatch_product_units.object_kind ==
          :pqs_outer_mismatch_product_doside_units_fixture
    @test outer_mismatch_product_units.status == :private_diagnostic_only
    @test outer_mismatch_product_units.fact.role ==
          :outer_mismatch_shared_molecular_shell
    @test outer_mismatch_product_units.fact.primitive_family ==
          :outer_mismatch_boundary_slab_set
    @test !outer_mismatch_product_units.fact.raw_product_box_operator_contract
    @test outer_mismatch_product_units.fact.product_box_construction_rule_available
    @test length(outer_mismatch_product_units.units) == 2
    @test map(unit -> unit.kind, outer_mismatch_product_units.units) ==
          (:product_doside, :product_doside)
    @test map(unit -> unit.role, outer_mismatch_product_units.units) ==
          (:outer_mismatch_z_low_slab, :outer_mismatch_z_high_slab)
    @test map(unit -> unit.column_range, outer_mismatch_product_units.units) ==
          (1:49, 50:98)
    @test map(unit -> length(unit.support_indices), outer_mismatch_product_units.units) ==
          (49, 49)
    outer_mismatch_piece_support_indices = vcat(
        outer_mismatch_product_units.units[1].support_indices,
        outer_mismatch_product_units.units[2].support_indices,
    )
    @test sort(outer_mismatch_piece_support_indices) ==
          sort(outer_mismatch_product_units.fact.support_indices)
    @test map(
        unit -> map(axis -> axis.kind, unit.axes),
        outer_mismatch_product_units.units,
    ) == ((:active, :active, :fixed), (:active, :active, :fixed))
    @test outer_mismatch_product_units.units[1].axes[1].interval == 1:7
    @test outer_mismatch_product_units.units[1].axes[2].interval == 1:7
    @test outer_mismatch_product_units.units[1].axes[3].fixed_index == 1
    @test outer_mismatch_product_units.units[2].axes[1].interval == 1:7
    @test outer_mismatch_product_units.units[2].axes[2].interval == 1:7
    @test outer_mismatch_product_units.units[2].axes[3].fixed_index == 15
    @test all(
        unit -> Matrix(unit.coefficient_matrix) == Matrix{Float64}(I, 49, 49),
        outer_mismatch_product_units.units,
    )
    @test all(
        unit -> unit.axis_function_indices ==
                GaussletBases._nested_product_axis_function_indices(3, 1, 7, 2, 7),
        outer_mismatch_product_units.units,
    )
    @test all(
        equivalence -> equivalence.support_indices_match,
        outer_mismatch_product_units.piece_equivalences,
    )
    @test all(
        equivalence -> equivalence.retained_count_match,
        outer_mismatch_product_units.piece_equivalences,
    )
    @test map(
        equivalence -> equivalence.column_range,
        outer_mismatch_product_units.piece_equivalences,
    ) == (1:49, 50:98)
    @test all(
        equivalence -> equivalence.coefficient_matrix_matches_direct_selector,
        outer_mismatch_product_units.piece_equivalences,
    )
    @test all(
        equivalence -> equivalence.max_parent_coefficient_error == 0.0,
        outer_mismatch_product_units.piece_equivalences,
    )
    @test outer_mismatch_product_units.aggregate_equivalence.support_indices_match
    @test outer_mismatch_product_units.aggregate_equivalence.audited_support_set_match
    @test outer_mismatch_product_units.aggregate_equivalence.retained_count_match
    @test outer_mismatch_product_units.aggregate_equivalence.column_range_partition
    @test outer_mismatch_product_units.aggregate_equivalence.coefficient_matrix_matches_direct_selector
    @test outer_mismatch_product_units.aggregate_equivalence.max_parent_coefficient_error == 0.0
    @test outer_mismatch_product_units.diagnostics.outer_mismatch_only
    @test outer_mismatch_product_units.diagnostics.boundary_slab_set
    @test outer_mismatch_product_units.diagnostics.product_doside_units_created
    @test outer_mismatch_product_units.diagnostics.unit_count == 2
    @test outer_mismatch_product_units.diagnostics.slab_piece_count == 2
    @test !outer_mismatch_product_units.diagnostics.route_descriptor_emitted
    @test !outer_mismatch_product_units.diagnostics.construction_mutated
    @test !outer_mismatch_product_units.diagnostics.sidecar_installation
    @test !outer_mismatch_product_units.diagnostics.packet_adoption
    @test !outer_mismatch_product_units.diagnostics.fixed_block_construction_changed
    @test !outer_mismatch_product_units.diagnostics.qwhamiltonian_changed
    @test !outer_mismatch_product_units.diagnostics.ida_weight_division_allowed
    @test outer_mismatch_product_units.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !outer_mismatch_product_units.diagnostics.input_fact_raw_product_box_operator_contract
    @test outer_mismatch_product_units.diagnostics.created_units_raw_product_box_operator_contract
    @test outer_mismatch_product_units.diagnostics.descriptor_piece_order_defines_columns
    @test outer_mismatch_product_units.diagnostics.audited_support_checked_as_set
    @test outer_mismatch_product_units.diagnostics.product_box_construction_rule_available
    @test !outer_mismatch_product_units.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    outer_mismatch_safe_term_metrics = _pqs_axis_metrics(bundles)
    outer_mismatch_safe_terms = (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    outer_mismatch_safe_term_comparison =
        CCPM._pqs_outer_mismatch_safe_term_operator_comparison(
            pqs_construction,
            outer_mismatch_safe_term_metrics,
        )
    @test outer_mismatch_safe_term_comparison.object_kind ==
          :pqs_outer_mismatch_safe_term_operator_comparison
    @test outer_mismatch_safe_term_comparison.status == :private_diagnostic_only
    @test outer_mismatch_safe_term_comparison.terms == outer_mismatch_safe_terms
    @test length(outer_mismatch_safe_term_comparison.fixture.units) == 2
    @test outer_mismatch_safe_term_comparison.max_block_error <= 1.0e-12
    @test outer_mismatch_safe_term_comparison.diagnostics.source ==
          :pqs_outer_mismatch_safe_term_operator_comparison
    @test outer_mismatch_safe_term_comparison.diagnostics.outer_mismatch_only
    @test outer_mismatch_safe_term_comparison.diagnostics.boundary_slab_set
    @test outer_mismatch_safe_term_comparison.diagnostics.private_diagnostic_only
    @test outer_mismatch_safe_term_comparison.diagnostics.terms_checked ==
          outer_mismatch_safe_terms
    @test outer_mismatch_safe_term_comparison.diagnostics.supported_terms ==
          outer_mismatch_safe_terms
    @test :weights in outer_mismatch_safe_term_comparison.diagnostics.unsupported_terms
    @test outer_mismatch_safe_term_comparison.diagnostics.product_path ==
          :_product_doside_source_box_reference_block
    @test outer_mismatch_safe_term_comparison.diagnostics.direct_oracle_path ==
          :support_local_direct_selector_contract_pair_block
    @test outer_mismatch_safe_term_comparison.diagnostics.product_doside_units_created
    @test outer_mismatch_safe_term_comparison.diagnostics.complete_slab_set_block_assembled
    @test outer_mismatch_safe_term_comparison.diagnostics.cross_slab_blocks_included
    @test outer_mismatch_safe_term_comparison.diagnostics.direct_support_oracle_compared
    @test !outer_mismatch_safe_term_comparison.diagnostics.route_descriptor_emitted
    @test !outer_mismatch_safe_term_comparison.diagnostics.construction_mutated
    @test !outer_mismatch_safe_term_comparison.diagnostics.sidecar_installation
    @test !outer_mismatch_safe_term_comparison.diagnostics.packet_adoption
    @test !outer_mismatch_safe_term_comparison.diagnostics.fixed_block_construction_changed
    @test !outer_mismatch_safe_term_comparison.diagnostics.qwhamiltonian_changed
    @test !outer_mismatch_safe_term_comparison.diagnostics.ida_weight_division_allowed
    @test outer_mismatch_safe_term_comparison.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !outer_mismatch_safe_term_comparison.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test outer_mismatch_safe_term_comparison.diagnostics.operator_factor_source ==
          :explicit_metric_operator_data
    @test outer_mismatch_safe_term_comparison.diagnostics.operator_metric_sources ==
          (:nested_pgdg_axis, :nested_pgdg_axis, :nested_pgdg_axis)
    @test !outer_mismatch_safe_term_comparison.diagnostics.input_metric_operator_data_pgdg_checked
    @test !outer_mismatch_safe_term_comparison.diagnostics.pgdg_analytic_operator_provenance_claimed
    @test !outer_mismatch_safe_term_comparison.diagnostics.numerical_reference_fallback
    @test outer_mismatch_safe_term_comparison.diagnostics.product_source_box_reference_compared
    @test outer_mismatch_safe_term_comparison.diagnostics.direct_support_oracle_entries_built
    @test outer_mismatch_safe_term_comparison.diagnostics.retained_count == 98
    @test outer_mismatch_safe_term_comparison.diagnostics.support_count == 98
    @test outer_mismatch_safe_term_comparison.diagnostics.unit_count == 2
    @test outer_mismatch_safe_term_comparison.diagnostics.max_block_error <= 1.0e-12
    @test outer_mismatch_safe_term_comparison.diagnostics.output_finite
    for term in outer_mismatch_safe_terms
        product_block = outer_mismatch_safe_term_comparison.product_blocks[term]
        oracle_block = outer_mismatch_safe_term_comparison.direct_oracle_blocks[term]
        @test size(product_block) == (98, 98)
        @test size(oracle_block) == (98, 98)
        @test all(isfinite, product_block)
        @test all(isfinite, oracle_block)
        @test product_block ≈ oracle_block atol = 1.0e-12 rtol = 0.0
        @test outer_mismatch_safe_term_comparison.term_errors[term] <= 1.0e-12
        @test outer_mismatch_safe_term_comparison.diagnostics.pair_block_counts[term] == 4
        @test outer_mismatch_safe_term_comparison.diagnostics.cross_slab_pair_counts[term] == 2
        @test length(outer_mismatch_safe_term_comparison.product_references[term].pair_references) == 4
    end
    @test_throws ArgumentError CCPM._pqs_outer_mismatch_safe_term_operator_comparison(
        pqs_construction,
        outer_mismatch_safe_term_metrics;
        terms = (:weights,),
    )
    left_atom_fact = pqs_retained_facts_by_role[:left_atom_box]
    right_atom_fact = pqs_retained_facts_by_role[:right_atom_box]
    @test left_atom_fact.classification ==
          :needs_direct_support_retained_unit_kind
    @test right_atom_fact.classification ==
          :needs_direct_support_retained_unit_kind
    @test !left_atom_fact.product_doside_unit
    @test !right_atom_fact.product_doside_unit
    @test !left_atom_fact.raw_product_box_operator_contract
    @test !right_atom_fact.raw_product_box_operator_contract
    @test !left_atom_fact.product_box_construction_rule_available
    @test !right_atom_fact.product_box_construction_rule_available
    @test left_atom_fact.safe_term_capability == :support_local_reference_only
    @test right_atom_fact.safe_term_capability == :support_local_reference_only
    atom_box_support_dense_units =
        CCPM._pqs_atom_box_support_dense_units(pqs_construction)
    @test atom_box_support_dense_units.object_kind ==
          :pqs_atom_box_support_dense_units_fixture
    @test atom_box_support_dense_units.status == :private_diagnostic_only
    @test length(atom_box_support_dense_units.units) == 2
    @test map(unit -> unit.role, atom_box_support_dense_units.units) ==
          (:left_atom_box, :right_atom_box)
    @test map(unit -> unit.kind, atom_box_support_dense_units.units) ==
          (:support_dense, :support_dense)
    @test map(unit -> unit.column_range, atom_box_support_dense_units.units) ==
          (99:223, 224:348)
    @test map(unit -> length(unit.support_indices), atom_box_support_dense_units.units) ==
          (125, 125)
    @test map(unit -> size(unit.coefficient_matrix), atom_box_support_dense_units.units) ==
          ((125, 125), (125, 125))
    @test all(
        unit -> map(axis -> axis.kind, unit.axes) == (:fixed, :fixed, :fixed),
        atom_box_support_dense_units.units,
    )
    @test all(
        unit -> all(index -> index == (1, 1, 1), unit.axis_function_indices),
        atom_box_support_dense_units.units,
    )
    @test all(
        equivalence -> equivalence.support_indices_match,
        atom_box_support_dense_units.equivalences,
    )
    @test all(
        equivalence -> equivalence.audited_support_set_match,
        atom_box_support_dense_units.equivalences,
    )
    @test all(
        equivalence -> equivalence.retained_count_match,
        atom_box_support_dense_units.equivalences,
    )
    @test all(
        equivalence -> equivalence.coefficient_matrix_matches_direct_support,
        atom_box_support_dense_units.equivalences,
    )
    @test all(
        equivalence -> equivalence.max_parent_coefficient_error == 0.0,
        atom_box_support_dense_units.equivalences,
    )
    @test all(
        equivalence -> equivalence.local_identity_error <= 1.0e-12,
        atom_box_support_dense_units.equivalences,
    )
    @test atom_box_support_dense_units.diagnostics.atom_box_only
    @test atom_box_support_dense_units.diagnostics.support_dense_direct_support_units_created
    @test !atom_box_support_dense_units.diagnostics.product_doside_units_created
    @test !atom_box_support_dense_units.diagnostics.raw_product_box_operator_contract
    @test atom_box_support_dense_units.diagnostics.support_local_reference_only
    @test !atom_box_support_dense_units.diagnostics.product_box_construction_rule_available
    @test !atom_box_support_dense_units.diagnostics.route_descriptor_emitted
    @test !atom_box_support_dense_units.diagnostics.construction_mutated
    @test !atom_box_support_dense_units.diagnostics.sidecar_installation
    @test !atom_box_support_dense_units.diagnostics.packet_adoption
    @test !atom_box_support_dense_units.diagnostics.fixed_block_construction_changed
    @test !atom_box_support_dense_units.diagnostics.qwhamiltonian_changed
    @test !atom_box_support_dense_units.diagnostics.ida_weight_division_allowed
    @test atom_box_support_dense_units.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !atom_box_support_dense_units.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test atom_box_support_dense_units.diagnostics.max_parent_coefficient_error == 0.0
    @test !atom_box_support_dense_units.diagnostics.local_identity_is_product_box_claim
    @test !atom_box_support_dense_units.diagnostics.safe_term_operator_comparison_added
    atom_box_safe_term_metrics = _pqs_axis_metrics(bundles)
    atom_box_safe_terms = (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    atom_box_safe_term_comparison =
        CCPM._pqs_atom_box_safe_term_operator_comparison(
            pqs_construction,
            atom_box_safe_term_metrics,
        )
    @test atom_box_safe_term_comparison.object_kind ==
          :pqs_atom_box_safe_term_operator_comparison
    @test atom_box_safe_term_comparison.status == :private_diagnostic_only
    @test atom_box_safe_term_comparison.terms == atom_box_safe_terms
    @test length(atom_box_safe_term_comparison.fixture.units) == 2
    @test atom_box_safe_term_comparison.max_block_error <= 1.0e-12
    @test atom_box_safe_term_comparison.diagnostics.source ==
          :pqs_atom_box_safe_term_operator_comparison
    @test atom_box_safe_term_comparison.diagnostics.atom_box_only
    @test atom_box_safe_term_comparison.diagnostics.support_dense_direct_support_units_created
    @test atom_box_safe_term_comparison.diagnostics.support_local_fallback_operator_comparison
    @test !atom_box_safe_term_comparison.diagnostics.product_doside_units_created
    @test !atom_box_safe_term_comparison.diagnostics.raw_product_box_operator_contract
    @test !atom_box_safe_term_comparison.diagnostics.product_box_construction_rule_available
    @test atom_box_safe_term_comparison.diagnostics.complete_atom_box_block_assembled
    @test atom_box_safe_term_comparison.diagnostics.cross_atom_blocks_included
    @test atom_box_safe_term_comparison.diagnostics.direct_support_oracle_compared
    @test !atom_box_safe_term_comparison.diagnostics.route_descriptor_emitted
    @test !atom_box_safe_term_comparison.diagnostics.construction_mutated
    @test !atom_box_safe_term_comparison.diagnostics.sidecar_installation
    @test !atom_box_safe_term_comparison.diagnostics.packet_adoption
    @test !atom_box_safe_term_comparison.diagnostics.fixed_block_construction_changed
    @test !atom_box_safe_term_comparison.diagnostics.qwhamiltonian_changed
    @test !atom_box_safe_term_comparison.diagnostics.ida_weight_division_allowed
    @test atom_box_safe_term_comparison.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !atom_box_safe_term_comparison.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test atom_box_safe_term_comparison.diagnostics.operator_factor_source ==
          :explicit_metric_operator_data
    @test atom_box_safe_term_comparison.diagnostics.operator_metric_sources ==
          (:nested_pgdg_axis, :nested_pgdg_axis, :nested_pgdg_axis)
    @test !atom_box_safe_term_comparison.diagnostics.input_metric_operator_data_pgdg_checked
    @test !atom_box_safe_term_comparison.diagnostics.pgdg_analytic_operator_provenance_claimed
    @test !atom_box_safe_term_comparison.diagnostics.numerical_reference_fallback
    @test atom_box_safe_term_comparison.diagnostics.retained_count == 250
    @test atom_box_safe_term_comparison.diagnostics.support_count == 250
    @test atom_box_safe_term_comparison.diagnostics.unit_count == 2
    @test atom_box_safe_term_comparison.diagnostics.max_block_error <= 1.0e-12
    @test atom_box_safe_term_comparison.diagnostics.output_finite
    for term in atom_box_safe_terms
        support_dense_block = atom_box_safe_term_comparison.support_dense_blocks[term]
        oracle_block = atom_box_safe_term_comparison.direct_oracle_blocks[term]
        @test size(support_dense_block) == (250, 250)
        @test size(oracle_block) == (250, 250)
        @test all(isfinite, support_dense_block)
        @test all(isfinite, oracle_block)
        @test support_dense_block ≈ oracle_block atol = 1.0e-12 rtol = 0.0
        @test atom_box_safe_term_comparison.term_errors[term] <= 1.0e-12
        @test atom_box_safe_term_comparison.diagnostics.pair_block_counts[term] == 4
        @test atom_box_safe_term_comparison.diagnostics.cross_atom_pair_counts[term] == 2
    end
    @test_throws ArgumentError CCPM._pqs_atom_box_safe_term_operator_comparison(
        pqs_construction,
        atom_box_safe_term_metrics;
        terms = (:weights,),
    )
    shared_pqs_fact =
        pqs_retained_facts_by_role[:regular_shared_molecular_shell]
    @test shared_pqs_fact.classification == :out_of_scope
    @test shared_pqs_fact.primitive_family == :projected_q_shell
    @test !shared_pqs_fact.product_doside_unit
    @test !shared_pqs_fact.raw_product_box_operator_contract
    @test !shared_pqs_fact.product_box_construction_rule_available
    @test shared_pqs_fact.safe_term_capability == :not_body_retained_unit
    @test shared_pqs_fact.notes.current_single_pqs_descriptor
    @test pqs_retained_unit_audit.diagnostics.private_diagnostic_only
    @test !pqs_retained_unit_audit.diagnostics.descriptor_emitted
    @test !pqs_retained_unit_audit.diagnostics.packet_adoption
    @test !pqs_retained_unit_audit.diagnostics.fixed_block_construction_changed
    @test !pqs_retained_unit_audit.diagnostics.qwhamiltonian_changed
    @test !pqs_retained_unit_audit.diagnostics.sidecar_mutation
    @test !pqs_retained_unit_audit.diagnostics.sidecar_installation
    @test !pqs_retained_unit_audit.diagnostics.direct_support_reinterpreted_as_product_doside
    @test pqs_retained_unit_audit.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !pqs_retained_unit_audit.diagnostics.ida_weight_division_allowed
    @test !pqs_retained_unit_audit.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    contact_cap_product_unit =
        CCPM._pqs_contact_cap_product_doside_unit(pqs_construction)
    @test contact_cap_product_unit.object_kind ==
          :pqs_contact_cap_product_doside_unit_fixture
    @test contact_cap_product_unit.status == :private_diagnostic_only
    @test contact_cap_product_unit.fact.role == :contact_cap
    @test contact_cap_product_unit.fact.classification ==
          :product_box_constructible
    @test !contact_cap_product_unit.fact.raw_product_box_operator_contract
    @test contact_cap_product_unit.fact.product_box_construction_rule_available
    @test contact_cap_product_unit.unit.role == :contact_cap_slab
    @test contact_cap_product_unit.unit.kind == :product_doside
    @test contact_cap_product_unit.unit.column_range == 349:373
    @test contact_cap_product_unit.unit.support_indices ==
          contact_cap_product_unit.fact.support_indices
    @test contact_cap_product_unit.unit.support_states == [
        GaussletBases._cartesian_unflat_index(index, (7, 7, 15))
        for index in contact_cap_product_unit.unit.support_indices
    ]
    @test Matrix(contact_cap_product_unit.unit.coefficient_matrix) ==
          Matrix{Float64}(I, 25, 25)
    @test map(axis -> axis.kind, contact_cap_product_unit.unit.axes) ==
          (:active, :active, :fixed)
    @test contact_cap_product_unit.unit.axes[1].interval == 2:6
    @test contact_cap_product_unit.unit.axes[2].interval == 2:6
    @test contact_cap_product_unit.unit.axes[3].fixed_index == 8
    @test Matrix(contact_cap_product_unit.unit.axes[1].coefficient_matrix) ==
          Matrix{Float64}(I, 5, 5)
    @test Matrix(contact_cap_product_unit.unit.axes[2].coefficient_matrix) ==
          Matrix{Float64}(I, 5, 5)
    @test contact_cap_product_unit.unit.axis_function_indices ==
          GaussletBases._nested_product_axis_function_indices(3, 1, 5, 2, 5)
    @test contact_cap_product_unit.equivalence.support_indices_match
    @test contact_cap_product_unit.equivalence.support_states_match
    @test contact_cap_product_unit.equivalence.retained_count_match
    @test contact_cap_product_unit.equivalence.column_range_match
    @test contact_cap_product_unit.equivalence.coefficient_matrix_matches_direct_selector
    @test contact_cap_product_unit.equivalence.max_parent_coefficient_error == 0.0
    @test contact_cap_product_unit.diagnostics.contact_cap_only
    @test contact_cap_product_unit.diagnostics.product_doside_unit_created
    @test !contact_cap_product_unit.diagnostics.route_descriptor_emitted
    @test !contact_cap_product_unit.diagnostics.construction_mutated
    @test !contact_cap_product_unit.diagnostics.sidecar_installation
    @test !contact_cap_product_unit.diagnostics.packet_adoption
    @test !contact_cap_product_unit.diagnostics.fixed_block_construction_changed
    @test !contact_cap_product_unit.diagnostics.qwhamiltonian_changed
    @test !contact_cap_product_unit.diagnostics.ida_weight_division_allowed
    @test contact_cap_product_unit.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test contact_cap_product_unit.diagnostics.product_box_construction_rule_available
    @test !contact_cap_product_unit.diagnostics.input_fact_raw_product_box_operator_contract
    @test contact_cap_product_unit.diagnostics.created_unit_raw_product_box_operator_contract
    contact_safe_term_metrics = _pqs_axis_metrics(bundles)
    contact_safe_term_comparison =
        CCPM._pqs_contact_cap_safe_term_operator_comparison(
            pqs_construction,
            contact_safe_term_metrics,
        )
    contact_safe_terms = (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    @test contact_safe_term_comparison.object_kind ==
          :pqs_contact_cap_safe_term_operator_comparison
    @test contact_safe_term_comparison.status == :private_diagnostic_only
    @test contact_safe_term_comparison.terms == contact_safe_terms
    @test contact_safe_term_comparison.fixture.unit.kind == :product_doside
    @test contact_safe_term_comparison.fixture.unit.column_range == 349:373
    @test contact_safe_term_comparison.max_block_error <= 1.0e-12
    @test contact_safe_term_comparison.diagnostics.source ==
          :pqs_contact_cap_safe_term_operator_comparison
    @test contact_safe_term_comparison.diagnostics.contact_cap_only
    @test contact_safe_term_comparison.diagnostics.private_diagnostic_only
    @test contact_safe_term_comparison.diagnostics.terms_checked ==
          contact_safe_terms
    @test contact_safe_term_comparison.diagnostics.supported_terms ==
          contact_safe_terms
    @test :weights in contact_safe_term_comparison.diagnostics.unsupported_terms
    @test contact_safe_term_comparison.diagnostics.product_path ==
          :_product_doside_source_box_reference_block
    @test contact_safe_term_comparison.diagnostics.direct_oracle_path ==
          :support_local_direct_selector_contract_pair_block
    @test contact_safe_term_comparison.diagnostics.current_direct_support_selector_compared
    @test contact_safe_term_comparison.diagnostics.product_doside_unit_created
    @test !contact_safe_term_comparison.diagnostics.input_fact_raw_product_box_operator_contract
    @test contact_safe_term_comparison.diagnostics.created_unit_raw_product_box_operator_contract
    @test contact_safe_term_comparison.diagnostics.product_box_construction_rule_available
    @test !contact_safe_term_comparison.diagnostics.route_descriptor_emitted
    @test !contact_safe_term_comparison.diagnostics.construction_mutated
    @test !contact_safe_term_comparison.diagnostics.sidecar_installation
    @test !contact_safe_term_comparison.diagnostics.packet_adoption
    @test !contact_safe_term_comparison.diagnostics.fixed_block_construction_changed
    @test !contact_safe_term_comparison.diagnostics.qwhamiltonian_changed
    @test !contact_safe_term_comparison.diagnostics.ida_weight_division_allowed
    @test contact_safe_term_comparison.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !contact_safe_term_comparison.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test contact_safe_term_comparison.diagnostics.operator_factor_source ==
          :explicit_metric_operator_data
    @test contact_safe_term_comparison.diagnostics.operator_metric_sources ==
          (:nested_pgdg_axis, :nested_pgdg_axis, :nested_pgdg_axis)
    @test !contact_safe_term_comparison.diagnostics.input_metric_operator_data_pgdg_checked
    @test !contact_safe_term_comparison.diagnostics.pgdg_analytic_operator_provenance_claimed
    @test !contact_safe_term_comparison.diagnostics.numerical_reference_fallback
    @test contact_safe_term_comparison.diagnostics.product_source_box_reference_compared
    @test contact_safe_term_comparison.diagnostics.direct_support_oracle_entries_built
    @test contact_safe_term_comparison.diagnostics.retained_count == 25
    @test contact_safe_term_comparison.diagnostics.support_count == 25
    @test contact_safe_term_comparison.diagnostics.column_range == 349:373
    @test contact_safe_term_comparison.diagnostics.max_block_error <= 1.0e-12
    @test contact_safe_term_comparison.diagnostics.output_finite
    for term in contact_safe_terms
        product_block = contact_safe_term_comparison.product_blocks[term]
        oracle_block = contact_safe_term_comparison.direct_oracle_blocks[term]
        @test size(product_block) == (25, 25)
        @test size(oracle_block) == (25, 25)
        @test all(isfinite, product_block)
        @test all(isfinite, oracle_block)
        @test product_block ≈ oracle_block atol = 1.0e-12 rtol = 0.0
        @test contact_safe_term_comparison.term_errors[term] <= 1.0e-12
        @test contact_safe_term_comparison.product_references[term].diagnostics.authoritative_block_compared
    end
    @test_throws ArgumentError CCPM._pqs_contact_cap_safe_term_operator_comparison(
        pqs_construction,
        contact_safe_term_metrics;
        terms = (:weights,),
    )
    current_route_inventory =
        CCPM._pqs_current_route_retained_unit_inventory(pqs_construction)
    @test current_route_inventory.object_kind ==
          :pqs_current_route_retained_unit_inventory_fixture
    @test current_route_inventory.status == :private_diagnostic_only
    @test length(current_route_inventory.units) == 6
    @test current_route_inventory.coverage.ordered_roles == (
        :outer_mismatch_z_low_slab,
        :outer_mismatch_z_high_slab,
        :left_atom_box,
        :right_atom_box,
        :contact_cap_slab,
        :regular_shared_molecular_shell,
    )
    @test current_route_inventory.coverage.ordered_column_ranges == (
        1:49,
        50:98,
        99:223,
        224:348,
        349:373,
        374:487,
    )
    @test current_route_inventory.coverage.first_column == 1
    @test current_route_inventory.coverage.last_column == 487
    @test current_route_inventory.coverage.represented_count == 487
    @test current_route_inventory.coverage.contiguous
    @test current_route_inventory.coverage.non_overlapping
    @test current_route_inventory.coverage.covers_every_column_once
    @test map(unit -> unit.category, current_route_inventory.units) == (
        :product_doside,
        :product_doside,
        :support_dense,
        :support_dense,
        :product_doside,
        :shell_realized_pqs_fixture,
    )
    @test map(unit -> unit.retained_count, current_route_inventory.units) ==
          (49, 49, 125, 125, 25, 114)
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].kind ==
          :projected_q_shell
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].active_representation_stage ==
          :shell_realized_pqs_fixture
    @test !current_route_inventory.by_role[:regular_shared_molecular_shell].raw_product_box_operator_contract
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].safe_term_capability ==
          :support_local_oracle_for_shell_realization
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].support_count == 362
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].retained_count == 114
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].raw_box_auxiliary_metadata.available
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].raw_box_auxiliary_metadata.reference_only
    @test !current_route_inventory.by_role[:regular_shared_molecular_shell].raw_box_auxiliary_metadata.active_current_route_contract
    shared_pqs_unit = current_route_inventory.by_role[:regular_shared_molecular_shell]
    shared_transform_fact = shared_pqs_unit.shell_realization_transform_fact
    @test shared_transform_fact.object_kind ==
          :pqs_current_route_shell_realization_transform_fact
    @test shared_transform_fact.status == :metadata_precursor
    @test shared_transform_fact.representation_stage == :shell_realized_pqs_fixture
    @test shared_transform_fact.source_box.source_mode_dims ==
          shared_pqs_unit.raw_box_auxiliary_metadata.source_mode_dims
    @test shared_transform_fact.source_box.source_mode_count ==
          prod(shared_transform_fact.source_box.source_mode_dims)
    @test shared_transform_fact.boundary_selection.mode_count == shared_pqs_unit.retained_count
    @test shared_transform_fact.shell_projection.support_count == shared_pqs_unit.support_count
    @test shared_transform_fact.shell_projection.matrix_shape ==
          (shared_pqs_unit.support_count, shared_pqs_unit.retained_count)
    @test shared_transform_fact.lowdin_cleanup.transform_shape ==
          (shared_pqs_unit.retained_count, shared_pqs_unit.retained_count)
    @test shared_transform_fact.retained_columns.support_local_coefficient_shape ==
          size(shared_pqs_unit.support_local_coefficient_matrix)
    @test shared_transform_fact.retained_columns.coefficient_matches_descriptor_realization
    @test shared_transform_fact.retained_columns.max_support_local_coefficient_error <=
          1.0e-12
    @test shared_transform_fact.shell_realization.shell_projection_used
    @test shared_transform_fact.shell_realization.lowdin_cleanup_used
    @test !shared_transform_fact.compact_source_space_transform.available
    @test !shared_transform_fact.source_box_operator_application_ready
    @test shared_transform_fact.diagnostics.support_local_oracle_used
    @test shared_transform_fact.diagnostics.shell_row_oracle_only
    @test shared_transform_fact.diagnostics.metadata_precursor
    shared_transform_fact_checked =
        CCPM._pqs_current_route_shell_realization_transform_fact(
            shared_pqs_unit;
            metrics = contact_safe_term_metrics,
        )
    @test shared_transform_fact_checked.shell_realization.isometry_checked
    @test shared_transform_fact_checked.shell_realization.isometry_error <= 1.0e-8
    @test shared_transform_fact_checked.shell_realization.isometric
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].equivalence.coefficient_matrix_matches_active_shell
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].equivalence.max_parent_coefficient_error == 0.0
    @test current_route_inventory.by_role[:left_atom_box].category == :support_dense
    @test current_route_inventory.by_role[:left_atom_box].safe_term_capability ==
          :support_local_fallback_safe_terms
    @test current_route_inventory.by_role[:contact_cap_slab].category == :product_doside
    @test current_route_inventory.by_role[:outer_mismatch_z_low_slab].category ==
          :product_doside
    @test map(policy -> policy.pair_type, current_route_inventory.pair_policies) == (
        :product_product,
        :support_support,
        :support_product,
        :shell_realized_pqs_product,
        :shell_realized_pqs_support,
        :shell_realized_pqs_pqs,
        :raw_box_pqs_helpers,
    )
    @test current_route_inventory.pair_policies[1].policy ==
          :product_doside_source_box_path
    @test current_route_inventory.pair_policies[4].policy ==
          :support_local_oracle_for_shell_realization
    @test !current_route_inventory.pair_policies[4].active_current_route
    @test !current_route_inventory.pair_policies[4].active_algorithmic_policy
    @test !current_route_inventory.pair_policies[4].source_box_algorithm_available
    @test current_route_inventory.pair_policies[4].support_local_oracle_used
    @test current_route_inventory.pair_policies[4].shell_row_oracle_only
    @test !current_route_inventory.pair_policies[end].active_current_route
    @test current_route_inventory.diagnostics.private_diagnostic_only
    @test current_route_inventory.diagnostics.current_route_inventory
    @test !current_route_inventory.diagnostics.route_descriptor_emitted
    @test !current_route_inventory.diagnostics.construction_mutated
    @test !current_route_inventory.diagnostics.sidecar_installation
    @test !current_route_inventory.diagnostics.packet_adoption
    @test !current_route_inventory.diagnostics.fixed_block_construction_changed
    @test !current_route_inventory.diagnostics.qwhamiltonian_changed
    @test !current_route_inventory.diagnostics.ida_weight_division_allowed
    @test current_route_inventory.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !current_route_inventory.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test current_route_inventory.diagnostics.shared_pqs_active_representation ==
          :shell_realized_pqs_fixture
    @test !current_route_inventory.diagnostics.shared_pqs_raw_box_operator_contract
    @test current_route_inventory.diagnostics.shared_pqs_shell_realization_transform_fact_count == 1
    @test current_route_inventory.diagnostics.shared_pqs_source_box_operator_application_ready_count == 0
    @test current_route_inventory.diagnostics.raw_box_pqs_auxiliary_reference_available
    @test !current_route_inventory.diagnostics.whole_route_safe_term_matrix_consumer
    @test current_route_inventory.diagnostics.fixed_dimension == 487
    @test current_route_inventory.diagnostics.coverage_complete
    current_route_pair_inventory =
        CCPM._pqs_current_route_retained_pair_inventory(current_route_inventory)
    @test current_route_pair_inventory.object_kind ==
          :pqs_current_route_retained_pair_inventory_fixture
    @test current_route_pair_inventory.status == :private_diagnostic_only
    @test current_route_pair_inventory.unit_inventory === current_route_inventory
    @test length(current_route_pair_inventory.pairs) == 21
    expected_pair_roles = Tuple(
        (current_route_inventory.units[left].role, current_route_inventory.units[right].role)
        for left in 1:length(current_route_inventory.units)
        for right in left:length(current_route_inventory.units)
    )
    expected_pair_shapes = Tuple(
        (
            current_route_inventory.units[left].retained_count,
            current_route_inventory.units[right].retained_count,
        ) for left in 1:length(current_route_inventory.units)
        for right in left:length(current_route_inventory.units)
    )
    @test map(
        pair -> (pair.left_role, pair.right_role),
        current_route_pair_inventory.pairs,
    ) == expected_pair_roles
    @test map(pair -> pair.pair_shape, current_route_pair_inventory.pairs) ==
          expected_pair_shapes
    @test current_route_pair_inventory.counts.pair_count == 21
    @test current_route_pair_inventory.counts.product_product == 6
    @test current_route_pair_inventory.counts.support_support == 3
    @test current_route_pair_inventory.counts.support_product == 6
    @test current_route_pair_inventory.counts.shell_realized_pqs_product == 3
    @test current_route_pair_inventory.counts.shell_realized_pqs_support == 2
    @test current_route_pair_inventory.counts.shell_realized_pqs_pqs == 1
    @test current_route_pair_inventory.counts.raw_box_pqs_active == 0
    @test current_route_pair_inventory.counts.active_algorithmic_policy == 15
    @test current_route_pair_inventory.counts.source_box_algorithm_available == 6
    @test current_route_pair_inventory.counts.product_doside_source_box_path == 6
    @test current_route_pair_inventory.counts.support_local_fallback == 9
    @test current_route_pair_inventory.counts.support_local_oracle_for_shell_realization == 6
    @test current_route_pair_inventory.counts.support_local_oracle_pair_count == 6
    @test current_route_pair_inventory.counts.shell_row_oracle_pair_count == 6
    @test current_route_pair_inventory.pairs[1].pair_group == :product_product
    @test current_route_pair_inventory.pairs[1].policy ==
          :product_doside_source_box_path
    @test current_route_pair_inventory.pairs[3].pair_group == :support_product
    @test current_route_pair_inventory.pairs[3].policy == :support_local_fallback
    @test current_route_pair_inventory.pairs[6].pair_group ==
          :shell_realized_pqs_product
    @test current_route_pair_inventory.pairs[6].policy ==
          :support_local_oracle_for_shell_realization
    @test current_route_pair_inventory.pairs[6].shell_row_oracle_only
    @test current_route_pair_inventory.pairs[6].support_local_oracle_used
    @test !current_route_pair_inventory.pairs[6].active_algorithmic_policy
    @test !current_route_pair_inventory.pairs[6].source_box_algorithm_available
    @test current_route_pair_inventory.pairs[15].pair_group ==
          :shell_realized_pqs_support
    @test current_route_pair_inventory.pairs[15].policy ==
          :support_local_oracle_for_shell_realization
    @test current_route_pair_inventory.pairs[end].pair_group ==
          :shell_realized_pqs_pqs
    @test current_route_pair_inventory.pairs[end].policy ==
          :support_local_oracle_for_shell_realization
    @test current_route_pair_inventory.pairs[end].shell_row_oracle_only
    @test current_route_pair_inventory.pairs[end].support_local_oracle_used
    @test count(pair -> pair.active_current_route, current_route_pair_inventory.pairs) == 15
    @test count(pair -> pair.active_algorithmic_policy, current_route_pair_inventory.pairs) == 15
    @test count(pair -> pair.shell_row_oracle_only, current_route_pair_inventory.pairs) == 6
    @test all(
        pair -> !pair.raw_box_pqs_active_pair_policy,
        current_route_pair_inventory.pairs,
    )
    @test current_route_pair_inventory.diagnostics.private_diagnostic_only
    @test current_route_pair_inventory.diagnostics.current_route_pair_inventory
    @test current_route_pair_inventory.diagnostics.unit_inventory_complete
    @test current_route_pair_inventory.diagnostics.upper_triangular_pairs
    @test current_route_pair_inventory.diagnostics.pair_count == 21
    @test current_route_pair_inventory.diagnostics.raw_box_pqs_active_pair_policy_count == 0
    @test current_route_pair_inventory.diagnostics.active_algorithmic_policy_pair_count == 15
    @test current_route_pair_inventory.diagnostics.source_box_algorithm_available_pair_count == 6
    @test current_route_pair_inventory.diagnostics.support_local_oracle_for_shell_realization_pair_count == 6
    @test current_route_pair_inventory.diagnostics.shell_row_oracle_pair_count == 6
    @test current_route_pair_inventory.diagnostics.shell_realized_pqs_pairs_are_oracle_only
    @test !current_route_pair_inventory.diagnostics.route_descriptor_emitted
    @test !current_route_pair_inventory.diagnostics.construction_mutated
    @test !current_route_pair_inventory.diagnostics.sidecar_installation
    @test !current_route_pair_inventory.diagnostics.packet_adoption
    @test !current_route_pair_inventory.diagnostics.fixed_block_construction_changed
    @test !current_route_pair_inventory.diagnostics.qwhamiltonian_changed
    @test !current_route_pair_inventory.diagnostics.ida_weight_division_allowed
    @test current_route_pair_inventory.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !current_route_pair_inventory.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test !current_route_pair_inventory.diagnostics.whole_route_safe_term_matrix_consumer
    be2_inventory_timed = @timed begin
        be2_pqs_construction =
            GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
                be2_basis,
                be2_bundles,
                be2_policy;
                nside = 5,
                term_coefficients = Float64.(expansion.coefficients),
                packet_kernel = :support_reference,
                shared_shell_realization = :projected_q_shell,
            )
        be2_current_route_inventory =
            CCPM._pqs_current_route_retained_unit_inventory(be2_pqs_construction)
        be2_current_route_pair_inventory =
            CCPM._pqs_current_route_retained_pair_inventory(
                be2_current_route_inventory,
            )
        (
            construction = be2_pqs_construction,
            inventory = be2_current_route_inventory,
            pair_inventory = be2_current_route_pair_inventory,
        )
    end
    be2_inventory_payload = be2_inventory_timed.value
    be2_current_route_inventory = be2_inventory_payload.inventory
    be2_current_route_pair_inventory = be2_inventory_payload.pair_inventory
    be2_shared_pqs_units = Tuple(
        unit for unit in be2_current_route_inventory.units
        if unit.category == :shell_realized_pqs_fixture
    )
    @test be2_inventory_timed.time >= 0.0
    @test be2_inventory_timed.bytes >= 0
    @test be2_current_route_inventory.object_kind ==
          :pqs_current_route_retained_unit_inventory_fixture
    @test be2_current_route_inventory.status == :private_diagnostic_only
    @test length(be2_current_route_inventory.units) == 8
    @test be2_current_route_inventory.diagnostics.unit_count == 8
    @test be2_current_route_inventory.diagnostics.shared_pqs_unit_count == 3
    @test be2_current_route_inventory.diagnostics.shared_pqs_roles == (
        :regular_shared_molecular_shell_1,
        :regular_shared_molecular_shell_2,
        :regular_shared_molecular_shell_3,
    )
    @test be2_current_route_inventory.diagnostics.shared_pqs_original_roles ==
          (:regular_shared_molecular_shell, :regular_shared_molecular_shell, :regular_shared_molecular_shell)
    @test map(unit -> unit.role, be2_shared_pqs_units) ==
          be2_current_route_inventory.diagnostics.shared_pqs_roles
    @test map(unit -> unit.original_role, be2_shared_pqs_units) ==
          be2_current_route_inventory.diagnostics.shared_pqs_original_roles
    @test map(unit -> unit.retained_count, be2_shared_pqs_units) == (98, 98, 114)
    @test map(unit -> unit.support_count, be2_shared_pqs_units) ==
          (1738, 1346, 1002)
    @test map(unit -> unit.column_range, be2_shared_pqs_units) ==
          (1174:1271, 1272:1369, 1370:1483)
    @test map(unit -> unit.original_column_range, be2_shared_pqs_units) ==
          (1174:1271, 1272:1369, 1370:1483)
    @test be2_current_route_inventory.diagnostics.shared_pqs_column_ranges ==
          (1174:1271, 1272:1369, 1370:1483)
    @test be2_current_route_inventory.diagnostics.shared_pqs_original_column_ranges ==
          (1174:1271, 1272:1369, 1370:1483)
    @test be2_current_route_inventory.diagnostics.shared_pqs_shell_realization_transform_fact_count == 3
    @test be2_current_route_inventory.diagnostics.shared_pqs_source_box_operator_application_ready_count == 0
    @test be2_current_route_inventory.coverage.first_column == 1
    @test be2_current_route_inventory.coverage.last_column == 1483
    @test be2_current_route_inventory.coverage.represented_count == 1483
    @test be2_current_route_inventory.coverage.covers_every_column_once
    @test be2_current_route_inventory.diagnostics.fixed_dimension == 1483
    @test be2_current_route_inventory.diagnostics.coverage_complete
    @test !be2_current_route_inventory.diagnostics.q4_single_shared_role_order_preserved
    @test all(
        unit -> unit.active_representation_stage == :shell_realized_pqs_fixture,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> unit.safe_term_capability == :support_local_oracle_for_shell_realization,
        be2_shared_pqs_units,
    )
    @test all(unit -> !unit.raw_product_box_operator_contract, be2_shared_pqs_units)
    @test all(
        unit -> unit.raw_box_auxiliary_metadata.available,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> unit.raw_box_auxiliary_metadata.reference_only,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> !unit.raw_box_auxiliary_metadata.active_current_route_contract,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> unit.diagnostics.representation_stage == :shell_realized_pqs_fixture,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> unit.diagnostics.shell_projection_lowdin_realization,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> !unit.diagnostics.raw_product_box_operator_contract,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> !unit.diagnostics.ida_weight_division_allowed,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> unit.diagnostics.retained_weight_semantics ==
                :not_positive_quadrature_weights,
        be2_shared_pqs_units,
    )
    @test length(be2_current_route_inventory.source_fixtures.shared_pqs) == 3
    @test be2_current_route_pair_inventory.object_kind ==
          :pqs_current_route_retained_pair_inventory_fixture
    @test be2_current_route_pair_inventory.unit_inventory ===
          be2_current_route_inventory
    @test length(be2_current_route_pair_inventory.pairs) == 36
    @test be2_current_route_pair_inventory.counts.pair_count == 36
    @test be2_current_route_pair_inventory.diagnostics.pair_count == 36
    @test be2_current_route_pair_inventory.diagnostics.expected_pair_count == 36
    @test be2_current_route_pair_inventory.diagnostics.unit_count == 8
    @test be2_current_route_pair_inventory.counts.raw_box_pqs_active == 0
    @test be2_current_route_pair_inventory.diagnostics.raw_box_pqs_active_pair_policy_count == 0
    @test be2_current_route_pair_inventory.counts.support_local_oracle_for_shell_realization == 21
    @test be2_current_route_pair_inventory.diagnostics.support_local_oracle_for_shell_realization_pair_count == 21
    @test be2_current_route_pair_inventory.diagnostics.shell_row_oracle_pair_count == 21
    @test be2_current_route_pair_inventory.diagnostics.shell_realized_pqs_pairs_are_oracle_only
    @test all(
        pair -> !pair.raw_box_pqs_active_pair_policy,
        be2_current_route_pair_inventory.pairs,
    )
    @test count(pair -> pair.shell_row_oracle_only, be2_current_route_pair_inventory.pairs) == 21
    @test count(pair -> pair.active_algorithmic_policy, be2_current_route_pair_inventory.pairs) == 15
    @test !be2_current_route_pair_inventory.diagnostics.whole_route_safe_term_matrix_consumer
    current_route_safe_terms = CCPM._pqs_current_route_safe_term_matrices(
        pqs_construction,
        contact_safe_term_metrics;
        inventory = current_route_inventory,
        pair_inventory = current_route_pair_inventory,
    )
    @test current_route_safe_terms.object_kind ==
          :pqs_current_route_safe_term_matrices_fixture
    @test current_route_safe_terms.status == :private_diagnostic_only
    @test current_route_safe_terms.terms == (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    @test current_route_safe_terms.global_max_error <= 1.0e-12
    for term in current_route_safe_terms.terms
        @test size(current_route_safe_terms.matrices[term]) == (487, 487)
        @test size(current_route_safe_terms.oracle_matrices[term]) == (487, 487)
        @test all(isfinite, current_route_safe_terms.matrices[term])
        @test all(isfinite, current_route_safe_terms.oracle_matrices[term])
        @test current_route_safe_terms.matrices[term] ≈
              current_route_safe_terms.oracle_matrices[term] atol = 1.0e-12 rtol = 0.0
        @test current_route_safe_terms.term_errors[term] <= 1.0e-12
    end
    @test current_route_safe_terms.diagnostics.private_diagnostic_only
    @test current_route_safe_terms.diagnostics.whole_route_safe_term_matrix_consumer
    @test current_route_safe_terms.diagnostics.retained_dimension == 487
    @test current_route_safe_terms.diagnostics.pair_count == 21
    @test current_route_safe_terms.diagnostics.product_source_box_pair_count == 6
    @test current_route_safe_terms.diagnostics.support_local_fallback_pair_count == 15
    @test current_route_safe_terms.diagnostics.support_local_oracle_for_shell_realization_pair_count == 6
    @test current_route_safe_terms.diagnostics.support_local_oracle_is_debug_validation
    @test current_route_safe_terms.diagnostics.shell_realized_pqs_pairs_use_oracle_not_algorithm
    @test !current_route_safe_terms.diagnostics.source_box_algorithm_available_for_shell_realized_pqs
    @test current_route_safe_terms.diagnostics.raw_box_pqs_active_pair_policy_count == 0
    @test current_route_safe_terms.diagnostics.global_max_error <= 1.0e-12
    @test current_route_safe_terms.diagnostics.finite_output
    @test current_route_safe_terms.diagnostics.support_local_oracle_compared
    @test !current_route_safe_terms.diagnostics.raw_box_pqs_active_policy_used
    @test !current_route_safe_terms.diagnostics.route_descriptor_emitted
    @test !current_route_safe_terms.diagnostics.construction_mutated
    @test !current_route_safe_terms.diagnostics.sidecar_installation
    @test !current_route_safe_terms.diagnostics.packet_adoption
    @test !current_route_safe_terms.diagnostics.fixed_block_construction_changed
    @test !current_route_safe_terms.diagnostics.qwhamiltonian_changed
    @test !current_route_safe_terms.diagnostics.ida_weight_division_allowed
    @test current_route_safe_terms.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !current_route_safe_terms.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test current_route_safe_terms.diagnostics.elapsed_seconds >= 0.0
    @test current_route_safe_terms.diagnostics.allocated_bytes >= 0
    @test_throws ArgumentError CCPM._pqs_current_route_safe_term_matrices(
        pqs_construction,
        contact_safe_term_metrics;
        inventory = current_route_inventory,
        pair_inventory = current_route_pair_inventory,
        terms = (:weights,),
    )
    @test :product_doside_unit in pqs_source_descriptor.non_contracts
    @test :dense_full_parent_fallback in pqs_source_descriptor.non_contracts
    @test pqs_source_descriptor.diagnostics.metadata_only
    @test !pqs_source_descriptor.active_consumption.fixed_block_sidecar_installed
    @test !pqs_source_descriptor.active_consumption.metric_packet_consumes
    @test !pqs_source_descriptor.active_consumption.by_center_consumes
    @test pqs_diagnostics.region_builds[2].retained_count == 114
    @test pqs_diagnostics.fixed_dimension == 487
    @test !isnothing(pqs_construction.sequence.packet)
    @test all(isfinite, pqs_construction.sequence.packet.overlap)
    @test all(isfinite, pqs_construction.sequence.packet.kinetic)
    @test all(isfinite, pqs_construction.sequence.packet.weights)
    @test all(isfinite, pqs_construction.sequence.packet.gaussian_sum)
    @test all(isfinite, pqs_construction.sequence.packet.pair_sum)
    @test norm(pqs_construction.sequence.packet.overlap - I, Inf) < 1.0e-8
    pqs_fixed_block =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_fixed_block(
            pqs_construction,
        )
    @test size(pqs_fixed_block.coefficient_matrix) == (7 * 7 * 15, 487)
    @test all(isfinite, pqs_fixed_block.overlap)
    @test all(isfinite, pqs_fixed_block.kinetic)
    @test all(isfinite, pqs_fixed_block.weights)
    @test all(isfinite, pqs_fixed_block.gaussian_sum)
    @test all(isfinite, pqs_fixed_block.pair_sum)
    @test norm(pqs_fixed_block.overlap - I, Inf) < 1.0e-8
    @test pqs_fixed_block.staged_by_center_sidecar[] === nothing
    current_route_authority_comparison =
        CCPM._pqs_current_route_safe_term_authority_comparison(
            pqs_construction,
            contact_safe_term_metrics;
            inventory = current_route_inventory,
            pair_inventory = current_route_pair_inventory,
            safe_terms = current_route_safe_terms,
            fixed_block = pqs_fixed_block,
        )
    @test current_route_authority_comparison.object_kind ==
          :pqs_current_route_safe_term_authority_comparison_fixture
    @test current_route_authority_comparison.status == :private_diagnostic_only
    @test current_route_authority_comparison.terms == current_route_safe_terms.terms
    @test current_route_authority_comparison.compared_terms ==
          current_route_safe_terms.terms
    @test isempty(current_route_authority_comparison.unavailable_terms)
    @test current_route_authority_comparison.max_authority_error <= 1.0e-8
    for term in current_route_authority_comparison.compared_terms
        @test current_route_authority_comparison.term_errors[term] <= 1.0e-8
        @test current_route_authority_comparison.authority_sources[term] ==
              :fixed_block
        @test current_route_authority_comparison.authority_fields[term] == term
        @test current_route_authority_comparison.authority_shapes[term] == (487, 487)
    end
    @test current_route_authority_comparison.diagnostics.private_diagnostic_only
    @test current_route_authority_comparison.diagnostics.current_route_safe_term_authority_comparison
    @test current_route_authority_comparison.diagnostics.authority_fixed_block_or_sequence_packet_only
    @test current_route_authority_comparison.diagnostics.support_local_oracle_secondary
    @test current_route_authority_comparison.diagnostics.support_local_oracle_global_max_error <=
          1.0e-12
    @test current_route_authority_comparison.diagnostics.compared_term_count ==
          length(current_route_safe_terms.terms)
    @test current_route_authority_comparison.diagnostics.unavailable_term_count == 0
    @test current_route_authority_comparison.diagnostics.max_authority_error <= 1.0e-8
    @test current_route_authority_comparison.diagnostics.finite_output
    @test !current_route_authority_comparison.diagnostics.construction_mutated
    @test !current_route_authority_comparison.diagnostics.sidecar_installation
    @test !current_route_authority_comparison.diagnostics.packet_adoption
    @test !current_route_authority_comparison.diagnostics.fixed_block_construction_changed
    @test !current_route_authority_comparison.diagnostics.qwhamiltonian_changed
    @test !current_route_authority_comparison.diagnostics.ida_weight_division_allowed
    @test current_route_authority_comparison.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !current_route_authority_comparison.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test current_route_authority_comparison.diagnostics.elapsed_seconds >= 0.0
    @test current_route_authority_comparison.diagnostics.allocated_bytes >= 0
    QWCS = GaussletBases.CartesianQWOperatorCarriedSpaces
    pqs_receipt = @test_logs min_level = Logging.Warn QWCS.cartesian_qw_operator_construction_receipt(
        pqs_fixed_block;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
        expansion = expansion,
    )
    pqs_operators = QWCS.qw_operator_construction_receipt_operators(pqs_receipt)
    pqs_receipt_diagnostics = QWCS.qw_operator_construction_receipt_diagnostics(
        pqs_receipt,
    )
    pqs_record_diagnostics = QWCS.qw_operator_construction_record_diagnostics(
        QWCS.qw_operator_construction_receipt_record(pqs_receipt),
    )
    @test pqs_receipt_diagnostics.delegated_to_existing_builder
    @test pqs_receipt_diagnostics.builder == :ordinary_cartesian_qiu_white_operators
    @test pqs_receipt_diagnostics.source_sidecar_agree
    @test isempty(pqs_receipt_diagnostics.mismatch_fields)
    @test pqs_receipt_diagnostics.operator_built
    @test pqs_receipt_diagnostics.gausslet_backend == :pgdg_localized_experimental
    @test pqs_receipt_diagnostics.interaction_treatment == :ggt_nearest
    @test pqs_receipt_diagnostics.nuclear_term_storage == :total_only
    @test pqs_receipt_diagnostics.dense_parent_matrix_used == false
    @test pqs_receipt_diagnostics.heavy_metric_packet_built == false
    @test pqs_receipt_diagnostics.new_hamiltonian_kernel_used == false
    @test pqs_receipt_diagnostics.numerical_outputs_changed == false
    @test pqs_receipt_diagnostics.operator_residual_count == 0
    @test pqs_record_diagnostics.source_sidecar_agree
    @test :carried_has_staged_sidecar in pqs_record_diagnostics.compared_fields
    @test !pqs_record_diagnostics.comparisons.carried_has_staged_sidecar.source_value
    @test !pqs_record_diagnostics.comparisons.carried_has_staged_sidecar.operator_value
    @test isnothing(
        pqs_record_diagnostics.comparisons.carried_staged_by_center_path.source_value,
    )
    @test isnothing(
        pqs_record_diagnostics.comparisons.carried_staged_by_center_path.operator_value,
    )
    @test pqs_record_diagnostics.source_parent_dimension == 7 * 7 * 15
    @test pqs_record_diagnostics.sidecar_parent_dimension == 7 * 7 * 15
    @test pqs_record_diagnostics.source_carried_dimension == 487
    @test pqs_record_diagnostics.sidecar_carried_dimension == 487
    @test pqs_record_diagnostics.source_carried_space_kind == :nested_fixed_block
    @test pqs_record_diagnostics.sidecar_input_kind == :nested_fixed_block_operator
    @test pqs_operators.gausslet_backend == :pgdg_localized_experimental
    @test pqs_operators.interaction_treatment == :ggt_nearest
    @test pqs_operators.nuclear_term_storage == :total_only
    @test pqs_operators.gausslet_count == 487
    @test pqs_operators.residual_count == 0
    @test size(pqs_operators.overlap) == (487, 487)
    @test size(pqs_operators.one_body_hamiltonian) == (487, 487)
    @test size(pqs_operators.interaction_matrix) == (487, 487)
    @test all(isfinite, pqs_operators.overlap)
    @test all(isfinite, pqs_operators.one_body_hamiltonian)
    @test all(isfinite, pqs_operators.interaction_matrix)
    @test norm(pqs_operators.overlap - I, Inf) < 1.0e-8
    @test norm(
        pqs_operators.one_body_hamiltonian - transpose(pqs_operators.one_body_hamiltonian),
        Inf,
    ) < 1.0e-12
    @test norm(
        pqs_operators.interaction_matrix - transpose(pqs_operators.interaction_matrix),
        Inf,
    ) < 1.0e-12
    @test_throws ArgumentError GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
        basis,
        bundles,
        policy;
        nside = 5,
        term_coefficients = Float64.(expansion.coefficients),
        packet_kernel = :factorized_direct,
        shared_shell_realization = :projected_q_shell,
    )

    shared_q5_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        policy.construction_plan;
        q_min = 4,
        atom_q = 4,
        atom_order = 4,
        shared_q = 5,
        shared_order = 5,
        contact_q = 4,
        contact_order = 4,
        outer_mismatch_q = 4,
        outer_mismatch_order = 4,
    )
    shared_q5_construction =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            shared_q5_policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :factorized_direct,
            build_sequence_packet = false,
        )
    shared_q5_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction_diagnostics(
            shared_q5_construction,
        )
    shared_q5_region_builds = shared_q5_diagnostics.region_builds
    @test shared_q5_diagnostics.recipe_label == :mixed_atom_cubic_shared_endcap_panel
    @test shared_q5_diagnostics.active_builder_consumes
    @test shared_q5_diagnostics.parent_dimension == 7 * 7 * 15
    @test shared_q5_diagnostics.fixed_dimension == 523
    @test shared_q5_diagnostics.support_coverage.coverage_ok
    @test [build.q for build in shared_q5_region_builds] == [4, 5, 4, 4, 4]
    @test [build.order for build in shared_q5_region_builds] == [4, 5, 4, 4, 4]
    @test [build.retained_count for build in shared_q5_region_builds] == [98, 150, 125, 125, 25]
    @test [build.column_range for build in shared_q5_region_builds] ==
          [1:98, 374:523, 99:223, 224:348, 349:373]
    @test shared_q5_region_builds[2].metadata.coefficient_contract == :product_doside
    @test shared_q5_region_builds[5].metadata.descriptor_scope == :middle_contact_cap
    @test shared_q5_diagnostics.metadata.q_policy == :atom_growth_endcap_panel_shared_q_variable
    @test !shared_q5_diagnostics.metadata.q4_acceptance_fixture
    @test shared_q5_diagnostics.metadata.non_shared_q_policy == :fixed_q4_order4
    @test shared_q5_diagnostics.metadata.shared_q_values == (5,)
    @test shared_q5_diagnostics.metadata.shared_order_values == (5,)
    @test isnothing(shared_q5_construction.sequence.packet)

    no_packet_readiness =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness(
            construction,
        )
    no_packet_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness_diagnostics(
            no_packet_readiness,
        )
    @test no_packet_diagnostics.parent_dimension == 7 * 7 * 15
    @test no_packet_diagnostics.fixed_dimension == 469
    @test no_packet_diagnostics.support_coverage.coverage_ok
    @test !no_packet_diagnostics.sequence_packet_available
    @test !no_packet_diagnostics.can_produce_fixed_block
    @test no_packet_diagnostics.fixed_block_missing_fields == [:sequence_packet]
    @test !no_packet_diagnostics.can_produce_nested_source
    @test :split_geometry in no_packet_diagnostics.nested_source_missing_fields
    @test no_packet_diagnostics.default_builders_unchanged

    packeted_construction =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :factorized_direct,
            build_sequence_packet = true,
        )
    readiness =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness(
            packeted_construction;
            build_fixed_block = true,
        )
    readiness_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness_diagnostics(
            readiness,
        )
    @test readiness isa
          GaussletBases._BondAlignedDiatomicHighOrderRecipeSourceReadiness3D
    @test readiness_diagnostics.sequence_packet_available
    @test readiness_diagnostics.overlap_available
    @test readiness_diagnostics.weights_available
    @test readiness_diagnostics.overlap_error < 1.0e-8
    @test readiness_diagnostics.can_produce_fixed_block
    @test isempty(readiness_diagnostics.fixed_block_missing_fields)
    @test readiness_diagnostics.fixed_block_built
    @test readiness_diagnostics.fixed_block_backend == :unknown
    @test readiness_diagnostics.fixed_block_dimension == 469
    @test readiness_diagnostics.fixed_block_support_count == 7 * 7 * 15
    @test readiness_diagnostics.staged_by_center_sidecar_available
    @test !readiness_diagnostics.can_produce_nested_source
    @test :child_sequences in readiness_diagnostics.nested_source_missing_fields
    fixed_block =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_fixed_block(
            packeted_construction,
        )
    @test fixed_block isa GaussletBases._NestedFixedBlock3D
    @test size(fixed_block.coefficient_matrix) == (7 * 7 * 15, 469)
    @test norm(fixed_block.overlap - I, Inf) < 1.0e-8
    @test all(isfinite, fixed_block.weights)
    @test fixed_block.staged_by_center_sidecar[] isa
          GaussletBases._CartesianNestedProductStagedByCenterSidecar3D

    receipt = @test_logs min_level = Logging.Warn QWCS.cartesian_qw_operator_construction_receipt(
        fixed_block;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
        expansion = expansion,
    )
    operators = QWCS.qw_operator_construction_receipt_operators(receipt)
    receipt_diagnostics = QWCS.qw_operator_construction_receipt_diagnostics(receipt)
    record_diagnostics = QWCS.qw_operator_construction_record_diagnostics(
        QWCS.qw_operator_construction_receipt_record(receipt),
    )
    @test receipt_diagnostics.delegated_to_existing_builder
    @test receipt_diagnostics.builder == :ordinary_cartesian_qiu_white_operators
    @test receipt_diagnostics.source_sidecar_agree
    @test isempty(receipt_diagnostics.mismatch_fields)
    @test receipt_diagnostics.operator_built
    @test receipt_diagnostics.gausslet_backend == :pgdg_localized_experimental
    @test receipt_diagnostics.interaction_treatment == :ggt_nearest
    @test receipt_diagnostics.nuclear_term_storage == :total_only
    @test receipt_diagnostics.dense_parent_matrix_used == false
    @test receipt_diagnostics.heavy_metric_packet_built == false
    @test receipt_diagnostics.new_hamiltonian_kernel_used == false
    @test receipt_diagnostics.numerical_outputs_changed == false
    @test record_diagnostics.source_sidecar_agree
    @test record_diagnostics.source_parent_dimension == 7 * 7 * 15
    @test record_diagnostics.sidecar_parent_dimension == 7 * 7 * 15
    @test record_diagnostics.source_carried_dimension == 469
    @test record_diagnostics.sidecar_carried_dimension == 469
    @test record_diagnostics.source_carried_space_kind == :nested_fixed_block
    @test record_diagnostics.sidecar_input_kind == :nested_fixed_block_operator
    @test operators.gausslet_backend == :pgdg_localized_experimental
    @test operators.interaction_treatment == :ggt_nearest
    @test operators.gausslet_count == 469
    @test operators.residual_count == 0
    @test size(operators.overlap) == (469, 469)
    @test size(operators.one_body_hamiltonian) == (469, 469)
    @test size(operators.interaction_matrix) == (469, 469)
    @test all(isfinite, operators.overlap)
    @test all(isfinite, operators.one_body_hamiltonian)
    @test all(isfinite, operators.interaction_matrix)
    @test norm(operators.overlap - I, Inf) < 1.0e-8
    @test norm(
        operators.one_body_hamiltonian - transpose(operators.one_body_hamiltonian),
        Inf,
    ) < 1.0e-12
    @test norm(
        operators.interaction_matrix - transpose(operators.interaction_matrix),
        Inf,
    ) < 1.0e-12

    shared_q5_packeted_construction =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            shared_q5_policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :factorized_direct,
            build_sequence_packet = true,
        )
    shared_q5_readiness =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness(
            shared_q5_packeted_construction;
            build_fixed_block = true,
        )
    shared_q5_readiness_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness_diagnostics(
            shared_q5_readiness,
        )
    @test shared_q5_readiness_diagnostics.sequence_packet_available
    @test shared_q5_readiness_diagnostics.overlap_available
    @test shared_q5_readiness_diagnostics.weights_available
    @test shared_q5_readiness_diagnostics.overlap_error < 1.0e-8
    @test shared_q5_readiness_diagnostics.can_produce_fixed_block
    @test isempty(shared_q5_readiness_diagnostics.fixed_block_missing_fields)
    @test shared_q5_readiness_diagnostics.fixed_block_built
    @test shared_q5_readiness_diagnostics.fixed_block_dimension == 523
    @test shared_q5_readiness_diagnostics.fixed_block_support_count == 7 * 7 * 15
    @test shared_q5_readiness_diagnostics.staged_by_center_sidecar_available
    shared_q5_fixed_block =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_fixed_block(
            shared_q5_packeted_construction,
        )
    @test shared_q5_fixed_block isa GaussletBases._NestedFixedBlock3D
    @test size(shared_q5_fixed_block.coefficient_matrix) == (7 * 7 * 15, 523)
    @test norm(shared_q5_fixed_block.overlap - I, Inf) < 1.0e-8
    @test all(isfinite, shared_q5_fixed_block.weights)
    @test shared_q5_fixed_block.staged_by_center_sidecar[] isa
          GaussletBases._CartesianNestedProductStagedByCenterSidecar3D

    shared_q5_receipt = @test_logs min_level = Logging.Warn QWCS.cartesian_qw_operator_construction_receipt(
        shared_q5_fixed_block;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
        expansion = expansion,
    )
    shared_q5_operators = QWCS.qw_operator_construction_receipt_operators(shared_q5_receipt)
    shared_q5_receipt_diagnostics =
        QWCS.qw_operator_construction_receipt_diagnostics(shared_q5_receipt)
    shared_q5_record_diagnostics = QWCS.qw_operator_construction_record_diagnostics(
        QWCS.qw_operator_construction_receipt_record(shared_q5_receipt),
    )
    @test shared_q5_receipt_diagnostics.delegated_to_existing_builder
    @test shared_q5_receipt_diagnostics.source_sidecar_agree
    @test isempty(shared_q5_receipt_diagnostics.mismatch_fields)
    @test shared_q5_receipt_diagnostics.operator_built
    @test shared_q5_receipt_diagnostics.gausslet_backend == :pgdg_localized_experimental
    @test shared_q5_receipt_diagnostics.interaction_treatment == :ggt_nearest
    @test shared_q5_receipt_diagnostics.nuclear_term_storage == :total_only
    @test shared_q5_receipt_diagnostics.dense_parent_matrix_used == false
    @test shared_q5_receipt_diagnostics.heavy_metric_packet_built == false
    @test shared_q5_receipt_diagnostics.new_hamiltonian_kernel_used == false
    @test shared_q5_receipt_diagnostics.numerical_outputs_changed == false
    @test shared_q5_record_diagnostics.source_sidecar_agree
    @test isempty(shared_q5_record_diagnostics.mismatch_fields)
    @test shared_q5_record_diagnostics.source_parent_dimension == 7 * 7 * 15
    @test shared_q5_record_diagnostics.sidecar_parent_dimension == 7 * 7 * 15
    @test shared_q5_record_diagnostics.source_carried_dimension == 523
    @test shared_q5_record_diagnostics.sidecar_carried_dimension == 523
    @test shared_q5_operators.gausslet_backend == :pgdg_localized_experimental
    @test shared_q5_operators.interaction_treatment == :ggt_nearest
    @test shared_q5_operators.gausslet_count == 523
    @test shared_q5_operators.residual_count == 0
    @test size(shared_q5_operators.overlap) == (523, 523)
    @test size(shared_q5_operators.one_body_hamiltonian) == (523, 523)
    @test size(shared_q5_operators.interaction_matrix) == (523, 523)
    @test all(isfinite, shared_q5_operators.overlap)
    @test all(isfinite, shared_q5_operators.one_body_hamiltonian)
    @test all(isfinite, shared_q5_operators.interaction_matrix)
    @test norm(shared_q5_operators.overlap - I, Inf) < 1.0e-8
    @test norm(
        shared_q5_operators.one_body_hamiltonian -
        transpose(shared_q5_operators.one_body_hamiltonian),
        Inf,
    ) < 1.0e-12
    @test norm(
        shared_q5_operators.interaction_matrix -
        transpose(shared_q5_operators.interaction_matrix),
        Inf,
    ) < 1.0e-12

    shared_q6_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        policy.construction_plan;
        q_min = 4,
        atom_q = 4,
        atom_order = 4,
        shared_q = 6,
        shared_order = 6,
        contact_q = 4,
        contact_order = 4,
        outer_mismatch_q = 4,
        outer_mismatch_order = 4,
    )
    shared_q6_packeted_construction =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            shared_q6_policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :factorized_direct,
            build_sequence_packet = true,
        )
    shared_q6_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction_diagnostics(
            shared_q6_packeted_construction,
        )
    shared_q6_region_builds = shared_q6_diagnostics.region_builds
    @test shared_q6_diagnostics.recipe_label == :mixed_atom_cubic_shared_endcap_panel
    @test shared_q6_diagnostics.active_builder_consumes
    @test shared_q6_diagnostics.parent_dimension == 7 * 7 * 15
    @test shared_q6_diagnostics.fixed_dimension == 589
    @test shared_q6_diagnostics.support_coverage.coverage_ok
    @test [build.q for build in shared_q6_region_builds] == [4, 6, 4, 4, 4]
    @test [build.order for build in shared_q6_region_builds] == [4, 6, 4, 4, 4]
    @test [build.retained_count for build in shared_q6_region_builds] ==
          [98, 216, 125, 125, 25]
    @test [build.column_range for build in shared_q6_region_builds] ==
          [1:98, 374:589, 99:223, 224:348, 349:373]
    @test shared_q6_region_builds[2].metadata.coefficient_contract == :product_doside
    @test all(
        build.retained_count == build.built_support_count for
        build in shared_q6_region_builds[[1, 3, 4, 5]]
    )
    @test shared_q6_diagnostics.metadata.non_shared_q_policy == :fixed_q4_order4
    @test shared_q6_diagnostics.metadata.shared_q_values == (6,)
    @test shared_q6_diagnostics.metadata.shared_order_values == (6,)

    shared_q6_readiness =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness(
            shared_q6_packeted_construction;
            build_fixed_block = true,
        )
    shared_q6_readiness_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness_diagnostics(
            shared_q6_readiness,
        )
    @test shared_q6_readiness_diagnostics.sequence_packet_available
    @test shared_q6_readiness_diagnostics.overlap_available
    @test shared_q6_readiness_diagnostics.weights_available
    @test shared_q6_readiness_diagnostics.overlap_error < 1.0e-8
    @test shared_q6_readiness_diagnostics.can_produce_fixed_block
    @test isempty(shared_q6_readiness_diagnostics.fixed_block_missing_fields)
    @test shared_q6_readiness_diagnostics.fixed_block_built
    @test shared_q6_readiness_diagnostics.fixed_block_dimension == 589
    @test shared_q6_readiness_diagnostics.fixed_block_support_count == 7 * 7 * 15
    @test shared_q6_readiness_diagnostics.staged_by_center_sidecar_available
    shared_q6_fixed_block = shared_q6_readiness.fixed_block
    @test shared_q6_fixed_block isa GaussletBases._NestedFixedBlock3D
    @test size(shared_q6_fixed_block.coefficient_matrix) == (7 * 7 * 15, 589)
    @test norm(shared_q6_fixed_block.overlap - I, Inf) < 1.0e-8
    @test all(isfinite, shared_q6_fixed_block.weights)
    @test shared_q6_fixed_block.staged_by_center_sidecar[] isa
          GaussletBases._CartesianNestedProductStagedByCenterSidecar3D

    shared_q6_receipt = @test_logs min_level = Logging.Warn QWCS.cartesian_qw_operator_construction_receipt(
        shared_q6_fixed_block;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
        expansion = expansion,
    )
    shared_q6_operators = QWCS.qw_operator_construction_receipt_operators(shared_q6_receipt)
    shared_q6_receipt_diagnostics =
        QWCS.qw_operator_construction_receipt_diagnostics(shared_q6_receipt)
    shared_q6_record_diagnostics = QWCS.qw_operator_construction_record_diagnostics(
        QWCS.qw_operator_construction_receipt_record(shared_q6_receipt),
    )
    @test shared_q6_receipt_diagnostics.delegated_to_existing_builder
    @test shared_q6_receipt_diagnostics.source_sidecar_agree
    @test isempty(shared_q6_receipt_diagnostics.mismatch_fields)
    @test shared_q6_receipt_diagnostics.operator_built
    @test shared_q6_receipt_diagnostics.gausslet_backend == :pgdg_localized_experimental
    @test shared_q6_receipt_diagnostics.interaction_treatment == :ggt_nearest
    @test shared_q6_receipt_diagnostics.nuclear_term_storage == :total_only
    @test shared_q6_receipt_diagnostics.dense_parent_matrix_used == false
    @test shared_q6_receipt_diagnostics.heavy_metric_packet_built == false
    @test shared_q6_receipt_diagnostics.new_hamiltonian_kernel_used == false
    @test shared_q6_receipt_diagnostics.numerical_outputs_changed == false
    @test shared_q6_record_diagnostics.source_sidecar_agree
    @test isempty(shared_q6_record_diagnostics.mismatch_fields)
    @test shared_q6_record_diagnostics.source_parent_dimension == 7 * 7 * 15
    @test shared_q6_record_diagnostics.sidecar_parent_dimension == 7 * 7 * 15
    @test shared_q6_record_diagnostics.source_carried_dimension == 589
    @test shared_q6_record_diagnostics.sidecar_carried_dimension == 589
    @test shared_q6_operators.gausslet_backend == :pgdg_localized_experimental
    @test shared_q6_operators.interaction_treatment == :ggt_nearest
    @test shared_q6_operators.gausslet_count == 589
    @test shared_q6_operators.residual_count == 0
    @test size(shared_q6_operators.overlap) == (589, 589)
    @test size(shared_q6_operators.one_body_hamiltonian) == (589, 589)
    @test size(shared_q6_operators.interaction_matrix) == (589, 589)
    @test all(isfinite, shared_q6_operators.overlap)
    @test all(isfinite, shared_q6_operators.one_body_hamiltonian)
    @test all(isfinite, shared_q6_operators.interaction_matrix)
    @test norm(shared_q6_operators.overlap - I, Inf) < 1.0e-8
    @test norm(
        shared_q6_operators.one_body_hamiltonian -
        transpose(shared_q6_operators.one_body_hamiltonian),
        Inf,
    ) < 1.0e-12
    @test norm(
        shared_q6_operators.interaction_matrix -
        transpose(shared_q6_operators.interaction_matrix),
        Inf,
    ) < 1.0e-12

    shared_q7_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        policy.construction_plan;
        q_min = 4,
        atom_q = 4,
        atom_order = 4,
        shared_q = 7,
        shared_order = 7,
        contact_q = 4,
        contact_order = 4,
        outer_mismatch_q = 4,
        outer_mismatch_order = 4,
    )
    shared_q7_error = try
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            shared_q7_policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :factorized_direct,
            build_sequence_packet = false,
        )
        nothing
    catch err
        err
    end
    @test shared_q7_error isa ArgumentError
    @test occursin(
        "nested doside retained_count must not exceed the interval size",
        sprint(showerror, shared_q7_error),
    )

    q_row_expectations = (
        (q = 4, fixed_dimension = 469, retained = (98, 96, 125, 125, 25), shared = 96),
        (q = 5, fixed_dimension = 523, retained = (98, 150, 125, 125, 25), shared = 150),
        (q = 6, fixed_dimension = 589, retained = (98, 216, 125, 125, 25), shared = 216),
    )
    multi_shared_region_builds = (
        (role = :outer_mismatch_shared_molecular_shell, retained_count = 12),
        (role = :regular_shared_molecular_shell, retained_count = 24),
        (role = :regular_shared_molecular_shell, retained_count = 36),
        (role = :left_atom_box, retained_count = 64),
    )
    @test QWCS._nested_q_row_shared_retained_counts(multi_shared_region_builds) ==
          (24, 36)
    @test QWCS._nested_q_row_shared_retained_count(multi_shared_region_builds) ===
          nothing
    @test QWCS._nested_q_row_shared_retained_counts(()) == ()
    @test QWCS._nested_q_row_shared_retained_count(()) === nothing
    for expected in q_row_expectations
        q_row_receipt = @test_logs min_level = Logging.Warn QWCS._nested_bond_aligned_diatomic_high_order_q_row_route_receipt(
            basis;
            shared_q = expected.q,
            shared_order = expected.q,
            protected_atom_side_count = 5,
            q_min = 4,
            nside = 5,
            expansion = expansion,
            nuclear_charges = [1.0, 1.0],
            nuclear_term_storage = :total_only,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :pgdg_localized_experimental,
        )
        q_row_diagnostics =
            QWCS._nested_bond_aligned_diatomic_high_order_q_row_route_diagnostics(
                q_row_receipt,
            )
        @test q_row_diagnostics.route_label ==
              :bond_aligned_diatomic_high_order_q_row_route
        @test q_row_diagnostics.shared_q == expected.q
        @test q_row_diagnostics.shared_order == expected.q
        @test q_row_diagnostics.non_shared_q_policy == :fixed_q4_order4
        @test q_row_diagnostics.parent_dimension == 7 * 7 * 15
        @test q_row_diagnostics.fixed_dimension == expected.fixed_dimension
        @test q_row_diagnostics.retained_counts_by_region == expected.retained
        @test q_row_diagnostics.shared_retained_counts == (expected.shared,)
        @test q_row_diagnostics.shared_retained_count == expected.shared
        @test q_row_diagnostics.overlap_error < 1.0e-8
        @test q_row_diagnostics.staged_sidecar_available
        @test q_row_diagnostics.backend == :pgdg_localized_experimental
        @test q_row_diagnostics.residual_count == 0
        @test q_row_diagnostics.gausslet_count == expected.fixed_dimension
        @test q_row_diagnostics.source_sidecar_agree
        @test isempty(q_row_diagnostics.mismatch_fields)
        @test q_row_diagnostics.dense_parent_matrix_used == false
        @test q_row_diagnostics.heavy_metric_packet_built == false
        @test q_row_diagnostics.new_hamiltonian_kernel_used == false
        @test q_row_diagnostics.default_source_builder_changed == false
    end

    for expected in q_row_expectations
        fixture_receipt = @test_logs min_level = Logging.Warn QWCS._nested_bond_aligned_homonuclear_high_order_q_row_fixture_receipt(
            bond_length = 5.0,
            core_spacing = 0.7,
            xmax_parallel = 8.0,
            xmax_transverse = 4.0,
            shared_q = expected.q,
            shared_order = expected.q,
            family = :G10,
            bond_axis = :z,
            nuclear_charge = 1.0,
            reference_spacing = 1.0,
            tail_spacing = 10.0,
            protected_atom_side_count = 5,
            q_min = 4,
            nside = 5,
            expansion = expansion,
            nuclear_term_storage = :total_only,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :pgdg_localized_experimental,
        )
        fixture_diagnostics =
            QWCS._nested_bond_aligned_homonuclear_high_order_q_row_fixture_diagnostics(
                fixture_receipt,
            )
        fixture_provenance =
            QWCS._nested_bond_aligned_homonuclear_high_order_q_row_fixture_provenance(
                fixture_receipt,
            )
        route_diagnostics = fixture_diagnostics.q_row_route_diagnostics
        @test fixture_diagnostics.route_label ==
              :bond_aligned_homonuclear_high_order_q_row_fixture
        @test fixture_diagnostics.receipt_contract ==
              :construct_homonuclear_basis_then_delegate_q_row_route
        @test fixture_diagnostics.basis_constructor ==
              :bond_aligned_homonuclear_qw_basis
        @test fixture_diagnostics.family == :G10
        @test fixture_diagnostics.bond_length == 5.0
        @test fixture_diagnostics.core_spacing == 0.7
        @test fixture_diagnostics.xmax_parallel == 8.0
        @test fixture_diagnostics.xmax_transverse == 4.0
        @test fixture_diagnostics.bond_axis == :z
        @test fixture_diagnostics.reference_spacing == 1.0
        @test fixture_diagnostics.tail_spacing == 10.0
        @test fixture_diagnostics.nuclear_charge == 1.0
        @test fixture_diagnostics.basis_nuclear_charges == (1.0, 1.0)
        @test fixture_receipt.basis.nuclear_charges == [1.0, 1.0]
        @test fixture_diagnostics.basis_nuclei == ((0.0, 0.0, -2.5), (0.0, 0.0, 2.5))
        @test fixture_diagnostics.parent_axis_counts == (7, 7, 15)
        @test fixture_diagnostics.parent_dimension == 7 * 7 * 15
        @test fixture_diagnostics.flat_index_convention.one_based
        @test fixture_diagnostics.flat_index_convention.fastest_axis == :z
        @test fixture_diagnostics.flat_index_convention.slowest_axis == :x
        @test fixture_diagnostics.shared_q == expected.q
        @test fixture_diagnostics.shared_order == expected.q
        @test fixture_diagnostics.q_min == 4
        @test fixture_diagnostics.protected_atom_side_count == 5
        @test fixture_diagnostics.nside == 5
        @test fixture_diagnostics.fixed_dimension == expected.fixed_dimension
        @test fixture_diagnostics.retained_counts_by_region == expected.retained
        @test fixture_diagnostics.shared_retained_counts == (expected.shared,)
        @test fixture_diagnostics.shared_retained_count == expected.shared
        @test route_diagnostics.route_label ==
              :bond_aligned_diatomic_high_order_q_row_route
        @test route_diagnostics.parent_dimension == fixture_diagnostics.parent_dimension
        @test route_diagnostics.fixed_dimension == expected.fixed_dimension
        @test route_diagnostics.backend == :pgdg_localized_experimental
        @test route_diagnostics.source_sidecar_agree
        @test isempty(route_diagnostics.mismatch_fields)
        @test fixture_provenance.charge_policy == :basis_nuclear_charges_only
        @test fixture_provenance.homonuclear_only
        @test !fixture_provenance.heteronuclear_support
        @test !fixture_provenance.public_api
        @test !fixture_provenance.science_validation
    end

    @test_throws ArgumentError QWCS._nested_bond_aligned_diatomic_high_order_q_row_route_receipt(
        basis;
        shared_q = 4,
        expansion = expansion,
        gausslet_backend = :numerical_reference,
    )
    q_row_q7_error = try
        QWCS._nested_bond_aligned_diatomic_high_order_q_row_route_receipt(
            basis;
            shared_q = 7,
            shared_order = 7,
            expansion = expansion,
            gausslet_backend = :pgdg_localized_experimental,
        )
        nothing
    catch err
        err
    end
    @test q_row_q7_error isa ArgumentError
    @test occursin(
        "nested doside retained_count must not exceed the interval size",
        sprint(showerror, q_row_q7_error),
    )

    @test_throws ArgumentError QWCS._nested_bond_aligned_homonuclear_high_order_q_row_fixture_receipt(
        bond_length = 5.0,
        core_spacing = 0.7,
        xmax_parallel = 8.0,
        xmax_transverse = 4.0,
        shared_q = 4,
        expansion = expansion,
        gausslet_backend = :numerical_reference,
    )
    @test_throws UndefKeywordError QWCS._nested_bond_aligned_homonuclear_high_order_q_row_fixture_receipt(
        bond_length = 5.0,
        core_spacing = 0.7,
        xmax_transverse = 4.0,
        shared_q = 4,
        expansion = expansion,
    )
    fixture_q7_error = try
        QWCS._nested_bond_aligned_homonuclear_high_order_q_row_fixture_receipt(
            bond_length = 5.0,
            core_spacing = 0.7,
            xmax_parallel = 8.0,
            xmax_transverse = 4.0,
            shared_q = 7,
            shared_order = 7,
            expansion = expansion,
            gausslet_backend = :pgdg_localized_experimental,
        )
        nothing
    catch err
        err
    end
    @test fixture_q7_error isa ArgumentError
    @test occursin(
        "nested doside retained_count must not exceed the interval size",
        sprint(showerror, fixture_q7_error),
    )

    annulus_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        policy.construction_plan;
        q_min = 4,
        atom_q = 4,
        shared_q = 5,
        contact_q = 4,
        outer_mismatch_q = 4,
        shared_exterior_family = :transverse_annulus_exterior,
    )
    @test_throws ArgumentError GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
        basis,
        bundles,
        annulus_policy;
        nside = 5,
        term_coefficients = Float64.(expansion.coefficients),
        packet_kernel = :factorized_direct,
        build_sequence_packet = false,
    )

    atom_q5_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        policy.construction_plan;
        q_min = 4,
        atom_q = 5,
        atom_order = 5,
        shared_q = 5,
        shared_order = 5,
        contact_q = 4,
        outer_mismatch_q = 4,
    )
    @test_throws ArgumentError GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
        basis,
        bundles,
        atom_q5_policy;
        nside = 5,
        term_coefficients = Float64.(expansion.coefficients),
        packet_kernel = :factorized_direct,
        build_sequence_packet = false,
    )
end

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

@testset "Cartesian nested shell first packet" begin
    function _fixed_a_nested_shell_basis(count::Int; a::Float64 = 0.25, xmax::Float64 = 10.0, tail_spacing::Float64 = 10.0)
        endpoint = (count - 1) / 2
        s = asinh(xmax / a) / (endpoint - xmax / tail_spacing)
        basis = build_basis(MappedUniformBasisSpec(:G10;
            count,
            mapping = AsinhMapping(; a, s, tail_spacing),
            reference_spacing = 1.0,
        ))
        return basis, s
    end

    expansion = coulomb_gaussian_expansion(doacc = false)
    basis, s = _fixed_a_nested_shell_basis(13)
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    pgdg = bundle.pgdg_intermediate
    term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
    interval = 2:(length(basis) - 1)
    shell = GaussletBases._nested_xy_shell_pair(
        bundle,
        interval,
        interval;
        retain_x = 4,
        retain_y = 3,
        term_coefficients,
    )
    packet = shell.packet
    face_low, face_high = shell.faces
    nface = size(face_low.coefficient_matrix, 2)
    direct_overlap = transpose(shell.coefficient_matrix) * shell.coefficient_matrix
    low_z_mean = sum(diag(packet.position_z)[1:nface]) / nface
    high_z_mean = sum(diag(packet.position_z)[(nface + 1):end]) / nface

    @test s > 0.0
    @test shell isa GaussletBases._CartesianNestedXYShell3D
    @test shell.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test size(shell.coefficient_matrix) == (length(basis)^3, 2 * nface)
    @test nface == 9
    @test length(shell.support_indices) == 2 * length(interval)^2
    @test length(shell.support_states) == length(shell.support_indices)
    @test isempty(intersect(face_low.support_indices, face_high.support_indices))
    @test norm(packet.overlap - I, Inf) < 1.0e-10
    @test packet.overlap ≈ direct_overlap atol = 1.0e-12 rtol = 1.0e-12
    @test packet.kinetic ≈ transpose(packet.kinetic) atol = 1.0e-10 rtol = 1.0e-10
    @test packet.position_x ≈ transpose(packet.position_x) atol = 1.0e-10 rtol = 1.0e-10
    @test packet.position_y ≈ transpose(packet.position_y) atol = 1.0e-10 rtol = 1.0e-10
    @test packet.position_z ≈ transpose(packet.position_z) atol = 1.0e-10 rtol = 1.0e-10
    @test packet.x2_x ≈ transpose(packet.x2_x) atol = 1.0e-10 rtol = 1.0e-10
    @test packet.x2_y ≈ transpose(packet.x2_y) atol = 1.0e-10 rtol = 1.0e-10
    @test packet.x2_z ≈ transpose(packet.x2_z) atol = 1.0e-10 rtol = 1.0e-10
    @test !hasproperty(packet, :gaussian_terms)
    @test !hasproperty(packet, :pair_terms)
    @test !hasproperty(packet, :term_storage)
    @test !isnothing(packet.gaussian_sum)
    @test !isnothing(packet.pair_sum)
    support_coefficients = GaussletBases._nested_support_coefficient_slice(
        shell.coefficient_matrix,
        shell.support_indices,
    )
    @test support_coefficients isa SparseMatrixCSC{Float64,Int}
    support_workspace, contraction_scratch = GaussletBases._nested_support_reference_workspaces(
        support_coefficients,
        length(shell.support_indices),
        size(shell.coefficient_matrix, 2),
    )
    @test size(support_workspace) == (0, 0)
    @test size(contraction_scratch) == (0, 0)
    support_weights = GaussletBases._nested_support_weights(shell.support_states, pgdg.weights)
    fixed_weights = vec(transpose(support_coefficients) * support_weights)
    weighted_support_coefficients = support_coefficients .* reshape(1.0 ./ fixed_weights, 1, :)
    @test weighted_support_coefficients isa SparseMatrixCSC{Float64,Int}
    support_axes = GaussletBases._nested_support_axes(shell.support_states)
    gaussian_reference = GaussletBases._nested_support_reference_gaussian_sum(
        support_axes,
        support_coefficients,
        support_workspace,
        contraction_scratch,
        term_coefficients,
        pgdg.gaussian_factor_terms,
        pgdg.gaussian_factor_terms,
        pgdg.gaussian_factor_terms,
    )
    pair_reference = GaussletBases._nested_support_reference_pair_sum(
        support_axes,
        weighted_support_coefficients,
        support_workspace,
        contraction_scratch,
        term_coefficients,
        pgdg.pair_factor_terms_raw,
        pgdg.pair_factor_terms_raw,
        pgdg.pair_factor_terms_raw,
    )
    @test packet.gaussian_sum ≈ gaussian_reference atol = 1.0e-10 rtol = 1.0e-10
    @test packet.pair_sum ≈ pair_reference atol = 1.0e-10 rtol = 1.0e-10
    @test low_z_mean < 0.0
    @test high_z_mean > 0.0
end

@testset "Cartesian nested support immediate contraction helpers" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    basis = build_basis(MappedUniformBasisSpec(:G10;
        count = 13,
        mapping = AsinhMapping(a = 0.25, s = asinh(10.0 / 0.25) / (6.0 - 1.0), tail_spacing = 10.0),
        reference_spacing = 1.0,
    ))
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    pgdg = bundle.pgdg_intermediate
    term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
    interval = 2:(length(basis) - 1)
    shell = GaussletBases._nested_rectangular_shell(
        bundle,
        interval,
        interval,
        interval;
        retain_xy = (4, 3),
        retain_xz = (4, 3),
        retain_yz = (4, 3),
        term_coefficients,
    )
    support_states = shell.support_states
    support_coefficients = Matrix{Float64}(shell.coefficient_matrix[shell.support_indices, :])
    nsupport = length(support_states)
    nfixed = size(support_coefficients, 2)
    workspace = Matrix{Float64}(undef, nsupport, nsupport)
    scratch = Matrix{Float64}(undef, nfixed, nsupport)

    overlap_support = GaussletBases._nested_support_product_matrix(
        support_states,
        pgdg.overlap,
        pgdg.overlap,
        pgdg.overlap,
    )
    overlap_workspace = Matrix{Float64}(undef, nsupport, nsupport)
    GaussletBases._nested_fill_support_product_matrix!(
        overlap_workspace,
        support_states,
        pgdg.overlap,
        pgdg.overlap,
        pgdg.overlap,
    )
    overlap_reference = transpose(support_coefficients) * overlap_support * support_coefficients
    overlap_contracted = Matrix{Float64}(undef, nfixed, nfixed)
    GaussletBases._nested_contract_support_product!(
        overlap_contracted,
        workspace,
        scratch,
        support_states,
        support_coefficients,
        pgdg.overlap,
        pgdg.overlap,
        pgdg.overlap;
        beta = 0.0,
    )

    kinetic_reference_support = GaussletBases._nested_sum_of_support_products(
        support_states,
        (
            (pgdg.kinetic, pgdg.overlap, pgdg.overlap),
            (pgdg.overlap, pgdg.kinetic, pgdg.overlap),
            (pgdg.overlap, pgdg.overlap, pgdg.kinetic),
        ),
    )
    kinetic_reference = transpose(support_coefficients) * kinetic_reference_support * support_coefficients
    kinetic_contracted = Matrix{Float64}(undef, nfixed, nfixed)
    GaussletBases._nested_contract_sum_of_support_products!(
        kinetic_contracted,
        workspace,
        scratch,
        support_states,
        support_coefficients,
        (
            (pgdg.kinetic, pgdg.overlap, pgdg.overlap),
            (pgdg.overlap, pgdg.kinetic, pgdg.overlap),
            (pgdg.overlap, pgdg.overlap, pgdg.kinetic),
        );
        beta = 0.0,
    )

    @test overlap_workspace ≈ overlap_support atol = 0.0 rtol = 0.0
    @test overlap_contracted ≈ overlap_reference atol = 1.0e-10 rtol = 1.0e-10
    @test kinetic_contracted ≈ kinetic_reference atol = 1.0e-10 rtol = 1.0e-10
end

@testset "Cartesian nested shell interface" begin
    function _fixed_a_multi_face_basis(count::Int; a::Float64 = 0.25, xmax::Float64 = 10.0, tail_spacing::Float64 = 10.0)
        endpoint = (count - 1) / 2
        s = asinh(xmax / a) / (endpoint - xmax / tail_spacing)
        basis = build_basis(MappedUniformBasisSpec(:G10;
            count,
            mapping = AsinhMapping(; a, s, tail_spacing),
            reference_spacing = 1.0,
        ))
        return basis, s
    end

    expansion = coulomb_gaussian_expansion(doacc = false)
    basis, s = _fixed_a_multi_face_basis(13)
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    interval = 2:(length(basis) - 1)
    term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
    shell = GaussletBases._nested_rectangular_shell(
        bundle,
        interval,
        interval,
        interval;
        retain_xy = (4, 3),
        retain_xz = (4, 3),
        retain_yz = (4, 3),
        term_coefficients,
    )
    packet = shell.packet
    direct_overlap = transpose(shell.coefficient_matrix) * shell.coefficient_matrix

    @test s > 0.0
    @test shell isa GaussletBases._CartesianNestedShell3D
    @test shell.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test length(shell.faces) == 6
    @test length(shell.face_column_ranges) == 6
    @test size(shell.coefficient_matrix) == (length(basis)^3, 54)
    @test length(shell.support_indices) == 6 * length(interval)^2
    @test length(shell.support_states) == length(shell.support_indices)
    @test all(length(face.support_indices) == length(interval)^2 for face in shell.faces)
    for left in 1:length(shell.faces), right in (left + 1):length(shell.faces)
        @test isempty(intersect(shell.faces[left].support_indices, shell.faces[right].support_indices))
    end
    @test norm(packet.overlap - I, Inf) < 1.0e-10
    @test packet.overlap ≈ direct_overlap atol = 1.0e-12 rtol = 1.0e-12
    @test packet.kinetic ≈ transpose(packet.kinetic) atol = 1.0e-10 rtol = 1.0e-10
    @test packet.position_x ≈ transpose(packet.position_x) atol = 1.0e-10 rtol = 1.0e-10
    @test packet.position_y ≈ transpose(packet.position_y) atol = 1.0e-10 rtol = 1.0e-10
    @test packet.position_z ≈ transpose(packet.position_z) atol = 1.0e-10 rtol = 1.0e-10
    @test !hasproperty(packet, :gaussian_terms)
    @test !hasproperty(packet, :pair_terms)
    @test !hasproperty(packet, :term_storage)
    @test !isnothing(packet.gaussian_sum)
    @test !isnothing(packet.pair_sum)
    support_coefficients = GaussletBases._nested_support_coefficient_slice(
        shell.coefficient_matrix,
        shell.support_indices,
    )
    support_workspace, contraction_scratch = GaussletBases._nested_support_reference_workspaces(
        support_coefficients,
        length(shell.support_indices),
        size(shell.coefficient_matrix, 2),
    )
    support_axes = GaussletBases._nested_support_axes(shell.support_states)
    support_weights = GaussletBases._nested_support_weights(shell.support_states, bundle.pgdg_intermediate.weights)
    fixed_weights = vec(transpose(support_coefficients) * support_weights)
    weighted_support_coefficients = support_coefficients .* reshape(1.0 ./ fixed_weights, 1, :)
    gaussian_reference = GaussletBases._nested_support_reference_gaussian_sum(
        support_axes,
        support_coefficients,
        support_workspace,
        contraction_scratch,
        term_coefficients,
        bundle.pgdg_intermediate.gaussian_factor_terms,
        bundle.pgdg_intermediate.gaussian_factor_terms,
        bundle.pgdg_intermediate.gaussian_factor_terms,
    )
    pair_reference = GaussletBases._nested_support_reference_pair_sum(
        support_axes,
        weighted_support_coefficients,
        support_workspace,
        contraction_scratch,
        term_coefficients,
        bundle.pgdg_intermediate.pair_factor_terms_raw,
        bundle.pgdg_intermediate.pair_factor_terms_raw,
        bundle.pgdg_intermediate.pair_factor_terms_raw,
    )
    @test packet.gaussian_sum ≈ gaussian_reference atol = 1.0e-10 rtol = 1.0e-10
    @test packet.pair_sum ≈ pair_reference atol = 1.0e-10 rtol = 1.0e-10

    for (face, columns) in zip(shell.faces, shell.face_column_ranges)
        mean_value =
            face.fixed_axis == :x ? sum(diag(packet.position_x)[columns]) / length(columns) :
            face.fixed_axis == :y ? sum(diag(packet.position_y)[columns]) / length(columns) :
            sum(diag(packet.position_z)[columns]) / length(columns)
        if face.fixed_side == :low
            @test mean_value < 0.0
        else
            @test mean_value > 0.0
        end
    end
end

@testset "Cartesian nested fixed-block QW-PGDG adapter" begin
    CCS = GaussletBases.CartesianCarriedSpaces
    QWCS = GaussletBases.CartesianQWOperatorCarriedSpaces
    (
        basis,
        bundle,
        shell,
        fixed_block,
        shell_plus_core,
        fixed_block_shell_plus_core,
        legacy,
        baseline,
        nested,
        nested_shell_plus_core,
        baseline_check,
        nested_check,
        nested_shell_plus_core_check,
    ) = _nested_qiu_white_nearest_fixture()

    @test fixed_block isa GaussletBases._NestedFixedBlock3D
    @test fixed_block.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test fixed_block.parent_basis === basis
    @test fixed_block.shell === shell
    @test size(fixed_block.coefficient_matrix) == size(shell.coefficient_matrix)
    @test size(fixed_block.overlap) == (54, 54)
    @test size(fixed_block.kinetic) == (54, 54)
    @test size(fixed_block.fixed_centers) == (54, 3)
    @test length(fixed_block.support_indices) == length(shell.support_indices)
    @test norm(fixed_block.overlap - I, Inf) < 1.0e-10
    @test fixed_block.overlap ≈ shell.packet.overlap atol = 0.0 rtol = 0.0
    @test !hasproperty(fixed_block, :gaussian_terms)
    @test !hasproperty(fixed_block, :pair_terms)
    @test !hasproperty(fixed_block, :term_storage)
    @test fixed_block.gaussian_sum ≈ shell.packet.gaussian_sum atol = 0.0 rtol = 0.0
    @test fixed_block.pair_sum ≈ shell.packet.pair_sum atol = 0.0 rtol = 0.0

    @test baseline.interaction_treatment == :ggt_nearest
    @test nested.interaction_treatment == :ggt_nearest
    @test nested.basis === fixed_block
    @test nested.gausslet_count == size(fixed_block.overlap, 1)
    @test nested.residual_count >= 1
    @test baseline.gausslet_count == length(bundle.pgdg_intermediate.centers)^3
    @test norm(nested.overlap - I, Inf) < 1.0e-10
    @test nested_check.overlap_error < 1.0e-10
    @test isfinite(nested_check.orbital_energy)
    @test isfinite(nested_check.vee_expectation)
    @test nested_check.orbital_energy < 0.0
    @test nested_check.vee_expectation > 0.0
    @test any(orbital.kind == :nested_fixed for orbital in orbitals(nested))
    @test all(startswith(orbital.label, "nf") for orbital in orbitals(nested)[1:nested.gausslet_count])

    overlap_before = copy(nested.overlap)
    one_body_before = copy(nested.one_body_hamiltonian)
    interaction_before = copy(nested.interaction_matrix)
    nested_sidecar = QWCS.cartesian_qw_operator_carried_space_sidecar(nested)
    nested_carried = QWCS.qw_operator_carried_space(nested_sidecar)
    nested_representation = QWCS.qw_operator_basis_representation(nested_sidecar)
    nested_diagnostics = QWCS.qw_operator_carried_space_diagnostics(nested_sidecar)
    @test QWCS.qw_operator_carried_space_provenance(nested_sidecar).input_kind ==
        :nested_fixed_block_operator
    @test nested_carried isa CCS.CartesianCarriedSpace3D
    @test nested_representation isa CartesianBasisRepresentation3D
    @test nested_diagnostics.operator_dimension == size(nested.overlap, 1)
    @test nested_diagnostics.operator_gausslet_count == nested.gausslet_count
    @test nested_diagnostics.operator_residual_count == nested.residual_count
    @test nested_diagnostics.raw_parent_dimension == size(nested.raw_to_final, 1)
    @test nested_diagnostics.carried_dimension == size(fixed_block.overlap, 1)
    @test nested_diagnostics.carried_dimension_matches_operator_gausslet_count
    @test nested_diagnostics.operator_representation_matches_operator_dimension
    @test nested_diagnostics.carried_has_contracted_parent
    @test nested_diagnostics.carried_has_staged_sidecar == false
    @test nested_diagnostics.dense_parent_matrix_used == false
    @test nested_diagnostics.heavy_metric_packet_built == false
    @test nested.overlap == overlap_before
    @test nested.one_body_hamiltonian == one_body_before
    @test nested.interaction_matrix == interaction_before

    nested_build_source = QWCS.cartesian_qw_operator_build_source(
        fixed_block,
        legacy;
        Z = 2.0,
        interaction_treatment = nested.interaction_treatment,
        gausslet_backend = nested.gausslet_backend,
    )
    nested_build_diagnostics =
        QWCS.operator_build_source_diagnostics(nested_build_source)
    @test QWCS.operator_build_source_provenance(nested_build_source).input_kind ==
        :atomic_nested_fixed_block_input
    @test nested_build_source.basis_family == :one_center_atomic
    @test nested_build_source.carried_space_kind == :nested_fixed_block
    @test nested_build_source.gausslet_backend == nested.gausslet_backend
    @test nested_build_source.interaction_treatment == nested.interaction_treatment
    @test nested_build_source.nuclear_term_storage == nested.nuclear_term_storage
    @test nested_build_diagnostics.carried_dimension ==
        nested_diagnostics.carried_dimension
    @test nested_build_diagnostics.carried_has_contracted_parent ==
        nested_diagnostics.carried_has_contracted_parent
    @test nested_build_diagnostics.carried_has_staged_sidecar ==
        nested_diagnostics.carried_has_staged_sidecar
    @test nested_build_diagnostics.dense_parent_matrix_used == false
    @test nested_build_diagnostics.heavy_metric_packet_built == false
    @test nested_build_diagnostics.operator_built == false

    nested_record =
        QWCS.cartesian_qw_operator_construction_record(nested_build_source, nested)
    nested_record_diagnostics =
        QWCS.qw_operator_construction_record_diagnostics(nested_record)
    @test nested_record_diagnostics.source_sidecar_agree
    @test isempty(nested_record_diagnostics.mismatch_fields)
    @test isempty(nested_record_diagnostics.ambiguous_mismatch_fields)
    @test :gausslet_backend in nested_record_diagnostics.compared_fields
    @test :nuclear_term_storage in nested_record_diagnostics.compared_fields
    @test :nuclear_charges in nested_record_diagnostics.compared_fields
    @test :carried_dimension in nested_record_diagnostics.compared_fields
    @test :carried_parent_axis_counts in nested_record_diagnostics.compared_fields
    @test :carried_parent_dimension in nested_record_diagnostics.compared_fields
    @test :carried_representation_basis_kind in nested_record_diagnostics.compared_fields
    @test :carried_representation_parent_kind in nested_record_diagnostics.compared_fields
    @test :carried_representation_final_dimension in nested_record_diagnostics.compared_fields
    @test :carried_axis_sharing in nested_record_diagnostics.compared_fields
    @test :carried_provenance_input_kind in nested_record_diagnostics.compared_fields
    @test :carried_provenance_route_metadata in nested_record_diagnostics.compared_fields
    @test nested_record_diagnostics.source_basis_family == :one_center_atomic
    @test nested_record_diagnostics.source_carried_space_kind == :nested_fixed_block
    @test nested_record_diagnostics.sidecar_input_kind == :nested_fixed_block_operator
    @test nested_record_diagnostics.source_parent_axis_counts ==
        nested_record_diagnostics.sidecar_parent_axis_counts
    @test nested_record_diagnostics.source_parent_dimension ==
        nested_record_diagnostics.sidecar_parent_dimension
    @test :coefficient_matrix_values in
        nested_record_diagnostics.intentionally_not_compared
    @test :overlap_matrix_values in
        nested_record_diagnostics.intentionally_not_compared
    @test nested_record_diagnostics.numerical_outputs_changed == false
    @test nested_record_diagnostics.dense_parent_matrix_used == false
    @test nested_record_diagnostics.heavy_metric_packet_built == false
    @test nested_record_diagnostics.operator_built == false
    @test QWCS.qw_operator_construction_record_sidecar(nested_record) isa
        QWCS.CartesianQWOperatorCarriedSpaceSidecar
    @test QWCS.qw_operator_construction_record_provenance(nested_record).source ==
        :cartesian_qw_operator_construction_record
    @test nested.overlap == overlap_before
    @test nested.one_body_hamiltonian == one_body_before
    @test nested.interaction_matrix == interaction_before

    @test shell_plus_core isa GaussletBases._CartesianNestedShellPlusCore3D
    @test fixed_block_shell_plus_core isa GaussletBases._NestedFixedBlock3D
    @test fixed_block_shell_plus_core.parent_basis === basis
    @test fixed_block_shell_plus_core.shell === shell_plus_core
    inner_len = length(basis) - 2
    @test first(shell_plus_core.core_column_range) == 1
    @test last(shell_plus_core.core_column_range) == length(shell_plus_core.core_indices)
    @test length(shell_plus_core.core_indices) == inner_len^3
    @test isempty(intersect(shell_plus_core.core_indices, shell.support_indices))
    @test size(fixed_block_shell_plus_core.overlap, 1) == length(shell_plus_core.core_indices) + size(shell.coefficient_matrix, 2)
    @test norm(fixed_block_shell_plus_core.overlap - I, Inf) < 1.0e-10
    @test nested_shell_plus_core.interaction_treatment == :ggt_nearest
    @test nested_shell_plus_core.basis === fixed_block_shell_plus_core
    @test nested_shell_plus_core.gausslet_count == size(fixed_block_shell_plus_core.overlap, 1)
    @test nested_shell_plus_core.residual_count >= 1
    @test nested_shell_plus_core_check.overlap_error < 1.0e-10
    @test isfinite(nested_shell_plus_core_check.orbital_energy)
    @test isfinite(nested_shell_plus_core_check.vee_expectation)
    @test nested_shell_plus_core_check.orbital_energy < 0.0
    @test nested_shell_plus_core_check.vee_expectation > 0.0
    @test abs(nested_shell_plus_core_check.vee_expectation - baseline_check.vee_expectation) < abs(nested_check.vee_expectation - baseline_check.vee_expectation)
    @test abs(nested_shell_plus_core_check.orbital_energy - baseline_check.orbital_energy) < abs(nested_check.orbital_energy - baseline_check.orbital_energy)
    @test abs(nested_shell_plus_core_check.vee_expectation - baseline_check.vee_expectation) < 1.0e-4
    @test abs(nested_shell_plus_core_check.orbital_energy - baseline_check.orbital_energy) < 1.0e-4
end

function _one_center_atomic_full_parent_contract_fixture(;
    Z::Float64 = 2.0,
    d::Float64 = 0.15,
    count::Int = 19,
    nside::Int = 7,
    tail_spacing::Float64 = 10.0,
)
    key = Symbol(
        :one_center_atomic_full_parent_contract_fixture,
        round(Int, 1000 * Z),
        round(Int, 1000 * d),
        count,
        nside,
    )
    return _cached_fixture(key, () -> begin
        basis = build_basis(MappedUniformBasisSpec(:G10;
            count,
            mapping = white_lindsey_atomic_mapping(; Z, d, tail_spacing),
            reference_spacing = 1.0,
        ))
        expansion = coulomb_gaussian_expansion(doacc = false)
        sequence = build_one_center_atomic_full_parent_shell_sequence(
            basis;
            exponents = expansion.exponents,
            gausslet_backend = :numerical_reference,
            refinement_levels = 0,
            nside,
        )
        audit = GaussletBases._nested_shell_sequence_contract_audit(sequence, (count, count, count))
        (basis, sequence, audit)
    end)
end

function _one_center_atomic_legacy_profile_contract_fixture(;
    Z::Float64 = 2.0,
    d::Float64 = 0.2,
    count::Int = 15,
    nside::Int = 5,
    working_box::UnitRange{Int} = 2:14,
    tail_spacing::Float64 = 10.0,
)
    key = Symbol(
        :one_center_atomic_legacy_profile_contract_fixture,
        round(Int, 1000 * Z),
        round(Int, 1000 * d),
        count,
        nside,
        first(working_box),
        last(working_box),
    )
    return _cached_fixture(key, () -> begin
        basis = build_basis(MappedUniformBasisSpec(:G10;
            count,
            mapping = white_lindsey_atomic_mapping(; Z, d, tail_spacing),
            reference_spacing = 1.0,
        ))
        expansion = coulomb_gaussian_expansion(doacc = false)
        sequence = build_one_center_atomic_legacy_profile_shell_sequence(
            basis;
            exponents = expansion.exponents,
            gausslet_backend = :numerical_reference,
            refinement_levels = 0,
            working_box,
            nside,
        )
        diagnostics = one_center_atomic_nested_structure_diagnostics(
            sequence;
            parent_side_count = count,
            nside,
        )
        ownership = GaussletBases._nested_shell_sequence_piece_ownership_audit(sequence)
        (basis, sequence, diagnostics, ownership)
    end)
end

function _ne_repo_v6z_sp_basis_text()
    return "#BASIS SET: Ne repo-v6z-sp\n" *
           "Ne    S\n" *
           "      9.024000e+05           5.510000e-06\n" *
           "      1.351000e+05           4.282000e-05\n" *
           "      3.075000e+04           2.251400e-04\n" *
           "      8.710000e+03           9.501600e-04\n" *
           "      2.842000e+03           3.447190e-03\n" *
           "      1.026000e+03           1.112545e-02\n" *
           "      4.001000e+02           3.220568e-02\n" *
           "      1.659000e+02           8.259891e-02\n" *
           "      7.221000e+01           1.799056e-01\n" *
           "      3.266000e+01           3.060521e-01\n" *
           "      1.522000e+01           3.401256e-01\n" *
           "      7.149000e+00           1.761682e-01\n" *
           "      2.957000e+00           2.101527e-02\n" *
           "      1.335000e+00          -5.074500e-04\n" *
           "      5.816000e-01           1.057850e-03\n" *
           "      2.463000e-01          -5.988000e-05\n" *
           "Ne    S\n" *
           "      7.149000e+00           1.000000e+00\n" *
           "Ne    S\n" *
           "      2.957000e+00           1.000000e+00\n" *
           "Ne    S\n" *
           "      1.335000e+00           1.000000e+00\n" *
           "Ne    S\n" *
           "      9.024000e+05          -1.290000e-06\n" *
           "      1.351000e+05          -1.005000e-05\n" *
           "      3.075000e+04          -5.293000e-05\n" *
           "      8.710000e+03          -2.231200e-04\n" *
           "      2.842000e+03          -8.133800e-04\n" *
           "      1.026000e+03          -2.632300e-03\n" *
           "      4.001000e+02          -7.759100e-03\n" *
           "      1.659000e+02          -2.045277e-02\n" *
           "      7.221000e+01          -4.797505e-02\n" *
           "      3.266000e+01          -9.340086e-02\n" *
           "      1.522000e+01          -1.427721e-01\n" *
           "      7.149000e+00          -1.022908e-01\n" *
           "      2.957000e+00           1.587858e-01\n" *
           "      1.335000e+00           4.494079e-01\n" *
           "      5.816000e-01           4.334854e-01\n" *
           "      2.463000e-01           1.215757e-01\n" *
           "Ne    S\n" *
           "      5.816000e-01           1.000000e+00\n" *
           "Ne    S\n" *
           "      2.463000e-01           1.000000e+00\n" *
           "Ne    P\n" *
           "      4.281000e+00           1.000000e+00\n" *
           "Ne    P\n" *
           "      1.915000e+00           1.000000e+00\n" *
           "Ne    P\n" *
           "      8.156000e+02           1.837600e-04\n" *
           "      1.933000e+02           1.585090e-03\n" *
           "      6.260000e+01           8.414640e-03\n" *
           "      2.361000e+01           3.220033e-02\n" *
           "      9.762000e+00           9.396390e-02\n" *
           "      4.281000e+00           2.004808e-01\n" *
           "      1.915000e+00           3.031137e-01\n" *
           "      8.476000e-01           3.297578e-01\n" *
           "      3.660000e-01           2.366743e-01\n" *
           "      1.510000e-01           6.911689e-02\n" *
           "Ne    P\n" *
           "      8.476000e-01           1.000000e+00\n" *
           "Ne    P\n" *
           "      3.660000e-01           1.000000e+00\n" *
           "Ne    P\n" *
           "      1.510000e-01           1.000000e+00\n" *
           "END\n"
end

function _one_center_atomic_legacy_profile_ne_residual_completion_fixture()
    return _cached_fixture(:one_center_atomic_legacy_profile_ne_residual_completion_fixture, () -> begin
        overlap_only_expansion = CoulombGaussianExpansion(
            [0.0],
            [1.0];
            del = 1.0,
            s = 1.0,
            c = 1.0,
            maxu = 1.0,
        )
        mktemp() do path, io
            write(io, _ne_repo_v6z_sp_basis_text())
            close(io)

            basis = build_basis(MappedUniformBasisSpec(:G10;
                count = 29,
                mapping = white_lindsey_atomic_mapping(Z = 10.0, d = 0.03, tail_spacing = 10.0),
                reference_spacing = 1.0,
            ))
            bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
                basis;
                exponents = overlap_only_expansion.exponents,
                center = 0.0,
                backend = :numerical_reference,
            )
            fixed_block = one_center_atomic_legacy_profile_fixed_block(
                bundle;
                working_box = 2:28,
                nside = 7,
            )
            supplement = legacy_atomic_gaussian_supplement(
                "Ne",
                "repo-v6z-sp";
                lmax = 1,
                basisfile = path,
            )
            supplement3d = GaussletBases._atomic_cartesian_shell_supplement_3d(supplement)
            blocks = GaussletBases._qwrg_atomic_cartesian_blocks_3d(
                bundle,
                supplement3d,
                overlap_only_expansion,
            )
            overlap_fg = GaussletBases._qwrg_contract_parent_ga_matrix(
                fixed_block.coefficient_matrix,
                blocks.overlap_ga,
            )
            near_null = diagnose_qwrg_residual_space(
                fixed_block.overlap,
                overlap_fg,
                blocks.overlap_aa;
                keep_policy = :near_null_only,
                keep_abs_tol = GaussletBases._qwrg_atomic_residual_keep_tol(),
                accept_tol = GaussletBases._qwrg_atomic_residual_accept_tol(),
            )
            near_null_data = GaussletBases._qwrg_residual_space(
                fixed_block.overlap,
                overlap_fg,
                blocks.overlap_aa;
                keep_policy = :near_null_only,
                keep_abs_tol = GaussletBases._qwrg_atomic_residual_keep_tol(),
                accept_tol = GaussletBases._qwrg_atomic_residual_accept_tol(),
            )
            near_null_total_basis = size(near_null_data.raw_to_final, 2)
            legacy_alias = diagnose_qwrg_residual_space(
                fixed_block.overlap,
                overlap_fg,
                blocks.overlap_aa;
                keep_policy = :legacy_profile,
                keep_abs_tol = GaussletBases._qwrg_atomic_residual_keep_tol(),
                accept_tol = GaussletBases._qwrg_atomic_residual_accept_tol(),
            )
            return (
                fixed_gausslet_count = size(fixed_block.overlap, 1),
                supplement_count = length(supplement3d.orbitals),
                near_null = near_null,
                near_null_data = near_null_data,
                near_null_total_basis = near_null_total_basis,
                legacy_alias = legacy_alias,
            )
        end
    end)
end

function _one_center_atomic_ns9_legacy_profile_qw_fixture()
    return _cached_fixture(:one_center_atomic_ns9_legacy_profile_qw_fixture, () -> begin
        mktemp() do path, io
            write(io, _ne_repo_v6z_sp_basis_text())
            close(io)

            basis = build_basis(MappedUniformBasisSpec(
                :G10;
                count = 29,
                mapping = white_lindsey_atomic_mapping(Z = 10.0, d = 0.03, tail_spacing = 10.0),
                reference_spacing = 1.0,
            ))
            expansion = coulomb_gaussian_expansion(doacc = false)
            fixed_block = one_center_atomic_legacy_profile_fixed_block(
                basis;
                expansion = expansion,
                working_box = 2:28,
                nside = 9,
            )
            supplement = legacy_atomic_gaussian_supplement(
                "Ne",
                "repo-v6z-sp";
                lmax = 1,
                basisfile = path,
            )
            supplement3d = GaussletBases._atomic_cartesian_shell_supplement_3d(supplement)
            bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
                basis;
                exponents = expansion.exponents,
                center = 0.0,
                backend = :numerical_reference,
            )
            blocks = GaussletBases._qwrg_atomic_cartesian_blocks_3d(
                bundle,
                supplement3d,
                expansion,
            )
            overlap_fg = GaussletBases._qwrg_contract_parent_ga_matrix(
                fixed_block.coefficient_matrix,
                blocks.overlap_ga,
            )
            residual_data = GaussletBases._qwrg_residual_space(
                fixed_block.overlap,
                overlap_fg,
                blocks.overlap_aa;
                keep_policy = :near_null_only,
                keep_abs_tol = GaussletBases._qwrg_atomic_residual_keep_tol(),
                accept_tol = GaussletBases._qwrg_atomic_residual_accept_tol(),
            )
            operators = ordinary_cartesian_qiu_white_operators(
                fixed_block,
                supplement;
                expansion,
                Z = 10.0,
                interaction_treatment = :ggt_nearest,
                residual_keep_policy = :near_null_only,
            )
            return (
                fixed_block = fixed_block,
                residual_data = residual_data,
                operators = operators,
            )
        end
    end)
end

@testset "One-center atomic full-parent nested contract" begin
    basis, sequence, audit = _one_center_atomic_full_parent_contract_fixture()
    count = length(basis)
    diagnostics = one_center_atomic_nested_structure_diagnostics(
        sequence;
        parent_side_count = count,
        nside = 7,
    )
    common_contract = GaussletBases._nested_glass_box_contract(diagnostics)
    count_only_27 = one_center_atomic_nested_structure_diagnostics(27; nside = 7)
    count_only_29 = one_center_atomic_nested_structure_diagnostics(29; nside = 7)
    count_only_5 = one_center_atomic_nested_structure_diagnostics(15; nside = 5)

    @test sequence isa GaussletBases._CartesianNestedShellSequence3D
    @test sequence.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test sequence.working_box == (1:count, 1:count, 1:count)
    @test audit.full_parent_working_box
    @test audit.support_count == count^3
    @test audit.expected_support_count == count^3
    @test audit.missing_row_count == 0
    @test audit.ownership_group_count_min == 1
    @test audit.ownership_group_count_max == 1
    @test audit.ownership_unowned_row_count == 0
    @test audit.ownership_multi_owned_row_count == 0

    @test diagnostics.parent_side_count == count
    @test diagnostics.working_box_side_count == count
    @test diagnostics.nside == 7
    @test diagnostics.core_side_count == 7
    @test diagnostics.shell_layer_count == 6
    @test diagnostics.expected_shell_increment == 7^3 - 5^3
    @test diagnostics.expected_shell_increment == 218
    @test diagnostics.expected_face_retained_count == 6 * 5^2
    @test diagnostics.expected_edge_retained_count == 12 * 5
    @test diagnostics.expected_corner_retained_count == 8
    @test diagnostics.total_face_retained_count == 6 * (6 * 5^2)
    @test diagnostics.total_edge_retained_count == 6 * (12 * 5)
    @test diagnostics.total_corner_retained_count == 6 * 8
    @test diagnostics.total_expected_gausslet_count == 7^3 + 6 * 218
    @test diagnostics.total_actual_gausslet_count == 7^3 + 6 * 218
    @test diagnostics.total_actual_gausslet_count == size(sequence.coefficient_matrix, 2)
    @test diagnostics.layers_match_expected
    @test length(diagnostics.layer_structures) == 6
    @test all(layer.face_retained_count == 150 for layer in diagnostics.layer_structures)
    @test all(layer.edge_retained_count == 60 for layer in diagnostics.layer_structures)
    @test all(layer.corner_retained_count == 8 for layer in diagnostics.layer_structures)
    @test all(layer.retained_dimension == 218 for layer in diagnostics.layer_structures)
    @test common_contract.fixed_dimension == diagnostics.total_actual_gausslet_count
    @test common_contract.contract_audit === diagnostics
    @test common_contract.layer_dimensions ==
          Int[layer.retained_dimension for layer in diagnostics.layer_structures]
    @test common_contract.layer_provenance ==
          GaussletBases._CartesianNestedShellLayerProvenance3D[
              shell.provenance for shell in sequence.shell_layers
          ]
    @test common_contract.leaf_count === nothing
    @test common_contract.layer_provenance[1].source_box == sequence.working_box
    @test common_contract.layer_provenance[1].source_point_count == 19^3 - 17^3
    @test common_contract.layer_provenance[end].next_inner_box == ((7:13), (7:13), (7:13))
    report_text = one_center_atomic_nested_structure_report(diagnostics)
    @test occursin("layer_1_source_box = ", report_text)
    @test occursin("layer_1_next_inner_box = ", report_text)
    @test occursin("layer_1_source_point_count = ", report_text)
    @test !occursin("leaf_count", report_text)

    @test count_only_5.expected_shell_increment == 5^3 - 3^3
    @test count_only_5.expected_shell_increment == 98
    @test count_only_5.layer_structures[1].provenance.source_point_count == 15^3 - 13^3

    @test count_only_27.parent_side_count == 27
    @test count_only_27.working_box_side_count == 27
    @test count_only_27.nside == 7
    @test count_only_27.core_side_count == 7
    @test count_only_27.shell_layer_count == 10
    @test count_only_27.expected_shell_increment == 218
    @test count_only_27.total_expected_gausslet_count == 343 + 10 * 218
    @test count_only_27.total_actual_gausslet_count == 343 + 10 * 218
    @test count_only_27.total_actual_gausslet_count == 2523
    @test count_only_27.layers_match_expected

    @test count_only_29.parent_side_count == 29
    @test count_only_29.working_box_side_count == 29
    @test count_only_29.shell_layer_count == 11
    @test count_only_29.expected_shell_increment == 218
    @test count_only_29.total_actual_gausslet_count == 343 + 11 * 218
end

@testset "One-center atomic legacy-profile nested contract" begin
    basis, sequence, diagnostics, ownership = _one_center_atomic_legacy_profile_contract_fixture()
    common_contract = GaussletBases._nested_glass_box_contract(diagnostics)
    range_groups = UnitRange{Int}[sequence.core_column_range]
    append!(range_groups, sequence.layer_column_ranges)
    support_group_counts = Int[]
    for row in sequence.support_indices
        nzcols = findall(!iszero, @view sequence.coefficient_matrix[row, :])
        touched_groups = 0
        for range in range_groups
            any(col -> col in range, nzcols) && (touched_groups += 1)
        end
        push!(support_group_counts, touched_groups)
    end

    @test sequence isa GaussletBases._CartesianNestedShellSequence3D
    @test sequence.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test sequence.working_box == (2:14, 2:14, 2:14)
    @test length(sequence.support_indices) == 13^3
    @test ownership.min_group_count == 0
    @test ownership.max_group_count == 1
    @test ownership.unowned_row_count == length(basis)^3 - 13^3
    @test ownership.multi_owned_row_count == 0
    @test minimum(support_group_counts) == 1
    @test maximum(support_group_counts) == 1
    @test diagnostics.parent_side_count == length(basis)
    @test diagnostics.working_box_side_count == 13
    @test diagnostics.nside == 5
    @test diagnostics.core_side_count == 5
    @test diagnostics.shell_layer_count == 4
    @test diagnostics.expected_shell_increment == 98
    @test diagnostics.total_actual_gausslet_count == 5^3 + 4 * 98
    @test diagnostics.layers_match_expected
    @test common_contract.fixed_dimension == diagnostics.total_actual_gausslet_count
    @test common_contract.contract_audit === diagnostics
    @test common_contract.layer_dimensions == [98, 98, 98, 98]
    @test common_contract.layer_provenance ==
          GaussletBases._CartesianNestedShellLayerProvenance3D[
              shell.provenance for shell in sequence.shell_layers
          ]
    @test common_contract.leaf_count === nothing
    @test common_contract.layer_provenance[1].source_box == ((2:14), (2:14), (2:14))
    @test common_contract.layer_provenance[end].next_inner_box == ((6:10), (6:10), (6:10))

    count_only_legacy_ne = one_center_atomic_nested_structure_diagnostics(
        29;
        working_box_side_count = 27,
        nside = 7,
    )
    count_only_modern_ne = one_center_atomic_nested_structure_diagnostics(29; nside = 7)
    @test count_only_legacy_ne.parent_side_count == 29
    @test count_only_legacy_ne.working_box_side_count == 27
    @test count_only_legacy_ne.shell_layer_count == 10
    @test count_only_legacy_ne.expected_shell_increment == 218
    @test count_only_legacy_ne.total_actual_gausslet_count == 2523
    @test count_only_modern_ne.working_box_side_count == 29
    @test count_only_modern_ne.shell_layer_count == 11
    @test count_only_modern_ne.total_actual_gausslet_count == 2741
end

@testset "One-center atomic fixed-block timing surface" begin
    function _timing_labels(report::GaussletBases.TimeG.TimingReport)
        labels = String[]
        function _visit(node)
            push!(labels, node.label)
            foreach(_visit, node.children)
        end
        foreach(_visit, report.roots)
        return labels
    end

    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 13,
            mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)

    timed_full = one_center_atomic_full_parent_fixed_block(
        basis;
        exponents = expansion.exponents,
        nside = 5,
        timing = true,
    )
    timed_legacy = one_center_atomic_legacy_profile_fixed_block(
        basis;
        exponents = expansion.exponents,
        working_box = 2:12,
        nside = 5,
        timing = true,
    )

    @test timed_full isa GaussletBases.TimedNestedFixedBlockBuild
    @test timed_legacy isa GaussletBases.TimedNestedFixedBlockBuild
    @test timed_full.timings isa GaussletBases.TimeG.TimingReport
    @test timed_legacy.timings isa GaussletBases.TimeG.TimingReport
    @test timed_full.fixed_block.shell.working_box == (1:13, 1:13, 1:13)
    @test timed_legacy.fixed_block.shell.working_box == (2:12, 2:12, 2:12)
    @test !hasproperty(timed_full.fixed_block, :gaussian_terms)
    @test !hasproperty(timed_full.fixed_block, :pair_terms)
    @test !hasproperty(timed_full.fixed_block, :term_storage)
    @test !isnothing(timed_full.fixed_block.gaussian_sum)
    @test !isnothing(timed_full.fixed_block.pair_sum)
    @test !hasproperty(timed_legacy.fixed_block, :gaussian_terms)
    @test !hasproperty(timed_legacy.fixed_block, :pair_terms)
    @test !hasproperty(timed_legacy.fixed_block, :term_storage)
    @test !isnothing(timed_legacy.fixed_block.gaussian_sum)
    @test !isnothing(timed_legacy.fixed_block.pair_sum)
    @test norm(timed_full.fixed_block.overlap - I, Inf) < 1.0e-10
    @test norm(timed_legacy.fixed_block.overlap - I, Inf) < 1.0e-10
    full_labels = _timing_labels(timed_full.timings)
    legacy_labels = _timing_labels(timed_legacy.timings)
    @test "fixed_block.total" in full_labels
    @test "fixed_block.parent_bundle" in full_labels
    @test "fixed_block.sequence_build" in full_labels
    @test "fixed_block.adapter" in full_labels
    @test "diatomic.packet.total" in full_labels
    @test "diatomic.packet.gaussian_terms" in full_labels
    @test "diatomic.packet.pair_terms" in full_labels
    @test "diatomic.packet.total" in legacy_labels

    full_report = nested_fixed_block_timing_report(timed_full)
    legacy_report = nested_fixed_block_timing_report(timed_legacy.timings)
    @test occursin("fixed_block.total", full_report)
    @test occursin("shell_layer.nonpacket", full_report)
    @test occursin("sequence_merge.nonpacket", full_report)
    @test occursin("diatomic.packet.gaussian_terms", full_report)
    @test occursin("diatomic.packet.pair_terms", full_report)
    @test occursin("diatomic.packet.total", legacy_report)
end

@testset "Global timing macro surface" begin
    old_config = GaussletBases.TimeG._TIMING_CONFIG[]
    try
        @test timing_enabled() == GaussletBases.TimeG.timing_enabled()
        @test timing_live_enabled() == GaussletBases.TimeG.timing_live_enabled()

        reset_timing_report!()
        set_timing!(true)
        set_timing_live!(false)
        set_timing_thresholds!(expand = 0.0, drop = 0.0)

        @timeg "outer" begin
            sleep(0.002)
            @timeg "inner" begin
                sleep(0.001)
            end
        end

        report = current_timing_report()
        @test report isa GaussletBases.TimeG.TimingReport
        @test length(report.roots) == 1
        root = only(report.roots)
        @test root.label == "outer"
        @test root.elapsed_seconds > 0.0
        @test root.self_seconds >= 0.0
        @test root.call_count == 1
        @test length(root.children) == 1
        child = only(root.children)
        @test child.label == "inner"
        @test child.elapsed_seconds > 0.0

        rendered = timing_report(report)
        @test occursin("GaussletBases timing report", rendered)
        @test occursin("outer", rendered)
        @test occursin("inner", rendered)

        live_path, live_io = mktemp()
        close(live_io)
        try
            open(live_path, "w") do io
                redirect_stdout(io) do
                    reset_timing_report!()
                    set_timing!(true)
                    set_timing_live!(true)
                    set_timing_thresholds!(expand = 0.0, drop = 0.0)
                    @timeg "live outer" begin
                        @timeg "live inner" begin
                            sleep(0.001)
                        end
                    end
                end
            end
            live_output = read(live_path, String)
            @test occursin("live outer: ", live_output)
            @test occursin("live inner: ", live_output)
            @test occursin("seconds", live_output)
        finally
            rm(live_path; force = true)
        end

        reset_timing_report!()
        set_timing!(false)
        set_timing_live!(false)
        @timeg "disabled" begin
            sleep(0.001)
        end
        disabled_report = current_timing_report()
        @test isempty(disabled_report.roots)
    finally
        GaussletBases.TimeG._TIMING_CONFIG[] = old_config
        reset_timing_report!()
    end
end

@testset "One-center atomic compact fixed-block term storage" begin
    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 13,
            mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)

    compact_full = one_center_atomic_full_parent_fixed_block(
        basis;
        expansion,
        nside = 5,
    )
    compact_legacy = one_center_atomic_legacy_profile_fixed_block(
        basis;
        expansion,
        working_box = 2:12,
        nside = 5,
    )

    @test !hasproperty(compact_full, :term_storage)
    @test !hasproperty(compact_full, :gaussian_terms)
    @test !hasproperty(compact_full, :pair_terms)
    @test !isnothing(compact_full.gaussian_sum)
    @test !isnothing(compact_full.pair_sum)

    @test !hasproperty(compact_legacy, :term_storage)
    @test !hasproperty(compact_legacy, :gaussian_terms)
    @test !hasproperty(compact_legacy, :pair_terms)
    @test !isnothing(compact_legacy.gaussian_sum)
    @test !isnothing(compact_legacy.pair_sum)

    full_contract = GaussletBases._nested_glass_box_contract(
        one_center_atomic_nested_structure_diagnostics(compact_full; nside = 5),
    )
    legacy_contract = GaussletBases._nested_glass_box_contract(
        one_center_atomic_nested_structure_diagnostics(compact_legacy; nside = 5),
    )
    @test full_contract.fixed_dimension == size(compact_full.overlap, 1)
    @test legacy_contract.fixed_dimension == size(compact_legacy.overlap, 1)
    @test full_contract.leaf_count === nothing
    @test legacy_contract.leaf_count === nothing

end

@testset "Cartesian basis representation for direct-product QW bases" begin
    CP = GaussletBases.CartesianParentGaussletBases
    CCS = GaussletBases.CartesianCarriedSpaces
    QWCS = GaussletBases.CartesianQWOperatorCarriedSpaces
    basis, operators, _check = _bond_aligned_diatomic_qw_fixture()
    representation = basis_representation(basis)
    metadata = basis_metadata(representation)
    chain_basis, _chain_ops, _chain_diagnostics = _bond_aligned_homonuclear_chain_qw_fixture()
    square_basis, _square_ops, _square_diagnostics, _square_check =
        _axis_aligned_homonuclear_square_lattice_qw_fixture()
    chain_representation = basis_representation(chain_basis)
    square_representation = basis_representation(square_basis)

    @test representation isa CartesianBasisRepresentation3D
    @test metadata.basis_kind == :direct_product
    @test metadata.parent_kind == :cartesian_product_basis
    @test metadata.parent_axis_counts == (length(basis.basis_x), length(basis.basis_y), length(basis.basis_z))
    @test metadata.parent_dimension == prod(metadata.parent_axis_counts)
    @test metadata.final_dimension == prod(metadata.parent_axis_counts)
    @test metadata.axis_sharing == :shared_xy
    @test metadata.route_metadata.basis_family == :bond_aligned_diatomic
    @test metadata.route_metadata.bond_axis == basis.bond_axis
    @test metadata.route_metadata.nuclei == basis.nuclei
    @test representation.contraction_kind == :identity
    @test isnothing(representation.coefficient_matrix)
    @test isnothing(representation.support_indices)
    @test isnothing(representation.support_states)
    @test metadata.basis_labels == representation.parent_labels
    @test metadata.basis_centers == representation.parent_centers
    @test size(metadata.basis_centers, 1) == metadata.final_dimension
    @test size(metadata.basis_centers, 2) == 3
    @test chain_representation.metadata.basis_kind == :direct_product
    @test chain_representation.metadata.route_metadata.basis_family ==
        :bond_aligned_homonuclear_chain
    @test square_representation.metadata.basis_kind == :direct_product
    @test square_representation.metadata.route_metadata.basis_family ==
        :axis_aligned_homonuclear_square_lattice

    carried = CCS.cartesian_carried_space(basis)
    chain_carried = CCS.cartesian_carried_space(chain_basis)
    square_carried = CCS.cartesian_carried_space(square_basis)
    @test CCS.carried_space_parent(carried) isa CP.CartesianParentGaussletBasis3D
    @test isnothing(CCS.carried_space_contracted_parent(carried))
    @test CCS.carried_space_representation(carried) isa CartesianBasisRepresentation3D
    @test CCS.carried_space_diagnostics(carried).parent_axis_counts ==
        metadata.parent_axis_counts
    @test CCS.carried_space_diagnostics(carried).parent_dimension ==
        metadata.parent_dimension
    @test CCS.carried_space_diagnostics(carried).representation_final_dimension ==
        metadata.final_dimension
    @test CCS.carried_space_diagnostics(carried).has_contracted_parent == false
    @test CCS.carried_space_diagnostics(carried).dense_parent_matrix_used == false
    @test CCS.carried_space_diagnostics(carried).heavy_metric_packet_built == false
    @test CCS.carried_space_provenance(carried).input_kind ==
        :bond_aligned_diatomic_qw_basis
    @test CCS.carried_space_diagnostics(chain_carried).parent_axis_counts ==
        chain_representation.metadata.parent_axis_counts
    @test CCS.carried_space_provenance(chain_carried).input_kind ==
        :bond_aligned_homonuclear_chain_qw_basis
    @test CCS.carried_space_diagnostics(square_carried).parent_axis_counts ==
        square_representation.metadata.parent_axis_counts
    @test CCS.carried_space_provenance(square_carried).input_kind ==
        :axis_aligned_homonuclear_square_lattice_qw_basis

    overlap_before = copy(operators.overlap)
    one_body_before = copy(operators.one_body_hamiltonian)
    interaction_before = copy(operators.interaction_matrix)
    operator_sidecar = QWCS.cartesian_qw_operator_carried_space_sidecar(operators)
    operator_carried = QWCS.qw_operator_carried_space(operator_sidecar)
    operator_representation = QWCS.qw_operator_basis_representation(operator_sidecar)
    operator_diagnostics = QWCS.qw_operator_carried_space_diagnostics(operator_sidecar)
    @test QWCS.qw_operator_carried_space_provenance(operator_sidecar).input_kind ==
        :bond_aligned_direct_product_operator
    @test operator_carried isa CCS.CartesianCarriedSpace3D
    @test operator_representation isa CartesianBasisRepresentation3D
    @test operator_diagnostics.operator_dimension == size(operators.overlap, 1)
    @test operator_diagnostics.operator_gausslet_count == operators.gausslet_count
    @test operator_diagnostics.operator_residual_count == operators.residual_count
    @test operator_diagnostics.carried_dimension == metadata.final_dimension
    @test operator_diagnostics.carried_dimension_matches_operator_gausslet_count
    @test operator_diagnostics.operator_representation_matches_operator_dimension
    @test operator_diagnostics.carried_has_contracted_parent == false
    @test operator_diagnostics.carried_has_staged_sidecar == false
    @test operator_diagnostics.dense_parent_matrix_used == false
    @test operator_diagnostics.heavy_metric_packet_built == false
    @test operators.overlap == overlap_before
    @test operators.one_body_hamiltonian == one_body_before
    @test operators.interaction_matrix == interaction_before

    build_source = QWCS.cartesian_qw_operator_build_source(
        basis;
        nuclear_charges = operators.nuclear_charges,
        nuclear_term_storage = :auto,
        interaction_treatment = operators.interaction_treatment,
        gausslet_backend = operators.gausslet_backend,
    )
    build_diagnostics = QWCS.operator_build_source_diagnostics(build_source)
    @test QWCS.operator_build_source_carried_space(build_source) isa
        CCS.CartesianCarriedSpace3D
    @test QWCS.operator_build_source_provenance(build_source).input_kind ==
        :bond_aligned_direct_product_input
    @test build_source.basis_family == :bond_aligned_diatomic
    @test build_source.carried_space_kind == :direct_product
    @test build_source.nuclear_charges == operators.nuclear_charges
    @test build_source.gausslet_backend == operators.gausslet_backend
    @test build_source.interaction_treatment == operators.interaction_treatment
    @test build_source.nuclear_term_storage == operators.nuclear_term_storage
    @test build_diagnostics.carried_dimension == operator_diagnostics.carried_dimension
    @test build_diagnostics.carried_has_contracted_parent ==
        operator_diagnostics.carried_has_contracted_parent
    @test build_diagnostics.carried_has_staged_sidecar ==
        operator_diagnostics.carried_has_staged_sidecar
    @test build_diagnostics.dense_parent_matrix_used == false
    @test build_diagnostics.heavy_metric_packet_built == false
    @test build_diagnostics.operator_built == false

    construction_record =
        QWCS.cartesian_qw_operator_construction_record(build_source, operators)
    record_diagnostics =
        QWCS.qw_operator_construction_record_diagnostics(construction_record)
    @test record_diagnostics.source_sidecar_agree
    @test isempty(record_diagnostics.mismatch_fields)
    @test isempty(record_diagnostics.ambiguous_mismatch_fields)
    @test :operator_input_kind in record_diagnostics.compared_fields
    @test :gausslet_backend in record_diagnostics.compared_fields
    @test :interaction_treatment in record_diagnostics.compared_fields
    @test :nuclear_charges in record_diagnostics.compared_fields
    @test :carried_parent_axis_counts in record_diagnostics.compared_fields
    @test :carried_parent_dimension in record_diagnostics.compared_fields
    @test :carried_representation_basis_kind in record_diagnostics.compared_fields
    @test :carried_representation_parent_kind in record_diagnostics.compared_fields
    @test :carried_representation_final_dimension in record_diagnostics.compared_fields
    @test :carried_axis_sharing in record_diagnostics.compared_fields
    @test :carried_provenance_input_kind in record_diagnostics.compared_fields
    @test :carried_provenance_route_metadata in record_diagnostics.compared_fields
    @test :carried_has_staged_sidecar in record_diagnostics.compared_fields
    @test record_diagnostics.source_basis_family == :bond_aligned_diatomic
    @test record_diagnostics.source_carried_space_kind == :direct_product
    @test record_diagnostics.sidecar_input_kind == :bond_aligned_direct_product_operator
    @test record_diagnostics.source_parent_axis_counts ==
        record_diagnostics.sidecar_parent_axis_counts
    @test record_diagnostics.source_parent_dimension ==
        record_diagnostics.sidecar_parent_dimension
    @test :coefficient_matrix_values in
        record_diagnostics.intentionally_not_compared
    @test :interaction_matrix_values in
        record_diagnostics.intentionally_not_compared
    @test record_diagnostics.numerical_outputs_changed == false
    @test record_diagnostics.dense_parent_matrix_used == false
    @test record_diagnostics.heavy_metric_packet_built == false
    @test record_diagnostics.operator_built == false
    @test QWCS.qw_operator_construction_record_sidecar(construction_record) isa
        QWCS.CartesianQWOperatorCarriedSpaceSidecar
    @test QWCS.qw_operator_construction_record_provenance(construction_record).source ==
        :cartesian_qw_operator_construction_record

    construction_receipt = QWCS.cartesian_qw_operator_construction_receipt(
        basis;
        nuclear_charges = operators.nuclear_charges,
        nuclear_term_storage = operators.nuclear_term_storage,
        interaction_treatment = operators.interaction_treatment,
        gausslet_backend = operators.gausslet_backend,
    )
    receipt_operators =
        QWCS.qw_operator_construction_receipt_operators(construction_receipt)
    receipt_record =
        QWCS.qw_operator_construction_receipt_record(construction_receipt)
    receipt_diagnostics =
        QWCS.qw_operator_construction_receipt_diagnostics(construction_receipt)
    @test QWCS.qw_operator_construction_receipt_source(construction_receipt) isa
        QWCS.CartesianOperatorBuildSource3D
    @test receipt_record isa QWCS.CartesianQWOperatorConstructionRecord3D
    @test receipt_diagnostics.delegated_to_existing_builder
    @test receipt_diagnostics.builder == :ordinary_cartesian_qiu_white_operators
    @test receipt_diagnostics.source_sidecar_agree
    @test isempty(receipt_diagnostics.mismatch_fields)
    @test receipt_diagnostics.operator_built
    @test receipt_diagnostics.new_hamiltonian_kernel_used == false
    @test receipt_diagnostics.dense_parent_matrix_used == false
    @test receipt_diagnostics.heavy_metric_packet_built == false
    @test receipt_diagnostics.numerical_outputs_changed == false
    @test QWCS.qw_operator_construction_receipt_provenance(construction_receipt).source ==
        :cartesian_qw_operator_construction_receipt
    @test receipt_operators.overlap == operators.overlap
    @test receipt_operators.one_body_hamiltonian == operators.one_body_hamiltonian
    @test receipt_operators.interaction_matrix == operators.interaction_matrix
    @test receipt_operators.gausslet_count == operators.gausslet_count
    @test receipt_operators.residual_count == operators.residual_count
    @test receipt_operators.gausslet_backend == operators.gausslet_backend
    @test receipt_operators.interaction_treatment == operators.interaction_treatment
    @test receipt_operators.nuclear_term_storage == operators.nuclear_term_storage

    mismatched_source = QWCS.cartesian_qw_operator_build_source(
        basis;
        nuclear_charges = [
            operators.nuclear_charges[1] + 0.5,
            operators.nuclear_charges[2],
        ],
        nuclear_term_storage = :auto,
        interaction_treatment = operators.interaction_treatment,
        gausslet_backend = operators.gausslet_backend,
    )
    mismatch_record = QWCS.cartesian_qw_operator_construction_record(
        mismatched_source,
        operators;
        throw_on_mismatch = false,
    )
    mismatch_diagnostics =
        QWCS.qw_operator_construction_record_diagnostics(mismatch_record)
    @test !mismatch_diagnostics.source_sidecar_agree
    @test :nuclear_charges in mismatch_diagnostics.mismatch_fields
    @test_throws ArgumentError QWCS.cartesian_qw_operator_construction_record(
        mismatched_source,
        operators,
    )
    mismatched_carried_source = QWCS.CartesianOperatorBuildSource3D(
        chain_carried,
        build_source.basis_family,
        build_source.carried_space_kind,
        build_source.nuclei,
        build_source.nuclear_charges,
        build_source.gausslet_backend,
        build_source.requested_gausslet_backend,
        build_source.interaction_treatment,
        build_source.nuclear_term_storage,
        build_source.requested_nuclear_term_storage,
        build_source.capabilities,
        build_source.diagnostics,
        build_source.provenance,
    )
    carried_mismatch_record = QWCS.cartesian_qw_operator_construction_record(
        mismatched_carried_source,
        operators;
        throw_on_mismatch = false,
    )
    carried_mismatch_diagnostics =
        QWCS.qw_operator_construction_record_diagnostics(carried_mismatch_record)
    @test !carried_mismatch_diagnostics.source_sidecar_agree
    @test !isempty(
        intersect(
            carried_mismatch_diagnostics.mismatch_fields,
            [
                :carried_parent_axis_counts,
                :carried_parent_dimension,
                :carried_representation_final_dimension,
                :carried_provenance_input_kind,
            ],
        ),
    )
    @test_throws ArgumentError QWCS.cartesian_qw_operator_construction_record(
        mismatched_carried_source,
        operators,
    )
end

@testset "Cartesian parent gausslet basis identity" begin
    CP = GaussletBases.CartesianParentGaussletBases
    atomic_axis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 7,
            mapping = IdentityMapping(),
            reference_spacing = 1.0,
        ),
    )
    atomic_parent = CP.CartesianParentGaussletBasis3D(atomic_axis)
    construction_allocated = @allocated CP.CartesianParentGaussletBasis3D(atomic_axis)

    diatomic_basis = bond_aligned_homonuclear_qw_basis(
        bond_length = 1.4;
        core_spacing = 0.5,
        xmax_parallel = 3.0,
        xmax_transverse = 2.0,
    )
    chain_basis = bond_aligned_homonuclear_chain_qw_basis(
        natoms = 3;
        spacing = 1.2,
        core_spacing = 0.5,
        xmax_parallel = 2.5,
        xmax_transverse = 2.0,
    )
    square_basis = axis_aligned_homonuclear_square_lattice_qw_basis(
        n = 2;
        spacing = 1.2,
        core_spacing = 0.5,
        xmax_in_plane = 2.5,
        xmax_transverse = 2.0,
    )
    diatomic_parent = CP.CartesianParentGaussletBasis3D(diatomic_basis)
    chain_parent = CP.CartesianParentGaussletBasis3D(chain_basis)
    square_parent = CP.CartesianParentGaussletBasis3D(square_basis)

    function _check_parent_index_contract(parent)
        dims = CP.parent_axis_counts(parent)
        states = (
            (1, 1, 1),
            (min(2, dims[1]), min(3, dims[2]), min(4, dims[3])),
            dims,
        )
        axes = CP.parent_axes(parent)
        for state in states
            flat = CP.parent_flat_index(parent, state...)
            @test flat == GaussletBases._cartesian_flat_index(state..., dims)
            @test CP.parent_unflat_index(parent, flat) ==
                GaussletBases._cartesian_unflat_index(flat, dims)
            @test CP.parent_center(parent, state) == (
                centers(axes.x)[state[1]],
                centers(axes.y)[state[2]],
                centers(axes.z)[state[3]],
            )
        end
    end

    @test CP.cartesian_parent_gausslet_basis(atomic_parent) === atomic_parent
    @test CP.parent_axis_counts(CP.cartesian_parent_gausslet_basis(atomic_axis)) ==
        CP.parent_axis_counts(atomic_parent)
    @test CP.parent_axis_counts(CP.cartesian_parent_gausslet_basis(diatomic_basis)) ==
        CP.parent_axis_counts(diatomic_parent)
    @test CP.parent_axis_counts(CP.cartesian_parent_gausslet_basis(chain_basis)) ==
        CP.parent_axis_counts(chain_parent)
    @test CP.parent_axis_counts(CP.cartesian_parent_gausslet_basis(square_basis)) ==
        CP.parent_axis_counts(square_parent)

    @test CP.parent_axes(atomic_parent).x === atomic_axis
    @test CP.axis_basis(atomic_parent, :y) === atomic_axis
    @test CP.parent_box(atomic_parent) == (1:7, 1:7, 1:7)
    @test CP.parent_axis_counts(atomic_parent) == (7, 7, 7)
    @test CP.parent_dimension(atomic_parent) == 7^3
    @test atomic_parent.axis_sharing == :shared_xyz
    @test atomic_parent.metadata.basis_family == :mapped_uniform_same_axis
    @test construction_allocated < 50_000
    @test fieldnames(typeof(atomic_parent)) == (:axes, :parent_box, :axis_sharing, :metadata)
    @test !hasproperty(atomic_parent, :gausslet_backend)
    @test !hasproperty(atomic_parent, :backend)
    @test !hasproperty(atomic_parent.metadata, :basis_centers)
    @test !hasproperty(atomic_parent.metadata, :parent_centers)

    @test CP.parent_axes(diatomic_parent).x === diatomic_basis.basis_x
    @test CP.parent_axes(diatomic_parent).z === diatomic_basis.basis_z
    @test diatomic_parent.axis_sharing == :shared_xy
    @test diatomic_parent.metadata.basis_family == :bond_aligned_diatomic
    @test diatomic_parent.metadata.bond_axis == diatomic_basis.bond_axis
    @test diatomic_parent.metadata.nuclei == diatomic_basis.nuclei
    @test diatomic_parent.metadata.nuclear_charges == diatomic_basis.nuclear_charges
    @test CP.parent_box(diatomic_parent) == (
        1:length(diatomic_basis.basis_x),
        1:length(diatomic_basis.basis_y),
        1:length(diatomic_basis.basis_z),
    )

    @test chain_parent.axis_sharing == :shared_xy
    @test chain_parent.metadata.basis_family == :bond_aligned_homonuclear_chain
    @test chain_parent.metadata.chain_axis == chain_basis.chain_axis
    @test chain_parent.metadata.chain_coordinates == chain_basis.chain_coordinates
    @test chain_parent.metadata.nuclei == chain_basis.nuclei

    @test square_parent.axis_sharing == :shared_xy
    @test square_parent.metadata.basis_family == :axis_aligned_homonuclear_square_lattice
    @test square_parent.metadata.lattice_size == square_basis.lattice_size
    @test square_parent.metadata.x_coordinates == square_basis.x_coordinates
    @test square_parent.metadata.y_coordinates == square_basis.y_coordinates

    for parent in (atomic_parent, diatomic_parent, chain_parent, square_parent)
        _check_parent_index_contract(parent)
        @test !hasproperty(parent, :gausslet_backend)
        @test !hasproperty(parent, :backend)
    end

    @test_throws ArgumentError CP.axis_basis(atomic_parent, :q)
    @test_throws ArgumentError CP.parent_flat_index(atomic_parent, 0, 1, 1)
    @test_throws ArgumentError CP.parent_unflat_index(atomic_parent, CP.parent_dimension(atomic_parent) + 1)
    @test_throws ArgumentError CP.cartesian_parent_gausslet_basis((not_a_parent = true,))
end

@testset "Cartesian contracted parent scaffold" begin
    CP = GaussletBases.CartesianParentGaussletBases
    CCP = GaussletBases.CartesianContractedParents
    axis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 3,
            mapping = IdentityMapping(),
            reference_spacing = 1.0,
        ),
    )
    parent = CP.cartesian_parent_gausslet_basis(axis)
    parent_dim = CP.parent_dimension(parent)

    coefficients = zeros(Float64, parent_dim, 4)
    coefficients[1, 1] = 1.0
    coefficients[2, 2] = 1.0
    coefficients[2, 3] = 2.0
    coefficients[5, 4] = 1.0
    coefficients[6, 4] = -1.0
    unit_a = CCP.CartesianContractionUnit3D(
        :cube_a,
        [1, 2, 3],
        1:2;
        metadata = (shape = :cube,),
    )
    unit_b = CCP.CartesianContractionUnit3D(
        :cube_b,
        [2, 5, 6],
        3:4;
        metadata = (shape = :cube,),
    )
    contracted = CCP.CartesianContractedParent3D(
        parent,
        coefficients;
        units = [unit_a, unit_b],
        metadata = (source = :synthetic,),
    )
    audit = CCP.contracted_parent_structural_audit(contracted)

    @test CCP.contracted_parent_basis(contracted) === parent
    @test CCP.contracted_parent_coefficients(contracted) == coefficients
    @test CCP.contracted_parent_parent_dimension(contracted) == parent_dim
    @test CCP.contracted_parent_dimension(contracted) == 4
    @test CCP.contracted_parent_metadata(contracted).source == :synthetic
    @test CCP.contracted_parent_units(contracted) == [unit_a, unit_b]
    @test CCP.contracted_parent_unit_column_ranges(contracted) == [1:2, 3:4]
    @test CCP.contracted_parent_unit_support_indices(contracted) == [[1, 2, 3], [2, 5, 6]]
    @test CCP.contracted_parent_support_indices(contracted) == [1, 2, 3, 2, 5, 6]
    @test CCP.contraction_unit_role(unit_a) == :cube_a
    @test CCP.contraction_unit_support_indices(unit_a) == [1, 2, 3]
    @test CCP.contraction_unit_column_range(unit_a) == 1:2
    @test CCP.contraction_unit_metadata(unit_a).shape == :cube

    @test coefficients[:, 3] == 2.0 .* coefficients[:, 2]
    @test audit.parent_dimension == parent_dim
    @test audit.contracted_dimension == 4
    @test audit.unit_count == 2
    @test audit.support_entry_count == 6
    @test audit.unique_support_count == 5
    @test audit.duplicate_support_count == 1
    @test audit.missing_support_count == parent_dim - 5
    @test audit.outside_support_count == 0
    @test !audit.support_complete
    @test audit.column_entry_count == 4
    @test audit.unique_column_count == 4
    @test audit.duplicate_column_count == 0
    @test audit.missing_column_count == 0
    @test audit.outside_column_count == 0
    @test audit.column_ranges_cover_contract
    @test audit.structural_ok
    @test !hasproperty(contracted, :gausslet_backend)
    @test !hasproperty(contracted, :backend)
    @test !hasproperty(contracted, :overlap)
    @test !hasproperty(contracted, :interaction_matrix)

    sparse_coefficients = sparse(coefficients)
    sparse_contracted = CCP.CartesianContractedParent3D(
        parent,
        sparse_coefficients;
        units = [unit_a, unit_b],
    )
    @test CCP.contracted_parent_coefficients(sparse_contracted) isa SparseMatrixCSC{Float64,Int}
    @test CCP.contracted_parent_coefficients(sparse_contracted) == sparse_coefficients

    outside_unit = CCP.CartesianContractionUnit3D(:outside, [1, parent_dim + 1], 1:1)
    outside = CCP.CartesianContractedParent3D(
        parent,
        coefficients[:, 1:1];
        units = [outside_unit],
    )
    outside_audit = CCP.contracted_parent_structural_audit(outside)
    @test outside_audit.outside_support_count == 1
    @test !outside_audit.support_complete
    @test !outside_audit.structural_ok

    overlapping_columns = CCP.CartesianContractedParent3D(
        parent,
        coefficients;
        units = [
            CCP.CartesianContractionUnit3D(:left, [1], 1:2),
            CCP.CartesianContractionUnit3D(:right, [2], 2:4),
        ],
    )
    overlapping_column_audit = CCP.contracted_parent_structural_audit(overlapping_columns)
    @test overlapping_column_audit.duplicate_column_count == 1
    @test !overlapping_column_audit.column_ranges_cover_contract
    @test !overlapping_column_audit.structural_ok

    missing_columns = CCP.CartesianContractedParent3D(
        parent,
        coefficients;
        units = [
            CCP.CartesianContractionUnit3D(:first, [1], 1:1),
            CCP.CartesianContractionUnit3D(:last, [2], 3:4),
        ],
    )
    missing_column_audit = CCP.contracted_parent_structural_audit(missing_columns)
    @test missing_column_audit.missing_column_count == 1
    @test !missing_column_audit.column_ranges_cover_contract
    @test !missing_column_audit.structural_ok

    @test_throws ArgumentError CCP.CartesianContractionUnit3D(:empty, [1], 1:0)
    @test_throws ArgumentError CCP.CartesianContractedParent3D(
        parent,
        coefficients;
        units = [CCP.CartesianContractionUnit3D(:bad_columns, [1], 4:5)],
    )
    @test_throws DimensionMismatch CCP.CartesianContractedParent3D(
        parent,
        coefficients[1:(end - 1), :],
    )
end

@testset "Cartesian contracted parent metric packet" begin
    CP = GaussletBases.CartesianParentGaussletBases
    CCP = GaussletBases.CartesianContractedParents
    CCPM = GaussletBases.CartesianContractedParentMetrics
    axis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 3,
            mapping = IdentityMapping(),
            reference_spacing = 1.0,
        ),
    )
    parent = CP.cartesian_parent_gausslet_basis(axis)
    parent_dim = CP.parent_dimension(parent)
    coefficients = zeros(Float64, parent_dim, 3)
    coefficients[1, 1] = 1.0
    coefficients[2, 1] = 0.25
    coefficients[5, 2] = 1.0
    coefficients[14, 2] = -0.5
    coefficients[parent_dim, 3] = 1.0
    units = [
        CCP.CartesianContractionUnit3D(:left, [1, 2], 1:1),
        CCP.CartesianContractionUnit3D(:middle, [5, 14], 2:2),
        CCP.CartesianContractionUnit3D(:right, [parent_dim], 3:3),
    ]
    contracted = CCP.CartesianContractedParent3D(parent, coefficients; units)
    packet = CCPM.cartesian_contracted_parent_metric_packet(contracted)
    reference = CCPM.cartesian_contracted_parent_metric_packet_dense_reference(contracted)

    @test packet isa CCPM.CartesianContractedParentMetricPacket3D
    @test CCPM.contracted_parent_metric_packet_parent(packet) === contracted
    @test packet.diagnostics.construction_path == :support_local_product
    @test packet.diagnostics.dense_parent_matrix_used == false
    @test reference.diagnostics.construction_path == :dense_reference_oracle
    @test reference.diagnostics.dense_parent_matrix_used == true
    @test size(packet.overlap) == (3, 3)
    @test size(packet.centers) == (3, 3)
    @test length(packet.weights) == 3
    @test packet.overlap ≈ reference.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test packet.position_x ≈ reference.position_x atol = 1.0e-12 rtol = 1.0e-12
    @test packet.position_y ≈ reference.position_y atol = 1.0e-12 rtol = 1.0e-12
    @test packet.position_z ≈ reference.position_z atol = 1.0e-12 rtol = 1.0e-12
    @test packet.weights ≈ reference.weights atol = 1.0e-12 rtol = 1.0e-12
    @test packet.first_moments ≈ reference.first_moments atol = 1.0e-12 rtol = 1.0e-12
    @test packet.centers ≈ reference.centers atol = 1.0e-12 rtol = 1.0e-12
    @test isfinite(packet.diagnostics.overlap_symmetry_error)
    @test isfinite(packet.diagnostics.overlap_identity_error)

    sparse_coefficients = sparse([1, 7, 13, 27], [1, 1, 2, 2], [0.5, -0.25, 1.0, 0.125], parent_dim, 2)
    sparse_contracted = CCP.CartesianContractedParent3D(parent, sparse_coefficients)
    sparse_packet = CCPM.cartesian_contracted_parent_metric_packet(sparse_contracted)
    sparse_reference = CCPM.cartesian_contracted_parent_metric_packet_dense_reference(sparse_contracted)
    @test CCP.contracted_parent_coefficients(sparse_contracted) isa SparseMatrixCSC{Float64,Int}
    @test sparse_packet.diagnostics.dense_parent_matrix_used == false
    @test sparse_packet.diagnostics.coefficient_storage == :SparseMatrixCSC
    @test sparse_packet.overlap ≈ sparse_reference.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test sparse_packet.weights ≈ sparse_reference.weights atol = 1.0e-12 rtol = 1.0e-12

    larger_axis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 5,
            mapping = IdentityMapping(),
            reference_spacing = 1.0,
        ),
    )
    larger_parent = CP.cartesian_parent_gausslet_basis(larger_axis)
    larger_coefficients = sparse([1, 125], [1, 2], [1.0, 1.0], 125, 2)
    larger_contracted = CCP.CartesianContractedParent3D(larger_parent, larger_coefficients)
    larger_packet = CCPM.cartesian_contracted_parent_metric_packet(larger_contracted)
    @test larger_packet.diagnostics.dense_parent_matrix_used == false
    @test_throws ArgumentError CCPM.cartesian_contracted_parent_metric_packet_dense_reference(
        larger_contracted;
        max_parent_dimension = 64,
    )

    distorted_axis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 3,
            mapping = AsinhMapping(a = 0.25, s = 0.5, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    distorted_parent = CP.cartesian_parent_gausslet_basis(distorted_axis)
    distorted_coefficients = sparse([1, 27], [1, 2], [1.0, 1.0], 27, 2)
    distorted_contracted = CCP.CartesianContractedParent3D(
        distorted_parent,
        distorted_coefficients,
    )
    explicit_axis_metric = (
        overlap = Matrix{Float64}(I, 3, 3),
        position = Matrix(Diagonal(centers(distorted_axis))),
        weights = ones(Float64, 3),
        centers = Float64.(centers(distorted_axis)),
        source = :test_explicit_no_quadrature,
    )
    explicit_packet = CCPM.cartesian_contracted_parent_metric_packet(
        distorted_contracted;
        axis_metrics = (
            x = explicit_axis_metric,
            y = explicit_axis_metric,
            z = explicit_axis_metric,
        ),
    )
    @test_throws ArgumentError CCPM.cartesian_contracted_parent_metric_packet(distorted_contracted)
    @test explicit_packet.diagnostics.dense_parent_matrix_used == false
    @test explicit_packet.diagnostics.axis_metric_sources.x == :test_explicit_no_quadrature
    @test_throws ArgumentError CCPM.cartesian_contracted_parent_metric_packet(
        sparse_contracted;
        construction_path = :product_staged_metric_contraction,
    )
end

@testset "Cartesian basis representation for nested fixed blocks" begin
    CP = GaussletBases.CartesianParentGaussletBases
    CCP = GaussletBases.CartesianContractedParents
    CCS = GaussletBases.CartesianCarriedSpaces
    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 13,
            mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    fixed_block = one_center_atomic_full_parent_fixed_block(
        basis;
        expansion,
        nside = 5,
    )
    fixed_parent = CP.cartesian_parent_gausslet_basis(fixed_block)
    fixed_contracted_parent = CCP.cartesian_contracted_parent(fixed_block)
    fixed_contracted_audit = CCP.contracted_parent_structural_audit(fixed_contracted_parent)
    representation = basis_representation(fixed_block)
    metadata = basis_metadata(representation)

    @test representation isa CartesianBasisRepresentation3D
    @test metadata.basis_kind == :nested_fixed_block
    @test metadata.parent_kind == :cartesian_product_basis
    @test metadata.parent_axis_counts == (13, 13, 13)
    @test metadata.parent_axis_counts == CP.parent_axis_counts(fixed_parent)
    @test metadata.parent_dimension == 13^3
    @test metadata.parent_dimension == CP.parent_dimension(fixed_parent)
    @test metadata.final_dimension == size(fixed_block.coefficient_matrix, 2)
    @test metadata.working_box == (1:13, 1:13, 1:13)
    @test metadata.route_metadata.shell_kind == :shell_sequence
    @test metadata.route_metadata.working_box_profile == :full_parent
    @test metadata.route_metadata.nside == 5
    @test metadata.route_metadata.support_count == length(fixed_block.support_indices)
    @test size(representation.coefficient_matrix) == size(fixed_block.coefficient_matrix)
    @test representation.support_indices == fixed_block.support_indices
    @test length(representation.support_states) == length(fixed_block.support_indices)
    @test size(metadata.basis_centers) == size(fixed_block.fixed_centers)
    @test CP.parent_center(fixed_parent, (1, 1, 1)) == (
        centers(basis)[1],
        centers(basis)[1],
        centers(basis)[1],
    )
    @test !hasproperty(fixed_parent, :gausslet_backend)
    @test !hasproperty(fixed_parent, :backend)
    @test CCP.contracted_parent_basis(fixed_contracted_parent).parent_box ==
        fixed_parent.parent_box
    @test CCP.contracted_parent_parent_dimension(fixed_contracted_parent) ==
        metadata.parent_dimension
    @test CCP.contracted_parent_dimension(fixed_contracted_parent) ==
        size(fixed_block.coefficient_matrix, 2)
    @test CCP.contracted_parent_coefficients(fixed_contracted_parent) isa
        SparseMatrixCSC{Float64,Int}
    @test CCP.contracted_parent_coefficients(fixed_contracted_parent) ==
        Matrix{Float64}(fixed_block.coefficient_matrix)
    @test only(CCP.contracted_parent_units(fixed_contracted_parent)).role ==
        :nested_fixed_block
    @test only(CCP.contracted_parent_unit_column_ranges(fixed_contracted_parent)) ==
        1:size(fixed_block.coefficient_matrix, 2)
    @test fixed_contracted_audit.outside_support_count == 0
    @test fixed_contracted_audit.column_ranges_cover_contract
    @test fixed_contracted_audit.structural_ok
    @test !hasproperty(fixed_contracted_parent, :gausslet_backend)
    @test !hasproperty(fixed_contracted_parent, :backend)
    @test !hasproperty(fixed_contracted_parent, :interaction_matrix)
    @test !isnothing(fixed_block.factorized_cartesian_parent_basis[])
    @test hasproperty(representation.parent_data, :factorized_cartesian_parent_basis)
    @test representation.parent_data.factorized_cartesian_parent_basis ===
          fixed_block.factorized_cartesian_parent_basis[]
    fixed_carried = CCS.cartesian_carried_space(fixed_block)
    @test CCS.carried_space_parent(fixed_carried) isa CP.CartesianParentGaussletBasis3D
    @test CCS.carried_space_contracted_parent(fixed_carried) isa CCP.CartesianContractedParent3D
    @test CCS.carried_space_representation(fixed_carried) isa CartesianBasisRepresentation3D
    @test CCS.carried_space_diagnostics(fixed_carried).parent_dimension ==
        metadata.parent_dimension
    @test CCS.carried_space_diagnostics(fixed_carried).contracted_dimension ==
        metadata.final_dimension
    @test CCS.carried_space_diagnostics(fixed_carried).contracted_dimension_matches_representation
    @test CCS.carried_space_diagnostics(fixed_carried).contracted_parent_dimension_matches_parent
    @test CCS.carried_space_provenance(fixed_carried).input_kind == :nested_fixed_block

    square_basis, _source, square_fixed_block, _diagnostics =
        _axis_aligned_homonuclear_square_lattice_nested_fixture()
    square_representation = basis_representation(square_fixed_block)
    square_metadata = basis_metadata(square_representation)
    @test square_metadata.basis_kind == :nested_fixed_block
    @test square_metadata.parent_axis_counts == (
        length(square_basis.basis_x),
        length(square_basis.basis_y),
        length(square_basis.basis_z),
    )
    @test square_metadata.parent_dimension == prod(square_metadata.parent_axis_counts)
    @test square_metadata.final_dimension == size(square_fixed_block.coefficient_matrix, 2)
    @test size(square_representation.coefficient_matrix) == size(square_fixed_block.coefficient_matrix)
    @test square_metadata.working_box == square_fixed_block.shell.working_box
    @test square_metadata.route_metadata.support_count == length(square_fixed_block.support_indices)
    @test !isnothing(square_fixed_block.factorized_cartesian_parent_basis[])
    @test hasproperty(square_representation.parent_data, :factorized_cartesian_parent_basis)
    @test square_representation.parent_data.factorized_cartesian_parent_basis ===
          square_fixed_block.factorized_cartesian_parent_basis[]
    square_carried = CCS.cartesian_carried_space(square_fixed_block)
    @test CCS.carried_space_parent(square_carried) isa CP.CartesianParentGaussletBasis3D
    @test CCS.carried_space_contracted_parent(square_carried) isa CCP.CartesianContractedParent3D
    @test CCS.carried_space_diagnostics(square_carried).parent_axis_counts ==
        square_metadata.parent_axis_counts
    @test CCS.carried_space_diagnostics(square_carried).contracted_dimension ==
        square_metadata.final_dimension
    @test CCS.carried_space_diagnostics(square_carried).contracted_dimension_matches_representation
    @test CCS.carried_space_provenance(square_carried).input_kind == :nested_fixed_block
end

function _with_sparse_nested_coefficients(fixed_block::GaussletBases._NestedFixedBlock3D)
    return GaussletBases._NestedFixedBlock3D(
        fixed_block.parent_basis,
        fixed_block.shell,
        fixed_block.gausslet_backend,
        sparse(fixed_block.coefficient_matrix),
        fixed_block.support_indices,
        fixed_block.overlap,
        fixed_block.kinetic,
        fixed_block.position_x,
        fixed_block.position_y,
        fixed_block.position_z,
        fixed_block.x2_x,
        fixed_block.x2_y,
        fixed_block.x2_z,
        fixed_block.weights,
        fixed_block.gaussian_sum,
        fixed_block.pair_sum,
        fixed_block.fixed_centers,
        GaussletBases._nested_factorized_basis_cache(
            fixed_block.factorized_cartesian_parent_basis[],
        ),
        GaussletBases._nested_staged_by_center_sidecar_cache(
            fixed_block.staged_by_center_sidecar[],
        ),
    )
end

@testset "Nested coefficient maps support sparse storage" begin
    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 13,
            mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    direct_fixed_block = one_center_atomic_full_parent_fixed_block(
        basis;
        expansion,
        nside = 5,
    )
    sparse_fixed_block = _with_sparse_nested_coefficients(direct_fixed_block)

    direct_representation = basis_representation(direct_fixed_block)
    sparse_representation = basis_representation(sparse_fixed_block)
    support_coefficients = GaussletBases._nested_support_coefficient_slice(
        direct_fixed_block.shell.coefficient_matrix,
        direct_fixed_block.shell.support_indices,
    )
    support_workspace, contraction_scratch = GaussletBases._nested_support_reference_workspaces(
        support_coefficients,
        length(direct_fixed_block.shell.support_indices),
        size(direct_fixed_block.shell.coefficient_matrix, 2),
    )

    @test direct_fixed_block.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test direct_representation.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test sparse_fixed_block.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test sparse_representation.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test support_coefficients isa SparseMatrixCSC{Float64,Int}
    @test size(support_workspace) == (0, 0)
    @test size(contraction_scratch) == (0, 0)
    @test size(support_coefficients) == (
        length(direct_fixed_block.shell.support_indices),
        size(direct_fixed_block.shell.coefficient_matrix, 2),
    )
    @test Matrix(sparse_representation.coefficient_matrix) ≈
        Matrix(direct_representation.coefficient_matrix) atol = 1.0e-12 rtol = 1.0e-12
    @test cross_overlap(sparse_representation, sparse_representation) ≈
        cross_overlap(direct_representation, direct_representation) atol = 1.0e-10 rtol = 1.0e-10

    mktemp() do sparse_path, sparse_io
        close(sparse_io)
        sparse_matrix = sparse(direct_fixed_block.coefficient_matrix)
        jldopen(sparse_path, "w") do file
            file["matrix"] = sparse_matrix
        end
        restored = jldopen(sparse_path, "r") do file
            file["matrix"]
        end
        @test restored isa SparseMatrixCSC{Float64,Int}
        @test restored == sparse_matrix
    end
end

function _atomic_hybrid_cartesian_representation_fixture()
    return _cached_fixture(:atomic_hybrid_cartesian_representation_fixture, () -> begin
        basis = build_basis(
            MappedUniformBasisSpec(
                :G10;
                count = 13,
                mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
                reference_spacing = 1.0,
            ),
        )
        expansion = coulomb_gaussian_expansion(doacc = false)
        fixed_full = one_center_atomic_full_parent_fixed_block(
            basis;
            expansion,
            nside = 5,
        )
        fixed_legacy = one_center_atomic_legacy_profile_fixed_block(
            basis;
            expansion,
            working_box = 2:12,
            nside = 5,
        )
        supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 1)
        full_ops = ordinary_cartesian_qiu_white_operators(
            fixed_full,
            supplement;
            expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
            residual_keep_policy = :near_null_only,
        )
        legacy_ops = ordinary_cartesian_qiu_white_operators(
            fixed_legacy,
            supplement;
            expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
            residual_keep_policy = :near_null_only,
        )
        (
            basis = basis,
            expansion = expansion,
            fixed_full = fixed_full,
            fixed_legacy = fixed_legacy,
            fixed_full_rep = basis_representation(fixed_full),
            fixed_legacy_rep = basis_representation(fixed_legacy),
            supplement = supplement,
            full_ops = full_ops,
            legacy_ops = legacy_ops,
            full_rep = basis_representation(full_ops),
            legacy_rep = basis_representation(legacy_ops),
        )
    end)
end

function _metric_normalize_orbital(
    coefficients::AbstractVector,
    overlap::AbstractMatrix{<:Real},
)
    orbital = Float64[Float64(real(value)) for value in coefficients]
    norm2 = Float64(real(dot(orbital, overlap * orbital)))
    norm2 > 0.0 || throw(ArgumentError("orbital must have nonzero target-metric norm"))
    return orbital ./ sqrt(norm2)
end

function _metric_orbital_overlap(
    left::AbstractVector,
    right::AbstractVector,
    overlap::AbstractMatrix{<:Real},
)
    normalized_left = _metric_normalize_orbital(left, overlap)
    normalized_right = _metric_normalize_orbital(right, overlap)
    return Float64(real(dot(normalized_left, overlap * normalized_right)))
end

function _ordinary_cartesian_hybrid_orbital_observables(
    operators::OrdinaryCartesianOperators3D,
    orbital::AbstractVector;
    overlap_tol::Real = 1.0e-7,
)
    overlap = Matrix{Float64}(operators.overlap)
    normalized = _metric_normalize_orbital(orbital, overlap)
    one_body = Float64(real(dot(normalized, operators.one_body_hamiltonian * normalized)))
    vee = GaussletBases.ordinary_cartesian_vee_expectation(
        operators,
        normalized;
        overlap_tol = overlap_tol,
    )
    return (
        orbital = normalized,
        metric_norm_error = abs(Float64(real(dot(normalized, overlap * normalized))) - 1.0),
        one_body = one_body,
        vee = vee,
        total = 2.0 * one_body + vee,
    )
end

function _atomic_direct_product_he_extent_change_contract_fixture(;
    source_count::Int = 3,
    target_count::Int = 5,
)
    key = Symbol(:atomic_direct_product_he_extent_change_contract, source_count, target_count)
    return _cached_fixture(key, () -> begin
        mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0)
        source_basis = build_basis(
            MappedUniformBasisSpec(
                :G10;
                count = source_count,
                mapping,
                reference_spacing = 1.0,
            ),
        )
        target_basis = build_basis(
            MappedUniformBasisSpec(
                :G10;
                count = target_count,
                mapping,
                reference_spacing = 1.0,
            ),
        )

        source_rep = basis_representation(source_basis)
        target_rep = basis_representation(target_basis)
        offset = (target_count - source_count) ÷ 2
        shared_slice = (offset + 1):(offset + source_count)

        return (
            source_count = source_count,
            target_count = target_count,
            shared_slice = shared_slice,
            source_rep = source_rep,
            target_rep = target_rep,
            centers_subset =
                source_rep.metadata.center_data == target_rep.metadata.center_data[shared_slice],
            weights_subset =
                source_rep.metadata.integral_weight_data ==
                target_rep.metadata.integral_weight_data[shared_slice],
            coefficient_core_match =
                source_rep.coefficient_matrix ==
                target_rep.coefficient_matrix[shared_slice, shared_slice],
        )
    end)
end

function _atomic_hybrid_he_same_parent_stress_fixture(;
    parent_count::Int = 7,
    source_working_box::UnitRange{Int} = 2:6,
    supplement_lmax::Int = 1,
)
    key = Symbol(
        :atomic_hybrid_he_same_parent_stress_fixture,
        parent_count,
        first(source_working_box),
        last(source_working_box),
        supplement_lmax,
    )
    return _cached_fixture(key, () -> begin
        expansion = coulomb_gaussian_expansion(doacc = false)
        mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0)
        supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = supplement_lmax)

        parent_basis = build_basis(
            MappedUniformBasisSpec(
                :G10;
                count = parent_count,
                mapping,
                reference_spacing = 1.0,
            ),
        )

        source_fixed = one_center_atomic_legacy_profile_fixed_block(
            parent_basis;
            expansion,
            working_box = source_working_box,
            nside = 5,
        )
        target_fixed = one_center_atomic_full_parent_fixed_block(
            parent_basis;
            expansion,
            nside = 5,
        )

        source_ops = ordinary_cartesian_qiu_white_operators(
            source_fixed,
            supplement;
            expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
            residual_keep_policy = :near_null_only,
        )
        target_ops = ordinary_cartesian_qiu_white_operators(
            target_fixed,
            supplement;
            expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
            residual_keep_policy = :near_null_only,
        )

        source_rep = basis_representation(source_ops)
        target_rep = basis_representation(target_ops)
        source_check = GaussletBases.ordinary_cartesian_1s2_check(
            source_ops;
            overlap_tol = 1.0e-7,
        )
        target_check = GaussletBases.ordinary_cartesian_1s2_check(
            target_ops;
            overlap_tol = 1.0e-7,
        )
        source_observables = _ordinary_cartesian_hybrid_orbital_observables(
            source_ops,
            source_check.orbital;
            overlap_tol = 1.0e-7,
        )
        target_observables = _ordinary_cartesian_hybrid_orbital_observables(
            target_ops,
            target_check.orbital;
            overlap_tol = 1.0e-7,
        )

        transfer = transfer_orbitals(source_observables.orbital, source_rep, target_rep)
        transferred_observables = _ordinary_cartesian_hybrid_orbital_observables(
            target_ops,
            transfer.coefficients;
            overlap_tol = 1.0e-7,
        )
        target_overlap = Matrix{Float64}(target_ops.overlap)
        overlap_with_target = _metric_orbital_overlap(
            transferred_observables.orbital,
            target_observables.orbital,
            target_overlap,
        )
        sign = overlap_with_target < 0.0 ? -1.0 : 1.0
        aligned_transferred_observables = _ordinary_cartesian_hybrid_orbital_observables(
            target_ops,
            sign .* transferred_observables.orbital;
            overlap_tol = 1.0e-7,
        )
        aligned_overlap_to_target = abs(
            _metric_orbital_overlap(
                aligned_transferred_observables.orbital,
                target_observables.orbital,
                target_overlap,
            ),
        )

        return (
            parent_basis = parent_basis,
            source_fixed = source_fixed,
            target_fixed = target_fixed,
            supplement = supplement,
            source_working_box = source_working_box,
            target_working_box = target_fixed.shell.working_box,
            source_ops = source_ops,
            target_ops = target_ops,
            source_rep = source_rep,
            target_rep = target_rep,
            source_check = source_check,
            target_check = target_check,
            source_observables = source_observables,
            target_observables = target_observables,
            transfer = transfer,
            transferred_observables = transferred_observables,
            aligned_transferred_observables = aligned_transferred_observables,
            aligned_overlap_to_target = aligned_overlap_to_target,
        )
    end)
end

@testset "Cartesian basis representation for atomic QW residual bases" begin
    fixture = _atomic_hybrid_cartesian_representation_fixture()
    operators = fixture.full_ops
    representation = fixture.full_rep
    metadata = basis_metadata(representation)
    supplement_representation = representation.parent_data.supplement_representation

    @test representation isa CartesianBasisRepresentation3D
    @test metadata.basis_kind == :hybrid_residual
    @test metadata.parent_kind == :cartesian_plus_supplement_raw
    @test metadata.final_dimension == length(operators.orbital_data)
    @test metadata.final_dimension == size(operators.raw_to_final, 2)
    @test metadata.parent_dimension == size(operators.raw_to_final, 1)
    @test metadata.route_metadata.gausslet_count == operators.gausslet_count
    @test metadata.route_metadata.residual_count == operators.residual_count
    @test metadata.route_metadata.supplement_kind == :atomic_cartesian_shell
    @test metadata.route_metadata.supplement_lmax == fixture.supplement.lmax
    @test size(representation.coefficient_matrix) == size(operators.raw_to_final)
    @test length(representation.parent_labels) == size(operators.raw_to_final, 1)
    @test size(representation.parent_centers, 1) == size(operators.raw_to_final, 1)
    @test hasproperty(representation.parent_data, :cartesian_parent_representation)
    @test representation.parent_data.cartesian_parent_representation.metadata.basis_kind ==
        :nested_fixed_block
    @test representation.parent_data.cartesian_parent_representation.metadata.final_dimension ==
        operators.gausslet_count
    @test hasproperty(representation.parent_data, :supplement_representation)
    @test hasproperty(representation.parent_data, :factorized_cartesian_parent_basis)
    @test hasproperty(representation.parent_data, :cartesian_supplement_axis_tables)
    @test supplement_representation isa CartesianGaussianShellSupplementRepresentation3D
    @test supplement_representation.supplement_kind == :atomic_cartesian_shell
    @test length(supplement_representation.orbitals) ==
        size(operators.raw_to_final, 1) - operators.gausslet_count
    @test size(representation.parent_data.cartesian_supplement_axis_tables.x, 2) ==
        length(supplement_representation.orbitals)
    @test size(representation.parent_data.cartesian_supplement_axis_tables.y, 2) ==
        length(supplement_representation.orbitals)
    @test size(representation.parent_data.cartesian_supplement_axis_tables.z, 2) ==
        length(supplement_representation.orbitals)
    @test any(
        orbital -> sum(orbital.angular_powers) > 0,
        supplement_representation.orbitals,
    )
end

@testset "Cartesian basis representation cross overlap" begin
    diatomic_basis14, diatomic_ops14, _check14 = _bond_aligned_diatomic_qw_fixture(; bond_length = 1.4)
    diatomic_basis20, _diatomic_ops20, _check20 = _bond_aligned_diatomic_qw_fixture(; bond_length = 2.0)
    diatomic_rep14 = basis_representation(diatomic_basis14)
    diatomic_rep20 = basis_representation(diatomic_basis20)
    direct_self = cross_overlap(diatomic_rep14, diatomic_rep14)
    direct_cross = cross_overlap(diatomic_rep14, diatomic_rep20)
    direct_cross_reverse = cross_overlap(diatomic_rep20, diatomic_rep14)

    @test size(direct_self) == size(diatomic_ops14.overlap)
    @test direct_self ≈ diatomic_ops14.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test size(direct_cross) == (
        diatomic_rep14.metadata.final_dimension,
        diatomic_rep20.metadata.final_dimension,
    )
    @test direct_cross ≈ transpose(direct_cross_reverse) atol = 1.0e-10 rtol = 1.0e-10

    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 13,
            mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    fixed_full = one_center_atomic_full_parent_fixed_block(
        basis;
        expansion = expansion,
        nside = 5,
    )
    fixed_legacy = one_center_atomic_legacy_profile_fixed_block(
        basis;
        expansion = expansion,
        working_box = 2:12,
        nside = 5,
    )
    fixed_full_rep = basis_representation(fixed_full)
    fixed_legacy_rep = basis_representation(fixed_legacy)
    S1d = basis_representation(basis).basis_matrices.overlap
    Sparent = kron(S1d, kron(S1d, S1d))

    fixed_self = cross_overlap(fixed_full_rep, fixed_full_rep)
    fixed_cross = cross_overlap(fixed_full_rep, fixed_legacy_rep)
    fixed_cross_reverse = cross_overlap(fixed_legacy_rep, fixed_full_rep)
    fixed_cross_expected =
        transpose(fixed_full.coefficient_matrix) * Sparent * fixed_legacy.coefficient_matrix

    @test size(fixed_self) == size(fixed_full.overlap)
    @test fixed_self ≈ fixed_full.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test size(fixed_cross) == (
        size(fixed_full.coefficient_matrix, 2),
        size(fixed_legacy.coefficient_matrix, 2),
    )
    @test fixed_cross ≈ fixed_cross_expected atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_cross ≈ transpose(fixed_cross_reverse) atol = 1.0e-10 rtol = 1.0e-10

    square_basis, _square_source, square_fixed_block, _square_diagnostics =
        _axis_aligned_homonuclear_square_lattice_nested_fixture()
    square_basis_rep = basis_representation(square_basis)
    square_fixed_rep = basis_representation(square_fixed_block)
    square_parent_x = basis_representation(square_basis.basis_x).basis_matrices.overlap
    square_parent_y = basis_representation(square_basis.basis_y).basis_matrices.overlap
    square_parent_z = basis_representation(square_basis.basis_z).basis_matrices.overlap
    square_parent_overlap = kron(square_parent_x, kron(square_parent_y, square_parent_z))
    square_cross = cross_overlap(square_basis_rep, square_fixed_rep)
    square_cross_expected = square_parent_overlap * square_fixed_block.coefficient_matrix

    @test size(square_cross) == (
        square_basis_rep.metadata.final_dimension,
        square_fixed_rep.metadata.final_dimension,
    )
    @test square_cross ≈ square_cross_expected atol = 1.0e-10 rtol = 1.0e-10

    hybrid_fixture = _atomic_hybrid_cartesian_representation_fixture()
    hybrid_self = cross_overlap(hybrid_fixture.full_rep, hybrid_fixture.full_rep)
    hybrid_cross = cross_overlap(hybrid_fixture.full_rep, hybrid_fixture.legacy_rep)
    hybrid_cross_reverse = cross_overlap(hybrid_fixture.legacy_rep, hybrid_fixture.full_rep)
    hybrid_parent_cross = cross_overlap(hybrid_fixture.fixed_full_rep, hybrid_fixture.full_rep)
    hybrid_parent_cross_reverse = cross_overlap(hybrid_fixture.full_rep, hybrid_fixture.fixed_full_rep)
    hybrid_self_dense = GaussletBases._cartesian_mixed_raw_cross_overlap_dense_reference(
        hybrid_fixture.full_rep,
        hybrid_fixture.full_rep,
    )

    @test size(hybrid_self) == size(hybrid_fixture.full_ops.overlap)
    @test hybrid_self ≈ hybrid_self_dense atol = 1.0e-10 rtol = 1.0e-10
    @test size(hybrid_cross) == (
        hybrid_fixture.full_rep.metadata.final_dimension,
        hybrid_fixture.legacy_rep.metadata.final_dimension,
    )
    @test hybrid_cross ≈ transpose(hybrid_cross_reverse) atol = 1.0e-10 rtol = 1.0e-10
    @test size(hybrid_parent_cross) == (
        hybrid_fixture.fixed_full_rep.metadata.final_dimension,
        hybrid_fixture.full_rep.metadata.final_dimension,
    )
    @test hybrid_parent_cross ≈ transpose(hybrid_parent_cross_reverse) atol = 1.0e-10 rtol = 1.0e-10
end

@testset "Cartesian basis projector and orbital transfer" begin
    square_basis, _square_source, square_fixed_block, _square_diagnostics =
        _axis_aligned_homonuclear_square_lattice_nested_fixture()
    square_basis_rep = basis_representation(square_basis)
    square_fixed_rep = basis_representation(square_fixed_block)

    direct_self_projector = basis_projector(square_basis_rep, square_basis_rep)
    direct_self_coefficients =
        reshape(sin.(Float64.(1:(2 * square_basis_rep.metadata.final_dimension))), :, 2)
    direct_self_transfer =
        transfer_orbitals(direct_self_coefficients, square_basis_rep, square_basis_rep)

    @test direct_self_projector.matrix ≈ Matrix{Float64}(
        LinearAlgebra.I,
        square_basis_rep.metadata.final_dimension,
        square_basis_rep.metadata.final_dimension,
    ) atol = 1.0e-10 rtol = 1.0e-10
    @test direct_self_transfer.coefficients ≈ direct_self_coefficients atol = 1.0e-10 rtol = 1.0e-10
    @test direct_self_transfer.diagnostics.transfer_path == :same_parent_cross_overlap_transfer
    @test direct_self_transfer.diagnostics.projector_residual_inf < 1.0e-10
    @test direct_self_transfer.diagnostics.transferred_residual_inf < 1.0e-10
    @test direct_self_transfer.diagnostics.source_metric_trace ≈ direct_self_transfer.diagnostics.target_metric_trace atol =
          1.0e-10 rtol = 1.0e-10

    nested_coefficients = cos.(Float64.(1:square_fixed_rep.metadata.final_dimension))
    nested_to_direct = transfer_orbitals(nested_coefficients, square_fixed_rep, square_basis_rep)
    nested_embedded = square_fixed_block.coefficient_matrix * nested_coefficients
    direct_back_to_nested =
        transfer_orbitals(nested_to_direct.coefficients, square_basis_rep, square_fixed_rep)

    @test nested_to_direct.coefficients ≈ nested_embedded atol = 1.0e-10 rtol = 1.0e-10
    @test direct_back_to_nested.coefficients ≈ nested_coefficients atol = 1.0e-10 rtol = 1.0e-10
    @test direct_back_to_nested.diagnostics.transferred_residual_inf < 1.0e-10

    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 13,
            mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    fixed_full = one_center_atomic_full_parent_fixed_block(
        basis;
        expansion = expansion,
        nside = 5,
    )
    fixed_legacy = one_center_atomic_legacy_profile_fixed_block(
        basis;
        expansion = expansion,
        working_box = 2:12,
        nside = 5,
    )
    fixed_full_rep = basis_representation(fixed_full)
    fixed_legacy_rep = basis_representation(fixed_legacy)

    full_to_legacy = basis_projector(fixed_full_rep, fixed_legacy_rep)
    legacy_to_full = basis_projector(fixed_legacy_rep, fixed_full_rep)
    full_coefficients =
        reshape(cos.(Float64.(1:(2 * fixed_full_rep.metadata.final_dimension))), :, 2)
    legacy_coefficients =
        reshape(sin.(Float64.(1:(2 * fixed_legacy_rep.metadata.final_dimension))), :, 2)
    transferred_full_to_legacy =
        transfer_orbitals(full_coefficients, fixed_full_rep, fixed_legacy_rep)
    transferred_legacy_to_full =
        transfer_orbitals(legacy_coefficients, fixed_legacy_rep, fixed_full_rep)

    @test full_to_legacy.matrix ≈ cross_overlap(fixed_legacy_rep, fixed_full_rep) atol =
          1.0e-10 rtol = 1.0e-10
    @test legacy_to_full.matrix ≈ cross_overlap(fixed_full_rep, fixed_legacy_rep) atol =
          1.0e-10 rtol = 1.0e-10
    @test transferred_full_to_legacy.diagnostics.transferred_residual_inf < 1.0e-10
    @test transferred_legacy_to_full.diagnostics.transferred_residual_inf < 1.0e-10

    hybrid_fixture = _atomic_hybrid_cartesian_representation_fixture()
    hybrid_self_projector = basis_projector(hybrid_fixture.full_rep, hybrid_fixture.full_rep)
    hybrid_self_dense = GaussletBases._cartesian_mixed_raw_cross_overlap_dense_reference(
        hybrid_fixture.full_rep,
        hybrid_fixture.full_rep,
    )
    hybrid_cross_dense = GaussletBases._cartesian_mixed_raw_cross_overlap_dense_reference(
        hybrid_fixture.legacy_rep,
        hybrid_fixture.full_rep,
    )
    hybrid_self_coefficients = reshape(
        sin.(Float64.(1:(2 * hybrid_fixture.full_rep.metadata.final_dimension))),
        :,
        2,
    )
    hybrid_self_transfer = transfer_orbitals(
        hybrid_self_coefficients,
        hybrid_fixture.full_rep,
        hybrid_fixture.full_rep,
    )
    hybrid_full_to_legacy = basis_projector(hybrid_fixture.full_rep, hybrid_fixture.legacy_rep)
    hybrid_parent_to_full =
        basis_projector(hybrid_fixture.fixed_full_rep, hybrid_fixture.full_rep)
    hybrid_full_to_parent =
        basis_projector(hybrid_fixture.full_rep, hybrid_fixture.fixed_full_rep)
    hybrid_transfer_from_projector =
        transfer_orbitals(hybrid_self_coefficients, hybrid_self_projector)

    @test hybrid_self_projector.matrix ≈ Matrix{Float64}(
        LinearAlgebra.I,
        hybrid_fixture.full_rep.metadata.final_dimension,
        hybrid_fixture.full_rep.metadata.final_dimension,
    ) atol = 1.0e-10 rtol = 1.0e-10
    @test cross_overlap(hybrid_fixture.full_rep, hybrid_fixture.full_rep) ≈
          hybrid_self_dense atol = 1.0e-10 rtol = 1.0e-10
    @test cross_overlap(hybrid_fixture.legacy_rep, hybrid_fixture.full_rep) ≈
          hybrid_cross_dense atol = 1.0e-10 rtol = 1.0e-10
    @test hybrid_full_to_legacy.matrix ≈ hybrid_cross_dense atol = 1.0e-10 rtol = 1.0e-10
    @test hybrid_self_transfer.coefficients ≈ hybrid_self_coefficients atol = 1.0e-10 rtol = 1.0e-10
    @test hybrid_self_transfer.projector !== nothing
    @test hybrid_self_transfer.projector.matrix ≈ hybrid_self_projector.matrix atol = 1.0e-10 rtol = 1.0e-10
    @test hybrid_self_transfer.diagnostics.transfer_path == :hybrid_mixed_raw_cross_overlap_transfer
    @test hybrid_self_transfer.diagnostics.projector_residual_inf < 1.0e-10
    @test hybrid_self_transfer.diagnostics.transferred_residual_inf < 1.0e-10
    @test hybrid_transfer_from_projector.coefficients ≈ hybrid_self_coefficients atol = 1.0e-10 rtol = 1.0e-10
    @test hybrid_transfer_from_projector.projector === hybrid_self_projector
    @test hybrid_full_to_legacy.matrix ≈
          cross_overlap(hybrid_fixture.legacy_rep, hybrid_fixture.full_rep) atol = 1.0e-10 rtol = 1.0e-10
    @test hybrid_parent_to_full.matrix ≈
          cross_overlap(hybrid_fixture.full_rep, hybrid_fixture.fixed_full_rep) atol = 1.0e-10 rtol = 1.0e-10
    @test hybrid_full_to_parent.matrix ≈
          cross_overlap(hybrid_fixture.fixed_full_rep, hybrid_fixture.full_rep) atol = 1.0e-10 rtol = 1.0e-10
end

@testset "Cartesian basis bundle export" begin
    square_basis, _square_source, square_fixed_block, _square_diagnostics =
        _axis_aligned_homonuclear_square_lattice_nested_fixture()
    square_basis_rep = basis_representation(square_basis)

    square_bundle = cartesian_basis_bundle_payload(
        square_basis;
        meta = (example = "test_cartesian_basis_bundle_basis_only",),
    )

    @test square_bundle.basis["format"] == "cartesian_basis_bundle_v1"
    @test square_bundle.basis["version"] == 1
    @test square_bundle.basis["basis_kind"] == "direct_product"
    @test square_bundle.basis["parent_kind"] == "cartesian_product_basis"
    @test square_bundle.basis["contraction_kind"] == "identity"
    @test size(square_bundle.basis["basis_centers"]) == size(square_basis_rep.metadata.basis_centers)
    @test length(square_bundle.basis["final_integral_weights"]) == square_basis_rep.metadata.final_dimension
    @test square_bundle.ham === nothing
    @test !square_bundle.meta["has_ham"]
    @test square_bundle.meta["example"] == "test_cartesian_basis_bundle_basis_only"

    fixed_bundle = cartesian_basis_bundle_payload(square_fixed_block)
    @test fixed_bundle.basis["basis_kind"] == "nested_fixed_block"
    @test fixed_bundle.basis["support_indices_present"]
    @test size(fixed_bundle.basis["support_states"], 2) == 3
    @test fixed_bundle.basis["final_integral_weights"] ≈ square_fixed_block.weights atol = 1.0e-12 rtol = 1.0e-12
    @test fixed_bundle.ham === nothing

    sparse_square_fixed_rep = basis_representation(_with_sparse_nested_coefficients(square_fixed_block))
    sparse_fixed_bundle = cartesian_basis_bundle_payload(sparse_square_fixed_rep)
    @test sparse_fixed_bundle.basis["coefficient_matrix"] isa SparseMatrixCSC{Float64,Int}
    @test Matrix(sparse_fixed_bundle.basis["coefficient_matrix"]) ≈
        Matrix(square_fixed_block.coefficient_matrix) atol = 1.0e-12 rtol = 1.0e-12

    diatomic_basis, diatomic_ops, _diatomic_check = _bond_aligned_diatomic_qw_fixture()
    operator_bundle = cartesian_basis_bundle_payload(
        diatomic_ops;
        meta = (example = "test_cartesian_basis_bundle_with_ham",),
    )
    operator_basis_only_bundle = cartesian_basis_bundle_payload(
        diatomic_ops;
        include_ham = false,
        meta = (example = "test_cartesian_basis_bundle_basis_only_from_ops",),
    )

    @test operator_bundle.basis["basis_kind"] == "direct_product"
    @test operator_bundle.ham !== nothing
    @test operator_bundle.ham["format"] == "cartesian_hamiltonian_bundle_v1"
    @test operator_bundle.ham["model_kind"] == "ordinary_cartesian_operators"
    @test size(operator_bundle.ham["overlap"]) == size(diatomic_ops.overlap)
    @test size(operator_bundle.ham["one_body_hamiltonian"]) == size(diatomic_ops.one_body_hamiltonian)
    @test size(operator_bundle.ham["interaction_matrix"]) == size(diatomic_ops.interaction_matrix)
    @test operator_bundle.ham["nuclear_term_storage"] == "by_center"
    @test operator_bundle.ham["default_nuclear_charges"] == [1.0, 1.0]
    @test operator_bundle.ham["nuclear_one_body_by_center/count"] == 2
    @test size(operator_bundle.ham["kinetic_one_body"]) == size(diatomic_ops.one_body_hamiltonian)
    @test operator_bundle.ham["basis_integral_weights"] == operator_bundle.basis["final_integral_weights"]
    @test operator_bundle.meta["has_ham"]

    mktempdir() do dir
        basis_only_path = joinpath(dir, "square_basis_only.jld2")
        sparse_fixed_path = joinpath(dir, "square_sparse_fixed.jld2")
        ops_path = joinpath(dir, "diatomic_ops_bundle.jld2")
        ops_basis_only_path = joinpath(dir, "diatomic_ops_basis_only_bundle.jld2")

        @test write_cartesian_basis_bundle_jld2(
            basis_only_path,
            square_basis;
            meta = (example = "test_cartesian_basis_bundle_basis_only",),
        ) == basis_only_path
        @test write_cartesian_basis_bundle_jld2(sparse_fixed_path, sparse_square_fixed_rep) ==
            sparse_fixed_path
        @test write_cartesian_basis_bundle_jld2(
            ops_path,
            diatomic_ops;
            meta = (example = "test_cartesian_basis_bundle_with_ham",),
        ) == ops_path
        @test write_cartesian_basis_bundle_jld2(
            ops_basis_only_path,
            diatomic_ops;
            include_ham = false,
            meta = (example = "test_cartesian_basis_bundle_basis_only_from_ops",),
        ) == ops_basis_only_path

        jldopen(basis_only_path, "r") do file
            top_keys = Set(
                key isa AbstractVector ? join(string.(key), "/") : string(key) for key in keys(file)
            )
            @test "basis" in top_keys
            @test "meta" in top_keys
            @test !("ham" in top_keys)
            @test String(file["basis/format"]) == "cartesian_basis_bundle_v1"
            @test Int(file["basis/version"]) == 1
            @test String(file["basis/basis_kind"]) == "direct_product"
            @test size(file["basis/basis_centers"]) == size(square_basis_rep.metadata.basis_centers)
            @test size(file["basis/final_integral_weights"]) == (square_basis_rep.metadata.final_dimension,)
            @test String(file["basis/axes/x/format"]) == "basis_representation_1d_v1"
            @test String(file["meta/producer"]) ==
                "GaussletBases.write_cartesian_basis_bundle_jld2"
        end

        jldopen(sparse_fixed_path, "r") do file
            basis_values = GaussletBases._cartesian_jld_group_values(file["basis"])
            meta_values = GaussletBases._cartesian_jld_group_values(file["meta"])
            @test file["basis/coefficient_matrix"] isa SparseMatrixCSC{Float64,Int}
            @test Set(keys(basis_values)) == Set(keys(sparse_fixed_bundle.basis))
            @test Set(keys(meta_values)) == Set(keys(sparse_fixed_bundle.meta))
            @test basis_values["final_integral_weights"] ≈
                sparse_fixed_bundle.basis["final_integral_weights"] atol = 1.0e-12 rtol = 1.0e-12
        end

        sparse_fixed_bundle_roundtrip = read_cartesian_basis_bundle(sparse_fixed_path)
        @test sparse_fixed_bundle_roundtrip.basis.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
        @test sparse_fixed_bundle_roundtrip.basis.coefficient_matrix ==
            sparse_square_fixed_rep.coefficient_matrix
        @test cross_overlap(sparse_fixed_bundle_roundtrip, sparse_fixed_bundle_roundtrip) ≈
            cross_overlap(sparse_square_fixed_rep, sparse_square_fixed_rep) atol = 1.0e-10 rtol = 1.0e-10

        jldopen(ops_path, "r") do file
            top_keys = Set(
                key isa AbstractVector ? join(string.(key), "/") : string(key) for key in keys(file)
            )
            ham_values = GaussletBases._cartesian_jld_group_values(file["ham"])
            @test "basis" in top_keys
            @test "ham" in top_keys
            @test "meta" in top_keys
            @test Set(keys(ham_values)) == Set(keys(operator_bundle.ham))
            @test String(file["ham/format"]) == "cartesian_hamiltonian_bundle_v1"
            @test String(file["ham/model_kind"]) == "ordinary_cartesian_operators"
            @test size(file["ham/overlap"]) == size(diatomic_ops.overlap)
            @test size(file["ham/one_body_hamiltonian"]) == size(diatomic_ops.one_body_hamiltonian)
            @test size(file["ham/interaction_matrix"]) == size(diatomic_ops.interaction_matrix)
            @test String(file["ham/nuclear_term_storage"]) == "by_center"
            @test Int(file["ham/nuclear_one_body_by_center/count"]) == 2
            @test size(file["ham/kinetic_one_body"]) == size(diatomic_ops.one_body_hamiltonian)
            @test String(file["meta/manifest/contract/format"]) == "cartesian_basis_bundle_v1"
            @test Bool(file["meta/has_ham"])
        end

        jldopen(ops_basis_only_path, "r") do file
            top_keys = Set(
                key isa AbstractVector ? join(string.(key), "/") : string(key) for key in keys(file)
            )
            basis_values = GaussletBases._cartesian_jld_group_values(file["basis"])
            meta_values = GaussletBases._cartesian_jld_group_values(file["meta"])
            @test "basis" in top_keys
            @test "meta" in top_keys
            @test !("ham" in top_keys)
            @test Set(keys(basis_values)) == Set(keys(operator_basis_only_bundle.basis))
            @test Set(keys(meta_values)) == Set(keys(operator_basis_only_bundle.meta))
            @test !Bool(file["meta/has_ham"])
            @test String(file["meta/example"]) == "test_cartesian_basis_bundle_basis_only_from_ops"
        end

        ops_basis_only_bundle_roundtrip = read_cartesian_basis_bundle(ops_basis_only_path)
        @test ops_basis_only_bundle_roundtrip.ham === nothing
        @test cross_overlap(ops_basis_only_bundle_roundtrip, ops_basis_only_bundle_roundtrip) ≈
            diatomic_ops.overlap atol = 1.0e-10 rtol = 1.0e-10
    end

    hybrid_fixture = _atomic_hybrid_cartesian_representation_fixture()
    hybrid_bundle = cartesian_basis_bundle_payload(
        hybrid_fixture.full_ops;
        meta = (example = "test_cartesian_hybrid_bundle",),
    )

    @test hybrid_bundle.basis["basis_kind"] == "hybrid_residual"
    @test hybrid_bundle.basis["parent_kind"] == "cartesian_plus_supplement_raw"
    @test hybrid_bundle.basis["parent/format"] == "cartesian_plus_supplement_raw_v1"
    @test hybrid_bundle.basis["parent/cartesian/format"] == "cartesian_basis_bundle_v1"
    @test hybrid_bundle.basis["parent/supplement/format"] ==
        "cartesian_gaussian_shell_supplement_v1"
    @test haskey(hybrid_bundle.basis, "parent/cartesian_supplement_axis_tables/x")
    @test haskey(hybrid_bundle.basis, "parent/cartesian_supplement_axis_tables/y")
    @test haskey(hybrid_bundle.basis, "parent/cartesian_supplement_axis_tables/z")
    @test haskey(hybrid_bundle.basis, "parent/exact_cartesian_supplement_overlap")
    @test haskey(hybrid_bundle.basis, "parent/exact_supplement_overlap")
    @test hybrid_bundle.basis["parent/supplement/orbital_count"] ==
        size(hybrid_fixture.full_ops.raw_to_final, 1) - hybrid_fixture.full_ops.gausslet_count
    @test hybrid_bundle.ham !== nothing
    @test hybrid_bundle.ham["model_kind"] == "ordinary_cartesian_operators"
    @test size(hybrid_bundle.ham["overlap"]) == size(hybrid_fixture.full_ops.overlap)
    @test hybrid_bundle.meta["example"] == "test_cartesian_hybrid_bundle"

    mktempdir() do dir
        hybrid_path = joinpath(dir, "atomic_hybrid_ops_bundle.jld2")

        @test write_cartesian_basis_bundle_jld2(
            hybrid_path,
            hybrid_fixture.full_ops;
            meta = (example = "test_cartesian_hybrid_bundle",),
        ) == hybrid_path

        jldopen(hybrid_path, "r") do file
            @test String(file["basis/parent_kind"]) == "cartesian_plus_supplement_raw"
            @test String(file["basis/parent/format"]) == "cartesian_plus_supplement_raw_v1"
            @test String(file["basis/parent/cartesian/format"]) == "cartesian_basis_bundle_v1"
            @test String(file["basis/parent/supplement/format"]) ==
                "cartesian_gaussian_shell_supplement_v1"
            @test size(file["basis/parent/cartesian_supplement_axis_tables/x"], 2) ==
                size(hybrid_fixture.full_ops.raw_to_final, 1) - hybrid_fixture.full_ops.gausslet_count
            @test size(file["basis/parent/cartesian_supplement_axis_tables/y"], 2) ==
                size(hybrid_fixture.full_ops.raw_to_final, 1) - hybrid_fixture.full_ops.gausslet_count
            @test size(file["basis/parent/cartesian_supplement_axis_tables/z"], 2) ==
                size(hybrid_fixture.full_ops.raw_to_final, 1) - hybrid_fixture.full_ops.gausslet_count
            @test size(file["basis/parent/exact_cartesian_supplement_overlap"]) ==
                (hybrid_fixture.full_ops.gausslet_count,
                 size(hybrid_fixture.full_ops.raw_to_final, 1) - hybrid_fixture.full_ops.gausslet_count)
            @test size(file["basis/parent/exact_supplement_overlap"]) ==
                (size(hybrid_fixture.full_ops.raw_to_final, 1) - hybrid_fixture.full_ops.gausslet_count,
                 size(hybrid_fixture.full_ops.raw_to_final, 1) - hybrid_fixture.full_ops.gausslet_count)
            @test Int(file["basis/parent/supplement/orbital_count"]) ==
                size(hybrid_fixture.full_ops.raw_to_final, 1) - hybrid_fixture.full_ops.gausslet_count
            @test size(file["ham/overlap"]) == size(hybrid_fixture.full_ops.overlap)
            @test size(file["ham/one_body_hamiltonian"]) ==
                size(hybrid_fixture.full_ops.one_body_hamiltonian)
        end
    end

    bond_aligned_hybrid_fixture = _bond_aligned_diatomic_nested_hybrid_bundle_fixture()
    bond_aligned_hybrid_trimmed_fixture =
        _bond_aligned_diatomic_nested_hybrid_bundle_fixture(; max_width = 1.0)
    bond_aligned_hybrid_supplement3d =
        GaussletBases._bond_aligned_diatomic_cartesian_shell_supplement_3d(
            bond_aligned_hybrid_fixture.supplement,
        )
    bond_aligned_hybrid_bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(
        bond_aligned_hybrid_fixture.basis,
        bond_aligned_hybrid_fixture.hybrid_ops.expansion;
        gausslet_backend = bond_aligned_hybrid_fixture.hybrid_ops.gausslet_backend,
    )
    bond_aligned_hybrid_overlap_blocks =
        GaussletBases._qwrg_diatomic_cartesian_shell_overlap_blocks_3d(
            bond_aligned_hybrid_bundles,
            bond_aligned_hybrid_supplement3d,
            bond_aligned_hybrid_fixture.basis,
            bond_aligned_hybrid_fixture.hybrid_ops.expansion,
        )
    bond_aligned_hybrid_full_blocks = GaussletBases._qwrg_diatomic_cartesian_shell_blocks_3d(
        bond_aligned_hybrid_bundles,
        bond_aligned_hybrid_supplement3d,
        bond_aligned_hybrid_fixture.basis,
        bond_aligned_hybrid_fixture.hybrid_ops.expansion,
        bond_aligned_hybrid_fixture.hybrid_ops.nuclear_charges,
    )
    bond_aligned_hybrid_bundle = cartesian_basis_bundle_payload(
        bond_aligned_hybrid_fixture.hybrid_ops;
        meta = (example = "test_cartesian_bond_aligned_diatomic_hybrid_bundle",),
    )
    bond_aligned_hybrid_trimmed_bundle = cartesian_basis_bundle_payload(
        bond_aligned_hybrid_trimmed_fixture.hybrid_ops;
        meta = (example = "test_cartesian_bond_aligned_diatomic_hybrid_bundle_trimmed",),
    )

    @test bond_aligned_hybrid_bundle.basis["basis_kind"] == "hybrid_residual"
    @test bond_aligned_hybrid_bundle.basis["parent_kind"] == "cartesian_plus_supplement_raw"
    @test bond_aligned_hybrid_bundle.basis["parent/format"] == "cartesian_plus_supplement_raw_v1"
    @test bond_aligned_hybrid_bundle.basis["parent/supplement/format"] ==
        "cartesian_gaussian_shell_supplement_v1"
    @test haskey(bond_aligned_hybrid_bundle.basis, "parent/cartesian_supplement_axis_tables/x")
    @test haskey(bond_aligned_hybrid_bundle.basis, "parent/exact_cartesian_supplement_overlap")
    @test haskey(bond_aligned_hybrid_bundle.basis, "parent/exact_supplement_overlap")
    @test bond_aligned_hybrid_overlap_blocks.overlap_ga ≈
        bond_aligned_hybrid_full_blocks.overlap_ga atol = 1.0e-12 rtol = 1.0e-12
    @test bond_aligned_hybrid_overlap_blocks.overlap_aa ≈
        bond_aligned_hybrid_full_blocks.overlap_aa atol = 1.0e-12 rtol = 1.0e-12
    @test bond_aligned_hybrid_trimmed_bundle.basis["parent/supplement/metadata/max_width"] == 1.0
    @test Int(bond_aligned_hybrid_trimmed_bundle.basis["parent/supplement/orbital_count"]) <
        Int(bond_aligned_hybrid_bundle.basis["parent/supplement/orbital_count"])

    mktempdir() do dir
        hybrid_path = joinpath(dir, "bond_aligned_diatomic_hybrid_ops_bundle.jld2")
        hybrid_trimmed_path =
            joinpath(dir, "bond_aligned_diatomic_hybrid_ops_bundle_trimmed.jld2")

        @test write_cartesian_basis_bundle_jld2(
            hybrid_path,
            bond_aligned_hybrid_fixture.hybrid_ops;
            meta = (example = "test_cartesian_bond_aligned_diatomic_hybrid_bundle",),
        ) == hybrid_path
        @test write_cartesian_basis_bundle_jld2(
            hybrid_trimmed_path,
            bond_aligned_hybrid_trimmed_fixture.hybrid_ops;
            meta = (example = "test_cartesian_bond_aligned_diatomic_hybrid_bundle_trimmed",),
        ) == hybrid_trimmed_path

        jldopen(hybrid_path, "r") do file
            @test String(file["basis/parent_kind"]) == "cartesian_plus_supplement_raw"
            @test String(file["basis/parent/format"]) == "cartesian_plus_supplement_raw_v1"
            @test String(file["basis/parent/supplement/format"]) ==
                "cartesian_gaussian_shell_supplement_v1"
            @test Int(file["basis/parent/supplement/orbital_count"]) ==
                size(bond_aligned_hybrid_fixture.hybrid_ops.raw_to_final, 1) -
                bond_aligned_hybrid_fixture.hybrid_ops.gausslet_count
            @test size(file["basis/parent/exact_cartesian_supplement_overlap"]) ==
                (bond_aligned_hybrid_fixture.hybrid_ops.gausslet_count,
                 size(bond_aligned_hybrid_fixture.hybrid_ops.raw_to_final, 1) -
                 bond_aligned_hybrid_fixture.hybrid_ops.gausslet_count)
            @test size(file["basis/parent/exact_supplement_overlap"]) ==
                (size(bond_aligned_hybrid_fixture.hybrid_ops.raw_to_final, 1) -
                 bond_aligned_hybrid_fixture.hybrid_ops.gausslet_count,
                 size(bond_aligned_hybrid_fixture.hybrid_ops.raw_to_final, 1) -
                 bond_aligned_hybrid_fixture.hybrid_ops.gausslet_count)
            @test String(file["meta/example"]) ==
                "test_cartesian_bond_aligned_diatomic_hybrid_bundle"
        end

        jldopen(hybrid_trimmed_path, "r") do file
            @test Float64(file["basis/parent/supplement/metadata/max_width"]) == 1.0
            @test Int(file["basis/parent/supplement/orbital_count"]) <
                Int(bond_aligned_hybrid_bundle.basis["parent/supplement/orbital_count"])
            @test String(file["meta/example"]) ==
                "test_cartesian_bond_aligned_diatomic_hybrid_bundle_trimmed"
        end
    end
end

@testset "Cartesian basis bundle overlap and projector" begin
    square_basis, _square_source, square_fixed_block, _square_diagnostics =
        _axis_aligned_homonuclear_square_lattice_nested_fixture()
    square_basis_rep = basis_representation(square_basis)
    square_fixed_rep = basis_representation(square_fixed_block)

    diatomic_basis14, diatomic_ops14, _diatomic_check14 =
        _bond_aligned_diatomic_qw_fixture(; bond_length = 1.4)
    diatomic_basis20, _diatomic_ops20, _diatomic_check20 =
        _bond_aligned_diatomic_qw_fixture(; bond_length = 2.0)
    diatomic_rep14 = basis_representation(diatomic_basis14)
    diatomic_rep20 = basis_representation(diatomic_basis20)
    bond_aligned_hybrid_fixture = _bond_aligned_diatomic_nested_hybrid_bundle_fixture()

    dir = mktempdir()
    try
        square_path = joinpath(dir, "square_basis.jld2")
        square_fixed_path = joinpath(dir, "square_fixed.jld2")
        diatomic14_path = joinpath(dir, "diatomic14.jld2")
        diatomic20_path = joinpath(dir, "diatomic20.jld2")
        diatomic_ops_path = joinpath(dir, "diatomic_ops.jld2")
        atomic_fixed_full_path = joinpath(dir, "atomic_fixed_full.jld2")
        atomic_hybrid_full_path = joinpath(dir, "atomic_hybrid_full.jld2")
        atomic_hybrid_legacy_path = joinpath(dir, "atomic_hybrid_legacy.jld2")
        bond_aligned_hybrid_fixed_path = joinpath(dir, "bond_aligned_hybrid_fixed.jld2")
        bond_aligned_hybrid_path = joinpath(dir, "bond_aligned_hybrid_ops.jld2")

        write_cartesian_basis_bundle_jld2(square_path, square_basis)
        write_cartesian_basis_bundle_jld2(square_fixed_path, square_fixed_block)
        write_cartesian_basis_bundle_jld2(diatomic14_path, diatomic_basis14)
        write_cartesian_basis_bundle_jld2(diatomic20_path, diatomic_basis20)
        write_cartesian_basis_bundle_jld2(diatomic_ops_path, diatomic_ops14)
        hybrid_fixture = _atomic_hybrid_cartesian_representation_fixture()
        write_cartesian_basis_bundle_jld2(atomic_fixed_full_path, hybrid_fixture.fixed_full)
        write_cartesian_basis_bundle_jld2(atomic_hybrid_full_path, hybrid_fixture.full_ops)
        write_cartesian_basis_bundle_jld2(atomic_hybrid_legacy_path, hybrid_fixture.legacy_ops)
        write_cartesian_basis_bundle_jld2(
            bond_aligned_hybrid_fixed_path,
            bond_aligned_hybrid_fixture.fixed_block,
        )
        write_cartesian_basis_bundle_jld2(
            bond_aligned_hybrid_path,
            bond_aligned_hybrid_fixture.hybrid_ops,
        )

        square_bundle = read_cartesian_basis_bundle(square_path)
        square_fixed_bundle = read_cartesian_basis_bundle(square_fixed_path)
        diatomic14_bundle = read_cartesian_basis_bundle(diatomic14_path)
        diatomic20_bundle = read_cartesian_basis_bundle(diatomic20_path)
        diatomic_ops_bundle = read_cartesian_basis_bundle(diatomic_ops_path)
        atomic_hybrid_full_bundle = read_cartesian_basis_bundle(atomic_hybrid_full_path)
        atomic_hybrid_legacy_bundle = read_cartesian_basis_bundle(atomic_hybrid_legacy_path)
        bond_aligned_hybrid_bundle = read_cartesian_basis_bundle(bond_aligned_hybrid_path)

        @test square_bundle.path == abspath(square_path)
        @test square_bundle.diagnostics.basis_kind == :direct_product
        @test square_bundle.diagnostics.final_dimension == square_basis_rep.metadata.final_dimension
        @test square_bundle.ham === nothing
        @test diatomic_ops_bundle.ham !== nothing
        @test diatomic_ops_bundle.ham["model_kind"] == "ordinary_cartesian_operators"
        @test diatomic_ops_bundle.diagnostics.has_ham

        diatomic_atom_a_ops = ordinary_cartesian_qiu_white_operators(
            diatomic_basis14;
            nuclear_charges = [1.0, 0.0],
            nuclear_term_storage = :total_only,
            interaction_treatment = :ggt_nearest,
        )
        diatomic_atom_b_ops = ordinary_cartesian_qiu_white_operators(
            diatomic_basis14;
            nuclear_charges = [0.0, 1.0],
            nuclear_term_storage = :total_only,
            interaction_treatment = :ggt_nearest,
        )
        @test assembled_one_body_hamiltonian(diatomic_ops14) ≈
              diatomic_ops14.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
        @test assembled_one_body_hamiltonian(diatomic_ops_bundle) ≈
              diatomic_ops14.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
        @test assembled_one_body_hamiltonian(diatomic_ops14; nuclear_charges = [1.0, 0.0]) ≈
              diatomic_atom_a_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(diatomic_ops_bundle; nuclear_charges = [1.0, 0.0]) ≈
              diatomic_atom_a_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(diatomic_ops14; nuclear_charges = [0.0, 1.0]) ≈
              diatomic_atom_b_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(diatomic_ops_bundle; nuclear_charges = [0.0, 1.0]) ≈
              diatomic_atom_b_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10

        hybrid_atom_a_ops = ordinary_cartesian_qiu_white_operators(
            bond_aligned_hybrid_fixture.fixed_block,
            bond_aligned_hybrid_fixture.supplement;
            nuclear_charges = [1.0, 0.0],
            nuclear_term_storage = :total_only,
            interaction_treatment = :ggt_nearest,
        )
        hybrid_atom_b_ops = ordinary_cartesian_qiu_white_operators(
            bond_aligned_hybrid_fixture.fixed_block,
            bond_aligned_hybrid_fixture.supplement;
            nuclear_charges = [0.0, 1.0],
            nuclear_term_storage = :total_only,
            interaction_treatment = :ggt_nearest,
        )
        @test assembled_one_body_hamiltonian(bond_aligned_hybrid_fixture.hybrid_ops) ≈
              bond_aligned_hybrid_fixture.hybrid_ops.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
        @test assembled_one_body_hamiltonian(bond_aligned_hybrid_bundle) ≈
              bond_aligned_hybrid_fixture.hybrid_ops.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
        @test assembled_one_body_hamiltonian(
                  bond_aligned_hybrid_fixture.hybrid_ops;
                  nuclear_charges = [1.0, 0.0],
              ) ≈ hybrid_atom_a_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(
                  bond_aligned_hybrid_bundle;
                  nuclear_charges = [1.0, 0.0],
              ) ≈ hybrid_atom_a_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(
                  bond_aligned_hybrid_fixture.hybrid_ops;
                  nuclear_charges = [0.0, 1.0],
              ) ≈ hybrid_atom_b_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(
                  bond_aligned_hybrid_bundle;
                  nuclear_charges = [0.0, 1.0],
              ) ≈ hybrid_atom_b_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10

        loaded_square_rep = load_cartesian_basis_representation(square_path)
        @test loaded_square_rep.metadata.basis_kind == square_basis_rep.metadata.basis_kind
        @test loaded_square_rep.metadata.final_dimension == square_basis_rep.metadata.final_dimension
        @test loaded_square_rep.metadata.parent_kind == square_basis_rep.metadata.parent_kind

        direct_self_disk = cross_overlap(square_bundle, square_bundle)
        direct_cross_disk = cross_overlap(diatomic14_bundle, diatomic20_bundle)
        nested_cross_disk = cross_overlap(square_path, square_fixed_path)

        @test direct_self_disk ≈ cross_overlap(square_basis_rep, square_basis_rep) atol = 1.0e-10 rtol = 1.0e-10
        @test direct_cross_disk ≈ cross_overlap(diatomic_rep14, diatomic_rep20) atol = 1.0e-10 rtol = 1.0e-10
        @test nested_cross_disk ≈ cross_overlap(square_basis_rep, square_fixed_rep) atol = 1.0e-10 rtol = 1.0e-10

        disk_projector = basis_projector(square_fixed_path, square_path)
        memory_projector = basis_projector(square_fixed_rep, square_basis_rep)
        @test disk_projector.matrix ≈ memory_projector.matrix atol = 1.0e-10 rtol = 1.0e-10
        @test disk_projector.diagnostics.transfer_path == memory_projector.diagnostics.transfer_path

        fixed_coefficients = cos.(Float64.(1:square_fixed_rep.metadata.final_dimension))
        disk_transfer = transfer_orbitals(fixed_coefficients, square_fixed_path, square_path)
        memory_transfer = transfer_orbitals(fixed_coefficients, square_fixed_rep, square_basis_rep)
        @test disk_transfer.coefficients ≈ memory_transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test disk_transfer.diagnostics.transferred_residual_inf < 1.0e-10

        @test atomic_hybrid_full_bundle.diagnostics.parent_kind == :cartesian_plus_supplement_raw
        @test atomic_hybrid_full_bundle.ham !== nothing
        @test hasproperty(
            atomic_hybrid_full_bundle.basis.parent_data,
            :cartesian_supplement_axis_tables,
        )
        @test bond_aligned_hybrid_bundle.diagnostics.parent_kind == :cartesian_plus_supplement_raw
        @test bond_aligned_hybrid_bundle.ham !== nothing
        @test hasproperty(
            bond_aligned_hybrid_bundle.basis.parent_data,
            :cartesian_supplement_axis_tables,
        )

        disk_hybrid_self = cross_overlap(atomic_hybrid_full_path, atomic_hybrid_full_path)
        memory_hybrid_self = cross_overlap(hybrid_fixture.full_rep, hybrid_fixture.full_rep)
        @test disk_hybrid_self ≈ memory_hybrid_self atol = 1.0e-10 rtol = 1.0e-10

        disk_hybrid_cross = cross_overlap(atomic_hybrid_full_path, atomic_hybrid_legacy_path)
        memory_hybrid_cross = cross_overlap(hybrid_fixture.full_rep, hybrid_fixture.legacy_rep)
        @test disk_hybrid_cross ≈ memory_hybrid_cross atol = 1.0e-10 rtol = 1.0e-10

        disk_hybrid_parent = cross_overlap(atomic_fixed_full_path, atomic_hybrid_full_path)
        memory_hybrid_parent = cross_overlap(hybrid_fixture.fixed_full_rep, hybrid_fixture.full_rep)
        @test disk_hybrid_parent ≈ memory_hybrid_parent atol = 1.0e-10 rtol = 1.0e-10

        disk_hybrid_projector =
            basis_projector(atomic_hybrid_full_path, atomic_hybrid_legacy_path)
        memory_hybrid_projector =
            basis_projector(hybrid_fixture.full_rep, hybrid_fixture.legacy_rep)
        @test disk_hybrid_projector.matrix ≈ memory_hybrid_projector.matrix atol = 1.0e-10 rtol = 1.0e-10
        @test disk_hybrid_projector.diagnostics.transfer_path ==
            memory_hybrid_projector.diagnostics.transfer_path

        hybrid_coefficients = cos.(Float64.(1:hybrid_fixture.full_rep.metadata.final_dimension))
        disk_hybrid_transfer = transfer_orbitals(
            hybrid_coefficients,
            atomic_hybrid_full_path,
            atomic_hybrid_legacy_path,
        )
        memory_hybrid_transfer = transfer_orbitals(
            hybrid_coefficients,
            hybrid_fixture.full_rep,
            hybrid_fixture.legacy_rep,
        )
        disk_hybrid_fast_transfer = transfer_orbitals(
            hybrid_coefficients,
            atomic_hybrid_full_path,
            atomic_hybrid_legacy_path;
            materialize_projector = false,
        )
        bundle_hybrid_fast_transfer = transfer_orbitals(
            hybrid_coefficients,
            atomic_hybrid_full_bundle,
            atomic_hybrid_legacy_bundle;
            materialize_projector = false,
        )
        memory_hybrid_fast_transfer = transfer_orbitals(
            hybrid_coefficients,
            hybrid_fixture.full_rep,
            hybrid_fixture.legacy_rep;
            materialize_projector = false,
        )
        @test disk_hybrid_transfer.coefficients ≈
            memory_hybrid_transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test disk_hybrid_transfer.projector !== nothing
        @test memory_hybrid_transfer.projector !== nothing
        @test disk_hybrid_transfer.projector.matrix ≈
            memory_hybrid_transfer.projector.matrix atol = 1.0e-10 rtol = 1.0e-10
        @test disk_hybrid_transfer.diagnostics.transferred_residual_inf < 1.0e-10
        @test memory_hybrid_fast_transfer.coefficients ≈
            memory_hybrid_transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test disk_hybrid_fast_transfer.coefficients ≈
            memory_hybrid_fast_transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test bundle_hybrid_fast_transfer.coefficients ≈
            memory_hybrid_fast_transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test disk_hybrid_fast_transfer.projector === nothing
        @test bundle_hybrid_fast_transfer.projector === nothing
        @test memory_hybrid_fast_transfer.projector === nothing
        @test disk_hybrid_fast_transfer.diagnostics.transfer_path == :hybrid_mixed_raw_cross_overlap_transfer
        @test bundle_hybrid_fast_transfer.diagnostics.transfer_path == :hybrid_mixed_raw_cross_overlap_transfer
        @test isnan(disk_hybrid_fast_transfer.diagnostics.projector_residual_inf)
        @test isnan(bundle_hybrid_fast_transfer.diagnostics.projector_residual_inf)
        @test isnan(disk_hybrid_fast_transfer.diagnostics.transferred_residual_inf)
        @test isnan(bundle_hybrid_fast_transfer.diagnostics.transferred_residual_inf)

        disk_bond_aligned_hybrid_self =
            cross_overlap(bond_aligned_hybrid_path, bond_aligned_hybrid_path)
        memory_bond_aligned_hybrid_self =
            cross_overlap(
                bond_aligned_hybrid_fixture.hybrid_rep,
                bond_aligned_hybrid_fixture.hybrid_rep,
            )
        @test disk_bond_aligned_hybrid_self ≈
            memory_bond_aligned_hybrid_self atol = 1.0e-10 rtol = 1.0e-10

        disk_bond_aligned_hybrid_cross =
            cross_overlap(bond_aligned_hybrid_fixed_path, bond_aligned_hybrid_path)
        memory_bond_aligned_hybrid_cross =
            cross_overlap(
                bond_aligned_hybrid_fixture.fixed_rep,
                bond_aligned_hybrid_fixture.hybrid_rep,
            )
        @test disk_bond_aligned_hybrid_cross ≈
            memory_bond_aligned_hybrid_cross atol = 1.0e-10 rtol = 1.0e-10

        disk_bond_aligned_hybrid_projector =
            basis_projector(bond_aligned_hybrid_fixed_path, bond_aligned_hybrid_path)
        memory_bond_aligned_hybrid_projector =
            basis_projector(
                bond_aligned_hybrid_fixture.fixed_rep,
                bond_aligned_hybrid_fixture.hybrid_rep,
            )
        @test disk_bond_aligned_hybrid_projector.matrix ≈
            memory_bond_aligned_hybrid_projector.matrix atol = 1.0e-10 rtol = 1.0e-10
        @test disk_bond_aligned_hybrid_projector.diagnostics.transfer_path ==
            memory_bond_aligned_hybrid_projector.diagnostics.transfer_path

        bond_aligned_coefficients =
            cos.(Float64.(1:bond_aligned_hybrid_fixture.fixed_rep.metadata.final_dimension))
        disk_bond_aligned_hybrid_transfer = transfer_orbitals(
            bond_aligned_coefficients,
            bond_aligned_hybrid_fixed_path,
            bond_aligned_hybrid_path,
        )
        memory_bond_aligned_hybrid_transfer = transfer_orbitals(
            bond_aligned_coefficients,
            bond_aligned_hybrid_fixture.fixed_rep,
            bond_aligned_hybrid_fixture.hybrid_rep,
        )
        @test disk_bond_aligned_hybrid_transfer.coefficients ≈
            memory_bond_aligned_hybrid_transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test disk_bond_aligned_hybrid_transfer.projector.matrix ≈
            memory_bond_aligned_hybrid_transfer.projector.matrix atol = 1.0e-10 rtol = 1.0e-10
        @test disk_bond_aligned_hybrid_transfer.diagnostics.transferred_residual_inf < 1.0e-10
    finally
        rm(dir; recursive = true, force = true)
    end
end

@testset "Atomic direct-product He extent change is not an outer-only identity" begin
    fixture = _atomic_direct_product_he_extent_change_contract_fixture()

    @test fixture.source_count == 3
    @test fixture.target_count == 5
    @test fixture.shared_slice == 2:4
    @test fixture.centers_subset
    @test !fixture.weights_subset
    @test !fixture.coefficient_core_match
end

@testset "Atomic hybrid He orbital transfer remains stable across same-parent different-final-contraction change" begin
    fixture = _atomic_hybrid_he_same_parent_stress_fixture()

    @test fixture.source_ops isa OrdinaryCartesianOperators3D
    @test fixture.target_ops isa OrdinaryCartesianOperators3D
    @test fixture.source_rep.metadata.parent_kind == :cartesian_plus_supplement_raw
    @test fixture.target_rep.metadata.parent_kind == :cartesian_plus_supplement_raw
    @test fixture.source_working_box == 2:6
    @test fixture.target_working_box == (1:7, 1:7, 1:7)
    @test length(fixture.source_ops.orbital_data) == 134
    @test length(fixture.target_ops.orbital_data) == 232
    @test fixture.transfer.diagnostics.transfer_path == :hybrid_mixed_raw_cross_overlap_transfer
    @test fixture.transfer.diagnostics.transferred_residual_inf < 1.0e-10

    @test fixture.source_observables.metric_norm_error < 1.0e-12
    @test fixture.target_observables.metric_norm_error < 1.0e-12
    @test fixture.aligned_transferred_observables.metric_norm_error < 1.0e-12

    source_self_overlap = cross_overlap(fixture.source_rep, fixture.source_rep)
    target_self_overlap = cross_overlap(fixture.target_rep, fixture.target_rep)
    cross_overlap_source_target = cross_overlap(fixture.source_rep, fixture.target_rep)
    @test source_self_overlap ≈ fixture.source_ops.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test target_self_overlap ≈ fixture.target_ops.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test maximum(abs, cross_overlap_source_target) <= 1.0 + 1.0e-10
    @test maximum(svdvals(cross_overlap_source_target)) <= 1.0 + 1.0e-10

    @test fixture.target_observables.total < fixture.source_observables.total
    @test fixture.aligned_overlap_to_target > 0.999995
    @test abs(
        fixture.aligned_transferred_observables.one_body - fixture.target_observables.one_body,
    ) < 1.0e-4
    @test abs(
        fixture.aligned_transferred_observables.vee - fixture.target_observables.vee,
    ) < 2.0e-4
    @test abs(
        fixture.aligned_transferred_observables.total - fixture.target_observables.total,
    ) < 5.0e-4

    mktempdir() do dir
        source_path = joinpath(dir, "he_source_hybrid.jld2")
        target_path = joinpath(dir, "he_target_hybrid.jld2")

        write_cartesian_basis_bundle_jld2(source_path, fixture.source_ops)
        write_cartesian_basis_bundle_jld2(target_path, fixture.target_ops)

        source_bundle = read_cartesian_basis_bundle(source_path)
        target_bundle = read_cartesian_basis_bundle(target_path)
        @test hasproperty(source_bundle.basis.parent_data, :exact_cartesian_supplement_overlap)
        @test hasproperty(source_bundle.basis.parent_data, :exact_supplement_overlap)
        @test hasproperty(target_bundle.basis.parent_data, :exact_cartesian_supplement_overlap)
        @test hasproperty(target_bundle.basis.parent_data, :exact_supplement_overlap)

        disk_source_self = cross_overlap(source_path, source_path)
        disk_target_self = cross_overlap(target_path, target_path)
        disk_cross = cross_overlap(source_path, target_path)
        @test disk_source_self ≈ fixture.source_ops.overlap atol = 1.0e-10 rtol = 1.0e-10
        @test disk_target_self ≈ fixture.target_ops.overlap atol = 1.0e-10 rtol = 1.0e-10
        @test disk_source_self ≈ source_self_overlap atol = 1.0e-10 rtol = 1.0e-10
        @test disk_target_self ≈ target_self_overlap atol = 1.0e-10 rtol = 1.0e-10
        @test disk_cross ≈ cross_overlap_source_target atol = 1.0e-10 rtol = 1.0e-10
        @test maximum(abs, disk_cross) <= 1.0 + 1.0e-10
        @test maximum(svdvals(disk_cross)) <= 1.0 + 1.0e-10

        disk_transfer = transfer_orbitals(
            fixture.source_observables.orbital,
            source_path,
            target_path,
        )

        @test disk_transfer.coefficients ≈ fixture.transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test disk_transfer.diagnostics.transfer_path ==
            fixture.transfer.diagnostics.transfer_path
    end
end

@testset "Mapped ordinary Cartesian 1D working representation uses localized Gaussian contract" begin
    mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0)
    basis_a = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 5,
            mapping = mapping,
            reference_spacing = 1.0,
        ),
    )
    basis_b = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 7,
            mapping = mapping,
            reference_spacing = 1.0,
        ),
    )

    rep_a = GaussletBases._mapped_ordinary_working_basis_representation(basis_a)
    rep_b = GaussletBases._mapped_ordinary_working_basis_representation(basis_b)
    S_AA = cross_overlap(rep_a, rep_a)
    S_BB = cross_overlap(rep_b, rep_b)
    S_AB = cross_overlap(rep_a, rep_b)
    I_A = Matrix{Float64}(I, size(S_AA, 1), size(S_AA, 2))
    I_B = Matrix{Float64}(I, size(S_BB, 1), size(S_BB, 2))

    @test all(primitive -> primitive isa Gaussian, primitives(primitive_set(rep_a)))
    @test all(primitive -> primitive isa Gaussian, primitives(primitive_set(rep_b)))
    @test norm(S_AA - I_A, Inf) < 1.0e-12
    @test norm(S_BB - I_B, Inf) < 1.0e-12
    @test maximum(svdvals(S_AB)) <= 1.0 + 1.0e-10

    fixture = _atomic_hybrid_he_same_parent_stress_fixture()
    @test all(
        primitive -> primitive isa Gaussian,
        primitives(primitive_set(fixture.source_rep.axis_representations.x)),
    )
    @test all(
        primitive -> primitive isa Gaussian,
        primitives(primitive_set(fixture.target_rep.axis_representations.x)),
    )
end

@testset "One-center atomic factorized direct packet kernel" begin
    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 13,
            mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )

    shell = GaussletBases._nested_rectangular_shell(
        bundle,
        2:12,
        2:12,
        2:12;
        retain_xy = (4, 3),
        retain_xz = (4, 3),
        retain_yz = (4, 3),
        packet_kernel = :factorized_direct,
        term_coefficients = term_coefficients,
    )
    factorized = GaussletBases._nested_extract_factorized_basis(
        shell.coefficient_matrix,
        (13, 13, 13),
    )
    reconstructed = GaussletBases._nested_reconstruct_factorized_coefficients(factorized)
    @test factorized.reconstruction_max_error < 1.0e-10
    @test reconstructed ≈ shell.coefficient_matrix atol = 1.0e-10 rtol = 1.0e-10

    dims = (3, 2, 2)
    x1 = [1.0, 0.5, 0.0]
    x2 = [1.0, -0.25, 0.2]
    y1 = [1.0, 0.0]
    y2 = [1.0, 0.3]
    z1 = [1.0, -0.4]
    z2 = [1.0, 0.25]
    amplitudes = [2.0, -1.5, 0.75]
    factorable = zeros(Float64, prod(dims), 3)
    factors = ((x1, y1, z1), (x1, y2, z1), (x2, y1, z2))
    for column in 1:3
        xvec, yvec, zvec = factors[column]
        amplitude = amplitudes[column]
        for ix in 1:dims[1], iy in 1:dims[2], iz in 1:dims[3]
            factorable[GaussletBases._cartesian_flat_index(ix, iy, iz, dims), column] =
                amplitude * xvec[ix] * yvec[iy] * zvec[iz]
        end
    end
    hand_factorized = GaussletBases._nested_extract_factorized_basis(factorable, dims)
    @test hand_factorized.reconstruction_max_error < 1.0e-12
    @test hand_factorized.basis_triplets == [(1, 1, 1), (1, 2, 1), (2, 1, 2)]
    @test hand_factorized.basis_amplitudes ≈ amplitudes atol = 1.0e-12 rtol = 1.0e-12
    @test hand_factorized.x_functions[:, 1] ≈ x1 atol = 1.0e-12 rtol = 1.0e-12
    @test hand_factorized.x_functions[:, 2] ≈ x2 atol = 1.0e-12 rtol = 1.0e-12
    @test hand_factorized.y_functions[:, 1] ≈ y1 atol = 1.0e-12 rtol = 1.0e-12
    @test hand_factorized.y_functions[:, 2] ≈ y2 atol = 1.0e-12 rtol = 1.0e-12
    @test hand_factorized.z_functions[:, 1] ≈ z1 atol = 1.0e-12 rtol = 1.0e-12
    @test hand_factorized.z_functions[:, 2] ≈ z2 atol = 1.0e-12 rtol = 1.0e-12
    @test GaussletBases._nested_reconstruct_factorized_coefficients(hand_factorized) ≈
          factorable atol = 1.0e-12 rtol = 1.0e-12

    broken_factorable = copy(factorable)
    broken_factorable[GaussletBases._cartesian_flat_index(2, 2, 2, dims), 2] += 1.0e-4
    @test_throws ArgumentError GaussletBases._nested_extract_factorized_basis(
        broken_factorable,
        dims,
    )

    full_reference = GaussletBases._build_one_center_atomic_shell_sequence(
        bundle,
        (1:13, 1:13, 1:13);
        nside = 5,
        packet_kernel = :support_reference,
        term_coefficients = term_coefficients,
    )
    full_direct = GaussletBases._build_one_center_atomic_shell_sequence(
        bundle,
        (1:13, 1:13, 1:13);
        nside = 5,
        packet_kernel = :factorized_direct,
        term_coefficients = term_coefficients,
    )
    legacy_reference = GaussletBases._build_one_center_atomic_shell_sequence(
        bundle,
        (2:12, 2:12, 2:12);
        nside = 5,
        packet_kernel = :support_reference,
        term_coefficients = term_coefficients,
    )
    legacy_direct = GaussletBases._build_one_center_atomic_shell_sequence(
        bundle,
        (2:12, 2:12, 2:12);
        nside = 5,
        packet_kernel = :factorized_direct,
        term_coefficients = term_coefficients,
    )

    fixed_full_reference = GaussletBases._nested_fixed_block(full_reference, bundle)
    fixed_full_direct = GaussletBases._nested_fixed_block(full_direct, bundle)
    fixed_legacy_reference = GaussletBases._nested_fixed_block(legacy_reference, bundle)
    fixed_legacy_direct = GaussletBases._nested_fixed_block(legacy_direct, bundle)

    carried_legacy = fixed_legacy_reference.factorized_cartesian_parent_basis[]
    @test !isnothing(carried_legacy)
    extracted_legacy = GaussletBases._nested_extract_factorized_basis(
        fixed_legacy_reference.coefficient_matrix,
        (length(basis), length(basis), length(basis)),
    )
    @test GaussletBases._nested_factorized_parent_basis(fixed_legacy_reference) === carried_legacy
    @test carried_legacy.basis_triplets == extracted_legacy.basis_triplets
    @test carried_legacy.basis_amplitudes ≈
          extracted_legacy.basis_amplitudes atol = 1.0e-12 rtol = 1.0e-12
    @test GaussletBases._nested_reconstruct_factorized_coefficients(carried_legacy) ≈
          fixed_legacy_reference.coefficient_matrix atol = 1.0e-10 rtol = 1.0e-10
    fixed_legacy_reference.factorized_cartesian_parent_basis[] = nothing
    lazy_legacy = GaussletBases._nested_factorized_parent_basis(fixed_legacy_reference)
    @test fixed_legacy_reference.factorized_cartesian_parent_basis[] === lazy_legacy
    @test lazy_legacy.basis_triplets == carried_legacy.basis_triplets
    @test lazy_legacy.basis_amplitudes ≈
          carried_legacy.basis_amplitudes atol = 1.0e-12 rtol = 1.0e-12

    @test fixed_full_direct.overlap ≈ fixed_full_reference.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.kinetic ≈ fixed_full_reference.kinetic atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.position_x ≈ fixed_full_reference.position_x atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.position_y ≈ fixed_full_reference.position_y atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.position_z ≈ fixed_full_reference.position_z atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.x2_x ≈ fixed_full_reference.x2_x atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.x2_y ≈ fixed_full_reference.x2_y atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.x2_z ≈ fixed_full_reference.x2_z atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.weights ≈ fixed_full_reference.weights atol = 1.0e-10 rtol = 1.0e-10
    @test !hasproperty(fixed_full_direct, :gaussian_terms)
    @test !hasproperty(fixed_full_direct, :pair_terms)
    @test !hasproperty(fixed_full_direct, :term_storage)
    @test !hasproperty(fixed_full_reference, :gaussian_terms)
    @test !hasproperty(fixed_full_reference, :pair_terms)
    @test !hasproperty(fixed_full_reference, :term_storage)
    @test fixed_full_direct.gaussian_sum ≈ fixed_full_reference.gaussian_sum atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.pair_sum ≈ fixed_full_reference.pair_sum atol = 1.0e-10 rtol = 1.0e-10

    @test fixed_legacy_direct.overlap ≈ fixed_legacy_reference.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_legacy_direct.kinetic ≈ fixed_legacy_reference.kinetic atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_legacy_direct.position_x ≈ fixed_legacy_reference.position_x atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_legacy_direct.position_y ≈ fixed_legacy_reference.position_y atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_legacy_direct.position_z ≈ fixed_legacy_reference.position_z atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_legacy_direct.x2_x ≈ fixed_legacy_reference.x2_x atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_legacy_direct.x2_y ≈ fixed_legacy_reference.x2_y atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_legacy_direct.x2_z ≈ fixed_legacy_reference.x2_z atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_legacy_direct.weights ≈ fixed_legacy_reference.weights atol = 1.0e-10 rtol = 1.0e-10
    @test !hasproperty(fixed_legacy_direct, :gaussian_terms)
    @test !hasproperty(fixed_legacy_direct, :pair_terms)
    @test !hasproperty(fixed_legacy_direct, :term_storage)
    @test !hasproperty(fixed_legacy_reference, :gaussian_terms)
    @test !hasproperty(fixed_legacy_reference, :pair_terms)
    @test !hasproperty(fixed_legacy_reference, :term_storage)
    @test fixed_legacy_direct.gaussian_sum ≈ fixed_legacy_reference.gaussian_sum atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_legacy_direct.pair_sum ≈ fixed_legacy_reference.pair_sum atol = 1.0e-10 rtol = 1.0e-10
end

@testset "QW residual-space keep policy is near-null-only and stabilized" begin
    # Literal residual-overlap spectrum observed on the anchored one-center
    # Ne legacy-profile case:
    # parent side = 29, working box = 2:28, nside = 7, supplement lmax = 1.
    residual_overlap_eigenvalues = Float64[
        6.486197469e-08,
        3.165964397e-06,
        3.165964398e-06,
        3.165964398e-06,
        5.681904400e-06,
        1.681965647e-05,
        3.337404514e-05,
        5.805312472e-05,
        5.805312472e-05,
        5.805312472e-05,
        7.256172691e-05,
        1.406818079e-04,
        1.406818079e-04,
        1.406818079e-04,
        1.927015773e-04,
        1.927015773e-04,
        1.927015773e-04,
        1.995510583e-04,
        4.261857498e-04,
        4.261857498e-04,
        4.261857498e-04,
        5.359433116e-04,
        1.945893481e-03,
        1.945893481e-03,
        1.945893481e-03,
    ]
    gausslet_overlap = Matrix{Float64}(I, 1, 1)
    overlap_ga = zeros(Float64, 1, length(residual_overlap_eigenvalues))
    overlap_aa = Matrix(Diagonal(residual_overlap_eigenvalues))
    near_null_diagnostics = diagnose_qwrg_residual_space(
        gausslet_overlap,
        overlap_ga,
        overlap_aa;
        keep_policy = :near_null_only,
        keep_abs_tol = GaussletBases._qwrg_atomic_residual_keep_tol(),
        accept_tol = GaussletBases._qwrg_atomic_residual_accept_tol(),
    )
    near_null_data = GaussletBases._qwrg_residual_space(
        gausslet_overlap,
        overlap_ga,
        overlap_aa;
        keep_policy = :near_null_only,
        keep_abs_tol = GaussletBases._qwrg_atomic_residual_keep_tol(),
        accept_tol = GaussletBases._qwrg_atomic_residual_accept_tol(),
    )
    legacy_alias_diagnostics = diagnose_qwrg_residual_space(
        gausslet_overlap,
        overlap_ga,
        overlap_aa;
        keep_policy = :legacy_profile,
        keep_abs_tol = GaussletBases._qwrg_atomic_residual_keep_tol(),
        accept_tol = GaussletBases._qwrg_atomic_residual_accept_tol(),
    )

    @test near_null_diagnostics.keep_policy == :near_null_only
    @test near_null_diagnostics.gaussian_count == 25
    @test near_null_diagnostics.supplement_numerical_rank == 25
    @test near_null_diagnostics.residual_numerical_rank == 25
    @test near_null_diagnostics.kept_count == 24
    @test near_null_diagnostics.discarded_count == 1
    @test near_null_diagnostics.keep_tol ≈ 1.0e-7 atol = 0.0 rtol = 0.0
    @test near_null_diagnostics.accept_tol ≈ 1.0e-7 atol = 0.0 rtol = 0.0
    @test near_null_diagnostics.kept_block_stabilization_null_tol ≈ 1.0e-12 atol = 1.0e-15 rtol = 0.0
    @test near_null_diagnostics.kept_block_stabilization_correction_passes >= 1
    @test near_null_diagnostics.kept_block_stabilization_clipped_count == 0
    @test near_null_diagnostics.kept_block_stabilization_dropped_count == 0
    @test near_null_diagnostics.kept_block_pre_stabilization_overlap_error < 1.0e-12
    @test near_null_diagnostics.kept_block_post_stabilization_overlap_error < 1.0e-12
    @test near_null_diagnostics.kept_block_pre_stabilization_symmetry_defect < 1.0e-12
    @test near_null_diagnostics.kept_block_post_stabilization_symmetry_defect < 1.0e-12
    @test near_null_diagnostics.kept_block_pre_stabilization_negative_count == 0
    @test near_null_diagnostics.kept_block_post_stabilization_negative_count == 0
    @test near_null_diagnostics.kept_block_pre_stabilization_near_null_count == 0
    @test near_null_diagnostics.kept_block_post_stabilization_near_null_count == 0
    @test norm(near_null_data.final_overlap - I, Inf) < 1.0e-10
    @test legacy_alias_diagnostics.keep_policy == :near_null_only
    @test legacy_alias_diagnostics.kept_count == near_null_diagnostics.kept_count
    @test legacy_alias_diagnostics.keep_tol == near_null_diagnostics.keep_tol
    @test legacy_alias_diagnostics.kept_block_post_stabilization_overlap_error ==
        near_null_diagnostics.kept_block_post_stabilization_overlap_error

    nsynthetic = 69
    synthetic_raw_overlap = Matrix{Float64}(I, nsynthetic, nsynthetic)
    synthetic_coefficients = Matrix{Float64}(I, nsynthetic, nsynthetic)
    @inbounds for i in 1:nsynthetic, j in 1:nsynthetic
        synthetic_coefficients[i, j] += 8.0e-9 * sin(Float64(i + 2 * j))
    end
    synthetic_stabilization = GaussletBases._qwrg_stabilize_residual_coefficients(
        synthetic_raw_overlap,
        synthetic_coefficients,
    )
    @test synthetic_stabilization.pre_error > 1.0e-8
    @test synthetic_stabilization.post_error < 1.0e-10
    @test synthetic_stabilization.post_symmetry_defect < 1.0e-12
    @test synthetic_stabilization.pre_negative_count == 0
    @test synthetic_stabilization.post_negative_count == 0
    @test synthetic_stabilization.dropped_count == 0
    @test synthetic_stabilization.correction_passes >= 1
end

@testset "One-center atomic legacy-profile residual completion contract" begin
    if !_RUN_SLOW_TESTS
        @test true
    else
        data = _one_center_atomic_legacy_profile_ne_residual_completion_fixture()

        @test data.fixed_gausslet_count == 2523
        @test data.supplement_count == 25

        @test data.near_null.keep_policy == :near_null_only
        @test data.near_null.residual_numerical_rank == 25
        @test data.near_null.kept_count == 24
        @test data.near_null.discarded_count == 1
        @test data.near_null.keep_tol ≈ 1.0e-7 atol = 0.0 rtol = 0.0
        @test data.near_null.accept_tol ≈ 1.0e-7 atol = 0.0 rtol = 0.0
        @test data.near_null.kept_block_pre_stabilization_overlap_error > 0.0
        @test data.near_null.kept_block_post_stabilization_overlap_error <
            data.near_null.kept_block_pre_stabilization_overlap_error
        @test data.near_null.kept_block_post_stabilization_overlap_error < 1.0e-9
        @test data.near_null.kept_block_post_stabilization_symmetry_defect < 1.0e-9
        @test data.near_null.kept_block_stabilization_dropped_count == 0
        @test norm(data.near_null_data.final_overlap - I, Inf) < 1.0e-7
        @test data.near_null_total_basis == 2547
        @test data.legacy_alias.keep_policy == :near_null_only
        @test data.legacy_alias.kept_count == data.near_null.kept_count
    end
end

@testset "Atomic residual keep policy rejects relative_case_scale on public QW routes" begin
    if !_legacy_basisfile_available()
        @test true
    else
        source_basis_qw, _legacy_qw, _ordinary_l0, _ordinary_l0_check = _qiu_white_full_nearest_fixture()
        supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 1)
        err = try
            ordinary_cartesian_qiu_white_operators(
                source_basis_qw,
                supplement;
                expansion = coulomb_gaussian_expansion(doacc = false),
                Z = 2.0,
                interaction_treatment = :ggt_nearest,
                residual_keep_policy = :relative_case_scale,
            )
            nothing
        catch caught
            caught
        end
        @test err isa ArgumentError
        @test occursin(":near_null_only", sprint(showerror, err))
        @test occursin(":legacy_profile", sprint(showerror, err))
    end
end

@testset "One-center atomic ns=9 legacy-profile residual stabilization closes center-extraction failure" begin
    if !_RUN_SLOW_TESTS
        @test true
    else
        data = _one_center_atomic_ns9_legacy_profile_qw_fixture()
        @test data.residual_data.diagnostics.kept_count == 24
        @test data.residual_data.diagnostics.keep_policy == :near_null_only
        @test data.residual_data.diagnostics.keep_tol ≈ 1.0e-7 atol = 0.0 rtol = 0.0
        @test data.residual_data.diagnostics.accept_tol ≈ 1.0e-7 atol = 0.0 rtol = 0.0
        @test data.residual_data.diagnostics.kept_block_pre_stabilization_overlap_error > 0.0
        @test data.residual_data.diagnostics.kept_block_post_stabilization_overlap_error <=
            data.residual_data.diagnostics.kept_block_pre_stabilization_overlap_error
        @test data.residual_data.diagnostics.kept_block_post_stabilization_overlap_error < 1.0e-9
        @test data.residual_data.diagnostics.kept_block_post_stabilization_symmetry_defect < 1.0e-9
        @test data.residual_data.diagnostics.kept_block_stabilization_dropped_count == 0
        @test norm(data.residual_data.final_overlap - I, Inf) < 1.0e-7
        @test data.operators.residual_count == 24
        @test norm(data.operators.overlap - I, Inf) < 1.0e-7
        check = GaussletBases.ordinary_cartesian_1s2_check(data.operators)
        @test isfinite(check.orbital_energy)
        @test check.overlap_error < 1.0e-7
    end
end

@testset "Cartesian nested shell sequence fixed-block" begin
    (
        basis,
        bundle,
        shell1,
        shell2,
        shell_plus_core,
        shell_sequence,
        fixed_shell_plus_core,
        fixed_sequence,
        legacy,
        baseline,
        shell_plus_core_ops,
        shell_sequence_ops,
        baseline_check,
        shell_plus_core_check,
        shell_sequence_check,
    ) = _nested_qiu_white_shell_sequence_fixture()

    @test shell_sequence isa GaussletBases._CartesianNestedShellSequence3D
    @test length(shell_sequence.shell_layers) == 2
    @test shell_sequence.shell_layers[1] === shell1
    @test shell_sequence.shell_layers[2] === shell2
    @test first(shell_sequence.core_column_range) == 1
    @test last(shell_sequence.core_column_range) == length(shell_sequence.core_indices)
    @test length(shell_sequence.layer_column_ranges) == 2
    @test length(shell_sequence.core_indices) == 11^3
    @test isempty(intersect(shell_sequence.core_indices, shell1.support_indices))
    @test isempty(intersect(shell_sequence.core_indices, shell2.support_indices))
    @test isempty(intersect(shell1.support_indices, shell2.support_indices))

    @test fixed_sequence isa GaussletBases._NestedFixedBlock3D
    @test fixed_sequence.parent_basis === basis
    @test fixed_sequence.shell === shell_sequence
    @test shell_plus_core_ops.gausslet_count == 1385
    @test shell_sequence_ops.gausslet_count == 1439
    @test baseline.gausslet_count == 17^3
    @test norm(fixed_shell_plus_core.overlap - I, Inf) < 1.0e-10
    @test norm(fixed_sequence.overlap - I, Inf) < 1.0e-10
    @test shell_plus_core_check.overlap_error < 1.0e-10
    @test shell_sequence_check.overlap_error < 1.0e-10
    @test shell_plus_core_check.orbital_energy < 0.0
    @test shell_sequence_check.orbital_energy < 0.0
    @test shell_plus_core_check.vee_expectation > 0.0
    @test shell_sequence_check.vee_expectation > 0.0
    @test abs(shell_plus_core_check.orbital_energy - baseline_check.orbital_energy) < 1.0e-4
    @test abs(shell_plus_core_check.vee_expectation - baseline_check.vee_expectation) < 1.0e-4
    @test abs(shell_sequence_check.orbital_energy - baseline_check.orbital_energy) < 1.0e-4
    @test abs(shell_sequence_check.vee_expectation - baseline_check.vee_expectation) < 1.0e-4
    @test abs(shell_sequence_check.orbital_energy - shell_plus_core_check.orbital_energy) < 1.0e-4
    @test abs(shell_sequence_check.vee_expectation - shell_plus_core_check.vee_expectation) < 1.0e-4
end

@testset "Cartesian nested fixed-nside compression policy" begin
    (
        basis,
        bundle,
        shell1,
        shell2,
        shell3,
        grow_sequence,
        shrinking_sequence,
        fixed_grow,
        fixed_shrink,
        legacy,
        baseline,
        grow_ops,
        shrink_ops,
        baseline_check,
        grow_check,
        shrink_check,
    ) = _nested_qiu_white_nside_sequence_fixture()

    @test GaussletBases._nested_shrunk_interval(4:14, 0; nside = 5) == 4:14
    @test GaussletBases._nested_shrunk_interval(4:14, 1; nside = 5) == 5:13
    @test GaussletBases._nested_shrunk_interval(4:14, 2; nside = 5) == 6:12
    @test GaussletBases._nested_shrunk_interval(4:14, 3; nside = 5) == 7:11
    @test GaussletBases._nested_shrunk_interval(4:14, 4; nside = 5) == 7:11

    @test shrinking_sequence isa GaussletBases._CartesianNestedShellSequence3D
    @test length(shrinking_sequence.shell_layers) == 3
    @test shrinking_sequence.shell_layers[1] === shell1
    @test shrinking_sequence.shell_layers[2] === shell2
    @test shrinking_sequence.shell_layers[3] === shell3
    @test length(grow_sequence.core_indices) == 5^3
    @test length(shrinking_sequence.core_indices) == 5^3
    @test first(shrinking_sequence.core_column_range) == 1
    @test last(shrinking_sequence.core_column_range) == 3^3
    @test isempty(intersect(shrinking_sequence.core_indices, shell1.support_indices))
    @test isempty(intersect(shrinking_sequence.core_indices, shell2.support_indices))
    @test isempty(intersect(shrinking_sequence.core_indices, shell3.support_indices))
    @test isempty(intersect(shell1.support_indices, shell2.support_indices))
    @test isempty(intersect(shell1.support_indices, shell3.support_indices))
    @test isempty(intersect(shell2.support_indices, shell3.support_indices))

    @test fixed_grow isa GaussletBases._NestedFixedBlock3D
    @test fixed_shrink isa GaussletBases._NestedFixedBlock3D
    @test grow_ops.gausslet_count == 287
    @test shrink_ops.gausslet_count == 189
    @test baseline.gausslet_count == 17^3
    @test norm(fixed_grow.overlap - I, Inf) < 1.0e-10
    @test norm(fixed_shrink.overlap - I, Inf) < 1.0e-10
    @test grow_check.overlap_error < 1.0e-10
    @test shrink_check.overlap_error < 1.0e-10
    @test isfinite(grow_check.orbital_energy)
    @test isfinite(grow_check.vee_expectation)
    @test isfinite(shrink_check.orbital_energy)
    @test isfinite(shrink_check.vee_expectation)
    @test grow_check.orbital_energy < 0.0
    @test shrink_check.orbital_energy < 0.0
    @test grow_check.vee_expectation > 0.0
    @test shrink_check.vee_expectation > 0.0
    @test shrink_ops.gausslet_count < grow_ops.gausslet_count
end

@testset "Cartesian nested complete shell layer" begin
    (
        basis,
        bundle,
        shell1_complete,
        shell2_complete,
        shell3_complete,
        shell4_complete,
        interval1,
        interval2,
        interval3,
        interval4,
        core5,
        complete_sequence,
        fixed_complete_sequence,
        legacy,
        baseline,
        complete_sequence_ops,
        baseline_check,
        complete_sequence_check,
    ) = _nested_qiu_white_complete_shell_sequence_fixture()

    @test shell1_complete isa GaussletBases._CartesianNestedCompleteShell3D
    @test shell2_complete isa GaussletBases._CartesianNestedCompleteShell3D
    @test shell3_complete isa GaussletBases._CartesianNestedCompleteShell3D
    @test shell4_complete isa GaussletBases._CartesianNestedCompleteShell3D
    @test length(shell1_complete.faces) == 6
    @test length(shell1_complete.edges) == 12
    @test length(shell1_complete.corners) == 8
    @test length(shell2_complete.faces) == 6
    @test length(shell2_complete.edges) == 12
    @test length(shell2_complete.corners) == 8
    @test length(shell3_complete.faces) == 6
    @test length(shell3_complete.edges) == 12
    @test length(shell3_complete.corners) == 8
    @test length(shell4_complete.faces) == 6
    @test length(shell4_complete.edges) == 12
    @test length(shell4_complete.corners) == 8

    @test length(shell1_complete.support_indices) == 13^3 - 11^3
    @test length(shell2_complete.support_indices) == 11^3 - 9^3
    @test length(shell3_complete.support_indices) == 9^3 - 7^3
    @test length(shell4_complete.support_indices) == 7^3 - 5^3
    @test shell1_complete.provenance.source_box == (
        (first(interval1) - 1):(last(interval1) + 1),
        (first(interval1) - 1):(last(interval1) + 1),
        (first(interval1) - 1):(last(interval1) + 1),
    )
    @test shell1_complete.provenance.next_inner_box == (interval1, interval1, interval1)
    @test shell1_complete.provenance.source_point_count == 13^3 - 11^3
    @test shell1_complete.provenance.retained_fixed_count == size(shell1_complete.coefficient_matrix, 2)
    @test shell2_complete.provenance.source_box == (
        (first(interval2) - 1):(last(interval2) + 1),
        (first(interval2) - 1):(last(interval2) + 1),
        (first(interval2) - 1):(last(interval2) + 1),
    )
    @test shell2_complete.provenance.next_inner_box == (interval2, interval2, interval2)
    @test shell2_complete.provenance.source_point_count == 11^3 - 9^3
    @test shell2_complete.provenance.retained_fixed_count == size(shell2_complete.coefficient_matrix, 2)
    @test shell3_complete.provenance.source_box == (
        (first(interval3) - 1):(last(interval3) + 1),
        (first(interval3) - 1):(last(interval3) + 1),
        (first(interval3) - 1):(last(interval3) + 1),
    )
    @test shell3_complete.provenance.next_inner_box == (interval3, interval3, interval3)
    @test shell3_complete.provenance.source_point_count == 9^3 - 7^3
    @test shell3_complete.provenance.retained_fixed_count == size(shell3_complete.coefficient_matrix, 2)
    @test shell4_complete.provenance.source_box == (
        (first(interval4) - 1):(last(interval4) + 1),
        (first(interval4) - 1):(last(interval4) + 1),
        (first(interval4) - 1):(last(interval4) + 1),
    )
    @test shell4_complete.provenance.next_inner_box == (interval4, interval4, interval4)
    @test shell4_complete.provenance.source_point_count == 7^3 - 5^3
    @test shell4_complete.provenance.retained_fixed_count == size(shell4_complete.coefficient_matrix, 2)
    @test sum(length(face.support_indices) for face in shell1_complete.faces) == 6 * 11^2
    @test sum(length(edge.support_indices) for edge in shell1_complete.edges) == 12 * 11
    @test sum(length(corner.support_indices) for corner in shell1_complete.corners) == 8

    @test complete_sequence isa GaussletBases._CartesianNestedShellSequence3D
    @test length(complete_sequence.core_indices) == 5^3
    @test complete_sequence.working_box == (3:15, 3:15, 3:15)
    @test baseline.gausslet_count == 17^3
    @test norm(fixed_complete_sequence.overlap - I, Inf) < 1.0e-10
    @test all(isfinite, fixed_complete_sequence.weights)
    @test minimum(fixed_complete_sequence.weights) > 0.0
    @test complete_sequence_check.overlap_error < 1.0e-10
    @test isfinite(complete_sequence_check.orbital_energy)
    @test isfinite(complete_sequence_check.vee_expectation)
    # Post-hardening residual-space route check for the complete-shell candidate only.
    @test abs(complete_sequence_check.vee_expectation - baseline_check.vee_expectation) < 3.0e-4

    expansion = coulomb_gaussian_expansion(doacc = false)
    overlap_parent, one_body_parent, interaction_parent = _nested_parent_fixed_problem(bundle, expansion; Z = 2.0)
    parent_modes = eigen(Hermitian(one_body_parent), Hermitian(overlap_parent))
    parent_ground = parent_modes.vectors[:, 1]
    parent_ground_vee = _nested_vee_from_orbital(interaction_parent, parent_ground)
    projected_complete = _nested_fixed_projected_orbital(overlap_parent, fixed_complete_sequence, parent_ground)
    projected_complete_vee = _nested_vee_from_orbital(
        GaussletBases._qwrg_fixed_block_interaction_matrix(fixed_complete_sequence, expansion),
        projected_complete,
    )
    term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
    @test abs(projected_complete_vee - parent_ground_vee) < 5.0e-4
    @test_throws ArgumentError GaussletBases._nested_shell_sequence(
        bundle,
        core5,
        core5,
        core5,
        [shell1_complete, shell2_complete, shell3_complete],
        term_coefficients = term_coefficients,
    )
end
