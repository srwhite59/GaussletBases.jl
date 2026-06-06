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

include("bond_aligned_diatomic_atom_growth_anatomy_policy_runtests.jl")
include("bond_aligned_diatomic_atom_growth_construction_plan_runtests.jl")
include("bond_aligned_diatomic_high_order_recipe_policy_metadata_runtests.jl")
include("bond_aligned_diatomic_high_order_recipe_realization_audit_runtests.jl")

include("bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl")

include("bond_aligned_diatomic_endcap_panel_shared_shell_source_policy_runtests.jl")

include("cartesian_nested_shell_first_packet_runtests.jl")
include("cartesian_nested_support_immediate_contraction_helpers_runtests.jl")
include("cartesian_nested_shell_interface_runtests.jl")
include("cartesian_nested_fixed_block_qw_pgdg_adapter_runtests.jl")

include("one_center_atomic_full_parent_nested_contract_runtests.jl")
include("one_center_atomic_legacy_profile_nested_contract_runtests.jl")
include("one_center_atomic_fixed_block_timing_surface_runtests.jl")
include("global_timing_macro_surface_runtests.jl")
include("one_center_atomic_compact_fixed_block_term_storage_runtests.jl")

include("cartesian_basis_representation_direct_product_qw_bases_runtests.jl")

include("cartesian_parent_gausslet_basis_identity_runtests.jl")

include("cartesian_contracted_parent_scaffold_runtests.jl")

include("cartesian_contracted_parent_metric_packet_runtests.jl")

include("cartesian_basis_representation_nested_fixed_blocks_runtests.jl")

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

include("nested_coefficient_maps_sparse_storage_runtests.jl")

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
