@testset "Ordinary QW internal hydrogenic/ESOI corrections" begin
    Z = 2.0
    basis = build_basis(MappedUniformBasisSpec(
        :G10;
        count = 9,
        mapping = white_lindsey_atomic_mapping(Z = Z, d = 0.2, tail_spacing = 10.0),
        reference_spacing = 1.0,
    ))
    supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 0)
    operators = ordinary_cartesian_qiu_white_operators(
        basis,
        supplement;
        Z = Z,
        interaction_treatment = :ggt_nearest,
        residual_keep_policy = :near_null_only,
    )
    target_one_body = -0.5 * Z^2
    target_closed_shell = -Z^2 + 5.0 * Z / 8.0

    @test_throws ArgumentError GaussletBases.HydrogenicCoreCorrectionSpec(;
        Z = Z,
        one_body_mode = :local_first_order,
        two_body_mode = :esoi_local,
    )

    for one_body_mode in (:projector, :local_exact)
        spec = GaussletBases.HydrogenicCoreCorrectionSpec(;
            Z = Z,
            one_body_mode = one_body_mode,
            two_body_mode = :esoi_local,
        )
        result = GaussletBases._apply_ordinary_cartesian_corrections(operators, spec)
        corrected = result.operators
        diagnostics = result.diagnostics
        check = GaussletBases.ordinary_cartesian_1s2_check(corrected)

        @test result isa GaussletBases.OrdinaryCartesianCorrectionResult
        @test diagnostics.one_body_mode == one_body_mode
        @test diagnostics.two_body_mode == :esoi_local
        @test diagnostics.initial_lowest_core_eigenvalue > target_one_body
        @test diagnostics.corrected_lowest_core_eigenvalue ≈ target_one_body atol = 1.0e-9 rtol = 0.0
        @test diagnostics.corrected_1s_coulomb ≈ 5.0 * Z / 8.0 atol = 1.0e-10 rtol = 0.0
        @test diagnostics.closed_shell_corrected_energy ≈ target_closed_shell atol = 1.0e-8 rtol = 0.0
        @test check.orbital_energy ≈ target_one_body atol = 1.0e-9 rtol = 0.0
        @test 2.0 * check.orbital_energy + check.vee_expectation ≈ target_closed_shell atol = 1.0e-8 rtol = 0.0
        @test corrected.kinetic_one_body === nothing
        @test corrected.nuclear_one_body_by_center === nothing
        @test corrected.nuclear_term_storage == :total_only
    end
end

@testset "Ordinary QW public hydrogenic projector corrections" begin
    Z = 2.0
    basis = build_basis(MappedUniformBasisSpec(
        :G10;
        count = 9,
        mapping = white_lindsey_atomic_mapping(Z = Z, d = 0.2, tail_spacing = 10.0),
        reference_spacing = 1.0,
    ))
    supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 0)
    operators = ordinary_cartesian_qiu_white_operators(
        basis,
        supplement;
        Z = Z,
        interaction_treatment = :ggt_nearest,
        residual_keep_policy = :near_null_only,
    )
    target_one_body = -0.5 * Z^2
    target_closed_shell = -Z^2 + 5.0 * Z / 8.0

    projector_only = apply_ordinary_cartesian_corrections(
        operators,
        HydrogenicCoreProjectorCorrectionSpec(; Z = Z),
    )
    projector_esoi = apply_ordinary_cartesian_corrections(
        operators,
        HydrogenicCoreProjectorCorrectionSpec(; Z = Z, include_esoi = true),
    )
    keyword_result = apply_ordinary_cartesian_corrections(
        operators;
        Z = Z,
        include_esoi = true,
    )

    @test projector_only isa OrdinaryCartesianCorrectionResult
    @test projector_only.diagnostics.one_body_mode == :projector
    @test projector_only.diagnostics.two_body_mode == :none
    @test projector_only.diagnostics.corrected_lowest_core_eigenvalue ≈ target_one_body atol = 1.0e-9 rtol = 0.0
    @test projector_only.diagnostics.closed_shell_corrected_energy - target_closed_shell ≈
          projector_only.diagnostics.corrected_1s_coulomb - 5.0 * Z / 8.0 atol = 1.0e-8 rtol = 0.0

    @test projector_esoi.diagnostics.one_body_mode == :projector
    @test projector_esoi.diagnostics.two_body_mode == :esoi_local
    @test projector_esoi.diagnostics.corrected_lowest_core_eigenvalue ≈ target_one_body atol = 1.0e-9 rtol = 0.0
    @test projector_esoi.diagnostics.corrected_1s_coulomb ≈ 5.0 * Z / 8.0 atol = 1.0e-10 rtol = 0.0
    @test projector_esoi.diagnostics.closed_shell_corrected_energy ≈ target_closed_shell atol = 1.0e-8 rtol = 0.0
    @test keyword_result.diagnostics.closed_shell_corrected_energy ≈
          projector_esoi.diagnostics.closed_shell_corrected_energy atol = 1.0e-12 rtol = 0.0

    for result in (projector_only, projector_esoi)
        @test result.operators.kinetic_one_body === nothing
        @test result.operators.nuclear_one_body_by_center === nothing
        @test result.operators.nuclear_term_storage == :total_only
    end
end

include(joinpath(@__DIR__, "high_order_doside_experimental_runtests.jl"))

@testset "Bond-aligned homonuclear chain ordinary QW reference path" begin
    basis2, operators2, diagnostics2 = _bond_aligned_homonuclear_chain_qw_fixture(; natoms = 2, spacing = 1.4)
    basis3, operators3, diagnostics3 = _bond_aligned_homonuclear_chain_qw_fixture(; natoms = 3, spacing = 1.2)

    @test basis2 isa BondAlignedHomonuclearChainQWBasis3D
    @test basis3 isa BondAlignedHomonuclearChainQWBasis3D
    @test operators2 isa QiuWhiteResidualGaussianOperators
    @test operators3 isa QiuWhiteResidualGaussianOperators
    @test operators2.gaussian_data === nothing
    @test operators3.gaussian_data === nothing
    @test operators2.residual_count == 0
    @test operators3.residual_count == 0
    @test operators2.gausslet_count == length(basis2.basis_x) * length(basis2.basis_y) * length(basis2.basis_z)
    @test operators3.gausslet_count == length(basis3.basis_x) * length(basis3.basis_y) * length(basis3.basis_z)
    @test size(operators2.overlap, 1) < 1000
    @test size(operators3.overlap, 1) < 1000
    @test norm(operators2.overlap - I, Inf) < 1.0e-8
    @test norm(operators3.overlap - I, Inf) < 1.0e-8
    @test operators2.one_body_hamiltonian ≈ transpose(operators2.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test operators3.one_body_hamiltonian ≈ transpose(operators3.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test operators2.interaction_matrix ≈ transpose(operators2.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test operators3.interaction_matrix ≈ transpose(operators3.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test all(isfinite, operators2.one_body_hamiltonian)
    @test all(isfinite, operators3.one_body_hamiltonian)
    @test all(isfinite, operators2.interaction_matrix)
    @test all(isfinite, operators3.interaction_matrix)
    @test minimum(diag(operators2.interaction_matrix)) > 0.0
    @test minimum(diag(operators3.interaction_matrix)) > 0.0
    @test diagnostics2.axis_monotone
    @test diagnostics3.axis_monotone
    @test all(diagnostics2.local_spacings_at_midpoints .> 0.45)
    @test all(diagnostics3.local_spacings_at_midpoints .> 0.45)
    @test length(basis2.basis_z) >= length(basis2.basis_x)
    @test length(basis3.basis_z) > length(basis3.basis_x)
end

@testset "Axis-aligned homonuclear square-lattice ordinary QW reference path" begin
    basis2, operators2, diagnostics2, check2 =
        _axis_aligned_homonuclear_square_lattice_qw_fixture(; n = 2, spacing = 1.4)
    basis3, operators3, diagnostics3, check3 =
        _axis_aligned_homonuclear_square_lattice_qw_fixture(; n = 3, spacing = 1.2)

    @test basis2 isa AxisAlignedHomonuclearSquareLatticeQWBasis3D
    @test basis3 isa AxisAlignedHomonuclearSquareLatticeQWBasis3D
    @test operators2 isa QiuWhiteResidualGaussianOperators
    @test operators3 isa QiuWhiteResidualGaussianOperators
    @test operators2.gaussian_data === nothing
    @test operators3.gaussian_data === nothing
    @test operators2.residual_count == 0
    @test operators3.residual_count == 0
    @test operators2.gausslet_count == length(basis2.basis_x) * length(basis2.basis_y) * length(basis2.basis_z)
    @test operators3.gausslet_count == length(basis3.basis_x) * length(basis3.basis_y) * length(basis3.basis_z)
    @test size(operators2.overlap, 1) < 1500
    @test size(operators3.overlap, 1) < 4000
    @test norm(operators2.overlap - I, Inf) < 1.0e-8
    @test norm(operators3.overlap - I, Inf) < 1.0e-8
    @test operators2.one_body_hamiltonian ≈ transpose(operators2.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test operators3.one_body_hamiltonian ≈ transpose(operators3.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test operators2.interaction_matrix ≈ transpose(operators2.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test operators3.interaction_matrix ≈ transpose(operators3.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test all(isfinite, operators2.one_body_hamiltonian)
    @test all(isfinite, operators3.one_body_hamiltonian)
    @test all(isfinite, operators2.interaction_matrix)
    @test all(isfinite, operators3.interaction_matrix)
    @test minimum(diag(operators2.interaction_matrix)) > 0.0
    @test minimum(diag(operators3.interaction_matrix)) > 0.0
    @test diagnostics2.x_axis_monotone
    @test diagnostics2.y_axis_monotone
    @test diagnostics3.x_axis_monotone
    @test diagnostics3.y_axis_monotone
    @test diagnostics2.xy_axis_center_match_error < 1.0e-12
    @test diagnostics3.xy_axis_center_match_error < 1.0e-12
    @test diagnostics2.local_spacing_at_plane_center_x > 0.45
    @test diagnostics2.local_spacing_at_plane_center_y > 0.45
    @test all(diagnostics3.representative_midpoint_spacings_x .> 0.45)
    @test all(diagnostics3.representative_midpoint_spacings_y .> 0.45)
    @test isfinite(check2.orbital_energy)
    @test isfinite(check3.orbital_energy)
    @test isfinite(check2.vee_expectation)
    @test isfinite(check3.vee_expectation)
end

@testset "Axis-aligned homonuclear square-lattice nested geometry diagnostics" begin
    basis2, source2, fixed2, diagnostics2 =
        _axis_aligned_homonuclear_square_lattice_nested_fixture(; n = 2, spacing = 1.4)
    basis3, source3, fixed3, diagnostics3 =
        _axis_aligned_homonuclear_square_lattice_nested_fixture(; n = 3, spacing = 1.2)

    @test basis2 isa AxisAlignedHomonuclearSquareLatticeQWBasis3D
    @test basis3 isa AxisAlignedHomonuclearSquareLatticeQWBasis3D
    @test diagnostics2.retention_contract.retain_xy == (3, 3)
    @test diagnostics3.retention_contract.retain_xy == (3, 3)
    @test diagnostics2.retention_contract.retain_x_edge == 3
    @test diagnostics3.retention_contract.retain_x_edge == 3
    @test diagnostics2.retention_contract.shell_increment == 98
    @test diagnostics3.retention_contract.shell_increment == 98
    @test diagnostics2.retention_contract.matches_nside_default
    @test diagnostics3.retention_contract.matches_nside_default
    @test diagnostics2.contract_audit.full_parent_working_box
    @test diagnostics3.contract_audit.full_parent_working_box
    @test diagnostics2.contract_audit.missing_row_count == 0
    @test diagnostics3.contract_audit.missing_row_count == 0
    @test diagnostics2.contract_audit.ownership_unowned_row_count == 0
    @test diagnostics3.contract_audit.ownership_unowned_row_count == 0
    @test diagnostics2.contract_audit.ownership_multi_owned_row_count == 0
    @test diagnostics3.contract_audit.ownership_multi_owned_row_count == 0
    @test diagnostics2.shared_shells_match_contract
    @test diagnostics3.shared_shells_match_contract
    @test Int[shell.retained_fixed_count for shell in diagnostics2.shared_shell_provenance] ==
          diagnostics2.shared_shell_dimensions
    @test Int[shell.retained_fixed_count for shell in diagnostics3.shared_shell_provenance] ==
          diagnostics3.shared_shell_dimensions
    @test all(
        shell.source_point_count ==
        prod(length.(shell.source_box)) - prod(length.(shell.next_inner_box)) for
        shell in diagnostics2.shared_shell_provenance
    )
    @test all(
        shell.source_point_count ==
        prod(length.(shell.source_box)) - prod(length.(shell.next_inner_box)) for
        shell in diagnostics3.shared_shell_provenance
    )
    @test diagnostics2.leaf_count == 4
    @test diagnostics3.leaf_count == 9
    @test diagnostics2.root_node.did_split
    @test diagnostics3.root_node.did_split
    @test diagnostics2.root_node.accepted_candidate_index == 1
    @test diagnostics3.root_node.accepted_candidate_index == 1
    @test diagnostics2.root_node.child_count == 2
    @test diagnostics3.root_node.child_count == 3
    @test !diagnostics2.root_node.local_resolution_warning
    @test !diagnostics3.root_node.local_resolution_warning
    @test diagnostics2.root_node.candidate_summaries[1].split_family == :split_x_binary
    @test diagnostics2.root_node.candidate_summaries[2].split_family == :split_y_binary
    @test diagnostics3.root_node.candidate_summaries[1].split_family == :split_x_ternary
    @test diagnostics3.root_node.candidate_summaries[2].split_family == :split_y_ternary
    @test diagnostics2.root_node.candidate_summaries[1].symmetry_preserving
    @test diagnostics2.root_node.candidate_summaries[2].symmetry_preserving
    @test diagnostics3.root_node.candidate_summaries[1].symmetry_preserving
    @test diagnostics3.root_node.candidate_summaries[2].symmetry_preserving
    @test all(
        ratio >= diagnostics2.root_node.min_in_plane_aspect_ratio for
        ratio in diagnostics2.root_node.candidate_summaries[1].child_in_plane_aspect_ratios
    )
    @test all(
        ratio >= diagnostics3.root_node.min_in_plane_aspect_ratio for
        ratio in diagnostics3.root_node.candidate_summaries[1].child_in_plane_aspect_ratios
    )
    @test length(source2.root_geometry.child_nodes) == 2
    @test length(source3.root_geometry.child_nodes) == 3
    @test source2.root_geometry.child_nodes[1].accepted_candidate_index == 1
    @test source2.root_geometry.child_nodes[2].accepted_candidate_index == 1
    @test all(child.accepted_candidate_index == 1 for child in source3.root_geometry.child_nodes)
    @test size(fixed2.overlap, 1) == diagnostics2.fixed_dimension
    @test size(fixed3.overlap, 1) == diagnostics3.fixed_dimension
    @test norm(fixed2.overlap - I, Inf) < 1.0e-8
    @test norm(fixed3.overlap - I, Inf) < 1.0e-8
    @test all(isfinite, fixed2.weights)
    @test all(isfinite, fixed3.weights)
    @test minimum(fixed2.weights) > 0.0
    @test minimum(fixed3.weights) > 0.0
    @test all(isfinite, fixed2.fixed_centers)
    @test all(isfinite, fixed3.fixed_centers)

    mktempdir() do dir
        report_path = joinpath(dir, "square_lattice_nested_geometry_report.txt")
        written = write_axis_aligned_homonuclear_square_lattice_nested_geometry_report(
            report_path,
            basis3;
            nside = 5,
            min_in_plane_aspect_ratio = 0.15,
        )
        @test written.leaf_count == diagnostics3.leaf_count
        report_text = read(report_path, String)
        @test occursin("# retain_xy = (3, 3)", report_text)
        @test occursin("# shell_increment = 98", report_text)
        @test occursin("# ownership_multi_owned_row_count = 0", report_text)
        @test occursin("shared_shell[1].source_box =", report_text)
        @test occursin("shared_shell[1].next_inner_box =", report_text)
        @test occursin("shared_shell[1].source_point_count =", report_text)
        @test occursin("shared_shell[1].retained_fixed_count =", report_text)
        @test occursin("candidate[1].split_family = split_x_ternary", report_text)
        @test occursin("candidate[2].split_family = split_y_ternary", report_text)
        @test occursin("[node root_1]", report_text)
    end
end

@testset "Experimental axis-aligned homonuclear square-lattice nested QW consumer path" begin
    basis2, path2, check2 =
        _axis_aligned_homonuclear_square_lattice_nested_qw_fixture(; n = 2, spacing = 1.4)
    basis3, path3, check3 =
        _axis_aligned_homonuclear_square_lattice_nested_qw_fixture(; n = 3, spacing = 1.2)

    for (basis, path, check, n, expected_dim) in (
        (basis2, path2, check2, 2, 773),
        (basis3, path3, check3, 3, 703),
    )
        operators = path.operators
        fixed_block = path.fixed_block
        diagnostics = path.diagnostics
        @test basis isa AxisAlignedHomonuclearSquareLatticeQWBasis3D
        @test path isa GaussletBases.ExperimentalAxisAlignedHomonuclearSquareLatticeNestedQWPath
        @test path.basis === basis
        @test path.source === diagnostics.source
        @test path.fixed_block === fixed_block
        @test path.min_in_plane_aspect_ratio == 0.15
        @test operators isa GaussletBases.QiuWhiteResidualGaussianOperators
        @test operators.basis === fixed_block
        @test operators.gaussian_data === nothing
        @test operators.interaction_treatment == :ggt_nearest
        @test operators.residual_count == 0
        @test operators.gausslet_count == size(fixed_block.overlap, 1)
        @test size(operators.residual_centers) == (0, 3)
        @test size(operators.residual_widths) == (0, 3)
        @test norm(operators.overlap - I, Inf) < 1.0e-8
        @test norm(operators.one_body_hamiltonian - transpose(operators.one_body_hamiltonian), Inf) < 1.0e-10
        @test norm(operators.interaction_matrix - transpose(operators.interaction_matrix), Inf) < 1.0e-10
        @test all(isfinite, operators.one_body_hamiltonian)
        @test all(isfinite, operators.interaction_matrix)
        @test all(isfinite, fixed_block.weights)
        @test minimum(fixed_block.weights) > 0.0
        @test all(isfinite, fixed_block.fixed_centers)
        @test size(fixed_block.overlap, 1) == diagnostics.fixed_dimension
        @test size(fixed_block.overlap, 1) == expected_dim
        @test check.overlap_error < 1.0e-8
        @test isfinite(check.orbital_energy)
        @test isfinite(check.vee_expectation)
        @test length(basis.nuclei) == n * n
        @test Int[shell.retained_fixed_count for shell in diagnostics.shared_shell_provenance] ==
              diagnostics.shared_shell_dimensions
        common_contract = GaussletBases._nested_source_common_contract(path.source)
        @test common_contract.fixed_dimension == diagnostics.fixed_dimension
        @test common_contract.contract_audit == diagnostics.contract_audit
        @test common_contract.shared_shell_dimensions == diagnostics.shared_shell_dimensions
        @test common_contract.shared_shell_provenance == diagnostics.shared_shell_provenance
        @test common_contract.leaf_count == diagnostics.leaf_count
    end

    @test path2.diagnostics.root_node.did_split
    @test path2.diagnostics.root_node.accepted_candidate_index == 1
    @test path2.diagnostics.leaf_count == 4
    @test path2.diagnostics.root_node.child_count == 2
    @test path2.diagnostics.root_node.min_in_plane_aspect_ratio == 0.15
    @test path2.diagnostics.root_node.candidate_summaries[1].split_family == :split_x_binary
    @test path2.diagnostics.root_node.candidate_summaries[1].did_split

    @test path3.diagnostics.root_node.did_split
    @test path3.diagnostics.root_node.accepted_candidate_index == 1
    @test path3.diagnostics.leaf_count == 9
    @test path3.diagnostics.root_node.child_count == 3
    @test path3.diagnostics.root_node.min_in_plane_aspect_ratio == 0.15
    @test path3.diagnostics.root_node.candidate_summaries[1].split_family == :split_x_ternary
    @test path3.diagnostics.root_node.candidate_summaries[2].split_family == :split_y_ternary
    @test path3.diagnostics.root_node.candidate_summaries[1].did_split
    @test path3.diagnostics.root_node.candidate_summaries[1].child_in_plane_aspect_ratios[2] < 0.2
    @test path3.diagnostics.root_node.candidate_summaries[1].child_in_plane_aspect_ratios[2] >=
        path3.min_in_plane_aspect_ratio
    @test !path3.diagnostics.root_node.local_resolution_warning

    square_context = GaussletBases._normalized_nested_source_frontend_context(
        basis3;
        nside = 5,
        min_in_plane_aspect_ratio = 0.15,
    )
    square_source = GaussletBases._nested_source_frontend_source(square_context)
    square_fixed = GaussletBases._nested_source_fixed_block(square_source)
    square_diagnostics = GaussletBases._nested_source_geometry_diagnostics(square_source)
    square_common_contract = GaussletBases._nested_source_common_contract(square_source)
    @test size(square_fixed.fixed_block.overlap, 1) == path3.diagnostics.fixed_dimension
    @test square_diagnostics.fixed_dimension == path3.diagnostics.fixed_dimension
    @test square_diagnostics.root_node.accepted_candidate_index ==
          path3.diagnostics.root_node.accepted_candidate_index
    @test square_common_contract.fixed_dimension == path3.diagnostics.fixed_dimension
    @test square_common_contract.leaf_count == path3.diagnostics.leaf_count
end

@testset "Experimental homonuclear square-lattice nested dense export" begin
    basis2, path2, _check2 =
        _axis_aligned_homonuclear_square_lattice_nested_qw_fixture(; n = 2, spacing = 1.4)
    basis3, path3, _check3 =
        _axis_aligned_homonuclear_square_lattice_nested_qw_fixture(; n = 3, spacing = 1.2)

    for (basis, path, n, natoms, expected_family) in (
        (basis2, path2, 2, 4, "split_x_binary"),
        (basis3, path3, 3, 9, "split_x_ternary"),
    )
        payload_data = experimental_homonuclear_square_lattice_nested_dense_payload(
            path;
            meta = (example = "test_experimental_square_lattice_nested_dense_export",),
        )
        @test payload_data.bridge_meta["format"] ==
            "experimental_homonuclear_square_lattice_nested_dense_v1"
        @test payload_data.bridge_meta["experimental"]
        @test payload_data.bridge_meta["lattice_size"] == n
        @test payload_data.bridge_meta["fixed_dimension"] == size(path.fixed_block.overlap, 1)
        @test payload_data.bridge_meta["leaf_count"] == path.diagnostics.leaf_count
        @test payload_data.bridge_meta["root_accepted_candidate_index"] ==
            something(path.diagnostics.root_node.accepted_candidate_index, 0)
        @test payload_data.bridge_meta["root_accepted_split_family"] == expected_family
        @test payload_data.bridge_meta["root_child_count"] == path.diagnostics.root_node.child_count
        @test payload_data.bridge_meta["min_in_plane_aspect_ratio"] == path.min_in_plane_aspect_ratio
        @test payload_data.bridge_meta["residual_sector_empty"]
        @test payload_data.bridge_meta["coordinate_provenance"] == "uniform_square_spacing"
        @test size(payload_data.payload["S"]) == size(path.operators.overlap)
        @test size(payload_data.payload["H1"]) == size(path.operators.one_body_hamiltonian)
        @test size(payload_data.payload["Vee"]) == size(path.operators.interaction_matrix)
        @test size(payload_data.payload["basis_centers"]) == size(path.fixed_block.fixed_centers)
        @test size(payload_data.payload["nuclear_coordinates_xyz"]) == (natoms, 3)
        @test payload_data.payload["nuclear_charges"] == fill(1.0, natoms)
        @test payload_data.payload["lattice_x_coordinates"] == basis.x_coordinates
        @test payload_data.payload["lattice_y_coordinates"] == basis.y_coordinates
        @test payload_data.payload["orbital_labels"] ==
            String[orbital.label for orbital in orbitals(path.operators)]
        @test payload_data.payload["geometry_report_text"] isa String
        @test occursin("min_in_plane_aspect_ratio = 0.15", payload_data.payload["geometry_report_text"])
        @test occursin("[node root]", payload_data.payload["geometry_report_text"])
        @test size(payload_data.payload["root_accepted_child_planar_counts"], 2) == 2
        @test size(payload_data.payload["root_accepted_child_physical_widths"], 2) == 3
        @test payload_data.meta_values["manifest/contract/status"] == "experimental"
        @test payload_data.meta_values["manifest/source/lattice_size"] == n
        @test payload_data.meta_values["manifest/source/min_in_plane_aspect_ratio"] ==
            path.min_in_plane_aspect_ratio
        @test payload_data.meta_values["manifest/source/fixed_dimension"] ==
            path.diagnostics.fixed_dimension
    end

    payload3 = experimental_homonuclear_square_lattice_nested_dense_payload(
        path3;
        meta = (example = "test_experimental_square_lattice_nested_dense_export",),
    )
    @test payload3.payload["root_accepted_child_in_plane_aspect_ratios"][2] < 0.2

    mktempdir() do dir
        export_path = joinpath(dir, "h3_square_lattice_nested_dense.jld2")
        @test write_experimental_homonuclear_square_lattice_nested_dense_jld2(
            export_path,
            path3;
            meta = (example = "test_experimental_square_lattice_nested_dense_export",),
        ) == export_path
        jldopen(export_path, "r") do file
            @test String(file["bridge/format"]) ==
                "experimental_homonuclear_square_lattice_nested_dense_v1"
            @test Bool(file["bridge/experimental"])
            @test Int(file["bridge/lattice_size"]) == 3
            @test String(file["bridge/root_accepted_split_family"]) == "split_x_ternary"
            @test String(file["bridge/root_accepted_split_axis"]) == "x"
            @test Int(file["bridge/root_accepted_candidate_index"]) == 1
            @test Int(file["bridge/root_child_count"]) == 3
            @test Float64(file["bridge/min_in_plane_aspect_ratio"]) == 0.15
            @test Bool(file["bridge/residual_sector_empty"])
            @test size(file["S"]) == size(path3.operators.overlap)
            @test size(file["H1"]) == size(path3.operators.one_body_hamiltonian)
            @test size(file["Vee"]) == size(path3.operators.interaction_matrix)
            @test size(file["nuclear_coordinates_xyz"]) == (9, 3)
            @test size(file["root_accepted_child_planar_counts"]) == (3, 2)
            @test Float64(file["root_accepted_child_in_plane_aspect_ratios"][2]) < 0.2
            @test String(file["meta/producer"]) ==
                "GaussletBases.write_experimental_homonuclear_square_lattice_nested_dense_jld2"
            @test String(file["meta/manifest/contract/status"]) == "experimental"
            @test Int(file["meta/manifest/source/lattice_size"]) == 3
            @test Float64(file["meta/manifest/source/min_in_plane_aspect_ratio"]) == 0.15
            @test String(file["meta/example"]) ==
                "test_experimental_square_lattice_nested_dense_export"
        end
    end
end

@testset "Bond-aligned homonuclear chain nested geometry diagnostics" begin
    basis3, source3, fixed3, diagnostics3 = _bond_aligned_homonuclear_chain_nested_fixture(; natoms = 3, spacing = 1.2)
    basis4, source4, fixed4, diagnostics4 = _bond_aligned_homonuclear_chain_nested_fixture(; natoms = 4, spacing = 1.2)
    basis5, source5, fixed5, diagnostics5 = _bond_aligned_homonuclear_chain_nested_fixture(; natoms = 5, spacing = 1.2)
    basis3r, source3r, fixed3r, diagnostics3r = _bond_aligned_homonuclear_chain_nested_fixture(; natoms = 3, spacing = 1.2, odd_chain_policy = :central_ternary_relaxed)
    basis5r, source5r, fixed5r, diagnostics5r = _bond_aligned_homonuclear_chain_nested_fixture(; natoms = 5, spacing = 1.2, odd_chain_policy = :central_ternary_relaxed)

    @test basis3 isa BondAlignedHomonuclearChainQWBasis3D
    @test basis4 isa BondAlignedHomonuclearChainQWBasis3D
    @test basis5 isa BondAlignedHomonuclearChainQWBasis3D
    @test diagnostics3.retention_contract.retain_xy == (3, 3)
    @test diagnostics4.retention_contract.retain_xy == (3, 3)
    @test diagnostics5.retention_contract.retain_xy == (3, 3)
    @test diagnostics3.retention_contract.retain_x_edge == 3
    @test diagnostics4.retention_contract.retain_x_edge == 3
    @test diagnostics5.retention_contract.retain_x_edge == 3
    @test diagnostics3.retention_contract.shell_increment == 98
    @test diagnostics4.retention_contract.shell_increment == 98
    @test diagnostics5.retention_contract.shell_increment == 98
    @test diagnostics3.retention_contract.matches_nside_default
    @test diagnostics4.retention_contract.matches_nside_default
    @test diagnostics5.retention_contract.matches_nside_default
    @test diagnostics3.contract_audit.full_parent_working_box
    @test diagnostics4.contract_audit.full_parent_working_box
    @test diagnostics5.contract_audit.full_parent_working_box
    @test diagnostics3.contract_audit.missing_row_count == 0
    @test diagnostics4.contract_audit.missing_row_count == 0
    @test diagnostics5.contract_audit.missing_row_count == 0
    @test diagnostics3.contract_audit.ownership_unowned_row_count == 0
    @test diagnostics4.contract_audit.ownership_unowned_row_count == 0
    @test diagnostics5.contract_audit.ownership_unowned_row_count == 0
    @test diagnostics3.contract_audit.ownership_multi_owned_row_count == 0
    @test diagnostics4.contract_audit.ownership_multi_owned_row_count == 0
    @test diagnostics5.contract_audit.ownership_multi_owned_row_count == 0
    @test diagnostics3.shared_shells_match_contract
    @test diagnostics4.shared_shells_match_contract
    @test diagnostics5.shared_shells_match_contract
    @test Int[shell.retained_fixed_count for shell in diagnostics3.shared_shell_provenance] ==
          diagnostics3.shared_shell_dimensions
    @test Int[shell.retained_fixed_count for shell in diagnostics4.shared_shell_provenance] ==
          diagnostics4.shared_shell_dimensions
    @test Int[shell.retained_fixed_count for shell in diagnostics5.shared_shell_provenance] ==
          diagnostics5.shared_shell_dimensions
    @test all(
        shell.source_point_count ==
        prod(length.(shell.source_box)) - prod(length.(shell.next_inner_box)) for
        shell in diagnostics4.shared_shell_provenance
    )
    @test diagnostics3.leaf_count == 1
    @test diagnostics4.leaf_count == 2
    @test diagnostics5.leaf_count == 1
    @test !diagnostics3.root_node.did_split
    @test diagnostics4.root_node.did_split
    @test !diagnostics5.root_node.did_split
    @test diagnostics3.root_node.child_count == 0
    @test diagnostics4.root_node.child_count == 2
    @test diagnostics5.root_node.child_count == 0
    @test diagnostics3.root_node.odd_chain_policy == :strict_current
    @test diagnostics5.root_node.odd_chain_policy == :strict_current
    @test diagnostics3.root_node.accepted_candidate_index === nothing
    @test diagnostics4.root_node.accepted_candidate_index == 2
    @test diagnostics5.root_node.accepted_candidate_index === nothing
    @test diagnostics3.root_node.candidate_summaries[1].split_kind == :ternary
    @test diagnostics4.root_node.candidate_summaries[2].split_kind == :binary
    @test diagnostics5.root_node.candidate_summaries[2].split_kind == :ternary
    @test diagnostics3.root_node.candidate_summaries[1].nucleus_ranges == [1:1, 2:2, 3:3]
    @test diagnostics4.root_node.candidate_summaries[2].nucleus_ranges == [1:2, 3:4]
    @test diagnostics5.root_node.candidate_summaries[2].nucleus_ranges == [1:2, 3:3, 4:5]
    @test diagnostics3.root_node.local_resolution_warning
    @test !diagnostics4.root_node.local_resolution_warning
    @test diagnostics5.root_node.local_resolution_warning
    @test length(source3.root_geometry.child_nodes) == 0
    @test length(source4.root_geometry.child_nodes) == 2
    @test length(source5.root_geometry.child_nodes) == 0
    @test !diagnostics3.root_node.candidate_summaries[1].did_split
    @test diagnostics4.root_node.candidate_summaries[2].did_split
    @test diagnostics5.root_node.candidate_summaries[2].count_eligible
    @test !diagnostics5.root_node.candidate_summaries[2].shape_eligible
    @test diagnostics5.root_node.candidate_summaries[2].child_parallel_counts == [6, 3, 6]
    @test diagnostics5.root_node.candidate_summaries[2].child_parallel_to_transverse_ratios[2] <
        diagnostics5.root_node.odd_chain_policy_thresholds.center_parallel_to_transverse_ratio_min
    @test diagnostics4.root_node.candidate_summaries[2].midpoint_values == [0.0]
    @test abs(
        length(diagnostics4.root_node.candidate_summaries[2].child_boxes[1][3]) -
        length(diagnostics4.root_node.candidate_summaries[2].child_boxes[2][3]),
    ) <= 1
    @test all(
        widths[3] > 0.0 for
        widths in diagnostics4.root_node.candidate_summaries[2].child_physical_widths
    )
    @test size(fixed3.overlap, 1) == diagnostics3.fixed_dimension
    @test size(fixed4.overlap, 1) == diagnostics4.fixed_dimension
    @test size(fixed5.overlap, 1) == diagnostics5.fixed_dimension
    @test norm(fixed3.overlap - I, Inf) < 1.0e-8
    @test norm(fixed4.overlap - I, Inf) < 1.0e-8
    @test norm(fixed5.overlap - I, Inf) < 1.0e-8
    @test all(isfinite, fixed3.weights)
    @test all(isfinite, fixed4.weights)
    @test all(isfinite, fixed5.weights)
    @test all(isfinite, fixed3.fixed_centers)
    @test all(isfinite, fixed4.fixed_centers)
    @test all(isfinite, fixed5.fixed_centers)

    @test diagnostics3r.root_node.odd_chain_policy == :central_ternary_relaxed
    @test diagnostics5r.root_node.odd_chain_policy == :central_ternary_relaxed
    @test diagnostics3r.root_node.did_split
    @test diagnostics5r.root_node.did_split
    @test diagnostics3r.root_node.accepted_candidate_index == 1
    @test diagnostics5r.root_node.accepted_candidate_index == 2
    @test diagnostics3r.root_node.child_count == 3
    @test diagnostics5r.root_node.child_count == 3
    @test diagnostics3r.root_node.candidate_summaries[1].did_split
    @test diagnostics5r.root_node.candidate_summaries[2].did_split
    @test diagnostics3r.root_node.candidate_summaries[1].child_parallel_counts == [4, 3, 4]
    @test diagnostics5r.root_node.candidate_summaries[2].child_parallel_counts == [6, 3, 6]
    @test diagnostics3r.root_node.odd_chain_policy_thresholds.center_parallel_count_min == 3
    @test diagnostics5r.root_node.odd_chain_policy_thresholds.center_parallel_to_transverse_ratio_min == 0.35
    @test size(fixed3r.overlap, 1) == diagnostics3r.fixed_dimension
    @test size(fixed5r.overlap, 1) == diagnostics5r.fixed_dimension
    @test norm(fixed3r.overlap - I, Inf) < 1.0e-8
    @test norm(fixed5r.overlap - I, Inf) < 1.0e-8
    @test all(isfinite, fixed3r.weights)
    @test all(isfinite, fixed5r.weights)
    @test all(isfinite, fixed3r.fixed_centers)
    @test all(isfinite, fixed5r.fixed_centers)
    @test length(source3r.root_geometry.child_nodes) == 3
    @test length(source5r.root_geometry.child_nodes) == 3

    mktempdir() do dir
        report_path = joinpath(dir, "chain_nested_report.txt")
        written = write_bond_aligned_homonuclear_chain_nested_geometry_report(
            report_path,
            basis5r;
            nside = 5,
            odd_chain_policy = :central_ternary_relaxed,
        )
        report_text = read(report_path, String)
        @test written.leaf_count == diagnostics5r.leaf_count
        @test occursin("[node root]", report_text)
        @test occursin("# retain_xy = (3, 3)", report_text)
        @test occursin("# shell_increment = 98", report_text)
        @test occursin("# ownership_multi_owned_row_count = 0", report_text)
        @test occursin("shared_shell[1].source_box =", report_text)
        @test occursin("shared_shell[1].next_inner_box =", report_text)
        @test occursin("shared_shell[1].source_point_count =", report_text)
        @test occursin("shared_shell[1].retained_fixed_count =", report_text)
        @test occursin("odd_chain_policy = central_ternary_relaxed", report_text)
        @test occursin("candidate[2].child_parallel_counts = [6, 3, 6]", report_text)
        @test occursin("candidate[2].accepted = true", report_text)
    end
end

@testset "Non-atomic nested routes default to nside-driven complete-shell retention" begin
    diatomic_basis = bond_aligned_homonuclear_qw_basis(
        bond_length = 1.4,
        core_spacing = 0.5,
        xmax_parallel = 4.0,
        xmax_transverse = 3.0,
        bond_axis = :z,
    )
    diatomic_diagnostics = bond_aligned_diatomic_nested_geometry_diagnostics(
        diatomic_basis;
        nside = 7,
    )
    @test diatomic_diagnostics.nside == 7
    @test diatomic_diagnostics.child_shell_retention_contract.retain_xy == (5, 5)
    @test diatomic_diagnostics.child_shell_retention_contract.retain_x_edge == 5
    @test diatomic_diagnostics.child_shell_retention_contract.shell_increment == 218
    @test diatomic_diagnostics.child_shell_retention_contract.matches_nside_default
    @test diatomic_diagnostics.shared_shell_retention_contract.matches_nside_default
    @test diatomic_diagnostics.contract_audit.support_count ==
        diatomic_diagnostics.contract_audit.expected_support_count
    @test diatomic_diagnostics.contract_audit.missing_row_count == 0
    @test diatomic_diagnostics.contract_audit.ownership_unowned_row_count == 0
    @test diatomic_diagnostics.contract_audit.ownership_multi_owned_row_count == 0

    chain_basis = bond_aligned_homonuclear_chain_qw_basis(
        natoms = 4,
        spacing = 1.2,
        core_spacing = 0.5,
        xmax_parallel = 2.0,
        xmax_transverse = 3.0,
        chain_axis = :z,
    )
    chain_diagnostics = bond_aligned_homonuclear_chain_nested_geometry_diagnostics(
        chain_basis;
        nside = 7,
    )
    @test chain_diagnostics.nside == 7
    @test chain_diagnostics.retention_contract.retain_xy == (5, 5)
    @test chain_diagnostics.retention_contract.retain_x_edge == 5
    @test chain_diagnostics.retention_contract.shell_increment == 218
    @test chain_diagnostics.retention_contract.matches_nside_default
    @test chain_diagnostics.shared_shell_dimensions == [218]
    @test length(chain_diagnostics.shared_shell_provenance) == 1
    @test chain_diagnostics.shared_shell_provenance[1].retained_fixed_count == 218
    @test chain_diagnostics.shared_shells_match_contract
    @test chain_diagnostics.contract_audit.full_parent_working_box
    @test chain_diagnostics.contract_audit.support_count ==
        chain_diagnostics.contract_audit.expected_support_count
    @test chain_diagnostics.contract_audit.missing_row_count == 0
    @test chain_diagnostics.contract_audit.ownership_unowned_row_count == 0
    @test chain_diagnostics.contract_audit.ownership_multi_owned_row_count == 0

    square_basis = axis_aligned_homonuclear_square_lattice_qw_basis(
        n = 3,
        spacing = 1.2,
        core_spacing = 0.5,
        xmax_in_plane = 3.0,
        xmax_transverse = 3.0,
    )
    square_diagnostics = axis_aligned_homonuclear_square_lattice_nested_geometry_diagnostics(
        square_basis;
        nside = 7,
    )
    @test square_diagnostics.nside == 7
    @test square_diagnostics.retention_contract.retain_xy == (5, 5)
    @test square_diagnostics.retention_contract.retain_x_edge == 5
    @test square_diagnostics.retention_contract.shell_increment == 218
    @test square_diagnostics.retention_contract.matches_nside_default
    @test square_diagnostics.shared_shell_dimensions == [218]
    @test length(square_diagnostics.shared_shell_provenance) == 1
    @test square_diagnostics.shared_shell_provenance[1].retained_fixed_count == 218
    @test square_diagnostics.shared_shells_match_contract
    @test square_diagnostics.contract_audit.full_parent_working_box
    @test square_diagnostics.contract_audit.support_count ==
        square_diagnostics.contract_audit.expected_support_count
    @test square_diagnostics.contract_audit.missing_row_count == 0
    @test square_diagnostics.contract_audit.ownership_unowned_row_count == 0
    @test square_diagnostics.contract_audit.ownership_multi_owned_row_count == 0
end

@testset "Bond-aligned diatomic nested source reuse path" begin
    basis = bond_aligned_homonuclear_qw_basis(
        bond_length = 1.4,
        core_spacing = 0.5,
        xmax_parallel = 4.0,
        xmax_transverse = 3.0,
        bond_axis = :z,
    )

    diagnostics_via_basis = bond_aligned_diatomic_nested_geometry_diagnostics(
        basis;
        nside = 5,
    )
    source = bond_aligned_diatomic_nested_fixed_source(
        basis;
        nside = 5,
    )
    diagnostics_via_source = bond_aligned_diatomic_nested_geometry_diagnostics(source)
    context = GaussletBases._normalized_nested_source_frontend_context(
        basis;
        nside = 5,
    )
    source_via_context = GaussletBases._nested_source_frontend_source(context)
    common_contract = GaussletBases._nested_source_common_contract(source)
    common_contract_via_context = GaussletBases._nested_source_common_contract(source_via_context)
    diagnostics_via_context = GaussletBases._nested_source_geometry_diagnostics(source_via_context)
    fixed_via_basis = bond_aligned_diatomic_nested_fixed_block(
        basis;
        nside = 5,
    )
    fixed_via_source = bond_aligned_diatomic_nested_fixed_block(source)
    fixed_via_context = GaussletBases._nested_source_fixed_block(source_via_context)
    source_payload = bond_aligned_diatomic_source_geometry_payload(source)
    source_slice = bond_aligned_diatomic_plane_slice(
        source_payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-5,
    )

    @test diagnostics_via_source.nside == diagnostics_via_basis.nside
    @test diagnostics_via_source.geometry.did_split == diagnostics_via_basis.geometry.did_split
    @test diagnostics_via_source.geometry.count_eligible ==
        diagnostics_via_basis.geometry.count_eligible
    @test diagnostics_via_source.geometry.unsplit_aspect_eligible ==
        diagnostics_via_basis.geometry.unsplit_aspect_eligible
    @test diagnostics_via_source.geometry.shape_eligible ==
        diagnostics_via_basis.geometry.shape_eligible
    @test diagnostics_via_source.geometry.split_index ==
        diagnostics_via_basis.geometry.split_index
    @test diagnostics_via_source.geometry.working_box ==
        diagnostics_via_basis.geometry.working_box
    @test diagnostics_via_source.geometry.shared_midpoint_box ==
        diagnostics_via_basis.geometry.shared_midpoint_box
    @test diagnostics_via_source.geometry.child_boxes ==
        diagnostics_via_basis.geometry.child_boxes
    @test maximum(
        abs,
        reduce(vcat, (
            collect(widths_source .- widths_basis) for
            (widths_source, widths_basis) in zip(
                diagnostics_via_source.geometry.child_physical_widths,
                diagnostics_via_basis.geometry.child_physical_widths,
            )
        );
        init = Float64[]),
    ) < 1.0e-12
    @test diagnostics_via_source.child_shell_retention_contract ==
        diagnostics_via_basis.child_shell_retention_contract
    @test diagnostics_via_source.shared_shell_retention_contract ==
        diagnostics_via_basis.shared_shell_retention_contract
    @test diagnostics_via_source.shared_shell_dimensions ==
        diagnostics_via_basis.shared_shell_dimensions
    @test diagnostics_via_source.shared_shell_provenance ==
        diagnostics_via_basis.shared_shell_provenance
    generic_contract = GaussletBases._nested_glass_box_contract(source)
    @test common_contract.fixed_dimension == diagnostics_via_source.fixed_dimension
    @test common_contract.contract_audit == diagnostics_via_source.contract_audit
    @test common_contract.shared_shell_dimensions == diagnostics_via_source.shared_shell_dimensions
    @test common_contract.shared_shell_provenance == diagnostics_via_source.shared_shell_provenance
    @test common_contract.leaf_count == diagnostics_via_source.child_sequence_count
    @test generic_contract.fixed_dimension == common_contract.fixed_dimension
    @test generic_contract.contract_audit == common_contract.contract_audit
    @test generic_contract.layer_dimensions == common_contract.shared_shell_dimensions
    @test generic_contract.layer_provenance == common_contract.shared_shell_provenance
    @test generic_contract.leaf_count == common_contract.leaf_count
    @test common_contract_via_context.fixed_dimension == common_contract.fixed_dimension
    @test common_contract_via_context.contract_audit == common_contract.contract_audit
    @test common_contract_via_context.shared_shell_dimensions ==
          common_contract.shared_shell_dimensions
    @test common_contract_via_context.shared_shell_provenance ==
          common_contract.shared_shell_provenance
    @test common_contract_via_context.leaf_count == common_contract.leaf_count
    @test diagnostics_via_context.fixed_dimension == diagnostics_via_source.fixed_dimension
    @test diagnostics_via_context.child_sequence_dimensions ==
          diagnostics_via_source.child_sequence_dimensions
    @test diagnostics_via_source.child_sequence_dimensions ==
        diagnostics_via_basis.child_sequence_dimensions
    @test diagnostics_via_source.fixed_dimension == diagnostics_via_basis.fixed_dimension
    @test diagnostics_via_source.contract_audit.support_count ==
        diagnostics_via_basis.contract_audit.support_count
    @test diagnostics_via_source.contract_audit.expected_support_count ==
        diagnostics_via_basis.contract_audit.expected_support_count
    @test diagnostics_via_source.contract_audit.missing_row_count ==
        diagnostics_via_basis.contract_audit.missing_row_count
    @test source_payload.bond_axis == basis.bond_axis
    @test length(source_payload.points) == prod(length.(source.geometry.parent_box))
    @test source_slice.selected_count > 0
    @test !hasproperty(source.sequence.packet, :term_storage)
    @test !hasproperty(source.sequence.packet, :gaussian_terms)
    @test !hasproperty(source.sequence.packet, :pair_terms)
    @test !isnothing(source.sequence.packet.gaussian_sum)
    @test !isnothing(source.sequence.packet.pair_sum)
    @test all(isnothing(sequence.support_states) for sequence in source.child_sequences)
    @test all(isnothing(sequence.packet) for sequence in source.child_sequences)
    @test !isnothing(source.sequence.support_states)
    @test !isnothing(source.sequence.packet)

    @test fixed_via_source.source === source
    @test norm(fixed_via_source.fixed_block.overlap - fixed_via_basis.fixed_block.overlap, Inf) <
        1.0e-12
    @test norm(
        fixed_via_source.fixed_block.coefficient_matrix -
        fixed_via_basis.fixed_block.coefficient_matrix,
        Inf,
    ) < 1.0e-12
    @test fixed_via_source.fixed_block.support_indices == fixed_via_basis.fixed_block.support_indices
    @test norm(
        source_via_context.sequence.coefficient_matrix - source.sequence.coefficient_matrix,
        Inf,
    ) < 1.0e-12
    @test norm(
        fixed_via_context.fixed_block.coefficient_matrix -
        fixed_via_source.fixed_block.coefficient_matrix,
        Inf,
    ) < 1.0e-12
    @test fixed_via_source.fixed_block.working_box == fixed_via_basis.fixed_block.working_box
end

@testset "Experimental bond-aligned homonuclear chain nested QW consumer path" begin
    basis3, path3, check3 = _bond_aligned_homonuclear_chain_nested_qw_fixture(; natoms = 3, spacing = 1.2)
    basis4, path4, check4 = _bond_aligned_homonuclear_chain_nested_qw_fixture(; natoms = 4, spacing = 1.2)
    basis5, path5, check5 = _bond_aligned_homonuclear_chain_nested_qw_fixture(; natoms = 5, spacing = 1.2)

    for (basis, path, check, natoms) in ((basis3, path3, check3, 3), (basis4, path4, check4, 4), (basis5, path5, check5, 5))
        operators = path.operators
        fixed_block = path.fixed_block
        diagnostics = path.diagnostics
        @test basis isa BondAlignedHomonuclearChainQWBasis3D
        @test path isa GaussletBases.ExperimentalBondAlignedHomonuclearChainNestedQWPath
        @test path.basis === basis
        @test path.source === diagnostics.source
        @test path.fixed_block === fixed_block
        @test operators isa GaussletBases.QiuWhiteResidualGaussianOperators
        @test operators.basis === fixed_block
        @test operators.gaussian_data === nothing
        @test operators.interaction_treatment == :ggt_nearest
        @test operators.residual_count == 0
        @test operators.gausslet_count == size(fixed_block.overlap, 1)
        @test size(operators.residual_centers) == (0, 3)
        @test size(operators.residual_widths) == (0, 3)
        @test norm(operators.overlap - I, Inf) < 1.0e-8
        @test norm(operators.one_body_hamiltonian - transpose(operators.one_body_hamiltonian), Inf) < 1.0e-10
        @test norm(operators.interaction_matrix - transpose(operators.interaction_matrix), Inf) < 1.0e-10
        @test all(isfinite, operators.one_body_hamiltonian)
        @test all(isfinite, operators.interaction_matrix)
        @test all(isfinite, fixed_block.weights)
        @test minimum(fixed_block.weights) > 0.0
        @test all(isfinite, fixed_block.fixed_centers)
        @test diagnostics.root_node.odd_chain_policy == :central_ternary_relaxed
        @test path.odd_chain_policy == :central_ternary_relaxed
        @test size(fixed_block.overlap, 1) == diagnostics.fixed_dimension
        @test size(fixed_block.overlap, 1) < 1000
        @test check.overlap_error < 1.0e-8
        @test isfinite(check.orbital_energy)
        @test isfinite(check.vee_expectation)
        @test length(basis.nuclei) == natoms
        @test Int[shell.retained_fixed_count for shell in diagnostics.shared_shell_provenance] ==
              diagnostics.shared_shell_dimensions
        common_contract = GaussletBases._nested_source_common_contract(path.source)
        @test common_contract.fixed_dimension == diagnostics.fixed_dimension
        @test common_contract.contract_audit == diagnostics.contract_audit
        @test common_contract.shared_shell_dimensions == diagnostics.shared_shell_dimensions
        @test common_contract.shared_shell_provenance == diagnostics.shared_shell_provenance
        @test common_contract.leaf_count == diagnostics.leaf_count
    end

    @test path3.diagnostics.root_node.did_split
    @test path3.diagnostics.root_node.accepted_candidate_index == 1
    @test path3.diagnostics.leaf_count == 3
    @test path3.diagnostics.root_node.child_count == 3

    @test path4.diagnostics.root_node.did_split
    @test path4.diagnostics.root_node.accepted_candidate_index == 2
    @test path4.diagnostics.leaf_count == 2
    @test path4.diagnostics.root_node.child_count == 2

    @test path5.diagnostics.root_node.did_split
    @test path5.diagnostics.root_node.accepted_candidate_index == 2
    @test path5.diagnostics.leaf_count == 3
    @test path5.diagnostics.root_node.child_count == 3
    @test path5.diagnostics.root_node.candidate_summaries[2].child_parallel_counts == [6, 3, 6]
    @test path5.diagnostics.root_node.candidate_summaries[2].child_parallel_to_transverse_ratios[2] > 0.35

    chain_context = GaussletBases._normalized_nested_source_frontend_context(
        basis5;
        nside = 5,
        odd_chain_policy = :central_ternary_relaxed,
    )
    chain_source = GaussletBases._nested_source_frontend_source(chain_context)
    chain_fixed = GaussletBases._nested_source_fixed_block(chain_source)
    chain_diagnostics = GaussletBases._nested_source_geometry_diagnostics(chain_source)
    chain_common_contract = GaussletBases._nested_source_common_contract(chain_source)
    @test size(chain_fixed.fixed_block.overlap, 1) == path5.diagnostics.fixed_dimension
    @test chain_diagnostics.fixed_dimension == path5.diagnostics.fixed_dimension
    @test chain_diagnostics.root_node.accepted_candidate_index ==
          path5.diagnostics.root_node.accepted_candidate_index
    @test chain_common_contract.fixed_dimension == path5.diagnostics.fixed_dimension
    @test chain_common_contract.leaf_count == path5.diagnostics.leaf_count
end

@testset "Experimental homonuclear chain nested dense export" begin
    basis3, path3, _check3 = _bond_aligned_homonuclear_chain_nested_qw_fixture(; natoms = 3, spacing = 1.2)
    basis4, path4, _check4 = _bond_aligned_homonuclear_chain_nested_qw_fixture(; natoms = 4, spacing = 1.2)
    basis5, path5, _check5 = _bond_aligned_homonuclear_chain_nested_qw_fixture(; natoms = 5, spacing = 1.2)

    for (basis, path, natoms) in ((basis3, path3, 3), (basis4, path4, 4), (basis5, path5, 5))
        payload_data = experimental_homonuclear_chain_nested_dense_payload(
            path;
            meta = (example = "test_experimental_chain_nested_dense_export",),
        )
        @test payload_data.bridge_meta["format"] == "experimental_homonuclear_chain_nested_dense_v1"
        @test payload_data.bridge_meta["experimental"]
        @test payload_data.bridge_meta["odd_chain_policy"] == "central_ternary_relaxed"
        @test payload_data.bridge_meta["fixed_dimension"] == size(path.fixed_block.overlap, 1)
        @test payload_data.bridge_meta["leaf_count"] == path.diagnostics.leaf_count
        @test payload_data.bridge_meta["root_accepted_candidate_index"] ==
            something(path.diagnostics.root_node.accepted_candidate_index, 0)
        @test payload_data.bridge_meta["residual_sector_empty"]
        @test size(payload_data.payload["S"]) == size(path.operators.overlap)
        @test size(payload_data.payload["H1"]) == size(path.operators.one_body_hamiltonian)
        @test size(payload_data.payload["Vee"]) == size(path.operators.interaction_matrix)
        @test size(payload_data.payload["basis_centers"]) == size(path.fixed_block.fixed_centers)
        @test size(payload_data.payload["nuclear_coordinates_xyz"]) == (natoms, 3)
        @test payload_data.payload["nuclear_charges"] == fill(1.0, natoms)
        @test payload_data.payload["orbital_labels"] ==
            String[orbital.label for orbital in orbitals(path.operators)]
        @test payload_data.payload["geometry_report_text"] isa String
        @test occursin("odd_chain_policy = central_ternary_relaxed", payload_data.payload["geometry_report_text"])
        @test occursin("[node root]", payload_data.payload["geometry_report_text"])
        @test payload_data.meta_values["manifest/contract/status"] == "experimental"
        @test payload_data.meta_values["manifest/source/odd_chain_policy"] == "central_ternary_relaxed"
        @test payload_data.meta_values["manifest/source/fixed_dimension"] == path.diagnostics.fixed_dimension
    end

    mktempdir() do dir
        export_path = joinpath(dir, "h5_chain_nested_dense.jld2")
        @test write_experimental_homonuclear_chain_nested_dense_jld2(
            export_path,
            path5;
            meta = (example = "test_experimental_chain_nested_dense_export",),
        ) == export_path
        jldopen(export_path, "r") do file
            @test String(file["bridge/format"]) == "experimental_homonuclear_chain_nested_dense_v1"
            @test Bool(file["bridge/experimental"])
            @test String(file["bridge/odd_chain_policy"]) == "central_ternary_relaxed"
            @test Bool(file["bridge/root_did_split"])
            @test Int(file["bridge/root_accepted_candidate_index"]) == 2
            @test Int(file["bridge/fixed_dimension"]) == size(path5.fixed_block.overlap, 1)
            @test Bool(file["bridge/residual_sector_empty"])
            @test size(file["S"]) == size(path5.operators.overlap)
            @test size(file["H1"]) == size(path5.operators.one_body_hamiltonian)
            @test size(file["Vee"]) == size(path5.operators.interaction_matrix)
            @test size(file["nuclear_coordinates_xyz"]) == (5, 3)
            @test String(file["meta/producer"]) ==
                "GaussletBases.write_experimental_homonuclear_chain_nested_dense_jld2"
            @test String(file["meta/manifest/contract/status"]) == "experimental"
            @test String(file["meta/manifest/source/odd_chain_policy"]) == "central_ternary_relaxed"
            @test String(file["meta/example"]) == "test_experimental_chain_nested_dense_export"
        end
    end
end

@testset "Atomic hybrid anchor comparison" begin
    if !_legacy_basisfile_available()
        @test true
    else
        (
            _source_basis,
            bundle,
            _shell1_complete,
            _shell2_complete,
            _shell3_complete,
            _shell4_complete,
            _interval1,
            _interval2,
            _interval3,
            _interval4,
            _core5,
            complete_sequence,
            fixed_complete_sequence,
            legacy,
            baseline,
            complete_sequence_ops,
            baseline_check,
            complete_sequence_check,
        ) = _nested_qiu_white_complete_shell_sequence_fixture(count = 15)

        expansion = coulomb_gaussian_expansion(doacc = false)
        overlap_parent, one_body_parent, interaction_parent =
            _nested_parent_fixed_problem(bundle, expansion; Z = 2.0)
        parent_modes = eigen(Hermitian(one_body_parent), Hermitian(overlap_parent))
        parent_ground = parent_modes.vectors[:, 1]
        parent_ground_vee = _nested_vee_from_orbital(interaction_parent, parent_ground)
        projected_complete = _nested_fixed_projected_orbital(
            overlap_parent,
            fixed_complete_sequence,
            parent_ground,
        )
        projected_complete_vee = _nested_vee_from_orbital(
            GaussletBases._qwrg_fixed_block_interaction_matrix(fixed_complete_sequence, expansion),
            projected_complete,
        )
        complete_ground_capture, complete_ground_energy = _nested_projector_stats(
            overlap_parent,
            one_body_parent,
            fixed_complete_sequence,
            parent_ground,
        )
        complete_average4_capture =
            sum(
                _nested_projector_stats(
                    overlap_parent,
                    one_body_parent,
                    fixed_complete_sequence,
                    parent_modes.vectors[:, index],
                )[1] for index in 1:4
            ) / 4

        @test legacy isa LegacyAtomicGaussianSupplement
        @test legacy.lmax == 0
        @test complete_sequence.working_box == (2:14, 2:14, 2:14)
        @test baseline.gausslet_count == 15^3
        @test norm(fixed_complete_sequence.overlap - I, Inf) < 1.0e-10
        @test complete_sequence_check.overlap_error < 1.0e-10
        @test all(isfinite, fixed_complete_sequence.weights)
        @test minimum(fixed_complete_sequence.weights) > 0.0
        @test maximum(fixed_complete_sequence.weights) < 10.0
        @test abs(complete_sequence_check.orbital_energy - baseline_check.orbital_energy) < 2.0e-4
        @test abs(complete_sequence_check.vee_expectation - baseline_check.vee_expectation) < 1.0e-4
        @test abs(projected_complete_vee - parent_ground_vee) < 1.5e-4
        @test complete_ground_capture > 0.99999
        @test complete_average4_capture > 0.999
        @test complete_ground_energy - parent_modes.values[1] < 1.0e-4
    end
end

@testset "Atomic hierarchical core-only refinement" begin
    if !_legacy_basisfile_available()
        @test true
    else
        for count in (17, 15)
            (
                _source_basis,
                bundle,
                core5,
                complete_sequence,
                fixed_complete_sequence,
                complete_sequence_ops,
                complete_sequence_check,
                refined_core,
                refined_sequence,
                fixed_refined_sequence,
                refined_sequence_ops,
                refined_sequence_check,
                legacy,
                _baseline,
                _baseline_check,
            ) = _nested_qiu_white_hierarchical_core_fixture(count = count)

            expansion = coulomb_gaussian_expansion(doacc = false)
            overlap_parent, one_body_parent, interaction_parent =
                _nested_parent_fixed_problem(bundle, expansion; Z = 2.0)
            parent_modes = eigen(Hermitian(one_body_parent), Hermitian(overlap_parent))
            parent_ground = parent_modes.vectors[:, 1]
            parent_ground_vee = _nested_vee_from_orbital(interaction_parent, parent_ground)
            projected_refined = _nested_fixed_projected_orbital(
                overlap_parent,
                fixed_refined_sequence,
                parent_ground,
            )
            projected_refined_vee = _nested_vee_from_orbital(
                GaussletBases._qwrg_fixed_block_interaction_matrix(fixed_refined_sequence, expansion),
                projected_refined,
            )
            refined_ground_capture, refined_ground_energy = _nested_projector_stats(
                overlap_parent,
                one_body_parent,
                fixed_refined_sequence,
                parent_ground,
            )

            inner_core = (first(core5) + 1):(last(core5) - 1)
            @test legacy isa LegacyAtomicGaussianSupplement
            @test legacy.lmax == 0
            @test refined_core isa GaussletBases._CartesianNestedShellSequence3D
            @test refined_sequence isa GaussletBases._CartesianNestedShellSequence3D
            @test refined_core.working_box == (core5, core5, core5)
            @test refined_sequence.working_box == complete_sequence.working_box
            @test length(refined_core.core_indices) == length(inner_core)^3
            @test size(refined_core.coefficient_matrix, 2) == 53
            @test refined_sequence_ops.gausslet_count == 445
            @test refined_sequence_ops.gausslet_count < complete_sequence_ops.gausslet_count
            @test norm(fixed_refined_sequence.overlap - I, Inf) < 1.0e-10
            @test refined_sequence_check.overlap_error < 1.0e-10
            @test all(isfinite, fixed_refined_sequence.weights)
            @test minimum(fixed_refined_sequence.weights) > 0.0
            @test refined_ground_capture > 0.999
            @test abs(refined_sequence_check.vee_expectation - complete_sequence_check.vee_expectation) > 5.0e-4
            @test abs(projected_refined_vee - parent_ground_vee) > 1.0e-3
            @test refined_ground_energy - parent_modes.values[1] > 1.0e-2
        end
    end
end

@testset "Mapped ordinary one-body backends" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    mild_basis = build_basis(MappedUniformBasisSpec(:G10;
        count = 5,
        mapping = fit_asinh_mapping_for_strength(s = 0.5, npoints = 5, xmax = 6.0),
        reference_spacing = 1.0,
    ))
    mild_reference = mapped_ordinary_one_body_operators(
        mild_basis;
        exponents = expansion.exponents[1:3],
        backend = :numerical_reference,
    )
    mild_analytic = mapped_ordinary_one_body_operators(
        mild_basis;
        exponents = expansion.exponents[1:3],
        backend = :pgdg_experimental,
    )
    mild_localized = mapped_ordinary_one_body_operators(
        mild_basis;
        exponents = expansion.exponents[1:3],
        backend = :pgdg_localized_experimental,
    )
    mild_oracle = GaussletBases._mapped_ordinary_localized_oracle_operators(
        mild_basis;
        exponents = expansion.exponents[1:3],
    )
    (_, _, overlap_reference_localized, kinetic_reference_localized, _) =
        _localized_numerical_reference_1d(mild_basis, expansion.exponents[1:3])
    raw_localized = mapped_pgdg_localized(GaussletBases.mapped_pgdg_logfit_prototype(mild_basis))
    raw_localized_kinetic = kinetic_matrix(raw_localized)

    @test mild_reference isa MappedOrdinaryOneBody1D
    @test mild_analytic isa MappedOrdinaryOneBody1D
    @test mild_localized isa MappedOrdinaryOneBody1D
    @test mild_oracle isa MappedOrdinaryOneBody1D
    @test mild_reference.backend == :numerical_reference
    @test mild_analytic.backend == :pgdg_experimental
    @test mild_localized.backend == :pgdg_localized_experimental
    @test mild_oracle.backend == :pgdg_localized_oracle
    @test occursin("experimental=true", sprint(show, mild_analytic))
    @test occursin("experimental=true", sprint(show, mild_localized))
    @test occursin("experimental=true", sprint(show, mild_oracle))
    @test !occursin("experimental=true", sprint(show, mild_reference))
    @test mild_reference.overlap ≈ transpose(mild_reference.overlap) atol = 1.0e-10 rtol = 1.0e-10
    @test mild_analytic.overlap ≈ transpose(mild_analytic.overlap) atol = 1.0e-10 rtol = 1.0e-10
    @test mild_localized.overlap ≈ transpose(mild_localized.overlap) atol = 1.0e-10 rtol = 1.0e-10
    @test mild_reference.kinetic ≈ transpose(mild_reference.kinetic) atol = 1.0e-10 rtol = 1.0e-10
    @test mild_analytic.kinetic ≈ transpose(mild_analytic.kinetic) atol = 1.0e-10 rtol = 1.0e-10
    @test mild_localized.kinetic ≈ transpose(mild_localized.kinetic) atol = 1.0e-10 rtol = 1.0e-10
    @test length(mild_reference.gaussian_factors) == 3
    @test length(mild_analytic.gaussian_factors) == 3
    @test length(mild_localized.gaussian_factors) == 3
    @test mild_basis isa MappedUniformBasis
    @test norm(mild_reference.overlap - mild_analytic.overlap, Inf) < 0.05
    @test norm(mild_reference.kinetic - mild_analytic.kinetic, Inf) < 0.05
    @test norm(mild_reference.gaussian_factors[1] - mild_analytic.gaussian_factors[1], Inf) < 0.05
    @test norm(mild_localized.overlap - I, Inf) < 1.0e-10
    @test norm(mild_localized.overlap - I, Inf) < norm(mild_analytic.overlap - I, Inf)
    @test norm(mild_localized.overlap - overlap_reference_localized, Inf) < 1.0e-10
    @test norm(mild_localized.kinetic - kinetic_reference_localized, Inf) <
          norm(raw_localized_kinetic - kinetic_reference_localized, Inf)
    @test norm(mild_oracle.kinetic - kinetic_reference_localized, Inf) <
          norm(mild_localized.kinetic - kinetic_reference_localized, Inf)
end

@testset "Public ordinary and nested backend contract" begin
    function _argument_error_text(route_builder)
        err = try
            route_builder()
            nothing
        catch err
            err
        end
        @test err isa ArgumentError
        return sprint(showerror, err)
    end

    function _reference_only_backend_error(route_builder)
        text = _argument_error_text(route_builder)
        @test occursin("numerical-reference-only route", text)
        @test occursin("PGDG production-contract support is not yet implemented here", text)
        return text
    end

    @testset "Structured final one-body mix matches dense congruence" begin
        carried_one_body = [
            1.1 0.2
            0.3 1.7
        ]
        one_body_ga = [
            -0.4 0.6
            0.5 -0.2
        ]
        one_body_aa = [
            0.9 -0.1
            0.25 1.3
        ]
        raw_to_final = [
            1.0 0.0 0.2 -0.1
            0.0 1.0 0.05 0.3
            0.0 0.0 0.7 0.1
            0.0 0.0 -0.2 0.6
        ]
        raw_one_body = [
            carried_one_body one_body_ga
            transpose(one_body_ga) one_body_aa
        ]
        dense_final = Matrix{Float64}(transpose(raw_to_final) * raw_one_body * raw_to_final)
        dense_final = 0.5 .* (dense_final .+ transpose(dense_final))

        structured_raw, structured_final = GaussletBases._qwrg_structured_final_one_body_matrices(
            carried_one_body,
            one_body_ga,
            one_body_aa,
            raw_to_final,
        )
        @test isnothing(structured_raw)
        @test structured_final ≈ dense_final atol = 1.0e-12 rtol = 1.0e-12

        mixed_raw, mixed_final = GaussletBases._qwrg_one_body_matrices(
            carried_one_body,
            one_body_ga,
            one_body_aa,
            raw_to_final,
        )
        @test isnothing(mixed_raw)
        @test mixed_final ≈ dense_final atol = 1.0e-12 rtol = 1.0e-12

        broken_raw_to_final = copy(raw_to_final)
        broken_raw_to_final[3, 1] = 1.0e-3
        @test isnothing(
            GaussletBases._qwrg_structured_final_one_body_matrices(
                carried_one_body,
                one_body_ga,
                one_body_aa,
                broken_raw_to_final,
            ),
        )

        fallback_raw, fallback_final = GaussletBases._qwrg_one_body_matrices(
            carried_one_body,
            one_body_ga,
            one_body_aa,
            broken_raw_to_final,
        )
        fallback_dense_final = Matrix{Float64}(transpose(broken_raw_to_final) * raw_one_body * broken_raw_to_final)
        fallback_dense_final = 0.5 .* (fallback_dense_final .+ transpose(fallback_dense_final))
        @test !isnothing(fallback_raw)
        @test fallback_final ≈ fallback_dense_final atol = 1.0e-12 rtol = 1.0e-12
    end

    function _check_direct_product_backend_pair(
        basis,
        nuclear_charges::AbstractVector{<:Real},
    )
        reference = ordinary_cartesian_qiu_white_operators(
            basis;
            nuclear_charges = nuclear_charges,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :numerical_reference,
        )
        localized = ordinary_cartesian_qiu_white_operators(
            basis;
            nuclear_charges = nuclear_charges,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :pgdg_localized_experimental,
        )
        reference_check = GaussletBases.ordinary_cartesian_1s2_check(reference)
        localized_check = GaussletBases.ordinary_cartesian_1s2_check(localized)

        @test reference.gausslet_backend == :numerical_reference
        @test localized.gausslet_backend == :pgdg_localized_experimental
        @test reference.gaussian_data === nothing
        @test localized.gaussian_data === nothing
        @test reference.residual_count == 0
        @test localized.residual_count == 0
        @test norm(reference.overlap - I, Inf) < 1.0e-8
        @test norm(localized.overlap - I, Inf) < 1.0e-8
        @test localized.one_body_hamiltonian ≈ transpose(localized.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
        @test localized.interaction_matrix ≈ transpose(localized.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
        @test norm(localized.overlap - reference.overlap, Inf) < 0.05
        @test norm(localized.one_body_hamiltonian - reference.one_body_hamiltonian, Inf) < 0.35
        @test norm(localized.interaction_matrix - reference.interaction_matrix, Inf) < 0.3
        @test abs(localized_check.orbital_energy - reference_check.orbital_energy) < 0.05
        @test abs(localized_check.vee_expectation - reference_check.vee_expectation) < 0.05
        return reference, localized, reference_check, localized_check
    end

    function _check_nested_fixed_block_backend_pair(
        fixed_block,
        nuclear_charges::AbstractVector{<:Real},
        expansion::CoulombGaussianExpansion,
    )
        @test basis_representation(fixed_block).metadata.parent_kind == :cartesian_product_basis

        reference = ordinary_cartesian_qiu_white_operators(
            fixed_block;
            nuclear_charges = nuclear_charges,
            expansion = expansion,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :numerical_reference,
        )
        localized = ordinary_cartesian_qiu_white_operators(
            fixed_block;
            nuclear_charges = nuclear_charges,
            expansion = expansion,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :pgdg_localized_experimental,
        )
        reference_check = GaussletBases.ordinary_cartesian_1s2_check(reference)
        localized_check = GaussletBases.ordinary_cartesian_1s2_check(localized)

        @test reference.gausslet_backend == :numerical_reference
        @test localized.gausslet_backend == :pgdg_localized_experimental
        @test reference.gaussian_data === nothing
        @test localized.gaussian_data === nothing
        @test reference.residual_count == 0
        @test localized.residual_count == 0
        @test norm(reference.overlap - I, Inf) < 1.0e-8
        @test norm(localized.overlap - I, Inf) < 1.0e-8
        @test localized.one_body_hamiltonian ≈ transpose(localized.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
        @test localized.interaction_matrix ≈ transpose(localized.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
        @test norm(localized.overlap - reference.overlap, Inf) < 0.05
        @test norm(localized.one_body_hamiltonian - reference.one_body_hamiltonian, Inf) < 0.4
        @test norm(localized.interaction_matrix - reference.interaction_matrix, Inf) < 0.35
        @test abs(localized_check.orbital_energy - reference_check.orbital_energy) < 0.06
        @test abs(localized_check.vee_expectation - reference_check.vee_expectation) < 0.06
        return reference, localized, reference_check, localized_check
    end

    function _check_diatomic_molecular_backend_pair(
        basis,
        supplement,
        nuclear_charges::AbstractVector{<:Real},
    )
        reference = ordinary_cartesian_qiu_white_operators(
            basis,
            supplement;
            nuclear_charges = nuclear_charges,
            nuclear_term_storage = :by_center,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :numerical_reference,
        )
        localized = ordinary_cartesian_qiu_white_operators(
            basis,
            supplement;
            nuclear_charges = nuclear_charges,
            nuclear_term_storage = :by_center,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :pgdg_localized_experimental,
        )
        atom_a_localized = ordinary_cartesian_qiu_white_operators(
            basis,
            supplement;
            nuclear_charges = [1.0, 0.0],
            nuclear_term_storage = :total_only,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :pgdg_localized_experimental,
        )
        atom_b_localized = ordinary_cartesian_qiu_white_operators(
            basis,
            supplement;
            nuclear_charges = [0.0, 1.0],
            nuclear_term_storage = :total_only,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :pgdg_localized_experimental,
        )
        reference_check = GaussletBases.ordinary_cartesian_1s2_check(reference)
        localized_check = GaussletBases.ordinary_cartesian_1s2_check(localized)

        @test reference.gausslet_backend == :numerical_reference
        @test localized.gausslet_backend == :pgdg_localized_experimental
        @test reference.gaussian_data === supplement
        @test localized.gaussian_data === supplement
        @test reference.residual_count > 0
        @test localized.residual_count == reference.residual_count
        @test size(localized.one_body_hamiltonian, 1) == localized.gausslet_count + localized.residual_count
        @test size(localized.one_body_hamiltonian) == size(reference.one_body_hamiltonian)
        @test reference.nuclear_term_storage == :by_center
        @test localized.nuclear_term_storage == :by_center
        @test !isnothing(localized.kinetic_one_body)
        @test !isnothing(localized.nuclear_one_body_by_center)
        @test length(localized.nuclear_one_body_by_center) == length(nuclear_charges)
        @test norm(reference.overlap - I, Inf) < 1.0e-8
        @test norm(localized.overlap - I, Inf) < 1.0e-8
        @test localized.one_body_hamiltonian ≈ transpose(localized.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
        @test localized.interaction_matrix ≈ transpose(localized.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(localized) ≈
              localized.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
        @test assembled_one_body_hamiltonian(localized; nuclear_charges = [1.0, 0.0]) ≈
              atom_a_localized.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(localized; nuclear_charges = [0.0, 1.0]) ≈
              atom_b_localized.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test norm(localized.overlap - reference.overlap, Inf) < 0.05
        @test norm(localized.interaction_matrix - reference.interaction_matrix, Inf) < 0.35
        @test abs(localized_check.orbital_energy - reference_check.orbital_energy) < 0.02
        @test abs(localized_check.vee_expectation - reference_check.vee_expectation) < 0.02
        return reference, localized, reference_check, localized_check
    end

    function _check_nested_diatomic_molecular_backend_pair(
        fixed_block,
        supplement,
        nuclear_charges::AbstractVector{<:Real},
        expansion::CoulombGaussianExpansion,
    )
        @test basis_representation(fixed_block).metadata.parent_kind == :cartesian_product_basis

        reference = ordinary_cartesian_qiu_white_operators(
            fixed_block,
            supplement;
            nuclear_charges = nuclear_charges,
            nuclear_term_storage = :by_center,
            expansion = expansion,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :numerical_reference,
        )
        localized = ordinary_cartesian_qiu_white_operators(
            fixed_block,
            supplement;
            nuclear_charges = nuclear_charges,
            nuclear_term_storage = :by_center,
            expansion = expansion,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :pgdg_localized_experimental,
        )
        atom_a_localized = ordinary_cartesian_qiu_white_operators(
            fixed_block,
            supplement;
            nuclear_charges = [1.0, 0.0],
            nuclear_term_storage = :total_only,
            expansion = expansion,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :pgdg_localized_experimental,
        )
        atom_b_localized = ordinary_cartesian_qiu_white_operators(
            fixed_block,
            supplement;
            nuclear_charges = [0.0, 1.0],
            nuclear_term_storage = :total_only,
            expansion = expansion,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :pgdg_localized_experimental,
        )
        reference_check = GaussletBases.ordinary_cartesian_1s2_check(reference)
        localized_check = GaussletBases.ordinary_cartesian_1s2_check(localized)

        @test reference.gausslet_backend == :numerical_reference
        @test localized.gausslet_backend == :pgdg_localized_experimental
        @test reference.gaussian_data === supplement
        @test localized.gaussian_data === supplement
        @test reference.residual_count > 0
        @test localized.residual_count == reference.residual_count
        @test size(localized.one_body_hamiltonian, 1) == localized.gausslet_count + localized.residual_count
        @test size(localized.one_body_hamiltonian) == size(reference.one_body_hamiltonian)
        @test reference.nuclear_term_storage == :by_center
        @test localized.nuclear_term_storage == :by_center
        @test !isnothing(localized.kinetic_one_body)
        @test !isnothing(localized.nuclear_one_body_by_center)
        @test length(localized.nuclear_one_body_by_center) == length(nuclear_charges)
        @test norm(reference.overlap - I, Inf) < 1.0e-8
        @test norm(localized.overlap - I, Inf) < 1.0e-8
        @test localized.one_body_hamiltonian ≈ transpose(localized.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
        @test localized.interaction_matrix ≈ transpose(localized.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(localized) ≈
              localized.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
        @test assembled_one_body_hamiltonian(localized; nuclear_charges = [1.0, 0.0]) ≈
              atom_a_localized.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(localized; nuclear_charges = [0.0, 1.0]) ≈
              atom_b_localized.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test norm(localized.overlap - reference.overlap, Inf) < 0.05
        @test norm(localized.interaction_matrix - reference.interaction_matrix, Inf) < 0.35
        @test abs(localized_check.orbital_energy - reference_check.orbital_energy) < 0.02
        @test abs(localized_check.vee_expectation - reference_check.vee_expectation) < 0.02
        return reference, localized, reference_check, localized_check
    end

    expansion = coulomb_gaussian_expansion(doacc = false)
    mapped_basis = build_basis(MappedUniformBasisSpec(:G10;
        count = 5,
        mapping = fit_asinh_mapping_for_strength(s = 0.5, npoints = 5, xmax = 6.0),
        reference_spacing = 1.0,
    ))
    diatomic_basis = bond_aligned_homonuclear_qw_basis(
        bond_length = 1.4,
        core_spacing = 0.5,
        xmax_parallel = 6.0,
        xmax_transverse = 4.0,
        bond_axis = :z,
    )
    chain_basis = bond_aligned_homonuclear_chain_qw_basis(
        natoms = 3,
        spacing = 1.2,
        core_spacing = 0.5,
        xmax_parallel = 2.0,
        xmax_transverse = 2.0,
        chain_axis = :z,
    )
    square_basis = axis_aligned_homonuclear_square_lattice_qw_basis(
        n = 2,
        spacing = 1.4,
        core_spacing = 0.5,
        xmax_in_plane = 2.0,
        xmax_transverse = 2.0,
    )
    (
        _nested_diatomic_basis,
        _nested_parent_ops,
        _nested_parent_check,
        nested_expansion,
        _nested_source,
        nested_fixed_block,
        _parent_modes,
        _parent_ground,
        _projected,
        _projected_vee,
        _capture,
        _projected_energy,
    ) = _bond_aligned_diatomic_nested_fixed_block_fixture(; bond_length = 1.4)
    (
        chain_nested_basis,
        _chain_nested_source,
        chain_nested_fixed_block,
        _chain_nested_diagnostics,
    ) = _bond_aligned_homonuclear_chain_nested_fixture(;
        natoms = 3,
        odd_chain_policy = :central_ternary_relaxed,
    )
    (
        square_nested_basis,
        _square_nested_source,
        square_nested_fixed_block,
        _square_nested_diagnostics,
    ) = _axis_aligned_homonuclear_square_lattice_nested_fixture(; n = 2)

    @test mapped_ordinary_one_body_operators(
        mapped_basis;
        exponents = expansion.exponents[1:3],
        backend = :pgdg_experimental,
    ).backend == :pgdg_experimental
    @test mapped_ordinary_one_body_operators(
        mapped_basis;
        exponents = expansion.exponents[1:3],
        backend = :pgdg_localized_experimental,
    ).backend == :pgdg_localized_experimental

    @testset "Bond-aligned normalized build context" begin
        diatomic_context = GaussletBases._normalized_bond_aligned_build_context(diatomic_basis)
        chain_context = GaussletBases._normalized_bond_aligned_build_context(chain_basis)
        square_context = GaussletBases._normalized_bond_aligned_build_context(square_basis)
        nested_diatomic_context = GaussletBases._normalized_bond_aligned_build_context(nested_fixed_block)
        nested_chain_context = GaussletBases._normalized_bond_aligned_build_context(chain_nested_fixed_block)
        nested_square_context = GaussletBases._normalized_bond_aligned_build_context(square_nested_fixed_block)

        @test diatomic_context.basis_family == :bond_aligned_diatomic
        @test chain_context.basis_family == :bond_aligned_homonuclear_chain
        @test square_context.basis_family == :axis_aligned_homonuclear_square_lattice
        @test diatomic_context.carried_space_kind == :direct_product
        @test chain_context.carried_space_kind == :direct_product
        @test square_context.carried_space_kind == :direct_product
        @test diatomic_context.parent_basis === diatomic_basis
        @test diatomic_context.carried === diatomic_basis
        @test diatomic_context.contraction === nothing
        @test diatomic_context.default_nuclear_charges == [1.0, 1.0]
        @test diatomic_context.route_metadata.basis_family == :bond_aligned_diatomic
        @test diatomic_context.capabilities.allowed_interaction_treatments == (:ggt_nearest, :mwg)
        @test diatomic_context.capabilities.allowed_gausslet_backends ==
              (:numerical_reference, :pgdg_localized_experimental)
        @test diatomic_context.capabilities.timing_label == "qwrg.bond_aligned_ordinary.total"

        @test nested_diatomic_context.basis_family == :bond_aligned_diatomic
        @test nested_chain_context.basis_family == :bond_aligned_homonuclear_chain
        @test nested_square_context.basis_family == :axis_aligned_homonuclear_square_lattice
        @test nested_diatomic_context.carried_space_kind == :nested_fixed_block
        @test nested_chain_context.carried_space_kind == :nested_fixed_block
        @test nested_square_context.carried_space_kind == :nested_fixed_block
        @test nested_diatomic_context.parent_basis === nested_fixed_block.parent_basis
        @test nested_diatomic_context.carried === nested_fixed_block
        @test nested_diatomic_context.contraction === nested_fixed_block.coefficient_matrix
        @test nested_diatomic_context.parent_route_metadata.basis_family == :bond_aligned_diatomic
        @test nested_diatomic_context.capabilities.allowed_interaction_treatments == (:ggt_nearest,)
        @test nested_diatomic_context.capabilities.allowed_gausslet_backends ==
              (:numerical_reference, :pgdg_localized_experimental)
        @test nested_diatomic_context.capabilities.localized_parent_kind == :cartesian_product_basis
        @test nested_diatomic_context.capabilities.timing_label == "qwrg.bond_aligned_nested_fixed.total"
    end

    diatomic_reference, diatomic_localized, diatomic_reference_check, diatomic_localized_check =
        _check_direct_product_backend_pair(diatomic_basis, [1.0, 1.0])
    chain_reference, chain_localized, chain_reference_check, chain_localized_check =
        _check_direct_product_backend_pair(
            chain_basis,
            fill(1.0, length(chain_basis.nuclei)),
        )
    square_reference, square_localized, square_reference_check, square_localized_check =
        _check_direct_product_backend_pair(
            square_basis,
            fill(1.0, length(square_basis.nuclei)),
        )

    @test diatomic_localized.gausslet_count == diatomic_reference.gausslet_count
    @test chain_localized.gausslet_count == chain_reference.gausslet_count
    @test square_localized.gausslet_count == square_reference.gausslet_count
    @test isfinite(diatomic_localized_check.orbital_energy)
    @test isfinite(chain_localized_check.orbital_energy)
    @test isfinite(square_localized_check.orbital_energy)
    @test isfinite(diatomic_localized_check.vee_expectation)
    @test isfinite(chain_localized_check.vee_expectation)
    @test isfinite(square_localized_check.vee_expectation)

    nested_diatomic_reference,
    nested_diatomic_localized,
    nested_diatomic_reference_check,
    nested_diatomic_localized_check = _check_nested_fixed_block_backend_pair(
        nested_fixed_block,
        [1.0, 1.0],
        nested_expansion,
    )
    nested_chain_reference,
    nested_chain_localized,
    nested_chain_reference_check,
    nested_chain_localized_check = _check_nested_fixed_block_backend_pair(
        chain_nested_fixed_block,
        fill(1.0, length(chain_nested_basis.nuclei)),
        expansion,
    )
    nested_square_reference,
    nested_square_localized,
    nested_square_reference_check,
    nested_square_localized_check = _check_nested_fixed_block_backend_pair(
        square_nested_fixed_block,
        fill(1.0, length(square_nested_basis.nuclei)),
        expansion,
    )

    @test nested_diatomic_localized.gausslet_count == nested_diatomic_reference.gausslet_count
    @test nested_chain_localized.gausslet_count == nested_chain_reference.gausslet_count
    @test nested_square_localized.gausslet_count == nested_square_reference.gausslet_count
    @test isfinite(nested_diatomic_localized_check.orbital_energy)
    @test isfinite(nested_chain_localized_check.orbital_energy)
    @test isfinite(nested_square_localized_check.orbital_energy)
    @test isfinite(nested_diatomic_localized_check.vee_expectation)
    @test isfinite(nested_chain_localized_check.vee_expectation)
    @test isfinite(nested_square_localized_check.vee_expectation)

    if !_legacy_basisfile_available()
        @test true
    else
        hybrid_basis = diatomic_basis
        hybrid_supplement = legacy_bond_aligned_diatomic_gaussian_supplement(
            "H",
            "cc-pVTZ",
            hybrid_basis.nuclei;
            lmax = 0,
            max_width = 1.0,
        )
        hybrid_fixed_block = nested_fixed_block
        hybrid_direct_context = GaussletBases._normalized_bond_aligned_build_context(
            hybrid_basis,
            hybrid_supplement,
        )
        hybrid_nested_context = GaussletBases._normalized_bond_aligned_build_context(
            hybrid_fixed_block,
            hybrid_supplement,
        )

        @test hybrid_direct_context.basis_family == :bond_aligned_diatomic
        @test hybrid_direct_context.carried_space_kind == :direct_product
        @test hybrid_direct_context.parent_basis === hybrid_basis
        @test hybrid_direct_context.carried === hybrid_basis
        @test hybrid_direct_context.gaussian_data === hybrid_supplement
        @test hybrid_direct_context.contraction === nothing
        @test hybrid_direct_context.capabilities.allowed_interaction_treatments == (:ggt_nearest,)
        @test hybrid_direct_context.capabilities.allowed_gausslet_backends ==
              (:numerical_reference, :pgdg_localized_experimental)
        @test hybrid_direct_context.capabilities.timing_label == "qwrg.diatomic_shell.total"

        @test hybrid_nested_context.basis_family == :bond_aligned_diatomic
        @test hybrid_nested_context.carried_space_kind == :nested_fixed_block
        @test hybrid_nested_context.parent_basis === hybrid_fixed_block.parent_basis
        @test hybrid_nested_context.carried === hybrid_fixed_block
        @test hybrid_nested_context.gaussian_data === hybrid_supplement
        @test hybrid_nested_context.contraction === hybrid_fixed_block.coefficient_matrix
        @test hybrid_nested_context.capabilities.allowed_interaction_treatments == (:ggt_nearest,)
        @test hybrid_nested_context.capabilities.allowed_gausslet_backends ==
              (:numerical_reference, :pgdg_localized_experimental)
        @test hybrid_nested_context.capabilities.localized_parent_kind == :cartesian_product_basis
        @test hybrid_nested_context.capabilities.timing_label == "qwrg.nested_diatomic_shell.total"

        hybrid_reference, hybrid_localized, hybrid_reference_check, hybrid_localized_check =
            _check_diatomic_molecular_backend_pair(
                hybrid_basis,
                hybrid_supplement,
                [1.0, 1.0],
            )
        hybrid_nested_reference,
        hybrid_nested_localized,
        hybrid_nested_reference_check,
        hybrid_nested_localized_check = _check_nested_diatomic_molecular_backend_pair(
            hybrid_fixed_block,
            hybrid_supplement,
            [1.0, 1.0],
            nested_expansion,
        )

        @test hybrid_localized.gausslet_backend == :pgdg_localized_experimental
        @test hybrid_localized.gaussian_data === hybrid_supplement
        @test hybrid_nested_localized.gausslet_backend == :pgdg_localized_experimental
        @test hybrid_nested_localized.gaussian_data === hybrid_supplement
        @test hybrid_localized.residual_count == hybrid_reference.residual_count
        @test hybrid_nested_localized.residual_count == hybrid_nested_reference.residual_count
        @test size(hybrid_localized.raw_to_final, 2) ==
              hybrid_localized.gausslet_count + hybrid_localized.residual_count
        @test size(hybrid_nested_localized.raw_to_final, 2) ==
              hybrid_nested_localized.gausslet_count + hybrid_nested_localized.residual_count
        @test isfinite(hybrid_localized_check.orbital_energy)
        @test isfinite(hybrid_localized_check.vee_expectation)
        @test isfinite(hybrid_nested_localized_check.orbital_energy)
        @test isfinite(hybrid_nested_localized_check.vee_expectation)

        hybrid_direct_error_text = _argument_error_text(() ->
            ordinary_cartesian_qiu_white_operators(
                hybrid_basis,
                hybrid_supplement;
                nuclear_charges = [1.0, 1.0],
                interaction_treatment = :ggt_nearest,
                gausslet_backend = :pgdg_experimental,
            )
        )
        @test occursin("bond-aligned diatomic molecular QW path", hybrid_direct_error_text)
        @test occursin(":pgdg_localized_experimental", hybrid_direct_error_text)

        hybrid_nested_error_text = _argument_error_text(() ->
            ordinary_cartesian_qiu_white_operators(
                hybrid_fixed_block,
                hybrid_supplement;
                nuclear_charges = [1.0, 1.0],
                expansion = nested_expansion,
                interaction_treatment = :ggt_nearest,
                gausslet_backend = :pgdg_experimental,
            )
        )
        @test occursin("bond-aligned diatomic nested molecular QW path", hybrid_nested_error_text)
        @test occursin(":pgdg_localized_experimental", hybrid_nested_error_text)

        mismatched_nuclei = [
            (hybrid_basis.nuclei[1][1] - 0.1, hybrid_basis.nuclei[1][2], hybrid_basis.nuclei[1][3]),
            (hybrid_basis.nuclei[2][1] + 0.1, hybrid_basis.nuclei[2][2], hybrid_basis.nuclei[2][3]),
        ]
        mismatched_supplement = legacy_bond_aligned_diatomic_gaussian_supplement(
            "H",
            "cc-pVTZ",
            mismatched_nuclei;
            lmax = 0,
            max_width = 1.0,
        )
        mismatched_nuclei_text = _argument_error_text(() ->
            ordinary_cartesian_qiu_white_operators(
                hybrid_basis,
                mismatched_supplement;
                nuclear_charges = [1.0, 1.0],
                interaction_treatment = :ggt_nearest,
                gausslet_backend = :numerical_reference,
            )
        )
        @test occursin(
            "bond-aligned diatomic molecular supplement nuclei must match the bond-aligned basis nuclei",
            mismatched_nuclei_text,
        )

        function _with_parent_kind(
            context::typeof(hybrid_nested_context),
            parent_kind::Symbol,
        )
            representation = context.carried_representation
            metadata = representation.metadata
            shifted_metadata = GaussletBases.CartesianBasisMetadata3D(
                metadata.basis_kind,
                metadata.axis_sharing,
                metadata.axis_metadata,
                parent_kind,
                metadata.parent_axis_counts,
                metadata.parent_dimension,
                metadata.final_dimension,
                metadata.working_box,
                metadata.basis_labels,
                metadata.basis_centers,
                metadata.route_metadata,
            )
            shifted_representation = GaussletBases.CartesianBasisRepresentation3D(
                shifted_metadata,
                representation.axis_representations,
                representation.contraction_kind,
                representation.coefficient_matrix,
                representation.parent_labels,
                representation.parent_centers,
                representation.support_indices,
                representation.support_states,
                representation.parent_data,
            )
            return typeof(context)(
                context.basis_family,
                context.carried_space_kind,
                context.parent_basis,
                context.carried,
                shifted_representation,
                context.parent_representation,
                context.nuclei,
                context.default_nuclear_charges,
                context.contraction,
                context.gaussian_data,
                context.route_metadata,
                context.parent_route_metadata,
                context.capabilities,
            )
        end

        transformed_parent_kind_text = _argument_error_text(() ->
            GaussletBases._validate_operator_route_backend(
                _with_parent_kind(hybrid_nested_context, :cartesian_plus_supplement_raw),
                :pgdg_localized_experimental,
            )
        )
        @test occursin("bond-aligned diatomic nested molecular QW path", transformed_parent_kind_text)
        @test occursin("parent_kind = :cartesian_product_basis", transformed_parent_kind_text)
        @test occursin(
            "transformed parent spaces remain numerical-reference-only",
            transformed_parent_kind_text,
        )
    end

    direct_product_error_text = _argument_error_text(() ->
        ordinary_cartesian_qiu_white_operators(
            diatomic_basis;
            nuclear_charges = [1.0, 1.0],
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :pgdg_experimental,
        )
    )
    @test occursin("bond-aligned ordinary_cartesian_qiu_white_operators", direct_product_error_text)
    @test occursin(":pgdg_localized_experimental", direct_product_error_text)

    direct_product_interaction_text = _argument_error_text(() ->
        ordinary_cartesian_qiu_white_operators(
            diatomic_basis;
            nuclear_charges = [1.0, 1.0],
            interaction_treatment = :bogus,
            gausslet_backend = :numerical_reference,
        )
    )
    @test occursin(
        "bond-aligned ordinary_cartesian_qiu_white_operators",
        direct_product_interaction_text,
    )
    @test occursin(":ggt_nearest or :mwg", direct_product_interaction_text)

    nested_fixed_block_error_text = _argument_error_text(() ->
        ordinary_cartesian_qiu_white_operators(
            nested_fixed_block;
            nuclear_charges = [1.0, 1.0],
            expansion = nested_expansion,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :pgdg_experimental,
        )
    )
    @test occursin("bond-aligned diatomic nested ordinary_cartesian_qiu_white_operators", nested_fixed_block_error_text)
    @test occursin(":pgdg_localized_experimental", nested_fixed_block_error_text)
    @test occursin("pure Cartesian-parent nested fixed blocks", nested_fixed_block_error_text)

    nested_fixed_block_interaction_text = _argument_error_text(() ->
        ordinary_cartesian_qiu_white_operators(
            nested_fixed_block;
            nuclear_charges = [1.0, 1.0],
            expansion = nested_expansion,
            interaction_treatment = :mwg,
            gausslet_backend = :numerical_reference,
        )
    )
    @test occursin(
        "bond-aligned diatomic nested ordinary_cartesian_qiu_white_operators",
        nested_fixed_block_interaction_text,
    )
    @test occursin(":ggt_nearest", nested_fixed_block_interaction_text)

    diatomic_nested_source_text = _reference_only_backend_error(() ->
        bond_aligned_diatomic_nested_fixed_source(
            diatomic_basis;
            expansion = expansion,
            gausslet_backend = :pgdg_localized_experimental,
        )
    )
    @test occursin("bond-aligned diatomic nested fixed source", diatomic_nested_source_text)

    diatomic_nested_fixed_text = _reference_only_backend_error(() ->
        bond_aligned_diatomic_nested_fixed_block(
            diatomic_basis;
            expansion = expansion,
            gausslet_backend = :pgdg_experimental,
        )
    )
    @test occursin("bond-aligned diatomic nested fixed source", diatomic_nested_fixed_text)

    chain_text = _reference_only_backend_error(() ->
        experimental_bond_aligned_homonuclear_chain_nested_qw_operators(
            chain_basis;
            nuclear_charges = fill(1.0, length(chain_basis.nuclei)),
            gausslet_backend = :pgdg_experimental,
        )
    )
    @test occursin("bond-aligned homonuclear chain nested fixed source", chain_text)

    square_text = _reference_only_backend_error(() ->
        experimental_axis_aligned_homonuclear_square_lattice_nested_qw_operators(
            square_basis;
            nuclear_charges = fill(1.0, length(square_basis.nuclei)),
            gausslet_backend = :pgdg_localized_experimental,
        )
    )
    @test occursin("axis-aligned homonuclear square-lattice nested fixed source", square_text)

    diatomic_nested_context = GaussletBases._normalized_nested_source_frontend_context(
        diatomic_basis;
        expansion = expansion,
        nside = 5,
    )
    chain_nested_context = GaussletBases._normalized_nested_source_frontend_context(
        chain_basis;
        expansion = expansion,
        nside = 5,
        odd_chain_policy = :central_ternary_relaxed,
    )
    square_nested_context = GaussletBases._normalized_nested_source_frontend_context(
        square_basis;
        expansion = expansion,
        nside = 5,
        min_in_plane_aspect_ratio = 0.15,
    )
    @test diatomic_nested_context.capabilities.route_label ==
          "bond-aligned diatomic nested fixed source"
    @test diatomic_nested_context.capabilities.total_timing_label == "diatomic.fixed_source.total"
    @test diatomic_nested_context.build_options.nside == 5
    @test chain_nested_context.capabilities.route_label ==
          "bond-aligned homonuclear chain nested fixed source"
    @test chain_nested_context.build_options.odd_chain_policy == :central_ternary_relaxed
    @test isnothing(chain_nested_context.capabilities.total_timing_label)
    @test square_nested_context.capabilities.route_label ==
          "axis-aligned homonuclear square-lattice nested fixed source"
    @test square_nested_context.build_options.min_in_plane_aspect_ratio == 0.15
    @test isnothing(square_nested_context.capabilities.total_timing_label)

    if !_legacy_basisfile_available()
        @test true
    else
        supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 0)
        (
            _source_basis_nested,
            _bundle_nested,
            _shell_nested,
            _fixed_block_nested,
            _shell_plus_core_nested,
            fixed_block_shell_plus_core,
            _legacy_nested,
            _baseline_nested,
            _nested_shell_only,
            _nested_shell_plus_core,
            _baseline_nested_check,
            _nested_shell_only_check,
            _nested_shell_plus_core_check,
        ) = _nested_qiu_white_nearest_fixture()
        atomic_direct_context = GaussletBases._normalized_atomic_build_context(
            mapped_basis,
            supplement,
        )
        atomic_nested_context = GaussletBases._normalized_atomic_build_context(
            fixed_block_shell_plus_core,
            supplement,
        )
        @test atomic_direct_context.carried_space_kind == :direct_product
        @test atomic_direct_context.parent_basis === mapped_basis
        @test atomic_direct_context.carried === mapped_basis
        @test atomic_direct_context.gaussian_data === supplement
        @test atomic_direct_context.contraction === nothing
        @test atomic_direct_context.capabilities.allowed_gausslet_backends ==
              (:numerical_reference,)
        @test atomic_direct_context.capabilities.allowed_interaction_treatments == (:ggt_nearest, :mwg)
        @test atomic_direct_context.capabilities.timing_label == "qwrg.atomic_shell.total"
        @test atomic_nested_context.carried_space_kind == :nested_fixed_block
        @test atomic_nested_context.parent_basis === fixed_block_shell_plus_core.parent_basis
        @test atomic_nested_context.carried === fixed_block_shell_plus_core
        @test atomic_nested_context.gaussian_data === supplement
        @test atomic_nested_context.contraction === fixed_block_shell_plus_core.coefficient_matrix
        @test atomic_nested_context.capabilities.allowed_gausslet_backends ==
              (:numerical_reference,)
        @test atomic_nested_context.capabilities.allowed_interaction_treatments == (:ggt_nearest,)
        @test atomic_nested_context.capabilities.timing_label == "qwrg.nested_atomic_shell.total"
        supplement_text = _reference_only_backend_error(() ->
            ordinary_cartesian_qiu_white_operators(
                mapped_basis,
                supplement;
                expansion = expansion,
                Z = 2.0,
                interaction_treatment = :ggt_nearest,
                gausslet_backend = :pgdg_localized_experimental,
            )
        )
        @test occursin("ordinary_cartesian_qiu_white_operators", supplement_text)
        direct_atomic_bad_interaction_text = _argument_error_text(() ->
            ordinary_cartesian_qiu_white_operators(
                mapped_basis,
                supplement;
                expansion = expansion,
                Z = 2.0,
                interaction_treatment = :bogus,
            )
        )
        @test occursin(
            "Qiu-White interaction_treatment must be :ggt_nearest or :mwg",
            direct_atomic_bad_interaction_text,
        )
        nested_supplement_text = _reference_only_backend_error(() ->
            ordinary_cartesian_qiu_white_operators(
                fixed_block_shell_plus_core,
                supplement;
                expansion = expansion,
                Z = 2.0,
                interaction_treatment = :ggt_nearest,
                gausslet_backend = :pgdg_localized_experimental,
            )
        )
        @test occursin("nested ordinary_cartesian_qiu_white_operators", nested_supplement_text)
        nested_atomic_mwg_text = _argument_error_text(() ->
            ordinary_cartesian_qiu_white_operators(
                fixed_block_shell_plus_core,
                supplement;
                expansion = expansion,
                Z = 2.0,
                interaction_treatment = :mwg,
            )
        )
        @test occursin(
            "nested ordinary_cartesian_qiu_white_operators currently supports only interaction_treatment = :ggt_nearest",
            nested_atomic_mwg_text,
        )
    end
end

@testset "Diatomic molecular one-body timing labels" begin
    function _capture_stdout_text(f::Function)
        path, io = mktemp()
        close(io)
        try
            open(path, "w") do stream
                redirect_stdout(stream) do
                    f()
                end
            end
            return read(path, String)
        finally
            rm(path; force = true)
        end
    end

    if !_legacy_basisfile_available()
        @test true
    else
        (
            basis,
            _operators,
            _check,
            expansion,
            _source,
            fixed_block,
            _parent_modes,
            _parent_ground,
            _projected,
            _projected_vee,
            _capture,
            _projected_energy,
        ) = _bond_aligned_diatomic_nested_fixed_block_fixture(; bond_length = 1.4)

        supplement = legacy_bond_aligned_diatomic_gaussian_supplement(
            "H",
            "cc-pVTZ",
            basis.nuclei;
            lmax = 0,
            max_width = 1.0,
        )

        direct_output = _capture_stdout_text(() ->
            ordinary_cartesian_qiu_white_operators(
                basis,
                supplement;
                nuclear_charges = [1.0, 1.0],
                nuclear_term_storage = :by_center,
                interaction_treatment = :ggt_nearest,
                gausslet_backend = :pgdg_localized_experimental,
                timing = true,
            ),
        )
        @test occursin("qwrg.diatomic_shell.one_body.carried", direct_output)
        @test occursin("qwrg.diatomic_shell.one_body.coupling", direct_output)
        @test occursin("qwrg.diatomic_shell.one_body.supplement", direct_output)
        @test occursin("qwrg.diatomic_shell.one_body.final_mix", direct_output)
        @test occursin("qwrg.diatomic_shell.one_body.by_center_final_mix", direct_output)

        eager_factorized_basis = fixed_block.factorized_cartesian_parent_basis[]
        @test !isnothing(eager_factorized_basis)
        eager_representation = basis_representation(fixed_block)
        @test hasproperty(eager_representation.parent_data, :factorized_cartesian_parent_basis)
        @test eager_representation.parent_data.factorized_cartesian_parent_basis ===
              eager_factorized_basis

        nested_output = _capture_stdout_text(() -> begin
            fixed_block.factorized_cartesian_parent_basis[] = nothing
            ordinary_cartesian_qiu_white_operators(
                fixed_block,
                supplement;
                nuclear_charges = [1.0, 1.0],
                nuclear_term_storage = :by_center,
                expansion = expansion,
                interaction_treatment = :ggt_nearest,
                gausslet_backend = :pgdg_localized_experimental,
                timing = true,
            )
        end)
        @test occursin("qwrg.nested_diatomic_shell.one_body.carried", nested_output)
        @test occursin("qwrg.nested_diatomic_shell.one_body.carried.kinetic", nested_output)
        @test occursin("qwrg.nested_diatomic_shell.one_body.carried.factorized_basis", nested_output)
        @test occursin("qwrg.nested_diatomic_shell.one_body.carried.nuclear_setup", nested_output)
        @test occursin("qwrg.nested_diatomic_shell.one_body.carried.nuclear_contract", nested_output)
        @test occursin("qwrg.nested_diatomic_shell.one_body.coupling", nested_output)
        @test occursin("qwrg.nested_diatomic_shell.one_body.supplement", nested_output)
        @test occursin("qwrg.nested_diatomic_shell.one_body.final_mix", nested_output)
        @test occursin("qwrg.nested_diatomic_shell.one_body.by_center_final_mix", nested_output)
        @test !isnothing(fixed_block.factorized_cartesian_parent_basis[])
        @test fixed_block.factorized_cartesian_parent_basis[].basis_triplets ==
              eager_factorized_basis.basis_triplets
        @test fixed_block.factorized_cartesian_parent_basis[].basis_amplitudes ≈
              eager_factorized_basis.basis_amplitudes atol = 1.0e-12 rtol = 1.0e-12
        nested_representation = basis_representation(fixed_block)
        @test hasproperty(nested_representation.parent_data, :factorized_cartesian_parent_basis)
        @test nested_representation.parent_data.factorized_cartesian_parent_basis ===
              fixed_block.factorized_cartesian_parent_basis[]
    end
end

@testset "Ordinary Cartesian IDA operators" begin
    mild_basis, mild_expansion, mild_analytic = _quick_ordinary_cartesian_ida_fixture(
        backend = :pgdg_experimental,
        mapped = true,
        s = 0.5,
    )
    (_, _, identity_analytic) = _quick_ordinary_cartesian_ida_fixture(
        backend = :pgdg_experimental,
        mapped = false,
    )

    function reconstruct_interaction(expansion::CoulombGaussianExpansion, factors::AbstractVector{<:AbstractMatrix})
        interaction = zeros(Float64, size(first(factors), 1)^3, size(first(factors), 1)^3)
        for term in eachindex(expansion.coefficients)
            factor = factors[term]
            interaction .+= expansion.coefficients[term] .* kron(factor, kron(factor, factor))
        end
        return interaction
    end

    @test mild_basis isa MappedUniformBasis
    @test mild_analytic isa OrdinaryCartesianIDAOperators
    @test identity_analytic isa OrdinaryCartesianIDAOperators
    @test mild_analytic.backend == :pgdg_experimental
    @test identity_analytic.backend == :pgdg_experimental
    @test occursin("experimental=true", sprint(show, mild_analytic))
    @test occursin("experimental=true", sprint(show, identity_analytic))
    @test length(orbitals(mild_analytic)) == length(mild_basis)^3
    @test orbitals(mild_analytic)[1].ix == 1
    @test orbitals(mild_analytic)[1].iy == 1
    @test orbitals(mild_analytic)[1].iz == 1
    @test orbitals(mild_analytic)[end].ix == length(mild_basis)
    @test orbitals(mild_analytic)[end].iy == length(mild_basis)
    @test orbitals(mild_analytic)[end].iz == length(mild_basis)
    @test mild_analytic.overlap_3d ≈ transpose(mild_analytic.overlap_3d) atol = 1.0e-10 rtol = 1.0e-10
    @test mild_analytic.one_body_hamiltonian ≈ transpose(mild_analytic.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test mild_analytic.interaction_matrix ≈ transpose(mild_analytic.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test all(
        factor -> isapprox(factor, transpose(factor); atol = 1.0e-10, rtol = 1.0e-10),
        mild_analytic.pair_factors_1d,
    )
    @test mild_analytic.interaction_matrix ≈ reconstruct_interaction(mild_expansion, mild_analytic.pair_factors_1d) atol = 1.0e-10 rtol = 1.0e-10
    @test minimum(diag(mild_analytic.interaction_matrix)) > 0.0
    @test minimum(diag(identity_analytic.interaction_matrix)) > 0.0
    @test size(identity_analytic.one_body_hamiltonian) == (length(identity_analytic.orbital_data), length(identity_analytic.orbital_data))
    @test opnorm(mild_analytic.interaction_matrix, Inf) > 0.0
end

@testset "Ordinary Cartesian 1s^2 Vee check" begin
    basis, operators, orbital_energy, orbital, vee = _quick_ordinary_cartesian_1s2_vee_fixture()
    reference_value = 1.25

    @test basis isa MappedUniformBasis
    @test operators isa OrdinaryCartesianIDAOperators
    @test operators.backend == :numerical_reference
    @test norm(operators.overlap_3d - I, Inf) < 1.0e-10
    @test operators.interaction_matrix ≈ transpose(operators.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test isfinite(orbital_energy)
    @test orbital_energy < 0.0
    @test isfinite(vee)
    @test vee > 0.0
    @test abs(sum(abs2, orbital) - 1.0) < 1.0e-10
    @test abs(vee - reference_value) < 0.02
    @test ordinary_cartesian_vee_expectation(operators, 2.0 .* orbital) ≈ vee atol = 1.0e-12 rtol = 0.0
end

@testset "Legacy 1D hybrid ordinary Cartesian 1s^2 Vee check" begin
    (
        source_basis,
        core_gaussians,
        pure_operators,
        hybrid_basis,
        hybrid_operators,
        pure_check,
        hybrid_check,
    ) = _quick_hybrid_cartesian_1s2_vee_fixture()
    reference_value = 1.25

    @test source_basis isa MappedUniformBasis
    @test hybrid_basis isa HybridMappedOrdinaryBasis1D
    @test hybrid_operators isa OrdinaryCartesianIDAOperators
    @test hybrid_operators.backend == :pgdg_localized_experimental
    @test length(core_gaussians) == 2
    @test length(hybrid_basis) == length(source_basis) + length(core_gaussians)
    @test length(orbitals(hybrid_operators)) == length(hybrid_basis)^3
    @test norm(hybrid_operators.overlap_3d - I, Inf) < 1.0e-10
    @test hybrid_operators.one_body_hamiltonian ≈ transpose(hybrid_operators.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test hybrid_operators.interaction_matrix ≈ transpose(hybrid_operators.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test isfinite(pure_check.orbital_energy)
    @test isfinite(hybrid_check.orbital_energy)
    @test isfinite(pure_check.vee_expectation)
    @test isfinite(hybrid_check.vee_expectation)
    @test pure_check.vee_expectation > 0.0
    @test hybrid_check.vee_expectation > 0.0
    @test hybrid_check.vee_expectation > pure_check.vee_expectation
    @test abs(hybrid_check.vee_expectation - reference_value) <
          abs(pure_check.vee_expectation - reference_value)
    @test abs(hybrid_check.orbital_energy + 2.0) < abs(pure_check.orbital_energy + 2.0)
@test ordinary_cartesian_vee_expectation(hybrid_operators, 3.0 .* hybrid_check.orbital) ≈
          hybrid_check.vee_expectation atol = 1.0e-12 rtol = 0.0
end

@testset "Legacy 1D hybrid residual Gaussian interaction treatment" begin
    (
        source_basis,
        hybrid_basis,
        pure_operators,
        combined_operators,
        residual_operators,
        pure_check,
        combined_check,
        residual_check,
    ) = _friendly_hybrid_residual_vee_fixture(11, 0.6)
    reference_value = 1.25

    @test source_basis isa MappedUniformBasis
    @test hybrid_basis isa HybridMappedOrdinaryBasis1D
    @test combined_operators.interaction_treatment == :combined_basis
    @test residual_operators.interaction_treatment == :residual_gaussian_nearest
    @test norm(combined_operators.overlap_3d - residual_operators.overlap_3d, Inf) < 1.0e-12
    @test norm(combined_operators.one_body_hamiltonian - residual_operators.one_body_hamiltonian, Inf) < 1.0e-12
    @test combined_operators.interaction_matrix ≈ transpose(combined_operators.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test residual_operators.interaction_matrix ≈ transpose(residual_operators.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test pure_check.vee_expectation > 0.0
    @test combined_check.vee_expectation > 0.0
    @test residual_check.vee_expectation > 0.0
    @test abs(combined_check.vee_expectation - residual_check.vee_expectation) > 1.0e-2
    @test abs(combined_check.vee_expectation - reference_value) <
          abs(residual_check.vee_expectation - reference_value)
    @test abs(combined_check.orbital_energy - residual_check.orbital_energy) < 1.0e-12
    @test abs(pure_check.vee_expectation - reference_value) < 0.02
    @test abs(combined_check.orbital_energy + 2.0) < abs(pure_check.orbital_energy + 2.0)
end

@testset "Legacy atomic Gaussian supplement" begin
    if !_legacy_basisfile_available()
        @test true
    else
        vtz0 = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 0)
        vtz1 = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 1)
        vtz = legacy_s_gaussian_data("He", "cc-pVTZ")
        vqz = legacy_s_gaussian_data("He", "cc-pVQZ")
        vtz_uncontracted = legacy_s_gaussian_data("He", "cc-pVTZ"; uncontracted = true)

        @test vtz0 isa LegacyAtomicGaussianSupplement
        @test vtz isa LegacySGaussianData
        @test vqz isa LegacySGaussianData
        @test vtz0.lmax == 0
        @test vtz1.lmax == 1
        @test length(vtz0.shells) == 3
        @test any(shell -> shell.l == 1, vtz1.shells)
        @test length(vtz1.shells) > length(vtz0.shells)
        @test vtz.primitive_exponents == vtz0.primitive_exponents
        @test vtz.primitive_widths == vtz0.primitive_widths
        @test vtz.contraction_matrix ≈ vtz0.contraction_matrix atol = 0.0 rtol = 0.0
        @test vtz.widths ≈ vtz0.widths atol = 0.0 rtol = 0.0
        @test vtz1.primitive_exponents == vtz0.primitive_exponents
        @test vtz1.primitive_widths == vtz0.primitive_widths
        @test vtz1.contraction_matrix ≈ vtz0.contraction_matrix atol = 0.0 rtol = 0.0
        @test vtz1.widths ≈ vtz0.widths atol = 0.0 rtol = 0.0
        @test !vtz.uncontracted
        @test !vqz.uncontracted
        @test vtz.max_width === nothing
        @test vqz.max_width === nothing
        @test vtz.primitive_exponents == [234.0, 35.16, 7.989, 2.212, 0.6669, 0.2089]
        @test vqz.primitive_exponents == [528.5, 79.31, 18.05, 5.085, 1.609, 0.5363, 0.1833]
        @test length(vtz.primitive_gaussians) == 6
        @test length(vqz.primitive_gaussians) == 7
        @test size(vtz.contraction_matrix) == (6, 3)
        @test size(vqz.contraction_matrix) == (7, 4)
        @test length(vtz.gaussians) == 3
        @test length(vqz.gaussians) == 4
        @test all(isfinite, vtz.widths)
        @test all(isfinite, vqz.widths)
        @test all(>(0.0), vtz.widths)
        @test all(>(0.0), vqz.widths)
        @test all(abs(gaussian.center_value) < 1.0e-12 for gaussian in vtz.gaussians)
        @test all(abs(gaussian.center_value) < 1.0e-12 for gaussian in vqz.gaussians)
        @test vtz.primitive_widths[1] ≈ inv(sqrt(2 * 234.0)) atol = 1.0e-12 rtol = 0.0
        @test vqz.primitive_widths[end] ≈ inv(sqrt(2 * 0.1833)) atol = 1.0e-12 rtol = 0.0
        @test vtz_uncontracted.uncontracted
        @test size(vtz_uncontracted.contraction_matrix) == (6, 6)
        @test vtz_uncontracted.contraction_matrix ≈ Matrix{Float64}(I, 6, 6) atol = 0.0 rtol = 0.0
        @test length(vtz_uncontracted.gaussians) == 6
    end
end

@testset "Legacy bond-aligned diatomic Gaussian supplement width cutoff" begin
    if !_legacy_basisfile_available()
        @test true
    else
        nuclei = [(-0.7, 0.0, 0.0), (0.7, 0.0, 0.0)]
        full = legacy_bond_aligned_diatomic_gaussian_supplement(
            "H",
            "cc-pVTZ",
            nuclei;
            lmax = 1,
        )
        trimmed = legacy_bond_aligned_diatomic_gaussian_supplement(
            "H",
            "cc-pVTZ",
            nuclei;
            lmax = 1,
            max_width = 1.0,
        )
        hetero_full = legacy_bond_aligned_heteronuclear_gaussian_supplement(
            "He",
            "cc-pVTZ",
            "H",
            "cc-pVTZ",
            nuclei;
            lmax = 1,
        )
        hetero_trimmed = legacy_bond_aligned_heteronuclear_gaussian_supplement(
            "He",
            "cc-pVTZ",
            "H",
            "cc-pVTZ",
            nuclei;
            lmax = 1,
            max_width = 1.0,
        )
        full_cart = GaussletBases._bond_aligned_diatomic_cartesian_shell_supplement_3d(full)
        trimmed_cart = GaussletBases._bond_aligned_diatomic_cartesian_shell_supplement_3d(trimmed)
        hetero_full_cart =
            GaussletBases._bond_aligned_diatomic_cartesian_shell_supplement_3d(hetero_full)
        hetero_trimmed_cart =
            GaussletBases._bond_aligned_diatomic_cartesian_shell_supplement_3d(hetero_trimmed)

        @test full.max_width === nothing
        @test trimmed.max_width == 1.0
        @test hetero_full.max_width === nothing
        @test hetero_trimmed.max_width == 1.0
        @test full.atomic_source.max_width === nothing
        @test trimmed.atomic_source.max_width === nothing
        @test length(full_cart.orbitals) == 18
        @test length(trimmed_cart.orbitals) == 8
        @test length(hetero_trimmed_cart.orbitals) < length(hetero_full_cart.orbitals)
        @test all(
            all(width <= 1.0 for width in 1.0 ./ sqrt.(2.0 .* orbital.exponents))
            for orbital in trimmed_cart.orbitals
        )
        @test all(
            all(width <= 1.0 for width in 1.0 ./ sqrt.(2.0 .* orbital.exponents))
            for orbital in hetero_trimmed_cart.orbitals
        )
    end

    mktemp() do path, io
        write(
            io,
            "#BASIS SET: He repo-wide-contraction\n" *
            "He    S\n" *
            "      8.0000000              1.0000000\n" *
            "      0.1250000              1.0000000\n" *
            "END\n",
        )
        close(io)

        nuclei = [(-0.7, 0.0, 0.0), (0.7, 0.0, 0.0)]
        trimmed = legacy_bond_aligned_diatomic_gaussian_supplement(
            "He",
            "repo-wide-contraction",
            nuclei;
            lmax = 0,
            basisfile = path,
            max_width = 1.0,
        )
        removed = legacy_bond_aligned_diatomic_gaussian_supplement(
            "He",
            "repo-wide-contraction",
            nuclei;
            lmax = 0,
            basisfile = path,
            max_width = 0.2,
        )
        trimmed_cart = GaussletBases._bond_aligned_diatomic_cartesian_shell_supplement_3d(trimmed)
        removed_cart = GaussletBases._bond_aligned_diatomic_cartesian_shell_supplement_3d(removed)

        @test length(trimmed_cart.orbitals) == 2
        @test all(orbital.exponents == [8.0] for orbital in trimmed_cart.orbitals)
        @test all(orbital.coefficients == [1.0] for orbital in trimmed_cart.orbitals)
        @test isempty(removed_cart.orbitals)
    end
end

@testset "Public GTO overlap and occupancy matrix helpers" begin
    mktemp() do path, io
        write(
            io,
            "#BASIS SET: He repo-sp\n" *
            "He    S\n" *
            "      1.0000000              1.0000000\n" *
            "He    P\n" *
            "      0.8000000              1.0000000\n" *
            "END\n" *
            "#BASIS SET: H repo-sp\n" *
            "H    S\n" *
            "      1.1000000              1.0000000\n" *
            "H    P\n" *
            "      0.7000000              1.0000000\n" *
            "END\n",
        )
        close(io)

        basis = build_basis(MappedUniformBasisSpec(:G10;
            count = 5,
            mapping = fit_asinh_mapping_for_strength(s = 0.5, npoints = 5, xmax = 4.0),
            reference_spacing = 1.0,
        ))
        supplement = legacy_atomic_gaussian_supplement(
            "He",
            "repo-sp";
            lmax = 1,
            basisfile = path,
        )
        probe_representation = basis_representation(supplement)
        working_representation = GaussletBases._cartesian_direct_product_representation(basis)
        overlap = gto_overlap_matrix(basis, supplement)
        reference_overlap = GaussletBases._cartesian_basis_supplement_cross(
            working_representation,
            probe_representation,
        )

        @test probe_representation isa CartesianGaussianShellSupplementRepresentation3D
        @test basis_representation(probe_representation) === probe_representation
        @test basis_metadata(probe_representation) == probe_representation.metadata
        @test size(overlap) == (length(basis)^3, length(probe_representation.orbitals))
        @test norm(overlap - reference_overlap, Inf) < 1.0e-12
        @test gto_overlap_matrix(working_representation, probe_representation) ≈ overlap atol = 1.0e-12 rtol = 1.0e-12

        block_indices = [1, 7, size(overlap, 1)]
        block_overlap = gto_overlap_matrix(basis, supplement; block_indices = block_indices)
        @test block_overlap == overlap[block_indices, :]
        @test gto_overlap_matrix(basis, supplement, block_indices) == block_overlap

        uniform_occupancy = gto_occupancy_matrix(basis, supplement)
        @test uniform_occupancy ≈ overlap * transpose(overlap) atol = 1.0e-12 rtol = 1.0e-12
        @test uniform_occupancy ≈ transpose(uniform_occupancy) atol = 1.0e-12 rtol = 1.0e-12

        weights = collect(range(0.25, 1.25; length = size(overlap, 2)))
        weighted_block_occupancy = gto_occupancy_matrix(
            basis,
            supplement;
            weights = weights,
            block_indices = block_indices,
        )
        expected_weighted_block =
            (block_overlap .* reshape(weights, 1, :)) * transpose(block_overlap)
        @test weighted_block_occupancy ≈ expected_weighted_block atol = 1.0e-12 rtol = 1.0e-12
        @test gto_occupancy_matrix(basis, supplement, block_indices; weights = weights) ≈
              expected_weighted_block atol = 1.0e-12 rtol = 1.0e-12

        shell_equalized_error = try
            gto_occupancy_matrix(basis, supplement; weights = :shell_equalized)
            nothing
        catch err
            err
        end
        @test shell_equalized_error isa ArgumentError
        @test occursin("not implemented yet", sprint(showerror, shell_equalized_error))

        expansion = _truncate_coulomb_expansion(coulomb_gaussian_expansion(doacc = false), 3)
        s_supplement = legacy_atomic_gaussian_supplement(
            "He",
            "repo-sp";
            lmax = 0,
            basisfile = path,
        )
        ordinary_ops = ordinary_cartesian_qiu_white_operators(
            basis,
            s_supplement;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        hybrid_representation = basis_representation(ordinary_ops)
        hybrid_raw = GaussletBases._cartesian_raw_components(hybrid_representation)
        hybrid_overlap = gto_overlap_matrix(ordinary_ops, s_supplement)
        hybrid_reference = transpose(hybrid_raw.raw_to_final) * [
            hybrid_raw.exact_cartesian_supplement_overlap
            hybrid_raw.exact_supplement_overlap
        ]
        @test size(hybrid_overlap, 1) == ordinary_ops.gausslet_count + ordinary_ops.residual_count
        @test norm(hybrid_overlap - hybrid_reference, Inf) < 1.0e-12

        nuclei = [(0.0, 0.0, -0.7), (0.0, 0.0, 0.7)]
        diatomic_basis = bond_aligned_homonuclear_qw_basis(
            bond_length = 1.4,
            core_spacing = 0.7,
            xmax_parallel = 2.5,
            xmax_transverse = 2.0,
            bond_axis = :z,
        )
        diatomic_supplement = legacy_bond_aligned_diatomic_gaussian_supplement(
            "H",
            "repo-sp",
            nuclei;
            lmax = 1,
            basisfile = path,
        )
        diatomic_probe = basis_representation(diatomic_supplement)
        diatomic_overlap = gto_overlap_matrix(diatomic_basis, diatomic_supplement)
        diatomic_reference = GaussletBases._cartesian_basis_supplement_cross(
            basis_representation(diatomic_basis),
            diatomic_probe,
        )
        @test diatomic_probe isa CartesianGaussianShellSupplementRepresentation3D
        @test norm(diatomic_overlap - diatomic_reference, Inf) < 1.0e-12

        heteronuclear_supplement = legacy_bond_aligned_heteronuclear_gaussian_supplement(
            "He",
            "repo-sp",
            "H",
            "repo-sp",
            nuclei;
            lmax = 0,
            basisfile = path,
        )
        @test basis_representation(heteronuclear_supplement) isa
              CartesianGaussianShellSupplementRepresentation3D
    end
end

@testset "Vendored legacy BasisSets lookup and overrides" begin
    vendored = GaussletBases._vendored_legacy_basisfile_path()
    @test isfile(vendored)
    @test occursin(joinpath("data", "legacy", "BasisSets"), vendored)

    withenv("GAUSSLETBASES_BASISSETS_PATH" => nothing) do
        @test GaussletBases._legacy_basisfile_path() == vendored
        he_supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 0)
        h_supplement = legacy_atomic_gaussian_supplement("H", "cc-pVTZ"; lmax = 1)
        @test he_supplement.basisfile == vendored
        @test h_supplement.basisfile == vendored
    end

    mktemp() do path, io
        write(
            io,
            "#BASIS SET: He repo-test\n" *
            "He    S\n" *
            "      1.0000000              1.0000000\n" *
            "END\n",
        )
        close(io)

        withenv("GAUSSLETBASES_BASISSETS_PATH" => path) do
            @test GaussletBases._legacy_basisfile_path() == path
            supplement = legacy_atomic_gaussian_supplement("He", "repo-test"; lmax = 0)
            @test supplement.basisfile == path
            @test supplement.primitive_exponents == [1.0]
        end

        supplement = legacy_atomic_gaussian_supplement("He", "repo-test"; lmax = 0, basisfile = path)
        @test supplement.basisfile == path
        @test supplement.primitive_exponents == [1.0]
    end
end

@testset "Atomic lmax=0 supplement uses the explicit 3D shell route in QW consumers" begin
    if !_legacy_basisfile_available()
        @test true
    else
        source_basis, _legacy_old, baseline_ops, baseline_check = _qiu_white_full_nearest_fixture()
        supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 0)
        @test !GaussletBases._legacy_atomic_has_nonseparable_shells(supplement)
        ordinary_ops = ordinary_cartesian_qiu_white_operators(
            source_basis,
            supplement;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        explicit_ordinary_context = GaussletBases._normalized_atomic_build_context(
            source_basis,
            supplement,
        )
        explicit_ordinary = GaussletBases._ordinary_cartesian_qiu_white_operators_atomic(
            explicit_ordinary_context;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :numerical_reference,
            timing = false,
        )
        ordinary_check = GaussletBases.ordinary_cartesian_1s2_check(ordinary_ops)
        explicit_ordinary_check = GaussletBases.ordinary_cartesian_1s2_check(explicit_ordinary)

        (
            _source_basis_nested,
            _bundle_nested,
            _shell_nested,
            _fixed_block_nested,
            _shell_plus_core_nested,
            fixed_block_shell_plus_core,
            _legacy_nested,
            _baseline_nested,
            _nested_shell_only,
            nested_shell_plus_core,
            _baseline_nested_check,
            _nested_shell_only_check,
            nested_shell_plus_core_check,
        ) = _nested_qiu_white_nearest_fixture()
        nested_ops = ordinary_cartesian_qiu_white_operators(
            fixed_block_shell_plus_core,
            supplement;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        explicit_nested_context = GaussletBases._normalized_atomic_build_context(
            fixed_block_shell_plus_core,
            supplement,
        )
        explicit_nested = GaussletBases._ordinary_cartesian_qiu_white_operators_atomic(
            explicit_nested_context;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :numerical_reference,
            timing = false,
        )
        nested_check = GaussletBases.ordinary_cartesian_1s2_check(nested_ops)
        explicit_nested_check = GaussletBases.ordinary_cartesian_1s2_check(explicit_nested)

        @test ordinary_ops.gaussian_data isa LegacyAtomicGaussianSupplement
        @test nested_ops.gaussian_data isa LegacyAtomicGaussianSupplement
        @test ordinary_ops.residual_count > 0
        @test nested_ops.residual_count > 0
        @test size(ordinary_ops.raw_to_final, 2) ==
              ordinary_ops.gausslet_count + ordinary_ops.residual_count
        @test size(nested_ops.raw_to_final, 2) ==
              nested_ops.gausslet_count + nested_ops.residual_count
        @test norm(ordinary_ops.overlap - I, Inf) < 1.0e-8
        @test norm(nested_ops.overlap - I, Inf) < 1.0e-8
        @test ordinary_check.overlap_error < 1.0e-8
        @test nested_check.overlap_error < 1.0e-8
        @test ordinary_check.orbital_energy ≈ explicit_ordinary_check.orbital_energy atol = 1.0e-12 rtol = 1.0e-12
        @test ordinary_check.vee_expectation ≈ explicit_ordinary_check.vee_expectation atol = 1.0e-12 rtol = 1.0e-12
        @test nested_check.orbital_energy ≈ explicit_nested_check.orbital_energy atol = 1.0e-12 rtol = 1.0e-12
        @test nested_check.vee_expectation ≈ explicit_nested_check.vee_expectation atol = 1.0e-12 rtol = 1.0e-12
        @test norm(ordinary_ops.overlap - explicit_ordinary.overlap, Inf) < 1.0e-12
        @test norm(ordinary_ops.one_body_hamiltonian - explicit_ordinary.one_body_hamiltonian, Inf) < 1.0e-12
        @test norm(nested_ops.overlap - explicit_nested.overlap, Inf) < 1.0e-12
        @test norm(nested_ops.one_body_hamiltonian - explicit_nested.one_body_hamiltonian, Inf) < 1.0e-12
        @test ordinary_check.orbital_energy ≈ baseline_check.orbital_energy atol = 1.0e-12 rtol = 1.0e-12
        @test nested_check.orbital_energy ≈ nested_shell_plus_core_check.orbital_energy atol = 1.0e-12 rtol = 1.0e-12
    end
end

@testset "Active atomic lmax=1 supplement is explicit and physical in QW routes" begin
    if !_legacy_basisfile_available()
        @test true
    else
        supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 1)
        supplement3d = GaussletBases._atomic_cartesian_shell_supplement_3d(supplement)

        @test GaussletBases._legacy_atomic_has_nonseparable_shells(supplement)
        @test any(orbital -> orbital.label == "px1", supplement3d.orbitals)
        @test any(orbital -> orbital.label == "py1", supplement3d.orbitals)
        @test any(orbital -> orbital.label == "pz1", supplement3d.orbitals)

        source_basis_hybrid, _legacy_old, _pure_check, _toy_check, _legacy_check =
            _legacy_he_s_hybrid_fixture("cc-pVTZ")
        hybrid_err = try
            hybrid_mapped_ordinary_basis(
                source_basis_hybrid;
                core_gaussians = supplement,
                backend = :pgdg_localized_experimental,
            )
            nothing
        catch err
            err
        end
        @test hybrid_err isa ArgumentError
        @test occursin("l > 0", sprint(showerror, hybrid_err))
        @test occursin("explicit 3D", sprint(showerror, hybrid_err))

        source_basis_qw, _legacy_qw, ordinary_l0, ordinary_l0_check = _qiu_white_full_nearest_fixture()
        ordinary_l1 = ordinary_cartesian_qiu_white_operators(
            source_basis_qw,
            supplement;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        ordinary_l1_check = GaussletBases.ordinary_cartesian_1s2_check(ordinary_l1)
        @test ordinary_l1.gaussian_data isa LegacyAtomicGaussianSupplement
        @test ordinary_l1.residual_count > 0
        @test ordinary_l1_check.overlap_error < 1.0e-8
        @test isfinite(ordinary_l1_check.orbital_energy)
        @test isfinite(ordinary_l1_check.vee_expectation)
        @test ordinary_l1_check.vee_expectation > 0.0
        @test abs(ordinary_l1_check.orbital_energy - ordinary_l0_check.orbital_energy) > 1.0e-6
        @test abs(ordinary_l1_check.vee_expectation - ordinary_l0_check.vee_expectation) > 1.0e-4
        @test size(ordinary_l1.one_body_hamiltonian) != size(ordinary_l0.one_body_hamiltonian) ||
              norm(ordinary_l1.one_body_hamiltonian - ordinary_l0.one_body_hamiltonian, Inf) > 1.0e-6

        (
            _source_basis_nested,
            _bundle_nested,
            _shell_nested,
            _fixed_block_nested,
            _shell_plus_core_nested,
            fixed_block_shell_plus_core,
            _legacy_nested,
            _baseline_nested,
            _nested_shell_only,
            nested_l0,
            _baseline_nested_check,
            _nested_shell_only_check,
            nested_l0_check,
        ) = _nested_qiu_white_nearest_fixture()
        nested_l1 = ordinary_cartesian_qiu_white_operators(
            fixed_block_shell_plus_core,
            supplement;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        nested_l1_check = GaussletBases.ordinary_cartesian_1s2_check(nested_l1)
        @test nested_l1.gaussian_data isa LegacyAtomicGaussianSupplement
        @test nested_l1.residual_count > 0
        @test nested_l1_check.overlap_error < 1.0e-8
        @test isfinite(nested_l1_check.orbital_energy)
        @test isfinite(nested_l1_check.vee_expectation)
        @test nested_l1_check.vee_expectation > 0.0
        @test abs(nested_l1_check.orbital_energy - nested_l0_check.orbital_energy) > 1.0e-6
        @test abs(nested_l1_check.vee_expectation - nested_l0_check.vee_expectation) > 1.0e-4
        @test size(nested_l1.one_body_hamiltonian) != size(nested_l0.one_body_hamiltonian) ||
              norm(nested_l1.one_body_hamiltonian - nested_l0.one_body_hamiltonian, Inf) > 1.0e-6
    end
end

@testset "Ne cc-pV6Z atomic shell referee for lmax=0 and lmax=1" begin
    mktemp() do path, io
        write(io, _ne_repo_v6z_sp_basis_text())
        close(io)

        basis = build_basis(MappedUniformBasisSpec(:G10;
            count = 9,
            mapping = fit_asinh_mapping_for_strength(s = 0.8, npoints = 9, xmax = 6.0),
            reference_spacing = 1.0,
        ))
        expansion = coulomb_gaussian_expansion(doacc = false)
        bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
            basis,
            exponents = expansion.exponents,
            center = 0.0,
            backend = :numerical_reference,
        )

        for lmax in (0, 1)
            supplement = legacy_atomic_gaussian_supplement("Ne", "repo-v6z-sp"; lmax = lmax, basisfile = path)
            supplement3d = GaussletBases._atomic_cartesian_shell_supplement_3d(supplement)
            public_ops = ordinary_cartesian_qiu_white_operators(
                basis,
                supplement;
                expansion = expansion,
                Z = 10.0,
                interaction_treatment = :ggt_nearest,
            )
            explicit_context = GaussletBases._normalized_atomic_build_context(
                basis,
                supplement,
            )
            explicit_ops = GaussletBases._ordinary_cartesian_qiu_white_operators_atomic(
                explicit_context;
                expansion = expansion,
                Z = 10.0,
                interaction_treatment = :ggt_nearest,
                gausslet_backend = :numerical_reference,
                timing = false,
            )
            supplement3d = GaussletBases._atomic_cartesian_shell_supplement_3d(supplement)
            blocks = GaussletBases._qwrg_atomic_cartesian_blocks_3d(bundle, supplement3d, expansion)
            shell_h1 = GaussletBases._qwrg_atomic_cartesian_one_body_aa(blocks, expansion; Z = 10.0)
            shell_vals = eigen(Hermitian(shell_h1), Hermitian(blocks.overlap_aa)).values

            @test norm(public_ops.overlap - explicit_ops.overlap, Inf) < 1.0e-12
            @test norm(public_ops.one_body_hamiltonian - explicit_ops.one_body_hamiltonian, Inf) < 1.0e-12
            @test abs(shell_vals[1] + 50.0) < 2.0e-3
            @test abs(shell_vals[2] + 12.5) < 2.0e-3

            if lmax == 0
                @test abs(shell_vals[3] + (50.0 / 9.0)) < 2.0e-3
            else
                @test length(supplement3d.orbitals) == 25
                @test maximum(abs.(shell_vals[3:5] .+ 12.5)) < 2.0e-3
                @test maximum(abs.(shell_vals[3:5] .- shell_vals[3])) < 1.0e-10
                @test abs(shell_vals[6] + (50.0 / 9.0)) < 2.0e-3
            end
        end
    end
end

@testset "Atomic lmax=2 supplement is explicit in QW routes but not yet molecular" begin
    mktemp() do path, io
        write(
            io,
            "#BASIS SET: He repo-spd\n" *
            "He    S\n" *
            "      1.0000000              1.0000000\n" *
            "He    P\n" *
            "      0.8000000              1.0000000\n" *
            "He    D\n" *
            "      0.6000000              1.0000000\n" *
            "END\n",
        )
        close(io)

        supplement = legacy_atomic_gaussian_supplement("He", "repo-spd"; lmax = 2, basisfile = path)
        supplement3d = GaussletBases._atomic_cartesian_shell_supplement_3d(supplement)

        @test any(orbital -> orbital.label == "dxx1", supplement3d.orbitals)
        @test any(orbital -> orbital.label == "dyy1", supplement3d.orbitals)
        @test any(orbital -> orbital.label == "dzz1", supplement3d.orbitals)
        @test any(orbital -> orbital.label == "dxy1", supplement3d.orbitals)
        @test any(orbital -> orbital.label == "dxz1", supplement3d.orbitals)
        @test any(orbital -> orbital.label == "dyz1", supplement3d.orbitals)

        source_basis = build_basis(MappedUniformBasisSpec(:G10;
            count = 9,
            mapping = fit_asinh_mapping_for_strength(s = 0.8, npoints = 9, xmax = 6.0),
            reference_spacing = 1.0,
        ))
        ordinary_ops = ordinary_cartesian_qiu_white_operators(
            source_basis,
            supplement;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        ordinary_check = GaussletBases.ordinary_cartesian_1s2_check(ordinary_ops)
        @test ordinary_ops.gaussian_data isa LegacyAtomicGaussianSupplement
        @test ordinary_ops.residual_count > 0
        @test ordinary_check.overlap_error < 1.0e-8
        @test isfinite(ordinary_check.orbital_energy)
        @test isfinite(ordinary_check.vee_expectation)
        @test ordinary_check.vee_expectation > 0.0

        diatomic = legacy_bond_aligned_diatomic_gaussian_supplement(
            "He",
            "repo-spd",
            [(-0.7, 0.0, 0.0), (0.7, 0.0, 0.0)];
            lmax = 2,
            basisfile = path,
        )
        diatomic_err = try
            GaussletBases._bond_aligned_diatomic_cartesian_shell_supplement_3d(diatomic)
            nothing
        catch err
            err
        end
        @test diatomic_err isa ArgumentError
        @test occursin("lmax <= 1", sprint(showerror, diatomic_err))
    end
end

@testset "Legacy He s hybrid supplement check" begin
    if !_legacy_basisfile_available()
        @test true
    else
        (
            source_basis,
            legacy,
            pure_check,
            toy_check,
            legacy_check,
        ) = _legacy_he_s_hybrid_fixture("cc-pVTZ")

        @test source_basis isa MappedUniformBasis
        @test legacy isa LegacySGaussianData
        @test isfinite(pure_check.orbital_energy)
        @test isfinite(toy_check.orbital_energy)
        @test isfinite(legacy_check.orbital_energy)
        @test isfinite(pure_check.vee_expectation)
        @test isfinite(toy_check.vee_expectation)
        @test isfinite(legacy_check.vee_expectation)
        @test legacy_check.vee_expectation > 0.0
        @test abs(legacy_check.orbital_energy + 2.0) < abs(toy_check.orbital_energy + 2.0)
        @test abs(legacy_check.orbital_energy + 2.0) < abs(pure_check.orbital_energy + 2.0)
    end
end

@testset "Legacy He s MWG residual interaction" begin
    if !_legacy_basisfile_available()
        @test true
    else
        (
            source_basis,
            legacy,
            hybrid_basis,
            pure_operators,
            combined_operators,
            nearest_operators,
            mwg_operators,
            pure_check,
            combined_check,
            nearest_check,
            mwg_check,
            mwg_data,
        ) = _legacy_he_s_mwg_fixture("cc-pVTZ")

        @test source_basis isa MappedUniformBasis
        @test legacy isa LegacySGaussianData
        @test hybrid_basis isa HybridMappedOrdinaryBasis1D
        @test combined_operators.interaction_treatment == :combined_basis
        @test nearest_operators.interaction_treatment == :residual_gaussian_nearest
        @test mwg_operators.interaction_treatment == :residual_gaussian_mwg
        @test norm(combined_operators.overlap_3d - nearest_operators.overlap_3d, Inf) < 1.0e-12
        @test norm(combined_operators.overlap_3d - mwg_operators.overlap_3d, Inf) < 1.0e-12
        @test norm(combined_operators.one_body_hamiltonian - nearest_operators.one_body_hamiltonian, Inf) < 1.0e-12
        @test norm(combined_operators.one_body_hamiltonian - mwg_operators.one_body_hamiltonian, Inf) < 1.0e-12
        @test mwg_operators.interaction_matrix ≈ transpose(mwg_operators.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
        @test all(isfinite, mwg_data.residual_centers)
        @test all(isfinite, mwg_data.residual_widths)
        @test all(>(0.0), mwg_data.residual_widths)
        @test maximum(abs.(mwg_data.residual_centers)) < 1.0e-8
        @test pure_check.vee_expectation > 0.0
        @test combined_check.vee_expectation > 0.0
        @test nearest_check.vee_expectation > 0.0
        @test mwg_check.vee_expectation > 0.0
        @test abs(mwg_check.vee_expectation - nearest_check.vee_expectation) > 1.0e-3
        @test abs(mwg_check.orbital_energy - combined_check.orbital_energy) < 1.0e-12
    end
end

@testset "Qiu-White residual Gaussian reference path" begin
    if !_legacy_basisfile_available() || !_RUN_SLOW_TESTS
        @test true
    else
        (
            source_basis,
            legacy,
            surrogate_mwg,
            qiu_nearest,
            qiu_mwg,
            surrogate_check,
            nearest_check,
            mwg_check,
        ) = _qiu_white_reference_fixture()

        @test source_basis isa MappedUniformBasis
        @test legacy isa LegacySGaussianData
        @test surrogate_mwg isa OrdinaryCartesianIDAOperators
        @test qiu_nearest isa QiuWhiteResidualGaussianOperators
        @test qiu_mwg isa QiuWhiteResidualGaussianOperators
        @test qiu_nearest.interaction_treatment == :ggt_nearest
        @test qiu_mwg.interaction_treatment == :mwg
        @test qiu_mwg.gausslet_backend == :numerical_reference
        @test qiu_nearest.residual_count > 0
        @test qiu_mwg.gausslet_count == length(source_basis)^3
        @test qiu_mwg.residual_count > 0
        @test size(qiu_mwg.raw_to_final, 1) == qiu_mwg.gausslet_count + length(legacy.gaussians)
        @test size(qiu_mwg.raw_to_final, 2) == qiu_mwg.gausslet_count + qiu_mwg.residual_count
        @test norm(qiu_nearest.overlap - I, Inf) < 1.0e-8
        @test norm(qiu_mwg.overlap - I, Inf) < 1.0e-8
        @test qiu_nearest.one_body_hamiltonian ≈ transpose(qiu_nearest.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
        @test qiu_mwg.one_body_hamiltonian ≈ transpose(qiu_mwg.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
        @test qiu_nearest.interaction_matrix ≈ transpose(qiu_nearest.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
        @test qiu_mwg.interaction_matrix ≈ transpose(qiu_mwg.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
        @test size(qiu_mwg.interaction_matrix) == (qiu_mwg.gausslet_count + qiu_mwg.residual_count, qiu_mwg.gausslet_count + qiu_mwg.residual_count)
        @test all(isfinite, qiu_nearest.residual_centers)
        @test all(isnan, qiu_nearest.residual_widths)
        @test all(isfinite, qiu_mwg.residual_centers)
        @test all(isfinite, qiu_mwg.residual_widths)
        @test all(>(0.0), vec(qiu_mwg.residual_widths))
        @test maximum(abs.(qiu_mwg.residual_centers)) < 1.0e-6
        @test isfinite(surrogate_check.orbital_energy)
        @test isfinite(nearest_check.orbital_energy)
        @test isfinite(mwg_check.orbital_energy)
        @test isfinite(surrogate_check.vee_expectation)
        @test isfinite(nearest_check.vee_expectation)
        @test isfinite(mwg_check.vee_expectation)
        @test surrogate_check.vee_expectation > 0.0
        @test nearest_check.vee_expectation > 0.0
        @test mwg_check.vee_expectation > 0.0

        (
            _full_source_basis,
            _full_legacy,
            qiu_full_nearest,
            qiu_full_nearest_check,
        ) = _qiu_white_full_nearest_fixture()

        @test norm(qiu_full_nearest.overlap - I, Inf) < 1.0e-8
        @test qiu_full_nearest.residual_count > 0
        @test -3.5 < qiu_full_nearest_check.orbital_energy < -1.8
        @test 1.0 < qiu_full_nearest_check.vee_expectation < 1.4
        @test abs(qiu_full_nearest_check.vee_expectation - 1.25) < 0.05
    end
end

@testset "Ordinary Cartesian localized backend" begin
    expansion = _truncate_coulomb_expansion(coulomb_gaussian_expansion(doacc = false), 3)
    basis = build_basis(MappedUniformBasisSpec(:G10;
        count = 5,
        mapping = fit_asinh_mapping_for_strength(s = 0.5, npoints = 5, xmax = 6.0),
        reference_spacing = 1.0,
    ))

    proxy = ordinary_cartesian_ida_operators(
        basis;
        expansion = expansion,
        Z = 2.0,
        backend = :pgdg_experimental,
    )
    localized = ordinary_cartesian_ida_operators(
        basis;
        expansion = expansion,
        Z = 2.0,
        backend = :pgdg_localized_experimental,
    )
    (_, _, overlap_reference_localized, kinetic_reference_localized, gaussian_reference_localized) =
        _localized_numerical_reference_1d(basis, expansion.exponents)
    raw_localized = mapped_pgdg_localized(GaussletBases.mapped_pgdg_logfit_prototype(basis))
    oracle_localized = GaussletBases._mapped_ordinary_localized_oracle_operators(
        basis;
        exponents = expansion.exponents,
        center = 0.0,
    )
    raw_h1, _ = _cartesian_hydrogen_energy(
        overlap_matrix(raw_localized),
        kinetic_matrix(raw_localized),
        gaussian_factor_matrices(raw_localized; exponents = expansion.exponents, center = 0.0),
        expansion;
        Z = 2.0,
    )
    corrected_h1, _ = _cartesian_hydrogen_energy(
        localized.one_body_1d.overlap,
        localized.one_body_1d.kinetic,
        localized.one_body_1d.gaussian_factors,
        expansion;
        Z = 2.0,
    )
    oracle_h1, _ = _cartesian_hydrogen_energy(
        oracle_localized.overlap,
        oracle_localized.kinetic,
        oracle_localized.gaussian_factors,
        expansion;
        Z = 2.0,
    )
    reference_h1, _ = _cartesian_hydrogen_energy(
        overlap_reference_localized,
        kinetic_reference_localized,
        gaussian_reference_localized,
        expansion;
        Z = 2.0,
    )

    @test localized isa OrdinaryCartesianIDAOperators
    @test localized.backend == :pgdg_localized_experimental
    @test occursin("experimental=true", sprint(show, localized))
    @test norm(localized.one_body_1d.overlap - I, Inf) < 1.0e-10
    @test norm(localized.overlap_3d - I, Inf) < 1.0e-9
    @test norm(localized.one_body_1d.overlap - I, Inf) < norm(proxy.one_body_1d.overlap - I, Inf)
    @test norm(localized.overlap_3d - I, Inf) < norm(proxy.overlap_3d - I, Inf)
    @test norm(localized.one_body_1d.overlap - overlap_reference_localized, Inf) < 1.0e-10
    @test norm(localized.one_body_1d.kinetic - kinetic_reference_localized, Inf) <
          norm(kinetic_matrix(raw_localized) - kinetic_reference_localized, Inf)
    @test norm(corrected_h1 - reference_h1, Inf) < norm(raw_h1 - reference_h1, Inf)
    @test norm(oracle_h1 - reference_h1, Inf) < norm(corrected_h1 - reference_h1, Inf)
    @test localized.one_body_hamiltonian ≈ transpose(localized.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test localized.interaction_matrix ≈ transpose(localized.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test minimum(diag(localized.interaction_matrix)) > 0.0
end

@testset "Legacy 1D hybrid mapped ordinary basis" begin
    (
        basis,
        core_gaussians,
        expansion,
        hybrid_reference,
        hybrid_analytic,
        one_body_reference,
        one_body_analytic,
        hard_reference,
        hard_analytic,
        energy_reference,
        energy_analytic,
        hard_energy_reference,
        hard_energy_analytic,
    ) = _quick_hybrid_mapped_ordinary_fixture()

    @test hybrid_reference isa HybridMappedOrdinaryBasis1D
    @test hybrid_analytic isa HybridMappedOrdinaryBasis1D
    @test hybrid_reference.backend == :numerical_reference
    @test hybrid_analytic.backend == :pgdg_localized_experimental
    @test length(core_gaussians) == 2
    @test length(hybrid_reference) == length(basis) + length(core_gaussians)
    @test length(hybrid_analytic) == length(basis) + length(core_gaussians)
    @test occursin("experimental=true", sprint(show, hybrid_analytic))
    @test !occursin("experimental=true", sprint(show, hybrid_reference))
    @test norm(one_body_reference.overlap - I, Inf) < 1.0e-10
    @test norm(one_body_analytic.overlap - I, Inf) < 1.0e-10
    @test one_body_analytic.kinetic ≈ transpose(one_body_analytic.kinetic) atol = 1.0e-10 rtol = 1.0e-10
    hybrid_factor_diff = maximum(
        norm(one_body_analytic.gaussian_factors[index] - one_body_reference.gaussian_factors[index], Inf)
        for index in eachindex(one_body_analytic.gaussian_factors)
    )
    hard_factor_diff = maximum(
        norm(hard_analytic.gaussian_factors[index] - hard_reference.gaussian_factors[index], Inf)
        for index in eachindex(hard_analytic.gaussian_factors)
    )
    @test hybrid_factor_diff < hard_factor_diff
    @test abs(energy_analytic - energy_reference) <
          abs(hard_energy_analytic - hard_energy_reference)
end

@testset "Legacy 1D hybrid ordinary mapped SHO smoke" begin
    (
        hybrid_analytic,
        centered_analytic,
        centered_hamiltonian,
    ) = _quick_ordinary_sho_smoke_fixture()

    @test occursin("experimental=true", sprint(show, hybrid_analytic))
    @test centered_hamiltonian.backend == :pgdg_localized_experimental
    @test centered_hamiltonian.overlap ≈ transpose(centered_hamiltonian.overlap) atol = 1.0e-10 rtol = 1.0e-10
    @test centered_hamiltonian.hamiltonian ≈ transpose(centered_hamiltonian.hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test norm(centered_hamiltonian.overlap - I, Inf) < 1.0e-10
    @test issorted(centered_analytic.eigenvalues)
    @test centered_analytic.eigenvalues[1] > 0.0
    @test abs(centered_analytic.eigenvalues[1] - centered_analytic.exact[1]) < 1.0e-2
    @test centered_analytic.kinetic_expectation > 0.0
    @test centered_analytic.displacement2_expectation > 0.0
end

if _RUN_SLOW_TESTS
    @testset "Ordinary mapped SHO spectra" begin
        (
            centered_reference,
            centered_analytic,
            shifted_reference,
            shifted_analytic,
            stress_reference,
            stress_analytic,
        ) = _slow_ordinary_sho_fixture()

        centered_diff = maximum(abs.(centered_reference.eigenvalues .- centered_analytic.eigenvalues))
        shifted_diff = maximum(abs.(shifted_reference.eigenvalues .- shifted_analytic.eigenvalues))
        stress_diff = maximum(abs.(stress_reference.eigenvalues .- stress_analytic.eigenvalues))

        @test centered_diff < 2.0e-4
        @test shifted_diff < 2.0e-4
        @test stress_diff > 1.0e-2
        @test maximum(abs.(centered_analytic.eigenvalues .- centered_analytic.exact)) <
              maximum(abs.(stress_analytic.eigenvalues .- stress_analytic.exact))
        @test maximum(
            abs.(
                abs.(centered_analytic.eigenvalues .- centered_analytic.exact) .-
                abs.(centered_reference.eigenvalues .- centered_reference.exact),
            ),
        ) < 5.0e-4
        @test maximum(
            abs.(
                abs.(shifted_analytic.eigenvalues .- shifted_analytic.exact) .-
                abs.(shifted_reference.eigenvalues .- shifted_reference.exact),
            ),
        ) < 5.0e-4
    end

    @testset "Mapped Cartesian hydrogen" begin
        basis, mapping_value, representation, overlap_1d, kinetic_1d, expansion, hamiltonian, energy =
            _quick_mapped_cartesian_hydrogen_fixture()
        (_, _, _, _, _, _, _, _, _, unmapped_energy) = _quick_cartesian_hydrogen_fixture()

        @test mapping_value isa AsinhMapping
        @test overlap_1d ≈ transpose(overlap_1d) atol = 1.0e-10 rtol = 1.0e-10
        @test kinetic_1d ≈ transpose(kinetic_1d) atol = 1.0e-10 rtol = 1.0e-10
        @test length(expansion) == 45
        @test size(hamiltonian) == (length(basis)^3, length(basis)^3)
        @test energy < unmapped_energy - 0.02
        @test energy < -0.46
        @test energy > -0.5
    end

    @testset "Mapped Coulomb expansion calibration" begin
        expansion_115 = coulomb_gaussian_expansion(doacc = true, maxu = 115.0)
        expansion_135 = coulomb_gaussian_expansion(doacc = true, maxu = 135.0)

        (_, _, _, _, _, energy_115, _, _, _, _) = _mapped_cartesian_hydrogen_comparison(expansion_115)
        (_, _, _, _, _, energy_135, _, _, _, _) = _mapped_cartesian_hydrogen_comparison(expansion_135)

        @test length(expansion_115) == 115
        @test length(expansion_135) == 135
        @test abs(energy_135 - energy_115) < 1.0e-8
    end

    @testset "Mapped PGDG hydrogen prototype" begin
        expansion = coulomb_gaussian_expansion(doacc = false)
        (
            basis,
            prototype,
            localized,
            refined,
            refined_localized,
            overlap_numeric,
            kinetic_numeric,
            hamiltonian_numeric,
            energy_numeric,
            overlap_pgdg,
            position_pgdg,
            kinetic_pgdg,
            hamiltonian_pgdg,
            energy_pgdg,
            overlap_localized,
            position_localized,
            kinetic_localized,
            hamiltonian_localized,
            energy_localized,
            overlap_refined,
            position_refined,
            kinetic_refined,
            hamiltonian_refined,
            energy_refined,
            overlap_refined_localized,
            position_refined_localized,
            kinetic_refined_localized,
            hamiltonian_refined_localized,
            energy_refined_localized,
        ) = _mapped_cartesian_hydrogen_comparison(expansion)

        representation_numeric = basis_representation(basis; operators = (:overlap, :position, :kinetic))
        transform_numeric, centers_numeric = _cleanup_comx_transform(
            representation_numeric.basis_matrices.overlap,
            representation_numeric.basis_matrices.position,
            integral_weights(basis),
        )
        overlap_numeric_localized = transform_numeric' * representation_numeric.basis_matrices.overlap * transform_numeric
        kinetic_numeric_localized = transform_numeric' * representation_numeric.basis_matrices.kinetic * transform_numeric
        position_numeric_localized = transform_numeric' * representation_numeric.basis_matrices.position * transform_numeric
        factor_numeric = gaussian_factor_matrix(basis; exponent = 0.35, center = 0.0)
        factor_pgdg = gaussian_factor_matrix(prototype; exponent = 0.35, center = 0.0)
        factor_refined = gaussian_factor_matrix(refined; exponent = 0.35, center = 0.0)
        factor_numeric_localized = transform_numeric' * factor_numeric * transform_numeric
        factor_localized = gaussian_factor_matrix(localized; exponent = 0.35, center = 0.0)
        factor_refined_localized = gaussian_factor_matrix(refined_localized; exponent = 0.35, center = 0.0)
        span_pre = _subspace_overlap_metric(basis, prototype)
        span_refined = _subspace_overlap_metric(basis, refined)
        projector_pre = _projector_difference_metric(basis, prototype)
        projector_refined = _projector_difference_metric(basis, refined)

        @test size(hamiltonian_numeric) == size(hamiltonian_pgdg)
        @test size(hamiltonian_numeric) == size(hamiltonian_localized)
        @test size(hamiltonian_numeric) == size(hamiltonian_refined)
        @test size(hamiltonian_numeric) == size(hamiltonian_refined_localized)
        @test overlap_pgdg ≈ transpose(overlap_pgdg) atol = 1.0e-10 rtol = 1.0e-10
        @test kinetic_pgdg ≈ transpose(kinetic_pgdg) atol = 1.0e-10 rtol = 1.0e-10
        @test overlap_localized ≈ I atol = 1.0e-10 rtol = 1.0e-10
        @test position_localized ≈ Diagonal(centers(localized)) atol = 1.0e-10 rtol = 1.0e-10
        @test overlap_numeric_localized ≈ I atol = 1.0e-10 rtol = 1.0e-10
        @test position_numeric_localized ≈ Diagonal(centers_numeric) atol = 1.0e-10 rtol = 1.0e-10
        @test overlap_refined ≈ transpose(overlap_refined) atol = 1.0e-10 rtol = 1.0e-10
        @test kinetic_refined ≈ transpose(kinetic_refined) atol = 1.0e-10 rtol = 1.0e-10
        @test overlap_refined_localized ≈ I atol = 1.0e-10 rtol = 1.0e-10
        @test position_refined_localized ≈ Diagonal(centers(refined_localized)) atol = 1.0e-10 rtol = 1.0e-10
        @test norm(kinetic_numeric_localized - kinetic_localized, Inf) < norm(kinetic_numeric - kinetic_pgdg, Inf)
        @test norm(factor_numeric_localized - factor_localized, Inf) < norm(factor_numeric - factor_pgdg, Inf)
        @test energy_numeric < -0.47
        @test energy_pgdg < -0.47
        @test energy_localized < -0.47
        @test energy_refined < -0.47
        @test energy_refined_localized < -0.47
        @test abs(energy_numeric - energy_pgdg) < 0.01
        @test abs(energy_localized - energy_pgdg) < 1.0e-10
        @test abs(energy_numeric - energy_localized) ≤ abs(energy_numeric - energy_pgdg) + 1.0e-10
        @test span_refined.min_sv > span_pre.min_sv
        @test projector_refined.frob < projector_pre.frob
        @test projector_refined.op < projector_pre.op
        @test abs(energy_numeric - energy_refined) < abs(energy_numeric - energy_pgdg)
        @test abs(energy_numeric - energy_refined_localized) ≤ abs(energy_numeric - energy_refined) + 1.0e-10
        @test norm(factor_numeric - factor_refined, Inf) ≤ norm(factor_numeric - factor_pgdg, Inf) + 0.02
        @test norm(factor_numeric_localized - factor_refined_localized, Inf) ≤ norm(factor_numeric_localized - factor_localized, Inf) + 0.02
    end

    @testset "Mapped ordinary backend hydrogen regimes" begin
        expansion = coulomb_gaussian_expansion(doacc = false)

        function backend_energy(s_value)
            basis = build_basis(MappedUniformBasisSpec(:G10;
                count = 5,
                mapping = fit_asinh_mapping_for_strength(s = s_value, npoints = 5, xmax = 6.0),
                reference_spacing = 1.0,
            ))
            reference = mapped_ordinary_one_body_operators(
                basis;
                exponents = expansion.exponents,
                backend = :numerical_reference,
            )
            analytic = mapped_ordinary_one_body_operators(
                basis;
                exponents = expansion.exponents,
                backend = :pgdg_experimental,
            )
            return (
                basis,
                mapped_cartesian_hydrogen_energy(reference, expansion; Z = 1.0),
                mapped_cartesian_hydrogen_energy(analytic, expansion; Z = 1.0),
            )
        end

        mild_basis, mild_energy_reference, mild_energy_analytic = backend_energy(0.5)
        moderate_basis, moderate_energy_reference, moderate_energy_analytic = backend_energy(1.0)
        stress_basis, stress_energy_reference, stress_energy_analytic = backend_energy(2.0)

        mild_diff = abs(mild_energy_reference - mild_energy_analytic)
        moderate_diff = abs(moderate_energy_reference - moderate_energy_analytic)
        stress_diff = abs(stress_energy_reference - stress_energy_analytic)

        @test mild_basis isa MappedUniformBasis
        @test moderate_basis isa MappedUniformBasis
        @test stress_basis isa MappedUniformBasis
        @test mild_diff < 3.0e-4
        @test moderate_diff < 1.0e-3
        @test stress_diff > moderate_diff
        @test stress_diff < 0.02
        @test stress_energy_analytic < -0.45
        @test stress_energy_reference < -0.45
    end

    @testset "Ordinary Cartesian IDA backend regimes" begin
        expansion = _truncate_coulomb_expansion(coulomb_gaussian_expansion(doacc = false), 3)

        function ida_difference(s_value)
            basis = build_basis(MappedUniformBasisSpec(:G10;
                count = 5,
                mapping = fit_asinh_mapping_for_strength(s = s_value, npoints = 5, xmax = 6.0),
                reference_spacing = 1.0,
            ))
            reference = ordinary_cartesian_ida_operators(
                basis;
                expansion = expansion,
                Z = 2.0,
                backend = :numerical_reference,
            )
            analytic = ordinary_cartesian_ida_operators(
                basis;
                expansion = expansion,
                Z = 2.0,
                backend = :pgdg_experimental,
            )
            return (
                norm(reference.one_body_hamiltonian - analytic.one_body_hamiltonian, Inf),
                norm(reference.interaction_matrix - analytic.interaction_matrix, Inf),
            )
        end

        mild_h1_diff, mild_vee_diff = ida_difference(0.5)
        stress_h1_diff, stress_vee_diff = ida_difference(2.0)

        @test mild_h1_diff < 0.05
        @test mild_vee_diff < 0.1
        @test stress_h1_diff > mild_h1_diff
        @test stress_vee_diff > mild_vee_diff
        @test stress_vee_diff < 1.0
    end

    @testset "Ordinary Cartesian IDA reference agreement" begin
        mild_basis, mild_expansion, mild_reference = _quick_ordinary_cartesian_ida_fixture(
            backend = :numerical_reference,
            mapped = true,
            s = 0.5,
        )
        (_, _, mild_analytic) = _quick_ordinary_cartesian_ida_fixture(
            backend = :pgdg_experimental,
            mapped = true,
            s = 0.5,
        )
        (_, _, identity_reference) = _quick_ordinary_cartesian_ida_fixture(
            backend = :numerical_reference,
            mapped = false,
        )
        (_, _, identity_analytic) = _quick_ordinary_cartesian_ida_fixture(
            backend = :pgdg_experimental,
            mapped = false,
        )

        @test mild_basis isa MappedUniformBasis
        @test mild_reference.backend == :numerical_reference
        @test !occursin("experimental=true", sprint(show, mild_reference))
        @test norm(identity_reference.one_body_hamiltonian - identity_analytic.one_body_hamiltonian, Inf) < 1.0e-6
        @test norm(identity_reference.interaction_matrix - identity_analytic.interaction_matrix, Inf) < 1.0e-6
        @test norm(mild_reference.one_body_hamiltonian - mild_analytic.one_body_hamiltonian, Inf) < 0.01
        @test norm(mild_reference.interaction_matrix - mild_analytic.interaction_matrix, Inf) < 1.0e-3
    end

    @testset "Ordinary Cartesian localized backend agreement" begin
        expansion = _truncate_coulomb_expansion(coulomb_gaussian_expansion(doacc = false), 3)

        function build_pair(s_value)
            basis = build_basis(MappedUniformBasisSpec(:G10;
                count = 5,
                mapping = fit_asinh_mapping_for_strength(s = s_value, npoints = 5, xmax = 6.0),
                reference_spacing = 1.0,
            ))
            reference = ordinary_cartesian_ida_operators(
                basis;
                expansion = expansion,
                Z = 2.0,
                backend = :numerical_reference,
            )
            localized = ordinary_cartesian_ida_operators(
                basis;
                expansion = expansion,
                Z = 2.0,
                backend = :pgdg_localized_experimental,
            )
            return reference, localized
        end

        mild_reference, mild_localized = build_pair(0.5)
        stress_reference, stress_localized = build_pair(2.0)

        mild_h1_diff = norm(mild_reference.one_body_hamiltonian - mild_localized.one_body_hamiltonian, Inf)
        mild_vee_diff = norm(mild_reference.interaction_matrix - mild_localized.interaction_matrix, Inf)
        stress_h1_diff = norm(stress_reference.one_body_hamiltonian - stress_localized.one_body_hamiltonian, Inf)
        stress_vee_diff = norm(stress_reference.interaction_matrix - stress_localized.interaction_matrix, Inf)

        @test norm(mild_localized.one_body_1d.overlap - I, Inf) < 1.0e-10
        @test norm(mild_localized.overlap_3d - I, Inf) < 1.0e-9
        @test mild_h1_diff < 0.02
        @test mild_vee_diff < 1.0e-2
        @test stress_h1_diff > mild_h1_diff
        @test stress_vee_diff > mild_vee_diff
        @test stress_h1_diff < 0.2
        @test stress_vee_diff < 0.2
    end
end

@testset "Cartesian hydrogen via Coulomb expansion" begin
    basis, representation, expansion, overlap_1d, kinetic_1d, gaussian_factors, overlap_3d, nuclear_3d, hamiltonian, energy =
        _quick_cartesian_hydrogen_fixture()

    @test length(basis) == 5
    @test size(overlap_1d) == (length(basis), length(basis))
    @test norm(overlap_1d - I, Inf) ≤ 1.0e-10
    @test size(overlap_3d) == (length(basis)^3, length(basis)^3)
    @test norm(overlap_3d - I, Inf) ≤ 1.0e-8
    @test hamiltonian ≈ transpose(hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test size(nuclear_3d) == size(hamiltonian)
    @test minimum(diag(nuclear_3d)) < 0.0
    @test -0.7 < energy < -0.3
end
