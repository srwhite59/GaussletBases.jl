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
