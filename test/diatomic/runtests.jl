@testset "Bond-aligned diatomic QW reference path" begin
    basis14, operators14, check14 = _bond_aligned_diatomic_qw_fixture(; bond_length = 1.4)
    basis20, operators20, check20 = _bond_aligned_diatomic_qw_fixture(; bond_length = 2.0)

    @test basis14 isa BondAlignedDiatomicQWBasis3D
    @test operators14 isa QiuWhiteResidualGaussianOperators
    @test operators14.gaussian_data === nothing
    @test operators14.residual_count == 0
    @test operators14.gausslet_count == length(basis14.basis_x) * length(basis14.basis_y) * length(basis14.basis_z)
    @test norm(operators14.overlap - I, Inf) < 1.0e-8
    @test norm(operators20.overlap - I, Inf) < 1.0e-8
    @test operators14.one_body_hamiltonian ≈ transpose(operators14.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test operators14.interaction_matrix ≈ transpose(operators14.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test isfinite(check14.orbital_energy)
    @test isfinite(check14.vee_expectation)
    @test isfinite(check20.orbital_energy)
    @test isfinite(check20.vee_expectation)
    @test check14.orbital_energy < -1.0
    @test check20.orbital_energy < -1.0
    @test 0.5 < check14.vee_expectation < 1.0
    @test 0.5 < check20.vee_expectation < 1.0
    @test check14.orbital_energy < check20.orbital_energy
    @test length(basis14.basis_x) == length(basis14.basis_y)
    @test length(basis14.basis_z) > length(basis14.basis_x)
    @test length(basis20.basis_z) >= length(basis14.basis_z)
end

@testset "Per-center nuclear one-body reassembly on diatomic routes" begin
    basis, _operators, _check = _bond_aligned_diatomic_qw_fixture(; bond_length = 1.4)
    full_ops = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :by_center,
        interaction_treatment = :ggt_nearest,
    )
    atom_a_ops = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = [1.0, 0.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
    )
    atom_b_ops = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = [0.0, 1.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
    )

    @test full_ops.nuclear_term_storage == :by_center
    @test full_ops.nuclear_charges == [1.0, 1.0]
    @test !isnothing(full_ops.kinetic_one_body)
    @test !isnothing(full_ops.nuclear_one_body_by_center)
    @test length(full_ops.nuclear_one_body_by_center) == 2
    @test assembled_one_body_hamiltonian(full_ops) ≈
          full_ops.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test assembled_one_body_hamiltonian(full_ops; nuclear_charges = [1.0, 0.0]) ≈
          atom_a_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
    @test assembled_one_body_hamiltonian(full_ops; nuclear_charges = [0.0, 1.0]) ≈
          atom_b_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10

    full_ops_localized = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :by_center,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
    )
    atom_a_ops_localized = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = [1.0, 0.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
    )
    atom_b_ops_localized = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = [0.0, 1.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
    )

    @test full_ops_localized.gausslet_backend == :pgdg_localized_experimental
    @test full_ops_localized.nuclear_term_storage == :by_center
    @test !isnothing(full_ops_localized.kinetic_one_body)
    @test !isnothing(full_ops_localized.nuclear_one_body_by_center)
    @test length(full_ops_localized.nuclear_one_body_by_center) == 2
    @test assembled_one_body_hamiltonian(full_ops_localized) ≈
          full_ops_localized.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test assembled_one_body_hamiltonian(full_ops_localized; nuclear_charges = [1.0, 0.0]) ≈
          atom_a_ops_localized.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
    @test assembled_one_body_hamiltonian(full_ops_localized; nuclear_charges = [0.0, 1.0]) ≈
          atom_b_ops_localized.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10

    (
        _nested_basis,
        _nested_parent_ops,
        _nested_parent_check,
        expansion,
        _nested_source,
        fixed_block,
        _parent_modes,
        _parent_ground,
        _projected,
        _projected_vee,
        _capture,
        _projected_energy,
    ) = _bond_aligned_diatomic_nested_fixed_block_fixture(; bond_length = 1.4)

    full_nested_ops_localized = ordinary_cartesian_qiu_white_operators(
        fixed_block;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :by_center,
        expansion,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
    )
    # Keep this comparison on the same nested by-center kinetic contract.
    # The nested fixed-block kinetic payload is packet-level data and is not the
    # same contract as rebuilding a total-only one-body matrix through a
    # separate parent-space contraction path.
    atom_a_nested_ops_localized = ordinary_cartesian_qiu_white_operators(
        fixed_block;
        nuclear_charges = [1.0, 0.0],
        nuclear_term_storage = :by_center,
        expansion,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
    )
    atom_b_nested_ops_localized = ordinary_cartesian_qiu_white_operators(
        fixed_block;
        nuclear_charges = [0.0, 1.0],
        nuclear_term_storage = :by_center,
        expansion,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
    )

    @test full_nested_ops_localized.gausslet_backend == :pgdg_localized_experimental
    @test full_nested_ops_localized.nuclear_term_storage == :by_center
    @test !isnothing(full_nested_ops_localized.kinetic_one_body)
    @test !isnothing(full_nested_ops_localized.nuclear_one_body_by_center)
    @test length(full_nested_ops_localized.nuclear_one_body_by_center) == 2
    nested_factorized_basis = GaussletBases._nested_factorized_parent_basis(fixed_block)
    @test GaussletBases._nested_by_center_sidecar_path(fixed_block) == :factorized_final
    for backend in (:numerical_reference, :pgdg_localized_experimental)
        nested_bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(
            fixed_block.parent_basis,
            expansion;
            gausslet_backend = backend,
        )
        direct_contracted_nested_nuclear =
            GaussletBases._qwrg_bond_aligned_direct_contracted_nuclear_one_body_by_center(
                fixed_block.parent_basis,
                nested_factorized_basis,
                nested_bundles.bundle_x,
                nested_bundles.bundle_y,
                nested_bundles.bundle_z,
                expansion,
            )
        reference_nested_parent_nuclear = GaussletBases._qwrg_diatomic_nuclear_one_body_by_center(
            fixed_block.parent_basis,
            nested_bundles.bundle_x,
            nested_bundles.bundle_y,
            nested_bundles.bundle_z,
            expansion,
        )
        reference_nested_contracted_nuclear = [
            GaussletBases._qwrg_contract_parent_symmetric_matrix(
                fixed_block.coefficient_matrix,
                matrix,
            ) for matrix in reference_nested_parent_nuclear
        ]
        @test length(direct_contracted_nested_nuclear) == length(reference_nested_contracted_nuclear)
        for nucleus_index in eachindex(
            direct_contracted_nested_nuclear,
            reference_nested_contracted_nuclear,
        )
            @test direct_contracted_nested_nuclear[nucleus_index] ≈
                  reference_nested_contracted_nuclear[nucleus_index] atol = 1.0e-10 rtol = 1.0e-10
        end
    end

    nonfactorized_coefficients = Matrix{Float64}(fixed_block.coefficient_matrix)
    triplets = nested_factorized_basis.basis_triplets
    mixed_columns = nothing
    for left in eachindex(triplets), right in (left + 1):length(triplets)
        sum(triplets[left][axis] != triplets[right][axis] for axis in 1:3) >= 2 || continue
        mixed_columns = (left, right)
        break
    end
    @test !isnothing(mixed_columns)
    left_column, right_column = mixed_columns
    column_one = copy(nonfactorized_coefficients[:, left_column])
    column_two = copy(nonfactorized_coefficients[:, right_column])
    nonfactorized_coefficients[:, left_column] .= (column_one .+ column_two) ./ sqrt(2.0)
    nonfactorized_coefficients[:, right_column] .= (column_one .- column_two) ./ sqrt(2.0)
    nonfactorized_cache = GaussletBases._nested_eager_factorized_basis_cache(
        fixed_block.parent_basis,
        nonfactorized_coefficients,
    )
    @test isnothing(nonfactorized_cache[])
    @test_throws ArgumentError GaussletBases._nested_eager_factorized_basis_cache(
        (;),
        nonfactorized_coefficients,
    )

    nonfactorized_fixed_block = GaussletBases._NestedFixedBlock3D(
        fixed_block.parent_basis,
        fixed_block.shell,
        fixed_block.gausslet_backend,
        nonfactorized_coefficients,
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
        nonfactorized_cache,
        GaussletBases._nested_staged_by_center_sidecar_cache(),
    )
    @test isnothing(
        GaussletBases._qwrg_try_nested_factorized_parent_basis(nonfactorized_fixed_block),
    )
    @test GaussletBases._nested_by_center_sidecar_path(nonfactorized_fixed_block) ==
          :general_parent_dense
    nested_bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(
        fixed_block.parent_basis,
        expansion;
        gausslet_backend = :numerical_reference,
    )
    fallback_nuclear = GaussletBases._qwrg_bond_aligned_nested_fixed_block_nuclear_one_body_by_center(
        fixed_block.parent_basis,
        nonfactorized_fixed_block,
        nested_bundles.bundle_x,
        nested_bundles.bundle_y,
        nested_bundles.bundle_z,
        expansion,
    )
    reference_nuclear = GaussletBases._qwrg_bond_aligned_general_contracted_nuclear_one_body_by_center(
        fixed_block.parent_basis,
        nonfactorized_coefficients,
        nested_bundles.bundle_x,
        nested_bundles.bundle_y,
        nested_bundles.bundle_z,
        expansion,
    )
    @test length(fallback_nuclear) == 2
    for nucleus_index in eachindex(fallback_nuclear, reference_nuclear)
        @test all(isfinite, fallback_nuclear[nucleus_index])
        @test fallback_nuclear[nucleus_index] ≈
              transpose(fallback_nuclear[nucleus_index]) atol = 1.0e-12 rtol = 1.0e-12
        @test fallback_nuclear[nucleus_index] ≈
              reference_nuclear[nucleus_index] atol = 1.0e-12 rtol = 1.0e-12
    end

    representatives_before_sidecar = Matrix{Float64}(nonfactorized_fixed_block.coefficient_matrix)
    staged_sidecar = GaussletBases._nested_attach_staged_by_center_sidecar!(
        nonfactorized_fixed_block;
        provenance = (; test = :diatomic_nonfactorized_by_center),
    )
    @test staged_sidecar isa GaussletBases._CartesianNestedStagedByCenterSidecar3D
    @test GaussletBases._nested_by_center_sidecar_path(nonfactorized_fixed_block) ==
          :staged_factorized
    @test isnothing(
        GaussletBases._qwrg_try_nested_factorized_parent_basis(nonfactorized_fixed_block),
    )
    @test Matrix{Float64}(nonfactorized_fixed_block.coefficient_matrix) ==
          representatives_before_sidecar
    @test size(nonfactorized_fixed_block.coefficient_matrix, 2) ==
          staged_sidecar.diagnostics.final_dimension
    @test staged_sidecar.diagnostics.block_count >= 1
    @test staged_sidecar.diagnostics.max_support_count <=
          size(nonfactorized_fixed_block.coefficient_matrix, 1)

    staged_nuclear = GaussletBases._qwrg_bond_aligned_nested_fixed_block_nuclear_one_body_by_center(
        fixed_block.parent_basis,
        nonfactorized_fixed_block,
        nested_bundles.bundle_x,
        nested_bundles.bundle_y,
        nested_bundles.bundle_z,
        expansion,
    )
    @test length(staged_nuclear) == length(reference_nuclear)
    for nucleus_index in eachindex(staged_nuclear, reference_nuclear)
        @test all(isfinite, staged_nuclear[nucleus_index])
        @test staged_nuclear[nucleus_index] ≈
              transpose(staged_nuclear[nucleus_index]) atol = 1.0e-12 rtol = 1.0e-12
        @test staged_nuclear[nucleus_index] ≈
              reference_nuclear[nucleus_index] atol = 1.0e-10 rtol = 1.0e-10
    end

    function _clone_fixed_block_with_sidecars(template, factorized_cache, staged_cache)
        return GaussletBases._NestedFixedBlock3D(
            template.parent_basis,
            template.shell,
            template.gausslet_backend,
            template.coefficient_matrix,
            template.support_indices,
            template.overlap,
            template.kinetic,
            template.position_x,
            template.position_y,
            template.position_z,
            template.x2_x,
            template.x2_y,
            template.x2_z,
            template.weights,
            template.gaussian_sum,
            template.pair_sum,
            template.fixed_centers,
            factorized_cache,
            staged_cache,
        )
    end

    endcap_nested = bond_aligned_diatomic_nested_fixed_block(
        basis;
        expansion,
        shared_shell_layer_policy = :endcap_panel_owned,
        shared_shell_endcap_panel_q = 4,
        shared_shell_endcap_panel_L = 4,
    )
    endcap_fixed_block = endcap_nested.fixed_block
    product_sidecar = endcap_fixed_block.staged_by_center_sidecar[]
    @test isnothing(fixed_block.staged_by_center_sidecar[])
    @test product_sidecar isa GaussletBases._CartesianNestedProductStagedByCenterSidecar3D
    @test product_sidecar.diagnostics.product_unit_count >= 6
    @test product_sidecar.diagnostics.product_unit_count % 6 == 0
    @test product_sidecar.diagnostics.generic_unit_count >= 1
    @test product_sidecar.diagnostics.final_dimension == size(endcap_fixed_block.coefficient_matrix, 2)
    if endcap_fixed_block.factorized_cartesian_parent_basis[] isa
       GaussletBases._CartesianNestedFactorizedBasis3D
        @test GaussletBases._nested_by_center_sidecar_path(endcap_fixed_block) == :factorized_final
    end

    endcap_bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(
        endcap_fixed_block.parent_basis,
        expansion;
        gausslet_backend = :numerical_reference,
    )
    endcap_dense_fixed_block = _clone_fixed_block_with_sidecars(
        endcap_fixed_block,
        GaussletBases._nested_factorized_basis_cache(),
        GaussletBases._nested_staged_by_center_sidecar_cache(),
    )
    @test GaussletBases._nested_by_center_sidecar_path(endcap_dense_fixed_block) ==
          :general_parent_dense
    endcap_dense_nuclear =
        GaussletBases._qwrg_bond_aligned_general_contracted_nuclear_one_body_by_center(
            endcap_fixed_block.parent_basis,
            endcap_fixed_block.coefficient_matrix,
            endcap_bundles.bundle_x,
            endcap_bundles.bundle_y,
            endcap_bundles.bundle_z,
            expansion,
        )

    endcap_generic_fixed_block = _clone_fixed_block_with_sidecars(
        endcap_fixed_block,
        GaussletBases._nested_factorized_basis_cache(),
        GaussletBases._nested_staged_by_center_sidecar_cache(),
    )
    endcap_generic_sidecar = GaussletBases._nested_attach_staged_by_center_sidecar!(
        endcap_generic_fixed_block;
        provenance = (; test = :diatomic_endcap_generic_staged_by_center),
    )
    @test endcap_generic_sidecar isa GaussletBases._CartesianNestedStagedByCenterSidecar3D
    @test GaussletBases._nested_by_center_sidecar_path(endcap_generic_fixed_block) ==
          :staged_factorized
    endcap_generic_nuclear =
        GaussletBases._qwrg_bond_aligned_nested_fixed_block_nuclear_one_body_by_center(
            endcap_fixed_block.parent_basis,
            endcap_generic_fixed_block,
            endcap_bundles.bundle_x,
            endcap_bundles.bundle_y,
            endcap_bundles.bundle_z,
            expansion,
        )

    endcap_product_fixed_block = _clone_fixed_block_with_sidecars(
        endcap_fixed_block,
        GaussletBases._nested_factorized_basis_cache(),
        GaussletBases._nested_staged_by_center_sidecar_cache(product_sidecar),
    )
    @test GaussletBases._nested_by_center_sidecar_path(endcap_product_fixed_block) ==
          :product_staged_factorized
    endcap_product_nuclear =
        GaussletBases._qwrg_bond_aligned_nested_fixed_block_nuclear_one_body_by_center(
            endcap_fixed_block.parent_basis,
            endcap_product_fixed_block,
            endcap_bundles.bundle_x,
            endcap_bundles.bundle_y,
            endcap_bundles.bundle_z,
            expansion,
        )

    for nucleus_index in eachindex(endcap_dense_nuclear)
        @test endcap_generic_nuclear[nucleus_index] ≈
              endcap_dense_nuclear[nucleus_index] atol = 1.0e-10 rtol = 1.0e-10
        @test endcap_product_nuclear[nucleus_index] ≈
              endcap_dense_nuclear[nucleus_index] atol = 1.0e-10 rtol = 1.0e-10
    end
    bad_product_sidecar = GaussletBases._CartesianNestedProductStagedByCenterSidecar3D(
        (1, 1, 1),
        product_sidecar.units,
        product_sidecar.provenance,
        product_sidecar.diagnostics,
    )
    bad_product_fixed_block = _clone_fixed_block_with_sidecars(
        endcap_fixed_block,
        GaussletBases._nested_factorized_basis_cache(),
        GaussletBases._nested_staged_by_center_sidecar_cache(bad_product_sidecar),
    )
    @test_throws ArgumentError GaussletBases._qwrg_bond_aligned_nested_fixed_block_nuclear_one_body_by_center(
        endcap_fixed_block.parent_basis,
        bad_product_fixed_block,
        endcap_bundles.bundle_x,
        endcap_bundles.bundle_y,
        endcap_bundles.bundle_z,
        expansion,
    )
    @test_throws ArgumentError GaussletBases._nested_attach_staged_by_center_sidecar!(
        nonfactorized_fixed_block;
        block_column_ranges = [1:1],
        replace = true,
    )
    @test_throws ArgumentError GaussletBases._nested_attach_staged_by_center_sidecar!(
        nonfactorized_fixed_block;
        block_column_ranges = [1:(size(nonfactorized_coefficients, 2) + 1)],
        replace = true,
    )
    saved_sidecar = nonfactorized_fixed_block.staged_by_center_sidecar[]
    nonfactorized_fixed_block.staged_by_center_sidecar[] = (; invalid = true)
    @test_throws ArgumentError GaussletBases._nested_staged_by_center_sidecar(
        nonfactorized_fixed_block,
    )
    @test_throws ArgumentError GaussletBases._qwrg_bond_aligned_nested_fixed_block_nuclear_one_body_by_center(
        fixed_block.parent_basis,
        nonfactorized_fixed_block,
        nested_bundles.bundle_x,
        nested_bundles.bundle_y,
        nested_bundles.bundle_z,
        expansion,
    )
    nonfactorized_fixed_block.staged_by_center_sidecar[] = saved_sidecar

    @test assembled_one_body_hamiltonian(full_nested_ops_localized) ≈
          full_nested_ops_localized.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test assembled_one_body_hamiltonian(full_nested_ops_localized; nuclear_charges = [1.0, 0.0]) ≈
          atom_a_nested_ops_localized.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
    @test assembled_one_body_hamiltonian(full_nested_ops_localized; nuclear_charges = [0.0, 1.0]) ≈
          atom_b_nested_ops_localized.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
end

@testset "Nuclear term storage auto stays lightweight on longer-center routes" begin
    basis, _operators, _diagnostics = _bond_aligned_homonuclear_chain_qw_fixture(; natoms = 3)
    auto_ops = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = fill(1.0, length(basis.nuclei)),
        nuclear_term_storage = :auto,
        interaction_treatment = :ggt_nearest,
    )

    @test auto_ops.nuclear_term_storage == :total_only
    @test isnothing(auto_ops.kinetic_one_body)
    @test isnothing(auto_ops.nuclear_one_body_by_center)
    @test_throws ArgumentError assembled_one_body_hamiltonian(
        auto_ops;
        nuclear_charges = [1.0, 0.0, 1.0],
    )
end

@testset "Ordinary Cartesian naming surface distinguishes geometry from supplement" begin
    @test OrdinaryCartesianOrbital3D === QiuWhiteHybridOrbital3D
    @test OrdinaryCartesianOperators3D === QiuWhiteResidualGaussianOperators

    (
        direct_basis,
        _direct_parent_ops,
        _direct_parent_check,
        direct_supplement,
        direct_hybrid_ops,
        _direct_hybrid_check,
    ) = _bond_aligned_diatomic_hybrid_qw_fixture(; bond_length = 1.4)

    direct_via_clear_name = ordinary_cartesian_product_operators(
        direct_basis,
        direct_supplement;
        nuclear_charges = [1.0, 1.0],
        interaction_treatment = :ggt_nearest,
    )

    @test direct_via_clear_name isa OrdinaryCartesianOperators3D
    @test direct_via_clear_name.overlap ≈ direct_hybrid_ops.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test direct_via_clear_name.one_body_hamiltonian ≈ direct_hybrid_ops.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test direct_via_clear_name.interaction_matrix ≈ direct_hybrid_ops.interaction_matrix atol = 1.0e-12 rtol = 1.0e-12
    @test direct_via_clear_name.residual_count == direct_hybrid_ops.residual_count

    (
        _nested_basis,
        _nested_parent_ops,
        _nested_parent_check,
        nested_expansion,
        _nested_source,
        nested_fixed_block,
        nested_pure_ops,
        _nested_pure_check,
    ) = _bond_aligned_diatomic_nested_qw_fixture(; bond_length = 1.4)

    nested_via_clear_name = nested_cartesian_operators(
        nested_fixed_block;
        nuclear_charges = [1.0, 1.0],
        expansion = nested_expansion,
        interaction_treatment = :ggt_nearest,
    )

    @test nested_via_clear_name isa OrdinaryCartesianOperators3D
    @test nested_via_clear_name.gaussian_data === nothing
    @test nested_via_clear_name.overlap ≈ nested_pure_ops.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test nested_via_clear_name.one_body_hamiltonian ≈ nested_pure_ops.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test nested_via_clear_name.interaction_matrix ≈ nested_pure_ops.interaction_matrix atol = 1.0e-12 rtol = 1.0e-12
    @test nested_via_clear_name.residual_count == 0
end

@testset "Bond-aligned diatomic split geometry" begin
    basis, _operators, _check = _bond_aligned_diatomic_qw_fixture(; bond_length = 1.4)
    expansion = coulomb_gaussian_expansion(doacc = false)
    bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(basis, expansion)
    parent_box = (
        1:length(basis.basis_x),
        1:length(basis.basis_y),
        1:length(basis.basis_z),
    )
    working_box = (2:8, 2:8, 2:12)
    midpoint = 0.0

    geometry = GaussletBases._nested_bond_aligned_diatomic_split_geometry(
        bundles,
        parent_box,
        working_box;
        bond_axis = :z,
        midpoint,
        nside = 5,
        min_parallel_to_transverse_ratio = 0.4,
    )
    sliver_geometry = GaussletBases._nested_bond_aligned_diatomic_split_geometry(
        bundles,
        parent_box,
        working_box;
        bond_axis = :z,
        midpoint,
        nside = 5,
        min_parallel_to_transverse_ratio = 0.75,
    )
    short_geometry = GaussletBases._nested_bond_aligned_diatomic_split_geometry(
        bundles,
        parent_box,
        (3:7, 3:7, 4:12);
        bond_axis = :z,
        midpoint,
        nside = 5,
        min_parallel_to_transverse_ratio = 0.4,
    )
    split_basis = bond_aligned_homonuclear_qw_basis(
        family = :G10,
        bond_length = 1.4,
        core_spacing = 0.3,
        xmax_parallel = 8.0,
        xmax_transverse = 6.0,
        bond_axis = :z,
        nuclear_charge = 2.0,
    )
    split_diagnostics = bond_aligned_diatomic_nested_geometry_diagnostics(
        split_basis;
        nside = 5,
    )
    split_bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(split_basis, expansion)
    split_parent_box = (
        1:length(split_basis.basis_x),
        1:length(split_basis.basis_y),
        1:length(split_basis.basis_z),
    )
    split_geometry = GaussletBases._nested_bond_aligned_diatomic_split_geometry(
        split_bundles,
        split_parent_box,
        split_diagnostics.geometry.working_box;
        bond_axis = :z,
        midpoint,
        nside = 5,
        min_parallel_to_transverse_ratio = 0.4,
    )

    @test geometry.count_eligible
    @test !geometry.unsplit_aspect_eligible
    @test geometry.shape_eligible
    @test !geometry.did_split
    @test geometry.split_index == 7
    @test geometry.working_box == working_box
    @test isnothing(geometry.shared_midpoint_box)
    @test geometry.child_boxes == [(2:8, 2:8, 2:6), (2:8, 2:8, 8:12)]
    @test all(widths[3] >= 0.4 * max(widths[1], widths[2]) for widths in geometry.child_physical_widths)

    @test split_geometry.count_eligible
    @test split_geometry.unsplit_aspect_eligible
    @test split_geometry.shape_eligible
    @test split_geometry.did_split
    @test !isnothing(split_geometry.shared_midpoint_box)
    @test length(split_geometry.child_boxes) == 2

    @test sliver_geometry.count_eligible
    @test !sliver_geometry.shape_eligible
    @test !sliver_geometry.did_split

    @test !short_geometry.count_eligible
    @test !short_geometry.did_split
end

function _bond_aligned_diatomic_shared_shell_policy_basis(nside::Int)
    return bond_aligned_homonuclear_qw_basis(
        family = :G10,
        bond_length = 4.0,
        core_spacing = 1.2 / (4.0 * (nside - 3)),
        xmax_parallel = 20.0,
        xmax_transverse = 20.0,
        bond_axis = :z,
        nuclear_charge = 4.0,
        reference_spacing = 1.0,
        tail_spacing = 10.0,
    )
end

function _bond_aligned_diatomic_shared_shell_ns_local_counts(
    nside::Int;
    angular_resolution_scale::Float64 = 1.4,
    max_layers::Union{Nothing,Int} = nothing,
)
    basis = _bond_aligned_diatomic_shared_shell_policy_basis(nside)
    expansion = coulomb_gaussian_expansion(doacc = false)
    bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(basis, expansion)
    dims = GaussletBases._nested_axis_lengths(bundles)
    current_box = (1:dims[1], 1:dims[2], 1:dims[3])
    retention = GaussletBases._nested_resolve_complete_shell_retention(nside)
    counts = NTuple{3,Int}[]
    while minimum(length.(current_box)) > nside && GaussletBases._nested_can_shrink_box(current_box)
        inner_box = GaussletBases._nested_inner_box(current_box)
        adaptive = GaussletBases._nested_diatomic_adaptive_shell_retention(
            basis,
            bundles,
            current_box,
            inner_box,
            retention;
            nside,
            shared_shell_angular_resolution_scale = angular_resolution_scale,
        )
        push!(
            counts,
            (
                adaptive.chosen_x.retain + 2,
                adaptive.chosen_y.retain + 2,
                adaptive.chosen_z.retain + 2,
            ),
        )
        current_box = inner_box
        if !isnothing(max_layers) && length(counts) >= max_layers
            break
        end
    end
    return counts
end

@testset "Bond-aligned diatomic shared-shell angular policy defaults to scaled reference" begin
    basis6 = _bond_aligned_diatomic_shared_shell_policy_basis(6)
    expansion6 = coulomb_gaussian_expansion(doacc = false)
    bundles6 = GaussletBases._qwrg_bond_aligned_axis_bundles(basis6, expansion6)
    dims6 = GaussletBases._nested_axis_lengths(bundles6)
    parent_box6 = (1:dims6[1], 1:dims6[2], 1:dims6[3])

    reference6 = GaussletBases._nested_diatomic_reference_band(
        basis6,
        bundles6,
        parent_box6;
        nside = 6,
        reference_fudge_factor = 1.0,
    )
    scaled_reference6 = GaussletBases._nested_diatomic_shared_shell_reference_band(
        basis6,
        bundles6,
        parent_box6;
        nside = 6,
        angular_resolution_scale = 1.4,
    )

    @test scaled_reference6.theta_min ≈ 1.4 * reference6.ideal_theta_min atol = 1.0e-12 rtol = 1.0e-12
    @test scaled_reference6.theta_max ≈ 1.4 * reference6.ideal_theta_max atol = 1.0e-12 rtol = 1.0e-12
    @test scaled_reference6.reference_bounds == reference6.reference_bounds
    @test scaled_reference6.reference_retain == reference6.reference_retain

    for nside in 5:8
        @test only(
            _bond_aligned_diatomic_shared_shell_ns_local_counts(
                nside;
                angular_resolution_scale = 1.4,
                max_layers = 1,
            ),
        ) == (nside, nside, nside)
    end

    @test only(
        _bond_aligned_diatomic_shared_shell_ns_local_counts(
            6;
            angular_resolution_scale = 1.0,
            max_layers = 1,
        ),
    ) == (8, 8, 8)
    @test only(
        _bond_aligned_diatomic_shared_shell_ns_local_counts(
            6;
            angular_resolution_scale = 1.2,
            max_layers = 1,
        ),
    ) == (7, 7, 7)

    nside6_counts = _bond_aligned_diatomic_shared_shell_ns_local_counts(
        6;
        angular_resolution_scale = 1.4,
        max_layers = 6,
    )
    @test length(nside6_counts) >= 2
    @test nside6_counts[1] == (6, 6, 6)
    @test all(count[1] == 6 && count[2] == 6 for count in nside6_counts)
    @test issorted([count[3] for count in nside6_counts])
    @test any(count[3] > 6 for count in nside6_counts[2:end])
end

@testset "Bond-aligned diatomic adaptive retained counts are parity-neutral" begin
    @test GaussletBases._nested_diatomic_candidate_counts(11, 4) == collect(4:11)
    @test GaussletBases._nested_diatomic_candidate_counts(11, 5) == collect(5:11)

    symmetric_centers = [-2.0, -1.0, 0.0, 1.0, 2.0]
    @test GaussletBases._nested_doside_retained_count(symmetric_centers, 4) == 3
    @test GaussletBases._nested_doside_retained_count(
        symmetric_centers,
        4;
        enforce_symmetric_odd = false,
    ) == 4

end

@testset "Bond-aligned heteronuclear split geometry" begin
    centers_axis = collect(-5.0:1.0:5.0)
    @test GaussletBases._nested_diatomic_split_plane_index(
        centers_axis,
        1:11,
        0.0;
        prefer_midpoint_tie_side = :left,
    ) == 6
    @test GaussletBases._nested_diatomic_split_plane_index(
        centers_axis,
        1:11,
        0.0;
        prefer_midpoint_tie_side = :right,
    ) == 5

    basis = bond_aligned_heteronuclear_qw_basis(
        atoms = ("He", "H"),
        bond_length = 1.45,
        core_spacings = (0.25, 0.5),
        nuclear_charges = (2.0, 1.0),
        xmax_parallel = 6.0,
        xmax_transverse = 4.0,
        bond_axis = :z,
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(basis, expansion)
    parent_box = (
        1:length(basis.basis_x),
        1:length(basis.basis_y),
        1:length(basis.basis_z),
    )
    working_box = (2:(length(basis.basis_x) - 1), 2:(length(basis.basis_y) - 1), 2:(length(basis.basis_z) - 1))
    midpoint = sum(GaussletBases._qwrg_axis_coordinate(nucleus, :z) for nucleus in basis.nuclei) / 2
    geometry = GaussletBases._nested_bond_aligned_diatomic_split_geometry(
        bundles,
        parent_box,
        working_box;
        bond_axis = :z,
        midpoint,
        nside = 5,
        min_parallel_to_transverse_ratio = 0.4,
        use_midpoint_slab = false,
        prefer_midpoint_tie_side = :left,
    )

    @test geometry.count_eligible
    @test !geometry.unsplit_aspect_eligible
    @test geometry.shape_eligible
    @test !geometry.did_split
    @test geometry.shared_midpoint_box === nothing
    @test length(geometry.child_boxes) == 2
    @test sum(length(box[3]) for box in geometry.child_boxes) == length(working_box[3])
    @test all(widths[3] >= 0.4 * max(widths[1], widths[2]) for widths in geometry.child_physical_widths)
end

@testset "Bond-aligned diatomic nested fixed block" begin
    (
        basis,
        operators,
        check,
        expansion,
        source,
        fixed_block,
        parent_modes,
        _parent_ground,
        _projected,
        projected_vee,
        capture,
        projected_energy,
    ) = _bond_aligned_diatomic_nested_fixed_block_fixture(; bond_length = 1.4)

    @test source isa GaussletBases._CartesianNestedBondAlignedDiatomicSource3D
    @test fixed_block isa GaussletBases._NestedFixedBlock3D
    @test !source.geometry.did_split
    @test length(source.shared_shell_layers) >= 1
    @test length(source.child_sequences) == 1
    @test isnothing(source.geometry.shared_midpoint_box)
    @test isnothing(source.midpoint_slab_column_range)
    @test length(source.child_column_ranges) == 1
    @test length(source.child_column_ranges[1]) == size(source.child_sequences[1].coefficient_matrix, 2)
    @test size(fixed_block.overlap, 1) == size(source.sequence.coefficient_matrix, 2)
    @test source.sequence.working_box == (
        1:length(basis.basis_x),
        1:length(basis.basis_y),
        1:length(basis.basis_z),
    )
    @test length(source.sequence.support_indices) ==
        length(basis.basis_x) * length(basis.basis_y) * length(basis.basis_z)
    @test size(fixed_block.coefficient_matrix, 1) ==
        length(basis.basis_x) * length(basis.basis_y) * length(basis.basis_z)
    @test size(fixed_block.coefficient_matrix, 2) < size(fixed_block.coefficient_matrix, 1)
    @test !hasproperty(source.sequence.packet, :term_storage)
    @test !hasproperty(source.sequence.packet, :gaussian_terms)
    @test !hasproperty(source.sequence.packet, :pair_terms)
    @test !hasproperty(fixed_block, :term_storage)
    @test !hasproperty(fixed_block, :gaussian_terms)
    @test !hasproperty(fixed_block, :pair_terms)
    @test !isnothing(fixed_block.gaussian_sum)
    @test !isnothing(fixed_block.pair_sum)
    @test norm(fixed_block.overlap - I, Inf) < 1.0e-10
    @test all(isfinite, fixed_block.weights)
    @test minimum(fixed_block.weights) > 0.0
    @test all(isfinite, fixed_block.fixed_centers)
    @test capture > 0.998
    @test projected_energy < -1.2
    @test abs(projected_energy - parent_modes.values[1]) < 0.03
    @test 0.7 < projected_vee < 0.8
    @test abs(projected_vee - check.vee_expectation) < 5.0e-4
end

@testset "Bond-aligned diatomic compact nested fixed-block contract" begin
    (
        basis,
        _operators,
        _check,
        expansion,
        default_source,
        default_fixed_block,
        _parent_modes,
        _parent_ground,
        _projected,
        _projected_vee,
        _capture,
        _projected_energy,
    ) = _bond_aligned_diatomic_nested_fixed_block_fixture(; bond_length = 1.4)

    term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
    explicit_nested = bond_aligned_diatomic_nested_fixed_block(
        basis;
        expansion,
        term_coefficients,
    )
    explicit_fixed_block = explicit_nested.fixed_block

    @test !hasproperty(default_source.sequence.packet, :term_storage)
    @test !hasproperty(default_source.sequence.packet, :gaussian_terms)
    @test !hasproperty(default_source.sequence.packet, :pair_terms)
    @test !hasproperty(default_fixed_block, :term_storage)
    @test !hasproperty(default_fixed_block, :gaussian_terms)
    @test !hasproperty(default_fixed_block, :pair_terms)
    @test !isnothing(default_fixed_block.gaussian_sum)
    @test !isnothing(default_fixed_block.pair_sum)

    @test !hasproperty(explicit_fixed_block, :term_storage)
    @test !hasproperty(explicit_fixed_block, :gaussian_terms)
    @test !hasproperty(explicit_fixed_block, :pair_terms)
    @test !isnothing(explicit_fixed_block.gaussian_sum)
    @test !isnothing(explicit_fixed_block.pair_sum)

    @test !hasproperty(explicit_nested.source.sequence.packet, :term_storage)
    @test !hasproperty(explicit_nested.source.sequence.packet, :gaussian_terms)
    @test !hasproperty(explicit_nested.source.sequence.packet, :pair_terms)

    @test norm(default_fixed_block.overlap - explicit_fixed_block.overlap, Inf) < 1.0e-12
    @test norm(default_fixed_block.coefficient_matrix - explicit_fixed_block.coefficient_matrix, Inf) < 1.0e-12
    @test norm(default_fixed_block.gaussian_sum - explicit_fixed_block.gaussian_sum, Inf) < 1.0e-12
    @test norm(default_fixed_block.pair_sum - explicit_fixed_block.pair_sum, Inf) < 1.0e-12

    @test GaussletBases._qwrg_fixed_block_one_body_matrix(default_fixed_block, expansion; Z = 1.0) ≈
        GaussletBases._qwrg_fixed_block_one_body_matrix(explicit_fixed_block, expansion; Z = 1.0) atol =
        1.0e-10 rtol = 1.0e-10
    @test GaussletBases._qwrg_fixed_block_interaction_matrix(default_fixed_block, expansion) ≈
        GaussletBases._qwrg_fixed_block_interaction_matrix(explicit_fixed_block, expansion) atol =
        1.0e-10 rtol = 1.0e-10

end

@testset "Bond-aligned diatomic packet-kernel parity" begin
    basis = bond_aligned_homonuclear_qw_basis(
        bond_length = 1.4,
        core_spacing = 0.5,
        xmax_parallel = 4.0,
        xmax_transverse = 3.0,
        bond_axis = :z,
    )

    reference = bond_aligned_diatomic_nested_fixed_block(
        basis;
        nside = 5,
        packet_kernel = :support_reference,
    )
    direct = bond_aligned_diatomic_nested_fixed_block(
        basis;
        nside = 5,
        packet_kernel = :factorized_direct,
    )

    reference_source = reference.source
    direct_source = direct.source
    reference_block = reference.fixed_block
    direct_block = direct.fixed_block
    atol = 1.0e-10
    rtol = 1.0e-10

    @test size(direct_block.coefficient_matrix, 2) == size(reference_block.coefficient_matrix, 2)
    @test direct_source.sequence.support_indices == reference_source.sequence.support_indices
    @test direct_source.sequence.coefficient_matrix ≈
          reference_source.sequence.coefficient_matrix atol = atol rtol = rtol
    @test direct_block.coefficient_matrix ≈
          reference_block.coefficient_matrix atol = atol rtol = rtol

    @test direct_block.overlap ≈ reference_block.overlap atol = atol rtol = rtol
    @test direct_block.kinetic ≈ reference_block.kinetic atol = atol rtol = rtol
    @test direct_block.position_x ≈ reference_block.position_x atol = atol rtol = rtol
    @test direct_block.position_y ≈ reference_block.position_y atol = atol rtol = rtol
    @test direct_block.position_z ≈ reference_block.position_z atol = atol rtol = rtol
    @test direct_block.x2_x ≈ reference_block.x2_x atol = atol rtol = rtol
    @test direct_block.x2_y ≈ reference_block.x2_y atol = atol rtol = rtol
    @test direct_block.x2_z ≈ reference_block.x2_z atol = atol rtol = rtol
    @test direct_block.gaussian_sum ≈ reference_block.gaussian_sum atol = atol rtol = rtol
    @test direct_block.pair_sum ≈ reference_block.pair_sum atol = atol rtol = rtol
end

@testset "Bond-aligned diatomic packet-kernel default" begin
    basis = bond_aligned_homonuclear_qw_basis(
        bond_length = 1.4,
        core_spacing = 0.5,
        xmax_parallel = 4.0,
        xmax_transverse = 3.0,
        bond_axis = :z,
    )

    default_context = GaussletBases._normalized_nested_source_frontend_context(
        basis;
        nside = 5,
    )
    @test default_context.build_options.packet_kernel == :factorized_direct

    default_source = bond_aligned_diatomic_nested_fixed_source(
        basis;
        nside = 5,
    )
    direct_source = bond_aligned_diatomic_nested_fixed_source(
        basis;
        nside = 5,
        packet_kernel = :factorized_direct,
    )
    support_source = bond_aligned_diatomic_nested_fixed_source(
        basis;
        nside = 5,
        packet_kernel = :support_reference,
    )

    @test size(default_source.sequence.coefficient_matrix, 2) ==
          size(direct_source.sequence.coefficient_matrix, 2)
    @test default_source.sequence.support_indices == direct_source.sequence.support_indices
    @test default_source.sequence.coefficient_matrix ≈
          direct_source.sequence.coefficient_matrix atol = 1.0e-12 rtol = 1.0e-12

    @test !isnothing(support_source.sequence.packet)
    @test size(support_source.sequence.coefficient_matrix, 2) ==
          size(default_source.sequence.coefficient_matrix, 2)
end

@testset "Bond-aligned diatomic nested QW consumer path" begin
    (
        _basis,
        parent_ops,
        parent_check,
        _expansion,
        _source,
        fixed_block,
        nested_ops,
        nested_check,
    ) = _bond_aligned_diatomic_nested_qw_fixture(; bond_length = 1.4)

    @test nested_ops isa GaussletBases.QiuWhiteResidualGaussianOperators
    @test nested_ops.basis === fixed_block
    @test nested_ops.gaussian_data === nothing
    @test nested_ops.interaction_treatment == :ggt_nearest
    @test nested_ops.residual_count == 0
    @test nested_ops.gausslet_count == size(fixed_block.overlap, 1)
    @test size(nested_ops.residual_centers) == (0, 3)
    @test size(nested_ops.residual_widths) == (0, 3)
    @test norm(nested_ops.overlap - I, Inf) < 1.0e-10
    @test norm(nested_ops.one_body_hamiltonian - transpose(nested_ops.one_body_hamiltonian), Inf) < 1.0e-12
    @test norm(nested_ops.interaction_matrix - transpose(nested_ops.interaction_matrix), Inf) < 1.0e-12
    @test all(orbital.kind == :nested_fixed for orbital in orbitals(nested_ops))
    @test all(isfinite, fixed_block.weights)
    @test minimum(fixed_block.weights) > 0.0
    @test nested_check.overlap_error < 1.0e-10
    @test nested_check.orbital_energy < -1.2
    @test 0.7 < nested_check.vee_expectation < 0.8
    @test abs(nested_check.orbital_energy - parent_check.orbital_energy) < 0.03
    @test abs(nested_check.vee_expectation - parent_check.vee_expectation) < 0.03
    @test parent_ops.residual_count == 0
end

@testset "Bond-aligned diatomic molecular supplement ordinary QW path" begin
    fixture = _bond_aligned_diatomic_hybrid_qw_fixture(; bond_length = 1.4)
    if fixture === nothing
        @test true
    else
        (
            _basis,
            parent_ops,
            parent_check,
            supplement,
            ordinary_ops,
            ordinary_check,
        ) = fixture

        @test supplement isa LegacyBondAlignedDiatomicGaussianSupplement
        @test supplement.atomic_source.lmax == 1
        @test ordinary_ops.gaussian_data === supplement
        @test ordinary_ops.interaction_treatment == :ggt_nearest
        @test ordinary_ops.residual_count > 0
        @test size(ordinary_ops.residual_centers, 1) == ordinary_ops.residual_count
        @test size(ordinary_ops.residual_widths, 1) == ordinary_ops.residual_count
        @test norm(ordinary_ops.overlap - I, Inf) < 1.0e-10
        @test ordinary_check.overlap_error < 1.0e-10
        @test ordinary_check.orbital_energy < -1.25
        @test 0.75 < ordinary_check.vee_expectation < 0.81
        @test abs(ordinary_check.orbital_energy - parent_check.orbital_energy) < 0.01
        @test abs(ordinary_check.vee_expectation - parent_check.vee_expectation) < 0.01
        @test parent_ops.residual_count == 0
    end
end

@testset "Bond-aligned heteronuclear molecular supplement ordinary QW path" begin
    (
        basis,
        parent_ops,
        parent_check,
        supplement,
        ordinary_ops,
        ordinary_check,
    ) = _bond_aligned_heteronuclear_hybrid_qw_fixture(; bond_length = 1.45)

    @test supplement isa LegacyBondAlignedHeteronuclearGaussianSupplement
    @test supplement.atomic_sources[1].atom == "He"
    @test supplement.atomic_sources[2].atom == "H"
    @test supplement.atomic_sources[1].basis_name == "cc-pVTZ"
    @test supplement.atomic_sources[2].basis_name == "cc-pVTZ"
    @test supplement.atomic_sources[1].lmax == 1
    @test supplement.atomic_sources[2].lmax == 1
    @test ordinary_ops.gaussian_data === supplement
    @test ordinary_ops.interaction_treatment == :ggt_nearest
    @test ordinary_ops.residual_count > 0
    @test size(ordinary_ops.residual_centers, 1) == ordinary_ops.residual_count
    @test size(ordinary_ops.residual_widths, 1) == ordinary_ops.residual_count
    @test norm(ordinary_ops.overlap - I, Inf) < 1.0e-10
    @test ordinary_check.overlap_error < 1.0e-10
    @test ordinary_check.orbital_energy < -2.0
    @test 0.8 < ordinary_check.vee_expectation < 1.5
    @test ordinary_check.orbital_energy < parent_check.orbital_energy
    @test ordinary_check.vee_expectation > parent_check.vee_expectation

    payload = bond_aligned_diatomic_geometry_payload(ordinary_ops)
    slice = bond_aligned_diatomic_plane_slice(
        payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-12,
    )
    @test payload.bond_axis == :z
    @test length(payload.nuclei) == 2
    @test count(point -> point.group_kind == :residual_gaussian, payload.points) == ordinary_ops.residual_count
    @test slice.selected_count > 0
    @test all(abs(point.y) <= slice.plane_tol for point in slice.points)
    @test parent_ops.residual_count == 0
end

@testset "Bond-aligned heteronuclear nested fixed block" begin
    (
        basis,
        parent_ops,
        parent_check,
        _supplement,
        ordinary_ops,
        ordinary_check,
        _expansion,
        source,
        fixed_block,
        _parent_modes,
        _parent_ground,
        _projected,
        projected_vee,
        capture,
        projected_energy,
    ) = _bond_aligned_heteronuclear_nested_fixed_block_fixture(; bond_length = 1.45)

    @test source isa GaussletBases._CartesianNestedBondAlignedDiatomicSource3D
    @test fixed_block isa GaussletBases._NestedFixedBlock3D
    @test !source.geometry.did_split
    @test source.geometry.shared_midpoint_box === nothing
    @test isnothing(source.midpoint_slab_column_range)
    @test length(source.shared_shell_layers) >= 1
    @test length(source.child_sequences) == 1
    @test length(source.child_column_ranges) == 1
    @test source.sequence.working_box == (
        1:length(basis.basis_x),
        1:length(basis.basis_y),
        1:length(basis.basis_z),
    )
    @test norm(fixed_block.overlap - I, Inf) < 1.0e-10
    @test all(isfinite, fixed_block.weights)
    @test minimum(fixed_block.weights) > 0.0
    @test all(isfinite, fixed_block.fixed_centers)
    @test parent_check.overlap_error < 1.0e-10
    @test ordinary_check.overlap_error < 1.0e-10
    @test capture > 0.99
    @test projected_energy < -2.0
    @test projected_vee > 0.8
    @test size(fixed_block.overlap, 1) < ordinary_ops.gausslet_count
    @test size(fixed_block.overlap, 1) < parent_ops.gausslet_count
    @test source.geometry.child_boxes[1][3] != source.geometry.child_boxes[2][3]
end

@testset "Bond-aligned heteronuclear nested QW consumer path" begin
    (
        _basis,
        parent_ops,
        parent_check,
        source,
        fixed_block,
        supplement,
        nested_ops,
        nested_check,
    ) = _bond_aligned_heteronuclear_nested_hybrid_qw_fixture(; bond_length = 1.45)

    @test nested_ops isa GaussletBases.QiuWhiteResidualGaussianOperators
    @test nested_ops.basis === fixed_block
    @test nested_ops.gaussian_data === supplement
    @test nested_ops.interaction_treatment == :ggt_nearest
    @test nested_ops.residual_count > 0
    @test nested_check.overlap_error < 1.0e-8
    @test nested_check.orbital_energy < -2.0
    @test nested_check.vee_expectation > 0.8
    @test abs(nested_check.orbital_energy - parent_check.orbital_energy) < 0.01
    @test nested_check.vee_expectation > parent_check.vee_expectation
    payload = bond_aligned_diatomic_geometry_payload(nested_ops, source)
    slice = bond_aligned_diatomic_plane_slice(
        payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-5,
    )
    @test count(point -> point.group_kind == :shared_midpoint_slab, payload.points) == 0
    @test count(point -> point.group_kind == :residual_gaussian, payload.points) == nested_ops.residual_count
    @test slice.selected_count > 0
    @test parent_ops.residual_count == 0
end

@testset "Bond-aligned diatomic molecular supplement nested QW path" begin
    fixture = _bond_aligned_diatomic_nested_hybrid_qw_fixture(; bond_length = 1.4)
    if fixture === nothing
        @test true
    else
        (
            _basis,
            parent_ops,
            parent_check,
            _source,
            fixed_block,
            supplement,
            nested_ops,
            nested_check,
        ) = fixture

        @test supplement isa LegacyBondAlignedDiatomicGaussianSupplement
        @test nested_ops.basis === fixed_block
        @test nested_ops.gaussian_data === supplement
        @test nested_ops.interaction_treatment == :ggt_nearest
        @test nested_ops.residual_count > 0
        @test size(nested_ops.residual_centers, 1) == nested_ops.residual_count
        @test size(nested_ops.residual_widths, 1) == nested_ops.residual_count
        @test norm(nested_ops.overlap - I, Inf) < 1.0e-10
        @test nested_check.overlap_error < 1.0e-10
        @test all(isfinite, fixed_block.weights)
        @test minimum(fixed_block.weights) > 0.0
        @test nested_check.orbital_energy < -1.25
        @test 0.75 < nested_check.vee_expectation < 0.81
        @test abs(nested_check.orbital_energy - parent_check.orbital_energy) < 0.03
        @test abs(nested_check.vee_expectation - parent_check.vee_expectation) < 0.01
        @test parent_ops.residual_count == 0
    end
end

@testset "Bond-aligned diatomic geometry payloads and plane slices" begin
    basis, ordinary_ops, _check = _bond_aligned_diatomic_qw_fixture(; bond_length = 1.4)
    basis_payload = bond_aligned_diatomic_geometry_payload(basis)
    ordinary_payload = bond_aligned_diatomic_geometry_payload(ordinary_ops)
    basis_slice = bond_aligned_diatomic_plane_slice(
        basis_payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-12,
    )

    @test basis_payload isa GaussletBases.BondAlignedDiatomicGeometryPayload3D
    @test length(basis_payload.nuclei) == 2
    @test basis_payload.bond_axis == :z
    @test length(basis_payload.points) ==
        length(basis.basis_x) * length(basis.basis_y) * length(basis.basis_z)
    @test Set(point.group_kind for point in basis_payload.points) == Set([:gausslet_product])
    @test length(basis_payload.box_outlines) == 1
    @test basis_payload.box_outlines[1].group_kind == :basis_box
    @test ordinary_payload.bond_axis == :z
    @test length(ordinary_payload.points) == ordinary_ops.gausslet_count
    @test Set(point.group_kind for point in ordinary_payload.points) == Set([:gausslet_product])
    @test basis_slice.plane_axis == :y
    @test basis_slice.plane_value == 0.0
    @test basis_slice.plane_tol == 1.0e-12
    @test basis_slice.total_count == length(basis_payload.points)
    expected_y_count = count(y -> abs(y) <= 1.0e-12, centers(basis.basis_y))
    @test basis_slice.selected_count == length(basis.basis_x) * expected_y_count * length(basis.basis_z)
    @test all(abs(point.y) <= basis_slice.plane_tol for point in basis_slice.points)
    @test length(basis_slice.nuclei) == 2

    (
        _nested_basis,
        _parent_ops,
        _parent_check,
        _expansion,
        source,
        fixed_block,
        _nested_ops,
        _nested_check,
    ) = _bond_aligned_diatomic_nested_qw_fixture(; bond_length = 1.4)
    nested_payload = bond_aligned_diatomic_geometry_payload(fixed_block, source)
    expected_nested_groups =
        source.geometry.did_split ?
        Set([:shared_shell_layer, :left_child, :shared_midpoint_slab, :right_child]) :
        Set([:shared_shell_layer, :shared_child])
    expected_nested_box_count =
        2 + length(source.geometry.child_boxes) + (isnothing(source.geometry.shared_midpoint_box) ? 0 : 1)

    @test nested_payload isa GaussletBases.BondAlignedDiatomicGeometryPayload3D
    @test length(nested_payload.points) == size(fixed_block.fixed_centers, 1)
    @test Set(point.group_kind for point in nested_payload.points) == expected_nested_groups
    @test length(nested_payload.box_outlines) == expected_nested_box_count
    @test length(nested_payload.shell_provenance) == length(source.shared_shell_layers)
    @test nested_payload.shell_provenance[1].source_box == source.shared_shell_layers[1].provenance.source_box
    @test nested_payload.shell_provenance[1].next_inner_box ==
        source.shared_shell_layers[1].provenance.next_inner_box
    @test nested_payload.box_outlines[1].group_kind == :parent_box
    @test nested_payload.box_outlines[2].group_kind == :working_box
    @test count(box -> box.group_kind == :child_box, nested_payload.box_outlines) == length(source.geometry.child_boxes)
    @test any(box -> box.group_kind == :shared_midpoint_slab_box, nested_payload.box_outlines) ==
        !isnothing(source.geometry.shared_midpoint_box)

    hybrid_fixture = _bond_aligned_diatomic_hybrid_qw_fixture(; bond_length = 1.4)
    @test hybrid_fixture !== nothing
    if hybrid_fixture !== nothing
        (
            _hybrid_basis,
            _hybrid_parent_ops,
            _hybrid_parent_check,
            _supplement,
            hybrid_ops,
            _hybrid_check,
        ) = hybrid_fixture
        hybrid_payload = bond_aligned_diatomic_geometry_payload(hybrid_ops)
        @test count(point -> point.group_kind == :residual_gaussian, hybrid_payload.points) == hybrid_ops.residual_count
    end

    nested_hybrid_fixture = _bond_aligned_diatomic_nested_hybrid_qw_fixture(; bond_length = 1.4)
    @test nested_hybrid_fixture !== nothing
    if nested_hybrid_fixture !== nothing
        (
            _basis2,
            _parent_ops2,
            _parent_check2,
            hybrid_source,
            _hybrid_fixed_block,
            _hybrid_supplement,
            hybrid_nested_ops,
            _hybrid_nested_check,
        ) = nested_hybrid_fixture
        hybrid_nested_payload = bond_aligned_diatomic_geometry_payload(hybrid_nested_ops, hybrid_source)
        expected_hybrid_nested_groups =
            hybrid_source.geometry.did_split ?
            Set([:shared_shell_layer, :left_child, :shared_midpoint_slab, :right_child, :residual_gaussian]) :
            Set([:shared_shell_layer, :shared_child, :residual_gaussian])
        @test count(point -> point.group_kind == :residual_gaussian, hybrid_nested_payload.points) ==
            hybrid_nested_ops.residual_count
        @test Set(point.group_kind for point in hybrid_nested_payload.points) ==
            expected_hybrid_nested_groups
    end
end

@testset "Bond-aligned diatomic raw source geometry and 3d export" begin
    nested_hybrid_fixture = _bond_aligned_diatomic_nested_hybrid_qw_fixture(; bond_length = 1.4)
    @test nested_hybrid_fixture !== nothing
    if nested_hybrid_fixture !== nothing
        (
            _basis,
            _parent_ops,
            _parent_check,
            source,
            _fixed_block,
            _supplement,
            hybrid_nested_ops,
            _hybrid_nested_check,
        ) = nested_hybrid_fixture

        fixed_payload = bond_aligned_diatomic_geometry_payload(hybrid_nested_ops, source)
        source_payload = bond_aligned_diatomic_source_geometry_payload(source)
        expected_source_groups =
            source.geometry.did_split ?
            Set([:shared_shell_region, :left_child_region, :shared_midpoint_slab_region, :right_child_region]) :
            Set([:shared_shell_region, :shared_child_region])
        expected_child_box_count = length(source.geometry.child_boxes)
        expected_shared_midpoint_box = !isnothing(source.geometry.shared_midpoint_box)
        expected_shared_shell_points = sum(length(layer.support_indices) for layer in source.shared_shell_layers)
        expected_left_child_points = source.geometry.did_split ? length(source.child_sequences[1].support_indices) : 0
        expected_right_child_points = source.geometry.did_split ? length(source.child_sequences[2].support_indices) : 0
        expected_shared_child_points = source.geometry.did_split ? 0 : length(source.child_sequences[1].support_indices)
        expected_midpoint_points =
            isnothing(source.geometry.shared_midpoint_box) ?
            0 :
            prod(length.(source.geometry.shared_midpoint_box))
        first_shell = source.shared_shell_layers[1].provenance
        first_payload_shell = source_payload.shell_provenance[1]
        first_shell_source_dims = string(
            length(first_shell.source_box[1]), "x",
            length(first_shell.source_box[2]), "x",
            length(first_shell.source_box[3]),
        )
        first_shell_next_inner_dims = string(
            length(first_shell.next_inner_box[1]), "x",
            length(first_shell.next_inner_box[2]), "x",
            length(first_shell.next_inner_box[3]),
        )

        @test all(isnothing(sequence.support_states) for sequence in source.child_sequences)
        @test Set(point.group_kind for point in source_payload.points) == expected_source_groups
        @test length(source_payload.points) == prod(length.(source.geometry.parent_box))
        @test length(source_payload.shell_provenance) == length(source.shared_shell_layers)
        @test first_shell.source_box == source.geometry.parent_box
        @test first_shell.source_point_count ==
            prod(length.(first_shell.source_box)) - prod(length.(first_shell.next_inner_box))
        @test first_shell.retained_fixed_count == length(source.sequence.layer_column_ranges[1])
        @test first_payload_shell.source_box == first_shell.source_box
        @test first_payload_shell.next_inner_box == first_shell.next_inner_box
        @test first_payload_shell.source_point_count == first_shell.source_point_count
        @test first_payload_shell.retained_fixed_count == first_shell.retained_fixed_count
        @test count(point -> point.group_kind == :shared_shell_region, source_payload.points) ==
            expected_shared_shell_points
        @test count(point -> point.group_kind == :left_child_region, source_payload.points) ==
            expected_left_child_points
        @test count(point -> point.group_kind == :right_child_region, source_payload.points) ==
            expected_right_child_points
        @test count(point -> point.group_kind == :shared_child_region, source_payload.points) ==
            expected_shared_child_points
        @test count(point -> point.group_kind == :shared_midpoint_slab_region, source_payload.points) ==
            expected_midpoint_points
        @test any(box.group_kind == :parent_box for box in source_payload.box_outlines)
        @test any(box.group_kind == :working_box for box in source_payload.box_outlines)
        @test count(box -> box.group_kind == :child_box, source_payload.box_outlines) == expected_child_box_count
        @test any(box -> box.group_kind == :shared_midpoint_slab_box, source_payload.box_outlines) ==
            expected_shared_midpoint_box

        fixed_slice = bond_aligned_diatomic_plane_slice(
            fixed_payload;
            plane_axis = :y,
            plane_value = 0.0,
            plane_tol = 1.0e-12,
        )
        debug_slice = bond_aligned_diatomic_plane_slice(
            fixed_payload;
            plane_axis = :y,
            plane_value = 0.0,
            plane_tol = 1.0e-5,
        )
        source_slice = bond_aligned_diatomic_plane_slice(
            source_payload;
            plane_axis = :y,
            plane_value = 0.0,
            plane_tol = 1.0e-5,
        )
        @test debug_slice.selected_count >= fixed_slice.selected_count
        @test source_slice.selected_count > debug_slice.selected_count
        @test source_slice.shell_provenance == source_payload.shell_provenance

    end
end
