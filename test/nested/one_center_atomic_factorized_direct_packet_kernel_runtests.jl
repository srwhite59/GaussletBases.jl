# Integration/slow test. Do not include in default nested runner.

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
