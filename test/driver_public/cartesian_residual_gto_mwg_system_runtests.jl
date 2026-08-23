using GaussletBases
using LinearAlgebra
using Test

@testset "Public residual-GTO/MWG system" begin
    nuclei = NTuple{3,Float64}[(0.0, 0.0, -2.0), (0.0, 0.0, 2.0)]
    system = (;
        atom_symbols = ["H", "H"], nuclear_charges = [1.0, 1.0],
        atom_locations = nuclei, nup = 1, ndn = 1)
    basis = (;
        q = 5, core_spacing = 0.5,
        xmax_parallel = 6.0, xmax_transverse = 4.0)
    supplement = (;
        basis_by_center = ["cc-pVTZ", "cc-pVTZ"], lmax = 1,
        uncontracted = false, width_filtering = nothing)

    result = cartesian_residual_gto_mwg_system(
        system; basis, supplement)
    ham = result.hamiltonian
    nfinal = 487 + 18
    matrices = [ham.kinetic, ham.electron_electron_ida,
        ham.nuclear_attraction_unit_by_center...]
    @test ham isa CartesianIDAHamiltonian{Float64}
    @test size(ham.kinetic) == size(ham.electron_electron_ida) == (nfinal, nfinal)
    @test all(matrix -> all(isfinite, matrix), matrices)
    @test all(matrix -> norm(matrix - transpose(matrix), Inf) <= 1.0e-10, matrices)
    @test ham.nup == ham.ndn == 1
    @test ham.nuclear_charges == [1.0, 1.0]
    @test ham.nuclear_positions == [nucleus[axis] for nucleus in nuclei, axis in 1:3]

    h1 = one_body_hamiltonian(ham)
    orbital = eigen(Symmetric(h1)).vectors[:, 1]
    density = abs2.(orbital)
    @test dot(density, ham.electron_electron_ida * density) ≈
        0.4574161883692301 atol = 1.0e-10

    probe_source = legacy_bond_aligned_diatomic_gaussian_supplement(
        "H", "cc-pVTZ", nuclei; lmax = 1)
    source_representation = basis_representation(probe_source)
    probe_orbitals = source_representation.orbitals[[4, 5]]
    @test getproperty.(probe_orbitals, :angular_powers) == [(1, 0, 0), (0, 1, 0)]
    @test all(orbital -> orbital.center == first(nuclei), probe_orbitals)
    @test probe_orbitals[1].exponents == probe_orbitals[2].exponents
    @test all(orbital -> length(orbital.exponents) == 1 &&
        orbital.coefficients == [1.0] &&
        orbital.primitive_normalization == :axiswise_normalized_cartesian_gaussian,
        probe_orbitals)
    probes = CartesianGaussianShellSupplementRepresentation3D(
        :public_ci_probe, probe_orbitals,
        basis_metadata(source_representation))
    overlap = gto_overlap_matrix(result, probes)
    rows = [nfinal, 1, 488, 1]
    @test size(overlap) == (nfinal, 2)
    @test all(isfinite, overlap)
    @test gto_overlap_matrix(result, probes; block_indices = rows) == overlap[rows, :]
    @test_throws BoundsError gto_overlap_matrix(result, probes; block_indices = [0])

    S_GG = Matrix{Float64}(I, 2, 2)
    coefficients = Matrix{Float64}(I, 2, 2)
    source_block = ExternalGTOOrbitalSpinBlock(:restricted, coefficients, [1.0, 1.0])
    packet = ExternalGTOOrbitalPacket(probes, S_GG, source_block)
    imported = import_external_gto_orbitals(result, packet)
    @test imported.cross_overlap_size == size(overlap)
    @test imported.cross_overlap_finite
    @test imported.ordering_fingerprint_valid
    @test imported.S_GG_fingerprint_valid
    @test imported.S_GG_symmetry_error <= 1.0e-12
    @test imported.S_GG_expected_error <= 1.0e-12
    @test imported.alpha.imported_coefficients == overlap
    @test imported.alpha.source_orthogonality_error <= 1.0e-12
    @test imported.alpha.capture_matrix ≈ transpose(overlap) * overlap
    @test imported.alpha.density_trace_source == 2.0
    @test imported.alpha.density_trace_capture ≈
        dot(source_block.occupations, imported.alpha.orbital_captures)
    @test abs(imported.alpha.worst_orbital_capture - 1.0) <= 1.0e-8

    theta = 0.37
    rotation = [cos(theta) -sin(theta); sin(theta) cos(theta)]
    rotated_block = ExternalGTOOrbitalSpinBlock(
        :restricted, coefficients * rotation, [1.0, 1.0])
    rotated_packet = ExternalGTOOrbitalPacket(probes, S_GG, rotated_block)
    rotated = import_external_gto_orbitals(result, rotated_packet)
    @test rotated.alpha.imported_coefficients ≈ overlap * rotation
    @test rotated.alpha.source_orthogonality_error <= 1.0e-12
    @test rotated.alpha.density_trace_capture ≈
        imported.alpha.density_trace_capture atol = 1.0e-12 rtol = 1.0e-12

    alpha = ExternalGTOOrbitalSpinBlock(:alpha, coefficients[:, 1:1], [1.0])
    beta = ExternalGTOOrbitalSpinBlock(:beta, coefficients[:, 2:2], [1.0])
    spin_packet = ExternalGTOOrbitalPacket(probes, S_GG, alpha; beta)
    spin_import = import_external_gto_orbitals(result, spin_packet)
    @test spin_import.alpha.spin == :alpha
    @test spin_import.beta !== nothing
    @test spin_import.beta.spin == :beta
    @test spin_import.alpha.imported_coefficients == overlap[:, 1:1]
    @test spin_import.beta.imported_coefficients == overlap[:, 2:2]

    stale_order = ExternalGTOOrbitalPacket(
        probes, S_GG, source_block; ordering_fingerprint = "stale")
    @test_throws ArgumentError import_external_gto_orbitals(result, stale_order)
    stale_overlap = ExternalGTOOrbitalPacket(
        probes, S_GG, source_block; S_GG_fingerprint = "stale")
    @test_throws ArgumentError import_external_gto_orbitals(result, stale_overlap)

    wrong_metric = copy(S_GG)
    wrong_metric[1, 1] += 0.05
    wrong_metric_packet = ExternalGTOOrbitalPacket(
        probes, wrong_metric, source_block;
        S_GG_fingerprint = external_gto_overlap_fingerprint(wrong_metric))
    @test_throws ArgumentError import_external_gto_orbitals(result, wrong_metric_packet)
    nonsymmetric_metric = copy(S_GG)
    nonsymmetric_metric[1, 2] += 1.0e-3
    nonsymmetric_packet = ExternalGTOOrbitalPacket(
        probes, nonsymmetric_metric, source_block;
        S_GG_fingerprint = external_gto_overlap_fingerprint(nonsymmetric_metric))
    @test_throws ArgumentError import_external_gto_orbitals(result, nonsymmetric_packet)

    beta_only = ExternalGTOOrbitalSpinBlock(:beta, coefficients[:, 2:2], [1.0])
    @test_throws ArgumentError ExternalGTOOrbitalPacket(probes, S_GG, beta_only)
    @test_throws ArgumentError ExternalGTOOrbitalPacket(
        probes, S_GG, source_block; beta)
    wrong_beta = ExternalGTOOrbitalSpinBlock(:alpha, coefficients[:, 2:2], [1.0])
    @test_throws ArgumentError ExternalGTOOrbitalPacket(probes, S_GG, alpha; beta = wrong_beta)

    nonorthogonal = ExternalGTOOrbitalSpinBlock(
        :restricted, 2.0 .* coefficients, [1.0, 1.0])
    nonorthogonal_packet = ExternalGTOOrbitalPacket(probes, S_GG, nonorthogonal)
    @test_throws ArgumentError import_external_gto_orbitals(result, nonorthogonal_packet)

    @test_throws ArgumentError cartesian_residual_gto_mwg_system(
        merge(system, (; extra = true)); basis, supplement)
    @test_throws ArgumentError cartesian_residual_gto_mwg_system(
        system; basis = merge(basis, (; radius = 4.0)), supplement)
    @test_throws ArgumentError cartesian_residual_gto_mwg_system(
        system; basis, supplement = merge(supplement, (; extra = true)))
end
