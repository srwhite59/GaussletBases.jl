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
    probes = CartesianGaussianShellSupplementRepresentation3D(
        :public_ci_probe, [source_representation.orbitals[2]],
        basis_metadata(source_representation))
    overlap = gto_overlap_matrix(result, probes)
    rows = [nfinal, 1, 488, 1]
    @test size(overlap) == (nfinal, 1)
    @test all(isfinite, overlap)
    @test gto_overlap_matrix(result, probes; block_indices = rows) == overlap[rows, :]
    @test_throws BoundsError gto_overlap_matrix(result, probes; block_indices = [0])

    source_block = ExternalGTOOrbitalSpinBlock(:restricted, ones(1, 1), [2.0])
    packet = ExternalGTOOrbitalPacket(probes, ones(1, 1), source_block)
    imported = import_external_gto_orbitals(result, packet)
    @test imported.cross_overlap_size == size(overlap)
    @test imported.cross_overlap_finite
    @test imported.alpha.imported_coefficients == overlap
    @test imported.alpha.source_orthogonality_error <= 1.0e-12
    @test abs(imported.alpha.worst_orbital_capture - 1.0) <= 1.0e-8
    stale_packet = ExternalGTOOrbitalPacket(
        probes, ones(1, 1), source_block; ordering_fingerprint = "stale")
    @test_throws ArgumentError import_external_gto_orbitals(result, stale_packet)

    @test_throws ArgumentError cartesian_residual_gto_mwg_system(
        merge(system, (; extra = true)); basis, supplement)
    @test_throws ArgumentError cartesian_residual_gto_mwg_system(
        system; basis = merge(basis, (; radius = 4.0)), supplement)
    @test_throws ArgumentError cartesian_residual_gto_mwg_system(
        system; basis, supplement = merge(supplement, (; extra = true)))
end
