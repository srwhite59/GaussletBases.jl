using LinearAlgebra
using Test

using GaussletBases

function _synthetic_external_gto_probes()
    orbitals = CartesianGaussianShellOrbitalRepresentation3D[
        CartesianGaussianShellOrbitalRepresentation3D(
            "s_core",
            (0, 0, 0),
            (0.0, 0.0, 0.0),
            [0.7],
            [1.0],
            :axiswise_normalized_cartesian_gaussian,
        ),
        CartesianGaussianShellOrbitalRepresentation3D(
            "pz_valence",
            (0, 0, 1),
            (0.0, 0.0, 0.0),
            [0.45],
            [1.0],
            :axiswise_normalized_cartesian_gaussian,
        ),
        CartesianGaussianShellOrbitalRepresentation3D(
            "s_shifted",
            (0, 0, 0),
            (0.35, 0.0, 0.0),
            [0.9],
            [1.0],
            :axiswise_normalized_cartesian_gaussian,
        ),
    ]
    return CartesianGaussianShellSupplementRepresentation3D(
        :external_gto_import_test,
        orbitals,
        (; source = :cartesian_external_gto_import_test),
    )
end

function _synthetic_working_basis()
    return build_basis(MappedUniformBasisSpec(:G10;
        count = 5,
        mapping = fit_asinh_mapping_for_strength(s = 0.5, npoints = 5, xmax = 4.0),
        reference_spacing = 1.0,
    ))
end

function _orthonormal_source_coefficients(S_GG)
    return inv(cholesky(Symmetric(S_GG)).U)
end

function _source_packet(; alpha_spin = :restricted, beta = nothing, fingerprint = nothing)
    probes = _synthetic_external_gto_probes()
    S_GG = Matrix{Float64}(GaussletBases._cartesian_supplement_cross_overlap(probes, probes))
    C = _orthonormal_source_coefficients(S_GG)
    alpha = ExternalGTOOrbitalSpinBlock(alpha_spin, C[:, 1:2], [1.0, 1.0])
    return ExternalGTOOrbitalPacket(
        probes,
        S_GG,
        alpha;
        beta = beta,
        S_GG_fingerprint = something(fingerprint, external_gto_overlap_fingerprint(S_GG)),
        provenance = (; source = :synthetic_test_packet),
    )
end

@testset "External GTO orbital import" begin
    working = _synthetic_working_basis()
    packet = _source_packet()
    result = import_external_gto_orbitals(working, packet)
    S_FG = gto_overlap_matrix(working, packet.probes)

    @test result.cross_overlap_size == (length(working)^3, length(packet.probes.orbitals))
    @test result.cross_overlap_finite
    @test result.ordering_fingerprint_valid
    @test result.S_GG_fingerprint_valid
    @test result.S_GG_symmetry_error < 1.0e-12
    @test result.S_GG_expected_error < 1.0e-12
    @test size(result.alpha.imported_coefficients) == (length(working)^3, 2)
    @test result.alpha.imported_coefficients ≈ S_FG * packet.alpha.coefficients
    @test result.alpha.source_orthogonality_error < 1.0e-12
    @test result.alpha.capture_matrix ≈
        transpose(result.alpha.imported_coefficients) * result.alpha.imported_coefficients
    @test result.alpha.density_trace_source == 2.0
    @test result.alpha.density_trace_capture ≈
        dot(packet.alpha.occupations, result.alpha.orbital_captures)
    @test result.alpha.density_trace_loss ≈
        result.alpha.density_trace_source - result.alpha.density_trace_capture
    @test result.alpha.worst_orbital_capture == minimum(result.alpha.orbital_captures)

    theta = 0.37
    rotation = [cos(theta) -sin(theta); sin(theta) cos(theta)]
    rotated_alpha = ExternalGTOOrbitalSpinBlock(
        :restricted,
        packet.alpha.coefficients * rotation,
        [1.0, 1.0],
    )
    rotated_packet = ExternalGTOOrbitalPacket(
        packet.probes,
        packet.S_GG,
        rotated_alpha;
        provenance = (; source = :rotated_synthetic_test_packet),
    )
    rotated_result = import_external_gto_orbitals(working, rotated_packet)
    @test rotated_result.alpha.source_orthogonality_error < 1.0e-12
    @test rotated_result.alpha.density_trace_capture ≈
        result.alpha.density_trace_capture atol = 1.0e-12 rtol = 1.0e-12

    C = _orthonormal_source_coefficients(packet.S_GG)
    beta = ExternalGTOOrbitalSpinBlock(:beta, C[:, 2:2], [1.0])
    spin_packet = _source_packet(alpha_spin = :alpha, beta = beta)
    spin_result = import_external_gto_orbitals(working, spin_packet)
    @test spin_result.alpha.spin == :alpha
    @test spin_result.beta !== nothing
    @test spin_result.beta.spin == :beta
    @test size(spin_result.beta.imported_coefficients, 2) == 1
    @test spin_result.beta.source_orthogonality_error < 1.0e-12

    bad_fingerprint_packet = _source_packet(fingerprint = "not-the-source-overlap")
    @test_throws ArgumentError import_external_gto_orbitals(working, bad_fingerprint_packet)

    wrong_metric = copy(packet.S_GG)
    wrong_metric[1, 1] += 0.05
    wrong_metric_alpha = ExternalGTOOrbitalSpinBlock(
        :restricted,
        _orthonormal_source_coefficients(wrong_metric)[:, 1:2],
        [1.0, 1.0],
    )
    wrong_metric_packet = ExternalGTOOrbitalPacket(
        packet.probes,
        wrong_metric,
        wrong_metric_alpha;
        S_GG_fingerprint = external_gto_overlap_fingerprint(wrong_metric),
        provenance = (; source = :wrong_but_refingerprinted_metric_test_packet),
    )
    @test_throws ArgumentError import_external_gto_orbitals(working, wrong_metric_packet)

    nonsymmetric_metric = copy(packet.S_GG)
    nonsymmetric_metric[1, 2] += 1.0e-3
    nonsymmetric_packet = ExternalGTOOrbitalPacket(
        packet.probes,
        nonsymmetric_metric,
        packet.alpha;
        S_GG_fingerprint = external_gto_overlap_fingerprint(nonsymmetric_metric),
        provenance = (; source = :nonsymmetric_metric_test_packet),
    )
    @test_throws ArgumentError import_external_gto_orbitals(working, nonsymmetric_packet)

    bad_order_packet = ExternalGTOOrbitalPacket(
        packet.probes,
        packet.S_GG,
        packet.alpha;
        ordering_fingerprint = "not-the-ao-order",
        provenance = (; source = :bad_ordering_test_packet),
    )
    @test_throws ArgumentError import_external_gto_orbitals(working, bad_order_packet)

    @test_throws ArgumentError _source_packet(alpha_spin = :beta)
    @test_throws ArgumentError _source_packet(alpha_spin = :restricted, beta = beta)
    bad_beta = ExternalGTOOrbitalSpinBlock(:alpha, C[:, 2:2], [1.0])
    @test_throws ArgumentError _source_packet(alpha_spin = :alpha, beta = bad_beta)

    bad_alpha = ExternalGTOOrbitalSpinBlock(
        :restricted,
        packet.alpha.coefficients .* 2.0,
        [1.0, 1.0],
    )
    bad_orthogonality_packet = ExternalGTOOrbitalPacket(
        packet.probes,
        packet.S_GG,
        bad_alpha;
        provenance = (; source = :bad_orthogonality_test_packet),
    )
    @test_throws ArgumentError import_external_gto_orbitals(working, bad_orthogonality_packet)
end
