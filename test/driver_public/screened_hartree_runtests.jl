using GaussletBases
using LinearAlgebra
using Test

@testset "Public supplied-field screened Hartree" begin
    orbitals = [
        CartesianGaussianShellOrbitalRepresentation3D(
            "1s", (0, 0, 0), (0.0, 0.0, 0.0), [0.8], [1.0],
            :axiswise_normalized_cartesian_gaussian),
        CartesianGaussianShellOrbitalRepresentation3D(
            "2px", (1, 0, 0), (0.0, 0.0, 0.0), [0.35], [1.0],
            :axiswise_normalized_cartesian_gaussian),
    ]
    eri = gaussian_coulomb_pair_matrix(
        orbitals; expansion = coulomb_gaussian_expansion(doacc = true))
    pair(i, j) = gaussian_coulomb_pair_index(i, j, 2)
    C = [1.0 1.0; 1.0 -1.0] ./ sqrt(2)
    occupations = [1.25, 0.75]
    P0 = C * Diagonal(occupations) * transpose(C)
    J0 = [sum(P0[k, l] * eri[pair(i, j), pair(k, l)]
        for k in 1:2, l in 1:2) for i in 1:2, j in 1:2]
    E0 = sum(P0 .* J0)
    V = [eri[pair(i, i), pair(j, j)] for i in 1:2, j in 1:2]
    ham = CartesianIDAHamiltonian([0.8 0.05; 0.05 1.1],
        [[-0.9 -0.02; -0.02 -0.6]], V, 1, 1;
        nuclear_charges = [2.0], nuclear_positions = [(0.0, 0.0, 0.0)])
    h1 = one_body_hamiltonian(ham)

    supplied = copy(J0)
    exact = ExactRepresentedHartreeField(supplied, E0;
        provenance = "bounded two-orbital Gaussian oracle")
    supplied[1, 1] += 1.0
    correction = screened_hartree_correction(
        ham.electron_electron_ida, exact, C, occupations)
    delta = screened_hartree_delta_one_body(correction)
    q0 = diag(P0)

    @test exact.matrix == J0
    @test correction isa ScreenedHartreeCorrection
    @test screened_hartree_field_kind(correction) === :exact_represented
    @test size(delta) == (2, 2)
    @test norm(delta - transpose(delta), Inf) <= 1.0e-10
    @test norm(delta, Inf) > 1.0e-8
    @test sum(q0) ≈ sum(occupations) atol = 1.0e-12
    @test screened_hartree_consistency_error(correction) ≈ 0.0 atol = 1.0e-11
    @test 0.5 * dot(q0, V * q0) + sum(P0 .* delta) +
        screened_hartree_energy_constant(correction) ≈ 0.5 * E0 atol = 1.0e-11
    @test (Diagonal(V * q0) + delta) * C ≈ J0 * C atol = 1.0e-11
    @test one_body_hamiltonian(ham) == h1
    delta[1, 1] += 1.0
    @test screened_hartree_delta_one_body(correction)[1, 1] != delta[1, 1]

    fitted_error = -2.0e-4
    fitted = FittedReferenceHartreeField(J0, E0 - fitted_error;
        provenance = "synthetic fitted-field control")
    fitted_correction = screened_hartree_correction(V, fitted, C, occupations)
    @test screened_hartree_field_kind(fitted_correction) === :fitted_reference
    @test screened_hartree_consistency_error(fitted_correction) ≈
        fitted_error atol = 1.0e-12

    @test_throws ArgumentError ExactRepresentedHartreeField(
        J0, E0; provenance = "")
    @test_throws ArgumentError ExactRepresentedHartreeField(
        [J0[1, 1] J0[1, 2] + 1.0e-5; J0[2, 1] J0[2, 2]], E0;
        provenance = "asymmetric")
    bad_finite = copy(J0)
    bad_finite[1, 1] = NaN
    @test_throws ArgumentError ExactRepresentedHartreeField(
        bad_finite, E0; provenance = "nonfinite")
    @test_throws DimensionMismatch screened_hartree_correction(
        ones(3, 3), exact, C, occupations)
    @test_throws ArgumentError screened_hartree_correction(
        [V[1, 1] V[1, 2] + 1.0e-5; V[2, 1] V[2, 2]], exact, C, occupations)
    @test_throws ArgumentError screened_hartree_correction(
        V, exact, C, [1.0, -0.1])
    @test_throws ArgumentError screened_hartree_correction(
        V, exact, 0.9 .* C, occupations)
    inconsistent = ExactRepresentedHartreeField(J0, E0 + 1.0e-4;
        provenance = "inconsistent exact control")
    @test_throws ArgumentError screened_hartree_correction(
        V, inconsistent, C, occupations)
end
