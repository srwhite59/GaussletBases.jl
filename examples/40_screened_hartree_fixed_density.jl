using GaussletBases
using LinearAlgebra

orbitals = [
    CartesianGaussianShellOrbitalRepresentation3D(
        "1s", (0, 0, 0), (0.0, 0.0, 0.0), [0.8], [1.0],
        :axiswise_normalized_cartesian_gaussian),
    CartesianGaussianShellOrbitalRepresentation3D(
        "2px", (1, 0, 0), (0.0, 0.0, 0.0), [0.35], [1.0],
        :axiswise_normalized_cartesian_gaussian),
]
expansion = coulomb_gaussian_expansion(doacc = true)
eri = gaussian_coulomb_pair_matrix(orbitals; expansion)
pair(i, j) = gaussian_coulomb_pair_index(i, j, length(orbitals))

reference_coefficients = [1.0 1.0; 1.0 -1.0] ./ sqrt(2)
occupations = [1.25, 0.75]
P0 = reference_coefficients * Diagonal(occupations) *
    transpose(reference_coefficients)
J0 = [sum(P0[k, l] * eri[pair(i, j), pair(k, l)]
    for k in 1:2, l in 1:2) for i in 1:2, j in 1:2]
reference_coulomb_self_integral = sum(P0 .* J0)
V_IDA = [eri[pair(i, i), pair(j, j)] for i in 1:2, j in 1:2]

ham = CartesianIDAHamiltonian(
    [0.8 0.05; 0.05 1.1],
    [[-0.9 -0.02; -0.02 -0.6]],
    V_IDA,
    1,
    1;
    nuclear_charges = [2.0],
    nuclear_positions = [(0.0, 0.0, 0.0)],
)
physical_h1 = one_body_hamiltonian(ham)
field = ExactRepresentedHartreeField(
    J0,
    reference_coulomb_self_integral;
    provenance = "two-orbital analytic Gaussian Coulomb reference",
)
correction = screened_hartree_correction(
    ham.electron_electron_ida,
    field,
    reference_coefficients,
    occupations,
)

q0 = diag(P0)
delta = screened_hartree_delta_one_body(correction)
corrected_direct = 0.5 * dot(q0, V_IDA * q0) + sum(P0 .* delta) +
    screened_hartree_energy_constant(correction)
@assert corrected_direct ≈ 0.5 * reference_coulomb_self_integral atol = 1.0e-11
@assert (Diagonal(V_IDA * q0) + delta) * reference_coefficients ≈
    J0 * reference_coefficients atol = 1.0e-11
@assert one_body_hamiltonian(ham) == physical_h1

println("field kind: ", screened_hartree_field_kind(correction))
println("reference charge: ", sum(occupations))
println("screened energy constant: ", screened_hartree_energy_constant(correction))
println("consistency error: ", screened_hartree_consistency_error(correction))
