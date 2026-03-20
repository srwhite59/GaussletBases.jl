using LinearAlgebra
using GaussletBases

Z = 2.0
basis = build_basis(MappedUniformBasisSpec(:G10;
    count = 9,
    mapping = IdentityMapping(),
    reference_spacing = 0.5,
))

ops = ordinary_cartesian_ida_operators(
    basis;
    expansion = coulomb_gaussian_expansion(doacc = false),
    Z = Z,
    backend = :numerical_reference,
)

decomposition = eigen(Hermitian(ops.one_body_hamiltonian))
orbital_energy = decomposition.values[1]
orbital = decomposition.vectors[:, 1]
vee_expectation = ordinary_cartesian_vee_expectation(ops, orbital)
reference_value = (5.0 / 8.0) * Z

println("Ordinary Cartesian IDA 1s^2 Vee check")
println("  basis: ", basis)
println("  operators: ", ops)
println("  orbital dimension: ", length(orbitals(ops)))
println("  3D overlap error: ", norm(ops.overlap_3d - I, Inf))
println("  lowest one-body orbital energy: ", orbital_energy)
println("  doubly occupied 1s-state IDA Vee expectation: ", vee_expectation)
println("  hydrogenic reference (5/8)Z: ", reference_value)
println("  difference: ", vee_expectation - reference_value)
