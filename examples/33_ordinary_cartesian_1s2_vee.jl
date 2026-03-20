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

check = GaussletBases.ordinary_cartesian_1s2_check(ops)
reference_value = (5.0 / 8.0) * Z

println("Ordinary Cartesian IDA 1s^2 Vee check")
println("  basis: ", basis)
println("  operators: ", ops)
println("  orbital dimension: ", length(orbitals(ops)))
println("  3D overlap error: ", check.overlap_error)
println("  lowest one-body orbital energy: ", check.orbital_energy)
println("  doubly occupied 1s-state IDA Vee expectation: ", check.vee_expectation)
println("  hydrogenic reference (5/8)Z: ", reference_value)
println("  difference: ", check.vee_expectation - reference_value)
