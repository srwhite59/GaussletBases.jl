using LinearAlgebra
using GaussletBases

Z = 1.0
basis = build_basis(UniformBasisSpec(:G10; xmin = -2.0, xmax = 2.0, spacing = 1.0))
representation = basis_representation(basis; operators = (:overlap, :kinetic))
overlap_1d = representation.basis_matrices.overlap
kinetic_1d = representation.basis_matrices.kinetic
expansion = coulomb_gaussian_expansion()

gaussian_factors = [
    gaussian_factor_matrix(basis; exponent = exponent, center = 0.0)
    for exponent in expansion.exponents
]

overlap_3d = kron(overlap_1d, kron(overlap_1d, overlap_1d))
kinetic_3d =
    kron(kinetic_1d, kron(overlap_1d, overlap_1d)) +
    kron(overlap_1d, kron(kinetic_1d, overlap_1d)) +
    kron(overlap_1d, kron(overlap_1d, kinetic_1d))

nuclear_3d = zeros(Float64, size(overlap_3d))
for term in eachindex(expansion.coefficients)
    factor = gaussian_factors[term]
    nuclear_3d .-= Z * expansion.coefficients[term] .* kron(factor, kron(factor, factor))
end

hamiltonian = kinetic_3d + nuclear_3d
eigenvalues = eigen(Hermitian(hamiltonian)).values
ground_energy = minimum(eigenvalues)

println("Ordinary Cartesian hydrogen through Coulomb expansion")
println("  1D basis functions: ", length(basis))
println("  3D product dimension: ", size(hamiltonian, 1))
println("  Coulomb Gaussian terms: ", length(expansion))
println("  1D overlap error: ", norm(overlap_1d - I, Inf))
println("  3D overlap error: ", norm(overlap_3d - I, Inf))
println("  ground-state energy: ", ground_energy)
println("  exact hydrogen energy: ", -0.5)
