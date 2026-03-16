using LinearAlgebra
using GaussletBases

function orthonormalize_1d(overlap, operators)
    decomposition = eigen(Symmetric(overlap))
    invhalf =
        decomposition.vectors *
        Diagonal(1.0 ./ sqrt.(decomposition.values)) *
        transpose(decomposition.vectors)
    return invhalf * overlap * invhalf, [invhalf * operator * invhalf for operator in operators]
end

Z = 1.0
npoints = 5
rmax = 6.0

mapping = fit_asinh_mapping_for_extent(npoints = npoints, xmax = rmax)
basis = build_basis(MappedUniformBasisSpec(:G10;
    count = npoints,
    mapping = mapping,
    reference_spacing = 1.0,
))

representation = basis_representation(basis; operators = (:overlap, :kinetic))
overlap_1d = representation.basis_matrices.overlap
kinetic_1d = representation.basis_matrices.kinetic

expansion = coulomb_gaussian_expansion(doacc = false)
gaussian_factors = gaussian_factor_matrices(
    basis;
    exponents = expansion.exponents,
    center = 0.0,
)

_, transformed = orthonormalize_1d(overlap_1d, [kinetic_1d, gaussian_factors...])
kinetic_orth = first(transformed)
gaussian_orth = transformed[2:end]

dimension = npoints^3
identity_1d = Matrix{Float64}(I, npoints, npoints)
hamiltonian =
    kron(kinetic_orth, kron(identity_1d, identity_1d)) +
    kron(identity_1d, kron(kinetic_orth, identity_1d)) +
    kron(identity_1d, kron(identity_1d, kinetic_orth))

for term in eachindex(expansion.coefficients)
    factor = gaussian_orth[term]
    hamiltonian .-= Z * expansion.coefficients[term] .* kron(factor, kron(factor, factor))
end

ground_energy = minimum(eigen(Hermitian(hamiltonian)).values)

println("Mapped Cartesian hydrogen through Coulomb expansion")
println("  basis: ", basis)
println("  mapping: ", mapping)
println("  outer centers: ", first(centers(basis)), "  ", last(centers(basis)))
println("  1D overlap error: ", norm(overlap_1d - I, Inf))
println("  3D product dimension: ", dimension)
println("  Coulomb Gaussian terms: ", length(expansion))
println("  ground-state energy: ", ground_energy)
println("  exact hydrogen energy: ", -0.5)
