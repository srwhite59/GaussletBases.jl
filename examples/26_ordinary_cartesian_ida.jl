using LinearAlgebra
using GaussletBases

Z = 2.0
npoints = 5
rmax = 6.0

mapping = fit_asinh_mapping_for_strength(s = 0.5, npoints = npoints, xmax = rmax)
basis = build_basis(MappedUniformBasisSpec(:G10;
    count = npoints,
    mapping = mapping,
    reference_spacing = 1.0,
))

ops = ordinary_cartesian_ida_operators(
    basis;
    expansion = coulomb_gaussian_expansion(doacc = false),
    Z = Z,
    backend = :pgdg_experimental,
)

println("Ordinary Cartesian He-style IDA ingredients")
println("  basis: ", basis)
println("  mapping: ", mapping)
println("  operators: ", ops)
println("  number of product orbitals: ", length(orbitals(ops)))
println("  one-body matrix size: ", size(ops.one_body_hamiltonian))
println("  interaction matrix size: ", size(ops.interaction_matrix))
println("  3D overlap infinity-norm error: ", norm(ops.overlap_3d - I, Inf))
println("  ||H1||_inf: ", opnorm(ops.one_body_hamiltonian, Inf))
println("  ||Vee||_inf: ", opnorm(ops.interaction_matrix, Inf))
println("  first orbital: ", orbitals(ops)[1])
println("  last orbital: ", orbitals(ops)[end])
