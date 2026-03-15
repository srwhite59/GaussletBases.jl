using LinearAlgebra
using GaussletBases

ub = build_basis(UniformBasisSpec(:G10; xmin = -2.0, xmax = 2.0, spacing = 1.0))
rep = basis_representation(ub)
hierarchy = refine_partition(hierarchical_partition(rep, [-2.5, -0.5, 0.5, 2.5]), 1)

base_pgdg = build_leaf_pgdg(hierarchy; primitives_per_leaf = 2, width_scale = 0.75)
augmented_pgdg = augment_leaf_pgdg(
    base_pgdg;
    by_leaf = Dict(
        4 => [LeafGaussianSpec1D(relative_position = 0.5, width_scale = 0.2)],
    ),
)

base_rep = basis_representation(base_pgdg)
augmented_rep = basis_representation(augmented_pgdg)

println(base_pgdg)
println("base leaf 4 primitives: ", collect(leaf_primitive_indices(base_pgdg, 4)))
println("base leaf 5 primitives: ", collect(leaf_primitive_indices(base_pgdg, 5)))

println(augmented_pgdg)
println("augmented leaf 4 primitives: ", collect(leaf_primitive_indices(augmented_pgdg, 4)))
println("augmented leaf 5 primitives: ", collect(leaf_primitive_indices(augmented_pgdg, 5)))
println("primitive origins: ", primitive_origins(augmented_pgdg))
println("primitive leaf boxes: ", primitive_leaf_boxes(augmented_pgdg))
println("augmented overlap size: ", size(augmented_rep.basis_matrices.overlap))
println("augmented kinetic inf-norm: ", norm(augmented_rep.basis_matrices.kinetic, Inf))
println("augmented position inf-norm: ", norm(augmented_rep.basis_matrices.position, Inf))
