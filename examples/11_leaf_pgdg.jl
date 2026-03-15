using LinearAlgebra
using GaussletBases

ub = build_basis(UniformBasisSpec(:G10; xmin = -2.0, xmax = 2.0, spacing = 1.0))
rep = basis_representation(ub)

coarse_hierarchy = hierarchical_partition(rep, [-2.5, -0.5, 0.5, 2.5])
refined_hierarchy = refine_partition(coarse_hierarchy, 1)

coarse_pgdg = build_leaf_pgdg(coarse_hierarchy; primitives_per_leaf = 2, width_scale = 0.75)
refined_pgdg = build_leaf_pgdg(refined_hierarchy; primitives_per_leaf = 2, width_scale = 0.75)

coarse_rep = basis_representation(coarse_pgdg)
refined_rep = basis_representation(refined_pgdg)

println(coarse_pgdg)
for leaf in leaf_boxes(coarse_hierarchy)
    println("coarse leaf ", leaf.index, " -> ", collect(leaf_primitive_indices(coarse_pgdg, leaf.index)))
end

println(refined_pgdg)
for leaf in leaf_boxes(refined_hierarchy)
    println("refined leaf ", leaf.index, " -> ", collect(leaf_primitive_indices(refined_pgdg, leaf.index)))
end

refined_region_count =
    length(leaf_primitive_indices(refined_pgdg, 4)) + length(leaf_primitive_indices(refined_pgdg, 5))
println("coarse primitives in box 1: ", length(leaf_primitive_indices(coarse_pgdg, 1)))
println("refined primitives in former box 1 region: ", refined_region_count)

println("coarse overlap size: ", size(coarse_rep.basis_matrices.overlap))
println("refined overlap size: ", size(refined_rep.basis_matrices.overlap))
println("refined kinetic inf-norm: ", norm(refined_rep.basis_matrices.kinetic, Inf))
