using GaussletBases

ub = build_basis(UniformBasisSpec(:G10; xmin = -2.0, xmax = 2.0, spacing = 1.0))
rep = basis_representation(ub)

hierarchy = hierarchical_partition(rep, [-2.5, -0.5, 0.5, 2.5])
hierarchy = refine_partition(hierarchy, 1)

println(hierarchy)
for box in boxes(hierarchy)
    println(box)
    println("  parent: ", box_parent(hierarchy, box.index))
    println("  children: ", box_children(hierarchy, box.index))
    println("  basis indices: ", box_indices(hierarchy, box.index))
end

println("leaf boxes:")
for box in leaf_boxes(hierarchy)
    println("  ", box)
end

println("overlap block for refined leaf 4:")
println(box_block(rep, hierarchy, :overlap, 4))

println("kinetic coupling between refined leaf 4 and coarse box 2:")
println(box_coupling(rep, hierarchy, :kinetic, 4, 2))
