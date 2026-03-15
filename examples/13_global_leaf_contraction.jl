using GaussletBases

mapping = AsinhMapping(c = 0.15, s = 0.15)
global_layer = build_global_mapped_primitive_layer(
    xmin = -2.0,
    xmax = 2.0,
    mapping = mapping,
    reference_spacing = 0.5,
    width_scale = 1.0,
)
global_rep = basis_representation(global_layer)

hierarchy = refine_partition(hierarchical_partition(global_layer, [-2.5, -0.5, 0.5, 2.5]), 1)
contracted = contract_leaf_boxes(global_layer, hierarchy; retained_per_leaf = 1)
contracted_rep = basis_representation(contracted)

println(global_layer)
println("global overlap size: ", size(global_rep.basis_matrices.overlap))
println("global kinetic size: ", size(global_rep.basis_matrices.kinetic))

println(contracted)
for contraction in leaf_contractions(contracted)
    println(contraction)
end

println("contracted overlap size: ", size(contracted_rep.basis_matrices.overlap))
println("contracted kinetic size: ", size(contracted_rep.basis_matrices.kinetic))
println("global primitive count: ", length(primitive_set(global_layer)))
println("contracted basis count: ", size(stencil_matrix(contracted), 2))
