using GaussletBases

ub = build_basis(UniformBasisSpec(:G10; xmin = -2.0, xmax = 2.0, spacing = 1.0))
rep = basis_representation(ub)
partition = basis_partition(rep, [-2.5, -0.5, 0.5, 2.5])

println(partition)
for (ibox, box) in enumerate(boxes(partition))
    println("box ", ibox, ": ", box)
    println("  basis indices: ", box_indices(partition, ibox))
end

println("overlap block (box 1):")
println(box_block(rep, partition, :overlap, 1))

println("kinetic coupling (box 1 -> box 2):")
println(box_coupling(rep, partition, :kinetic, 1, 2))
