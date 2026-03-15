using LinearAlgebra
using GaussletBases

ub = build_basis(UniformBasisSpec(:G10; xmin = -1.0, xmax = 1.0, spacing = 1.0))
rep = basis_representation(ub)

println(rep)
println("basis kind: ", rep.metadata.basis_kind)
println("basis functions: ", size(rep.coefficient_matrix, 2))
println("primitive functions: ", length(rep.primitive_set))
println("available matrices: ", collect(keys(rep.basis_matrices)))
println("overlap deviation: ", norm(rep.basis_matrices.overlap - I, Inf))
println("position diagonal: ", diag(rep.basis_matrices.position))
