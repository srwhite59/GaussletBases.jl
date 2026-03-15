using LinearAlgebra
using GaussletBases

function midpoint_reference_matrices(basis; xmin = -20.0, xmax = 20.0, h = 0.02)
    points = collect((xmin + 0.5 * h):h:(xmax - 0.5 * h))
    weights = fill(h, length(points))
    values = [basis[j](x) for x in points, j in 1:length(basis)]
    derivatives = [derivative(basis[j], x) for x in points, j in 1:length(basis)]
    overlap = transpose(values) * (weights .* values)
    kinetic = 0.5 .* (transpose(derivatives) * (weights .* derivatives))
    return overlap, kinetic
end

ub = build_basis(UniformBasisSpec(:G10; xmin = -1.0, xmax = 1.0, spacing = 1.0))
P = primitive_set(ub)

Smu = overlap_matrix(P)
Tmu = kinetic_matrix(P)

Sb = contract_primitive_matrix(ub, Smu)
Tb = contract_primitive_matrix(ub, Tmu)

Sref, Tref = midpoint_reference_matrices(ub)

println("uniform basis length: ", length(ub))
println("primitive set length: ", length(P))
println("contracted overlap error: ", norm(Sb - Sref, Inf))
println("contracted kinetic error: ", norm(Tb - Tref, Inf))
