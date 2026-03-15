using LinearAlgebra
using GaussletBases

function midpoint_reference_position_matrix(basis; xmin = -20.0, xmax = 20.0, h = 0.02)
    points = collect((xmin + 0.5 * h):h:(xmax - 0.5 * h))
    weights = fill(h, length(points))
    values = [basis[j](x) for x in points, j in 1:length(basis)]
    return transpose(values) * (((weights .* points)) .* values)
end

ub = build_basis(UniformBasisSpec(:G10; xmin = -1.0, xmax = 1.0, spacing = 1.0))
P = primitive_set(ub)

Xmu = position_matrix(P)
Xb = contract_primitive_matrix(ub, Xmu)
Xref = midpoint_reference_position_matrix(ub)

println("uniform basis length: ", length(ub))
println("primitive set length: ", length(P))
println("contracted position error: ", norm(Xb - Xref, Inf))
