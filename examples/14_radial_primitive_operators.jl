using LinearAlgebra
using GaussletBases

Z = 2.0
rb = build_basis(RadialBasisSpec(:G10;
    count = 6,
    mapping = AsinhMapping(c = 0.15, s = 0.15),
    reference_spacing = 1.0,
    tails = 3,
    odd_even_kmax = 2,
    xgaussians = [XGaussian(alpha = 0.2)],
))
grid = radial_quadrature(rb; accuracy = :medium)
P = primitive_set(rb)

Smu = overlap_matrix(P, grid)
Tmu = kinetic_matrix(P, grid)
Vmu = nuclear_matrix(P, grid; Z = Z)
C2mu = centrifugal_matrix(P, grid; l = 2)

S = contract_primitive_matrix(rb, Smu)
T = contract_primitive_matrix(rb, Tmu)
V = contract_primitive_matrix(rb, Vmu)
C2 = contract_primitive_matrix(rb, C2mu)

Sdirect = overlap_matrix(rb, grid)
Tdirect = kinetic_matrix(rb, grid)
Vdirect = nuclear_matrix(rb, grid; Z = Z)
C2direct = centrifugal_matrix(rb, grid; l = 2)

relative_difference(A, B) = norm(A - B, Inf) / max(norm(B, Inf), 1.0)

println("primitive count: ", length(P))
println("basis count: ", length(rb))
println("overlap difference: ", norm(S - Sdirect, Inf))
println("overlap relative difference: ", relative_difference(S, Sdirect))
println("kinetic difference: ", norm(T - Tdirect, Inf))
println("kinetic relative difference: ", relative_difference(T, Tdirect))
println("nuclear difference: ", norm(V - Vdirect, Inf))
println("nuclear relative difference: ", relative_difference(V, Vdirect))
println("centrifugal difference: ", norm(C2 - C2direct, Inf))
println("centrifugal relative difference: ", relative_difference(C2, C2direct))
