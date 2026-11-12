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
grid = radial_quadrature(rb)
P = primitive_set(rb)

Smu = overlap_matrix(P, grid)
Tmu = kinetic_matrix(P, grid)
Vmu = nuclear_matrix(P, grid; Z = Z)

S = contract_primitive_matrix(rb, Smu)
T = contract_primitive_matrix(rb, Tmu)
V = contract_primitive_matrix(rb, Vmu)

Sdirect = overlap_matrix(rb, grid)
Tdirect = kinetic_matrix(rb, grid)
Vdirect = nuclear_matrix(rb, grid; Z = Z)
C2direct = centrifugal_matrix(rb, grid; l = 2)
# Nonvanishing primitive origin values make continuum inverse-square contraction ill-conditioned.
h = 12/512
finite_grid = RadialQuadratureGrid(h .* ((1:512) .- 0.5), fill(h,512); mapping=mapping(rb))
C2finite = centrifugal_matrix(rb, finite_grid; l = 2)
C2 = contract_primitive_matrix(rb, centrifugal_matrix(P, finite_grid; l = 2))

relative_difference(A, B) = norm(A - B, Inf) / max(norm(B, Inf), 1.0)

println("primitive count: ", length(P))
println("basis count: ", length(rb))
println("overlap difference: ", norm(S - Sdirect, Inf))
println("overlap relative difference: ", relative_difference(S, Sdirect))
println("kinetic difference: ", norm(T - Tdirect, Inf))
println("kinetic relative difference: ", relative_difference(T, Tdirect))
println("nuclear difference: ", norm(V - Vdirect, Inf))
println("nuclear relative difference: ", relative_difference(V, Vdirect))
println("default direct-basis centrifugal norm: ", norm(C2direct, Inf))
println("finite-grid centrifugal difference (512 midpoints, 0..12 bohr): ", norm(C2 - C2finite, Inf))
println("finite-grid centrifugal relative difference: ", relative_difference(C2, C2finite))
