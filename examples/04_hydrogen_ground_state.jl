using LinearAlgebra
using GaussletBases

Z = 1.0
s = 0.2

map = AsinhMapping(c = s / (2Z), s = s)

rb = build_basis(RadialBasisSpec(:G10;
    rmax = 30.0,
    mapping = map,
))

grid = radial_quadrature(rb)

S = overlap_matrix(rb, grid)
H = kinetic_matrix(rb, grid) +
    nuclear_matrix(rb, grid; Z = Z) +
    centrifugal_matrix(rb, grid; l = 0)

E0 = minimum(real(eigen(Hermitian(H)).values))

println("Hydrogen ground-state energy: ", E0)
println("Error relative to -0.5 Ha: ", E0 + 0.5)
