using GaussletBases

map = AsinhMapping(c = 0.15, s = 0.15)
rb = build_basis(RadialBasisSpec(:G10;
    count = 9,
    mapping = map,
))
grid = radial_quadrature(rb)
ops = atomic_operators(rb, grid; Z = 2.0, lmax = 2)

println(ops)
println("overlap size: ", size(ops.overlap))
println("kinetic size: ", size(ops.kinetic))
println("nuclear size: ", size(ops.nuclear))
println("centrifugal l=2 size: ", size(centrifugal(ops, 2)))
println("multipole L=1 size: ", size(multipole(ops, 1)))
