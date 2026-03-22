using GaussletBases

Z = 2.0
s = 0.2
lmax = 2

rb = build_basis(RadialBasisSpec(:G10;
    rmax = 30.0,
    mapping = AsinhMapping(c = s / (2Z), s = s),
))

grid = radial_quadrature(rb)
radial_ops = atomic_operators(rb, grid; Z = Z, lmax = lmax)
ida = atomic_ida_operators(radial_ops; lmax = lmax)

println("channel set: ", ida.one_body.channels)
println("radial basis functions per channel: ", length(rb))
println("total orbitals: ", length(orbitals(ida)))
println("expected channel count for lmax = ", lmax, ": ", sum(2l + 1 for l in 0:lmax))
println("one-body matrix size: ", size(ida.one_body.hamiltonian))
println("radial multipole count: ", length(radial_ops.multipole_data))
println("radial multipole L=1 size: ", size(radial_multipole(ida, 1)))
println("Gaunt tensor L=1 size: ", size(gaunt_tensor(ida, 1)))
println("angular kernel L=1 size: ", size(angular_kernel(ida, 1)))
println("L=1 radial multipoles are radial-only blocks; the angular kernel carries the channel coupling.")
println("Next step: examples/19_atomic_ida_direct.jl uses these ingredients to build a direct matrix from a trial density.")
println("first three orbitals:")
for orbital in orbitals(ida)[1:3]
    println("  ", orbital)
end
