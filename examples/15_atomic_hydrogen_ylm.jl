using LinearAlgebra
using GaussletBases

Z = 1.0
s = 0.2
lmax = 2

rb = build_basis(RadialBasisSpec(:G10;
    rmax = 30.0,
    mapping = AsinhMapping(c = s / (2Z), s = s),
    reference_spacing = 1.0,
    tails = 6,
    odd_even_kmax = 6,
    xgaussians = XGaussian[],
))

grid = radial_quadrature(rb)
radial_ops = atomic_operators(rb, grid; Z = Z, lmax = lmax)
atom = atomic_one_body_operators(radial_ops; lmax = lmax)

energies = sort(real(eigen(Hermitian(atom.hamiltonian)).values))

println("channel set: ", atom.channels)
println("basis functions per channel: ", length(rb))
println("lowest global energies:")
for energy in energies[1:min(end, 10)]
    println("  ", energy)
end

println()
println("lowest energy in each (l,m) channel:")
for channel in atom.channels
    channel_eig = eigen(Hermitian(channel_hamiltonian(atom, channel)))
    println("  (l, m) = (", channel.l, ", ", channel.m, "): ", minimum(real(channel_eig.values)))
end
