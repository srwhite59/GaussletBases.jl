using GaussletBases
using LinearAlgebra

Z = 2.0
s = 0.5
lmax = 1

rb = build_basis(RadialBasisSpec(:G10;
    rmax = 8.0,
    mapping = AsinhMapping(c = s / (2Z), s = s),
    reference_spacing = 1.0,
    tails = 6,
    odd_even_kmax = 6,
))

grid = radial_quadrature(rb; accuracy = :high)
radial_ops = atomic_operators(rb, grid; Z = Z, lmax = lmax)
ida = atomic_ida_operators(radial_ops; lmax = lmax)

radial_dim = length(rb)
nchannels = length(ida.one_body.channels)
norbitals = length(orbitals(ida))
orbital_index(channel_index, radial_index) = (channel_index - 1) * radial_dim + radial_index

channel_labels = [(channel.l, channel.m) for channel in ida.one_body.channels]

# Positive-semidefinite trial spatial density built from two inspectable orbitals.
orbital1 = zeros(Float64, norbitals)
orbital1[orbital_index(1, 1)] = sqrt(0.75)  # channel 1 = (0,0), first radial function
orbital1[orbital_index(3, 1)] = sqrt(0.25)  # channel 3 = (1,0), first radial function

orbital2 = zeros(Float64, norbitals)
orbital2[orbital_index(4, 2)] = 1.0         # channel 4 = (1,1), second radial function

density = 1.4 .* (orbital1 * orbital1') .+ 0.4 .* (orbital2 * orbital2')

direct = direct_matrix(ida, density)
exchange = exchange_matrix(ida, density)
fock = fock_matrix(ida, density)
density_eigmin = minimum(real(eigen(Hermitian(density)).values))

println("channel set: ", ida.one_body.channels)
println("channel index map:")
for (index, label) in enumerate(channel_labels)
    println("  ", index, " => (l, m) = ", label)
end
println("spatial orbitals: ", norbitals)
println("trial density minimum eigenvalue: ", density_eigmin)
println("one-body size: ", size(ida.one_body.hamiltonian))
println("direct size: ", size(direct))
println("exchange size: ", size(exchange))
println("fock size: ", size(fock))
println("||J||_inf: ", norm(direct, Inf))
println("||K||_inf: ", norm(exchange, Inf))
println("||F - (h + J - K)||_inf: ", norm(fock - (ida.one_body.hamiltonian + direct - exchange), Inf))
println("Fock Hermiticity error: ", norm(fock - transpose(fock), Inf))
println("first radial-pair Fock block (p=1, q=1):")
println(fock[[orbital_index(channel, 1) for channel in 1:nchannels], [orbital_index(channel, 1) for channel in 1:nchannels]])
