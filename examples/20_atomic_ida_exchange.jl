using GaussletBases
using LinearAlgebra

function dense_exchange_reference(ops::AtomicIDAOperators, density::AbstractMatrix{<:Real})
    radial_dim = size(radial_multipole(ops, 0), 1)
    nchannels = length(ops.one_body.channels)
    norbitals = length(orbitals(ops))
    reference = zeros(Float64, norbitals, norbitals)
    radial_multipoles = [radial_multipole(ops, L) for L in 0:GaussletBases.gaunt_Lmax(ops.gaunt_table)]
    angular_kernels = [angular_kernel(ops, L) for L in 0:GaussletBases.gaunt_Lmax(ops.gaunt_table)]

    orbital_index(channel_index, radial_index) = (channel_index - 1) * radial_dim + radial_index

    for left_radial in 1:radial_dim, right_radial in 1:radial_dim, left_channel in 1:nchannels, right_channel in 1:nchannels
        total = 0.0
        for (level_index, multipole) in enumerate(radial_multipoles), left_source in 1:nchannels, right_source in 1:nchannels
            total += multipole[left_radial, right_radial] *
                     angular_kernels[level_index][left_channel, left_source, right_channel, right_source] *
                     density[orbital_index(left_source, left_radial), orbital_index(right_source, right_radial)]
        end
        reference[orbital_index(left_channel, left_radial), orbital_index(right_channel, right_radial)] = total
    end

    return 0.5 .* (reference .+ transpose(reference))
end

Z = 2.0
s = 0.5
lmax = 1

rb = build_basis(RadialBasisSpec(:G10;
    rmax = 8.0,
    mapping = AsinhMapping(c = s / (2Z), s = s),
    reference_spacing = 1.0,
    tails = 6,
    odd_even_kmax = 6,
    xgaussians = XGaussian[],
))

grid = radial_quadrature(rb; accuracy = :high)
radial_ops = atomic_operators(rb, grid; Z = Z, lmax = lmax)
ida = atomic_ida_operators(radial_ops; lmax = lmax)

radial_dim = length(rb)
nchannels = length(ida.one_body.channels)
norbitals = length(orbitals(ida))
orbital_index(channel_index, radial_index) = (channel_index - 1) * radial_dim + radial_index

density = zeros(Float64, norbitals, norbitals)
density[orbital_index(1, 1), orbital_index(1, 1)] = 0.8
density[orbital_index(3, 1), orbital_index(1, 1)] = 0.25
density[orbital_index(1, 1), orbital_index(3, 1)] = 0.25
density[orbital_index(1, 1), orbital_index(3, 2)] = -0.2
density[orbital_index(3, 2), orbital_index(1, 1)] = -0.2
density[orbital_index(4, 2), orbital_index(4, 2)] = 0.4

exchange = exchange_matrix(ida, density)
reference = dense_exchange_reference(ida, density)

hermiticity_error = norm(exchange - transpose(exchange), Inf)
comparison_error = norm(exchange - reference, Inf)

println("channel set: ", ida.one_body.channels)
println("radial basis functions per channel: ", radial_dim)
println("spatial orbitals: ", norbitals)
println("exchange matrix size: ", size(exchange))
println("exchange Hermiticity error: ", hermiticity_error)
println("sectorized vs dense difference: ", comparison_error)
println("first radial-pair exchange block (p=1, q=1):")
println(exchange[[orbital_index(channel, 1) for channel in 1:nchannels], [orbital_index(channel, 1) for channel in 1:nchannels]])
