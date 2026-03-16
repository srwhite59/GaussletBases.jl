using GaussletBases
using LinearAlgebra

function dense_direct_reference(ops::AtomicIDAOperators, density::AbstractMatrix{<:Real})
    radial_dim = size(radial_multipole(ops, 0), 1)
    nchannels = length(ops.one_body.channels)
    norbitals = length(orbitals(ops))
    reference = zeros(Float64, norbitals, norbitals)

    orbital_index(channel_index, radial_index) = (channel_index - 1) * radial_dim + radial_index

    for left_radial in 1:radial_dim, left_channel in 1:nchannels, right_channel in 1:nchannels
        total = 0.0
        for right_radial in 1:radial_dim, L in 0:GaussletBases.gaunt_Lmax(ops.gaunt_table), beta in 1:nchannels, betap in 1:nchannels
            density_value = density[orbital_index(beta, right_radial), orbital_index(betap, right_radial)]
            total += radial_multipole(ops, L)[left_radial, right_radial] *
                     angular_kernel(ops, L)[left_channel, right_channel, beta, betap] *
                     density_value
        end
        reference[orbital_index(left_channel, left_radial), orbital_index(right_channel, left_radial)] = total
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

# One simple inspectable trial density: a rank-one projector with mild s/p mixing.
trial_orbital = zeros(Float64, norbitals)
trial_orbital[orbital_index(1, 1)] = sqrt(0.8)   # (l,m) = (0,0), first radial function
trial_orbital[orbital_index(3, 1)] = sqrt(0.2)   # (l,m) = (1,0), first radial function
density = 2.0 .* (trial_orbital * trial_orbital')

direct = direct_matrix(ida, density)
reference = dense_direct_reference(ida, density)

hermiticity_error = norm(direct - transpose(direct), Inf)
comparison_error = norm(direct - reference, Inf)
off_radial_norm = let maximum_off_radial = 0.0
    for left_channel in 1:nchannels, right_channel in 1:nchannels, left_radial in 1:radial_dim, right_radial in 1:radial_dim
        left_radial == right_radial && continue
        maximum_off_radial = max(
            maximum_off_radial,
            abs(direct[orbital_index(left_channel, left_radial), orbital_index(right_channel, right_radial)]),
        )
    end
    maximum_off_radial
end

println("channel set: ", ida.one_body.channels)
println("radial basis functions per channel: ", radial_dim)
println("spatial orbitals: ", norbitals)
println("density rank: ", rank(density))
println("direct matrix size: ", size(direct))
println("direct Hermiticity error: ", hermiticity_error)
println("sectorized vs dense difference: ", comparison_error)
println("largest off-radial entry: ", off_radial_norm)
println("first radial direct block:")
println(direct[[orbital_index(channel, 1) for channel in 1:nchannels], [orbital_index(channel, 1) for channel in 1:nchannels]])
