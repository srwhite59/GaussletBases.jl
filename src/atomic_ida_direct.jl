@inline function _atomic_orbital_index(radial_dim::Int, channel_index::Int, radial_index::Int)
    return (channel_index - 1) * radial_dim + radial_index
end

function _radial_diagonal_density_blocks(ops::AtomicIDAOperators, density::AbstractMatrix{<:Real})
    norbitals = length(ops.orbital_data)
    size(density) == (norbitals, norbitals) ||
        throw(DimensionMismatch("density matrix must have size ($(norbitals), $(norbitals))"))

    radial_dim = size(radial_multipole(ops, 0), 1)
    nchannels = length(ops.one_body.channels)
    blocks = [zeros(Float64, nchannels, nchannels) for _ in 1:radial_dim]

    for radial_index in 1:radial_dim
        block = blocks[radial_index]
        for left_channel in 1:nchannels, right_channel in 1:nchannels
            left = _atomic_orbital_index(radial_dim, left_channel, radial_index)
            right = _atomic_orbital_index(radial_dim, right_channel, radial_index)
            block[left_channel, right_channel] = density[left, right]
        end
    end

    return blocks
end

function _sectorized_direct_blocks(
    ops::AtomicIDAOperators,
    density_blocks::AbstractVector{<:AbstractMatrix{<:Real}},
)
    radial_dim = size(radial_multipole(ops, 0), 1)
    nchannels = length(ops.one_body.channels)
    length(density_blocks) == radial_dim ||
        throw(DimensionMismatch("density block count must match the radial dimension"))
    all(size(block) == (nchannels, nchannels) for block in density_blocks) ||
        throw(DimensionMismatch("each radial density block must have size ($(nchannels), $(nchannels))"))

    direct_blocks = [zeros(Float64, nchannels, nchannels) for _ in 1:radial_dim]

    for L in 0:gaunt_Lmax(ops.gaunt_table)
        entries_by_M = _gaunt_entries_by_M(ops.gaunt_table, ops.one_body.channels, ops.channel_lookup, L)
        density_coefficients = zeros(Float64, radial_dim, 2 * L + 1)

        for slot in eachindex(entries_by_M)
            entries = entries_by_M[slot]
            for (left_channel, right_channel, value) in entries
                for radial_index in 1:radial_dim
                    density_coefficients[radial_index, slot] +=
                        value * density_blocks[radial_index][left_channel, right_channel]
                end
            end
        end

        weighted_coefficients = radial_multipole(ops, L) * density_coefficients
        prefactor = 4 * pi / (2 * L + 1)

        for M in -L:L
            slot = M + L + 1
            partner_slot = -M + L + 1
            scale = prefactor * (isodd(M) ? -1.0 : 1.0)

            for (left_channel, right_channel, value) in entries_by_M[slot]
                for radial_index in 1:radial_dim
                    direct_blocks[radial_index][left_channel, right_channel] +=
                        scale * value * weighted_coefficients[radial_index, partner_slot]
                end
            end
        end
    end

    return direct_blocks
end

"""
    direct_matrix(ops::AtomicIDAOperators, density::AbstractMatrix)

Build the direct/Hartree one-body matrix in the current channel-major orbital
basis from a spatial one-particle density matrix.

The input density is interpreted in the spatial-orbital basis exposed by
`orbitals(ops)`. Within the present local diagonal approximation, only the
radial-diagonal density blocks contribute:

- `density[(p, α), (p, β)]` contributes
- `density[(p, α), (q, β)]` with `p != q` is ignored

The result is an effective one-body matrix on that same spatial-orbital basis.
It is block diagonal in the radial index and uses the sparse Gaunt table
directly rather than the dense `angular_kernel` view.
"""
function direct_matrix(ops::AtomicIDAOperators, density::AbstractMatrix{<:Real})
    radial_dim = size(radial_multipole(ops, 0), 1)
    nchannels = length(ops.one_body.channels)
    density_blocks = _radial_diagonal_density_blocks(ops, density)
    direct_blocks = _sectorized_direct_blocks(ops, density_blocks)

    norbitals = length(ops.orbital_data)
    direct = zeros(Float64, norbitals, norbitals)

    for radial_index in 1:radial_dim
        block = direct_blocks[radial_index]
        for left_channel in 1:nchannels, right_channel in 1:nchannels
            left = _atomic_orbital_index(radial_dim, left_channel, radial_index)
            right = _atomic_orbital_index(radial_dim, right_channel, radial_index)
            direct[left, right] = block[left_channel, right_channel]
        end
    end

    return 0.5 .* (direct .+ transpose(direct))
end
