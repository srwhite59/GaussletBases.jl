function _radial_pair_density_blocks(ops::AtomicIDAOperators, density::AbstractMatrix{<:Real})
    norbitals = length(ops.orbital_data)
    size(density) == (norbitals, norbitals) ||
        throw(DimensionMismatch("density matrix must have size ($(norbitals), $(norbitals))"))

    radial_dim = size(radial_multipole(ops, 0), 1)
    nchannels = length(ops.one_body.channels)
    blocks = Matrix{Matrix{Float64}}(undef, radial_dim, radial_dim)

    for left_radial in 1:radial_dim, right_radial in 1:radial_dim
        block = zeros(Float64, nchannels, nchannels)
        for left_channel in 1:nchannels, right_channel in 1:nchannels
            left = _atomic_orbital_index(radial_dim, left_channel, left_radial)
            right = _atomic_orbital_index(radial_dim, right_channel, right_radial)
            block[left_channel, right_channel] = density[left, right]
        end
        blocks[left_radial, right_radial] = block
    end

    return blocks
end

function _sectorized_exchange_blocks(
    ops::AtomicIDAOperators,
    density_blocks::AbstractMatrix{<:AbstractMatrix{<:Real}},
)
    radial_dim = size(radial_multipole(ops, 0), 1)
    nchannels = length(ops.one_body.channels)
    size(density_blocks) == (radial_dim, radial_dim) ||
        throw(DimensionMismatch("density block grid must have size ($(radial_dim), $(radial_dim))"))
    all(size(block) == (nchannels, nchannels) for block in density_blocks) ||
        throw(DimensionMismatch("each radial-pair density block must have size ($(nchannels), $(nchannels))"))

    exchange_blocks = Matrix{Matrix{Float64}}(undef, radial_dim, radial_dim)
    sectors = ops.angular_sectors
    radial_multipoles = [radial_multipole(ops, L) for L in 0:gaunt_Lmax(ops.gaunt_table)]

    for left_radial in 1:radial_dim, right_radial in 1:radial_dim
        density_block = density_blocks[left_radial, right_radial]
        exchange_block = zeros(Float64, nchannels, nchannels)

        for (sector_index, sector) in enumerate(sectors.sectors)
            x = Vector{Float64}(undef, length(sector.pair_indices))
            for local_index in eachindex(sector.pair_indices)
                left_channel = sector.left_channel_indices[local_index]
                right_channel = sector.right_channel_indices[local_index]
                x[local_index] = density_block[left_channel, right_channel]
            end

            y = zeros(Float64, length(x))
            for (level_index, multipole) in enumerate(radial_multipoles)
                coefficient = multipole[left_radial, right_radial]
                coefficient == 0.0 && continue
                y .+= coefficient .* (sectors.sector_matrices[level_index][sector_index] * x)
            end

            for local_index in eachindex(sector.pair_indices)
                left_channel = sector.left_channel_indices[local_index]
                right_channel = sector.right_channel_indices[local_index]
                exchange_block[left_channel, right_channel] = y[local_index]
            end
        end

        exchange_blocks[left_radial, right_radial] = exchange_block
    end

    return exchange_blocks
end

"""
    exchange_matrix(ops::AtomicIDAOperators, density::AbstractMatrix)

Build the exchange/Fock-style one-body matrix in the current atomic IDA/local
diagonal approximation from a spatial one-particle density matrix.

The input density is interpreted in the spatial-orbital basis exposed by
`orbitals(ops)`. Unlike `direct_matrix`, this exchange term uses the full
radial-pair density blocks

```math
\rho_{(p,\alpha),(q,\beta)}
```

and produces a full two-site one-body matrix on that same spatial-orbital
basis. The current approximation still comes from the local diagonal radial
multipole tables in `AtomicIDAOperators`; it is not a fully general Coulomb
contraction beyond the present atomic IDA structure.
"""
function exchange_matrix(ops::AtomicIDAOperators, density::AbstractMatrix{<:Real})
    radial_dim = size(radial_multipole(ops, 0), 1)
    nchannels = length(ops.one_body.channels)
    density_blocks = _radial_pair_density_blocks(ops, density)
    exchange_blocks = _sectorized_exchange_blocks(ops, density_blocks)

    norbitals = length(ops.orbital_data)
    exchange = zeros(Float64, norbitals, norbitals)

    for left_radial in 1:radial_dim, right_radial in 1:radial_dim
        block = exchange_blocks[left_radial, right_radial]
        for left_channel in 1:nchannels, right_channel in 1:nchannels
            left = _atomic_orbital_index(radial_dim, left_channel, left_radial)
            right = _atomic_orbital_index(radial_dim, right_channel, right_radial)
            exchange[left, right] = block[left_channel, right_channel]
        end
    end

    return 0.5 .* (exchange .+ transpose(exchange))
end
