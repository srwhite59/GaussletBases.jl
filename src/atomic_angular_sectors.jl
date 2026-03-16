"""
    AngularPairSector

One conserved-`m` pair sector in the complex `Y_{lm}` channel basis.

The sector groups channel pairs `(α, β)` by the conserved sum `m_α + m_β`.
"""
struct AngularPairSector
    msum::Int
    pair_indices::Vector{Int}
    left_channel_indices::Vector{Int}
    right_channel_indices::Vector{Int}
end

"""
    AngularKernelSectors

Sectorized angular-kernel preparation built from the sparse Gaunt table.

For each multipole `L`, the object stores one dense matrix per conserved-`m`
pair sector rather than one dense kernel over the full channel-pair space.
"""
struct AngularKernelSectors
    nchannels::Int
    sectors::Vector{AngularPairSector}
    pair_to_sector::Vector{Int}
    pair_to_local::Vector{Int}
    sector_matrices::Vector{Vector{Matrix{Float64}}}
end

@inline _channel_pair_index(nchannels::Int, left::Int, right::Int) = (left - 1) * nchannels + right

function _pair_groups_by_msum(channels::YlmChannelSet)
    nchannels = length(channels)
    grouped = Dict{Int, Vector{NTuple{3,Int}}}()

    for left in 1:nchannels, right in 1:nchannels
        msum = channels[left].m + channels[right].m
        pair_index = _channel_pair_index(nchannels, left, right)
        push!(get!(grouped, msum, NTuple{3,Int}[]), (left, right, pair_index))
    end

    sectors = AngularPairSector[]
    pair_to_sector = zeros(Int, nchannels^2)
    pair_to_local = zeros(Int, nchannels^2)

    for (sector_index, msum) in enumerate(sort(collect(keys(grouped))))
        tuples = grouped[msum]
        sort!(tuples; by = tuple -> tuple[3])

        pair_indices = Int[]
        left_indices = Int[]
        right_indices = Int[]

        for (local_index, (left, right, pair_index)) in enumerate(tuples)
            push!(pair_indices, pair_index)
            push!(left_indices, left)
            push!(right_indices, right)
            pair_to_sector[pair_index] = sector_index
            pair_to_local[pair_index] = local_index
        end

        push!(sectors, AngularPairSector(msum, pair_indices, left_indices, right_indices))
    end

    return sectors, pair_to_sector, pair_to_local
end

function _gaunt_entries_by_M(
    table::GauntTable{Float64},
    channels::YlmChannelSet,
    lookup::Dict{Tuple{Int,Int},Int},
    L::Int,
)
    entries_by_M = [Tuple{Int,Int,Float64}[] for _ in 1:(2 * L + 1)]

    for (l1, l2, entries) in gaunt_each_block(table, L)
        for entry in entries
            left = lookup[(l1, entry.m1)]
            right = lookup[(l2, entry.m2)]
            push!(entries_by_M[entry.M + L + 1], (left, right, entry.val))
        end
    end

    return entries_by_M
end

function _build_angular_kernel_sectors(
    table::GauntTable{Float64},
    channels::YlmChannelSet,
    lookup::Dict{Tuple{Int,Int},Int},
)
    nchannels = length(channels)
    sectors, pair_to_sector, pair_to_local = _pair_groups_by_msum(channels)
    sector_matrices = Vector{Vector{Matrix{Float64}}}(undef, gaunt_Lmax(table) + 1)

    for L in 0:gaunt_Lmax(table)
        prefactor = 4 * pi / (2 * L + 1)
        entries_by_M = _gaunt_entries_by_M(table, channels, lookup, L)
        level = [zeros(Float64, length(sector.pair_indices), length(sector.pair_indices)) for sector in sectors]

        for M in -L:L
            left_entries = entries_by_M[M + L + 1]
            right_entries = entries_by_M[-M + L + 1]
            scale = prefactor * (isodd(M) ? -1.0 : 1.0)

            for (alpha, alphap, left_value) in left_entries
                for (beta, betap, right_value) in right_entries
                    left_pair = _channel_pair_index(nchannels, alpha, beta)
                    right_pair = _channel_pair_index(nchannels, alphap, betap)
                    sector_index = pair_to_sector[left_pair]
                    pair_to_sector[right_pair] == sector_index || continue
                    local_left = pair_to_local[left_pair]
                    local_right = pair_to_local[right_pair]
                    level[sector_index][local_left, local_right] += scale * left_value * right_value
                end
            end
        end

        sector_matrices[L + 1] = level
    end

    return AngularKernelSectors(nchannels, sectors, pair_to_sector, pair_to_local, sector_matrices)
end

function _dense_angular_kernel_from_sectors(sectors::AngularKernelSectors, L::Int)
    L >= 0 || throw(ArgumentError("_dense_angular_kernel_from_sectors requires L >= 0"))
    L < length(sectors.sector_matrices) || throw(BoundsError(sectors.sector_matrices, L + 1))

    kernel = zeros(Float64, sectors.nchannels, sectors.nchannels, sectors.nchannels, sectors.nchannels)
    level = sectors.sector_matrices[L + 1]

    for (sector_index, sector) in enumerate(sectors.sectors)
        matrix = level[sector_index]
        for left_local in eachindex(sector.pair_indices)
            alpha = sector.left_channel_indices[left_local]
            beta = sector.right_channel_indices[left_local]
            for right_local in eachindex(sector.pair_indices)
                alphap = sector.left_channel_indices[right_local]
                betap = sector.right_channel_indices[right_local]
                kernel[alpha, alphap, beta, betap] = matrix[left_local, right_local]
            end
        end
    end

    return kernel
end

function _angular_kernel_sector_value(
    sectors::AngularKernelSectors,
    L::Int,
    alpha::Int,
    beta::Int,
    alphap::Int,
    betap::Int,
)
    left_pair = _channel_pair_index(sectors.nchannels, alpha, beta)
    right_pair = _channel_pair_index(sectors.nchannels, alphap, betap)
    left_sector = sectors.pair_to_sector[left_pair]
    right_sector = sectors.pair_to_sector[right_pair]
    left_sector == right_sector || return 0.0
    left_local = sectors.pair_to_local[left_pair]
    right_local = sectors.pair_to_local[right_pair]
    return sectors.sector_matrices[L + 1][left_sector][left_local, right_local]
end

