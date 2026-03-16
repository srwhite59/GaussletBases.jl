"""
    AtomicOrbital

One atomic orbital index in the first static IDA assembly layer.

The ordering is channel-major: channels follow the `YlmChannelSet`, and within
each channel the radial basis index runs from `1:nradial`.
"""
struct AtomicOrbital
    index::Int
    channel::YlmChannel
    radial_index::Int

    function AtomicOrbital(index::Int, channel::YlmChannel, radial_index::Int)
        index >= 1 || throw(ArgumentError("AtomicOrbital requires index >= 1"))
        radial_index >= 1 || throw(ArgumentError("AtomicOrbital requires radial_index >= 1"))
        return new(index, channel, radial_index)
    end
end

function Base.show(io::IO, orbital::AtomicOrbital)
    print(
        io,
        "AtomicOrbital(index=",
        orbital.index,
        ", channel=",
        orbital.channel,
        ", radial_index=",
        orbital.radial_index,
        ")",
    )
end

"""
    AtomicIDAOperators

Static interacting atomic IDA ingredients built on top of the radial substrate
and the explicit `(l,m)` channel layer.

The object bundles:

- the one-body atomic blocks
- the underlying radial multipole tables
- a sparse/block Gaunt table
- sectorized M-summed angular-kernel preparation
- orbital indexing metadata

It does **not** solve the many-electron problem. It only assembles the static
ingredients needed for later He / IDA work.
"""
struct AtomicIDAOperators{A <: AbstractDiagonalApproximation}
    one_body::AtomicOneBodyOperators
    radial_operators::RadialAtomicOperators{A}
    gaunt_table::GauntTable{Float64}
    angular_sectors::AngularKernelSectors
    channel_lookup::Dict{Tuple{Int,Int},Int}
    orbital_data::Vector{AtomicOrbital}
end

function Base.show(io::IO, ops::AtomicIDAOperators)
    print(
        io,
        "AtomicIDAOperators(nchannels=",
        length(ops.one_body.channels),
        ", norbitals=",
        length(ops.orbital_data),
        ", Lmax=",
        gaunt_Lmax(ops.gaunt_table),
        ")",
    )
end

orbitals(ops::AtomicIDAOperators) = ops.orbital_data

channel_range(ops::AtomicIDAOperators, index::Int) = channel_range(ops.one_body, index)
channel_range(ops::AtomicIDAOperators, channel::YlmChannel) = channel_range(ops.one_body, channel)
channel_hamiltonian(ops::AtomicIDAOperators, index::Int) = channel_hamiltonian(ops.one_body, index)
channel_hamiltonian(ops::AtomicIDAOperators, channel::YlmChannel) = channel_hamiltonian(ops.one_body, channel)
channel_overlap(ops::AtomicIDAOperators, index::Int) = channel_overlap(ops.one_body, index)
channel_overlap(ops::AtomicIDAOperators, channel::YlmChannel) = channel_overlap(ops.one_body, channel)

function radial_multipole(ops::AtomicIDAOperators, L::Int)
    return multipole(ops.radial_operators, L)
end

function gaunt_tensor(ops::AtomicIDAOperators, L::Int)
    L >= 0 || throw(ArgumentError("gaunt_tensor requires L >= 0"))
    L <= gaunt_Lmax(ops.gaunt_table) || throw(BoundsError(ops.gaunt_table.blocks, L + 1))
    return _dense_gaunt_tensor(ops.gaunt_table, ops.one_body.channels, ops.channel_lookup, L)
end

function angular_kernel(ops::AtomicIDAOperators, L::Int)
    L >= 0 || throw(ArgumentError("angular_kernel requires L >= 0"))
    L <= gaunt_Lmax(ops.gaunt_table) || throw(BoundsError(ops.gaunt_table.blocks, L + 1))
    return _dense_angular_kernel_from_sectors(ops.angular_sectors, L)
end

function gaunt_coefficient(
    ops::AtomicIDAOperators,
    L::Int,
    M::Int,
    left::YlmChannel,
    right::YlmChannel,
)
    L >= 0 || throw(ArgumentError("gaunt_coefficient requires L >= 0"))
    abs(M) <= L || throw(ArgumentError("gaunt_coefficient requires |M| <= L"))
    return gaunt_value(ops.gaunt_table, L, left.l, left.m, right.l, right.m, M)
end

function _channel_lookup(channels::YlmChannelSet)
    lookup = Dict{Tuple{Int,Int},Int}()
    for index in eachindex(channels.channel_data)
        channel = channels[index]
        lookup[(channel.l, channel.m)] = index
    end
    return lookup
end

function _dense_gaunt_tensor(
    table::GauntTable{Float64},
    channels::YlmChannelSet,
    lookup::Dict{Tuple{Int,Int},Int},
    L::Int,
)
    nchannels = length(channels)
    tensor = zeros(Float64, nchannels, nchannels, 2 * L + 1)

    for (l1, l2, entries) in gaunt_each_block(table, L)
        for entry in entries
            alpha = lookup[(l1, entry.m1)]
            alphap = lookup[(l2, entry.m2)]
            tensor[alpha, alphap, entry.M + L + 1] = entry.val
        end
    end

    return tensor
end

function _atomic_orbitals(channels::YlmChannelSet, radial_dim::Int)
    orbital_data = AtomicOrbital[]
    index = 1
    for channel in channels
        for radial_index in 1:radial_dim
            push!(orbital_data, AtomicOrbital(index, channel, radial_index))
            index += 1
        end
    end
    return orbital_data
end

"""
    atomic_ida_operators(radial_ops::RadialAtomicOperators, channels::YlmChannelSet)
    atomic_ida_operators(radial_ops::RadialAtomicOperators; lmax::Int)

Build the first static interacting atomic IDA ingredient bundle on top of the
current radial and angular layers.

The result contains:

- the one-body atomic blocks
- the radial multipole tables from `radial_ops`
- sparse/block Gaunt data for the complex `Y_{lm}` channels
- sectorized M-summed angular kernels `Q_L`
- channel-major orbital indexing metadata

This object assembles the Hamiltonian ingredients but does not solve the
many-electron problem. The dense `gaunt_tensor` and `angular_kernel` accessors
remain available and reconstruct those objects on demand from the sparse and
sectorized internal data.
"""
function atomic_ida_operators(radial_ops::RadialAtomicOperators, channels::YlmChannelSet)
    maximum_l = _maximum_l(channels)
    maximum_l >= 0 || throw(ArgumentError("atomic_ida_operators requires at least one channel"))
    2 * maximum_l < length(radial_ops.multipole_data) ||
        throw(
            ArgumentError(
                "interacting atomic assembly requires radial multipoles through L = $(2 * maximum_l), " *
                "but the supplied RadialAtomicOperators only contain data through L = $(length(radial_ops.multipole_data) - 1)",
            ),
        )

    one_body = atomic_one_body_operators(radial_ops, channels)
    gaunt_table = build_gaunt_table(maximum_l; Lmax = 2 * maximum_l, atol = 1.0e-14, basis = :complex)
    lookup = _channel_lookup(channels)
    angular_sectors = _build_angular_kernel_sectors(gaunt_table, channels, lookup)
    orbital_data = _atomic_orbitals(channels, size(radial_ops.overlap, 1))
    return AtomicIDAOperators(one_body, radial_ops, gaunt_table, angular_sectors, lookup, orbital_data)
end

function atomic_ida_operators(radial_ops::RadialAtomicOperators; lmax::Int)
    return atomic_ida_operators(radial_ops, ylm_channels(lmax))
end
