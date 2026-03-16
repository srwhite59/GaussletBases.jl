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
- angular Gaunt tensors and M-summed kernels
- orbital indexing metadata

It does **not** solve the many-electron problem. It only assembles the static
ingredients needed for later He / IDA work.
"""
struct AtomicIDAOperators{A <: AbstractDiagonalApproximation}
    one_body::AtomicOneBodyOperators
    radial_operators::RadialAtomicOperators{A}
    gaunt_data::Vector{Array{Float64,3}}
    angular_kernel_data::Vector{Array{Float64,4}}
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
        length(ops.angular_kernel_data) - 1,
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
    L < length(ops.gaunt_data) || throw(BoundsError(ops.gaunt_data, L + 1))
    return ops.gaunt_data[L + 1]
end

function angular_kernel(ops::AtomicIDAOperators, L::Int)
    L >= 0 || throw(ArgumentError("angular_kernel requires L >= 0"))
    L < length(ops.angular_kernel_data) || throw(BoundsError(ops.angular_kernel_data, L + 1))
    return ops.angular_kernel_data[L + 1]
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
    left_index = _channel_index(ops.one_body.channels, left)
    right_index = _channel_index(ops.one_body.channels, right)
    return gaunt_tensor(ops, L)[left_index, right_index, M + L + 1]
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

function _gaunt_entries_by_M(
    table::GauntTable{Float64},
    channels::YlmChannelSet,
    lookup::Dict{Tuple{Int,Int},Int},
    L::Int,
)
    entries_by_M = [Tuple{Int,Int,Float64}[] for _ in 1:(2 * L + 1)]

    for (l1, l2, entries) in gaunt_each_block(table, L)
        for entry in entries
            alpha = lookup[(l1, entry.m1)]
            alphap = lookup[(l2, entry.m2)]
            push!(entries_by_M[entry.M + L + 1], (alpha, alphap, entry.val))
        end
    end

    return entries_by_M
end

function _build_gaunt_tensors(table::GauntTable{Float64}, channels::YlmChannelSet)
    lookup = _channel_lookup(channels)
    return [
        _dense_gaunt_tensor(table, channels, lookup, L)
        for L in 0:gaunt_Lmax(table)
    ]
end

function _build_angular_kernels(table::GauntTable{Float64}, channels::YlmChannelSet)
    nchannels = length(channels)
    lookup = _channel_lookup(channels)
    angular_kernel_data = Vector{Array{Float64,4}}(undef, gaunt_Lmax(table) + 1)

    for L in 0:gaunt_Lmax(table)
        kernel = zeros(Float64, nchannels, nchannels, nchannels, nchannels)
        prefactor = 4 * pi / (2 * L + 1)
        entries_by_M = _gaunt_entries_by_M(table, channels, lookup, L)

        for M in -L:L
            left_entries = entries_by_M[M + L + 1]
            right_entries = entries_by_M[-M + L + 1]
            phase = isodd(M) ? -1.0 : 1.0
            scale = prefactor * phase

            for (alpha, alphap, left_value) in left_entries
                for (beta, betap, right_value) in right_entries
                    kernel[alpha, alphap, beta, betap] += scale * left_value * right_value
                end
            end
        end
        angular_kernel_data[L + 1] = kernel
    end

    return angular_kernel_data
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
- Gaunt tensors for the complex `Y_{lm}` channels
- M-summed angular kernels `Q_L`
- channel-major orbital indexing metadata

This object assembles the Hamiltonian ingredients but does not solve the
many-electron problem.
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
    gaunt_data = _build_gaunt_tensors(gaunt_table, channels)
    angular_kernel_data = _build_angular_kernels(gaunt_table, channels)
    orbital_data = _atomic_orbitals(channels, size(radial_ops.overlap, 1))
    return AtomicIDAOperators(one_body, radial_ops, gaunt_data, angular_kernel_data, orbital_data)
end

function atomic_ida_operators(radial_ops::RadialAtomicOperators; lmax::Int)
    return atomic_ida_operators(radial_ops, ylm_channels(lmax))
end
