"""
    YlmChannel

Angular channel labeled by spherical-harmonic quantum numbers `(l, m)`.
"""
struct YlmChannel
    l::Int
    m::Int

    function YlmChannel(l::Int, m::Int)
        l >= 0 || throw(ArgumentError("YlmChannel requires l >= 0"))
        abs(m) <= l || throw(ArgumentError("YlmChannel requires |m| <= l"))
        return new(l, m)
    end
end

Base.:(==)(a::YlmChannel, b::YlmChannel) = a.l == b.l && a.m == b.m
Base.hash(channel::YlmChannel, h::UInt) = hash((channel.l, channel.m), h)

function Base.show(io::IO, channel::YlmChannel)
    print(io, "YlmChannel(l=", channel.l, ", m=", channel.m, ")")
end

"""
    YlmChannelSet

Ordered list of spherical-harmonic channels up to a fixed `lmax`.
"""
struct YlmChannelSet
    lmax::Int
    channel_data::Vector{YlmChannel}
end

Base.length(channels::YlmChannelSet) = length(channels.channel_data)
Base.firstindex(channels::YlmChannelSet) = firstindex(channels.channel_data)
Base.lastindex(channels::YlmChannelSet) = lastindex(channels.channel_data)
Base.getindex(channels::YlmChannelSet, index::Integer) = channels.channel_data[index]
Base.iterate(channels::YlmChannelSet, state...) = iterate(channels.channel_data, state...)

function Base.show(io::IO, channels::YlmChannelSet)
    print(io, "YlmChannelSet(lmax=", channels.lmax, ", nchannels=", length(channels), ")")
end

_maximum_l(channels::YlmChannelSet) = isempty(channels.channel_data) ? -1 : maximum(channel.l for channel in channels.channel_data)

"""
    ylm_channels(lmax::Int)

Build the ordered list of angular channels `(l, m)` with `0 <= l <= lmax` and
`m = -l:l`.

The ordering is by increasing `l`, and for fixed `l` by increasing `m`.
"""
function ylm_channels(lmax::Int)
    lmax >= 0 || throw(ArgumentError("ylm_channels requires lmax >= 0"))
    channel_data = YlmChannel[YlmChannel(l, m) for l in 0:lmax for m in -l:l]
    return YlmChannelSet(lmax, channel_data)
end

"""
    AtomicOneBodyOperators

Block-diagonal one-body atomic operator bundle built from a radial one-body
bundle and an explicit `YlmChannelSet`.
"""
struct AtomicOneBodyOperators
    channels::YlmChannelSet
    channel_ranges::Vector{UnitRange{Int}}
    overlap::Matrix{Float64}
    hamiltonian::Matrix{Float64}
end

function Base.show(io::IO, ops::AtomicOneBodyOperators)
    print(
        io,
        "AtomicOneBodyOperators(nchannels=",
        length(ops.channels),
        ", matrix_size=",
        size(ops.hamiltonian),
        ")",
    )
end

function _channel_index(channels::YlmChannelSet, channel::YlmChannel)
    index = findfirst(==(channel), channels.channel_data)
    index === nothing && throw(ArgumentError("channel $(channel) is not present in the supplied YlmChannelSet"))
    return index
end

"""
    channel_range(ops::AtomicOneBodyOperators, index::Int)
    channel_range(ops::AtomicOneBodyOperators, channel::YlmChannel)

Return the matrix-index range belonging to one angular channel.
"""
function channel_range(ops::AtomicOneBodyOperators, index::Int)
    index in eachindex(ops.channel_ranges) || throw(BoundsError(ops.channel_ranges, index))
    return ops.channel_ranges[index]
end

function channel_range(ops::AtomicOneBodyOperators, channel::YlmChannel)
    return ops.channel_ranges[_channel_index(ops.channels, channel)]
end

"""
    channel_hamiltonian(radial_ops::RadialAtomicOperators, channel::YlmChannel)

Return the radial one-body block for angular channel `channel`.
"""
function channel_hamiltonian(radial_ops::RadialAtomicOperators, channel::YlmChannel)
    return radial_ops.kinetic + radial_ops.nuclear + centrifugal(radial_ops, channel.l)
end

"""
    channel_hamiltonian(ops::AtomicOneBodyOperators, index::Int)
    channel_hamiltonian(ops::AtomicOneBodyOperators, channel::YlmChannel)

Return one channel block of the assembled atomic one-body Hamiltonian.
"""
function channel_hamiltonian(ops::AtomicOneBodyOperators, index::Int)
    block = channel_range(ops, index)
    return ops.hamiltonian[block, block]
end

function channel_hamiltonian(ops::AtomicOneBodyOperators, channel::YlmChannel)
    block = channel_range(ops, channel)
    return ops.hamiltonian[block, block]
end

"""
    channel_overlap(ops::AtomicOneBodyOperators, index::Int)
    channel_overlap(ops::AtomicOneBodyOperators, channel::YlmChannel)

Return one channel block of the assembled atomic overlap matrix.
"""
function channel_overlap(ops::AtomicOneBodyOperators, index::Int)
    block = channel_range(ops, index)
    return ops.overlap[block, block]
end

function channel_overlap(ops::AtomicOneBodyOperators, channel::YlmChannel)
    block = channel_range(ops, channel)
    return ops.overlap[block, block]
end

"""
    atomic_one_body_operators(radial_ops::RadialAtomicOperators, channels::YlmChannelSet)
    atomic_one_body_operators(radial_ops::RadialAtomicOperators; lmax::Int)

Build the first combined radial-angular one-body atomic operator bundle.

For the present central one-body atomic path, the Hamiltonian is block diagonal
in the angular channels. Each `(l, m)` channel gets the radial block

    T + V_nuc + C_l

with the common overlap block repeated on the diagonal.
"""
function atomic_one_body_operators(radial_ops::RadialAtomicOperators, channels::YlmChannelSet)
    radial_dim = size(radial_ops.overlap, 1)
    nchannels = length(channels)
    total_dim = radial_dim * nchannels
    maximum_l = _maximum_l(channels)
    maximum_l < length(radial_ops.centrifugal_data) ||
        throw(
            ArgumentError(
                "angular channels require radial centrifugal data through l = $(maximum_l), " *
                "but the supplied RadialAtomicOperators only contain data through l = $(length(radial_ops.centrifugal_data) - 1)",
            ),
        )

    overlap = zeros(Float64, total_dim, total_dim)
    hamiltonian = zeros(Float64, total_dim, total_dim)
    channel_ranges = Vector{UnitRange{Int}}(undef, nchannels)

    for i in eachindex(channels.channel_data)
        block = ((i - 1) * radial_dim + 1):(i * radial_dim)
        channel_ranges[i] = block
        overlap[block, block] = radial_ops.overlap
        hamiltonian[block, block] = channel_hamiltonian(radial_ops, channels[i])
    end

    return AtomicOneBodyOperators(channels, channel_ranges, overlap, hamiltonian)
end

function atomic_one_body_operators(radial_ops::RadialAtomicOperators; lmax::Int)
    return atomic_one_body_operators(radial_ops, ylm_channels(lmax))
end
