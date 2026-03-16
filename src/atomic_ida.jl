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

function _legal_angular_triple(l1::Int, l2::Int, L::Int)
    return abs(l1 - l2) <= L && l1 + l2 >= L && iseven(l1 + l2 + L)
end

function _factorial_table(nmax::Int)
    return [BigFloat(factorial(big(n))) for n in 0:nmax]
end

@inline _fact(table::AbstractVector{BigFloat}, n::Int) = table[n + 1]

function _wigner3j(
    fact_table::AbstractVector{BigFloat},
    j1::Int,
    j2::Int,
    j3::Int,
    m1::Int,
    m2::Int,
    m3::Int,
)
    m1 + m2 + m3 == 0 || return 0.0
    abs(m1) <= j1 || return 0.0
    abs(m2) <= j2 || return 0.0
    abs(m3) <= j3 || return 0.0
    _legal_angular_triple(j1, j2, j3) || return 0.0

    delta = sqrt(
        _fact(fact_table, j1 + j2 - j3) *
        _fact(fact_table, j1 - j2 + j3) *
        _fact(fact_table, -j1 + j2 + j3) /
        _fact(fact_table, j1 + j2 + j3 + 1),
    )

    pref = (isodd(j1 - j2 - m3) ? -one(BigFloat) : one(BigFloat)) * delta
    pref *= sqrt(
        _fact(fact_table, j1 + m1) *
        _fact(fact_table, j1 - m1) *
        _fact(fact_table, j2 + m2) *
        _fact(fact_table, j2 - m2) *
        _fact(fact_table, j3 + m3) *
        _fact(fact_table, j3 - m3),
    )

    tmin = max(0, j2 - j3 - m1, j1 - j3 + m2)
    tmax = min(j1 + j2 - j3, j1 - m1, j2 + m2)
    tmin <= tmax || return 0.0

    total = zero(BigFloat)
    for t in tmin:tmax
        denom =
            _fact(fact_table, t) *
            _fact(fact_table, j1 + j2 - j3 - t) *
            _fact(fact_table, j1 - m1 - t) *
            _fact(fact_table, j2 + m2 - t) *
            _fact(fact_table, j3 - j2 + m1 + t) *
            _fact(fact_table, j3 - j1 - m2 + t)
        total += (isodd(t) ? -one(BigFloat) : one(BigFloat)) / denom
    end

    return Float64(pref * total)
end

function _complex_gaunt(
    fact_table::AbstractVector{BigFloat},
    left::YlmChannel,
    L::Int,
    M::Int,
    right::YlmChannel,
)
    _legal_angular_triple(left.l, right.l, L) || return 0.0
    -left.m + M + right.m == 0 || return 0.0

    pref = (isodd(left.m) ? -1.0 : 1.0) *
           sqrt(((2 * left.l + 1) * (2 * L + 1) * (2 * right.l + 1)) / (4 * pi))
    return pref *
           _wigner3j(fact_table, left.l, L, right.l, 0, 0, 0) *
           _wigner3j(fact_table, left.l, L, right.l, -left.m, M, right.m)
end

function _build_gaunt_tensors(channels::YlmChannelSet)
    nchannels = length(channels)
    Lmax = 2 * _maximum_l(channels)
    fact_table = _factorial_table(max(1, 4 * max(channels.lmax, Lmax) + 1))
    gaunt_data = Vector{Array{Float64,3}}(undef, Lmax + 1)

    for L in 0:Lmax
        tensor = zeros(Float64, nchannels, nchannels, 2 * L + 1)
        for alpha in eachindex(channels.channel_data)
            left = channels[alpha]
            for alphap in eachindex(channels.channel_data)
                right = channels[alphap]
                M = left.m - right.m
                abs(M) <= L || continue
                tensor[alpha, alphap, M + L + 1] = _complex_gaunt(fact_table, left, L, M, right)
            end
        end
        gaunt_data[L + 1] = tensor
    end

    return gaunt_data
end

function _build_angular_kernels(gaunt_data::Vector{Array{Float64,3}})
    angular_kernel_data = Vector{Array{Float64,4}}(undef, length(gaunt_data))

    for L in 0:(length(gaunt_data) - 1)
        tensor = gaunt_data[L + 1]
        nchannels = size(tensor, 1)
        kernel = zeros(Float64, nchannels, nchannels, nchannels, nchannels)
        pref = 4 * pi / (2 * L + 1)

        for alpha in 1:nchannels, alphap in 1:nchannels, beta in 1:nchannels, betap in 1:nchannels
            total = 0.0
            for M in -L:L
                phase = isodd(M) ? -1.0 : 1.0
                total += tensor[alpha, alphap, M + L + 1] * phase * tensor[beta, betap, -M + L + 1]
            end
            kernel[alpha, alphap, beta, betap] = pref * total
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
    gaunt_data = _build_gaunt_tensors(channels)
    angular_kernel_data = _build_angular_kernels(gaunt_data)
    orbital_data = _atomic_orbitals(channels, size(radial_ops.overlap, 1))
    return AtomicIDAOperators(one_body, radial_ops, gaunt_data, angular_kernel_data, orbital_data)
end

function atomic_ida_operators(radial_ops::RadialAtomicOperators; lmax::Int)
    return atomic_ida_operators(radial_ops, ylm_channels(lmax))
end
