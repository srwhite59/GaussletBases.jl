_m_sign_rank(m::Integer) = m == 0 ? 0 : (m < 0 ? 1 : 2)

function _mzigzag_then_l_channel_permutation(channels::YlmChannelSet)
    idxs = collect(eachindex(channels.channel_data))
    return sortperm(idxs, by = i -> (
        abs(channels[i].m),
        _m_sign_rank(channels[i].m),
        channels[i].l,
        i,
    ))
end

function _l0_desc_mzigzag_channel_permutation(channels::YlmChannelSet)
    idxs = collect(eachindex(channels.channel_data))
    return sortperm(idxs, by = i -> (
        channels[i].l == 0 ? 0 : 1,
        channels[i].l == 0 ? 0 : -channels[i].l,
        abs(channels[i].m),
        _m_sign_rank(channels[i].m),
        i,
    ))
end

function _atomic_sliced_permutation(ops::AtomicIDAOperators)
    radial_dim = size(ops.radial_operators.overlap, 1)
    channel_perm = _l0_desc_mzigzag_channel_permutation(ops.one_body.channels)
    orbital_perm = Int[
        (channel - 1) * radial_dim + radial for radial in 1:radial_dim for channel in channel_perm
    ]
    return orbital_perm, channel_perm
end

function _atomic_hamv6_permutation(ops::AtomicIDAOperators)
    radial_dim = size(ops.radial_operators.overlap, 1)
    channel_perm = _mzigzag_then_l_channel_permutation(ops.one_body.channels)
    orbital_perm = Int[
        (channel - 1) * radial_dim + radial for radial in 1:radial_dim for channel in channel_perm
    ]
    return orbital_perm, channel_perm
end

function _dense_block_to_coo(block::AbstractMatrix; threshold::Real = 0.0)
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    for row in axes(block, 1), col in axes(block, 2)
        value = Float64(block[row, col])
        abs(value) > threshold || continue
        push!(rows, row)
        push!(cols, col)
        push!(vals, value)
    end
    return (rows = rows, cols = cols, vals = vals)
end

function _density_density_pair_block_to_coo(block::AbstractMatrix; threshold::Real = 0.0)
    dn, dm = size(block)
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    for p in 1:dn, q in 1:dm
        value = Float64(block[p, q])
        abs(value) > threshold || continue
        pair_index = (p - 1) * dm + q
        push!(rows, pair_index)
        push!(cols, pair_index)
        push!(vals, value)
    end
    return (rows = rows, cols = cols, vals = vals)
end

function _coo_blocks_to_dense(blocks, dims::AbstractVector{<:Integer})
    offs = _shell_offsets(dims)
    total_dim = sum(Int.(dims))
    matrix = zeros(Float64, total_dim, total_dim)
    for n in eachindex(dims), m in eachindex(dims)
        row_offset = offs[n] - 1
        col_offset = offs[m] - 1
        block = blocks[n][m]
        for k in eachindex(block.vals)
            matrix[row_offset + Int(block.rows[k]), col_offset + Int(block.cols[k])] = block.vals[k]
        end
    end
    return matrix
end

function _pairdiag_blocks_to_density_matrix(blocks, dims::AbstractVector{<:Integer})
    offs = _shell_offsets(dims)
    total_dim = sum(Int.(dims))
    matrix = zeros(Float64, total_dim, total_dim)
    for n in eachindex(dims), m in eachindex(dims)
        dm = Int(dims[m])
        row_offset = offs[n] - 1
        col_offset = offs[m] - 1
        block = blocks[n][m]
        for k in eachindex(block.vals)
            row = Int(block.rows[k])
            col = Int(block.cols[k])
            row == col || throw(ArgumentError("density-density pair block must be diagonal in pair space"))
            p = (row - 1) ÷ dm + 1
            q = (row - 1) % dm + 1
            matrix[row_offset + p, col_offset + q] = block.vals[k]
        end
    end
    return matrix
end

function _atomic_onebody_component_matrices(ops::AtomicIDAOperators)
    total_dim = size(ops.one_body.hamiltonian, 1)
    Tblocks = zeros(Float64, total_dim, total_dim)
    Vnucblocks = zeros(Float64, total_dim, total_dim)
    for (channel_index, channel) in enumerate(ops.one_body.channels)
        block = channel_range(ops.one_body, channel_index)
        Tblocks[block, block] = ops.radial_operators.kinetic + centrifugal(ops.radial_operators, channel.l)
        Vnucblocks[block, block] = ops.radial_operators.nuclear
    end
    return (
        H1 = Matrix{Float64}(ops.one_body.hamiltonian),
        T = Tblocks,
        Vnuc = Vnucblocks,
    )
end

function _atomic_sliced_payload_impl(
    ops::AtomicIDAOperators;
    orbital_perm::AbstractVector{<:Integer},
    channel_perm::AbstractVector{<:Integer},
    within_slice::AbstractString,
    ordering_description::AbstractString,
    orbital_ordering::AbstractString,
    producer_entrypoint::AbstractString,
    consumer_shape::AbstractString,
    format_name::AbstractString,
    threshold::Real = 0.0,
    nelec::Union{Nothing,Int} = nothing,
    meta = nothing,
)
    ordered_channels = ops.one_body.channels.channel_data[channel_perm]
    nchannels = length(ordered_channels)
    radial_dim = size(ops.radial_operators.overlap, 1)
    norb = nchannels * radial_dim
    shell_centers_r = Float64[Float64(value) for value in ops.radial_operators.shell_centers_r]
    dims = fill(nchannels, radial_dim)
    offs = _shell_offsets(dims)

    components = _atomic_onebody_component_matrices(ops)
    H1 = Matrix{Float64}(components.H1[orbital_perm, orbital_perm])
    T = Matrix{Float64}(components.T[orbital_perm, orbital_perm])
    Vnuc = Matrix{Float64}(components.Vnuc[orbital_perm, orbital_perm])
    Vee = atomic_ida_density_interaction_matrix(ops; ordering = collect(Int, orbital_perm))

    H1blocks = [
        [
            _dense_block_to_coo(
                view(H1, offs[n]:(offs[n + 1] - 1), offs[m]:(offs[m + 1] - 1));
                threshold = threshold,
            )
            for m in 1:radial_dim
        ]
        for n in 1:radial_dim
    ]
    Tblocks = [
        [
            _dense_block_to_coo(
                view(T, offs[n]:(offs[n + 1] - 1), offs[m]:(offs[m + 1] - 1));
                threshold = threshold,
            )
            for m in 1:radial_dim
        ]
        for n in 1:radial_dim
    ]
    Vnucblocks = [
        [
            _dense_block_to_coo(
                view(Vnuc, offs[n]:(offs[n + 1] - 1), offs[m]:(offs[m + 1] - 1));
                threshold = threshold,
            )
            for m in 1:radial_dim
        ]
        for n in 1:radial_dim
    ]
    Vblocks = [
        [
            _density_density_pair_block_to_coo(
                view(Vee, offs[n]:(offs[n + 1] - 1), offs[m]:(offs[m + 1] - 1));
                threshold = threshold,
            )
            for m in 1:radial_dim
        ]
        for n in 1:radial_dim
    ]

    m_template = Int[channel.m for channel in ordered_channels]
    l_template = Int[channel.l for channel in ordered_channels]
    m_by_slice = [copy(m_template) for _ in 1:radial_dim]
    l_by_slice = [copy(l_template) for _ in 1:radial_dim]
    labels_by_slice = [
        String["r=$(slice),l=$(channel.l),m=$(channel.m)" for channel in ordered_channels]
        for slice in 1:radial_dim
    ]

    layout_values = Dict{String,Any}(
        "nslices" => radial_dim,
        "dims" => dims,
        "offs" => offs,
        "slice_coord" => shell_centers_r,
        "slice_index" => collect(1:radial_dim),
    )
    basis_values = Dict{String,Any}(
        "m_by_slice" => m_by_slice,
        "l_by_slice" => l_by_slice,
        "m_flat" => vcat(m_by_slice...),
        "l_flat" => vcat(l_by_slice...),
        "labels_by_slice" => labels_by_slice,
    )
    ordering_values = Dict{String,Any}(
        "major" => "slice_major",
        "slice_meaning" => "radial_shell",
        "within_slice" => within_slice,
        "description" => ordering_description,
        "permutation_from_in_memory" => collect(Int, orbital_perm),
    )
    onebody_values = Dict{String,Any}(
        "stored" => "coo",
        "is_hermitian" => true,
        "decomposition" => "H1blocks = Tblocks + Vnucblocks, with Tblocks including the centrifugal contribution",
        "H1blocks" => H1blocks,
        "Tblocks" => Tblocks,
        "Vnucblocks" => Vnucblocks,
    )
    twobody_values = Dict{String,Any}(
        "stored" => "coo_all",
        "convention" => "density_density_pairdiag_v1",
        "symmetry" => "pair_diagonal_density_density",
        "description" => "pair-diagonal COO blocks encoding the current two-index IDA density-density interaction model",
        "Vblocks" => Vblocks,
    )

    meta_values = _atomic_export_meta_values(
        ops;
        producer_entrypoint = producer_entrypoint,
        interaction_model = "density_density_ida",
        interaction_detail = "two_index_ida_pair_diagonal",
        extra_values = Dict{String,Any}(
            "format" => format_name,
            "consumer_shape" => consumer_shape,
            "onebody_model" => "H1 = Tblocks + Vnucblocks, with centrifugal included in Tblocks",
            "twobody_encoding" => "pair_diagonal_density_density",
            "slice_kind" => "radial_shell",
            "slice_coord_kind" => "physical_radial_center",
            "slice_index_kind" => "radial_index",
            "orbital_ordering" => orbital_ordering,
            "norb" => norb,
            "nelec" => something(nelec, 0),
            "has_nelec" => nelec !== nothing,
            "permutation_from_in_memory" => collect(Int, orbital_perm),
        ),
        meta = meta,
    )

    return (
        layout_values = layout_values,
        basis_values = basis_values,
        ordering_values = ordering_values,
        onebody_values = onebody_values,
        twobody_values = twobody_values,
        meta_values = meta_values,
    )
end

"""
    sliced_ham_payload(ops::AtomicIDAOperators; nelec=nothing, meta=nothing, threshold=0.0)

Build the current sliced/block atomic Hamiltonian export payload in memory
without writing a JLD2 file.

The returned named tuple contains:

- `layout_values`
- `basis_values`
- `ordering_values`
- `onebody_values`
- `twobody_values`
- `meta_values`

and matches the data written by [`write_sliced_ham_jld2`](@ref).
"""
function sliced_ham_payload(
    ops::AtomicIDAOperators;
    nelec::Union{Nothing,Int} = nothing,
    meta = nothing,
    threshold::Real = 0.0,
)
    orbital_perm, channel_perm = _atomic_sliced_permutation(ops)
    return _atomic_sliced_payload_impl(
        ops;
        orbital_perm = orbital_perm,
        channel_perm = channel_perm,
        within_slice = "l0_desc_mzigzag",
        ordering_description = "slice-major by radial index; within slice l=0 first, then descending l with m-zigzag order",
        orbital_ordering = "slice_major_by_radial_index_then_l0_desc_mzigzag",
        producer_entrypoint = "GaussletBases.write_sliced_ham_jld2",
        consumer_shape = "slicedmrgutils.HamIO",
        format_name = "atomic_ida_sliced_v1",
        threshold = threshold,
        nelec = nelec,
        meta = meta,
    )
end

"""
    atomic_hamv6_payload(ops::AtomicIDAOperators; nelec=nothing, meta=nothing, threshold=0.0)

Build the current atomic IDA sliced/block payload in the explicit HamV6 /
`HamIO.jl` compatibility ordering used by the present downstream bridge stack.

This keeps the current atomic density-density IDA model and grouped COO block
storage, but adapts the within-slice orbital order to `mzigzag_then_l`.
"""
function atomic_hamv6_payload(
    ops::AtomicIDAOperators;
    nelec::Union{Nothing,Int} = nothing,
    meta = nothing,
    threshold::Real = 0.0,
)
    orbital_perm, channel_perm = _atomic_hamv6_permutation(ops)
    return _atomic_sliced_payload_impl(
        ops;
        orbital_perm = orbital_perm,
        channel_perm = channel_perm,
        within_slice = "mzigzag_then_l",
        ordering_description = "slice-major; within slice m=0,-1,+1,-2,+2,...; for each m, l increasing",
        orbital_ordering = "slice_major_by_radial_index_then_mzigzag_then_l",
        producer_entrypoint = "GaussletBases.write_atomic_hamv6_jld2",
        consumer_shape = "slicedmrgutils.HamIO/HamV6",
        format_name = "atomic_hamv6_v1",
        threshold = threshold,
        nelec = nelec,
        meta = meta,
    )
end

function _angular_benchmark_exact_hamv6_payload_impl(
    ops::AtomicIDAOperators;
    benchmark_type::AbstractString,
    exact_common_lmax::Int,
    shell_orders::AbstractVector{<:Integer},
    nelec::Union{Nothing,Int} = nothing,
    meta = nothing,
    threshold::Real = 0.0,
    producer_entrypoint::AbstractString,
)
    data = atomic_hamv6_payload(ops; nelec = nelec, meta = meta, threshold = threshold)
    meta_values = copy(data.meta_values)
    merge!(
        meta_values,
        Dict{String,Any}(
            "producer" => producer_entrypoint,
            "producer_type" => benchmark_type,
            "manifest/producer/entrypoint" => producer_entrypoint,
            "manifest/producer/object_type" => benchmark_type,
            "angular_bridge_kind" => "exact_common_low_l_reference",
            "angular_bridge_scope" => "exact_common_low_l_reference_only",
            "angular_bridge_consumer_language" => "slicedmrgutils.HamIO/HamV6",
            "angular_exact_common_lmax" => exact_common_lmax,
            "angular_shell_orders" => Int.(collect(shell_orders)),
            "angular_shell_count" => length(shell_orders),
            "angular_underlying_export" => "atomic_hamv6_payload",
            "angular_underlying_object_type" => "AtomicIDAOperators",
        ),
    )
    return (
        layout_values = data.layout_values,
        basis_values = data.basis_values,
        ordering_values = data.ordering_values,
        onebody_values = data.onebody_values,
        twobody_values = data.twobody_values,
        meta_values = meta_values,
    )
end

"""
    angular_benchmark_exact_hamv6_payload(
        benchmark::AtomicInjectedAngularHFStyleBenchmark;
        nelec=nothing,
        meta=nothing,
        threshold=0.0,
    )

Build an explicit HamV6 / `HamIO.jl` compatibility payload for the exact common
low-`l` reference embedded in the angular benchmark line.

This is intentionally narrow. It does **not** export the full mixed shell-local
angular basis, because that basis is not yet representable in the current
consumer language that assumes definite per-orbital `l,m` labels. Instead, it
exports the honest exact common reference already carried by the benchmark.
"""
function angular_benchmark_exact_hamv6_payload(
    benchmark::AtomicInjectedAngularHFStyleBenchmark;
    nelec::Union{Nothing,Int} = nothing,
    meta = nothing,
    threshold::Real = 0.0,
)
    return _angular_benchmark_exact_hamv6_payload_impl(
        benchmark.exact_ida_reference;
        benchmark_type = "AtomicInjectedAngularHFStyleBenchmark",
        exact_common_lmax = benchmark.one_body.exact_common_lmax,
        shell_orders = benchmark.one_body.angular_assembly.shell_orders,
        nelec = nelec,
        meta = meta,
        threshold = threshold,
        producer_entrypoint = "GaussletBases.write_angular_benchmark_exact_hamv6_jld2",
    )
end

"""
    angular_benchmark_exact_hamv6_payload(
        benchmark::AtomicInjectedAngularSmallEDBenchmark;
        nelec=nothing,
        meta=nothing,
        threshold=0.0,
    )

Build the same exact-common-reference HamV6 payload starting from the small-ED
benchmark wrapper.
"""
function angular_benchmark_exact_hamv6_payload(
    benchmark::AtomicInjectedAngularSmallEDBenchmark;
    nelec::Union{Nothing,Int} = nothing,
    meta = nothing,
    threshold::Real = 0.0,
)
    return _angular_benchmark_exact_hamv6_payload_impl(
        benchmark.hf_style.exact_ida_reference;
        benchmark_type = "AtomicInjectedAngularSmallEDBenchmark",
        exact_common_lmax = benchmark.hf_style.one_body.exact_common_lmax,
        shell_orders = benchmark.hf_style.one_body.angular_assembly.shell_orders,
        nelec = nelec,
        meta = meta,
        threshold = threshold,
        producer_entrypoint = "GaussletBases.write_angular_benchmark_exact_hamv6_jld2",
    )
end

"""
    write_atomic_hamv6_jld2(path, ops::AtomicIDAOperators; nelec=nothing, meta=nothing, threshold=0.0)

Write the explicit HamV6 / `HamIO.jl` compatibility export for the current
atomic density-density IDA model.
"""
function write_atomic_hamv6_jld2(
    path::AbstractString,
    ops::AtomicIDAOperators;
    nelec::Union{Nothing,Int} = nothing,
    meta = nothing,
    threshold::Real = 0.0,
)
    data = atomic_hamv6_payload(ops; nelec = nelec, meta = meta, threshold = threshold)
    jldopen(path, "w") do file
        _write_prefixed_values!(file, "layout", data.layout_values)
        _write_prefixed_values!(file, "basis", data.basis_values)
        _write_prefixed_values!(file, "ordering", data.ordering_values)
        _write_prefixed_values!(file, "onebody", data.onebody_values)
        _write_prefixed_values!(file, "twobody", data.twobody_values)
        _write_prefixed_values!(file, "meta", data.meta_values)
    end
    return path
end

"""
    write_angular_benchmark_exact_hamv6_jld2(
        path,
        benchmark::Union{
            AtomicInjectedAngularHFStyleBenchmark,
            AtomicInjectedAngularSmallEDBenchmark,
        };
        nelec=nothing,
        meta=nothing,
        threshold=0.0,
    )

Write the narrow angular bridge artifact that reuses the proven HamV6 grouped
shape for the benchmark's exact common low-`l` reference only.
"""
function write_angular_benchmark_exact_hamv6_jld2(
    path::AbstractString,
    benchmark::Union{
        AtomicInjectedAngularHFStyleBenchmark,
        AtomicInjectedAngularSmallEDBenchmark,
    };
    nelec::Union{Nothing,Int} = nothing,
    meta = nothing,
    threshold::Real = 0.0,
)
    data = angular_benchmark_exact_hamv6_payload(
        benchmark;
        nelec = nelec,
        meta = meta,
        threshold = threshold,
    )
    jldopen(path, "w") do file
        _write_prefixed_values!(file, "layout", data.layout_values)
        _write_prefixed_values!(file, "basis", data.basis_values)
        _write_prefixed_values!(file, "ordering", data.ordering_values)
        _write_prefixed_values!(file, "onebody", data.onebody_values)
        _write_prefixed_values!(file, "twobody", data.twobody_values)
        _write_prefixed_values!(file, "meta", data.meta_values)
    end
    return path
end

"""
    write_sliced_ham_jld2(path, ops::AtomicIDAOperators; nelec=nothing, meta=nothing, threshold=0.0)

Write the grouped sliced/block Hamiltonian export for the current atomic IDA model.

The written file follows the existing `HamIO.jl` grouped structure:

- `layout/*`
- `basis/*`
- `ordering/*`
- `onebody/*`
- `twobody/*`
- `meta/*`

The present export remains honest about the current model:

- atomic-only
- density-density / two-index IDA interaction
- not a full four-index Coulomb Hamiltonian
"""
function write_sliced_ham_jld2(
    path::AbstractString,
    ops::AtomicIDAOperators;
    nelec::Union{Nothing,Int} = nothing,
    meta = nothing,
    threshold::Real = 0.0,
)
    data = sliced_ham_payload(ops; nelec = nelec, meta = meta, threshold = threshold)
    jldopen(path, "w") do file
        _write_prefixed_values!(file, "layout", data.layout_values)
        _write_prefixed_values!(file, "basis", data.basis_values)
        _write_prefixed_values!(file, "ordering", data.ordering_values)
        _write_prefixed_values!(file, "onebody", data.onebody_values)
        _write_prefixed_values!(file, "twobody", data.twobody_values)
        _write_prefixed_values!(file, "meta", data.meta_values)
    end
    return path
end
