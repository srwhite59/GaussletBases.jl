function _shell_offsets(dims::AbstractVector{<:Integer})
    offsets = Vector{Int}(undef, length(dims) + 1)
    offsets[1] = 1
    for i in eachindex(dims)
        offsets[i + 1] = offsets[i] + Int(dims[i])
    end
    return offsets
end

function _shell_index(dims::AbstractVector{<:Integer})
    index = Vector{Int}(undef, sum(Int.(dims)))
    p = 1
    for shell in eachindex(dims)
        d = Int(dims[shell])
        if d > 0
            index[p:(p + d - 1)] .= shell
            p += d
        end
    end
    return index
end

function _normalize_meta_dict(meta)
    meta === nothing && return Dict{String,Any}()
    if meta isa NamedTuple
        return Dict{String,Any}(String(key) => value for (key, value) in pairs(meta))
    elseif meta isa AbstractDict
        return Dict{String,Any}(string(key) => value for (key, value) in pairs(meta))
    else
        throw(ArgumentError("meta must be nothing, a NamedTuple, or an AbstractDict"))
    end
end

function _write_prefixed_values!(file, prefix::AbstractString, values::AbstractDict{String,Any})
    for (key, value) in values
        file[string(prefix, "/", key)] = value
    end
    return nothing
end

function _atomic_shell_major_permutation(ops::AtomicIDAOperators)
    nchannels = length(ops.one_body.channels)
    radial_dim = size(ops.radial_operators.overlap, 1)
    return Int[(channel - 1) * radial_dim + radial for radial in 1:radial_dim for channel in 1:nchannels]
end

function _atomic_fullida_dense_payload(
    ops::AtomicIDAOperators;
    Vps = nothing,
    nelec::Union{Nothing,Int} = nothing,
    meta = nothing,
)
    perm = _atomic_shell_major_permutation(ops)
    orbital_data = ops.orbital_data[perm]
    nchannels = length(ops.one_body.channels)
    radial_dim = size(ops.radial_operators.overlap, 1)
    norb = length(orbital_data)

    H1 = Matrix{Float64}(ops.one_body.hamiltonian[perm, perm])
    Vee = _ida_density_interaction_matrix(ops, orbital_data)
    basis_centers = zeros(Float64, norb, 3)
    dims_per_shell = fill(nchannels, radial_dim)
    orders = collect(1:radial_dim)
    shell_offsets = _shell_offsets(dims_per_shell)
    shell_index = _shell_index(dims_per_shell)

    bridge_meta = Dict{String,Any}(
        "format" => "fullida_dense_v1",
        "version" => 1,
        "site_type" => "Electron",
        "interaction_model" => "density_density",
        "model_detail" => "two_index_ida",
        "source_branch" => "atomic_ida",
        "onebody_key" => "H1",
        "interaction_key" => "Vee",
        "optional_interaction_key" => "Vps",
        "has_optional_interaction" => Vps !== nothing,
        "norb" => norb,
        "nelec" => something(nelec, 0),
        "has_nelec" => nelec !== nothing,
        "order/convention" => "shell_major_by_radial_index",
        "order/within_shell" => "ylm_channel_order",
        "order/dims_per_shell" => dims_per_shell,
        "order/shell_offsets" => shell_offsets,
        "order/shell_index" => shell_index,
        "order/orders_per_shell" => orders,
        "order/basis_centers_kind" => "origin_only_atomic_orbitals",
        "order/shell_centers_r" => fill(NaN, radial_dim),
        "order/basis_radius" => NaN,
        "order/permutation_from_in_memory" => perm,
    )

    meta_values = Dict{String,Any}(
        "producer" => "GaussletBases.write_fullida_dense_jld2",
        "producer_type" => "AtomicIDAOperators",
        "source_branch" => "atomic_ida",
        "interaction_model" => "density_density_ida",
        "nchannels" => nchannels,
        "nradial" => radial_dim,
    )
    merge!(meta_values, _normalize_meta_dict(meta))

    payload = Dict{String,Any}(
        "H1" => H1,
        "Vee" => Vee,
        "dims_per_shell" => dims_per_shell,
        "orders" => orders,
        "basis_centers" => basis_centers,
    )
    Vps !== nothing && (payload["Vps"] = Matrix{Float64}(Vps))

    return (payload = payload, bridge_meta = bridge_meta, meta_values = meta_values)
end

"""
    write_fullida_dense_jld2(path, ops::AtomicIDAOperators; Vps=nothing, nelec=nothing, meta=nothing)

Write the current dense full-IDA bridge file for downstream solver consumers.

The written file follows the existing `fullida_dense_v1` bridge shape:

- `H1`
- `Vee` (and optional `Vps`)
- `dims_per_shell`
- `orders`
- `basis_centers`
- `bridge/...` metadata
- `meta/...` producer metadata

The present export is honest about the current model:

- one-body `H1`
- dense two-index IDA `Vee`
- density-density interaction model

It is not a full four-index Coulomb Hamiltonian export.
"""
function write_fullida_dense_jld2(
    path::AbstractString,
    ops::AtomicIDAOperators;
    Vps = nothing,
    nelec::Union{Nothing,Int} = nothing,
    meta = nothing,
)
    data = _atomic_fullida_dense_payload(ops; Vps = Vps, nelec = nelec, meta = meta)
    jldopen(path, "w") do file
        for (key, value) in data.payload
            file[key] = value
        end
        _write_prefixed_values!(file, "bridge", data.bridge_meta)
        _write_prefixed_values!(file, "meta", data.meta_values)
    end
    return path
end
