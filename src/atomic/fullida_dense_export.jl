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

function _atomic_mapping_manifest(mapping_value::IdentityMapping)
    return Dict{String,Any}(
        "manifest/source/mapping/type" => "IdentityMapping",
        "manifest/source/mapping/is_identity" => true,
        "manifest/source/mapping/a" => NaN,
        "manifest/source/mapping/c" => NaN,
        "manifest/source/mapping/s" => NaN,
        "manifest/source/mapping/tail_spacing" => NaN,
    )
end

function _atomic_mapping_manifest(mapping_value::AsinhMapping)
    return Dict{String,Any}(
        "manifest/source/mapping/type" => "AsinhMapping",
        "manifest/source/mapping/is_identity" => false,
        "manifest/source/mapping/a" => mapping_value.a,
        "manifest/source/mapping/c" => mapping_value.a * mapping_value.s,
        "manifest/source/mapping/s" => mapping_value.s,
        "manifest/source/mapping/tail_spacing" => mapping_value.tail_spacing,
    )
end

function _atomic_mapping_manifest(mapping_value::AbstractCoordinateMapping)
    return Dict{String,Any}(
        "manifest/source/mapping/type" => string(nameof(typeof(mapping_value))),
        "manifest/source/mapping/is_identity" => false,
        "manifest/source/mapping/a" => NaN,
        "manifest/source/mapping/c" => NaN,
        "manifest/source/mapping/s" => NaN,
        "manifest/source/mapping/tail_spacing" => NaN,
    )
end

function _atomic_source_manifest(ops::AtomicIDAOperators)
    source = ops.radial_operators.source_manifest
    spec = source.basis_spec
    channels = ops.one_body.channels
    xgaussian_alphas = Float64[xgaussian.alpha for xgaussian in spec.xgaussians]
    basis_description =
        spec.rmax === nothing ?
        "RadialBasisSpec(count=$(spec.count), family=$(spec.family_value.name))" :
        "RadialBasisSpec(rmax=$(spec.rmax), family=$(spec.family_value.name))"
    supplement_description =
        isempty(xgaussian_alphas) ?
        "none" :
        "XGaussian(alpha=$(join(string.(xgaussian_alphas), ",")))"

    values = Dict{String,Any}(
        "manifest/source/branch" => "atomic_ida",
        "manifest/source/model" => "radial_atomic_ida",
        "manifest/source/atomic_charge" => source.nuclear_charge,
        "manifest/source/basis_spec_type" => "RadialBasisSpec",
        "manifest/source/basis_family" => string(spec.family_value.name),
        "manifest/source/basis_description" => basis_description,
        "manifest/source/public_extent_kind" => spec.rmax === nothing ? "count" : "rmax",
        "manifest/source/public_rmax" => spec.rmax === nothing ? NaN : spec.rmax,
        "manifest/source/public_count" => spec.count === nothing ? 0 : spec.count,
        "manifest/source/has_public_rmax" => spec.rmax !== nothing,
        "manifest/source/has_public_count" => spec.count !== nothing,
        "manifest/source/reference_spacing" => spec.reference_spacing,
        "manifest/source/tails" => spec.tails,
        "manifest/source/odd_even_kmax" => spec.odd_even_kmax,
        "manifest/source/supplement_kind" => "xgaussian",
        "manifest/source/supplement_count" => length(spec.xgaussians),
        "manifest/source/supplement/xgaussian_alphas" => xgaussian_alphas,
        "manifest/source/supplement_description" => supplement_description,
        "manifest/source/radial_dimension" => size(ops.radial_operators.overlap, 1),
        "manifest/source/channel_count" => length(channels),
        "manifest/source/channel_lmax" => channels.lmax,
        "manifest/source/channel_l" => Int[channel.l for channel in channels],
        "manifest/source/channel_m" => Int[channel.m for channel in channels],
        "manifest/source/channel_convention" => "ylm_channels_increasing_l_then_m",
    )
    merge!(values, _atomic_mapping_manifest(spec.mapping_value))
    return values
end

function _atomic_export_meta_values(
    ops::AtomicIDAOperators;
    producer_entrypoint::AbstractString,
    interaction_model::AbstractString,
    interaction_detail::AbstractString,
    extra_values::AbstractDict{String,Any} = Dict{String,Any}(),
    meta = nothing,
)
    source = ops.radial_operators.source_manifest
    spec = source.basis_spec
    meta_values = _normalize_meta_dict(meta)
    merge!(
        meta_values,
        Dict{String,Any}(
            "Z" => source.nuclear_charge,
            "producer" => producer_entrypoint,
            "producer_type" => "AtomicIDAOperators",
            "source_branch" => "atomic_ida",
            "source_model" => "radial_atomic_ida",
            "interaction_model" => interaction_model,
            "interaction_detail" => interaction_detail,
            "nchannels" => length(ops.one_body.channels),
            "nradial" => size(ops.radial_operators.overlap, 1),
            "public_extent_kind" => spec.rmax === nothing ? "count" : "rmax",
            "public_rmax" => spec.rmax === nothing ? NaN : spec.rmax,
            "public_count" => spec.count === nothing ? 0 : spec.count,
            "has_public_rmax" => spec.rmax !== nothing,
            "has_public_count" => spec.count !== nothing,
            "manifest/producer/package" => "GaussletBases",
            "manifest/producer/version" => string(Base.pkgversion(@__MODULE__)),
            "manifest/producer/entrypoint" => producer_entrypoint,
            "manifest/producer/object_type" => "AtomicIDAOperators",
            "manifest/interaction/model" => interaction_model,
            "manifest/interaction/detail" => interaction_detail,
        ),
    )
    merge!(meta_values, _atomic_source_manifest(ops))
    merge!(meta_values, extra_values)
    return meta_values
end

function _atomic_shell_major_permutation(ops::AtomicIDAOperators)
    nchannels = length(ops.one_body.channels)
    radial_dim = size(ops.radial_operators.overlap, 1)
    return Int[(channel - 1) * radial_dim + radial for radial in 1:radial_dim for channel in 1:nchannels]
end

"""
    fullida_dense_payload(ops::AtomicIDAOperators; Vps=nothing, nelec=nothing, meta=nothing)

Build the current dense full-IDA export payload in memory without writing a
JLD2 file.

The returned named tuple contains:

- `payload`
- `bridge_meta`
- `meta_values`

and matches the data written by [`write_fullida_dense_jld2`](@ref).
"""
function fullida_dense_payload(
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
    shell_centers_r = Float64[Float64(value) for value in ops.radial_operators.shell_centers_r]

    H1 = Matrix{Float64}(ops.one_body.hamiltonian[perm, perm])
    Vee = atomic_ida_density_interaction_matrix(ops; ordering = perm)
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
        "order/shell_centers_r" => shell_centers_r,
        "order/basis_radius" => isempty(shell_centers_r) ? NaN : maximum(shell_centers_r),
        "order/permutation_from_in_memory" => perm,
    )

    meta_values = _atomic_export_meta_values(
        ops;
        producer_entrypoint = "GaussletBases.write_fullida_dense_jld2",
        interaction_model = "density_density_ida",
        interaction_detail = "two_index_ida",
        meta = meta,
    )

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
    data = fullida_dense_payload(ops; Vps = Vps, nelec = nelec, meta = meta)
    jldopen(path, "w") do file
        for (key, value) in data.payload
            file[key] = value
        end
        _write_prefixed_values!(file, "bridge", data.bridge_meta)
        _write_prefixed_values!(file, "meta", data.meta_values)
    end
    return path
end
