struct AtomicFixedRadialAngularSequenceLevel
    sequence_id::String
    level_id::String
    level_index::Int
    N_sph::Int
    radial_basis_id::String
    shell_ids::Vector{Int}
    shell_centers_r::Vector{Float64}
    shell_dimensions::Vector{Int}
    shell_orders::Vector{Int}
    profile::ShellLocalAngularProfile
    payload::AtomicInjectedAngularHFDMRGHFAdapter
end

struct AtomicFixedRadialAngularSequenceOverlapSidecar
    sequence_id::String
    source_level_index::Int
    target_level_index::Int
    source_N_sph::Int
    target_N_sph::Int
    source_profile_id::String
    target_profile_id::String
    source_labels::Vector{String}
    target_labels::Vector{String}
    overlap::Matrix{Float64}
    source_exact_count::Int
    source_mixed_count::Int
    target_exact_count::Int
    target_mixed_count::Int
    shell_independent::Bool
    diagnostics::NamedTuple
end

struct AtomicFixedRadialAngularSequence
    sequence_id::String
    radial_basis_id::String
    shell_ids::Vector{Int}
    shell_centers_r::Vector{Float64}
    N_sph_values::Vector{Int}
    profile_settings::NamedTuple
    levels::Vector{AtomicFixedRadialAngularSequenceLevel}
    adjacent_overlaps::Vector{AtomicFixedRadialAngularSequenceOverlapSidecar}
end

function Base.show(io::IO, level::AtomicFixedRadialAngularSequenceLevel)
    print(
        io,
        "AtomicFixedRadialAngularSequenceLevel(level_index=",
        level.level_index,
        ", N_sph=",
        level.N_sph,
        ", dim=",
        size(level.payload.hamiltonian, 1),
        ", profile_id=",
        repr(level.profile.profile_id),
        ")",
    )
end

function Base.show(io::IO, sidecar::AtomicFixedRadialAngularSequenceOverlapSidecar)
    print(
        io,
        "AtomicFixedRadialAngularSequenceOverlapSidecar(",
        sidecar.source_N_sph,
        "->",
        sidecar.target_N_sph,
        ", size=",
        size(sidecar.overlap),
        ", shell_independent=",
        sidecar.shell_independent,
        ")",
    )
end

function Base.show(io::IO, sequence::AtomicFixedRadialAngularSequence)
    print(
        io,
        "AtomicFixedRadialAngularSequence(levels=",
        length(sequence.levels),
        ", N_sph_values=",
        repr(sequence.N_sph_values),
        ", radial_basis_id=",
        repr(sequence.radial_basis_id),
        ")",
    )
end

function _fixed_radial_sequence_mapping_fields(mapping_value::IdentityMapping)
    return ("IdentityMapping", NaN, NaN, NaN, NaN)
end

function _fixed_radial_sequence_mapping_fields(mapping_value::AsinhMapping)
    return ("AsinhMapping", mapping_value.a, mapping_value.a * mapping_value.s, mapping_value.s, mapping_value.tail_spacing)
end

function _fixed_radial_sequence_mapping_fields(mapping_value::AbstractCoordinateMapping)
    return (string(nameof(typeof(mapping_value))), NaN, NaN, NaN, NaN)
end

function _fixed_radial_sequence_source_values(radial_ops::RadialAtomicOperators)
    source = radial_ops.source_manifest
    spec = source.basis_spec
    mapping_type, mapping_a, mapping_c, mapping_s, mapping_tail_spacing =
        _fixed_radial_sequence_mapping_fields(spec.mapping_value)
    xgaussian_alphas = Float64[xgaussian.alpha for xgaussian in spec.xgaussians]
    return Dict{String,Any}(
        "manifest/source/branch" => "atomic_fixed_radial_angular",
        "manifest/source/model" => "radial_atomic_substrate",
        "manifest/source/atomic_charge" => source.nuclear_charge,
        "manifest/source/basis_spec_type" => "RadialBasisSpec",
        "manifest/source/basis_family" => string(spec.family_value.name),
        "manifest/source/reference_spacing" => spec.reference_spacing,
        "manifest/source/public_extent_kind" => spec.rmax === nothing ? "count" : "rmax",
        "manifest/source/public_rmax" => spec.rmax === nothing ? NaN : spec.rmax,
        "manifest/source/public_count" => spec.count === nothing ? 0 : spec.count,
        "manifest/source/has_public_rmax" => spec.rmax !== nothing,
        "manifest/source/has_public_count" => spec.count !== nothing,
        "manifest/source/tails" => spec.tails,
        "manifest/source/odd_even_kmax" => spec.odd_even_kmax,
        "manifest/source/supplement_kind" => "xgaussian",
        "manifest/source/supplement_count" => length(spec.xgaussians),
        "manifest/source/supplement/xgaussian_alphas" => xgaussian_alphas,
        "manifest/source/mapping/type" => mapping_type,
        "manifest/source/mapping/a" => mapping_a,
        "manifest/source/mapping/c" => mapping_c,
        "manifest/source/mapping/s" => mapping_s,
        "manifest/source/mapping/tail_spacing" => mapping_tail_spacing,
        "manifest/source/radial_dimension" => size(radial_ops.overlap, 1),
        "manifest/source/centrifugal_lmax" => length(radial_ops.centrifugal_data) - 1,
        "manifest/source/multipole_Lmax" => length(radial_ops.multipole_data) - 1,
        "manifest/source/approximation" => string(nameof(typeof(radial_ops.approximation))),
    )
end

function _fixed_radial_basis_id(radial_ops::RadialAtomicOperators)
    source_values = _fixed_radial_sequence_source_values(radial_ops)
    lines = String[]
    for key in sort!(collect(keys(source_values)))
        value = source_values[key]
        if value isa AbstractVector
            push!(lines, string(key, "=", join(string.(value), ",")))
        else
            push!(lines, string(key, "=", value))
        end
    end
    append!(
        lines,
        (
            "shell_centers_checksum=" * _angular_digest_float_vector(radial_ops.shell_centers_r),
            "overlap_checksum=" * _angular_digest_float_matrix(radial_ops.overlap),
            "kinetic_checksum=" * _angular_digest_float_matrix(radial_ops.kinetic),
            "nuclear_checksum=" * _angular_digest_float_matrix(radial_ops.nuclear),
            "centrifugal_l0_checksum=" * _angular_digest_float_matrix(centrifugal(radial_ops, 0)),
            "multipole_l0_checksum=" * _angular_digest_float_matrix(multipole(radial_ops, 0)),
        ),
    )
    return bytes2hex(SHA.sha1(join(lines, "\n")))
end

function _fixed_radial_sequence_level_id(
    sequence_id::AbstractString,
    level_index::Int,
    N_sph::Int,
    profile_id::AbstractString,
)
    return bytes2hex(
        SHA.sha1(
            join(
                (
                    sequence_id,
                    string(level_index),
                    string(N_sph),
                    profile_id,
                ),
                "\n",
            ),
        ),
    )
end

function _fixed_radial_sequence_id(
    radial_basis_id::AbstractString,
    N_sph_values::AbstractVector{<:Integer},
    profile_settings::NamedTuple,
)
    lines = String[
        "radial_basis_id=" * radial_basis_id,
        "N_sph_values=" * join(string.(Int.(N_sph_values)), ","),
    ]
    for (key, value) in pairs(profile_settings)
        push!(lines, string(key, "=", value))
    end
    return bytes2hex(SHA.sha1(join(lines, "\n")))
end

function _build_atomic_fixed_radial_angular_sequence_level(
    sequence_id::AbstractString,
    radial_basis_id::AbstractString,
    radial_ops::RadialAtomicOperators,
    level_index::Int,
    N_sph::Int;
    beta::Real,
    l_inject::Union{Int,Symbol},
    tau::Real,
    whiten::Symbol,
    nelec::Int,
)
    shell_ids = collect(1:length(radial_ops.shell_centers_r))
    profile =
        shell_local_angular_profile(
            N_sph;
            beta = beta,
            l_inject = l_inject,
            tau = tau,
            whiten = whiten,
        )
    payload =
        build_atomic_injected_angular_hfdmrg_payload(
            radial_ops;
            shell_orders = fill(N_sph, length(shell_ids)),
            beta = beta,
            l_inject = l_inject,
            tau = tau,
            whiten = whiten,
            nelec = nelec,
        )
    shell_centers_r = Float64[Float64(value) for value in radial_ops.shell_centers_r]
    shell_dimensions = copy(payload.one_body.angular_assembly.shell_dimensions)
    shell_orders = copy(payload.one_body.angular_assembly.shell_orders)
    level_id = _fixed_radial_sequence_level_id(sequence_id, level_index, N_sph, profile.profile_id)
    return AtomicFixedRadialAngularSequenceLevel(
        String(sequence_id),
        level_id,
        level_index,
        N_sph,
        String(radial_basis_id),
        shell_ids,
        shell_centers_r,
        shell_dimensions,
        shell_orders,
        profile,
        payload,
    )
end

function _build_atomic_fixed_radial_overlap_sidecar(
    sequence_id::AbstractString,
    source::AtomicFixedRadialAngularSequenceLevel,
    target::AtomicFixedRadialAngularSequenceLevel,
)
    overlap =
        adjacent_shell_local_angular_profile_overlap(source.profile, target.profile)
    return AtomicFixedRadialAngularSequenceOverlapSidecar(
        String(sequence_id),
        source.level_index,
        target.level_index,
        source.N_sph,
        target.N_sph,
        overlap.source_profile_id,
        overlap.target_profile_id,
        copy(overlap.source_labels),
        copy(overlap.target_labels),
        Matrix{Float64}(overlap.overlap),
        overlap.source_exact_count,
        overlap.source_mixed_count,
        overlap.target_exact_count,
        overlap.target_mixed_count,
        overlap.shell_independent,
        overlap.diagnostics,
    )
end

"""
    build_atomic_fixed_radial_angular_sequence(
        radial_ops::RadialAtomicOperators,
        N_sph_values;
        beta=2.0,
        l_inject=:auto,
        tau=1e-12,
        whiten=:svd,
        nelec=round(Int, radial_ops.source_manifest.nuclear_charge),
    )

Build the fixed-radial-basis angular profile ladder for a user-specified list
of `N_sph` values.

The radial substrate is held fixed across every level. Each level reuses one
cached shell-local angular profile for its `N_sph` value and exports the dense
HF-facing bridge data independently.
"""
function build_atomic_fixed_radial_angular_sequence(
    radial_ops::RadialAtomicOperators,
    N_sph_values::AbstractVector{<:Integer};
    beta::Real = 2.0,
    l_inject::Union{Int,Symbol} = :auto,
    tau::Real = 1.0e-12,
    whiten::Symbol = :svd,
    nelec::Int = round(Int, radial_ops.source_manifest.nuclear_charge),
)
    isempty(N_sph_values) &&
        throw(ArgumentError("build_atomic_fixed_radial_angular_sequence requires at least one N_sph value"))
    resolved_N_sph = Int[Int(value) for value in N_sph_values]
    any(value -> value <= 0, resolved_N_sph) &&
        throw(ArgumentError("N_sph values must be positive"))

    radial_basis_id = _fixed_radial_basis_id(radial_ops)
    profile_settings = (
        beta = Float64(beta),
        l_inject_request = l_inject,
        tau = Float64(tau),
        whiten = whiten,
        gauge_version = _SHELL_LOCAL_ANGULAR_PROFILE_GAUGE_VERSION,
    )
    sequence_id = _fixed_radial_sequence_id(radial_basis_id, resolved_N_sph, profile_settings)
    levels = Vector{AtomicFixedRadialAngularSequenceLevel}(undef, length(resolved_N_sph))
    for i in eachindex(resolved_N_sph)
        levels[i] =
            _build_atomic_fixed_radial_angular_sequence_level(
                sequence_id,
                radial_basis_id,
                radial_ops,
                i,
                resolved_N_sph[i];
                beta = beta,
                l_inject = l_inject,
                tau = tau,
                whiten = whiten,
                nelec = nelec,
            )
    end

    adjacent_overlaps = AtomicFixedRadialAngularSequenceOverlapSidecar[]
    for i in 1:(length(levels) - 1)
        push!(
            adjacent_overlaps,
            _build_atomic_fixed_radial_overlap_sidecar(sequence_id, levels[i], levels[i + 1]),
        )
    end

    return AtomicFixedRadialAngularSequence(
        sequence_id,
        radial_basis_id,
        collect(1:length(radial_ops.shell_centers_r)),
        Float64[Float64(value) for value in radial_ops.shell_centers_r],
        resolved_N_sph,
        profile_settings,
        levels,
        adjacent_overlaps,
    )
end

function build_atomic_fixed_radial_angular_sequence(
    basis::RadialBasis,
    grid::RadialQuadratureGrid,
    N_sph_values::AbstractVector{<:Integer};
    Z::Real,
    lmax::Int = 0,
    approximation::AbstractDiagonalApproximation = IntegralDiagonal(),
    kwargs...,
)
    radial_ops = atomic_operators(
        basis,
        grid;
        Z = Z,
        lmax = lmax,
        approximation = approximation,
    )
    return build_atomic_fixed_radial_angular_sequence(radial_ops, N_sph_values; kwargs...)
end

function _fixed_radial_level_meta_values(
    level::AtomicFixedRadialAngularSequenceLevel;
    meta = nothing,
)
    meta_values = _normalize_meta_dict(meta)
    merge!(
        meta_values,
        Dict{String,Any}(
            "sequence_id" => level.sequence_id,
            "level_id" => level.level_id,
            "level_index" => level.level_index,
            "N_sph" => level.N_sph,
            "radial_basis_id" => level.radial_basis_id,
            "angular_profile_id" => level.profile.profile_id,
            "manifest/producer/package" => "GaussletBases",
            "manifest/producer/version" => string(Base.pkgversion(@__MODULE__)),
            "manifest/producer/entrypoint" => "GaussletBases.write_atomic_fixed_radial_angular_level_jld2",
            "manifest/producer/object_type" => "AtomicFixedRadialAngularSequenceLevel",
            "manifest/sequence/kind" => "fixed_radial_n_sph_ladder",
            "manifest/sequence/shell_independent_adjacent_overlap" => true,
        ),
    )
    merge!(meta_values, _fixed_radial_sequence_source_values(level.payload.one_body.radial_operators))
    return meta_values
end

"""
    atomic_fixed_radial_angular_level_dense_payload(level; meta=nothing)

Build the current dense per-level bridge payload for one fixed-radial angular
sequence level without writing a file.
"""
function atomic_fixed_radial_angular_level_dense_payload(
    level::AtomicFixedRadialAngularSequenceLevel;
    meta = nothing,
)
    profile = level.profile
    bridge_meta = Dict{String,Any}(
        "format" => "angular_fixed_radial_dense_v1",
        "version" => 1,
        "sequence_id" => level.sequence_id,
        "level_id" => level.level_id,
        "level_index" => level.level_index,
        "N_sph" => level.N_sph,
        "radial_basis_id" => level.radial_basis_id,
        "angular_profile_id" => profile.profile_id,
        "gauge_version" => string(profile.key.gauge_version),
        "shell_count" => length(level.shell_ids),
        "shell_ids" => copy(level.shell_ids),
        "shell_centers_r" => copy(level.shell_centers_r),
        "shell_dimensions" => copy(level.shell_dimensions),
        "shell_orders" => copy(level.shell_orders),
        "within_shell_label_count" => length(profile.labels),
        "within_shell_exact_count" => profile.basis.injected_count,
        "within_shell_mixed_count" => profile.basis.final_count - profile.basis.injected_count,
        "route" => string(level.payload.route),
        "solver_mode" => string(level.payload.solver_mode),
    )
    payload = Dict{String,Any}(
        "H1" => Matrix{Float64}(level.payload.hamiltonian),
        "Vee" => Matrix{Float64}(level.payload.interaction),
        "psiup0" => Matrix{Float64}(level.payload.psiup0),
        "psidn0" => Matrix{Float64}(level.payload.psidn0),
        "shell_ids" => copy(level.shell_ids),
        "shell_centers_r" => copy(level.shell_centers_r),
        "shell_dimensions" => copy(level.shell_dimensions),
        "shell_orders" => copy(level.shell_orders),
        "within_shell_labels" => copy(profile.labels),
        "within_shell_block_kinds" => String[string(kind) for kind in profile.block_kinds],
        "within_shell_exact_labels" => copy(profile.exact_labels),
        "within_shell_mixed_labels" => copy(profile.mixed_labels),
    )
    return (
        payload = payload,
        bridge_meta = bridge_meta,
        meta_values = _fixed_radial_level_meta_values(level; meta = meta),
    )
end

"""
    write_atomic_fixed_radial_angular_level_jld2(path, level; meta=nothing)

Write one native dense per-level artifact for the fixed-radial increasing-`N_sph`
sequence line.

This writer stores the current level Hamiltonian, interaction, seeds, stable
radial shell metadata, and stable within-shell labels from the cached angular
profile. It is the native producer-side level export for this sequence family.
"""
function write_atomic_fixed_radial_angular_level_jld2(
    path::AbstractString,
    level::AtomicFixedRadialAngularSequenceLevel;
    meta = nothing,
)
    data = atomic_fixed_radial_angular_level_dense_payload(level; meta = meta)
    jldopen(path, "w") do file
        for (key, value) in data.payload
            file[key] = value
        end
        _write_prefixed_values!(file, "bridge", data.bridge_meta)
        _write_prefixed_values!(file, "meta", data.meta_values)
    end
    return path
end

function _fixed_radial_overlap_meta_values(
    sidecar::AtomicFixedRadialAngularSequenceOverlapSidecar;
    meta = nothing,
)
    meta_values = _normalize_meta_dict(meta)
    merge!(
        meta_values,
        Dict{String,Any}(
            "sequence_id" => sidecar.sequence_id,
            "source_level_index" => sidecar.source_level_index,
            "target_level_index" => sidecar.target_level_index,
            "source_N_sph" => sidecar.source_N_sph,
            "target_N_sph" => sidecar.target_N_sph,
            "source_profile_id" => sidecar.source_profile_id,
            "target_profile_id" => sidecar.target_profile_id,
            "manifest/producer/package" => "GaussletBases",
            "manifest/producer/version" => string(Base.pkgversion(@__MODULE__)),
            "manifest/producer/entrypoint" => "GaussletBases.write_atomic_fixed_radial_angular_overlap_sidecar_jld2",
            "manifest/producer/object_type" => "AtomicFixedRadialAngularSequenceOverlapSidecar",
        ),
    )
    return meta_values
end

const _LEGACY_DMRGATOM_CENTER_POLICY_NAME =
    "representative_dominant_prototype_direction_v1"

function _legacy_dmrgatom_within_shell_dominant_prototype_indices(
    profile::ShellLocalAngularProfile,
)
    coefficients = profile.basis.prototype_coefficients
    size(coefficients, 1) == profile.basis.point_set.cardinality || throw(
        DimensionMismatch(
            "prototype coefficient rows must match point-set cardinality for legacy dmrgatom representative centers",
        ),
    )
    norbitals = size(coefficients, 2)
    indices = Vector{Int}(undef, norbitals)
    magnitudes = Vector{Float64}(undef, norbitals)
    for orbital in 1:norbitals
        dominant = _dominant_column_index(view(coefficients, :, orbital))
        magnitude = abs(coefficients[dominant, orbital])
        magnitude > 0.0 || throw(
            ArgumentError(
                "legacy dmrgatom representative centers require a nonzero dominant prototype coefficient for every orbital",
            ),
        )
        indices[orbital] = dominant
        magnitudes[orbital] = magnitude
    end
    return indices, magnitudes
end

function _legacy_dmrgatom_basis_centers(
    level::AtomicFixedRadialAngularSequenceLevel,
)
    profile = level.profile
    per_shell_dim = profile.basis.final_count
    all(dim -> dim == per_shell_dim, level.shell_dimensions) || throw(
        ArgumentError(
            "legacy dmrgatom export currently requires shell-independent profile dimensions matching the cached profile final count",
        ),
    )

    dominant_indices, dominant_magnitudes =
        _legacy_dmrgatom_within_shell_dominant_prototype_indices(profile)
    directions = profile.basis.point_set.coordinates[dominant_indices, :]
    nshells = length(level.shell_dimensions)
    norbitals = sum(level.shell_dimensions)
    basis_centers = Matrix{Float64}(undef, norbitals, 3)
    row = 1
    for shell in 1:nshells
        radius = level.shell_centers_r[shell]
        for orbital in 1:per_shell_dim
            basis_centers[row, :] .= radius .* view(directions, orbital, :)
            row += 1
        end
    end
    return basis_centers, dominant_indices, dominant_magnitudes
end

function _legacy_dmrgatom_meta_values(
    level::AtomicFixedRadialAngularSequenceLevel;
    meta = nothing,
)
    nuclear_charge = level.payload.one_body.radial_operators.source_manifest.nuclear_charge
    nuclear_charge_integer = round(Int, nuclear_charge)
    isapprox(nuclear_charge, nuclear_charge_integer; atol = 1.0e-12, rtol = 0.0) || throw(
        ArgumentError(
            "legacy dmrgatom export requires an integer-valued atomic charge for meta/Z; got $(nuclear_charge)",
        ),
    )
    meta_values = _fixed_radial_level_meta_values(level; meta = meta)
    merge!(
        meta_values,
        Dict{String,Any}(
            "Z" => nuclear_charge_integer,
            "producer" => "GaussletBases.write_atomic_fixed_radial_legacy_dmrgatom_jld2",
            "producer_type" => "AtomicFixedRadialAngularSequenceLevel",
            "manifest/producer/entrypoint" =>
                "GaussletBases.write_atomic_fixed_radial_legacy_dmrgatom_jld2",
            "manifest/producer/object_type" => "AtomicFixedRadialAngularSequenceLevel",
            "manifest/consumer/family" => "legacy_dmrgatom",
            "manifest/consumer/target" => "GaussletModules/bin/dmrgatom.jl",
            "manifest/geometry/basis_centers_kind" =>
                "representative_orbital_centers_for_legacy_ordering",
            "manifest/geometry/center_policy_name" =>
                _LEGACY_DMRGATOM_CENTER_POLICY_NAME,
        ),
    )
    return meta_values
end

"""
    atomic_fixed_radial_angular_overlap_sidecar_payload(sidecar; meta=nothing)

Build the adjacent shell-local profile-overlap sidecar for one neighboring
`N_sph[k] -> N_sph[k+1]` pair without writing a file.
"""
function atomic_fixed_radial_angular_overlap_sidecar_payload(
    sidecar::AtomicFixedRadialAngularSequenceOverlapSidecar;
    meta = nothing,
)
    bridge_meta = Dict{String,Any}(
        "format" => "angular_fixed_radial_profile_overlap_v1",
        "version" => 1,
        "sequence_id" => sidecar.sequence_id,
        "source_level_index" => sidecar.source_level_index,
        "target_level_index" => sidecar.target_level_index,
        "source_N_sph" => sidecar.source_N_sph,
        "target_N_sph" => sidecar.target_N_sph,
        "source_profile_id" => sidecar.source_profile_id,
        "target_profile_id" => sidecar.target_profile_id,
        "source_exact_count" => sidecar.source_exact_count,
        "source_mixed_count" => sidecar.source_mixed_count,
        "target_exact_count" => sidecar.target_exact_count,
        "target_mixed_count" => sidecar.target_mixed_count,
        "shell_independent" => sidecar.shell_independent,
        "min_singular_value" => sidecar.diagnostics.min_singular_value,
        "max_singular_value" => sidecar.diagnostics.max_singular_value,
        "overlap_inf_norm" => sidecar.diagnostics.inf_norm,
        "overlap_frobenius_norm" => sidecar.diagnostics.frobenius_norm,
        "exact_block_inf_norm" => sidecar.diagnostics.exact_block_inf_norm,
    )
    payload = Dict{String,Any}(
        "overlap" => Matrix{Float64}(sidecar.overlap),
        "source_labels" => copy(sidecar.source_labels),
        "target_labels" => copy(sidecar.target_labels),
    )
    return (
        payload = payload,
        bridge_meta = bridge_meta,
        meta_values = _fixed_radial_overlap_meta_values(sidecar; meta = meta),
    )
end

"""
    write_atomic_fixed_radial_angular_overlap_sidecar_jld2(path, sidecar; meta=nothing)

Write the native adjacent-overlap sidecar for one neighboring
`N_sph[k] -> N_sph[k+1]` pair in the fixed-radial sequence line.

This writer stores the compact shell-local cross-overlap, stable source/target
labels, and overlap diagnostics/provenance needed by external continuation
consumers.
"""
function write_atomic_fixed_radial_angular_overlap_sidecar_jld2(
    path::AbstractString,
    sidecar::AtomicFixedRadialAngularSequenceOverlapSidecar;
    meta = nothing,
)
    data = atomic_fixed_radial_angular_overlap_sidecar_payload(sidecar; meta = meta)
    jldopen(path, "w") do file
        for (key, value) in data.payload
            file[key] = value
        end
        _write_prefixed_values!(file, "bridge", data.bridge_meta)
        _write_prefixed_values!(file, "meta", data.meta_values)
    end
    return path
end

"""
    atomic_fixed_radial_legacy_dmrgatom_payload(level; meta=nothing)

Build the narrow legacy-direct dense payload for the `dmrgatom.jl` consumer
family from one fixed-radial angular sequence level without writing a file.

The required top-level datasets are:

- `H1`
- `Vee`
- `basis_centers`
- `dims_per_shell`

with `meta/Z` and additional provenance carried under `meta/` and `bridge/`.
The exported `basis_centers` use a representative-center policy for legacy
ordering heuristics, not exact physical orbital centers.
"""
function atomic_fixed_radial_legacy_dmrgatom_payload(
    level::AtomicFixedRadialAngularSequenceLevel;
    meta = nothing,
)
    profile = level.profile
    basis_centers, dominant_indices, dominant_magnitudes =
        _legacy_dmrgatom_basis_centers(level)
    bridge_meta = Dict{String,Any}(
        "format" => "legacy_dmrgatom_dense_v1",
        "version" => 1,
        "sequence_id" => level.sequence_id,
        "level_id" => level.level_id,
        "level_index" => level.level_index,
        "N_sph" => level.N_sph,
        "radial_basis_id" => level.radial_basis_id,
        "angular_profile_id" => profile.profile_id,
        "gauge_version" => string(profile.key.gauge_version),
        "shell_count" => length(level.shell_ids),
        "shell_ids" => copy(level.shell_ids),
        "shell_centers_r" => copy(level.shell_centers_r),
        "dims_per_shell" => copy(level.shell_dimensions),
        "within_shell_labels" => copy(profile.labels),
        "within_shell_block_kinds" => String[string(kind) for kind in profile.block_kinds],
        "within_shell_exact_labels" => copy(profile.exact_labels),
        "within_shell_mixed_labels" => copy(profile.mixed_labels),
        "within_shell_dominant_prototype_indices" => copy(dominant_indices),
        "within_shell_dominant_prototype_magnitudes" => copy(dominant_magnitudes),
        "basis_centers_kind" => "representative_orbital_centers_for_legacy_ordering",
        "basis_centers_center_policy_name" => _LEGACY_DMRGATOM_CENTER_POLICY_NAME,
        "basis_centers_are_representative" => true,
        "route" => string(level.payload.route),
        "solver_mode" => string(level.payload.solver_mode),
    )
    payload = Dict{String,Any}(
        "H1" => Matrix{Float64}(level.payload.hamiltonian),
        "Vee" => Matrix{Float64}(level.payload.interaction),
        "basis_centers" => basis_centers,
        "dims_per_shell" => copy(level.shell_dimensions),
    )
    return (
        payload = payload,
        bridge_meta = bridge_meta,
        meta_values = _legacy_dmrgatom_meta_values(level; meta = meta),
    )
end

"""
    write_atomic_fixed_radial_legacy_dmrgatom_jld2(path, level; meta=nothing)

Write the narrow legacy-direct dense export for the `dmrgatom.jl` consumer
family from one fixed-radial angular sequence level.

This writer is intentionally consumer-specific. It does not attempt to solve
the broader legacy chemistry-style JLD2 compatibility problem.
"""
function write_atomic_fixed_radial_legacy_dmrgatom_jld2(
    path::AbstractString,
    level::AtomicFixedRadialAngularSequenceLevel;
    meta = nothing,
)
    data = atomic_fixed_radial_legacy_dmrgatom_payload(level; meta = meta)
    jldopen(path, "w") do file
        for (key, value) in data.payload
            file[key] = value
        end
        _write_prefixed_values!(file, "bridge", data.bridge_meta)
        _write_prefixed_values!(file, "meta", data.meta_values)
    end
    return path
end
