using LinearAlgebra
using Printf
using GaussletBases

const DEFAULTS = Dict(
    "family" => "ns7",
    "mode" => "compare",
    "Z" => "2.0",
    "rmax" => "10.0",
    "tail_spacing" => "10.0",
    "backend" => "numerical_reference",
)

const FAMILY_SPECS = Dict(
    "ns7" => (nside = 7, d = 0.15),
    "ns9" => (nside = 9, d = 0.10),
)

function _parse_args(args)
    values = copy(DEFAULTS)
    for arg in args
        occursin('=', arg) || continue
        key, value = split(arg, '='; limit = 2)
        values[key] = value
    end
    return values
end

function _resolve_count_and_basis(mapping; target_rmax::Real)
    count = 2 * ceil(Int, uofx(mapping, Float64(target_rmax))) + 1
    basis = build_basis(MappedUniformBasisSpec(:G10;
        count = count,
        mapping = mapping,
        reference_spacing = 1.0,
    ))
    while maximum(abs.(Float64.(centers(basis)))) + 1.0e-12 < Float64(target_rmax)
        count += 2
        basis = build_basis(MappedUniformBasisSpec(:G10;
            count = count,
            mapping = mapping,
            reference_spacing = 1.0,
        ))
    end
    return count, basis
end

function _build_full_parent_sequence(family_label::String; Z::Float64, rmax::Float64, tail_spacing::Float64, backend::Symbol)
    spec = FAMILY_SPECS[family_label]
    basis_count, basis = _resolve_count_and_basis(
        white_lindsey_atomic_mapping(Z = Z, d = spec.d, tail_spacing = tail_spacing);
        target_rmax = rmax,
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    sequence = build_one_center_atomic_full_parent_shell_sequence(
        basis;
        exponents = expansion.exponents,
        gausslet_backend = backend,
        refinement_levels = 0,
        nside = spec.nside,
    )
    return (
        family = family_label,
        d = spec.d,
        nside = spec.nside,
        count = basis_count,
        basis = basis,
        expansion = expansion,
        sequence = sequence,
    )
end

function _reference_weight_aware_pair_terms(
    bundles,
    support_states::AbstractVector{<:NTuple{3,Int}},
    support_coefficients::AbstractMatrix{<:Real},
)
    pgdg_x = GaussletBases._nested_axis_pgdg(bundles, :x)
    pgdg_y = GaussletBases._nested_axis_pgdg(bundles, :y)
    pgdg_z = GaussletBases._nested_axis_pgdg(bundles, :z)
    nterms = size(pgdg_x.pair_factor_terms, 1)
    nfixed = size(support_coefficients, 2)
    parent_weight_outer_x = pgdg_x.weights * transpose(pgdg_x.weights)
    parent_weight_outer_y = pgdg_y.weights * transpose(pgdg_y.weights)
    parent_weight_outer_z = pgdg_z.weights * transpose(pgdg_z.weights)
    support_weights = GaussletBases._nested_support_weights(
        support_states,
        pgdg_x.weights,
        pgdg_y.weights,
        pgdg_z.weights,
    )
    fixed_weights = vec(transpose(support_coefficients) * support_weights)
    fixed_weight_outer = fixed_weights * transpose(fixed_weights)
    pair_terms = zeros(Float64, nterms, nfixed, nfixed)
    for term in 1:nterms
        raw_x = @view(pgdg_x.pair_factor_terms[term, :, :]) .* parent_weight_outer_x
        raw_y = @view(pgdg_y.pair_factor_terms[term, :, :]) .* parent_weight_outer_y
        raw_z = @view(pgdg_z.pair_factor_terms[term, :, :]) .* parent_weight_outer_z
        raw_support = GaussletBases._nested_support_product_matrix(
            support_states,
            raw_x,
            raw_y,
            raw_z,
        )
        raw_contracted = transpose(support_coefficients) * raw_support * support_coefficients
        matrix = Matrix{Float64}(raw_contracted ./ fixed_weight_outer)
        @views pair_terms[term, :, :] .= 0.5 .* (matrix .+ transpose(matrix))
    end
    return (weights = fixed_weights, pair_terms = pair_terms)
end

function _reference_shell_packet(
    bundles,
    coefficient_matrix::AbstractMatrix{<:Real},
    support_indices::AbstractVector{Int},
)
    dims = GaussletBases._nested_axis_lengths(bundles)
    pgdg_x = GaussletBases._nested_axis_pgdg(bundles, :x)
    pgdg_y = GaussletBases._nested_axis_pgdg(bundles, :y)
    pgdg_z = GaussletBases._nested_axis_pgdg(bundles, :z)
    support_states = [GaussletBases._cartesian_unflat_index(index, dims) for index in support_indices]
    support_coefficients = Matrix{Float64}(coefficient_matrix[support_indices, :])
    overlap_support = GaussletBases._nested_support_product_matrix(
        support_states,
        pgdg_x.overlap,
        pgdg_y.overlap,
        pgdg_z.overlap,
    )
    kinetic_support = GaussletBases._nested_sum_of_support_products(
        support_states,
        (
            (pgdg_x.kinetic, pgdg_y.overlap, pgdg_z.overlap),
            (pgdg_x.overlap, pgdg_y.kinetic, pgdg_z.overlap),
            (pgdg_x.overlap, pgdg_y.overlap, pgdg_z.kinetic),
        ),
    )
    position_x_support = GaussletBases._nested_support_product_matrix(
        support_states,
        pgdg_x.position,
        pgdg_y.overlap,
        pgdg_z.overlap,
    )
    position_y_support = GaussletBases._nested_support_product_matrix(
        support_states,
        pgdg_x.overlap,
        pgdg_y.position,
        pgdg_z.overlap,
    )
    position_z_support = GaussletBases._nested_support_product_matrix(
        support_states,
        pgdg_x.overlap,
        pgdg_y.overlap,
        pgdg_z.position,
    )
    x2_x_support = GaussletBases._nested_support_product_matrix(
        support_states,
        pgdg_x.x2,
        pgdg_y.overlap,
        pgdg_z.overlap,
    )
    x2_y_support = GaussletBases._nested_support_product_matrix(
        support_states,
        pgdg_x.overlap,
        pgdg_y.x2,
        pgdg_z.overlap,
    )
    x2_z_support = GaussletBases._nested_support_product_matrix(
        support_states,
        pgdg_x.overlap,
        pgdg_y.overlap,
        pgdg_z.x2,
    )
    nshell = size(coefficient_matrix, 2)
    nterms = size(pgdg_x.gaussian_factor_terms, 1)
    pair_data = _reference_weight_aware_pair_terms(bundles, support_states, support_coefficients)
    gaussian_terms = zeros(Float64, nterms, nshell, nshell)
    for term in 1:nterms
        factor_support = GaussletBases._nested_support_product_matrix(
            support_states,
            @view(pgdg_x.gaussian_factor_terms[term, :, :]),
            @view(pgdg_y.gaussian_factor_terms[term, :, :]),
            @view(pgdg_z.gaussian_factor_terms[term, :, :]),
        )
        @views gaussian_terms[term, :, :] .= transpose(support_coefficients) * factor_support * support_coefficients
    end
    return GaussletBases._CartesianNestedShellPacket3D(
        transpose(support_coefficients) * overlap_support * support_coefficients,
        transpose(support_coefficients) * kinetic_support * support_coefficients,
        transpose(support_coefficients) * position_x_support * support_coefficients,
        transpose(support_coefficients) * position_y_support * support_coefficients,
        transpose(support_coefficients) * position_z_support * support_coefficients,
        transpose(support_coefficients) * x2_x_support * support_coefficients,
        transpose(support_coefficients) * x2_y_support * support_coefficients,
        transpose(support_coefficients) * x2_z_support * support_coefficients,
        pair_data.weights,
        gaussian_terms,
        pair_data.pair_terms,
    )
end

function _metric(name::String, current, reference)
    diff = current .- reference
    return (
        name = name,
        inf = maximum(abs.(diff)),
        fro = norm(vec(diff)),
        ref_inf = maximum(abs.(reference)),
    )
end

function _compare_packets(family_label::String; Z::Float64, rmax::Float64, tail_spacing::Float64, backend::Symbol)
    data = _build_full_parent_sequence(family_label; Z = Z, rmax = rmax, tail_spacing = tail_spacing, backend = backend)
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        data.basis;
        exponents = data.expansion.exponents,
        backend = backend,
        refinement_levels = 0,
    )
    bundles = GaussletBases._CartesianNestedAxisBundles3D(bundle, bundle, bundle)
    sequence = data.sequence

    t0 = time_ns()
    reference_packet = _reference_shell_packet(bundles, sequence.coefficient_matrix, sequence.support_indices)
    reference_seconds = (time_ns() - t0) / 1.0e9

    t1 = time_ns()
    current_packet = GaussletBases._nested_shell_packet(bundles, sequence.coefficient_matrix, sequence.support_indices).packet
    current_seconds = (time_ns() - t1) / 1.0e9

    metrics = (
        _metric("overlap", current_packet.overlap, reference_packet.overlap),
        _metric("kinetic", current_packet.kinetic, reference_packet.kinetic),
        _metric("position_x", current_packet.position_x, reference_packet.position_x),
        _metric("position_y", current_packet.position_y, reference_packet.position_y),
        _metric("position_z", current_packet.position_z, reference_packet.position_z),
        _metric("x2_x", current_packet.x2_x, reference_packet.x2_x),
        _metric("x2_y", current_packet.x2_y, reference_packet.x2_y),
        _metric("x2_z", current_packet.x2_z, reference_packet.x2_z),
        _metric("weights", current_packet.weights, reference_packet.weights),
        _metric("gaussian_terms", current_packet.gaussian_terms, reference_packet.gaussian_terms),
        _metric("pair_terms", current_packet.pair_terms, reference_packet.pair_terms),
    )

    println("COMPARE")
    println("  family = ", family_label)
    println("  count = ", data.count)
    println("  nside = ", data.nside)
    println("  nshells = ", length(sequence.shell_layers))
    println("  fixed_dim = ", size(sequence.coefficient_matrix, 2))
    println("  support_count = ", length(sequence.support_indices))
    println("  reference_packet_seconds = ", @sprintf("%.6f", reference_seconds))
    println("  current_packet_seconds = ", @sprintf("%.6f", current_seconds))
    for metric in metrics
        println("  METRIC")
        println("    name = ", metric.name)
        println("    inf = ", @sprintf("%.3e", metric.inf))
        println("    fro = ", @sprintf("%.3e", metric.fro))
        println("    ref_inf = ", @sprintf("%.3e", metric.ref_inf))
    end
end

function _current_fixed_block(family_label::String; Z::Float64, rmax::Float64, tail_spacing::Float64, backend::Symbol)
    spec = FAMILY_SPECS[family_label]
    count, basis = _resolve_count_and_basis(
        white_lindsey_atomic_mapping(Z = Z, d = spec.d, tail_spacing = tail_spacing);
        target_rmax = rmax,
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    t0 = time_ns()
    fixed_block = one_center_atomic_full_parent_fixed_block(
        basis;
        exponents = expansion.exponents,
        gausslet_backend = backend,
        refinement_levels = 0,
        nside = spec.nside,
    )
    elapsed = (time_ns() - t0) / 1.0e9
    println("CURRENT_FIXED_BLOCK")
    println("  family = ", family_label)
    println("  count = ", count)
    println("  nside = ", spec.nside)
    println("  fixed_dim = ", size(fixed_block.overlap, 1))
    println("  support_count = ", length(fixed_block.support_indices))
    println("  wall_seconds = ", @sprintf("%.6f", elapsed))
    println("  overlap_error = ", @sprintf("%.3e", norm(fixed_block.overlap - I, Inf)))
end

function main(args)
    config = _parse_args(args)
    family_label = config["family"]
    mode = config["mode"]
    Z = parse(Float64, config["Z"])
    rmax = parse(Float64, config["rmax"])
    tail_spacing = parse(Float64, config["tail_spacing"])
    backend = Symbol(config["backend"])
    haskey(FAMILY_SPECS, family_label) || throw(ArgumentError("unknown family=$family_label"))

    if mode == "compare"
        _compare_packets(family_label; Z = Z, rmax = rmax, tail_spacing = tail_spacing, backend = backend)
    elseif mode == "current_fixed_block"
        _current_fixed_block(family_label; Z = Z, rmax = rmax, tail_spacing = tail_spacing, backend = backend)
    else
        throw(ArgumentError("unknown mode=$mode"))
    end
end

main(ARGS)
