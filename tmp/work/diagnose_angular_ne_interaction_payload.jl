using GaussletBases
using LinearAlgebra
using Printf

function parse_kv_args(args)
    values = Dict{String,String}()
    for arg in args
        parts = split(arg, "="; limit = 2)
        length(parts) == 2 || error("expected key=value argument, got $(repr(arg))")
        values[parts[1]] = parts[2]
    end
    return values
end

function get_arg(values, key, default)
    return haskey(values, key) ? parse(typeof(default), values[key]) : default
end

function get_int_list_arg(values, key, default::Vector{Int})
    haskey(values, key) || return default
    raw = strip(values[key])
    isempty(raw) && return default
    return [parse(Int, strip(part)) for part in split(raw, ",") if !isempty(strip(part))]
end

function paper_style_radial_ops(; Z::Float64, lmax::Int, s::Float64, rmax::Float64)
    rb = build_basis(
        RadialBasisSpec(
            :G10;
            rmax = rmax,
            mapping = AsinhMapping(c = s / (2 * Z), s = s),
            reference_spacing = 1.0,
            tails = 6,
            odd_even_kmax = 6,
            xgaussian_count = 2,
            rmax_count_policy = :legacy_strict_trim,
        ),
    )
    grid = radial_quadrature(rb)
    return atomic_operators(rb, grid; Z = Z, lmax = lmax)
end

function shell_angular_weights(moment_blocks::Dict{Int,Matrix{Float64}})
    haskey(moment_blocks, 0) ||
        throw(ArgumentError("shell-local moment blocks must contain the monopole L=0 block"))
    weights = sqrt(4 * pi) .* vec(moment_blocks[0][1, :])
    any(weight -> abs(weight) <= 1.0e-12, weights) &&
        throw(ArgumentError("shell-local angular weights are too small for density-density assembly"))
    return weights
end

function scaled_shell_moment_block(moment_blocks::Dict{Int,Matrix{Float64}}, L::Int)
    haskey(moment_blocks, L) || return nothing
    weights = shell_angular_weights(moment_blocks)
    return moment_blocks[L] ./ reshape(weights, 1, :)
end

function assemble_interaction_from_blocks(
    radial_ops::RadialAtomicOperators,
    assembly::AtomicShellLocalInjectedAngularAssembly,
    shell_blocks::Vector{Dict{Int,Matrix{Float64}}},
)
    nshells = length(assembly.shells)
    shell_ranges = [
        assembly.shell_offsets[i]:(assembly.shell_offsets[i] + assembly.shell_dimensions[i] - 1)
        for i in eachindex(assembly.shell_dimensions)
    ]
    interaction = zeros(Float64, size(assembly.overlap))
    Lmax = length(radial_ops.multipole_data) - 1

    for L in 0:Lmax
        prefactor = 4 * pi / (2 * L + 1)
        radial_block = multipole(radial_ops, L)
        for a in 1:nshells
            ia = shell_ranges[a]
            left = scaled_shell_moment_block(shell_blocks[a], L)
            left === nothing && continue
            for b in 1:nshells
                ib = shell_ranges[b]
                right = scaled_shell_moment_block(shell_blocks[b], L)
                right === nothing && continue
                interaction[ia, ib] .+= prefactor .* radial_block[a, b] .* (transpose(left) * right)
            end
        end
    end

    return 0.5 .* (interaction .+ transpose(interaction))
end

function restricted_energy_decomposition(
    one_body::AbstractMatrix{<:Real},
    interaction::AbstractMatrix{<:Real},
    density::AbstractMatrix{<:Real},
)
    H = Matrix{Float64}(one_body)
    V = Matrix{Float64}(interaction)
    rho = 0.5 .* (Matrix{Float64}(density) .+ transpose(Matrix{Float64}(density)))
    occupations = vec(diag(rho))
    one_body_energy = 2.0 * tr(rho * H)
    direct_energy = 2.0 * dot(occupations, V * occupations)
    exchange_energy = dot(vec(rho), vec(V .* rho))
    return (
        one_body = one_body_energy,
        direct = direct_energy,
        exchange = exchange_energy,
        total = one_body_energy + direct_energy - exchange_energy,
    )
end

function row_summary(mat::AbstractMatrix{<:Real}, n::Int)
    count = min(n, size(mat, 1))
    return [
        (
            row = i,
            diag = Float64(mat[i, i]),
            rowsum = Float64(sum(view(mat, i, :))),
            rowabs = Float64(sum(abs, view(mat, i, :))),
        )
        for i in 1:count
    ]
end

function low_shell_report(assembly::AtomicShellLocalInjectedAngularAssembly, count::Int = 6)
    nshells = min(count, length(assembly.shells))
    return [
        (
            shell = i,
            radius = assembly.shell_radii[i],
            order = assembly.shell_orders[i],
            l_inject = assembly.shells[i].l_inject,
            Ny = assembly.shells[i].injected_count,
            dim = assembly.shells[i].final_count,
            bare_moment_lmax = maximum(keys(assembly.shell_moment_blocks[i])),
            kinetic_lmax = maximum(keys(assembly.shell_kinetic_moment_blocks[i])),
            kinetic_lexpand = assembly.shell_kinetic_lexpand[i],
            interaction_lmax = maximum(keys(assembly.shell_interaction_moment_blocks[i])),
            interaction_lexpand = assembly.shell_interaction_lexpand[i],
            weight_min = minimum(shell_angular_weights(assembly.shell_interaction_moment_blocks[i])),
            weight_max = maximum(shell_angular_weights(assembly.shell_interaction_moment_blocks[i])),
        )
        for i in 1:nshells
    ]
end

function print_decomposition(io::IO, label::AbstractString, scf_result, one_body, interaction)
    decomp = restricted_energy_decomposition(one_body, interaction, scf_result.density)
    @printf(io, "  %s energy = %.12f\n", label, scf_result.energy)
    @printf(
        io,
        "    one_body = %.12f  direct = %.12f  exchange = %.12f\n",
        decomp.one_body,
        decomp.direct,
        decomp.exchange,
    )
end

function report_case(io::IO, order::Int; Z::Float64, lmax::Int, top::Int)
    radial_ops = paper_style_radial_ops(; Z = Z, lmax = lmax, s = 0.2, rmax = 30.0)
    one_body = build_atomic_injected_angular_one_body_benchmark(
        radial_ops;
        shell_orders = fill(order, length(radial_ops.shell_centers_r)),
    )
    assembly = one_body.angular_assembly
    truncated_interaction =
        assemble_interaction_from_blocks(radial_ops, assembly, assembly.shell_moment_blocks)
    assembled_interaction =
        GaussletBases._assemble_atomic_injected_angular_interaction(radial_ops, assembly)
    helper_interaction = assemble_interaction_from_blocks(
        radial_ops,
        assembly,
        assembly.shell_interaction_moment_blocks,
    )
    projected_truncated =
        one_body.exact_transform * truncated_interaction * transpose(one_body.exact_transform)
    projected_assembled =
        one_body.exact_transform * assembled_interaction * transpose(one_body.exact_transform)
    exact_ida = atomic_ida_operators(radial_ops; lmax = one_body.exact_common_lmax)
    exact_interaction =
        atomic_ida_density_interaction_matrix(exact_ida; ordering = :shell_major)

    exact_overlap = one_body.exact_overlap
    exact_h1 = one_body.exact_hamiltonian
    nelec = Int(round(Z))
    truncated_scf = GaussletBases._closed_shell_density_density_scf(
        one_body.hamiltonian,
        one_body.overlap,
        truncated_interaction;
        nelec = nelec,
        maxiter = 100,
        damping = 0.25,
        tol = 1.0e-10,
    )
    assembled_scf = GaussletBases._closed_shell_density_density_scf(
        one_body.hamiltonian,
        one_body.overlap,
        assembled_interaction;
        nelec = nelec,
        maxiter = 100,
        damping = 0.25,
        tol = 1.0e-10,
    )
    exact_scf = GaussletBases._closed_shell_density_density_scf(
        exact_h1,
        exact_overlap,
        exact_interaction;
        nelec = nelec,
        maxiter = 100,
        damping = 0.25,
        tol = 1.0e-10,
    )
    one_body_spectrum =
        GaussletBases._generalized_spectrum(one_body.hamiltonian, one_body.overlap)

    bare_moment_lmax = maximum(maximum(keys(blocks)) for blocks in assembly.shell_moment_blocks)
    interaction_moment_lmax =
        maximum(maximum(keys(blocks)) for blocks in assembly.shell_interaction_moment_blocks)
    required_product_lmax = 2 * one_body.exact_common_lmax

    @printf(io, "\n=== Ne interaction payload diagnostic: N = %d ===\n", order)
    @printf(io, "  nr = %d  total_dim = %d  exact_common_lmax = %d  required_product_lmax = %d\n",
        size(radial_ops.overlap, 1),
        size(one_body.hamiltonian, 1),
        one_body.exact_common_lmax,
        required_product_lmax)
    @printf(io, "  overlap_identity_error = %.6e\n", opnorm(one_body.overlap - I, Inf))
    @printf(io, "  one_body_symmetry_error = %.6e\n", opnorm(one_body.hamiltonian - transpose(one_body.hamiltonian), Inf))
    @printf(io, "  interaction_symmetry_error(truncated) = %.6e\n", opnorm(truncated_interaction - transpose(truncated_interaction), Inf))
    @printf(io, "  interaction_symmetry_error(assembled) = %.6e\n", opnorm(assembled_interaction - transpose(assembled_interaction), Inf))
    @printf(io, "  interaction_helper_consistency = %.6e\n", opnorm(assembled_interaction - helper_interaction, Inf))
    @printf(io, "  bare_moment_lmax = %d  interaction_moment_lmax = %d\n", bare_moment_lmax, interaction_moment_lmax)
    @printf(io, "  low one-body eigenvalues = %s\n", repr(one_body_spectrum[1:min(top, length(one_body_spectrum))]))
    @printf(io, "  projected_exact_interaction_error(truncated) = %.6e\n",
        opnorm(projected_truncated - exact_interaction, Inf))
    @printf(io, "  projected_exact_interaction_error(assembled) = %.6e\n",
        opnorm(projected_assembled - exact_interaction, Inf))
    @printf(io, "  interaction_diff(truncated->assembled) = %.6e\n", opnorm(assembled_interaction - truncated_interaction, Inf))
    @printf(io, "  projected_interaction_diff(truncated->assembled) = %.6e\n",
        opnorm(projected_assembled - projected_truncated, Inf))
    @printf(io, "  truncated V row summary = %s\n", repr(row_summary(truncated_interaction, top)))
    @printf(io, "  assembled V row summary = %s\n", repr(row_summary(assembled_interaction, top)))
    @printf(io, "  exact projected V row summary = %s\n", repr(row_summary(exact_interaction, min(top, size(exact_interaction, 1)))))
    @printf(io, "  shell weight/moment summary = %s\n", repr(low_shell_report(assembly)))
    print_decomposition(io, "truncated", truncated_scf, one_body.hamiltonian, truncated_interaction)
    print_decomposition(io, "assembled", assembled_scf, one_body.hamiltonian, assembled_interaction)
    print_decomposition(io, "exact_common", exact_scf, exact_h1, exact_interaction)
end

function main(args)
    values = parse_kv_args(args)
    orders = get_int_list_arg(values, "orders", [10, 15, 32])
    Z = get_arg(values, "Z", 10.0)
    lmax = get_arg(values, "lmax", 6)
    top = get_arg(values, "top", 8)

    @printf("Angular Ne interaction payload diagnostic\n")
    @printf("  orders = %s\n", repr(orders))
    @printf("  setup = legacy_strict_trim, uniform shell orders, direct repo payload\n")
    for order in orders
        report_case(stdout, order; Z = Z, lmax = lmax, top = top)
    end
end

main(ARGS)
