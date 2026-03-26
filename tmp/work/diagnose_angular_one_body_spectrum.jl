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

function shell_exact_projector(shell, exact_channels)
    lookup = Dict{YlmChannel,Int}(
        shell.injected_channels.channel_data[i] => i for i in eachindex(shell.injected_channels.channel_data)
    )
    rows = Int[]
    for channel in exact_channels.channel_data
        row = get(lookup, channel, 0)
        row == 0 || push!(rows, row)
    end
    isempty(rows) && return zeros(Float64, shell.final_count, shell.final_count)
    exact_rows = shell.injected_overlap[rows, :]
    return transpose(exact_rows) * exact_rows
end

function top_shell_weights(weights::AbstractVector{<:Real}; ntop::Int = 5)
    order = sortperm(weights; rev = true)
    return [(index = idx, weight = Float64(weights[idx])) for idx in order[1:min(ntop, length(order))]]
end

function main(args)
    values = parse_kv_args(args)
    Z = get_arg(values, "Z", 4.0)
    order = get_arg(values, "N", 15)
    lmax = get_arg(values, "lmax", 2)
    top = get_arg(values, "top", 8)
    s = get_arg(values, "s", 0.2)
    rmax = get_arg(values, "rmax", 30.0)

    radial_ops = paper_style_radial_ops(; Z = Z, lmax = lmax, s = s, rmax = rmax)
    benchmark = build_atomic_injected_angular_one_body_benchmark(
        radial_ops;
        shell_orders = fill(order, length(radial_ops.shell_centers_r)),
    )
    eig = GaussletBases._generalized_eigensystem(benchmark.hamiltonian, benchmark.overlap)
    full_spectrum = eig.values
    exact_spectrum =
        GaussletBases._generalized_spectrum(benchmark.exact_hamiltonian, benchmark.exact_overlap)
    noninteracting_closed_shell = 2.0 * (full_spectrum[1] + full_spectrum[2])

    shell_projectors = [
        shell_exact_projector(shell, benchmark.exact_channels) for shell in benchmark.angular_assembly.shells
    ]

    @printf("Angular one-body spectrum diagnostic\n")
    @printf("  Z = %.1f  lmax = %d  uniform_order = %d\n", Z, lmax, order)
    @printf("  shell_count = %d  total_dim = %d  exact_common_lmax = %d\n",
            length(benchmark.angular_assembly.shells),
            size(benchmark.hamiltonian, 1),
            benchmark.exact_common_lmax)
    @printf("  shell_radii = %s\n", repr(benchmark.angular_assembly.shell_radii))
    @printf("  shell_orders = %s\n", repr(benchmark.angular_assembly.shell_orders))
    @printf("  shell_kinetic_lcap = %s\n", repr(benchmark.angular_assembly.shell_kinetic_lcap))
    @printf("  shell_kinetic_lexpand = %s\n", repr(benchmark.angular_assembly.shell_kinetic_lexpand))
    @printf("  shell_kinetic_tail = %s\n", repr(benchmark.angular_assembly.shell_kinetic_tail))
    @printf("  low full eigenvalues = %s\n", repr(full_spectrum[1:min(top, length(full_spectrum))]))
    @printf("  low exact eigenvalues = %s\n", repr(exact_spectrum[1:min(top, length(exact_spectrum))]))
    @printf("  closed_shell_noninteracting = %.12f\n", noninteracting_closed_shell)

    for state in 1:min(top, length(full_spectrum))
        coeffs = eig.coefficients[:, state]
        total_norm = sum(abs2, coeffs)
        shell_weights = Float64[]
        exact_weight = 0.0
        for (shell_index, range) in enumerate(benchmark.shell_ranges)
            shell_coeffs = coeffs[range]
            push!(shell_weights, sum(abs2, shell_coeffs) / total_norm)
            exact_weight += real(dot(shell_coeffs, shell_projectors[shell_index] * shell_coeffs)) / total_norm
        end
        sector = exact_weight >= 0.5 ? "exact_common_dominant" : "mixed_remainder_dominant"
        @printf("\nstate %d\n", state)
        @printf("  eigenvalue = %.12f\n", full_spectrum[state])
        @printf("  exact_common_fraction = %.6f  mixed_fraction = %.6f  sector = %s\n",
                exact_weight, max(0.0, 1.0 - exact_weight), sector)
        @printf("  top_shell_weights = %s\n", repr(top_shell_weights(shell_weights)))
    end
end

main(ARGS)
