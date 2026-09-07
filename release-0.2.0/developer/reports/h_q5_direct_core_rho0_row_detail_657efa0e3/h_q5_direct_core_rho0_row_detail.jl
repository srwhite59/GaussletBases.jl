#!/usr/bin/env julia

using GaussletBases
using LinearAlgebra
using Printf
using SpecialFunctions: erf
using Statistics

const GB = GaussletBases
const CFR = GB.CartesianFinalBasisRealization
const CRG = GB.CartesianResidualGaussians
const HEAD_TAG = readchomp(`git rev-parse --short HEAD`)
const RUN_DIR = "/Users/srw/dmrgtmp/h_q5_direct_core_rho0_row_detail_$(HEAD_TAG)"
const PIDFILE = joinpath(RUN_DIR, "h_q5_direct_core_rho0_row_detail.pid")
const ALPHA = parse(Float64, get(ENV, "RHO0_DIRECT_CORE_ALPHA", "8.0"))

sym(A) = Matrix{Float64}(0.5 .* (A .+ transpose(A)))
gaussian_potential(beta, r) = r == 0.0 ? 2.0 * sqrt(beta / pi) : erf(sqrt(beta) * r) / r

function write_table(path, rows, fields)
    open(path, "w") do io
        println(io, join(String.(fields), '\t'))
        for row in rows
            println(io, join((getproperty(row, field) for field in fields), '\t'))
        end
    end
end

function screened_coulomb_expansion(alpha)
    expansion = GB.coulomb_gaussian_expansion(doacc = false)
    coefficients = similar(expansion.coefficients)
    exponents = similar(expansion.exponents)
    for index in eachindex(coefficients)
        zeta = expansion.exponents[index]
        coefficients[index] = expansion.coefficients[index] *
            (alpha / (alpha + zeta))^(3 / 2)
        exponents[index] = alpha * zeta / (alpha + zeta)
    end
    return GB.CoulombGaussianExpansion(coefficients, exponents;
        del = expansion.del, s = expansion.s, c = expansion.c, maxu = expansion.maxu)
end

function terminal_screen_center_matrix(base, screen_expansion, center)
    basis = base.terminal_basis
    bundles = base.parent.parent_axis_bundle_object
    pgdg = Tuple(GB._nested_axis_pgdg(bundles, axis) for axis in (:x, :y, :z))
    matrix = zeros(Float64, basis.final_dimension, basis.final_dimension)
    factors = ntuple(axis ->
        CFR._r3a_centered_factor_terms(pgdg[axis], screen_expansion, center[axis]), 3)
    CFR._accumulate_terminal_gaussian_sum!(
        matrix, basis, screen_expansion.coefficients, factors[1], factors[2], factors[3];
        scale = 1.0)
    return matrix
end

function exact_base_screen_unit(base, alpha)
    expansion = screened_coulomb_expansion(alpha)
    center = Tuple(Float64.(base.input.locations[1]))
    return sym(terminal_screen_center_matrix(base, expansion, center))
end

function u0_base_unit(base, alpha)
    centers = reshape(Float64[base.input.locations[1]...], 1, 3)
    widths = fill(inv(sqrt(alpha)), 1, 3)
    expansion = GB.coulomb_gaussian_expansion(doacc = false)
    pair_terms = CRG._mwg_axis_pairs(
        base.parent.parent_axis_bundle_object, expansion, centers, widths)
    cols = CRG._terminal_mwg_fixed_residual(
        base.terminal_basis, base.parent.parent_axis_bundle_object,
        pair_terms, expansion.coefficients)
    return vec(cols[:, 1])
end

axis_bundle(base, axis) = GB._nested_axis_pgdg(base.parent.parent_axis_bundle_object, axis)
axis_values(base, field) = ntuple(i -> getproperty(axis_bundle(base, (:x, :y, :z)[i]), field), 3)

function terminal_weights(base)
    weights = zeros(Float64, base.terminal_basis.final_dimension)
    wx, wy, wz = axis_values(base, :weights)
    for block in base.terminal_basis.blocks
        support_weights = [wx[s[1]] * wy[s[2]] * wz[s[3]] for s in block.support_states]
        if isnothing(block.coefficients)
            weights[block.column_range] .= support_weights
        else
            weights[block.column_range] .= vec(transpose(block.coefficients) * support_weights)
        end
    end
    return weights
end

function local_spacing(centers, i)
    values = Float64[]
    i > firstindex(centers) && push!(values, abs(centers[i] - centers[i - 1]))
    i < lastindex(centers) && push!(values, abs(centers[i + 1] - centers[i]))
    return isempty(values) ? NaN : minimum(values)
end

function main()
    mkpath(RUN_DIR)
    println("probe\t", @__FILE__)
    println("pid\t", getpid())
    open(PIDFILE, "w") do io
        println(io, getpid())
    end
    println("pidfile\t", PIDFILE)
    println("run_dir\t", RUN_DIR)
    println("scope\tH_q5_core_spacing_0.3_base_terminal_direct_core_only")

    system = (; atom_symbols = ["H"], nuclear_charges = [1.0],
        atom_locations = [(0.0, 0.0, 0.0)], nup = 1, ndn = 0)
    basis = (; q = 5, core_spacing = 0.3, radius = 5.0,
        reference_spacing = 1.0, d = 0.3)
    t_build = @elapsed base = GB.cartesian_base_working_basis(system; basis)
    t_j = @elapsed J = exact_base_screen_unit(base, ALPHA)
    t_u = @elapsed u = u0_base_unit(base, ALPHA)
    weights = terminal_weights(base)
    jw = J * weights
    row_action = jw ./ weights
    residual = jw .- u .* weights

    block = first(base.terminal_basis.blocks)
    rows = NamedTuple[]
    xs, ys, zs = axis_values(base, :centers)
    wx, wy, wz = axis_values(base, :weights)
    support_weight = [wx[s[1]] * wy[s[2]] * wz[s[3]] for s in block.support_states]
    direct_total = sum(abs2, residual[block.column_range])
    all_total = sum(abs2, residual)
    small_cut = quantile(weights[block.column_range], 0.25)
    dd = base.terminal_due_diligence.terminal_rows[1]
    lo = first.(dd.outer_box)
    hi = last.(dd.outer_box)
    for (local_index, col) in enumerate(block.column_range)
        state = block.support_states[local_index]
        coord = (xs[state[1]], ys[state[2]], zs[state[3]])
        r = sqrt(sum(abs2, coord))
        parent_weight = support_weight[local_index]
        point = gaussian_potential(ALPHA, r)
        smooth = gaussian_potential(ALPHA / 2.0, r)
        is_boundary = any(axis -> state[axis] == lo[axis] || state[axis] == hi[axis], 1:3)
        is_core_center = r <= 1.0e-10
        is_small_weight = weights[col] <= small_cut
        sx, sy, sz = local_spacing(xs, state[1]), local_spacing(ys, state[2]),
            local_spacing(zs, state[3])
        push!(rows, (; row_index = col, local_index, state_x = state[1],
            state_y = state[2], state_z = state[3], x = coord[1], y = coord[2],
            z = coord[3], nearest_center_distance = r,
            terminal_support_weight = weights[col],
            parent_product_weight = parent_weight,
            parent_minus_terminal_weight = parent_weight - weights[col],
            u_direct = u[col], row_action = row_action[col],
            ratio_row_action_over_u_direct = row_action[col] / u[col],
            row_residual = residual[col],
            residual_abs = abs(residual[col]),
            residual_sq_fraction_direct_core = abs2(residual[col]) / direct_total,
            residual_sq_fraction_all = abs2(residual[col]) / all_total,
            analytic_point_potential = point,
            analytic_equal_width_gaussian_average = smooth,
            row_action_minus_point = row_action[col] - point,
            u_direct_minus_point = u[col] - point,
            u_direct_minus_equal_width_average = u[col] - smooth,
            boundary = is_boundary,
            core_center = is_core_center,
            small_weight_q25 = is_small_weight,
            local_spacing_x = sx, local_spacing_y = sy, local_spacing_z = sz,
            local_spacing_min = minimum((sx, sy, sz)),
            local_spacing_geom = (sx * sy * sz)^(1 / 3)))
    end

    direct = collect(block.column_range)
    parent_weights = copy(weights)
    parent_weights[direct] .= support_weight
    parent_jw = J * parent_weights
    parent_residual = parent_jw .- u .* parent_weights
    point_error = row_action[direct] .- [row.analytic_point_potential for row in rows]
    smooth_error = u[direct] .- [row.analytic_equal_width_gaussian_average for row in rows]
    summary = [(; alpha = ALPHA, dimension = length(weights),
        direct_core_count = length(direct),
        build_elapsed_s = t_build, J_elapsed_s = t_j, u_elapsed_s = t_u,
        terminal_row_rel = norm(residual) / norm(jw),
        parent_product_row_rel = norm(parent_residual) / norm(parent_jw),
        parent_terminal_weight_maxabs = maximum(abs.(support_weight .- weights[direct])),
        direct_core_row_rel_global = norm(residual[direct]) / norm(jw),
        direct_core_residual_sq_fraction = direct_total / all_total,
        row_action_point_rel = norm(point_error) /
            norm([row.analytic_point_potential for row in rows]),
        row_action_point_maxabs = maximum(abs.(point_error)),
        u_equal_width_average_rel = norm(smooth_error) /
            norm([row.analytic_equal_width_gaussian_average for row in rows]),
        u_equal_width_average_maxabs = maximum(abs.(smooth_error)),
        worst_row = rows[argmax([row.residual_abs for row in rows])].row_index,
        worst_state = string((rows[argmax([row.residual_abs for row in rows])].state_x,
            rows[argmax([row.residual_abs for row in rows])].state_y,
            rows[argmax([row.residual_abs for row in rows])].state_z)),
        conclusion = "terminal_direct_core_weight_matches_parent_product_but_u_direct_is_gaussian_smoothed_not_row_action_point_potential")]

    write_table(joinpath(RUN_DIR, "direct_core_rows.tsv"), rows, propertynames(first(rows)))
    write_table(joinpath(RUN_DIR, "summary.tsv"), summary, propertynames(first(summary)))
    @printf("RESULT\trow_rel=%.12e\tparent_row_rel=%.12e\trow_action_point_rel=%.12e\tu_smooth_rel=%.12e\n",
        summary[1].terminal_row_rel, summary[1].parent_product_row_rel,
        summary[1].row_action_point_rel, summary[1].u_equal_width_average_rel)
    println("run_dir\t", RUN_DIR)
end

t = @elapsed main()
@printf("elapsed_s=%.6f\n", t)
