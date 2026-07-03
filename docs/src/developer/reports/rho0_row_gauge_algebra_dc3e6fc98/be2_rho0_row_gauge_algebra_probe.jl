#!/usr/bin/env julia

include(joinpath(@__DIR__, "protected_localized_rho0_galerkin_audit.jl"))

using LinearAlgebra
using Printf
using SpecialFunctions: erf
using Statistics

const BE2_RECIPE = (;
    label = "Be2_row_gauge",
    atom_symbols = ["Be", "Be"],
    nuclear_charges = [4.0, 4.0],
    atom_locations = [(0.0, 0.0, -2.5), (0.0, 0.0, 2.5)],
    nup = 4,
    ndn = 4,
    ns = 4,
    nesting = :pqs,
    core_spacing = 0.3,
    xmax_parallel = 7.5,
    xmax_transverse = 5.0,
    basisname = "cc-pVDZ",
    basisfile = "/Users/srw/Dropbox/GaussletModules/BasisSets",
    lmax = 1,
    uncontracted = false,
    width_filtering = nothing,
    residual_occupation_cutoff = 1.0e-6,
    tau_neg_abs = 1.0e-12,
    tau_neg_rel = 1.0e-12,
    tau_merge_abs = 1.0e-12,
    tau_merge_rel = 1.0e-12,
)

const ROW_ALPHA = parse(Float64, get(ENV, "RHO0_ROW_ALPHA", "8.0"))
const ROW_N_SCREENS = (0.0, 10.0, 18.0, 22.0)
const ROW_HEAD_TAG = readchomp(`git rev-parse --short HEAD`)
const ROW_RUN_DIR = "/Users/srw/dmrgtmp/rho0_row_gauge_algebra_$(ROW_HEAD_TAG)"
const ROW_PIDFILE = joinpath(ROW_RUN_DIR, "rho0_row_gauge_algebra.pid")

function row_write(path, rows, fields)
    open(path, "w") do io
        println(io, join(String.(fields), '\t'))
        for row in rows
            println(io, join((getproperty(row, field) for field in fields), '\t'))
        end
    end
end

function row_stage!(rows, label, f)
    GC.gc()
    timed = @timed f()
    push!(rows, (; label, seconds = timed.time, alloc_mib = timed.bytes / 2.0^20,
        gc_s = timed.gctime))
    @printf("stage\t%s\telapsed_s=%.6f\talloc_mib=%.3f\tgc_s=%.6f\n",
        label, timed.time, timed.bytes / 2.0^20, timed.gctime)
    flush(stdout)
    return timed.value
end

function augmented_screen_unit(inputs, residual, alpha)
    screen_expansion = screened_coulomb_expansion(alpha)
    raw = screen_raw_blocks(inputs, screen_expansion)
    n = residual.base_dimension + residual.residual_dimension
    screen = zeros(Float64, n, n)
    for (center_index, center) in pairs(center_tuples(inputs.base.input.locations))
        U_GG = terminal_screen_center_matrix(inputs.base, screen_expansion, center)
        screen .+= CRG.transform_augmented_operator(
            U_GG, -raw.mixed[center_index], -raw.self[center_index], residual)
    end
    return sym(screen)
end

function terminal_density_weights(base)
    basis = base.terminal_basis
    bundles = base.parent.parent_axis_bundle_object
    weights = zeros(Float64, basis.final_dimension)
    for block in basis.blocks
        support_weights = CRG._mwg_support_weights(block.support_states, bundles)
        if block.coefficients === nothing
            weights[block.column_range] .= support_weights
        else
            weights[block.column_range] .= vec(transpose(block.coefficients) * support_weights)
        end
    end
    return weights
end

function residual_density_weights(operators, residual)
    centers, widths = CRG.moment_matched_gaussians(operators, residual)
    effective = getfield(GB, :_qwrg_effective_gaussians)(centers, widths)
    axis_weights = ntuple(axis ->
        getfield(GB, :_qwrg_supplement_integral_weights)(effective[axis]), 3)
    return axis_weights[1] .* axis_weights[2] .* axis_weights[3]
end

m_density_weights(inputs, compact, residual) = vcat(
    terminal_density_weights(inputs.base), residual_density_weights(compact.products, residual))

function row_action_check(label, J, u, w)
    jw = J * w
    residual = jw .- u .* w
    offdiag = Matrix{Float64}(J)
    for i in axes(offdiag, 1)
        offdiag[i, i] = 0.0
    end
    return (; label, dimension = length(w), finite = all(isfinite, J) &&
            all(isfinite, u) && all(isfinite, w),
        weight_min = minimum(w), weight_max = maximum(w),
        weight_negative_count = count(<(0.0), w),
        row_abs = norm(residual),
        row_rel = norm(residual) / max(norm(jw), eps(Float64)),
        row_maxabs = maximum(abs.(residual)),
        diag_abs = norm(diag(J) .- u),
        diag_rel = norm(diag(J) .- u) / max(norm(u), eps(Float64)),
        offdiag_fro = norm(offdiag),
        offdiag_rel = norm(offdiag) / max(norm(J), eps(Float64)),
        jw_norm = norm(jw))
end

function direct_u_from_row_action(J, w)
    u = similar(w)
    for i in eachindex(w)
        u[i] = abs(w[i]) <= 1.0e-14 ? NaN : dot(view(J, i, :), w) / w[i]
    end
    return u
end

function occupied_shift_rows(label, H, delta, nocc)
    values, vectors = eigen(Symmetric(sym(H)))
    rows = NamedTuple[]
    total = 0.0
    for orbital in 1:nocc
        vector = vectors[:, orbital]
        shift = dot(vector, delta * vector)
        total += 2.0 * shift
        push!(rows, (; label, orbital, h1 = values[orbital], unit_shift = shift,
            closed_shell_unit_shift_running = total))
    end
    return rows
end

function linearity_rows(delta_unit)
    rows = NamedTuple[]
    zero = 0.0 .* delta_unit
    push!(rows, (; test = "N_screen_zero", error_inf = norm(zero, Inf),
        error_fro = norm(zero), reference_fro = 0.0))
    for n in ROW_N_SCREENS
        lhs = n .* delta_unit
        rhs = n .* delta_unit
        push!(rows, (; test = "explicit_unit_scale_N=$(n)",
            error_inf = norm(lhs .- rhs, Inf), error_fro = norm(lhs .- rhs),
            reference_fro = norm(rhs)))
    end
    push!(rows, (; test = "N18_vs_18_over_10_N10",
        error_inf = norm(18.0 .* delta_unit .- 1.8 .* (10.0 .* delta_unit), Inf),
        error_fro = norm(18.0 .* delta_unit .- 1.8 .* (10.0 .* delta_unit)),
        reference_fro = norm(18.0 .* delta_unit)))
    push!(rows, (; test = "N22_vs_22_over_10_N10",
        error_inf = norm(22.0 .* delta_unit .- 2.2 .* (10.0 .* delta_unit), Inf),
        error_fro = norm(22.0 .* delta_unit .- 2.2 .* (10.0 .* delta_unit)),
        reference_fro = norm(22.0 .* delta_unit)))
    return rows
end

function gaussian_normalization_checks(alpha)
    radius_max = 12.0 / sqrt(alpha)
    points = 200_000
    dr = radius_max / points
    charge = 0.0
    self = 0.0
    for i in 1:points
        r = (i - 0.5) * dr
        rho = (alpha / pi)^(3 / 2) * exp(-alpha * r^2)
        potential = erf(sqrt(alpha) * r) / r
        shell = 4.0 * pi * r^2 * dr
        charge += shell * rho
        self += shell * rho * potential
    end
    expansion_center = screened_coulomb_expansion(alpha)(0.0)
    return [(; alpha, numeric_charge = charge, analytic_charge = 1.0,
        charge_error = charge - 1.0,
        numeric_self_energy = self,
        analytic_self_energy = sqrt(2.0 * alpha / pi),
        self_energy_error = self - sqrt(2.0 * alpha / pi),
        expansion_center_potential = expansion_center,
        analytic_center_potential = 2.0 * sqrt(alpha / pi),
        center_potential_error = expansion_center - 2.0 * sqrt(alpha / pi))]
end

function main()
    mkpath(ROW_RUN_DIR)
    println("probe\t", @__FILE__)
    println("pid\t", getpid())
    open(ROW_PIDFILE, "w") do io
        println(io, getpid())
    end
    println("pidfile\t", ROW_PIDFILE)
    println("run_dir\t", ROW_RUN_DIR)
    println("authority\tHP-RG-RHO0-GAL-AUDIT-01")
    println("interpretation\talgebra_only_no_HF")
    println("rho0\tsingle normalized spherical Gaussian per Be center alpha=", ROW_ALPHA)

    stages = NamedTuple[]
    inputs = row_stage!(stages, "build Be2 inputs/base/supplement", () ->
        build_inputs(BE2_RECIPE))
    residual = row_stage!(stages, "compact residual parent M", () ->
        compact_residual_basis(inputs, BE2_RECIPE))
    compact = row_stage!(stages, "compact M operators and inherited Vee", () ->
        augmented_data(inputs, residual))
    geometry = row_stage!(stages, "staged protected geometry", () ->
        protected_geometry(inputs, BE2_RECIPE))
    loc = row_stage!(stages, "localized protected L transform", () ->
        localized_transform(geometry))
    J_M = row_stage!(stages, "exact Galerkin rho0 screen unit in M", () ->
        augmented_screen_unit(inputs, residual, ROW_ALPHA))
    u0_M = row_stage!(stages, "direct row-gauge u0 unit in M", () ->
        u0_m_unit(inputs, residual, compact, ROW_ALPHA))
    J_L = row_stage!(stages, "exact Galerkin rho0 screen unit in L", () ->
        protected_screen_unit(inputs, geometry, loc, ROW_ALPHA))
    H_M = row_stage!(stages, "compact M one-body", () ->
        h1_from_augmented(inputs, compact))
    H_F = row_stage!(stages, "exact one-body in F=[Z,MQperp]", () ->
        protected_onebody(inputs, geometry, compact).one_body_hamiltonian)
    H_L = row_stage!(stages, "exact one-body in L", () -> localized_operator(H_F, loc))

    w_M = row_stage!(stages, "M row density weights", () ->
        m_density_weights(inputs, compact, residual))
    ML = hcat(loc.C, loc.Qp) * loc.W
    w_L_projected = transpose(ML) * w_M
    u0_L_projected = direct_u_from_row_action(J_L, w_L_projected)

    checks = NamedTuple[]
    push!(checks, row_action_check("M_same_basis_direct_u0_M", J_M, u0_M, w_M))
    push!(checks, row_action_check("M_same_basis_u_from_Jw", J_M,
        direct_u_from_row_action(J_M, w_M), w_M))
    push!(checks, row_action_check("L_using_u0_M_and_w_M", J_L, u0_M, w_M))
    push!(checks, row_action_check("L_using_u0_L_projected_and_w_L_projected",
        J_L, u0_L_projected, w_L_projected))

    delta_M = Matrix{Float64}(J_M)
    for i in eachindex(u0_M)
        delta_M[i, i] -= u0_M[i]
    end
    delta_L_u0M = Matrix{Float64}(J_L)
    for i in eachindex(u0_M)
        delta_L_u0M[i, i] -= u0_M[i]
    end

    weight_rows = [(; label = "M_weights", count = length(w_M), min = minimum(w_M),
        median = median(w_M), max = maximum(w_M), sum = sum(w_M),
        negative_count = count(<(0.0), w_M), near_zero_count = count(<(1.0e-14), abs.(w_M))),
        (; label = "L_projected_weights", count = length(w_L_projected),
            min = minimum(w_L_projected), median = median(w_L_projected),
            max = maximum(w_L_projected), sum = sum(w_L_projected),
            negative_count = count(<(0.0), w_L_projected),
            near_zero_count = count(<(1.0e-14), abs.(w_L_projected)))]
    weight_compare = [(; test = "w_L_projected_minus_w_M",
        abs = norm(w_L_projected .- w_M),
        rel = norm(w_L_projected .- w_M) / max(norm(w_M), eps(Float64)),
        maxabs = maximum(abs.(w_L_projected .- w_M)))]

    occupied_rows = vcat(
        occupied_shift_rows("M_delta_direct_u0_M", H_M, delta_M, BE2_RECIPE.nup),
        occupied_shift_rows("L_delta_u0_M_wM", H_L, delta_L_u0M, BE2_RECIPE.nup))

    summary = [(; case_label = BE2_RECIPE.label,
        alpha = ROW_ALPHA,
        base_dimension = inputs.base.terminal_basis.final_dimension,
        compact_R_dimension = residual.residual_dimension,
        M_dimension = length(w_M),
        protected_originals = geometry.protected_original_count,
        broad_Z = size(geometry.Z_broad, 2),
        Z_dimension = size(geometry.Z, 2),
        L_dimension = size(J_L, 1),
        B_min = geometry.b_min,
        B_lt_0p99 = geometry.b_lt_0p99,
        M_row_rel = checks[1].row_rel,
        L_u0M_row_rel = checks[3].row_rel,
        L_projected_row_rel = checks[4].row_rel,
        weight_rel = weight_compare[1].rel,
        note = "No HF; row-gauge algebra only")]

    row_write(joinpath(ROW_RUN_DIR, "summary.tsv"), summary, propertynames(first(summary)))
    row_write(joinpath(ROW_RUN_DIR, "row_gauge_checks.tsv"), checks, propertynames(first(checks)))
    row_write(joinpath(ROW_RUN_DIR, "linearity.tsv"), linearity_rows(delta_M),
        (:test, :error_inf, :error_fro, :reference_fro))
    row_write(joinpath(ROW_RUN_DIR, "weights.tsv"), weight_rows, propertynames(first(weight_rows)))
    row_write(joinpath(ROW_RUN_DIR, "weight_compare.tsv"), weight_compare,
        propertynames(first(weight_compare)))
    row_write(joinpath(ROW_RUN_DIR, "occupied_orbital_shift_proxies.tsv"), occupied_rows,
        propertynames(first(occupied_rows)))
    row_write(joinpath(ROW_RUN_DIR, "analytic_gaussian_checks.tsv"),
        gaussian_normalization_checks(ROW_ALPHA),
        (:alpha, :numeric_charge, :analytic_charge, :charge_error,
            :numeric_self_energy, :analytic_self_energy, :self_energy_error,
            :expansion_center_potential, :analytic_center_potential,
            :center_potential_error))
    row_write(joinpath(ROW_RUN_DIR, "stages.tsv"), stages,
        (:label, :seconds, :alloc_mib, :gc_s))

    @printf("RESULT\tM_row_rel=%.12e\tL_u0M_row_rel=%.12e\tL_projected_row_rel=%.12e\tweight_rel=%.12e\n",
        checks[1].row_rel, checks[3].row_rel, checks[4].row_rel, weight_compare[1].rel)
    println("run_dir\t", ROW_RUN_DIR)
end

t = @elapsed main()
@printf("elapsed_s=%.6f\n", t)
