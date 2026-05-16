using Dates
using LinearAlgebra
using Printf

using GaussletBases

const GB = GaussletBases

set_timing!(false)

const R_BOHR = 4.0
const NUCLEAR_CHARGES = [1.0, 1.0]
const CORE_SPACING = 0.5
const XMAX_PARALLEL = 6.0
const XMAX_TRANSVERSE = 4.0
const BOND_AXIS = :z
const NSIDE = 5
const SUPPLEMENT_ATOM = "H"
const SUPPLEMENT_BASIS = "cc-pVTZ"
const SUPPLEMENT_LMAX = 1

# The high-order old-standard-aligned reference rows used the nearest-center
# route. Keep that convention here so the chemistry numbers are comparable.
const INTERACTION_TREATMENT = :ggt_nearest
const ENDCAP_Q = 4
const ENDCAP_L = 4

const BO_TOTAL_R4 = -1.0163902529471283
const HIGH_ORDER_OLD_SP_DIM_R4 = 481
const HIGH_ORDER_OLD_SP_ED_TOTAL_R4 = -1.0156138376908870
const HIGH_ORDER_ALL_SHARED_Q4_DIM_R4 = 461
const HIGH_ORDER_ALL_SHARED_Q4_ED_TOTAL_R4 = -1.015663743783

const HF_TOL = 1.0e-10
const HF_MAXIT = 250
const ED_TOL = 1.0e-10
const ED_MAXIT = 2000

const DIAG2PTLE_CANDIDATES = (
    "/Users/srwhite/Dropbox/GaussletModules/Diag2ptle.jl",
    "/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/Diag2ptle.jl",
)

function _first_existing_path(paths)
    for path in paths
        isfile(path) && return path
    end
    error("Diag2ptle.jl not found in candidates: $(join(paths, ", "))")
end

const DIAG2PTLE_PATH = _first_existing_path(DIAG2PTLE_CANDIDATES)
const DIAG2PTLE_DIR = dirname(DIAG2PTLE_PATH)
DIAG2PTLE_DIR in LOAD_PATH || pushfirst!(LOAD_PATH, DIAG2PTLE_DIR)
const DIAG2PTLE_LOAD_START = time()
include(DIAG2PTLE_PATH)
const DIAG2PTLE_LOAD_SECONDS = time() - DIAG2PTLE_LOAD_START
const DIAG2PTLE =
    isdefined(Main, :Diag2ptle) ? getfield(Main, :Diag2ptle) :
    isdefined(@__MODULE__, :Diag2ptle) ? getfield(@__MODULE__, :Diag2ptle) :
    error("Diag2ptle include did not define a visible Diag2ptle module at $(DIAG2PTLE_PATH)")

mutable struct StageTimer
    rows::Vector{NamedTuple}
end

StageTimer() = StageTimer(NamedTuple[])

function timed!(timer::StageTimer, label::AbstractString, f::Function)
    GC.gc()
    result = @timed f()
    push!(
        timer.rows,
        (
            label = String(label),
            time = result.time,
            bytes = result.bytes,
            gctime = result.gctime,
        ),
    )
    return result.value
end

timed!(f::Function, timer::StageTimer, label::AbstractString) = timed!(timer, label, f)

fmt(value) = value === missing ? "NA" : string(value)
fmt_float(value) = value === missing ? "NA" : @sprintf("%.12e", Float64(value))

symmetrize(matrix::AbstractMatrix{<:Real}) =
    Matrix{Float64}(0.5 .* (matrix .+ transpose(matrix)))

function symmetry_error(matrix::AbstractMatrix{<:Real})
    isempty(matrix) && return 0.0
    return norm(matrix - transpose(matrix), Inf)
end

function width_status(widths::AbstractMatrix)
    all(isnan, widths) && return :all_nan
    all(isfinite, widths) && all(>(0.0), widths) && return :finite_positive
    return :mixed_or_invalid
end

ensure_diag2ptle_loaded() = DIAG2PTLE

function fock_matrix(
    H::AbstractMatrix{<:Real},
    V::AbstractMatrix{<:Real},
    rho::AbstractMatrix{<:Real},
)
    density_diagonal = diag(rho)
    return symmetrize(H + 2.0 * Diagonal(V * density_diagonal) - rho .* V)
end

function hf_energy(
    H::AbstractMatrix{<:Real},
    V::AbstractMatrix{<:Real},
    rho::AbstractMatrix{<:Real},
)
    density_diagonal = diag(rho)
    return (
        2.0 * tr(rho * H) +
        2.0 * dot(density_diagonal, V * density_diagonal) -
        dot(vec(rho), vec(V .* rho))
    )
end

function orbital_residual(F::AbstractMatrix{<:Real}, c::AbstractVector{<:Real})
    Fc = F * c
    epsilon = dot(c, Fc)
    return norm(Fc - epsilon * c)
end

function closed_shell_hf(H::AbstractMatrix{<:Real}, V::AbstractMatrix{<:Real})
    h_eigen = eigen(Symmetric(symmetrize(H)))
    c = Vector{Float64}(h_eigen.vectors[:, 1])
    c ./= norm(c)

    previous_energy = Inf
    last_energy_delta = Inf
    last_density_delta = Inf
    iterations = 0
    converged = false
    F = fock_matrix(H, V, c * transpose(c))

    for iteration in 1:HF_MAXIT
        iterations = iteration
        rho = c * transpose(c)
        F = fock_matrix(H, V, rho)
        f_eigen = eigen(Symmetric(F))
        c_new = Vector{Float64}(f_eigen.vectors[:, 1])
        dot(c_new, c) < 0.0 && (c_new .*= -1.0)
        c_new ./= norm(c_new)
        rho_new = c_new * transpose(c_new)
        F_new = fock_matrix(H, V, rho_new)
        energy = hf_energy(H, V, rho_new)
        energy_delta = isfinite(previous_energy) ? abs(energy - previous_energy) : Inf
        density_delta = norm(rho_new - rho, Inf)
        residual = orbital_residual(F_new, c_new)
        c = c_new
        previous_energy = energy
        last_energy_delta = energy_delta
        last_density_delta = density_delta
        F = F_new
        if energy_delta < HF_TOL && density_delta < sqrt(HF_TOL) && residual < sqrt(HF_TOL)
            converged = true
            break
        end
    end

    rho = c * transpose(c)
    F = fock_matrix(H, V, rho)
    return (
        converged = converged,
        iterations = iterations,
        orbital = c,
        rho = rho,
        fock = F,
        electronic_energy = hf_energy(H, V, rho),
        total_energy = hf_energy(H, V, rho) + 1.0 / R_BOHR,
        last_energy_delta = last_energy_delta,
        last_density_delta = last_density_delta,
        orbital_residual = orbital_residual(F, c),
    )
end

function two_electron_apply(
    H::AbstractMatrix{<:Real},
    V::AbstractMatrix{<:Real},
    psi::AbstractMatrix{<:Real},
)
    return H * psi + psi * H + V .* psi
end

function two_electron_residual(H, V, psi, energy)
    psi_value = Matrix{Float64}(psi)
    psi_value ./= norm(psi_value)
    return norm(two_electron_apply(H, V, psi_value) - Float64(energy) * psi_value)
end

function davidson_ed(
    H::Matrix{Float64},
    V::Matrix{Float64},
    orbital::Vector{Float64},
    fock::Matrix{Float64},
)
    Diag2ptle = ensure_diag2ptle_loaded()
    psi0 = orbital * transpose(orbital)
    psi0 ./= norm(psi0)
    multV!(out, psi) = (out .= V .* psi)
    psi, energy, info = Diag2ptle.davidson_ground_sylv(
        H,
        multV!,
        psi0;
        Aprec = fock,
        tol = ED_TOL,
        maxit = ED_MAXIT,
        mmax = 192,
        keep = 24,
        block2 = true,
        verbose = false,
        warn_iter = 200,
        warn_every = 200,
    )
    return (
        route = :davidson_ground_sylv,
        electronic_energy = Float64(energy),
        total_energy = Float64(energy) + 1.0 / R_BOHR,
        residual = two_electron_residual(H, V, psi, energy),
        iterations = getproperty(info, :iters),
        solver_dimension = getproperty(info, :dim),
        reason = String(getproperty(info, :reason)),
        wavefunction_norm = norm(psi),
    )
end

function build_basis()
    return bond_aligned_homonuclear_qw_basis(
        bond_length = R_BOHR,
        core_spacing = CORE_SPACING,
        xmax_parallel = XMAX_PARALLEL,
        xmax_transverse = XMAX_TRANSVERSE,
        bond_axis = BOND_AXIS,
    )
end

function build_source_and_fixed(basis, expansion, policy::Symbol)
    bundles = GB._qwrg_bond_aligned_axis_bundles(basis, expansion)
    source = GB._nested_bond_aligned_diatomic_source(
        basis,
        bundles;
        bond_axis = BOND_AXIS,
        nside = NSIDE,
        term_coefficients = Float64.(expansion.coefficients),
        packet_kernel = :factorized_direct,
        shared_shell_layer_policy = policy,
        shared_shell_endcap_panel_q = ENDCAP_Q,
        shared_shell_endcap_panel_L = ENDCAP_L,
    )
    fixed = GB._nested_fixed_block(source)
    return source, fixed
end

function build_operators(fixed_block, supplement, expansion)
    operators = ordinary_cartesian_qiu_white_operators(
        fixed_block,
        supplement;
        nuclear_charges = NUCLEAR_CHARGES,
        expansion = expansion,
        interaction_treatment = INTERACTION_TREATMENT,
    )
    Hraw = Matrix{Float64}(operators.one_body_hamiltonian)
    Vraw = Matrix{Float64}(operators.interaction_matrix)
    H = symmetrize(Hraw)
    V = symmetrize(Vraw)
    return (
        operators = operators,
        H = H,
        V = V,
        dimension = size(H, 1),
        residual_count = operators.residual_count,
        overlap_error = norm(operators.overlap - I, Inf),
        h_symmetry_error = symmetry_error(Hraw),
        v_symmetry_error = symmetry_error(Vraw),
        width_status = width_status(operators.residual_widths),
        owner_set = sort!(collect(Set(operators.residual_nucleus_indices))),
    )
end

function route_metadata(label::Symbol, source, fixed_block, operators_diag)
    shared_columns = [size(layer.coefficient_matrix, 2) for layer in source.shared_shell_layers]
    shared_types = [nameof(typeof(layer)) for layer in source.shared_shell_layers]
    return (
        label = label,
        shared_layer_types = Tuple(shared_types),
        shared_layer_columns = Tuple(shared_columns),
        fixed_block_size = size(fixed_block.coefficient_matrix),
        fixed_overlap_error = norm(fixed_block.overlap - I, Inf),
        operator_dimension = operators_diag.dimension,
        residual_count = operators_diag.residual_count,
        operator_overlap_error = operators_diag.overlap_error,
        h_symmetry_error = operators_diag.h_symmetry_error,
        v_symmetry_error = operators_diag.v_symmetry_error,
        width_status = operators_diag.width_status,
        owner_set = Tuple(operators_diag.owner_set),
    )
end

function run_route!(timer::StageTimer, label::Symbol, policy::Symbol, basis, supplement, expansion)
    source, fixed = timed!(timer, "$(label) source and fixed block") do
        build_source_and_fixed(basis, expansion, policy)
    end
    operators_diag = timed!(timer, "$(label) ordinary QW operators") do
        build_operators(fixed, supplement, expansion)
    end
    hf = timed!(timer, "$(label) restricted closed-shell HF") do
        closed_shell_hf(operators_diag.H, operators_diag.V)
    end
    if !hf.converged
        error("$(label) HF did not converge after $(hf.iterations) iterations")
    end
    ed = timed!(timer, "$(label) Diag2ptle ED") do
        davidson_ed(
            operators_diag.H,
            operators_diag.V,
            hf.orbital,
            Matrix{Float64}(hf.fock),
        )
    end
    return (
        metadata = route_metadata(label, source, fixed, operators_diag),
        hf = hf,
        ed = ed,
    )
end

function print_route(label::AbstractString, row)
    meta = row.metadata
    hf = row.hf
    ed = row.ed
    println(label)
    println("  shared_layer_types = ", meta.shared_layer_types)
    println("  shared_layer_columns = ", meta.shared_layer_columns)
    println("  fixed_block_size = ", meta.fixed_block_size)
    @printf("  fixed_overlap_error = %.12e\n", meta.fixed_overlap_error)
    println("  operator_dimension = ", meta.operator_dimension)
    println("  residual_count = ", meta.residual_count)
    @printf("  operator_overlap_error = %.12e\n", meta.operator_overlap_error)
    @printf("  h_symmetry_error = %.12e\n", meta.h_symmetry_error)
    @printf("  v_symmetry_error = %.12e\n", meta.v_symmetry_error)
    println("  residual_width_status = ", meta.width_status)
    println("  residual_owner_set = ", meta.owner_set)
    println("  hf_converged = ", hf.converged)
    println("  hf_iterations = ", hf.iterations)
    @printf("  hf_total = %.12f\n", hf.total_energy)
    @printf("  hf_orbital_residual = %.12e\n", hf.orbital_residual)
    println("  ed_route = ", ed.route)
    println("  ed_reason = ", ed.reason)
    println("  ed_iterations = ", ed.iterations)
    println("  ed_solver_dimension = ", ed.solver_dimension)
    @printf("  ed_total = %.12f\n", ed.total_energy)
    @printf("  ed_residual = %.12e\n", ed.residual)
    @printf("  ed_hf_lowering_mha = %.12f\n", 1000.0 * (hf.total_energy - ed.total_energy))
    @printf("  bo_error_mha = %.12f\n", 1000.0 * (ed.total_energy - BO_TOTAL_R4))
end

function print_comparison(default, endcap)
    default_delta_ref = 1000.0 * (default.ed.total_energy - HIGH_ORDER_OLD_SP_ED_TOTAL_R4)
    endcap_delta_ref = 1000.0 * (endcap.ed.total_energy - HIGH_ORDER_ALL_SHARED_Q4_ED_TOTAL_R4)
    endcap_minus_default = 1000.0 * (endcap.ed.total_energy - default.ed.total_energy)
    dim_delta = endcap.metadata.operator_dimension - default.metadata.operator_dimension

    println("comparison")
    println("  high_order_old_sp_dim = ", HIGH_ORDER_OLD_SP_DIM_R4)
    @printf("  high_order_old_sp_ed_total = %.12f\n", HIGH_ORDER_OLD_SP_ED_TOTAL_R4)
    println("  high_order_all_shared_q4_dim = ", HIGH_ORDER_ALL_SHARED_Q4_DIM_R4)
    @printf("  high_order_all_shared_q4_ed_total = %.12f\n", HIGH_ORDER_ALL_SHARED_Q4_ED_TOTAL_R4)
    println("  mainline_dimension_delta_endcap_minus_default = ", dim_delta)
    @printf("  mainline_ed_delta_endcap_minus_default_mha = %.12f\n", endcap_minus_default)
    @printf("  mainline_default_minus_high_order_old_sp_mha = %.12f\n", default_delta_ref)
    @printf("  mainline_endcap_minus_high_order_all_shared_q4_mha = %.12f\n", endcap_delta_ref)
end

function print_timing(timer::StageTimer)
    println("timing")
    for row in timer.rows
        @printf(
            "  %s: %.6f s, %d bytes, gc %.6f s\n",
            row.label,
            row.time,
            row.bytes,
            row.gctime,
        )
    end
end

function run_study()
    timer = StageTimer()
    println("H2 endcap/panel HF/ED chemistry reproduction")
    println("  date = ", Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS"))
    println("  R_bohr = ", R_BOHR)
    println("  core_spacing = ", CORE_SPACING)
    println("  xmax_parallel = ", XMAX_PARALLEL)
    println("  xmax_transverse = ", XMAX_TRANSVERSE)
    println("  supplement = ", SUPPLEMENT_ATOM, "/", SUPPLEMENT_BASIS, " lmax=", SUPPLEMENT_LMAX)
    println("  interaction_treatment = ", INTERACTION_TREATMENT)
    println("  endcap_policy = :endcap_panel_owned, q=", ENDCAP_Q, ", L=", ENDCAP_L)
    println("  diag2ptle_path = ", DIAG2PTLE_PATH)
    @printf("  diag2ptle_top_level_load_seconds = %.6f\n", DIAG2PTLE_LOAD_SECONDS)
    println()

    expansion = timed!(timer, "coulomb expansion") do
        coulomb_gaussian_expansion(doacc = false)
    end
    basis = timed!(timer, "bond-aligned basis") do
        build_basis()
    end
    supplement = timed!(timer, "H cc-pVTZ S/P supplement") do
        legacy_bond_aligned_diatomic_gaussian_supplement(
            SUPPLEMENT_ATOM,
            SUPPLEMENT_BASIS,
            basis.nuclei;
            lmax = SUPPLEMENT_LMAX,
        )
    end

    default = run_route!(
        timer,
        :default_complete_rectangular,
        :complete_rectangular,
        basis,
        supplement,
        expansion,
    )
    endcap = run_route!(
        timer,
        :endcap_panel_owned,
        :endcap_panel_owned,
        basis,
        supplement,
        expansion,
    )

    print_route("default complete-rectangular route", default)
    print_route("endcap/panel owned route", endcap)
    print_comparison(default, endcap)
    print_timing(timer)

    return (default = default, endcap = endcap, timing = Tuple(timer.rows))
end

run_study()
