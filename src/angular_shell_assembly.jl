"""
    AtomicShellLocalInjectedAngularAssembly

Experimental atom-side angular assembly built from shell-local injected angular
bases across a list of shell radii.

This object is still angular-only. It does not yet include radial coupling,
Coulomb assembly, or any end-to-end atomic workflow.
"""
struct AtomicShellLocalInjectedAngularAssembly
    shell_radii::Vector{Float64}
    shell_ids::Vector{Int}
    shell_orders::Vector{Int}
    shell_offsets::Vector{Int}
    shell_dimensions::Vector{Int}
    shells::Vector{ShellLocalInjectedAngularBasis}
    shell_moment_blocks::Vector{Dict{Int,Matrix{Float64}}}
    shell_kinetic_moment_blocks::Vector{Dict{Int,Matrix{Float64}}}
    shell_exact_lmax::Vector{Int}
    shell_kinetic_lcap::Vector{Int}
    shell_kinetic_lexpand::Vector{Int}
    shell_kinetic_tail::Vector{Float64}
    overlap_blocks::Matrix{Matrix{Float64}}
    kinetic_blocks::Matrix{Matrix{Float64}}
    overlap::Matrix{Float64}
    kinetic::Matrix{Float64}
end

function Base.show(io::IO, assembly::AtomicShellLocalInjectedAngularAssembly)
    print(
        io,
        "AtomicShellLocalInjectedAngularAssembly(nshells=",
        length(assembly.shells),
        ", total_dim=",
        size(assembly.overlap, 1),
        ", shell_orders=",
        repr(assembly.shell_orders),
        ")",
    )
end

_smooth_step_logistic(x, width) = width <= 0 ? (x >= 0 ? 1.0 : 0.0) : inv(1 + exp(-x / width))

function _smooth_shell_order_of_radius(
    radius::Real;
    ord_min::Int,
    ord_max::Int,
    r_lo::Real,
    r_hi::Real,
    w_lo::Real,
    w_hi::Real,
)
    ord_min <= ord_max || throw(ArgumentError("ord_min must be <= ord_max"))
    radius_f = Float64(radius)
    if radius_f <= r_lo || radius_f >= r_hi
        return ord_min
    end

    left = _smooth_step_logistic(radius_f - r_lo, w_lo)
    right = _smooth_step_logistic(r_hi - radius_f, w_hi)
    window = clamp(left * right, 0.0, 1.0)
    fraction = window^2.5
    return ord_min + round(Int, fraction * (ord_max - ord_min))
end

function _promote_to_available_shell_order(order::Int, available_orders::Vector{Int})
    idx = searchsortedfirst(available_orders, order)
    idx <= length(available_orders) || throw(
        ArgumentError(
            "requested shell order $order is above the vendored point-set range; available orders are $(join(available_orders, ", "))",
        ),
    )
    return available_orders[idx]
end

"""
    assign_atomic_angular_shell_orders(shell_radii; ord_min, ord_max, r_lo=0.15, r_hi=7.0, w_lo=0.2, w_hi=0.7, available_orders=sphere_point_set_orders())

Assign experimental angular shell orders across a list of shell radii using the
first narrow smooth schedule, then promote each requested order to the next
available vendored point-set order.
"""
function assign_atomic_angular_shell_orders(
    shell_radii::AbstractVector{<:Real};
    ord_min::Int,
    ord_max::Int,
    r_lo::Real = 0.15,
    r_hi::Real = 7.0,
    w_lo::Real = 0.2,
    w_hi::Real = 0.7,
    available_orders::AbstractVector{<:Integer} = sphere_point_set_orders(),
)
    available = sort(Int[order for order in available_orders])
    isempty(available) && throw(ArgumentError("no vendored sphere-point-set orders are available"))
    ord_min ≥ available[1] || throw(ArgumentError("ord_min must be within the vendored order range"))
    ord_max ≤ available[end] || throw(ArgumentError("ord_max must be within the vendored order range"))

    assigned = Int[]
    for radius in shell_radii
        requested = _smooth_shell_order_of_radius(
            radius;
            ord_min = ord_min,
            ord_max = ord_max,
            r_lo = r_lo,
            r_hi = r_hi,
            w_lo = w_lo,
            w_hi = w_hi,
        )
        push!(assigned, _promote_to_available_shell_order(requested, available))
    end
    return assigned
end

function _ylm_rows_for_l(channels::YlmChannelSet, l::Int)
    return findall(channel -> channel.l == l, channels.channel_data)
end

function _shell_local_ylm_moment_block(shell::ShellLocalInjectedAngularBasis, l::Int)
    l >= 0 || throw(ArgumentError("_shell_local_ylm_moment_block requires l >= 0"))
    channels = YlmChannelSet(l, [YlmChannel(l, m) for m in -l:l])
    raw_overlap, _ = _ylm_prototype_coupling(
        shell.point_set.coordinates,
        shell.kappa,
        Matrix{Float64}(I, shell.prototype_count, shell.prototype_count),
        channels,
    )
    gaussian_part = raw_overlap * shell.prototype_coefficients

    y_rows = zeros(Float64, 2 * l + 1, shell.final_count)
    injected_rows = _ylm_rows_for_l(shell.injected_channels, l)
    if !isempty(injected_rows)
        y_rows[:, :] .= shell.ylm_coefficients[injected_rows, :]
    end
    return y_rows .+ gaussian_part
end

function _shell_raw_overlap(shell_a::ShellLocalInjectedAngularBasis, shell_b::ShellLocalInjectedAngularBasis)
    channels_a = shell_a.injected_channels
    channels_b = shell_b.injected_channels
    ny_a = shell_a.injected_count
    ny_b = shell_b.injected_count
    ng_a = shell_a.prototype_count
    ng_b = shell_b.prototype_count

    y_y = zeros(Float64, ny_a, ny_b)
    lookup_b = Dict{YlmChannel,Int}(channels_b.channel_data[i] => i for i in eachindex(channels_b.channel_data))
    for i in eachindex(channels_a.channel_data)
        j = get(lookup_b, channels_a[i], 0)
        j == 0 || (y_y[i, j] = 1.0)
    end

    y_g_a, _ = _ylm_prototype_coupling(
        shell_b.point_set.coordinates,
        shell_b.kappa,
        Matrix{Float64}(I, ng_b, ng_b),
        channels_a,
    )
    y_g_b, _ = _ylm_prototype_coupling(
        shell_a.point_set.coordinates,
        shell_a.kappa,
        Matrix{Float64}(I, ng_a, ng_a),
        channels_b,
    )

    g_g = Matrix{Float64}(undef, ng_a, ng_b)
    for j in 1:ng_b, i in 1:ng_a
        g_g[i, j] = _prototype_overlap(shell_a.point_set.coordinates, shell_a.kappa, i, i)
    end
    for j in 1:ng_b, i in 1:ng_a
        ki = shell_a.kappa[i]
        kj = shell_b.kappa[j]
        c = clamp(dot(view(shell_a.point_set.coordinates, i, :), view(shell_b.point_set.coordinates, j, :)), -1.0, 1.0)
        radius_squared = (ki - kj)^2 + 2.0 * ki * kj * (1.0 + c)
        radius = sqrt(max(radius_squared, 0.0))
        if radius < _angular_small_argument_cutoff()
            sinhc_series = 1.0 + radius_squared / 6.0 + (radius_squared^2) / 120.0
            g_g[i, j] = 4 * pi * exp(-(ki + kj)) * sinhc_series
        else
            tail = radius - (ki + kj)
            one_minus_tail = -expm1(-2 * radius)
            g_g[i, j] = 4 * pi * 0.5 * exp(tail) * (one_minus_tail / radius)
        end
    end

    raw = zeros(Float64, ny_a + ng_a, ny_b + ng_b)
    raw[1:ny_a, 1:ny_b] .= y_y
    raw[1:ny_a, ny_b+1:end] .= y_g_a
    raw[ny_a+1:end, 1:ny_b] .= transpose(y_g_b)
    raw[ny_a+1:end, ny_b+1:end] .= g_g
    return raw
end

function _shell_local_overlap(shell_a::ShellLocalInjectedAngularBasis, shell_b::ShellLocalInjectedAngularBasis)
    raw_overlap = _shell_raw_overlap(shell_a, shell_b)
    coeff_a = vcat(shell_a.ylm_coefficients, shell_a.prototype_coefficients)
    coeff_b = vcat(shell_b.ylm_coefficients, shell_b.prototype_coefficients)
    return transpose(coeff_a) * raw_overlap * coeff_b
end

function _shell_local_pair_kinetic(
    shell_a::ShellLocalInjectedAngularBasis,
    shell_b::ShellLocalInjectedAngularBasis,
    blocks_a::Dict{Int,Matrix{Float64}},
    blocks_b::Dict{Int,Matrix{Float64}},
    lcap_a::Int,
    lcap_b::Int,
)
    pair_lmax = min(lcap_a, lcap_b)
    kinetic = zeros(Float64, shell_a.final_count, shell_b.final_count)
    for l in 1:pair_lmax
        kinetic .+= _kinetic_eigenvalue(YlmChannel(l, 0)) .* (transpose(blocks_a[l]) * blocks_b[l])
    end
    return kinetic
end

_shell_ylm_residual_tol() = 1.0e-12
_shell_ylm_span_cap() = 256
_shell_ylm_small_streak() = 2
_shell_ylm_extra_required() = 4

function _ylm_increment!(captured::AbstractVector{<:Real}, mat::AbstractMatrix{<:Real})
    contrib = vec(sum(abs2, mat; dims = 1))
    captured .+= contrib
    return sqrt(maximum(contrib))
end

function _ylm_tail_max(diag_norms::AbstractVector{<:Real}, captured::AbstractVector{<:Real})
    tail = 0.0
    for i in eachindex(diag_norms)
        diff = Float64(diag_norms[i]) - Float64(captured[i])
        diff = diff > 0 ? diff : 0.0
        tail = max(tail, sqrt(diff))
    end
    return tail
end

function _build_shell_local_kinetic_moment_blocks(
    shell::ShellLocalInjectedAngularBasis;
    residual_tol::Real = _shell_ylm_residual_tol(),
    span_cap::Int = _shell_ylm_span_cap(),
    small_streak_required::Int = _shell_ylm_small_streak(),
    lextra::Int = _shell_ylm_extra_required(),
)
    residual_tol > 0 || throw(ArgumentError("residual_tol must be positive"))
    span_cap >= 0 || throw(ArgumentError("span_cap must be nonnegative"))
    small_streak_required >= 1 ||
        throw(ArgumentError("small_streak_required must be at least one"))
    lextra >= 0 || throw(ArgumentError("lextra must be nonnegative"))

    diag_norms = Float64[shell.final_overlap[i, i] for i in 1:shell.final_count]
    captured = zeros(Float64, shell.final_count)
    required_min = max(0, shell.l_inject + lextra)
    blocks = Dict{Int,Matrix{Float64}}()
    lcap = -1
    lexpand = 0
    small_streak = 0
    tail = _ylm_tail_max(diag_norms, captured)

    for l in 0:span_cap
        mat = _shell_local_ylm_moment_block(shell, l)
        blocks[l] = mat
        delta = _ylm_increment!(captured, mat)
        lcap = l
        if delta > residual_tol
            lexpand = l
        end
        tail = _ylm_tail_max(diag_norms, captured)
        if l >= required_min && delta <= residual_tol
            small_streak += 1
        else
            small_streak = 0
        end
        small_streak >= small_streak_required && break
    end

    lcap >= 0 || error("_build_shell_local_kinetic_moment_blocks failed to build any shell moments")
    return (blocks = blocks, lcap = lcap, lexpand = max(lexpand, required_min), tail = tail)
end

function _assemble_shell_block_matrix(blocks::Matrix{Matrix{Float64}}, offsets::Vector{Int}, dims::Vector{Int})
    total_dim = sum(dims)
    full = zeros(Float64, total_dim, total_dim)
    nshells = length(dims)
    for a in 1:nshells
        ia = offsets[a]:(offsets[a] + dims[a] - 1)
        for b in 1:nshells
            ib = offsets[b]:(offsets[b] + dims[b] - 1)
            full[ia, ib] .= blocks[a, b]
        end
    end
    return full
end

"""
    build_atomic_shell_local_angular_assembly(shell_radii; shell_orders=nothing, beta=2.0, l_inject=:auto, tau=1e-12, whiten=:svd, ord_min=minimum(sphere_point_set_orders()), ord_max=maximum(sphere_point_set_orders()), r_lo=0.15, r_hi=7.0, w_lo=0.2, w_hi=0.7)

Build the first experimental atom-side angular assembly layer by assigning one
shell-local injected angular basis to each shell radius and assembling exact
pairwise overlap plus injected-sector angular kinetic blocks.
"""
function build_atomic_shell_local_angular_assembly(
    shell_radii::AbstractVector{<:Real};
    shell_orders::Union{Nothing,AbstractVector{<:Integer}} = nothing,
    beta::Real = 2.0,
    l_inject::Union{Int,Symbol} = :auto,
    tau::Real = 1.0e-12,
    whiten::Symbol = :svd,
    ord_min::Int = minimum(sphere_point_set_orders()),
    ord_max::Int = maximum(sphere_point_set_orders()),
    r_lo::Real = 0.15,
    r_hi::Real = 7.0,
    w_lo::Real = 0.2,
    w_hi::Real = 0.7,
)
    radii = Float64[shell_radii[i] for i in eachindex(shell_radii)]
    nshells = length(radii)
    nshells ≥ 1 || throw(ArgumentError("build_atomic_shell_local_angular_assembly requires at least one shell radius"))

    assigned_orders =
        shell_orders === nothing ?
        assign_atomic_angular_shell_orders(
            radii;
            ord_min = ord_min,
            ord_max = ord_max,
            r_lo = r_lo,
            r_hi = r_hi,
            w_lo = w_lo,
            w_hi = w_hi,
        ) :
        Int[order for order in shell_orders]
    length(assigned_orders) == nshells ||
        throw(ArgumentError("shell_orders length must match the number of shell radii"))

    cached_shells = Dict{Tuple{Int,Float64,Union{Int,Symbol},Float64,Symbol},ShellLocalInjectedAngularBasis}()
    shells = Vector{ShellLocalInjectedAngularBasis}(undef, nshells)
    shell_moment_blocks = Vector{Dict{Int,Matrix{Float64}}}(undef, nshells)
    cached_kinetic_moment_blocks = Dict{
        Tuple{Int,Float64,Union{Int,Symbol},Float64,Symbol},
        NamedTuple{(:blocks, :lcap, :lexpand, :tail),Tuple{Dict{Int,Matrix{Float64}},Int,Int,Float64}},
    }()
    shell_kinetic_moment_blocks = Vector{Dict{Int,Matrix{Float64}}}(undef, nshells)
    shell_exact_lmax = Vector{Int}(undef, nshells)
    shell_kinetic_lcap = Vector{Int}(undef, nshells)
    shell_kinetic_lexpand = Vector{Int}(undef, nshells)
    shell_kinetic_tail = Vector{Float64}(undef, nshells)

    for i in 1:nshells
        key = (assigned_orders[i], Float64(beta), l_inject, Float64(tau), whiten)
        shell = get!(cached_shells, key) do
            build_shell_local_injected_angular_basis(
                assigned_orders[i];
                beta = beta,
                l_inject = l_inject,
                tau = tau,
                whiten = whiten,
            )
        end
        shells[i] = shell
        moment_blocks = Dict{Int,Matrix{Float64}}()
        for l in 0:shell.l_inject
            moment_blocks[l] = _shell_local_ylm_moment_block(shell, l)
        end
        shell_moment_blocks[i] = moment_blocks
        # The shell-local exact injection only fixes the low-l subspace. The
        # one-body angular kinetic needs a larger Y-moment span for the mixed
        # remainder, or spurious extra s-like states collapse onto the 1s branch.
        kinetic_moment_data = get!(cached_kinetic_moment_blocks, key) do
            _build_shell_local_kinetic_moment_blocks(shell)
        end
        shell_kinetic_moment_blocks[i] = kinetic_moment_data.blocks
        shell_exact_lmax[i] = shell.l_inject
        shell_kinetic_lcap[i] = kinetic_moment_data.lcap
        shell_kinetic_lexpand[i] = kinetic_moment_data.lexpand
        shell_kinetic_tail[i] = kinetic_moment_data.tail
    end

    shell_dimensions = [shell.final_count for shell in shells]
    shell_offsets = cumsum(vcat(1, shell_dimensions[1:end-1]))
    overlap_blocks = Matrix{Matrix{Float64}}(undef, nshells, nshells)
    kinetic_blocks = Matrix{Matrix{Float64}}(undef, nshells, nshells)

    for a in 1:nshells
        for b in 1:nshells
            overlap_blocks[a, b] = _shell_local_overlap(shells[a], shells[b])
            kinetic_blocks[a, b] = _shell_local_pair_kinetic(
                shells[a],
                shells[b],
                shell_kinetic_moment_blocks[a],
                shell_kinetic_moment_blocks[b],
                shell_kinetic_lcap[a],
                shell_kinetic_lcap[b],
            )
        end
    end

    overlap = _assemble_shell_block_matrix(overlap_blocks, shell_offsets, shell_dimensions)
    overlap = 0.5 .* (overlap .+ transpose(overlap))
    kinetic = _assemble_shell_block_matrix(kinetic_blocks, shell_offsets, shell_dimensions)
    kinetic = 0.5 .* (kinetic .+ transpose(kinetic))

    return AtomicShellLocalInjectedAngularAssembly(
        radii,
        collect(1:nshells),
        assigned_orders,
        shell_offsets,
        shell_dimensions,
        shells,
        shell_moment_blocks,
        shell_kinetic_moment_blocks,
        shell_exact_lmax,
        shell_kinetic_lcap,
        shell_kinetic_lexpand,
        shell_kinetic_tail,
        overlap_blocks,
        kinetic_blocks,
        overlap,
        kinetic,
    )
end

"""
    atomic_shell_local_angular_diagnostics(assembly::AtomicShellLocalInjectedAngularAssembly)

Return compact diagnostics for the first experimental shell-to-atom angular
assembly layer.
"""
function atomic_shell_local_angular_diagnostics(assembly::AtomicShellLocalInjectedAngularAssembly)
    nshells = length(assembly.shells)
    shell_diags = [shell_local_injected_angular_diagnostics(shell) for shell in assembly.shells]
    diagonal_overlap_errors = Float64[]
    diagonal_kinetic_errors = Float64[]
    pair_overlap_symmetry = Float64[]
    pair_kinetic_symmetry = Float64[]

    for a in 1:nshells
        ia = assembly.shell_offsets[a]:(assembly.shell_offsets[a] + assembly.shell_dimensions[a] - 1)
        push!(diagonal_overlap_errors, opnorm(assembly.overlap[ia, ia] - Matrix(I, assembly.shell_dimensions[a], assembly.shell_dimensions[a]), Inf))

        expected = Diagonal([_kinetic_eigenvalue(channel) for channel in assembly.shells[a].injected_channels])
        actual = assembly.shells[a].injected_kinetic
        push!(diagonal_kinetic_errors, opnorm(actual - expected, Inf))

        for b in a:nshells
            push!(pair_overlap_symmetry, opnorm(assembly.overlap_blocks[a, b] - transpose(assembly.overlap_blocks[b, a]), Inf))
            push!(pair_kinetic_symmetry, opnorm(assembly.kinetic_blocks[a, b] - transpose(assembly.kinetic_blocks[b, a]), Inf))
        end
    end

    return (
        nshells = nshells,
        total_dim = size(assembly.overlap, 1),
        shell_orders = copy(assembly.shell_orders),
        shell_dimensions = copy(assembly.shell_dimensions),
        shell_exact_lmax = copy(assembly.shell_exact_lmax),
        shell_kinetic_lcap = copy(assembly.shell_kinetic_lcap),
        shell_kinetic_lexpand = copy(assembly.shell_kinetic_lexpand),
        max_shell_kinetic_tail = maximum(assembly.shell_kinetic_tail),
        max_shell_overlap_error = maximum(diag.overlap_error for diag in shell_diags),
        max_shell_injected_exactness_error = maximum(diag.injected_exactness_error for diag in shell_diags),
        max_shell_injected_kinetic_error = maximum(diag.injected_kinetic_error for diag in shell_diags),
        max_diagonal_overlap_error = maximum(diagonal_overlap_errors),
        max_diagonal_injected_kinetic_error = maximum(diagonal_kinetic_errors),
        max_pair_overlap_symmetry_error = maximum(pair_overlap_symmetry),
        max_pair_kinetic_symmetry_error = maximum(pair_kinetic_symmetry),
    )
end
