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
    profiles::Vector{ShellLocalAngularProfile}
    shells::Vector{ShellLocalInjectedAngularBasis}
    shell_moment_blocks::Vector{Dict{Int,Matrix{Float64}}}
    shell_kinetic_moment_blocks::Vector{Dict{Int,Matrix{Float64}}}
    shell_interaction_moment_blocks::Vector{Dict{Int,Matrix{Float64}}}
    shell_exact_lmax::Vector{Int}
    shell_kinetic_lcap::Vector{Int}
    shell_kinetic_lexpand::Vector{Int}
    shell_kinetic_tail::Vector{Float64}
    shell_interaction_lcap::Vector{Int}
    shell_interaction_lexpand::Vector{Int}
    shell_interaction_tail::Vector{Float64}
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

"""
    atomic_contact_integrals(rb, assembly; radial_grid, angular_quadrature,
        supplements=nothing, supplement_labels=nothing, provenance=nothing)
    atomic_contact_integrals(rb, shells, ranges; labels, kwargs...)

Exact real-function contact geometry on explicit quadrature, not an IDA interaction.
Returns `(factors, site_map, diagnostics, provenance)` in the supplied final order.
The lower entrance uses actual angular bases and contiguous base-site ranges.
Supplement input `(grid, reduced_values, channels)` supplies final phased reduced
radial columns and one YlmChannel per appended site, with explicit labels.
No function reconstruction, normalization or phase alignment occurs. Volume
weights are squared signed integrals; zero weights are legal (and diagnostic
only for supplements). Accuracy and support coverage remain caller obligations.
"""
function atomic_contact_integrals(rb::RadialBasis, assembly::AtomicShellLocalInjectedAngularAssembly; kwargs...)
    ranges = [o:o+d-1 for (o, d) in zip(assembly.shell_offsets, assembly.shell_dimensions)]
    return atomic_contact_integrals(rb, assembly.shells, ranges;
        labels=reduce(vcat, [p.labels for p in assembly.profiles]; init=String[]), kwargs...)
end

function atomic_contact_integrals(rb::RadialBasis, shells::AbstractVector, ranges::AbstractVector;
    radial_grid::RadialQuadratureGrid, angular_quadrature, labels,
    supplements=nothing, supplement_labels=nothing, provenance=nothing)
    r, w = _validate_radial_operator_grid(rb, radial_grid)
    all(isfinite, r) && all(>(0), r) && all(isfinite, w) ||
        throw(ArgumentError("contact radial points and weights must be finite and positive"))
    xyz, wa = angular_quadrature.directions, angular_quadrature.weights
    wa isa AbstractVector{<:Real} && size(xyz, 1) == length(wa) > 0 && all(x -> isfinite(x) && x > 0, wa) ||
        throw(ArgumentError("sphere weights must be finite, positive and match directions"))
    length(shells) == length(ranges) == length(rb) > 0 ||
        throw(DimensionMismatch("one actual angular shell and range per radial function required"))
    offset = 1
    for (shell, range) in zip(shells, ranges)
        shell isa ShellLocalInjectedAngularBasis || throw(ArgumentError("actual angular bases required"))
        range == offset:offset+shell.final_count-1 || throw(ArgumentError("base ranges must be contiguous and ordered"))
        offset += shell.final_count
    end
    length(labels) == offset-1 && all(x -> x isa AbstractString, labels) ||
        throw(ArgumentError("explicit labels must match the base site order"))
    source = provenance === nothing ? (producer_commit="", basis_manifest_sha256="") : provenance
    source isa NamedTuple && propertynames(source) == (:producer_commit, :basis_manifest_sha256) &&
        all(x -> x isa AbstractString, Base.values(source)) || throw(ArgumentError("invalid contact provenance"))
    supplements === nothing || (supplements isa NamedTuple &&
        propertynames(supplements) == (:grid, :reduced_values, :channels)) || throw(ArgumentError("invalid supplement record"))
    ns = supplements === nothing ? 0 : length(supplements.channels)
    if supplements !== nothing
        quadrature_points(supplements.grid) == r && quadrature_weights(supplements.grid) == w ||
            throw(ArgumentError("supplement samples must use exactly the contact radial grid"))
        supplements.reduced_values isa AbstractMatrix{<:Real} &&
            size(supplements.reduced_values) == (length(r), ns) && all(isfinite, supplements.reduced_values) ||
            throw(ArgumentError("invalid final reduced-radial supplement samples"))
        all(ch -> ch isa YlmChannel, supplements.channels) || throw(ArgumentError("Ylm channels required"))
        supplement_labels !== nothing && length(supplement_labels) == ns &&
            all(x -> x isa AbstractString, supplement_labels) || throw(ArgumentError("explicit supplement labels required"))
    elseif supplement_labels !== nothing
        throw(ArgumentError("supplement labels require supplement functions"))
    end
    angular_values, unique_shells = Matrix{Float64}[], ShellLocalInjectedAngularBasis[]
    profile_index = Int[]
    values = @timeg "contact function evaluation" begin
        for shell in shells
            p = findfirst(s -> s === shell, unique_shells)
            if p === nothing
                push!(unique_shells, shell)
                push!(angular_values, evaluate_angular_basis(shell, xyz))
                p = length(angular_values)
            end
            push!(profile_index, p)
        end
        channels = Dict{YlmChannel,Int}()
        for ch in (supplements === nothing ? YlmChannel[] : supplements.channels)
            p = get!(channels, ch) do
                push!(angular_values, reshape([_real_spherical_harmonic(ch, view(xyz, i, :)) for i in axes(xyz, 1)], :, 1))
                length(angular_values)
            end
            push!(profile_index, p)
        end
        base = _basis_values_matrix(rb, r)
        ns == 0 ? base : hcat(base, supplements.reduced_values)
    end
    radial, radial_overlap = @timeg "contact radial contractions" begin
        squared = values.^2
        R = transpose(squared) * ((w ./ r.^2) .* squared)
        Tr = transpose(values.^3) * ((w ./ r.^2) .* values)
        ((R=R, Tr=Tr, wrad=(transpose(values) * (w .* r)).^2), transpose(values) * (w .* values))
    end
    A, Ta, angular_overlap = Matrix{Float64}[], Matrix{Float64}[], Matrix{Float64}[]
    np = length(angular_values)
    pair_index = zeros(Int, np, np)
    wang = Vector{Float64}[]
    @timeg "contact angular pair contractions" begin
        for F in angular_values
            push!(wang, (transpose(F) * wa).^2)
        end
        for q in 1:np, p in 1:np
            F, G = angular_values[p], angular_values[q]
            push!(A, transpose(F.^2) * (wa .* G.^2))
            push!(Ta, transpose(F.^3) * (wa .* G))
            push!(angular_overlap, transpose(F) * (wa .* G))
            pair_index[p, q] = length(A)
        end
    end
    all(isfinite, radial.R) && all(isfinite, radial.Tr) && all(isfinite, radial.wrad) &&
        all(M -> all(isfinite, M), A) && all(M -> all(isfinite, M), Ta) &&
        all(v -> all(isfinite, v), wang) || throw(ArgumentError("nonfinite contact integrals"))
    site_labels = ns == 0 ? String.(labels) : vcat(String.(labels), String.(supplement_labels))
    site_map = [(group=g, local_index=k, label="") for g in eachindex(profile_index)
        for k in axes(angular_values[profile_index[g]], 2)]
    site_map = [(group=s.group, local_index=s.local_index, label=site_labels[i]) for (i, s) in enumerate(site_map)]
    row_sums = zeros(length(site_map))
    full_error, cursor = 0.0, 1
    @timeg "contact overlap and row diagnostics" for g in eachindex(profile_index)
        p = profile_index[g]
        rows = zeros(length(wang[p]))
        for h in eachindex(profile_index)
            q = profile_index[h]
            pair = pair_index[p, q]
            block = radial_overlap[g, h] .* angular_overlap[pair]
            if g == h
                block -= I
            end
            full_error = max(full_error, maximum(abs, block))
            rows .+= radial.R[g, h] .* vec(sum(A[pair]; dims=2))
        end
        count = length(rows)
        row_sums[cursor:cursor+count-1] .= radial.wrad[g] .* wang[p] .* rows
        cursor += count
    end
    overlap = (radial=maximum(abs, radial_overlap[1:length(rb), 1:length(rb)]-I),
        angular=maximum(maximum(abs, angular_overlap[pair_index[p, p]]-I) for p in 1:np), full=full_error)
    factors = (radial=radial, angular=(A=A, Ta=Ta, wang=wang, profile_index=profile_index, pair_index=pair_index))
    diagnostics = (overlap=overlap, row_sums=row_sums,
        quadrature=(radial_nodes=length(r), sphere_nodes=length(wa), radial_extent=extrema(r)))
    return (factors=factors, site_map=site_map, diagnostics=diagnostics, provenance=source)
end

"""
    materialize_atomic_contact_integrals(result)

Derive dense `(U,T,w)` from atomic contact factors, including every appended
supplement. T remains signed and directed; no quadrature is repeated.
"""
function materialize_atomic_contact_integrals(result)
    radial, angular, sites = result.factors.radial, result.factors.angular, result.site_map
    n = length(sites)
    U, T, w = zeros(n, n), zeros(n, n), zeros(n)
    @timeg "contact dense materialization" for (i, a) in enumerate(sites)
        p = angular.profile_index[a.group]
        w[i] = radial.wrad[a.group] * angular.wang[p][a.local_index]
        for (j, b) in enumerate(sites)
            pair = angular.pair_index[p, angular.profile_index[b.group]]
            U[i, j] = radial.R[a.group, b.group] * angular.A[pair][a.local_index, b.local_index]
            T[i, j] = radial.Tr[a.group, b.group] * angular.Ta[pair][a.local_index, b.local_index]
        end
    end
    return (U=U, T=T, w=w)
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

"""
    assign_atomic_angular_shell_orders(radial_ops::RadialAtomicOperators;
                                      ord_max, estimated_energy_cutoff, required_l=1)

Experimental inner-shell point counts from one estimated angular-correlation
budget in Ha. A correlated radial s-sector reference and projected IDA couplings
estimate discarded channels; this is not a guaranteed reconstructed/continuum
energy bound. Only the first `min(12,n)` shells may shrink, retaining at least
`required_l` (s/p by default). Outer shells retain `ord_max`.

Supports 1:40 radial functions and vendored 10:100 point rules, with
`IntegralDiagonal` radial operators, inverse-r2 and all required multipoles
(stored or available from retained samples). The bounded dense s-reference
uses O(n^4) storage. Pass the returned `Vector{Int}` as `shell_orders` to the
normal builder; keep `interaction_lmax=:auto`. Radius/default paths are unchanged.
"""
function assign_atomic_angular_shell_orders(
    radial_ops::RadialAtomicOperators; ord_max, estimated_energy_cutoff, required_l=1,
)
    ord_max isa Integer && !(ord_max isa Bool) && 10 <= ord_max <= 100 &&
        ord_max in sphere_point_set_orders() || throw(ArgumentError("ord_max must be a vendored 10:100 point count"))
    estimated_energy_cutoff isa Real && !(estimated_energy_cutoff isa Bool) &&
        isfinite(estimated_energy_cutoff) && estimated_energy_cutoff > 0 ||
        throw(ArgumentError("estimated_energy_cutoff must be finite and positive (Ha)"))
    cutoff = Float64(estimated_energy_cutoff)
    isfinite(cutoff) && cutoff > 0 || throw(ArgumentError("cutoff must be representable as positive finite Float64"))
    lmax = _choose_l_inject(Int(ord_max))
    required_l isa Integer && !(required_l isa Bool) && 1 <= required_l <= lmax ||
        throw(ArgumentError("required_l must lie between 1 and the maximum rule's auto-injected l"))
    radii = radial_ops.shell_centers_r
    n = length(radii)
    1 <= n <= 40 || throw(ArgumentError("estimated allocation supports 1:40 radial functions"))
    all(r -> isfinite(r) && r > 0, radii) && all(i -> radii[i] < radii[i+1], 1:n-1) ||
        throw(ArgumentError("shell radii must be finite, positive and strictly ordered"))
    radial_ops.approximation isa IntegralDiagonal || throw(ArgumentError("estimated allocation requires IntegralDiagonal"))
    length(radial_ops.centrifugal_data) >= 2 || throw(ArgumentError("inverse-r2 data are required"))
    inverse_r2 = centrifugal(radial_ops, 1)
    for matrix in (radial_ops.overlap, radial_ops.kinetic, radial_ops.nuclear, inverse_r2)
        size(matrix) == (n, n) && all(isfinite, matrix) || throw(ArgumentError("radial matrices must be finite and match the shell count"))
    end
    shell = shell_local_angular_profile(Int(ord_max)).basis
    moments = _build_shell_local_interaction_moment_blocks(shell)
    moments.lcap <= _stored_multipole_lmax(radial_ops) || _radial_multipoles_extendable(radial_ops) ||
        throw(ArgumentError("radial multipoles must cover the maximum profile's moment span"))
    couplings = _angular_correlation_couplings(radial_ops, shell, moments.blocks)
    weights = _angular_correlation_weights(radial_ops.kinetic + radial_ops.nuclear, inverse_r2, couplings)
    keeps = _angular_correlation_keeps(weights, Int(required_l), cutoff)
    available = filter(order -> 10 <= order <= ord_max, sphere_point_set_orders())
    counts = [l == lmax ? Int(ord_max) : first(order for order in available if _choose_l_inject(order) >= l) for l in 1:lmax]
    return [i <= min(12, n) ? counts[keeps[i]] : Int(ord_max) for i in 1:n]
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
_shell_ylm_kinetic_extra_required() = 4
_shell_ylm_interaction_extra_required() = 4

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

function _build_shell_local_expanded_moment_blocks(
    shell::ShellLocalInjectedAngularBasis;
    residual_tol::Real = _shell_ylm_residual_tol(),
    span_cap::Int = _shell_ylm_span_cap(),
    small_streak_required::Int = _shell_ylm_small_streak(),
    required_min::Int,
)
    residual_tol > 0 || throw(ArgumentError("residual_tol must be positive"))
    span_cap >= 0 || throw(ArgumentError("span_cap must be nonnegative"))
    small_streak_required >= 1 ||
        throw(ArgumentError("small_streak_required must be at least one"))
    required_min >= 0 || throw(ArgumentError("required_min must be nonnegative"))

    diag_norms = Float64[shell.final_overlap[i, i] for i in 1:shell.final_count]
    captured = zeros(Float64, shell.final_count)
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

function _build_shell_local_kinetic_moment_blocks(
    shell::ShellLocalInjectedAngularBasis;
    residual_tol::Real = _shell_ylm_residual_tol(),
    span_cap::Int = _shell_ylm_span_cap(),
    small_streak_required::Int = _shell_ylm_small_streak(),
    lextra::Int = _shell_ylm_kinetic_extra_required(),
)
    lextra >= 0 || throw(ArgumentError("lextra must be nonnegative"))
    required_min = max(0, shell.l_inject + lextra)
    return _build_shell_local_expanded_moment_blocks(
        shell;
        residual_tol = residual_tol,
        span_cap = span_cap,
        small_streak_required = small_streak_required,
        required_min = required_min,
    )
end

function _build_shell_local_interaction_moment_blocks(
    shell::ShellLocalInjectedAngularBasis;
    residual_tol::Real = _shell_ylm_residual_tol(),
    span_cap::Int = _shell_ylm_span_cap(),
    small_streak_required::Int = _shell_ylm_small_streak(),
    lextra::Int = _shell_ylm_interaction_extra_required(),
)
    lextra >= 0 || throw(ArgumentError("lextra must be nonnegative"))
    required_min = max(2 * shell.l_inject, shell.l_inject + lextra)
    return _build_shell_local_expanded_moment_blocks(
        shell;
        residual_tol = residual_tol,
        span_cap = span_cap,
        small_streak_required = small_streak_required,
        required_min = required_min,
    )
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

    profiles = @timeg "angular profile preparation" [shell_local_angular_profile(
        order; beta=beta, l_inject=l_inject, tau=tau, whiten=whiten,
    ) for order in assigned_orders]
    return _assemble_shell_local_angular_profiles(radii, profiles; assigned_orders=assigned_orders)
end

"""
    build_atomic_shell_local_angular_assembly(shell_radii, profiles)

Assemble one explicit `ShellLocalAngularProfile` per radial shell, preserving
its actual basis, order and gauge. Input containers are copied; profile entries
and repeated-key moment/pair arrays are shared read-only. Reported shell orders
are final function counts. No family or radial schedule is inferred.

Reuses expanded tables per profile and overlap/kinetic blocks per ordered key
pair within this assembly. Final dense matrices still require O(D^2) storage.
"""
function build_atomic_shell_local_angular_assembly(
    shell_radii::AbstractVector{<:Real},
    profiles::AbstractVector{ShellLocalAngularProfile},
)
    prof_vec = ShellLocalAngularProfile[profiles[i] for i in eachindex(profiles)]
    return _assemble_shell_local_angular_profiles(shell_radii, prof_vec)
end

struct _CachedShellMomentData
    blocks::Dict{Int,Matrix{Float64}}
    lcap::Int
    lexpand::Int
    tail::Float64
end

function _assemble_shell_local_angular_profiles(
    shell_radii::AbstractVector{<:Real}, profiles::Vector{ShellLocalAngularProfile};
    assigned_orders::Vector{Int}=Int[p.basis.final_count for p in profiles],
)
    radii = Float64[shell_radii[i] for i in eachindex(shell_radii)]
    nshells = length(radii)
    nshells >= 1 || throw(ArgumentError("angular assembly requires at least one shell radius"))
    length(profiles) == length(assigned_orders) == nshells ||
        throw(ArgumentError("profiles and assigned_orders lengths must match shell_radii"))
    shells = Vector{ShellLocalInjectedAngularBasis}(undef, nshells)
    shell_moment_blocks = Vector{Dict{Int,Matrix{Float64}}}(undef, nshells)
    cached_kinetic_moment_blocks = Dict{ShellLocalAngularProfileKey,_CachedShellMomentData}()
    shell_kinetic_moment_blocks = Vector{Dict{Int,Matrix{Float64}}}(undef, nshells)
    cached_interaction_moment_blocks = Dict{ShellLocalAngularProfileKey,_CachedShellMomentData}()
    shell_interaction_moment_blocks = Vector{Dict{Int,Matrix{Float64}}}(undef, nshells)
    shell_exact_lmax = Vector{Int}(undef, nshells)
    shell_kinetic_lcap = Vector{Int}(undef, nshells)
    shell_kinetic_lexpand = Vector{Int}(undef, nshells)
    shell_kinetic_tail = Vector{Float64}(undef, nshells)
    shell_interaction_lcap = Vector{Int}(undef, nshells)
    shell_interaction_lexpand = Vector{Int}(undef, nshells)
    shell_interaction_tail = Vector{Float64}(undef, nshells)

    @timeg "angular moment tables" for i in 1:nshells
        profile = profiles[i]
        shell = profile.basis
        shells[i] = shell
        key = profile.key
        moment_blocks = Dict{Int,Matrix{Float64}}()
        for l in 0:shell.l_inject
            moment_blocks[l] = _shell_local_ylm_moment_block(shell, l)
        end
        shell_moment_blocks[i] = moment_blocks
        # The shell-local exact injection only fixes the low-l subspace. The
        # one-body angular kinetic needs a larger Y-moment span for the mixed
        # remainder, or spurious extra s-like states collapse onto the 1s branch.
        kinetic_moment_data = get!(cached_kinetic_moment_blocks, key) do
            raw = _build_shell_local_kinetic_moment_blocks(shell)
            _CachedShellMomentData(raw.blocks, raw.lcap, raw.lexpand, raw.tail)
        end
        # The interaction needs at least the product span of the exact injected
        # sector, but for the mixed shell-local basis this is only a lower
        # bound. Follow the legacy sphgatomps principle and grow the shell
        # Y-moment tables adaptively beyond 2*l_inject.
        interaction_moment_data = get!(cached_interaction_moment_blocks, key) do
            raw = _build_shell_local_interaction_moment_blocks(shell)
            _CachedShellMomentData(raw.blocks, raw.lcap, raw.lexpand, raw.tail)
        end
        shell_kinetic_moment_blocks[i] = kinetic_moment_data.blocks
        shell_interaction_moment_blocks[i] = interaction_moment_data.blocks
        shell_exact_lmax[i] = shell.l_inject
        shell_kinetic_lcap[i] = kinetic_moment_data.lcap
        shell_kinetic_lexpand[i] = kinetic_moment_data.lexpand
        shell_kinetic_tail[i] = kinetic_moment_data.tail
        shell_interaction_lcap[i] = interaction_moment_data.lcap
        shell_interaction_lexpand[i] = interaction_moment_data.lexpand
        shell_interaction_tail[i] = interaction_moment_data.tail
    end

    shell_dimensions = [shell.final_count for shell in shells]
    shell_offsets = cumsum(vcat(1, shell_dimensions[1:end-1]))
    overlap_blocks = Matrix{Matrix{Float64}}(undef, nshells, nshells)
    kinetic_blocks = Matrix{Matrix{Float64}}(undef, nshells, nshells)
    pair_ov_cache = Dict{Tuple{ShellLocalAngularProfileKey,ShellLocalAngularProfileKey},Matrix{Float64}}()
    pair_kin_cache = Dict{Tuple{ShellLocalAngularProfileKey,ShellLocalAngularProfileKey},Matrix{Float64}}()

    @timeg "angular ordered pair blocks" for a in 1:nshells, b in 1:nshells
        key = (profiles[a].key, profiles[b].key)
        overlap_blocks[a, b] = get!(pair_ov_cache, key) do
            _shell_local_overlap(shells[a], shells[b])
        end
        kinetic_blocks[a, b] = get!(pair_kin_cache, key) do
            _shell_local_pair_kinetic(
                shells[a],
                shells[b],
                shell_kinetic_moment_blocks[a],
                shell_kinetic_moment_blocks[b],
                shell_kinetic_lcap[a],
                shell_kinetic_lcap[b],
            )
        end
    end

    overlap, kinetic = @timeg "angular final dense assembly" begin
        overlap = _assemble_shell_block_matrix(overlap_blocks, shell_offsets, shell_dimensions)
        overlap = 0.5 .* (overlap .+ transpose(overlap))
        kinetic = _assemble_shell_block_matrix(kinetic_blocks, shell_offsets, shell_dimensions)
        kinetic = 0.5 .* (kinetic .+ transpose(kinetic))
        (overlap, kinetic)
    end

    return AtomicShellLocalInjectedAngularAssembly(
        radii,
        collect(1:nshells),
        assigned_orders,
        shell_offsets,
        shell_dimensions,
        profiles,
        shells,
        shell_moment_blocks,
        shell_kinetic_moment_blocks,
        shell_interaction_moment_blocks,
        shell_exact_lmax,
        shell_kinetic_lcap,
        shell_kinetic_lexpand,
        shell_kinetic_tail,
        shell_interaction_lcap,
        shell_interaction_lexpand,
        shell_interaction_tail,
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
        shell_interaction_lcap = copy(assembly.shell_interaction_lcap),
        shell_interaction_lexpand = copy(assembly.shell_interaction_lexpand),
        interaction_pair_plan = 2 * maximum(assembly.shell_interaction_lexpand),
        interaction_exact_lower_bound = 2 * maximum(assembly.shell_exact_lmax),
        max_shell_kinetic_tail = maximum(assembly.shell_kinetic_tail),
        max_shell_interaction_tail = maximum(assembly.shell_interaction_tail),
        max_shell_overlap_error = maximum(diag.overlap_error for diag in shell_diags),
        max_shell_injected_exactness_error = maximum(diag.injected_exactness_error for diag in shell_diags),
        max_shell_injected_kinetic_error = maximum(diag.injected_kinetic_error for diag in shell_diags),
        max_diagonal_overlap_error = maximum(diagonal_overlap_errors),
        max_diagonal_injected_kinetic_error = maximum(diagonal_kinetic_errors),
        max_pair_overlap_symmetry_error = maximum(pair_overlap_symmetry),
        max_pair_kinetic_symmetry_error = maximum(pair_kinetic_symmetry),
    )
end
