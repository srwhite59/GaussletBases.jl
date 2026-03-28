"""
    ShellLocalInjectedAngularBasis

Experimental shell-local injected angular basis built on one vendored sphere
point set.

This is the first narrow scientific import for the angular research track. It
is intentionally shell-local only: no radial coupling, no shell schedule, and
no atom-wide assembly are included here yet.
"""
struct ShellLocalInjectedAngularBasis
    point_set::SpherePointSet
    beta::Float64
    shell_weights::Vector{Float64}
    theta_nn::Vector{Float64}
    kappa::Vector{Float64}
    l_inject::Int
    injected_channels::YlmChannelSet
    prototype_count::Int
    injected_count::Int
    whitened_complement_count::Int
    final_count::Int
    prototype_overlap::Matrix{Float64}
    prototype_kinetic::Matrix{Float64}
    prototype_orthonormalizer::Matrix{Float64}
    ylm_prototype_coupling::Matrix{Float64}
    raw_to_grand::Matrix{Float64}
    grand_coefficients::Matrix{Float64}
    ylm_coefficients::Matrix{Float64}
    prototype_coefficients::Matrix{Float64}
    final_overlap::Matrix{Float64}
    final_kinetic::Matrix{Float64}
    injected_overlap::Matrix{Float64}
    injected_kinetic::Matrix{Float64}
end

struct ShellLocalAngularProfileKey
    point_set_source_tag::String
    point_set_checksum::String
    order::Int
    beta::Float64
    l_inject::Int
    tau::Float64
    whiten::Symbol
    gauge_version::Symbol
end

struct ShellLocalAngularProfile
    key::ShellLocalAngularProfileKey
    basis::ShellLocalInjectedAngularBasis
    labels::Vector{String}
    block_kinds::Vector{Symbol}
    exact_labels::Vector{String}
    mixed_labels::Vector{String}
    gauge_metadata::NamedTuple
    diagnostics::NamedTuple
    profile_id::String
end

struct ShellLocalAngularProfileOverlap
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

const _SHELL_LOCAL_ANGULAR_PROFILE_GAUGE_VERSION =
    :v1_seed_order_dominant_positive
const _SHELL_LOCAL_ANGULAR_PROFILE_CACHE =
    Ref{Union{Nothing,Dict{ShellLocalAngularProfileKey,ShellLocalAngularProfile}}}(nothing)
const _SHELL_LOCAL_ANGULAR_PROFILE_OVERLAP_CACHE =
    Ref{Union{Nothing,Dict{Tuple{String,String},ShellLocalAngularProfileOverlap}}}(nothing)

function Base.show(io::IO, shell::ShellLocalInjectedAngularBasis)
    print(
        io,
        "ShellLocalInjectedAngularBasis(order=",
        shell.point_set.order,
        ", l_inject=",
        shell.l_inject,
        ", injected_count=",
        shell.injected_count,
        ", whitened_complement_count=",
        shell.whitened_complement_count,
        ", final_count=",
        shell.final_count,
        ")",
    )
end

function Base.show(io::IO, profile::ShellLocalAngularProfile)
    print(
        io,
        "ShellLocalAngularProfile(order=",
        profile.key.order,
        ", l_inject=",
        profile.key.l_inject,
        ", dim=",
        profile.basis.final_count,
        ", profile_id=",
        repr(profile.profile_id),
        ")",
    )
end

function Base.show(io::IO, overlap::ShellLocalAngularProfileOverlap)
    print(
        io,
        "ShellLocalAngularProfileOverlap(source=",
        repr(overlap.source_profile_id),
        ", target=",
        repr(overlap.target_profile_id),
        ", size=",
        size(overlap.overlap),
        ", shell_independent=",
        overlap.shell_independent,
        ")",
    )
end

function _angular_digest_float_vector(values::AbstractVector{<:Real})
    text = join((bitstring(Float64(value)) for value in values), "\n")
    return bytes2hex(SHA.sha1(text))
end

function _angular_digest_float_matrix(values::AbstractMatrix{<:Real})
    text =
        join((bitstring(Float64(value)) for value in vec(Matrix{Float64}(values))), "\n")
    return bytes2hex(SHA.sha1(text))
end

function _angular_digest_strings(values::AbstractVector{<:AbstractString})
    return bytes2hex(SHA.sha1(join(values, "\n")))
end

function _shell_local_angular_profile_cache()
    cache = _SHELL_LOCAL_ANGULAR_PROFILE_CACHE[]
    if cache === nothing
        cache = Dict{ShellLocalAngularProfileKey,ShellLocalAngularProfile}()
        _SHELL_LOCAL_ANGULAR_PROFILE_CACHE[] = cache
    end
    return cache
end

function _shell_local_angular_profile_overlap_cache()
    cache = _SHELL_LOCAL_ANGULAR_PROFILE_OVERLAP_CACHE[]
    if cache === nothing
        cache = Dict{Tuple{String,String},ShellLocalAngularProfileOverlap}()
        _SHELL_LOCAL_ANGULAR_PROFILE_OVERLAP_CACHE[] = cache
    end
    return cache
end

_angular_small_argument_cutoff() = 1.0e-8

function _double_factorial_odd(n::Int)
    n >= -1 || throw(ArgumentError("_double_factorial_odd requires n >= -1"))
    value = 1.0
    k = n
    while k > 1
        value *= k
        k -= 2
    end
    return value
end

function _associated_legendre(l::Int, m::Int, x::Float64)
    l >= 0 || throw(ArgumentError("_associated_legendre requires l >= 0"))
    0 <= m <= l || throw(ArgumentError("_associated_legendre requires 0 <= m <= l"))
    x_clamped = clamp(x, -1.0, 1.0)

    pmm = 1.0
    if m > 0
        factor = sqrt(max(1.0 - x_clamped^2, 0.0))
        pmm = ((-1)^m) * _double_factorial_odd(2 * m - 1) * factor^m
    end
    l == m && return pmm

    pmmp1 = x_clamped * (2 * m + 1) * pmm
    l == m + 1 && return pmmp1

    p_prev = pmm
    p_curr = pmmp1
    for ell in (m + 2):l
        p_next = ((2 * ell - 1) * x_clamped * p_curr - (ell + m - 1) * p_prev) / (ell - m)
        p_prev = p_curr
        p_curr = p_next
    end
    return p_curr
end

function _real_spherical_harmonic(channel::YlmChannel, direction::AbstractVector{<:Real})
    x = Float64(direction[1])
    y = Float64(direction[2])
    z = Float64(direction[3])
    theta = acos(clamp(z, -1.0, 1.0))
    phi = atan(y, x)
    l = channel.l
    m = channel.m
    mm = abs(m)
    log_ratio =
        SpecialFunctions.loggamma(l - mm + 1) - SpecialFunctions.loggamma(l + mm + 1)
    prefactor = sqrt(((2 * l + 1) / (4 * pi)) * exp(log_ratio))
    plm = _associated_legendre(l, mm, cos(theta))

    if m == 0
        return prefactor * plm
    elseif m > 0
        return sqrt(2.0) * prefactor * plm * cos(mm * phi)
    else
        return sqrt(2.0) * prefactor * plm * sin(mm * phi)
    end
end

function _scaled_spherical_besseli(l::Int, k::Float64)
    k == 0.0 && return l == 0 ? 1.0 : 0.0
    return sqrt(pi / (2 * k)) * SpecialFunctions.besselix(l + 0.5, k)
end

function _prototype_overlap(coordinates::AbstractMatrix{<:Real}, kappa::AbstractVector{<:Real}, i::Int, j::Int)
    ki = Float64(kappa[i])
    kj = Float64(kappa[j])
    c = clamp(dot(view(coordinates, i, :), view(coordinates, j, :)), -1.0, 1.0)
    radius_squared = (ki - kj)^2 + 2.0 * ki * kj * (1.0 + c)
    radius = sqrt(max(radius_squared, 0.0))

    if radius < _angular_small_argument_cutoff()
        sinhc_series = 1.0 + radius_squared / 6.0 + (radius_squared^2) / 120.0
        return 4 * pi * exp(-(ki + kj)) * sinhc_series
    end

    tail = radius - (ki + kj)
    one_minus_tail = -expm1(-2 * radius)
    return 4 * pi * 0.5 * exp(tail) * (one_minus_tail / radius)
end

function _prototype_kinetic(coordinates::AbstractMatrix{<:Real}, kappa::AbstractVector{<:Real}, i::Int, j::Int)
    ki = Float64(kappa[i])
    kj = Float64(kappa[j])
    c = clamp(dot(view(coordinates, i, :), view(coordinates, j, :)), -1.0, 1.0)
    radius_squared = (ki - kj)^2 + 2.0 * ki * kj * (1.0 + c)
    radius = sqrt(max(radius_squared, 0.0))

    if radius < 1.0e-6
        ecut = exp(-(ki + kj))
        z0 = 4 * pi * ecut * (1 + radius_squared / 6 + (radius_squared^2) / 120)
        z1 = 4 * pi * ecut * (1 / 3 + radius_squared / 30 + (radius_squared^2) / 840)
        psi = 4 * pi * ecut / 15
        ua = ki + kj * c
        va = kj + ki * c
        return 0.5 * (ki * kj * (c * (z0 - z1) - (ua * va) * psi))
    end

    tail = radius - (ki + kj)
    scale = 0.5 * exp(tail)
    one_minus_tail = -expm1(-2 * radius)
    exp_minus_two_radius = 1 - one_minus_tail
    term1 = scale * (1 + exp_minus_two_radius)
    term2 = scale * (one_minus_tail / radius)
    z0 = 4 * pi * term2
    z1 = 4 * pi * (term1 - term2) / (radius^2)
    z2 = 4 * pi * (term2 - 2 * (term1 - term2) / (radius^2))
    psi = (z2 - z1) / (radius^2)
    ua = ki + kj * c
    va = kj + ki * c
    return 0.5 * (ki * kj * (c * (z0 - z1) - (ua * va) * psi))
end

function _kappa_from_nearest_neighbor(coordinates::AbstractMatrix{<:Real}; beta::Real = 2.0, eps::Real = 1.0e-12)
    npoints = size(coordinates, 1)
    dots = Matrix{Float64}(coordinates * transpose(coordinates))
    for i in 1:npoints
        dots[i, i] = -Inf
    end
    cmax = [maximum(view(dots, i, :)) for i in 1:npoints]
    theta_nn = acos.(clamp.(cmax, -1.0, 1.0))
    kappa = (Float64(beta)^2) ./ (theta_nn .^ 2 .+ Float64(eps))
    return kappa, theta_nn
end

function _prototype_overlap_matrix(coordinates::AbstractMatrix{<:Real}, kappa::AbstractVector{<:Real})
    npoints = length(kappa)
    overlap = Matrix{Float64}(undef, npoints, npoints)
    for j in 1:npoints, i in 1:j
        value = _prototype_overlap(coordinates, kappa, i, j)
        overlap[i, j] = value
        overlap[j, i] = value
    end
    return overlap
end

function _prototype_kinetic_matrix(coordinates::AbstractMatrix{<:Real}, kappa::AbstractVector{<:Real})
    npoints = length(kappa)
    kinetic = Matrix{Float64}(undef, npoints, npoints)
    for j in 1:npoints, i in 1:j
        value = _prototype_kinetic(coordinates, kappa, i, j)
        kinetic[i, j] = value
        kinetic[j, i] = value
    end
    return kinetic
end

function _symmetric_inverse_square_root(overlap::AbstractMatrix{<:Real}; tau::Real = 1.0e-12)
    eigen_data = eigen(Symmetric(Matrix{Float64}(overlap)))
    max_eval = max(1.0, maximum(abs.(eigen_data.values)))
    cutoff = Float64(tau) * max_eval
    repaired = similar(eigen_data.values)
    for i in eachindex(eigen_data.values)
        value = eigen_data.values[i]
        value < -cutoff &&
            throw(ArgumentError("prototype overlap has a negative eigenvalue $value beyond the stability cutoff $cutoff"))
        repaired[i] = max(value, cutoff)
    end
    return eigen_data.vectors * Diagonal(1.0 ./ sqrt.(repaired)) * transpose(eigen_data.vectors)
end

function _choose_l_inject(npoints::Int)
    ny = 0
    l_keep = 0
    for l in 0:200
        ny_new = ny + (2 * l + 1)
        if 2 * ny_new <= npoints
            ny = ny_new
            l_keep = l
        else
            break
        end
    end
    return l_keep
end

function _ylm_prototype_coupling(
    coordinates::AbstractMatrix{<:Real},
    kappa::AbstractVector{<:Real},
    orthonormalizer::AbstractMatrix{<:Real},
    channels::YlmChannelSet,
)
    npoints = length(kappa)
    nchannels = length(channels)
    ylm_to_prototype = zeros(Float64, nchannels, npoints)

    for row in eachindex(channels.channel_data), j in 1:npoints
        channel = channels[row]
        ylm_to_prototype[row, j] =
            4 * pi * _scaled_spherical_besseli(channel.l, Float64(kappa[j])) *
            _real_spherical_harmonic(channel, view(coordinates, j, :))
    end

    return ylm_to_prototype, ylm_to_prototype * orthonormalizer
end

function _grand_from_injected_coupling(coupling::AbstractMatrix{<:Real}; tau::Real = 1.0e-12, whiten::Symbol = :svd)
    nchannels, npoints = size(coupling)

    if whiten == :chol
        complement_overlap = Symmetric(Matrix(I, npoints, npoints) .- transpose(coupling) * coupling)
        eigen_data = eigen(complement_overlap)
        cutoff = Float64(tau) * max(1.0, maximum(eigen_data.values))
        keep = findall(i -> eigen_data.values[i] > cutoff, eachindex(eigen_data.values))
        complement_transform = eigen_data.vectors[:, keep] * Diagonal(1.0 ./ sqrt.(eigen_data.values[keep]))
    elseif whiten == :svd
        factor = svd(coupling; full = true)
        delta = ones(Float64, npoints)
        for i in eachindex(factor.S)
            delta[i] = max(1 - factor.S[i]^2, 0.0)
        end
        cutoff = Float64(tau) * max(1.0, maximum(delta))
        keep = findall(i -> delta[i] > cutoff, 1:npoints)
        complement_transform = factor.V[:, keep] * Diagonal(1.0 ./ sqrt.(delta[keep]))
    else
        throw(ArgumentError("unsupported whitening mode $(repr(whiten)); expected :svd or :chol"))
    end

    complement_count = size(complement_transform, 2)
    grand_transform = [
        Matrix(I, nchannels, nchannels)   -coupling * complement_transform;
        zeros(Float64, npoints, nchannels) complement_transform;
    ]
    return grand_transform, complement_transform, complement_count
end

function _orthogonal_complement(keep::AbstractMatrix{<:Real}, ambient::AbstractMatrix{<:Real})
    size(keep, 1) == size(ambient, 1) ||
        throw(ArgumentError("orthogonal complement requires matching row counts"))

    nrows = size(ambient, 1)
    nkeep = size(keep, 2)
    if nkeep == 0
        factor = qr(Matrix{Float64}(ambient))
        return Matrix(factor.Q)[:, 1:size(ambient, 2)]
    end

    keep_augmented = Matrix(I, nrows, nrows)
    keep_augmented[:, 1:nkeep] = keep
    keep_span = Matrix(qr(keep_augmented).Q)[:, 1:nkeep]
    ambient_residual = ambient - keep_span * (transpose(keep_span) * ambient)
    nret = size(ambient, 2) - nkeep
    factor = svd(ambient_residual)
    return factor.U[:, 1:nret]
end

function _put_in_orbitals(seed::AbstractMatrix{<:Real}, exact_orbitals::AbstractMatrix{<:Real})
    other_orbitals = _orthogonal_complement(exact_orbitals, seed)
    scaffold = hcat(exact_orbitals, other_orbitals)
    projected = scaffold * (transpose(scaffold) * seed)
    norm_square = Symmetric(transpose(projected) * projected)
    return projected * inv(sqrt(norm_square))
end

function _dominant_column_index(values::AbstractVector{<:Real})
    best_index = 1
    best_value = abs(Float64(values[1]))
    for i in 2:length(values)
        trial = abs(Float64(values[i]))
        if trial > best_value + 1.0e-15
            best_index = i
            best_value = trial
        elseif abs(trial - best_value) <= 1.0e-15 && i < best_index
            best_index = i
        end
    end
    return best_index
end

function _canonicalize_grand_orbitals(
    grand_coefficients::AbstractMatrix{<:Real},
    injected_count::Int,
)
    grand = Matrix{Float64}(grand_coefficients)
    nambient, norbitals = size(grand)
    injected_count <= norbitals ||
        throw(ArgumentError("canonical grand orbital split requires injected_count <= norbitals"))

    exact = copy(grand[:, 1:injected_count])
    exact_signs = ones(Float64, injected_count)
    for i in 1:injected_count
        exact[i, i] < 0.0 && (exact[:, i] .*= -1.0; exact_signs[i] = -1.0)
    end

    mixed_count = norbitals - injected_count
    mixed = zeros(Float64, nambient, mixed_count)
    dominant_indices = Int[]
    dominant_values = Float64[]
    mixed_signs = Float64[]

    if mixed_count > 0
        # Preserve the seed-order mixed complement produced by `_put_in_orbitals`
        # so downstream density-density benchmarks keep the historical working
        # basis. The explicit gauge layer here only fixes per-vector signs.
        mixed .= grand[:, injected_count+1:end]
        dominant_indices = Vector{Int}(undef, mixed_count)
        dominant_values = Vector{Float64}(undef, mixed_count)
        mixed_signs = ones(Float64, mixed_count)
        for j in 1:mixed_count
            dominant_index = _dominant_column_index(view(mixed, :, j))
            dominant_value = mixed[dominant_index, j]
            if dominant_value < 0.0
                mixed[:, j] .*= -1.0
                dominant_value = -dominant_value
                mixed_signs[j] = -1.0
            end
            dominant_indices[j] = dominant_index
            dominant_values[j] = dominant_value
        end
    end

    canonical = hcat(exact, mixed)
    metadata = (
        gauge_version = _SHELL_LOCAL_ANGULAR_PROFILE_GAUGE_VERSION,
        mixed_orientation_strategy = :seed_order,
        exact_signs = exact_signs,
        mixed_dominant_grand_indices = dominant_indices,
        mixed_dominant_grand_values = dominant_values,
        mixed_signs = mixed_signs,
    )
    return canonical, metadata
end

function _kinetic_eigenvalue(channel::YlmChannel)
    return 0.5 * channel.l * (channel.l + 1)
end

function _shell_local_raw_overlap(
    ylm_to_prototype::AbstractMatrix{<:Real},
    prototype_overlap::AbstractMatrix{<:Real},
)
    ny = size(ylm_to_prototype, 1)
    npoints = size(prototype_overlap, 1)
    overlap = zeros(Float64, ny + npoints, ny + npoints)
    overlap[1:ny, 1:ny] .= Matrix(I, ny, ny)
    overlap[1:ny, ny+1:end] .= ylm_to_prototype
    overlap[ny+1:end, 1:ny] .= transpose(ylm_to_prototype)
    overlap[ny+1:end, ny+1:end] .= prototype_overlap
    return overlap
end

function _shell_local_raw_kinetic(
    channels::YlmChannelSet,
    ylm_to_prototype::AbstractMatrix{<:Real},
    prototype_kinetic::AbstractMatrix{<:Real},
)
    ny = length(channels)
    npoints = size(prototype_kinetic, 1)
    kinetic = zeros(Float64, ny + npoints, ny + npoints)
    for i in 1:ny
        eigenvalue = _kinetic_eigenvalue(channels[i])
        kinetic[i, i] = eigenvalue
        kinetic[i, ny+1:end] .= eigenvalue .* view(ylm_to_prototype, i, :)
    end
    kinetic[ny+1:end, 1:ny] .= transpose(kinetic[1:ny, ny+1:end])
    kinetic[ny+1:end, ny+1:end] .= prototype_kinetic
    return kinetic
end

function _shell_local_point_set_checksum(point_set::SpherePointSet)
    return _angular_digest_float_matrix(point_set.coordinates)
end

function _shell_local_profile_key(
    point_set::SpherePointSet;
    beta::Real,
    l_inject::Union{Int,Symbol},
    tau::Real,
    whiten::Symbol,
)
    resolved_l_inject =
        l_inject === :auto ? _choose_l_inject(point_set.cardinality) : Int(l_inject)
    resolved_l_inject >= 0 || throw(ArgumentError("l_inject must be >= 0"))
    return ShellLocalAngularProfileKey(
        point_set.provenance.source_tag,
        _shell_local_point_set_checksum(point_set),
        point_set.order,
        Float64(beta),
        resolved_l_inject,
        Float64(tau),
        whiten,
        _SHELL_LOCAL_ANGULAR_PROFILE_GAUGE_VERSION,
    )
end

function _exact_profile_labels(channels::YlmChannelSet)
    return ["exact_l$(channel.l)_m$(channel.m)" for channel in channels.channel_data]
end

function _mixed_profile_labels(mixed_count::Int)
    return ["mixed_$(i)" for i in 1:mixed_count]
end

function _build_shell_local_angular_profile_uncached(
    point_set::SpherePointSet;
    beta::Real = 2.0,
    l_inject::Union{Int,Symbol} = :auto,
    tau::Real = 1.0e-12,
    whiten::Symbol = :svd,
)
    key = _shell_local_profile_key(
        point_set;
        beta = beta,
        l_inject = l_inject,
        tau = tau,
        whiten = whiten,
    )
    prototype_count = point_set.cardinality
    shell_weights = fill(4 * pi / prototype_count, prototype_count)
    coordinates = point_set.coordinates
    kappa, theta_nn = _kappa_from_nearest_neighbor(coordinates; beta = key.beta)
    prototype_overlap = _prototype_overlap_matrix(coordinates, kappa)
    prototype_kinetic = _prototype_kinetic_matrix(coordinates, kappa)
    prototype_orthonormalizer =
        _symmetric_inverse_square_root(prototype_overlap; tau = key.tau)

    injected_channels = ylm_channels(key.l_inject)
    injected_count = length(injected_channels)
    injected_count < prototype_count ||
        throw(
            ArgumentError(
                "shell-local injection requires injected_count < prototype_count; got $(injected_count) and $(prototype_count)",
            ),
        )

    ylm_to_prototype, ylm_to_orthogonalized =
        _ylm_prototype_coupling(
            coordinates,
            kappa,
            prototype_orthonormalizer,
            injected_channels,
        )
    raw_to_grand, complement_transform, complement_count =
        _grand_from_injected_coupling(ylm_to_orthogonalized; tau = key.tau, whiten = key.whiten)

    exact_orbitals =
        vcat(
            Matrix(I, injected_count, injected_count),
            zeros(Float64, complement_count, injected_count),
        )
    seed_complement =
        Symmetric(
            Matrix(I, prototype_count, prototype_count) .-
            transpose(ylm_to_orthogonalized) * ylm_to_orthogonalized,
        )
    seed =
        vcat(
            ylm_to_orthogonalized,
            transpose(complement_transform) * seed_complement,
        )
    grand_coefficients_raw = _put_in_orbitals(seed, exact_orbitals)
    grand_coefficients, gauge_metadata =
        _canonicalize_grand_orbitals(grand_coefficients_raw, injected_count)

    raw_coefficients_ylm_phi = raw_to_grand * grand_coefficients
    ylm_coefficients = raw_coefficients_ylm_phi[1:injected_count, :]
    phi_coefficients = raw_coefficients_ylm_phi[injected_count+1:end, :]
    prototype_coefficients = prototype_orthonormalizer * phi_coefficients

    raw_overlap = _shell_local_raw_overlap(ylm_to_prototype, prototype_overlap)
    raw_kinetic =
        _shell_local_raw_kinetic(
            injected_channels,
            ylm_to_prototype,
            prototype_kinetic,
        )
    raw_coefficients = vcat(ylm_coefficients, prototype_coefficients)

    final_overlap = transpose(raw_coefficients) * raw_overlap * raw_coefficients
    final_overlap = 0.5 .* (final_overlap .+ transpose(final_overlap))
    final_kinetic = transpose(raw_coefficients) * raw_kinetic * raw_coefficients
    final_kinetic = 0.5 .* (final_kinetic .+ transpose(final_kinetic))
    injected_overlap = raw_overlap[1:injected_count, :] * raw_coefficients
    injected_kinetic = injected_overlap * final_kinetic * transpose(injected_overlap)
    injected_kinetic = 0.5 .* (injected_kinetic .+ transpose(injected_kinetic))

    basis = ShellLocalInjectedAngularBasis(
        point_set,
        key.beta,
        shell_weights,
        theta_nn,
        kappa,
        key.l_inject,
        injected_channels,
        prototype_count,
        injected_count,
        complement_count,
        prototype_count,
        prototype_overlap,
        prototype_kinetic,
        prototype_orthonormalizer,
        ylm_to_prototype,
        raw_to_grand,
        grand_coefficients,
        ylm_coefficients,
        prototype_coefficients,
        final_overlap,
        final_kinetic,
        injected_overlap,
        injected_kinetic,
    )

    exact_labels = _exact_profile_labels(injected_channels)
    mixed_labels = _mixed_profile_labels(prototype_count - injected_count)
    labels = vcat(exact_labels, mixed_labels)
    block_kinds = vcat(
        fill(:exact, injected_count),
        fill(:mixed, prototype_count - injected_count),
    )

    diagnostics = (
        shell = shell_local_injected_angular_diagnostics(basis),
        final_overlap_checksum = _angular_digest_float_matrix(final_overlap),
        final_kinetic_checksum = _angular_digest_float_matrix(final_kinetic),
        grand_coefficients_checksum = _angular_digest_float_matrix(grand_coefficients),
        label_checksum = _angular_digest_strings(labels),
        exact_count = injected_count,
        mixed_count = prototype_count - injected_count,
    )
    profile_id = bytes2hex(
        SHA.sha1(
            join(
                (
                    key.point_set_source_tag,
                    key.point_set_checksum,
                    string(key.order),
                    bitstring(key.beta),
                    string(key.l_inject),
                    bitstring(key.tau),
                    string(key.whiten),
                    string(key.gauge_version),
                    diagnostics.final_overlap_checksum,
                    diagnostics.grand_coefficients_checksum,
                ),
                "\n",
            ),
        ),
    )

    return ShellLocalAngularProfile(
        key,
        basis,
        labels,
        block_kinds,
        exact_labels,
        mixed_labels,
        gauge_metadata,
        diagnostics,
        profile_id,
    )
end

function shell_local_angular_profile(
    order::Int;
    beta::Real = 2.0,
    l_inject::Union{Int,Symbol} = :auto,
    tau::Real = 1.0e-12,
    whiten::Symbol = :svd,
)
    return shell_local_angular_profile(
        sphere_point_set(order);
        beta = beta,
        l_inject = l_inject,
        tau = tau,
        whiten = whiten,
    )
end

function shell_local_angular_profile(
    point_set::SpherePointSet;
    beta::Real = 2.0,
    l_inject::Union{Int,Symbol} = :auto,
    tau::Real = 1.0e-12,
    whiten::Symbol = :svd,
)
    key = _shell_local_profile_key(
        point_set;
        beta = beta,
        l_inject = l_inject,
        tau = tau,
        whiten = whiten,
    )
    cache = _shell_local_angular_profile_cache()
    return get!(cache, key) do
        _build_shell_local_angular_profile_uncached(
            point_set;
            beta = beta,
            l_inject = key.l_inject,
            tau = tau,
            whiten = whiten,
        )
    end
end

function adjacent_shell_local_angular_profile_overlap(
    source::ShellLocalAngularProfile,
    target::ShellLocalAngularProfile,
)
    cache = _shell_local_angular_profile_overlap_cache()
    cache_key = (source.profile_id, target.profile_id)
    return get!(cache, cache_key) do
        overlap = _shell_local_overlap(source.basis, target.basis)
        singular_values = svdvals(overlap)
        diagnostics = (
            overlap_checksum = _angular_digest_float_matrix(overlap),
            inf_norm = opnorm(overlap, Inf),
            frobenius_norm = norm(overlap),
            min_singular_value = minimum(singular_values),
            max_singular_value = maximum(singular_values),
            exact_block_inf_norm =
                opnorm(overlap[1:source.basis.injected_count, 1:target.basis.injected_count], Inf),
        )
        ShellLocalAngularProfileOverlap(
            source.profile_id,
            target.profile_id,
            copy(source.labels),
            copy(target.labels),
            overlap,
            source.basis.injected_count,
            source.basis.final_count - source.basis.injected_count,
            target.basis.injected_count,
            target.basis.final_count - target.basis.injected_count,
            true,
            diagnostics,
        )
    end
end

function adjacent_shell_local_angular_profile_overlap(
    source_order::Int,
    target_order::Int;
    beta::Real = 2.0,
    l_inject::Union{Int,Symbol} = :auto,
    tau::Real = 1.0e-12,
    whiten::Symbol = :svd,
)
    source = shell_local_angular_profile(
        source_order;
        beta = beta,
        l_inject = l_inject,
        tau = tau,
        whiten = whiten,
    )
    target = shell_local_angular_profile(
        target_order;
        beta = beta,
        l_inject = l_inject,
        tau = tau,
        whiten = whiten,
    )
    return adjacent_shell_local_angular_profile_overlap(source, target)
end

"""
    build_shell_local_injected_angular_basis(order::Int; beta=2.0, l_inject=:auto, tau=1e-12, whiten=:svd)
    build_shell_local_injected_angular_basis(point_set::SpherePointSet; beta=2.0, l_inject=:auto, tau=1e-12, whiten=:svd)

Build the first experimental shell-local injected angular basis on one vendored
sphere point set.

This constructor intentionally stops at shell-local overlap and angular kinetic
diagnostics. It does not couple shells radially or assemble atom-wide angular
workflows yet.
"""
function build_shell_local_injected_angular_basis(
    order::Int;
    beta::Real = 2.0,
    l_inject::Union{Int,Symbol} = :auto,
    tau::Real = 1.0e-12,
    whiten::Symbol = :svd,
)
    return build_shell_local_injected_angular_basis(
        sphere_point_set(order);
        beta = beta,
        l_inject = l_inject,
        tau = tau,
        whiten = whiten,
    )
end

function build_shell_local_injected_angular_basis(
    point_set::SpherePointSet;
    beta::Real = 2.0,
    l_inject::Union{Int,Symbol} = :auto,
    tau::Real = 1.0e-12,
    whiten::Symbol = :svd,
)
    return shell_local_angular_profile(
        point_set;
        beta = beta,
        l_inject = l_inject,
        tau = tau,
        whiten = whiten,
    ).basis
end

"""
    shell_local_injected_angular_diagnostics(shell::ShellLocalInjectedAngularBasis)

Return a compact diagnostic bundle for one experimental shell-local injected
angular basis.
"""
function shell_local_injected_angular_diagnostics(shell::ShellLocalInjectedAngularBasis)
    identity_final = Matrix(I, shell.final_count, shell.final_count)
    identity_injected = Matrix(I, shell.injected_count, shell.injected_count)
    expected_kinetic = Diagonal([_kinetic_eigenvalue(channel) for channel in shell.injected_channels])
    injected_gram = shell.injected_overlap * transpose(shell.injected_overlap)

    return (
        order = shell.point_set.order,
        cardinality = shell.point_set.cardinality,
        l_inject = shell.l_inject,
        prototype_count = shell.prototype_count,
        injected_count = shell.injected_count,
        whitened_complement_count = shell.whitened_complement_count,
        final_count = shell.final_count,
        overlap_error = opnorm(shell.final_overlap - identity_final, Inf),
        injected_exactness_error = opnorm(injected_gram - identity_injected, Inf),
        injected_kinetic_error = opnorm(shell.injected_kinetic - expected_kinetic, Inf),
        injected_kinetic_eigenvalues = eigvals(Symmetric(shell.injected_kinetic)),
        expected_injected_kinetic_eigenvalues = diag(expected_kinetic),
        nn_ratio = shell.point_set.nn_ratio,
        nearest_neighbor_angles = copy(shell.theta_nn),
        kappa = copy(shell.kappa),
    )
end
