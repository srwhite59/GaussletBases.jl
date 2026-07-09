Base.@kwdef struct AtomicHFReferencePacketSpec
    atom::String
    nuclear_charge::Float64
    electron_count::Int
    basis_name::String
    lmax::Int
    center::NTuple{3,Float64} = (0.0, 0.0, 0.0)
    ns::Int = 5
    core_spacing::Float64
    radius::Float64 = 5.0
    reference_spacing::Float64 = 1.0
    s_factor::Float64 = 1.0
    basisfile::Union{Nothing,String} = nothing
    fill_shell_convention::String
    reference_energy::Union{Nothing,Float64} = nothing
end

Base.@kwdef struct AtomicDensityFitOptions
    base_width_min::Float64 = 0.002
    base_width_max::Float64 = 6.0
    base_ratio::Float64 = 1.14
    drop_low::Int = 0
    drop_high::Int = 7
    radial_grid_min::Float64 = 1.0e-6
    radial_grid_max::Float64 = 30.0
    radial_grid_count::Int = 4000
    svd_rtol::Float64 = 1.0e-12
    charge_weight::Float64 = 5.0e4
end

Base.@kwdef struct AtomicPotentialFitOptions
    fixed_broad_terms::Int = 5
    drop_tight_terms::Int = 12
    refit_after_trim::Bool = false
    radial_grid_core_min::Float64 = 1.0e-8
    radial_grid_core_max::Float64 = 0.2
    radial_grid_tail_max::Float64 = 100.0
    radial_grid_core_count::Int = 900
    radial_grid_mid_count::Int = 900
    radial_grid_tail_count::Int = 500
    svd_rtol::Float64 = 1.0e-12
end

struct AtomicReferenceDensityFit
    betas::Vector{Float64}
    widths::Vector{Float64}
    weights::Vector{Float64}
    radial_grid::Vector{Float64}
    radial_target::Vector{Float64}
    radial_fit::Vector{Float64}
    row::NamedTuple
end

struct AtomicReferencePotentialFit
    coefficients::Vector{Float64}
    exponents::Vector{Float64}
    radial_grid::Vector{Float64}
    radial_exact::Vector{Float64}
    radial_fit::Vector{Float64}
    radial_error::Vector{Float64}
    row::NamedTuple
end

struct AtomicHFReferencePacket
    spec::AtomicHFReferencePacketSpec
    supplement::Any
    overlap::Matrix{Float64}
    overlap_fingerprint::String
    occupied_coefficients::Matrix{Float64}
    occupations::Vector{Float64}
    orbital_energies::Vector{Float64}
    density_matrix::Matrix{Float64}
    rhf_diagnostics::NamedTuple
    density_fit::AtomicReferenceDensityFit
    potential_fit::AtomicReferencePotentialFit
    validation::NamedTuple
    provenance::NamedTuple
end

function atomic_hf_reference_packet_spec(;
    atom::AbstractString,
    nuclear_charge::Real,
    electron_count::Integer,
    basis_name::AbstractString = "cc-pV5Z",
    lmax::Integer = 1,
    center = (0.0, 0.0, 0.0),
    ns::Integer = 5,
    core_spacing::Real,
    radius::Real = 5.0,
    reference_spacing::Real = 1.0,
    s_factor::Real = 1.0,
    basisfile = nothing,
    fill_shell_convention::AbstractString,
    reference_energy = nothing,
)
    electron_count > 0 || throw(ArgumentError("electron_count must be positive"))
    iseven(electron_count) ||
        throw(ArgumentError("initial atomic reference packets are closed-shell RHF only"))
    lmax >= 0 || throw(ArgumentError("lmax must be nonnegative"))
    ns > 0 || throw(ArgumentError("ns must be positive"))
    core_spacing > 0 || throw(ArgumentError("core_spacing must be positive"))
    s_factor > 0 || throw(ArgumentError("s_factor must be positive"))
    center_tuple = ntuple(i -> Float64(center[i]), 3)
    return AtomicHFReferencePacketSpec(
        String(atom),
        Float64(nuclear_charge),
        Int(electron_count),
        String(basis_name),
        Int(lmax),
        center_tuple,
        Int(ns),
        Float64(core_spacing),
        Float64(radius),
        Float64(reference_spacing),
        Float64(s_factor),
        isnothing(basisfile) ? nothing : String(basisfile),
        String(fill_shell_convention),
        isnothing(reference_energy) ? nothing : Float64(reference_energy),
    )
end

be_core_reference_packet_spec(; ns::Integer = 5,
    core_spacing::Real = 1.2 / (4.0 * (Int(ns) - 1)),
    kwargs...) = atomic_hf_reference_packet_spec(;
        atom = "Be", nuclear_charge = 4.0, electron_count = 2,
        basis_name = "cc-pV5Z", lmax = 1, ns, core_spacing,
        fill_shell_convention = "Be Z=4 closed-shell 1s^2 core-screening reference",
        kwargs...)

ne_all_electron_reference_packet_spec(; ns::Integer = 5,
    core_spacing::Real = 1.2 / (10.0 * (Int(ns) - 1)),
    kwargs...) = atomic_hf_reference_packet_spec(;
        atom = "Ne", nuclear_charge = 10.0, electron_count = 10,
        basis_name = "cc-pV5Z", lmax = 1, ns, core_spacing,
        fill_shell_convention = "Ne closed-shell all-electron 1s^2 2s^2 2p^6 reference",
        reference_energy = -128.547098109,
        kwargs...)

_sym(A) = Matrix{Float64}(0.5 .* (Matrix{Float64}(A) .+ transpose(Matrix{Float64}(A))))

function _matrix_fingerprint(A)
    bytes = reinterpret(UInt8, vec(Matrix{Float64}(A)))
    return bytes2hex(sha256(bytes))
end

function _system_spec(spec::AtomicHFReferencePacketSpec)
    neutral_count = round(Int, spec.nuclear_charge)
    isapprox(spec.nuclear_charge, neutral_count; atol = 1.0e-12, rtol = 0.0) ||
        throw(ArgumentError("one-center basis construction requires integer explicit nuclear charge"))
    nup = cld(neutral_count, 2)
    ndn = neutral_count - nup
    return (; atom_symbols = [spec.atom],
        nuclear_charges = [spec.nuclear_charge],
        atom_locations = [spec.center],
        nup,
        ndn)
end

function _basis_spec(spec::AtomicHFReferencePacketSpec)
    return (; q = spec.ns,
        core_spacing = spec.core_spacing,
        radius = spec.radius,
        reference_spacing = spec.reference_spacing,
        d = spec.core_spacing,
        s_factor = spec.s_factor)
end

function _supplement_spec(spec::AtomicHFReferencePacketSpec)
    return (; basis_by_center = [spec.basis_name],
        lmax = spec.lmax,
        uncontracted = false,
        width_filtering = nothing,
        basisfile = spec.basisfile)
end

function build_atomic_reference_supplement(spec::AtomicHFReferencePacketSpec)
    parent = _GB_PARENT
    base = getfield(parent, :cartesian_base_working_basis)(
        _system_spec(spec); basis = _basis_spec(spec), supplemented = true)
    supplement = getfield(parent, :cartesian_residual_gto_supplement_basis)(
        base, _supplement_spec(spec))
    return supplement.basis
end

function _orbital_arrays(supplement)
    n = length(supplement.orbitals)
    labels = String[orb.label for orb in supplement.orbitals]
    powers = Matrix{Int}(undef, 3, n)
    centers = Matrix{Float64}(undef, 3, n)
    exponents = Vector{Vector{Float64}}(undef, n)
    coefficients = Vector{Vector{Float64}}(undef, n)
    primitive_counts = Vector{Int}(undef, n)
    for (i, orb) in pairs(supplement.orbitals)
        powers[:, i] .= Int.(orb.angular_powers)
        centers[:, i] .= Float64.(orb.center)
        exponents[i] = Float64.(orb.exponents)
        coefficients[i] = Float64.(orb.coefficients)
        primitive_counts[i] = length(orb.exponents)
    end
    return (; labels, powers, centers, exponents, coefficients, primitive_counts)
end

function _axis_prefactor(orbital, primitive, axis)
    return getfield(_GB_PARENT.GaussianAnalyticIntegrals,
        :polynomial_gaussian_shell_prefactor)(
        Float64(orbital.exponents[primitive]), orbital.angular_powers[axis])
end

function _axis_basic(left, li, right, ri, axis; extra_exponent = 0.0,
    extra_center = 0.0)
    return getfield(_GB_PARENT.GaussianAnalyticIntegrals,
        :polynomial_gaussian_basic_integral)(
        Float64(left.exponents[li]),
        Float64(left.center[axis]),
        left.angular_powers[axis],
        _axis_prefactor(left, li, axis),
        Float64(right.exponents[ri]),
        Float64(right.center[axis]),
        right.angular_powers[axis],
        _axis_prefactor(right, ri, axis);
        extra_exponent = Float64(extra_exponent),
        extra_center = Float64(extra_center))
end

function _contracted_pair_value(f, left, right)
    value = 0.0
    for ri in eachindex(right.exponents), li in eachindex(left.exponents)
        value += Float64(left.coefficients[li]) *
            Float64(right.coefficients[ri]) * f(left, li, right, ri)
    end
    return value
end

function _axis_kinetic(left, li, right, ri, axis)
    return getfield(_GB_PARENT.GaussianAnalyticIntegrals,
        :polynomial_gaussian_kinetic_integral)(
        Float64(left.exponents[li]),
        Float64(left.center[axis]),
        left.angular_powers[axis],
        _axis_prefactor(left, li, axis),
        Float64(right.exponents[ri]),
        Float64(right.center[axis]),
        right.angular_powers[axis],
        _axis_prefactor(right, ri, axis))
end

function _nuclear_matrix(supplement, charge, center, expansion)
    n = length(supplement.orbitals)
    V = zeros(Float64, n, n)
    for (coefficient, exponent) in zip(expansion.coefficients, expansion.exponents)
        for j in 1:n, i in 1:j
            left = supplement.orbitals[i]
            right = supplement.orbitals[j]
            value = _contracted_pair_value(left, right) do l, li, r, ri
                _axis_basic(l, li, r, ri, 1; extra_exponent = exponent,
                    extra_center = center[1]) *
                _axis_basic(l, li, r, ri, 2; extra_exponent = exponent,
                    extra_center = center[2]) *
                _axis_basic(l, li, r, ri, 3; extra_exponent = exponent,
                    extra_center = center[3])
            end
            V[i, j] -= Float64(charge) * Float64(coefficient) * value
            V[j, i] = V[i, j]
        end
    end
    return _sym(V)
end

function _supplement_overlap_kinetic(supplement)
    n = length(supplement.orbitals)
    S = Matrix{Float64}(getfield(_GB_PARENT, :_cartesian_supplement_cross_overlap)(
        supplement, supplement))
    K = zeros(Float64, n, n)
    for j in 1:n, i in 1:j
        left = supplement.orbitals[i]
        right = supplement.orbitals[j]
        kinetic = _contracted_pair_value(left, right) do l, li, r, ri
            _axis_kinetic(l, li, r, ri, 1) *
            _axis_basic(l, li, r, ri, 2) *
            _axis_basic(l, li, r, ri, 3) +
            _axis_basic(l, li, r, ri, 1) *
            _axis_kinetic(l, li, r, ri, 2) *
            _axis_basic(l, li, r, ri, 3) +
            _axis_basic(l, li, r, ri, 1) *
            _axis_basic(l, li, r, ri, 2) *
            _axis_kinetic(l, li, r, ri, 3)
        end
        K[i, j] = kinetic
        K[j, i] = kinetic
    end
    return S, _sym(K)
end

_pair_index(i, j, n) = getfield(_GB_PARENT, :gaussian_coulomb_pair_index)(i, j, n)

function _coulomb_exchange_matrices(eri, density)
    n = size(density, 1)
    J = zeros(Float64, n, n)
    K = zeros(Float64, n, n)
    for j in 1:n, i in 1:n
        ij = _pair_index(i, j, n)
        value_j = 0.0
        value_k = 0.0
        for l in 1:n, k in 1:n
            value_j += density[k, l] * eri[ij, _pair_index(k, l, n)]
            value_k += density[k, l] *
                eri[_pair_index(i, k, n), _pair_index(j, l, n)]
        end
        J[i, j] = value_j
        K[i, j] = value_k
    end
    return J, K
end

function _orthogonalizer(S; cutoff = 1.0e-12)
    eig = eigen(Symmetric(_sym(S)))
    minimum(eig.values) > cutoff ||
        throw(ArgumentError("atomic reference supplement overlap is rank deficient"))
    X = eig.vectors * Diagonal(1 ./ sqrt.(eig.values)) * transpose(eig.vectors)
    sqrtS = eig.vectors * Diagonal(sqrt.(eig.values)) * transpose(eig.vectors)
    return X, sqrtS, eig.values
end

function _generalized_lowest_orbitals(F, S)
    X, _, _ = _orthogonalizer(S)
    eig = eigen(Symmetric(transpose(X) * F * X))
    C = X * eig.vectors
    return eig.values, C
end

function solve_atomic_supplement_rhf(
    supplement,
    spec::AtomicHFReferencePacketSpec;
    maxiter::Int = 100,
    energy_tol::Float64 = 1.0e-11,
    density_tol::Float64 = 1.0e-10,
)
    S, T = _supplement_overlap_kinetic(supplement)
    expansion = getfield(_GB_PARENT, :coulomb_gaussian_expansion)(doacc = false)
    V = _nuclear_matrix(supplement, spec.nuclear_charge, spec.center, expansion)
    hcore = T + V
    eri = getfield(_GB_PARENT, :gaussian_coulomb_pair_matrix)(
        supplement; expansion, max_orbitals = nothing)
    values, C = _generalized_lowest_orbitals(hcore, S)
    occ = div(spec.electron_count, 2)
    density = 2.0 .* (C[:, 1:occ] * transpose(C[:, 1:occ]))
    energy = Inf
    converged = false
    iterations = 0
    for iteration in 1:maxiter
        J, K = _coulomb_exchange_matrices(eri, density)
        F = hcore + J - 0.5 .* K
        values, C = _generalized_lowest_orbitals(F, S)
        new_density = 2.0 .* (C[:, 1:occ] * transpose(C[:, 1:occ]))
        new_J, new_K = _coulomb_exchange_matrices(eri, new_density)
        new_F = hcore + new_J - 0.5 .* new_K
        new_energy = 0.5 * sum(new_density .* (hcore + new_F))
        density_change = norm(new_density - density, Inf)
        if isfinite(energy) && abs(new_energy - energy) < energy_tol &&
                density_change < density_tol
            density = new_density
            energy = new_energy
            converged = true
            iterations = iteration
            break
        end
        density = new_density
        energy = new_energy
        iterations = iteration
    end
    _, sqrtS, overlap_values = _orthogonalizer(S)
    density = _sym(density)
    density_eigs = eigvals(Symmetric(sqrtS * density * sqrtS))
    return (; method = :closed_shell_supplement_rhf,
        spin_convention = :closed_shell_spin_averaged,
        overlap = S,
        kinetic = T,
        nuclear = V,
        hcore,
        orbital_energies = values,
        occupied_orbitals = Matrix{Float64}(C[:, 1:occ]),
        occupations = fill(2.0, occ),
        occupied_count = occ,
        density_total = density,
        energy,
        converged,
        iterations,
        overlap_eig_min = minimum(overlap_values),
        overlap_eig_max = maximum(overlap_values),
        density_eig_min = minimum(density_eigs),
        density_eig_max = maximum(density_eigs),
        density_trace = tr(density * S),
        density_symmetry_error = norm(density - transpose(density), Inf))
end

function _supplement_self_energy(supplement, density)
    pair = getfield(_GB_PARENT, :gaussian_coulomb_pair_matrix)(
        supplement; max_orbitals = length(supplement.orbitals))
    n = size(density, 1)
    value = 0.0
    for a in 1:n, b in 1:n, c in 1:n, d in 1:n
        value += density[a, b] * density[c, d] *
            pair[_pair_index(a, b, n), _pair_index(c, d, n)]
    end
    return Float64(value)
end

function _cloud_supplement(betas, spec::AtomicHFReferencePacketSpec)
    orbitals = _GB_PARENT.CartesianGaussianShellOrbitalRepresentation3D[
        _GB_PARENT.CartesianGaussianShellOrbitalRepresentation3D(
            "fit_s$(i)", (0, 0, 0), spec.center,
            [0.5 * Float64(beta)], [1.0],
            :axiswise_normalized_cartesian_gaussian)
        for (i, beta) in pairs(betas)]
    metadata = (; source_kind = :atomic_hf_reference_density_fit_cloud,
        atom = spec.atom, basis_name = "$(spec.basis_name)_density_fit",
        lmax = 0, nuclei = NTuple{3,Float64}[spec.center])
    return _GB_PARENT.CartesianGaussianShellSupplementRepresentation3D(
        :atomic_hf_reference_density_fit_cloud, orbitals, metadata)
end

function _cloud_self_energy(betas, weights, spec::AtomicHFReferencePacketSpec)
    cloud = _cloud_supplement(betas, spec)
    pair = getfield(_GB_PARENT, :gaussian_coulomb_pair_matrix)(
        cloud; max_orbitals = length(betas))
    n = length(betas)
    value = 0.0
    for a in 1:n, b in 1:n
        value += weights[a] * weights[b] *
            pair[_pair_index(a, a, n), _pair_index(b, b, n)]
    end
    return Float64(value)
end

function _width_grid(width_min, width_max, ratio)
    widths = Float64[]
    value = Float64(width_min)
    while value < width_max * (1.0 - 1.0e-12)
        push!(widths, value)
        value *= ratio
    end
    if isempty(widths) || abs(last(widths) - width_max) / width_max > 1.0e-10
        push!(widths, Float64(width_max))
    end
    return widths
end

function _fit_widths(options::AtomicDensityFitOptions)
    widths = _width_grid(options.base_width_min, options.base_width_max,
        options.base_ratio)
    first_index = options.drop_low + 1
    last_index = length(widths) - options.drop_high
    1 <= first_index <= last_index <= length(widths) ||
        throw(ArgumentError("density fit drop range removed all widths"))
    return widths[first_index:last_index]
end

function _radial_grid(options::AtomicDensityFitOptions)
    return vcat(0.0, exp.(range(log(options.radial_grid_min),
        log(options.radial_grid_max), length = options.radial_grid_count)))
end

function _trapz_weights(r)
    w = zeros(Float64, length(r))
    length(r) <= 1 && return w
    w[1] = 0.5 * (r[2] - r[1])
    for i in 2:(length(r) - 1)
        w[i] = 0.5 * (r[i + 1] - r[i - 1])
    end
    w[end] = 0.5 * (r[end] - r[end - 1])
    return w
end

function _orbital_value_axis(orbital, x::Float64)
    dx = x - orbital.center[1]
    dy = -orbital.center[2]
    dz = -orbital.center[3]
    powers = orbital.angular_powers
    value = 0.0
    for primitive in eachindex(orbital.exponents)
        alpha = Float64(orbital.exponents[primitive])
        pref = _axis_prefactor(orbital, primitive, 1) *
            _axis_prefactor(orbital, primitive, 2) *
            _axis_prefactor(orbital, primitive, 3)
        value += Float64(orbital.coefficients[primitive]) * pref *
            dx^powers[1] * dy^powers[2] * dz^powers[3] *
            exp(-alpha * (dx^2 + dy^2 + dz^2))
    end
    return value
end

function _reference_radial_density(supplement, C_occ, occupations, r::Float64)
    vals = [_orbital_value_axis(orbital, r) for orbital in supplement.orbitals]
    density = 0.0
    for col in axes(C_occ, 2)
        amplitude = dot(vals, view(C_occ, :, col))
        density += occupations[col] * amplitude^2
    end
    return density
end

_normalized_s_density(beta, r) = (beta / pi)^1.5 * exp(-beta * r^2)

function _cumulative_charge(r, rho)
    q = zeros(Float64, length(r))
    for i in 2:length(r)
        f0 = 4pi * r[i - 1]^2 * rho[i - 1]
        f1 = 4pi * r[i]^2 * rho[i]
        q[i] = q[i - 1] + 0.5 * (f0 + f1) * (r[i] - r[i - 1])
    end
    return q
end

function _radial_potential_from_density(r, rho)
    dr = _trapz_weights(r)
    shell = 4pi .* r.^2 .* rho
    enclosed = cumsum(shell .* dr)
    tail_integrand = 4pi .* r .* rho
    tail = reverse(cumsum(reverse(tail_integrand .* dr)))
    potential = similar(r)
    for i in eachindex(r)
        potential[i] = r[i] == 0.0 ? tail[i] : enclosed[i] / r[i] + tail[i]
    end
    return potential
end

_gaussian_density_potential(beta, r) =
    r == 0.0 ? 2.0 * sqrt(beta / pi) : erf(sqrt(beta) * r) / r

function _solve_weighted_svd(A, b, weights; rtol = 1.0e-12)
    colnorms = [max(norm(weights .* view(A, :, j)), eps(Float64))
        for j in axes(A, 2)]
    Aw = (A .* weights) ./ reshape(colnorms, 1, :)
    bw = b .* weights
    F = svd(Aw; full = false)
    cutoff = max(rtol * maximum(F.S), eps(Float64))
    keep = findall(>=(cutoff), F.S)
    scaled = F.V[:, keep] * ((transpose(F.U[:, keep]) * bw) ./ F.S[keep])
    return scaled ./ colnorms, F.S, length(keep)
end

function fit_atomic_reference_density(
    supplement,
    rhf,
    spec::AtomicHFReferencePacketSpec;
    options::AtomicDensityFitOptions = AtomicDensityFitOptions(),
)
    widths = _fit_widths(options)
    betas = 1.0 ./ (widths .^ 2)
    r = _radial_grid(options)
    target = Float64[_reference_radial_density(
        supplement, rhf.occupied_orbitals, rhf.occupations, radius) for radius in r]
    target_charge = sum(rhf.occupations)
    exact_self = _supplement_self_energy(supplement, rhf.density_total)
    A0 = Matrix{Float64}(undef, length(r), length(betas))
    for (j, beta) in pairs(betas)
        A0[:, j] .= _normalized_s_density.(beta, r)
    end
    dr = _trapz_weights(r)
    radial_weights = sqrt.(4pi .* r.^2 .* dr)
    A = [A0 .* radial_weights; options.charge_weight .* ones(1, length(betas))]
    b = [target .* radial_weights; options.charge_weight * target_charge]
    colnorms = [max(norm(view(A, :, j)), eps(Float64)) for j in axes(A, 2)]
    As = A ./ reshape(colnorms, 1, :)
    F = svd(As; full = false)
    cutoff = max(options.svd_rtol * maximum(F.S), eps(Float64))
    keep = findall(>=(cutoff), F.S)
    scaled = F.V[:, keep] * ((transpose(F.U[:, keep]) * b) ./ F.S[keep])
    weights = scaled ./ colnorms
    fit = A0 * weights
    residual = fit .- target
    enclosed_target = _cumulative_charge(r, target)
    enclosed_fit = _cumulative_charge(r, fit)
    potential_target = _radial_potential_from_density(r, target)
    potential_fit = zeros(Float64, length(r))
    for (beta, weight) in zip(betas, weights)
        potential_fit .+= weight .* _gaussian_density_potential.(beta, r)
    end
    potential_residual = potential_fit .- potential_target
    fit_self = _cloud_self_energy(betas, weights, spec)
    row = (; fit_kind = :signed_svd_log_width_grid,
        base_width_min = options.base_width_min,
        base_width_max = options.base_width_max,
        base_ratio = options.base_ratio,
        drop_low = options.drop_low,
        drop_high = options.drop_high,
        term_count = length(betas),
        retained_rank = length(keep),
        svd_rtol = options.svd_rtol,
        singular_max = maximum(F.S),
        singular_min = minimum(F.S),
        retained_singular_min = isempty(keep) ? NaN : minimum(F.S[keep]),
        condition_full = maximum(F.S) / minimum(F.S),
        condition_retained = isempty(keep) ? NaN : maximum(F.S[keep]) / minimum(F.S[keep]),
        charge = sum(weights),
        charge_error = sum(weights) - target_charge,
        exact_self_energy = exact_self,
        fit_self_energy = fit_self,
        self_energy_relative_error = abs(fit_self - exact_self) /
            max(abs(exact_self), eps(Float64)),
        radial_relative_l2 = norm(residual .* radial_weights) /
            max(norm(target .* radial_weights), eps(Float64)),
        radial_relative_max = maximum(abs.(residual)) /
            max(maximum(abs.(target)), eps(Float64)),
        enclosed_charge_relative_max = maximum(abs.(enclosed_fit .- enclosed_target)) /
            max(abs(enclosed_target[end]), eps(Float64)),
        potential_relative_l2 = norm(potential_residual .* radial_weights) /
            max(norm(potential_target .* radial_weights), eps(Float64)),
        potential_relative_max = maximum(abs.(potential_residual)) /
            max(maximum(abs.(potential_target)), eps(Float64)),
        weight_min = minimum(weights),
        weight_max = maximum(weights),
        weight_abs_sum = sum(abs, weights),
        abs_sum_over_charge = sum(abs, weights) /
            max(abs(sum(weights)), eps(Float64)),
        negative_weight_count = count(<(-1.0e-12), weights))
    return AtomicReferenceDensityFit(Vector{Float64}(betas), Vector{Float64}(widths),
        Vector{Float64}(weights), Vector{Float64}(r), Vector{Float64}(target),
        Vector{Float64}(fit), row)
end

function _potential_fit_grid(options::AtomicPotentialFitOptions)
    core = exp.(range(log(options.radial_grid_core_min),
        log(options.radial_grid_core_max), length = options.radial_grid_core_count))
    mid = exp.(range(log(options.radial_grid_core_max), log(10.0),
        length = options.radial_grid_mid_count))
    tail = exp.(range(log(10.0), log(options.radial_grid_tail_max),
        length = options.radial_grid_tail_count))
    return unique!(sort!(vcat(0.0, core, mid, tail)))
end

function _density_fit_potential_value(fit::AtomicReferenceDensityFit, r::Float64)
    if r == 0.0
        return sum(fit.weights .* (2.0 .* sqrt.(fit.betas ./ pi)))
    end
    return sum(fit.weights .* erf.(sqrt.(fit.betas) .* r)) / r
end

function _expansion_value(expansion, coeffs, r::Float64)
    return sum(coeffs .* exp.(-expansion.exponents .* r^2))
end

function _radial_fit_weights(grid, exact)
    weights = 1.0 ./ max.(abs.(exact), 1.0e-10)
    weights[grid .<= 0.05] .*= 20.0
    weights[grid .<= 0.005] .*= 20.0
    weights[grid .>= 20.0] .*= 8.0
    return weights
end

function _potential_fit_errors(grid, exact, fit_values, charge)
    err = fit_values .- exact
    rel = abs.(err) ./ max.(abs.(exact), 1.0e-12)
    core = grid .<= 0.05
    mid = (grid .> 0.05) .& (grid .<= 5.0)
    tail = grid .>= 30.0
    return (; err,
        absmax = maximum(abs.(err)),
        relmax = maximum(rel),
        relrms = sqrt(sum(abs2, rel) / length(rel)),
        core_relmax = maximum(rel[core]),
        mid_relmax = maximum(rel[mid]),
        tail_relmax = maximum(rel[tail]),
        tail_charge_error = maximum(abs.(grid[tail] .* fit_values[tail] .- charge)))
end

function fit_atomic_reference_potential(
    density_fit::AtomicReferenceDensityFit;
    options::AtomicPotentialFitOptions = AtomicPotentialFitOptions(),
)
    expansion = getfield(_GB_PARENT, :coulomb_gaussian_expansion)(doacc = false)
    term_count = length(expansion.coefficients)
    fixed = options.fixed_broad_terms
    drop = options.drop_tight_terms
    0 < fixed < term_count || throw(ArgumentError("fixed_broad_terms must be in 1:$(term_count - 1)"))
    0 <= drop < term_count - fixed ||
        throw(ArgumentError("drop_tight_terms removes all optimized terms"))
    charge = density_fit.row.charge
    grid = _potential_fit_grid(options)
    exact = Float64[_density_fit_potential_value(density_fit, r) for r in grid]
    base_coeffs = charge .* expansion.coefficients
    fixed_sources = collect(1:fixed)
    opt_sources = collect((fixed + 1):term_count)
    fixed_values = zeros(Float64, length(grid))
    for idx in fixed_sources
        fixed_values .+= base_coeffs[idx] .* exp.(-expansion.exponents[idx] .* grid.^2)
    end
    A = Matrix{Float64}(undef, length(grid), length(opt_sources))
    for (j, idx) in pairs(opt_sources)
        A[:, j] .= exp.(-expansion.exponents[idx] .* grid.^2)
    end
    weights = _radial_fit_weights(grid, exact)
    opt_coeffs, singulars, rank = _solve_weighted_svd(
        A, exact .- fixed_values, weights; rtol = options.svd_rtol)
    all_coeffs = copy(base_coeffs)
    all_coeffs[(fixed + 1):term_count] .= opt_coeffs
    keep_sources = collect(1:(term_count - drop))
    if options.refit_after_trim
        opt_keep = collect((fixed + 1):(term_count - drop))
        Atrim = Matrix{Float64}(undef, length(grid), length(opt_keep))
        for (j, idx) in pairs(opt_keep)
            Atrim[:, j] .= exp.(-expansion.exponents[idx] .* grid.^2)
        end
        trim_coeffs, trim_singulars, trim_rank = _solve_weighted_svd(
            Atrim, exact .- fixed_values, weights; rtol = options.svd_rtol)
        total_coeffs = vcat(all_coeffs[fixed_sources], trim_coeffs)
        total_exps = expansion.exponents[keep_sources]
        fit_values = fixed_values .+ Atrim * trim_coeffs
        singulars = trim_singulars
        rank = trim_rank
    else
        total_coeffs = all_coeffs[keep_sources]
        total_exps = expansion.exponents[keep_sources]
        fit_values = zeros(Float64, length(grid))
        for (coeff, exponent) in zip(total_coeffs, total_exps)
            fit_values .+= coeff .* exp.(-exponent .* grid.^2)
        end
    end
    errors = _potential_fit_errors(grid, exact, fit_values, charge)
    kept_widths = inv.(sqrt.(total_exps))
    row = (; fit_kind = :fixed_broad_tail_refit_tight_then_trim,
        source_coulomb_terms = term_count,
        fixed_broad_terms = fixed,
        optimized_terms_before_trim = length(opt_sources),
        drop_tight_terms = drop,
        refit_after_trim = options.refit_after_trim,
        total_terms = length(total_coeffs),
        retained_rank = rank,
        singular_min = minimum(singulars),
        singular_max = maximum(singulars),
        charge,
        width_min = minimum(kept_widths),
        width_max = maximum(kept_widths),
        dropped_source_indices = drop == 0 ? "" :
            join(((term_count + 1 - drop):term_count), ","),
        coefficient_min = minimum(total_coeffs),
        coefficient_max = maximum(total_coeffs),
        negative_coefficient_count = count(<(0.0), total_coeffs),
        errors.absmax,
        errors.relmax,
        errors.relrms,
        errors.core_relmax,
        errors.mid_relmax,
        errors.tail_relmax,
        errors.tail_charge_error)
    return AtomicReferencePotentialFit(Vector{Float64}(total_coeffs),
        Vector{Float64}(total_exps), Vector{Float64}(grid),
        Vector{Float64}(exact), Vector{Float64}(fit_values),
        Vector{Float64}(errors.err), row)
end

function _packet_validation(spec, supplement, rhf, density_fit, potential_fit)
    S = Matrix{Float64}(rhf.overlap)
    C = Matrix{Float64}(rhf.occupied_orbitals)
    P = C * Diagonal(rhf.occupations) * transpose(C)
    return (;
        occupied_orthogonality_error = norm(transpose(C) * S * C -
            Matrix{Float64}(I, size(C, 2), size(C, 2)), Inf),
        density_trace = tr(P * S),
        density_trace_error = tr(P * S) - spec.electron_count,
        density_from_coefficients_error = norm(P - rhf.density_total, Inf),
        exact_self_energy_nohalf = density_fit.row.exact_self_energy,
        density_fit_charge_error = density_fit.row.charge_error,
        density_fit_self_energy_relative_error =
            density_fit.row.self_energy_relative_error,
        potential_fit_tail_charge_error = potential_fit.row.tail_charge_error,
        potential_fit_radial_relmax = potential_fit.row.relmax,
        potential_fit_matrix_check_status = :not_run,
        potential_fit_matrix_relative_fro = NaN,
        potential_fit_matrix_max_abs = NaN,
        potential_fit_anchor_status = :not_run,
        potential_fit_anchor_error = NaN)
end

function build_atomic_hf_reference_packet(
    spec::AtomicHFReferencePacketSpec;
    density_options::AtomicDensityFitOptions = AtomicDensityFitOptions(),
    potential_options::AtomicPotentialFitOptions = AtomicPotentialFitOptions(),
)
    supplement = build_atomic_reference_supplement(spec)
    rhf = solve_atomic_supplement_rhf(supplement, spec)
    density_fit = fit_atomic_reference_density(
        supplement, rhf, spec; options = density_options)
    potential_fit = fit_atomic_reference_potential(
        density_fit; options = potential_options)
    validation = _packet_validation(spec, supplement, rhf, density_fit, potential_fit)
    provenance = (; git_commit = _git_commit(),
        construction_timestamp = string(time()),
        code_version = :HP_PQS_ATOMREF_PACKET_FN_01,
        density_fit_tolerances = density_options,
        potential_fit_tolerances = potential_options)
    return AtomicHFReferencePacket(spec, supplement, Matrix{Float64}(rhf.overlap),
        _matrix_fingerprint(rhf.overlap), Matrix{Float64}(rhf.occupied_orbitals),
        Vector{Float64}(rhf.occupations), Vector{Float64}(rhf.orbital_energies),
        Matrix{Float64}(rhf.density_total),
        (; method = rhf.method,
            spin_convention = rhf.spin_convention,
            rhf_energy = rhf.energy,
            converged = rhf.converged,
            iterations = rhf.iterations,
            occupied_count = rhf.occupied_count,
            overlap_eig_min = rhf.overlap_eig_min,
            overlap_eig_max = rhf.overlap_eig_max,
            density_eig_min = rhf.density_eig_min,
            density_eig_max = rhf.density_eig_max,
            density_trace = rhf.density_trace,
            density_symmetry_error = rhf.density_symmetry_error),
        density_fit, potential_fit, validation, provenance)
end

function _git_commit()
    try
        return readchomp(`git rev-parse HEAD`)
    catch
        return "unknown"
    end
end

function _write_common_packet_fields(f, packet::AtomicHFReferencePacket)
    spec = packet.spec
    arrays = _orbital_arrays(packet.supplement)
    f["artifact_kind"] = :atomic_hf_reference_density_fit
    f["convention_id"] = :atomic_hf_reference_density_fit_v1
    f["system/atom"] = spec.atom
    f["system/element"] = spec.atom
    f["system/Z"] = spec.nuclear_charge
    f["system/nuclear_charge"] = spec.nuclear_charge
    f["system/charge"] = spec.nuclear_charge - spec.electron_count
    f["system/electron_count"] = spec.electron_count
    f["system/method"] = :RHF
    f["system/closed_shell"] = true
    f["system/spin_convention"] = :closed_shell_spin_averaged
    f["system/fill_shell_convention"] = spec.fill_shell_convention
    f["system/center"] = collect(spec.center)

    for prefix in ("basis", "supplement_basis")
        f["$(prefix)/name"] = spec.basis_name
        f["$(prefix)/basisfile"] = isnothing(spec.basisfile) ? "" : spec.basisfile
        f["$(prefix)/lmax"] = spec.lmax
        f["$(prefix)/supplement_count"] = length(packet.supplement.orbitals)
        f["$(prefix)/orbital_labels"] = arrays.labels
        f["$(prefix)/angular_powers"] = arrays.powers
        f["$(prefix)/centers"] = arrays.centers
        f["$(prefix)/exponents"] = arrays.exponents
        f["$(prefix)/coefficients"] = arrays.coefficients
        f["$(prefix)/primitive_counts"] = arrays.primitive_counts
        f["$(prefix)/overlap"] = packet.overlap
        f["$(prefix)/overlap_fingerprint_sha256"] = packet.overlap_fingerprint
    end
end

function write_atomic_hf_reference_packet(path::AbstractString,
    packet::AtomicHFReferencePacket)
    mkpath(dirname(path))
    jldopen(path, "w") do f
        _write_common_packet_fields(f, packet)
        hf = packet.rhf_diagnostics
        f["hf/occupied_coefficients_AA"] = packet.occupied_coefficients
        f["hf/occupations"] = packet.occupations
        f["hf/orbital_energies"] = packet.orbital_energies
        f["hf/density_matrix_AA"] = packet.density_matrix
        f["hf/rhf_energy"] = hf.rhf_energy
        f["hf/converged"] = hf.converged
        f["hf/iterations"] = hf.iterations
        f["hf/method"] = hf.method
        f["hf/spin_convention"] = hf.spin_convention
        f["hf/occupied_count"] = hf.occupied_count
        f["hf/S_AA_orthogonality_error"] =
            packet.validation.occupied_orthogonality_error
        f["hf/density_trace"] = hf.density_trace
        f["hf/density_symmetry_error"] = hf.density_symmetry_error
        f["hf/overlap_eig_min"] = hf.overlap_eig_min
        f["hf/overlap_eig_max"] = hf.overlap_eig_max
        f["hf/density_eig_min"] = hf.density_eig_min
        f["hf/density_eig_max"] = hf.density_eig_max

        df = packet.density_fit
        f["density_fit/fit_kind"] = df.row.fit_kind
        f["density_fit/betas"] = df.betas
        f["density_fit/widths"] = df.widths
        f["density_fit/weights"] = df.weights
        f["density_fit/radial_grid"] = df.radial_grid
        f["density_fit/radial_target"] = df.radial_target
        f["density_fit/radial_fit"] = df.radial_fit
        for name in propertynames(df.row)
            name == :fit_kind && continue
            f["density_fit/$(String(name))"] = getproperty(df.row, name)
        end

        pf = packet.potential_fit
        f["potential_fit/fit_kind"] = pf.row.fit_kind
        f["potential_fit/total_coefficients"] = pf.coefficients
        f["potential_fit/total_exponents"] = pf.exponents
        f["potential_fit/radial_grid"] = pf.radial_grid
        f["potential_fit/radial_exact"] = pf.radial_exact
        f["potential_fit/radial_fit"] = pf.radial_fit
        f["potential_fit/radial_error"] = pf.radial_error
        for name in propertynames(pf.row)
            name == :fit_kind && continue
            f["potential_fit/$(String(name))"] = getproperty(pf.row, name)
        end

        for name in propertynames(packet.validation)
            f["validation/$(String(name))"] = getproperty(packet.validation, name)
        end
        for name in propertynames(packet.provenance)
            value = getproperty(packet.provenance, name)
            f["provenance/$(String(name))"] = value isa AtomicDensityFitOptions ||
                value isa AtomicPotentialFitOptions ? repr(value) : value
        end
        f["provenance/input_ns"] = packet.spec.ns
        f["provenance/input_core_spacing"] = packet.spec.core_spacing
        f["provenance/input_s_factor"] = packet.spec.s_factor
        f["provenance/reference_energy"] =
            isnothing(packet.spec.reference_energy) ? NaN : packet.spec.reference_energy
    end
    return path
end

function read_atomic_hf_reference_packet(path::AbstractString)
    return jldopen(path, "r") do f
        fit_row = (; charge = Float64(f["density_fit/charge"]),
            charge_error = Float64(f["density_fit/charge_error"]),
            exact_self_energy = Float64(f["density_fit/exact_self_energy"]),
            fit_self_energy = Float64(f["density_fit/fit_self_energy"]),
            self_energy_relative_error =
                Float64(f["density_fit/self_energy_relative_error"]),
            radial_relative_l2 = Float64(f["density_fit/radial_relative_l2"]),
            radial_relative_max = Float64(f["density_fit/radial_relative_max"]),
            enclosed_charge_relative_max =
                Float64(f["density_fit/enclosed_charge_relative_max"]),
            potential_relative_l2 = Float64(f["density_fit/potential_relative_l2"]),
            potential_relative_max = Float64(f["density_fit/potential_relative_max"]),
            term_count = Int(f["density_fit/term_count"]),
            retained_rank = Int(f["density_fit/retained_rank"]),
            negative_weight_count = Int(f["density_fit/negative_weight_count"]))
        pot_row = (; charge = Float64(f["potential_fit/charge"]),
            total_terms = Int(f["potential_fit/total_terms"]),
            fixed_broad_terms = Int(f["potential_fit/fixed_broad_terms"]),
            drop_tight_terms = Int(f["potential_fit/drop_tight_terms"]),
            refit_after_trim = Bool(f["potential_fit/refit_after_trim"]),
            absmax = Float64(f["potential_fit/absmax"]),
            relmax = Float64(f["potential_fit/relmax"]),
            relrms = Float64(f["potential_fit/relrms"]),
            tail_charge_error = Float64(f["potential_fit/tail_charge_error"]))
        return (;
            path = String(path),
            artifact_kind = f["artifact_kind"],
            convention_id = f["convention_id"],
            atom = String(f["system/atom"]),
            nuclear_charge = Float64(f["system/nuclear_charge"]),
            electron_count = Int(f["system/electron_count"]),
            fill_shell_convention = String(f["system/fill_shell_convention"]),
            center = Tuple(Float64.(f["system/center"])),
            basis_name = String(f["supplement_basis/name"]),
            basisfile = String(f["supplement_basis/basisfile"]),
            lmax = Int(f["supplement_basis/lmax"]),
            labels = Vector{String}(f["supplement_basis/orbital_labels"]),
            angular_powers = Matrix{Int}(f["supplement_basis/angular_powers"]),
            overlap = Matrix{Float64}(f["supplement_basis/overlap"]),
            overlap_fingerprint = String(f["supplement_basis/overlap_fingerprint_sha256"]),
            C_occ = Matrix{Float64}(f["hf/occupied_coefficients_AA"]),
            occupations = Vector{Float64}(f["hf/occupations"]),
            density_AA = Matrix{Float64}(f["hf/density_matrix_AA"]),
            rhf_energy = Float64(f["hf/rhf_energy"]),
            rhf_converged = Bool(f["hf/converged"]),
            rhf_iterations = Int(f["hf/iterations"]),
            density_fit = (; betas = Vector{Float64}(f["density_fit/betas"]),
                widths = Vector{Float64}(f["density_fit/widths"]),
                weights = Vector{Float64}(f["density_fit/weights"]),
                row = fit_row),
            potential_fit = (; coefficients =
                    Vector{Float64}(f["potential_fit/total_coefficients"]),
                exponents = Vector{Float64}(f["potential_fit/total_exponents"]),
                row = pot_row))
    end
end

function validate_atomic_hf_reference_packet(packet)
    S = packet.overlap
    C = packet.occupied_coefficients
    P = C * Diagonal(packet.occupations) * transpose(C)
    return (;
        occupied_orthogonality_error = norm(transpose(C) * S * C -
            Matrix{Float64}(I, size(C, 2), size(C, 2)), Inf),
        density_trace = tr(P * S),
        density_trace_error = tr(P * S) - packet.spec.electron_count,
        density_matrix_error = norm(P - packet.density_matrix, Inf),
        density_fit_charge_error = packet.density_fit.row.charge_error,
        density_fit_self_energy_relative_error =
            packet.density_fit.row.self_energy_relative_error,
        potential_fit_tail_charge_error =
            packet.potential_fit.row.tail_charge_error,
        potential_fit_radial_relmax = packet.potential_fit.row.relmax)
end

function validate_atomic_hf_reference_packet(path::AbstractString)
    packet = read_atomic_hf_reference_packet(path)
    S = packet.overlap
    C = packet.C_occ
    P = C * Diagonal(packet.occupations) * transpose(C)
    return (;
        artifact_kind = packet.artifact_kind,
        convention_id = packet.convention_id,
        occupied_orthogonality_error = norm(transpose(C) * S * C -
            Matrix{Float64}(I, size(C, 2), size(C, 2)), Inf),
        density_trace = tr(P * S),
        density_trace_error = tr(P * S) - packet.electron_count,
        density_matrix_error = norm(P - packet.density_AA, Inf),
        density_fit_charge_error = packet.density_fit.row.charge_error,
        density_fit_self_energy_relative_error =
            packet.density_fit.row.self_energy_relative_error,
        potential_fit_tail_charge_error =
            packet.potential_fit.row.tail_charge_error,
        potential_fit_radial_relmax = packet.potential_fit.row.relmax)
end

function atomic_reference_packet_p0_q0(packet)
    C = hasproperty(packet, :C_occ) ? packet.C_occ : packet.occupied_coefficients
    occupations = packet.occupations
    P = _sym(C * Diagonal(occupations) * transpose(C))
    return (; P_AA = P, q_AA = diag(P), trace = tr(P * packet.overlap))
end

function _packet_center(packet)
    center =
        hasproperty(packet, :spec) ? packet.spec.center :
        hasproperty(packet, :center) ? packet.center :
        (0.0, 0.0, 0.0)
    return ntuple(axis -> Float64(center[axis]), 3)
end

function _packet_density_fit_spec(packet, center)
    spec = hasproperty(packet, :spec) ? packet.spec :
        atomic_hf_reference_packet_spec(;
            atom = packet.atom,
            nuclear_charge = packet.nuclear_charge,
            electron_count = packet.electron_count,
            basis_name = packet.basis_name,
            lmax = packet.lmax,
            core_spacing = 1.0,
            center,
            fill_shell_convention = packet.fill_shell_convention)
    if hasproperty(packet, :spec) && packet.spec.center != center
        spec = atomic_hf_reference_packet_spec(;
            atom = spec.atom,
            nuclear_charge = spec.nuclear_charge,
            electron_count = spec.electron_count,
            basis_name = spec.basis_name,
            lmax = spec.lmax,
            center,
            ns = spec.ns,
            core_spacing = spec.core_spacing,
            radius = spec.radius,
            reference_spacing = spec.reference_spacing,
            s_factor = spec.s_factor,
            basisfile = spec.basisfile,
            fill_shell_convention = spec.fill_shell_convention,
            reference_energy = spec.reference_energy)
    end
    return spec
end

function _packet_cloud_from_readback(packet; center = _packet_center(packet))
    spec = _packet_density_fit_spec(packet, center)
    cloud = _cloud_supplement(packet.density_fit.betas, spec)
    density = Diagonal(packet.density_fit.weights)
    return cloud, density
end

function atomic_reference_packet_terminal_hartree_gg(
    base,
    packet;
    source::Symbol,
    center = _packet_center(packet),
)
    placement = ntuple(axis -> Float64(center[axis]), 3)
    if source == :density_fit
        cloud, density = _packet_cloud_from_readback(packet; center = placement)
        result = _GB_PARENT.CartesianGaussianRawBlocks.atomic_reference_hartree_gg_block(
            base.terminal_basis, base.parent.parent_axis_bundle_object, cloud, density)
        return result.GG
    elseif source == :potential_fit
        template = getfield(_GB_PARENT, :coulomb_gaussian_expansion)(doacc = false)
        expansion = _GB_PARENT.CoulombGaussianExpansion(
            packet.potential_fit.coefficients,
            packet.potential_fit.exponents;
            del = template.del, s = template.s, c = template.c,
            maxu = template.maxu)
        basis = base.terminal_basis
        bundles = base.parent.parent_axis_bundle_object
        pgdg = Tuple(getfield(_GB_PARENT, :_nested_axis_pgdg)(
            bundles, axis) for axis in (:x, :y, :z))
        matrix = zeros(Float64, basis.final_dimension, basis.final_dimension)
        factors = ntuple(axis ->
            getfield(_GB_PARENT.CartesianFinalBasisRealization,
                :_r3a_centered_factor_terms)(
                pgdg[axis], expansion, placement[axis]), 3)
        getfield(_GB_PARENT.CartesianFinalBasisRealization,
            :_accumulate_terminal_gaussian_sum!)(
            matrix, basis, expansion.coefficients, factors[1], factors[2],
            factors[3]; scale = 1.0)
        return _sym(matrix)
    else
        throw(ArgumentError("source must be :density_fit or :potential_fit"))
    end
end

function atomic_reference_packet_matrix_check(base, packet)
    exact = atomic_reference_packet_terminal_hartree_gg(
        base, packet; source = :density_fit)
    fast = atomic_reference_packet_terminal_hartree_gg(
        base, packet; source = :potential_fit)
    diff = _sym(fast - exact)
    return (; dimension = size(exact, 1),
        relative_fro = norm(diff) / max(norm(exact), eps(Float64)),
        max_abs = maximum(abs, diff),
        diag_max_abs = maximum(abs.(diag(diff))),
        trace_difference = tr(diff),
        exact_trace = tr(exact),
        fast_trace = tr(fast),
        exact_finite = all(isfinite, exact),
        fast_finite = all(isfinite, fast),
        symmetry_error = norm(diff - transpose(diff), Inf))
end
