"""
    EGOIDensityDensityCorrectionResult

Matrix-level result returned by [`egoi_density_density_correction`](@ref).
`interaction_matrix` is the corrected density-density interaction matrix,
`interaction_delta` is the symmetric correction added to the input matrix, and
`diagnostics` records target-space residuals, correction sizes, and product
matrix conditioning.
"""
struct EGOIDensityDensityCorrectionResult
    interaction_matrix::Matrix{Float64}
    interaction_delta::Matrix{Float64}
    diagnostics::NamedTuple
end

"""
    StationaryFockCorrectionResult

Matrix-level result returned by [`stationary_fock_one_body_correction`](@ref).
`fock_matrix` is `F + one_body_delta`; the same symmetric `one_body_delta` is
intended to be added to the one-body Hamiltonian that generated `F`.
"""
struct StationaryFockCorrectionResult
    fock_matrix::Matrix{Float64}
    one_body_delta::Matrix{Float64}
    diagnostics::NamedTuple
end

"""
    HamiltonianCorrectionResult

Combined matrix-level result returned by
[`egoi_stationary_hamiltonian_correction`](@ref) and the ordinary convenience
adapter. `one_body_hamiltonian` and `interaction_matrix` are the final matrices;
`one_body_delta` and `interaction_delta` are final-minus-initial matrix deltas.
"""
struct HamiltonianCorrectionResult
    one_body_hamiltonian::Matrix{Float64}
    interaction_matrix::Matrix{Float64}
    one_body_delta::Matrix{Float64}
    interaction_delta::Matrix{Float64}
    diagnostics::NamedTuple
end

"""
    OrdinaryProjectedHamiltonianCorrectionTarget

Convenience target object produced by
[`ordinary_cartesian_projected_gaussian_target`](@ref). It carries Gaussian
target coefficients, their projection into an orthonormal ordinary working
basis, the dense Gaussian exact target Coulomb matrix, occupations, selected
Gaussian column indices, and projection diagnostics.
"""
struct OrdinaryProjectedHamiltonianCorrectionTarget
    projected_orbitals::Matrix{Float64}
    gaussian_coefficients::Matrix{Float64}
    exact_target::Matrix{Float64}
    occupations::Vector{Float64}
    indices::Vector{Int}
    diagnostics::NamedTuple
end

function _hc_symmetrize(matrix::AbstractMatrix{<:Real})
    dense = Matrix{Float64}(matrix)
    return Matrix{Float64}(0.5 .* (dense .+ transpose(dense)))
end

function _hc_max_abs(values)
    length(values) == 0 && return 0.0
    return Float64(maximum(abs, values))
end

function _hc_rms_abs(values)
    length(values) == 0 && return 0.0
    return Float64(sqrt(sum(abs2, values) / length(values)))
end

function _hc_fro_norm(values)
    return Float64(norm(values))
end

function _hc_relative(numerator::Real, denominator_values)
    denominator = _hc_fro_norm(denominator_values)
    denominator == 0.0 && return (Float64(numerator) == 0.0 ? 0.0 : Inf)
    return Float64(numerator) / denominator
end

function _hc_matrix_delta_diagnostics(delta::AbstractMatrix{<:Real}, baseline::AbstractMatrix{<:Real})
    delta_matrix = Matrix{Float64}(delta)
    baseline_matrix = Matrix{Float64}(baseline)
    delta_max = _hc_max_abs(delta_matrix)
    delta_fro = _hc_fro_norm(delta_matrix)
    return (
        max = delta_max,
        rms = _hc_rms_abs(delta_matrix),
        fro = delta_fro,
        relative_max = _hc_max_abs(baseline_matrix) == 0.0 ?
                       (delta_max == 0.0 ? 0.0 : Inf) :
                       delta_max / _hc_max_abs(baseline_matrix),
        relative_fro = _hc_relative(delta_fro, baseline_matrix),
    )
end

function _hc_residual_diagnostics(residual::AbstractMatrix{<:Real})
    residual_matrix = Matrix{Float64}(residual)
    return (
        fro = _hc_fro_norm(residual_matrix),
        max = _hc_max_abs(residual_matrix),
    )
end

function _hc_qr_orthonormalize_columns(matrix::AbstractMatrix{<:Real})
    dense = Matrix{Float64}(matrix)
    column_count = size(dense, 2)
    column_count > 0 || throw(ArgumentError("target orbital matrix must have at least one column"))
    size(dense, 1) >= column_count || throw(
        ArgumentError("target orbital matrix must have at least as many rows as columns"),
    )
    singular_values = Float64[Float64(value) for value in svdvals(dense)]
    rank_cutoff = max(size(dense)...) * eps(Float64) * maximum(singular_values)
    rank = count(value -> value > rank_cutoff, singular_values)
    rank == column_count || throw(
        ArgumentError("target orbital columns must be linearly independent for QR orthonormalization"),
    )
    qr_factor = qr(dense)
    q = Matrix(qr_factor.Q)[:, 1:column_count]
    return q, singular_values, rank
end

"""
    egoi_target_product_matrix(Qtarget)

Return the density-product matrix used by the EGOI correction. `Qtarget`
columns are target orbitals represented in an orthonormal working basis. Output
columns are pointwise products `Qtarget[:, i] .* Qtarget[:, j]`, with target
pair index `(j - 1) * ntarget + i`.
"""
function egoi_target_product_matrix(Qtarget::AbstractMatrix{<:Real})
    q = Matrix{Float64}(Qtarget)
    basis_count, target_count = size(q)
    product = zeros(Float64, basis_count, target_count^2)
    for j in 1:target_count, i in 1:target_count
        product[:, (j - 1) * target_count + i] .= q[:, i] .* q[:, j]
    end
    return product
end

function _egoi_svd_delta_matrix(
    product::AbstractMatrix{<:Real},
    delta_target::AbstractMatrix{<:Real};
    regularization::Real,
)
    regularization_value = Float64(regularization)
    regularization_value >= 0.0 ||
        throw(ArgumentError("regularization must be nonnegative"))
    decomposition = svd(Matrix{Float64}(product))
    singular_values = Float64[Float64(value) for value in decomposition.S]
    right_vectors = decomposition.V
    left_vectors = decomposition.U
    transformed_target = transpose(right_vectors) * Matrix{Float64}(delta_target) * right_vectors
    transformed_delta = zeros(Float64, length(singular_values), length(singular_values))
    for nu in eachindex(singular_values), mu in eachindex(singular_values)
        sigma_product = singular_values[mu] * singular_values[nu]
        denominator = regularization_value + sigma_product^2
        transformed_delta[mu, nu] =
            denominator == 0.0 ? 0.0 :
            transformed_target[mu, nu] * sigma_product / denominator
    end
    diagnostic_cutoff =
        isempty(singular_values) ? 0.0 : max(size(product)...) * eps(Float64) * maximum(singular_values)
    rank = count(value -> value > diagnostic_cutoff, singular_values)
    return (
        delta = _hc_symmetrize(left_vectors * transformed_delta * transpose(left_vectors)),
        singular_values = singular_values,
        rank = rank,
        diagnostic_cutoff = diagnostic_cutoff,
    )
end

function _egoi_gaussian_pair_transform(Ctarget::AbstractMatrix{<:Real})
    coefficients = Matrix{Float64}(Ctarget)
    gaussian_count, target_count = size(coefficients)
    transform = zeros(Float64, gaussian_count^2, target_count^2)
    for j in 1:target_count, i in 1:target_count
        column = (j - 1) * target_count + i
        for q in 1:gaussian_count, p in 1:gaussian_count
            row = gaussian_coulomb_pair_index(p, q, gaussian_count)
            transform[row, column] = coefficients[p, i] * coefficients[q, j]
        end
    end
    return transform
end

"""
    egoi_target_coulomb_matrix(pair_coulomb, Ctarget)

Transform a dense Gaussian pair Coulomb matrix into target-orbital pair space.
`pair_coulomb` uses the [`gaussian_coulomb_pair_index`](@ref) convention, and
`Ctarget` columns are target orbitals in that Gaussian basis.
"""
function egoi_target_coulomb_matrix(
    pair_coulomb::AbstractMatrix{<:Real},
    Ctarget::AbstractMatrix{<:Real},
)
    pair_matrix = _hc_symmetrize(pair_coulomb)
    pair_count = size(pair_matrix, 1)
    size(pair_matrix, 2) == pair_count ||
        throw(DimensionMismatch("pair_coulomb must be square"))
    gaussian_count = round(Int, sqrt(pair_count))
    gaussian_count^2 == pair_count || throw(
        DimensionMismatch("pair_coulomb dimension must be a perfect square Gaussian pair count"),
    )
    size(Ctarget, 1) == gaussian_count || throw(
        DimensionMismatch("Ctarget row count must match the Gaussian orbital count"),
    )
    transform = _egoi_gaussian_pair_transform(Ctarget)
    return _hc_symmetrize(transpose(transform) * pair_matrix * transform)
end

"""
    egoi_density_density_correction(V, Qtarget, exact_target; regularization=1e-18)

Compute a symmetric density-density correction `ΔV` so the target-orbital pair
space reproduces `exact_target` as closely as possible:
`P' * (V + ΔV) * P ≈ exact_target`, where `P =
egoi_target_product_matrix(Qtarget)`.

The correction follows the CR2 smooth SVD regularization in the product basis,
so rank-deficient product spaces are handled explicitly and reported in
diagnostics.
"""
function egoi_density_density_correction(
    V::AbstractMatrix{<:Real},
    Qtarget::AbstractMatrix{<:Real},
    exact_target::AbstractMatrix{<:Real};
    regularization::Real = 1.0e-18,
)
    interaction = _hc_symmetrize(V)
    size(interaction, 1) == size(Qtarget, 1) || throw(
        DimensionMismatch("V dimension must match the target working-basis row count"),
    )
    product = egoi_target_product_matrix(Qtarget)
    target = _hc_symmetrize(exact_target)
    size(target) == (size(product, 2), size(product, 2)) || throw(
        DimensionMismatch("exact_target must have one row/column per target orbital pair"),
    )
    initial_target = _hc_symmetrize(transpose(product) * interaction * product)
    target_residual_before = initial_target - target
    delta_target = target - initial_target
    svd_delta = _egoi_svd_delta_matrix(product, delta_target; regularization)
    delta_v = svd_delta.delta
    corrected = _hc_symmetrize(interaction + delta_v)
    corrected_target = _hc_symmetrize(transpose(product) * corrected * product)
    target_residual_after = corrected_target - target
    delta_diagnostics = _hc_matrix_delta_diagnostics(delta_v, interaction)
    diagnostics = (
        target_residual_fro_before = _hc_fro_norm(target_residual_before),
        target_residual_max_before = _hc_max_abs(target_residual_before),
        target_residual_fro_after = _hc_fro_norm(target_residual_after),
        target_residual_max_after = _hc_max_abs(target_residual_after),
        delta_v_max = delta_diagnostics.max,
        delta_v_rms = delta_diagnostics.rms,
        delta_v_fro = delta_diagnostics.fro,
        delta_v_relative_max = delta_diagnostics.relative_max,
        delta_v_relative_fro = delta_diagnostics.relative_fro,
        product_rank = svd_delta.rank,
        product_singular_values = svd_delta.singular_values,
        product_min_singular_value = isempty(svd_delta.singular_values) ? 0.0 : minimum(svd_delta.singular_values),
        product_max_singular_value = isempty(svd_delta.singular_values) ? 0.0 : maximum(svd_delta.singular_values),
        regularization = Float64(regularization),
        diagnostic_singular_cutoff = svd_delta.diagnostic_cutoff,
    )
    return EGOIDensityDensityCorrectionResult(corrected, delta_v, diagnostics)
end

function _hc_quantile_sorted(sorted_values, q::Real)
    isempty(sorted_values) && return 0.0
    length(sorted_values) == 1 && return Float64(only(sorted_values))
    pos = 1 + (length(sorted_values) - 1) * Float64(q)
    lo, hi = floor(Int, pos), ceil(Int, pos)
    lo == hi && return Float64(sorted_values[lo])
    return Float64((hi - pos) * sorted_values[lo] + (pos - lo) * sorted_values[hi])
end
_hc_median(values) = _hc_quantile_sorted(sort(Float64.(values)), 0.5)
_hc_p95(values) = _hc_quantile_sorted(sort(Float64.(values)), 0.95)
_hc_distance(a, b) = sqrt(sum((Float64(a[i]) - Float64(b[i]))^2 for i in 1:3))

function _hc_target_channel(label)
    text = string(label)
    endswith(text, "s1") && return :s1
    endswith(text, "s2") && return :s2
    return :other
end

function _hc_project_supplement_to_protected_localized(coeffs, X, S_AA, geometry, loc)
    crg = CartesianResidualGaussians
    components = crg.protected_original_fixed_sector_components(geometry)
    qZ = transpose(components.Z) * S_AA * coeffs
    qQ = transpose(components.G_perp) * X * coeffs +
        transpose(components.A_perp) * S_AA * coeffs
    return transpose(loc.W) * vcat(qZ, qQ)
end

function _protected_retained_original_gto_egoi_target(
    supplement,
    S_AA,
    X,
    geometry,
    candidate_owner_indices;
    channels = (:s1, :s2),
    loc = nothing,
)
    allowed_channels = Set((:s1, :s2))
    channel_list = Symbol.(collect(channels))
    all(channel -> channel in allowed_channels, channel_list) ||
        throw(ArgumentError("retained-original EGOI currently supports only s1/s2 targets"))
    source_indices = Int.(collect(geometry.protected_original_source_indices))
    labels = String[String(supplement.orbitals[index].label) for index in source_indices]
    source_channels = [_hc_target_channel(label) for label in labels]
    keep = [i for i in eachindex(source_indices) if source_channels[i] in channel_list]
    isempty(keep) && throw(ArgumentError("retained-original EGOI target selection is empty"))
    selected_sources = source_indices[keep]
    selected_owners = Int.(candidate_owner_indices[selected_sources])
    owners = sort(unique(selected_owners))
    owner_counts = [count(==(owner), selected_owners) for owner in owners]
    length(unique(owner_counts)) == 1 ||
        throw(ArgumentError("retained-original EGOI targets must be owner-balanced"))
    for owner in owners
        found = Set(source_channels[keep[selected_owners .== owner]])
        Set(channel_list) == found || throw(ArgumentError(
            "retained-original EGOI target owner $(owner) is missing requested channels"))
    end
    nA = length(supplement.orbitals)
    norms = sqrt.(diag(S_AA))
    localized = isnothing(loc) ?
        CartesianResidualGaussians.protected_localized_inherited_site_transform(geometry) : loc
    Qtarget = zeros(Float64, size(localized.W, 2), length(selected_sources))
    projection_norms = zeros(Float64, length(selected_sources))
    for (column, source_index) in pairs(selected_sources)
        coeffs = zeros(Float64, nA)
        coeffs[source_index] = inv(norms[source_index])
        projected = _hc_project_supplement_to_protected_localized(
            coeffs, X, S_AA, geometry, localized)
        projection_norms[column] = norm(projected)
        projection_norms[column] > 1.0e-12 ||
            throw(ArgumentError("retained-original EGOI target projection is near zero"))
        Qtarget[:, column] .= projected ./ projection_norms[column]
    end
    orbitals = [supplement.orbitals[index] for index in selected_sources]
    pair_coulomb = gaussian_coulomb_pair_matrix(
        orbitals; max_orbitals = length(selected_sources))
    Ctarget = Matrix(Diagonal(1.0 ./ norms[selected_sources]))
    exact_ordered = egoi_target_coulomb_matrix(pair_coulomb, Ctarget)
    target_rows = NamedTuple[]
    for (target_index, source_index) in pairs(selected_sources)
        push!(target_rows, (;
            target_index,
            source_index,
            label = String(supplement.orbitals[source_index].label),
            channel = _hc_target_channel(supplement.orbitals[source_index].label),
            owner = selected_owners[target_index],
            center = supplement.orbitals[source_index].center,
            raw_norm = norms[source_index],
            projected_norm = projection_norms[target_index],
            projection_loss = max(0.0, 1.0 - projection_norms[target_index]^2)))
    end
    return (; Qtarget, exact_target_ordered = exact_ordered, target_rows,
        owners = selected_owners, labels = labels[keep], source_indices = selected_sources,
        channels = source_channels[keep], loc = localized)
end

function _hc_owner_local_pairs(owners)
    pairs = Tuple{Int,Int}[]
    for owner in sort(unique(owners))
        local_indices = findall(==(owner), owners)
        for j in eachindex(local_indices), i in 1:j
            push!(pairs, (local_indices[i], local_indices[j]))
        end
    end
    return pairs
end

function _hc_product_matrix(Qtarget, pairs)
    product = Matrix{Float64}(undef, size(Qtarget, 1), length(pairs))
    for (column, (i, j)) in enumerate(pairs)
        product[:, column] .= Qtarget[:, i] .* Qtarget[:, j]
    end
    return product
end

function _hc_exact_pair_block(exact_ordered, ntarget::Int, pairs)
    exact = Matrix{Float64}(undef, length(pairs), length(pairs))
    for b in eachindex(pairs), a in eachindex(pairs)
        i, j = pairs[a]
        p, q = pairs[b]
        exact[a, b] = exact_ordered[(j - 1) * ntarget + i, (q - 1) * ntarget + p]
    end
    return _hc_symmetrize(exact)
end

function _hc_local_scales(row_centers, row_owners, core_spacing::Real; nearest_count::Integer = 6)
    n = length(row_centers)
    length(row_owners) == n || throw(DimensionMismatch("row owner count mismatch"))
    ell = zeros(Float64, n)
    for owner in unique(row_owners)
        rows = findall(==(owner), row_owners)
        fallback_distances = Float64[]
        for a in rows, b in rows
            a < b || continue
            d = _hc_distance(row_centers[a], row_centers[b])
            d > 1.0e-12 && push!(fallback_distances, d)
        end
        fallback = isempty(fallback_distances) ? 1.0 : _hc_median(fallback_distances)
        for row in rows
            distances = sort(Float64[_hc_distance(row_centers[row], row_centers[other])
                for other in rows if other != row &&
                    _hc_distance(row_centers[row], row_centers[other]) > 1.0e-12])
            local_scale = isempty(distances) ? fallback :
                _hc_median(distances[1:min(Int(nearest_count), length(distances))])
            ell[row] = max(local_scale, Float64(core_spacing))
        end
    end
    return (; centers = row_centers, owners = Int.(row_owners), ell,
        core_spacing = Float64(core_spacing), nearest_count = Int(nearest_count))
end

function _hc_m2_variables(scales; factor::Real = 1.75)
    n = length(scales.ell)
    left = Int[]; right = Int[]
    local_ratio = Float64[]; distance = Float64[]; local_scale = Float64[]
    variable_class = Symbol[]
    for j in 1:n, i in 1:j
        scales.owners[i] == scales.owners[j] || continue
        d = i == j ? 0.0 : _hc_distance(scales.centers[i], scales.centers[j])
        scale = max(scales.ell[i], scales.ell[j], eps(Float64))
        d <= Float64(factor) * scale + 1.0e-12 || continue
        ratio = i == j ? 0.0 : d / scale
        push!(left, i); push!(right, j); push!(distance, d)
        push!(local_scale, scale); push!(local_ratio, ratio)
        push!(variable_class,
            i == j ? :diag :
            ratio <= 1.25 + 1.0e-12 ? :nearest_local :
            :next_local)
    end
    return (; left, right, local_ratio, distance, local_scale, variable_class)
end

function _hc_upper_triangle_vector(matrix)
    values = Float64[]
    for j in axes(matrix, 2), i in 1:j
        push!(values, matrix[i, j])
    end
    return values
end

function _hc_upper_triangle_design(product, variables)
    np = size(product, 2)
    design = Matrix{Float64}(undef, div(np * (np + 1), 2), length(variables.left))
    for variable_index in eachindex(variables.left)
        i, j = variables.left[variable_index], variables.right[variable_index]
        row = 1
        @inbounds for b in 1:np, a in 1:b
            design[row, variable_index] = i == j ?
                product[i, a] * product[i, b] :
                product[i, a] * product[j, b] + product[j, a] * product[i, b]
            row += 1
        end
    end
    return design
end

function _hc_relative_denominators(V, variables)
    diag_scale = _hc_median(abs.(diag(V)))
    floor_value = max(1.0e-12, 1.0e-8 * diag_scale)
    return Float64[
        max(abs(V[variables.left[i], variables.right[i]]), floor_value)
        for i in eachindex(variables.left)]
end

function _hc_local_cap_values(V, variables)
    fractions = Float64[
        variables.variable_class[i] === :diag ? 0.20 :
        variables.variable_class[i] === :nearest_local ? 0.20 :
        variables.variable_class[i] === :next_local ? 0.10 : 0.0
        for i in eachindex(variables.left)]
    return fractions .* _hc_relative_denominators(V, variables)
end

function _hc_delta_from_variables(n::Int, variables, x)
    delta = zeros(Float64, n, n)
    for index in eachindex(x)
        i, j = variables.left[index], variables.right[index]
        delta[i, j] = x[index]
        delta[j, i] = x[index]
    end
    return delta
end

function _hc_solve_scaled_ridge(design, target_vector, scales, free, lambda_relative)
    isempty(free) && return zeros(Float64, 0)
    D = @view design[:, free]
    scaled = D .* reshape(scales[free], 1, :)
    gram = scaled * transpose(scaled)
    singulars = sqrt.(sort(max.(0.0, eigvals(Symmetric(gram))); rev = true))
    maxsv = isempty(singulars) ? 1.0 : max(maximum(singulars), eps(Float64))
    y = transpose(scaled) * ((gram + Float64(lambda_relative * maxsv^2) *
        Matrix{Float64}(I, size(gram, 1), size(gram, 1))) \ target_vector)
    return scales[free] .* y
end

function _hc_active_set_clipped_ridge(design, target_vector, V, variables)
    scales = _hc_relative_denominators(V, variables)
    caps = _hc_local_cap_values(V, variables)
    relative_lambdas = sort(unique(vcat(10.0 .^ collect(-16.0:1.0:8.0),
        collect(0.01:0.0025:0.10))))
    scaled = design .* reshape(scales, 1, :)
    gram = scaled * transpose(scaled)
    singulars = sqrt.(sort(max.(0.0, eigvals(Symmetric(gram))); rev = true))
    maxsv = isempty(singulars) ? 1.0 : max(maximum(singulars), eps(Float64))
    identity_rows = Matrix{Float64}(I, size(gram, 1), size(gram, 1))
    best = nothing
    scan_rows = NamedTuple[]
    for lambda_relative in relative_lambdas
        lambda = Float64(lambda_relative * maxsv^2)
        y = transpose(scaled) * ((gram + lambda .* identity_rows) \ target_vector)
        x = scales .* y
        residual = design * x - target_vector
        rel = abs.(x) ./ scales
        row = (; lambda, lambda_relative, residual_fro_after = norm(residual),
            residual_reduction_percent = norm(target_vector) == 0.0 ? 0.0 :
                100.0 * (1.0 - norm(residual) / norm(target_vector)),
            delta_v_relative_max = isempty(rel) ? 0.0 : maximum(rel),
            delta_v_relative_p95 = _hc_p95(rel),
            delta_v_relative_median = _hc_median(rel))
        push!(scan_rows, row)
        if row.delta_v_relative_max <= 0.20 && row.delta_v_relative_p95 <= 0.10 &&
                (isnothing(best) ||
                    row.residual_reduction_percent > best.row.residual_reduction_percent)
            best = (; row, x)
        end
    end
    if isnothing(best)
        order = sortperm(eachindex(scan_rows);
            by = i -> (scan_rows[i].delta_v_relative_p95,
                scan_rows[i].delta_v_relative_max,
                -scan_rows[i].residual_reduction_percent))
        row = scan_rows[first(order)]
        y = transpose(scaled) * ((gram + row.lambda .* identity_rows) \ target_vector)
        best = (; row, x = scales .* y)
    end
    x = clamp.(best.x, .-caps, caps)
    frozen = abs.(best.x) .> caps
    lambda_relative = max(best.row.lambda_relative, 0.01)
    for _ in 1:6
        free = findall(!, frozen)
        isempty(free) && break
        residual_target = target_vector - design * x
        dx = _hc_solve_scaled_ridge(design, residual_target, scales, free, lambda_relative)
        current = norm(design * x - target_vector)
        best_x, best_frozen, best_residual = copy(x), copy(frozen), current
        for step in (1.0, 0.5, 0.25, 0.125)
            proposal, proposal_frozen = copy(x), copy(frozen)
            for (local_index, variable_index) in pairs(free)
                value = proposal[variable_index] + step * dx[local_index]
                cap = caps[variable_index]
                if abs(value) > cap
                    proposal[variable_index] = sign(value) * cap
                    proposal_frozen[variable_index] = true
                else
                    proposal[variable_index] = value
                end
            end
            residual = norm(design * proposal - target_vector)
            if residual < best_residual
                best_x, best_frozen, best_residual = proposal, proposal_frozen, residual
            end
        end
        x, frozen = best_x, best_frozen
        best_residual >= current - 1.0e-14 && break
    end
    residual = design * x - target_vector
    rel = abs.(x) ./ scales
    fit_row = (; lambda = NaN, lambda_relative,
        residual_fro_after = norm(residual),
        residual_reduction_percent = norm(target_vector) == 0.0 ? 0.0 :
            100.0 * (1.0 - norm(residual) / norm(target_vector)),
        delta_v_max = _hc_max_abs(x),
        delta_v_fro = norm(x),
        delta_v_relative_max = isempty(rel) ? 0.0 : maximum(rel),
        delta_v_relative_p95 = _hc_p95(rel),
        delta_v_relative_median = _hc_median(rel))
    return (; x, caps, scan_rows, singulars, row = fit_row,
        saturated = abs.(x) .>= 0.999 .* max.(caps, eps(Float64)))
end

function _hc_delta_class_rows(V, variables, x, caps)
    denominators = _hc_relative_denominators(V, variables)
    rel = abs.(x) ./ max.(denominators, eps(Float64))
    saturated = abs.(x) .>= 0.999 .* max.(caps, eps(Float64))
    rows = NamedTuple[]
    for class in (:diag, :nearest_local, :next_local)
        indices = [i for i in eachindex(x) if variables.variable_class[i] === class]
        isempty(indices) && continue
        push!(rows, (; variable_class = class, variable_count = length(indices),
            saturated_count = count(saturated[indices]),
            delta_abs_max = maximum(abs.(x[indices])),
            delta_fro = norm(x[indices]),
            relative_max = maximum(rel[indices]),
            relative_p95 = _hc_p95(rel[indices]),
            relative_median = _hc_median(rel[indices])))
    end
    return rows
end

function _hc_product_class(owners, pair)
    owners[pair[1]] == owners[pair[2]] || return :AB
    return owners[pair[1]] == minimum(owners) ? :AA : :BB
end

function _hc_residual_block_rows(residual_before, residual_after, pairs, owners)
    classes = [_hc_product_class(owners, pair) for pair in pairs]
    rows = NamedTuple[]
    for block in (:AA_AA, :BB_BB, :AA_BB)
        values_before = Float64[]; values_after = Float64[]
        for b in eachindex(pairs), a in eachindex(pairs)
            class_a, class_b = classes[a], classes[b]
            keep = block === :AA_AA ? class_a === :AA && class_b === :AA :
                block === :BB_BB ? class_a === :BB && class_b === :BB :
                class_a != class_b
            keep || continue
            push!(values_before, residual_before[a, b])
            push!(values_after, residual_after[a, b])
        end
        push!(rows, (; block, entry_count = length(values_before),
            residual_fro_before = norm(values_before),
            residual_fro_after = norm(values_after),
            residual_max_before = isempty(values_before) ? 0.0 : maximum(abs, values_before),
            residual_max_after = isempty(values_after) ? 0.0 : maximum(abs, values_after),
            residual_reduction_percent = norm(values_before) == 0.0 ? 0.0 :
                100.0 * (1.0 - norm(values_after) / norm(values_before))))
    end
    return rows
end

function _hc_aa_bb_subblock_rows(residual_before, residual_after, pairs, owners)
    classes = [_hc_product_class(owners, pair) for pair in pairs]
    diagkind(pair) = pair[1] == pair[2] ? :diag : :offdiag
    rows = NamedTuple[]
    for sub in (:diag_diag, :diag_offdiag, :offdiag_offdiag)
        before = Float64[]; after = Float64[]
        for b in eachindex(pairs), a in eachindex(pairs)
            classes[a] != classes[b] || continue
            kinds = sort((diagkind(pairs[a]), diagkind(pairs[b])))
            row_sub = kinds == (:diag, :diag) ? :diag_diag :
                kinds == (:diag, :offdiag) ? :diag_offdiag : :offdiag_offdiag
            row_sub === sub || continue
            push!(before, residual_before[a, b])
            push!(after, residual_after[a, b])
        end
        push!(rows, (; aa_bb_subblock = sub, entry_count = length(before),
            residual_fro_before = norm(before), residual_fro_after = norm(after),
            residual_max_before = isempty(before) ? 0.0 : maximum(abs, before),
            residual_max_after = isempty(after) ? 0.0 : maximum(abs, after),
            residual_reduction_percent = norm(before) == 0.0 ? 0.0 :
                100.0 * (1.0 - norm(after) / norm(before))))
    end
    return rows
end

function _hc_low_fock_shift(H, V_before, V_after, Qtarget, occupations, dense_limit)
    isnothing(H) && return (; status = :not_requested, before_min = NaN, after_min = NaN,
        shift_min = NaN)
    n = size(V_before, 1)
    n <= Int(dense_limit) || return (; status = :skipped_dimension,
        before_min = NaN, after_min = NaN, shift_min = NaN)
    density = projected_orbital_density(Qtarget, occupations)
    before = eigvals(Symmetric(density_density_restricted_fock(H, V_before, density)))
    after = eigvals(Symmetric(density_density_restricted_fock(H, V_after, density)))
    return (; status = :dense, before_min = before[1], after_min = after[1],
        shift_min = after[1] - before[1])
end

function _retained_original_gto_local_product_egoi(
    V::AbstractMatrix{<:Real},
    target;
    row_centers,
    row_owners,
    core_spacing::Real,
    H = nothing,
    occupations = nothing,
    low_fock_dense_limit::Integer = 2500,
)
    interaction = _hc_symmetrize(V)
    Qtarget = Matrix{Float64}(target.Qtarget)
    exact_ordered = _hc_symmetrize(target.exact_target_ordered)
    size(interaction, 1) == size(Qtarget, 1) ||
        throw(DimensionMismatch("V dimension must match retained-GTO target rows"))
    target_owners = Int.(target.owners)
    pairs = _hc_owner_local_pairs(target_owners)
    product = _hc_product_matrix(Qtarget, pairs)
    exact = _hc_exact_pair_block(exact_ordered, size(Qtarget, 2), pairs)
    residual_before = _hc_symmetrize(transpose(product) * interaction * product) - exact
    variables = _hc_m2_variables(_hc_local_scales(row_centers, row_owners, core_spacing))
    design = _hc_upper_triangle_design(product, variables)
    fit = _hc_active_set_clipped_ridge(
        design, _hc_upper_triangle_vector(-residual_before), interaction, variables)
    delta = _hc_delta_from_variables(size(interaction, 1), variables, fit.x)
    corrected = _hc_symmetrize(interaction + delta)
    residual_after = _hc_symmetrize(transpose(product) * corrected * product) - exact
    product_singulars = svdvals(product)
    rank_cutoff = max(size(product)...) * eps(Float64) * maximum(product_singulars)
    occ = occupations === nothing ? ones(Float64, size(Qtarget, 2)) : Float64.(occupations)
    low = _hc_low_fock_shift(H, interaction, corrected, Qtarget, occ, low_fock_dense_limit)
    diagnostics = (;
        convention = :retained_original_gto_local_product_m2,
        target_count = size(Qtarget, 2),
        product_count = size(product, 2),
        product_pairs = pairs,
        product_rank = count(>(rank_cutoff), product_singulars),
        product_singular_values = product_singulars,
        product_min_singular_value = minimum(product_singulars),
        product_median_singular_value = _hc_median(product_singulars),
        product_max_singular_value = maximum(product_singulars),
        allowed_variable_count = length(variables.left),
        design_singular_values = fit.singulars,
        residual_fro_before = norm(residual_before),
        residual_fro_after = norm(residual_after),
        residual_max_before = maximum(abs, residual_before),
        residual_max_after = maximum(abs, residual_after),
        residual_reduction_percent = norm(residual_before) == 0.0 ? 0.0 :
            100.0 * (1.0 - norm(residual_after) / norm(residual_before)),
        delta_v_matrix_relative_fro = norm(delta) / norm(interaction),
        delta_v_relative_max = fit.row.delta_v_relative_max,
        delta_v_relative_p95 = fit.row.delta_v_relative_p95,
        delta_v_relative_median = fit.row.delta_v_relative_median,
        saturated_count = count(fit.saturated),
        max_disallowed_delta_v = 0.0,
        corrected_v_finite = all(isfinite, corrected),
        corrected_v_symmetry_error = norm(corrected - transpose(corrected), Inf),
        target_rows = target.target_rows,
        residual_block_rows = _hc_residual_block_rows(
            residual_before, residual_after, pairs, target_owners),
        aa_bb_subblock_rows = _hc_aa_bb_subblock_rows(
            residual_before, residual_after, pairs, target_owners),
        delta_class_rows = _hc_delta_class_rows(interaction, variables, fit.x, fit.caps),
        low_fock = low)
    return (; interaction_matrix = corrected, interaction_delta = delta,
        Qtarget, exact_target = exact, local_product_matrix = product,
        variables, diagnostics)
end

function _protected_retained_original_gto_local_product_egoi(
    V::AbstractMatrix{<:Real},
    supplement,
    S_AA,
    X,
    geometry,
    candidate_owner_indices;
    row_centers,
    row_owners,
    core_spacing::Real,
    H = nothing,
    channels = (:s1, :s2),
    loc = nothing,
    occupations = nothing,
    low_fock_dense_limit::Integer = 2500,
)
    target = _protected_retained_original_gto_egoi_target(
        supplement, S_AA, X, geometry, candidate_owner_indices; channels, loc)
    return merge(
        _retained_original_gto_local_product_egoi(
            V, target; row_centers, row_owners, core_spacing, H, occupations,
            low_fock_dense_limit),
        (; target))
end

"""
    projected_orbital_density(Qtarget, occupations)

Return the full projected density matrix `sum_a occupations[a] * q_a*q_a'` in
an orthonormal working basis.
"""
function projected_orbital_density(
    Qtarget::AbstractMatrix{<:Real},
    occupations::AbstractVector{<:Real},
)
    q = Matrix{Float64}(Qtarget)
    occupation_values = Float64[Float64(value) for value in occupations]
    size(q, 2) == length(occupation_values) || throw(
        DimensionMismatch("occupation count must match target orbital count"),
    )
    density = zeros(Float64, size(q, 1), size(q, 1))
    for column in axes(q, 2)
        orbital = @view q[:, column]
        density .+= occupation_values[column] .* (orbital * transpose(orbital))
    end
    return _hc_symmetrize(density)
end

"""
    density_density_restricted_fock(H, V, D)

Build the density-density restricted Fock matrix
`H + Diagonal(V * diag(D)) - 0.5 .* (D .* V)` in the repo two-index
interaction convention. `D` is the full projected working-basis density matrix
from [`projected_orbital_density`](@ref).
"""
function density_density_restricted_fock(
    H::AbstractMatrix{<:Real},
    V::AbstractMatrix{<:Real},
    D::AbstractMatrix{<:Real},
)
    h = _hc_symmetrize(H)
    interaction = _hc_symmetrize(V)
    density = _hc_symmetrize(D)
    size(h) == size(interaction) == size(density) || throw(
        DimensionMismatch("H, V, and D dimensions must match"),
    )
    hartree = Diagonal(interaction * diag(density))
    exchange_like = 0.5 .* (density .* interaction)
    return _hc_symmetrize(h + hartree - exchange_like)
end

"""
    occupied_virtual_fock_residual(F, Qocc)

Return `(I - P) * F * P` in the Euclidean working-basis metric, where columns
of `Qocc` are QR-orthonormalized before forming `P = Q * Q'`.
"""
function occupied_virtual_fock_residual(
    F::AbstractMatrix{<:Real},
    Qocc::AbstractMatrix{<:Real},
)
    fock = _hc_symmetrize(F)
    q, _singular_values, _rank = _hc_qr_orthonormalize_columns(Qocc)
    size(fock, 1) == size(q, 1) || throw(
        DimensionMismatch("F dimension must match occupied orbital row count"),
    )
    projector = q * transpose(q)
    return Matrix{Float64}(fock * projector - projector * fock * projector)
end

"""
    stationary_fock_one_body_correction(F, Qocc)

Return a symmetric one-body correction that cancels the occupied-virtual Fock
residual for the QR-orthonormalized occupied subspace. This routine assumes an
orthonormal Euclidean working basis and intentionally has no general
metric-overlap API.
"""
function stationary_fock_one_body_correction(
    F::AbstractMatrix{<:Real},
    Qocc::AbstractMatrix{<:Real};
    reference_matrix::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
)
    fock = _hc_symmetrize(F)
    q, singular_values, target_rank = _hc_qr_orthonormalize_columns(Qocc)
    size(fock, 1) == size(q, 1) || throw(
        DimensionMismatch("F dimension must match occupied orbital row count"),
    )
    projector = q * transpose(q)
    residual_before = Matrix{Float64}(fock * projector - projector * fock * projector)
    delta_h = _hc_symmetrize(-(residual_before + transpose(residual_before)))
    corrected_fock = _hc_symmetrize(fock + delta_h)
    residual_after = Matrix{Float64}(
        corrected_fock * projector - projector * corrected_fock * projector,
    )
    baseline = reference_matrix === nothing ? fock : Matrix{Float64}(reference_matrix)
    delta_diagnostics = _hc_matrix_delta_diagnostics(delta_h, baseline)
    diagnostics = (
        occupied_virtual_residual_fro_before = _hc_fro_norm(residual_before),
        occupied_virtual_residual_max_before = _hc_max_abs(residual_before),
        occupied_virtual_residual_fro_after = _hc_fro_norm(residual_after),
        occupied_virtual_residual_max_after = _hc_max_abs(residual_after),
        delta_h_max = delta_diagnostics.max,
        delta_h_rms = delta_diagnostics.rms,
        delta_h_fro = delta_diagnostics.fro,
        delta_h_relative_max = delta_diagnostics.relative_max,
        delta_h_relative_fro = delta_diagnostics.relative_fro,
        target_rank = target_rank,
        target_orthonormality_error = norm(transpose(q) * q - I, Inf),
        target_singular_values = singular_values,
        target_min_singular_value = isempty(singular_values) ? 0.0 : minimum(singular_values),
        target_max_singular_value = isempty(singular_values) ? 0.0 : maximum(singular_values),
    )
    return StationaryFockCorrectionResult(corrected_fock, delta_h, diagnostics)
end

function _hc_occupied_target_columns(
    Qtarget::AbstractMatrix{<:Real},
    occupations::AbstractVector{<:Real},
)
    occupation_values = Float64[Float64(value) for value in occupations]
    occupied_indices = Int[
        index for (index, occupation) in pairs(occupation_values) if occupation > 0.0
    ]
    isempty(occupied_indices) &&
        throw(ArgumentError("at least one occupation must be positive for stationary correction"))
    return Matrix{Float64}(Qtarget)[:, occupied_indices]
end

"""
    egoi_stationary_hamiltonian_correction(H, V, Qtarget, exact_target, occupations;
        include_egoi=true, include_stationary=true, regularization=1e-18)

Compose the EGOI density-density interaction correction and the stationary-Fock
one-body correction on dense matrices in an orthonormal working basis.
"""
function egoi_stationary_hamiltonian_correction(
    H::AbstractMatrix{<:Real},
    V::AbstractMatrix{<:Real},
    Qtarget::AbstractMatrix{<:Real},
    exact_target::AbstractMatrix{<:Real},
    occupations::AbstractVector{<:Real};
    include_egoi::Bool = true,
    include_stationary::Bool = true,
    regularization::Real = 1.0e-18,
)
    h = _hc_symmetrize(H)
    interaction = _hc_symmetrize(V)
    size(h) == size(interaction) ||
        throw(DimensionMismatch("H and V must have matching dimensions"))
    size(h, 1) == size(Qtarget, 1) ||
        throw(DimensionMismatch("Qtarget row count must match H/V dimension"))
    length(occupations) == size(Qtarget, 2) ||
        throw(DimensionMismatch("occupation count must match target orbital count"))

    egoi_result =
        include_egoi ?
        egoi_density_density_correction(
            interaction,
            Qtarget,
            exact_target;
            regularization,
        ) :
        nothing
    corrected_v = egoi_result === nothing ? interaction : egoi_result.interaction_matrix
    interaction_delta =
        egoi_result === nothing ? zeros(Float64, size(interaction)) : egoi_result.interaction_delta

    density = projected_orbital_density(Qtarget, occupations)
    fock = density_density_restricted_fock(h, corrected_v, density)
    stationary_result =
        include_stationary ?
        stationary_fock_one_body_correction(
            fock,
            _hc_occupied_target_columns(Qtarget, occupations);
            reference_matrix = h,
        ) :
        nothing
    one_body_delta =
        stationary_result === nothing ? zeros(Float64, size(h)) : stationary_result.one_body_delta
    corrected_h = _hc_symmetrize(h + one_body_delta)
    diagnostics = (
        include_egoi = include_egoi,
        include_stationary = include_stationary,
        density = density,
        egoi = egoi_result === nothing ? nothing : egoi_result.diagnostics,
        stationary = stationary_result === nothing ? nothing : stationary_result.diagnostics,
    )
    return HamiltonianCorrectionResult(
        corrected_h,
        corrected_v,
        one_body_delta,
        interaction_delta,
        diagnostics,
    )
end

"""
    ordinary_cartesian_projected_gaussian_target(operators, supplement, coefficients;
        indices=nothing, occupations=nothing, expansion=coulomb_gaussian_expansion(doacc=false),
        max_orbitals=64)

Project selected Gaussian target orbitals into an ordinary Cartesian working
basis, normalize each projected column, and build their dense Gaussian exact
Coulomb target matrix. Diagnostics retain the raw projected Gram/column norms
and the normalized Gram so projection quality remains visible. This is a
convenience adapter for the shared matrix correction layer, not a separate
ordinary-specific correction model.
"""
function ordinary_cartesian_projected_gaussian_target(
    operators::OrdinaryCartesianOperators3D,
    supplement,
    coefficients::AbstractMatrix{<:Real};
    indices = nothing,
    occupations = nothing,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    max_orbitals = 64,
)
    coefficient_matrix = Matrix{Float64}(coefficients)
    selected_indices =
        indices === nothing ? collect(axes(coefficient_matrix, 2)) : Int[Int(index) for index in indices]
    isempty(selected_indices) &&
        throw(ArgumentError("ordinary projected correction target requires at least one selected index"))
    all(index -> 1 <= index <= size(coefficient_matrix, 2), selected_indices) || throw(
        ArgumentError("target coefficient indices must be valid coefficient matrix columns"),
    )
    target_coefficients = coefficient_matrix[:, selected_indices]
    overlap = gto_overlap_matrix(operators, supplement)
    size(overlap, 2) == size(target_coefficients, 1) || throw(
        DimensionMismatch("Gaussian target coefficient row count must match supplement orbital count"),
    )
    raw_projected_orbitals = Matrix{Float64}(overlap * target_coefficients)
    raw_projected_gram = transpose(raw_projected_orbitals) * raw_projected_orbitals
    raw_projected_norms = Float64[
        norm(view(raw_projected_orbitals, :, column)) for column in axes(raw_projected_orbitals, 2)
    ]
    projected_orbitals = copy(raw_projected_orbitals)
    for column in axes(projected_orbitals, 2)
        raw_projected_norms[column] > 1.0e-14 || throw(
            ArgumentError("projected Gaussian target column $(column) has near-zero norm"),
        )
        projected_orbitals[:, column] ./= raw_projected_norms[column]
    end
    occupation_values =
        occupations === nothing ?
        ones(Float64, length(selected_indices)) :
        Float64[Float64(value) for value in occupations]
    length(occupation_values) == length(selected_indices) || throw(
        DimensionMismatch("occupation count must match selected target count"),
    )
    pair_coulomb = gaussian_coulomb_pair_matrix(
        supplement;
        expansion,
        max_orbitals,
    )
    exact_target = egoi_target_coulomb_matrix(pair_coulomb, target_coefficients)
    normalized_projected_gram = transpose(projected_orbitals) * projected_orbitals
    normalized_projected_norms = Float64[
        norm(view(projected_orbitals, :, i)) for i in axes(projected_orbitals, 2)
    ]
    diagnostics = (
        selected_indices = Tuple(selected_indices),
        gaussian_orbital_count = size(coefficient_matrix, 1),
        target_count = length(selected_indices),
        raw_projected_gram = raw_projected_gram,
        raw_projected_column_norms = raw_projected_norms,
        normalized_projected_gram = normalized_projected_gram,
        normalized_projected_column_norms = normalized_projected_norms,
        projected_column_norms = normalized_projected_norms,
        projected_overlap_error = norm(normalized_projected_gram - I, Inf),
        dense_pair_dimension = size(pair_coulomb, 1),
    )
    return OrdinaryProjectedHamiltonianCorrectionTarget(
        projected_orbitals,
        target_coefficients,
        exact_target,
        occupation_values,
        selected_indices,
        diagnostics,
    )
end

"""
    ordinary_cartesian_egoi_stationary_correction(operators; target,
        nuclear_charges=nothing, include_egoi=true, include_stationary=true,
        regularization=1e-18, overlap_tol=1e-8)

Apply the shared matrix-level EGOI/stationary correction machinery to an
ordinary Cartesian QW operator payload. The result is a
[`HamiltonianCorrectionResult`](@ref), not a new operator payload.
"""
function ordinary_cartesian_egoi_stationary_correction(
    operators::OrdinaryCartesianOperators3D;
    target::OrdinaryProjectedHamiltonianCorrectionTarget,
    nuclear_charges = nothing,
    include_egoi::Bool = true,
    include_stationary::Bool = true,
    regularization::Real = 1.0e-18,
    overlap_tol::Real = 1.0e-8,
)
    overlap_error = norm(operators.overlap - I, Inf)
    overlap_error <= Float64(overlap_tol) || throw(
        ArgumentError(
            "ordinary_cartesian_egoi_stationary_correction requires an orthonormal final ordinary basis; got overlap error $(overlap_error)",
        ),
    )
    charges = nuclear_charges === nothing ? operators.nuclear_charges : nuclear_charges
    h = assembled_one_body_hamiltonian(operators; nuclear_charges = charges)
    result = egoi_stationary_hamiltonian_correction(
        h,
        operators.interaction_matrix,
        target.projected_orbitals,
        target.exact_target,
        target.occupations;
        include_egoi,
        include_stationary,
        regularization,
    )
    diagnostics = merge(
        result.diagnostics,
        (
            ordinary_adapter = (
                overlap_error = overlap_error,
                nuclear_charges = charges === nothing ? nothing : Tuple(Float64.(collect(charges))),
            ),
            target = target.diagnostics,
        ),
    )
    return HamiltonianCorrectionResult(
        result.one_body_hamiltonian,
        result.interaction_matrix,
        result.one_body_delta,
        result.interaction_delta,
        diagnostics,
    )
end
