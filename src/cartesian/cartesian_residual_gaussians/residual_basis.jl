struct CartesianResidualGaussianBasis
    base_dimension::Int
    candidate_count::Int
    residual_dimension::Int
    candidate_labels::Vector{String}
    candidate_owner_indices::Vector{Int}
    candidate_centers::Vector{NTuple{3,Float64}}
    compact_source_candidate_indices::Union{Nothing,Vector{Int}}
    residual_source_owner_indices::Vector{Int}
    residual_occupations::Vector{Float64}
    owner_retained_counts::Vector{Int}
    residual_labels::Vector{String}
    T_G::Matrix{Float64}
    T_A::Matrix{Float64}
    injected_G::Union{Nothing,Matrix{Float64}}
    injected_A::Union{Nothing,Matrix{Float64}}
    occupation_cutoff::Float64
    residual_injection_cutoff::Float64
    tau_neg_abs::Float64
    tau_neg_rel::Float64
    tau_merge_abs::Float64
    tau_merge_rel::Float64
    selection_rule::Symbol
    orientation::Symbol
    sign_rule::Symbol
end
residual_gaussian_center(center) = (Float64(center[1]), Float64(center[2]), Float64(center[3]))
residual_gaussian_float_centers(centers) = NTuple{3,Float64}[residual_gaussian_center(center) for center in centers]
residual_gaussian_candidate_labels(supplement) = String[String(orbital.label) for orbital in supplement.orbitals]
residual_gaussian_candidate_centers(supplement) = NTuple{3,Float64}[residual_gaussian_center(orbital.center) for orbital in supplement.orbitals]
terminal_residual_mixed_overlap(basis, bundles, supplement, block_builder, expansion) = block_builder(basis, bundles, supplement, (), expansion).mixed.overlap
_rg_sym(matrix) = Matrix{Float64}(0.5 .* (matrix .+ transpose(matrix)))
_rg_distance(a, b) = sqrt(sum((a[i] - b[i])^2 for i in 1:3))
function _rg_median(values)
    sorted = sort(values)
    mid = length(sorted) ÷ 2
    return isodd(length(sorted)) ? sorted[mid + 1] : 0.5 * (sorted[mid] + sorted[mid + 1])
end
function _rg_weighted_contracted_tail(orbital, distance)
    coefficients = abs.(Float64.(orbital.coefficients))
    weights = coefficients ./ max(sum(coefficients), eps(Float64))
    return sum(weights .* exp.(-Float64.(orbital.exponents) .* distance^2))
end
function _rg_weighted_width(orbital)
    coefficients = abs.(Float64.(orbital.coefficients))
    weights = coefficients ./ max(sum(coefficients), eps(Float64))
    return sum(weights .* inv.(sqrt.(Float64.(orbital.exponents))))
end
function residual_candidate_owner(center, nuclei)
    matches = findall(==(center), nuclei)
    length(matches) == 1 ||
        throw(ArgumentError("residual-Gaussian candidate center must exactly match one nucleus"))
    return first(matches)
end
function residual_gaussian_owner_centers(candidate_centers, candidate_owner_indices)
    owner_centers = Vector{Union{Nothing,NTuple{3,Float64}}}(nothing, maximum(candidate_owner_indices))
    for (center, owner) in zip(candidate_centers, candidate_owner_indices)
        if isnothing(owner_centers[owner])
            owner_centers[owner] = center
        else
            owner_centers[owner] == center ||
                throw(ArgumentError("residual-Gaussian compactness owner centers are inconsistent"))
        end
    end
    return NTuple{3,Float64}[center::NTuple{3,Float64} for center in owner_centers]
end
function residual_gaussian_candidate_tail_compactness_values(supplement, candidate_centers, candidate_owner_indices, metric)
    metric_symbol = Symbol(metric)
    owner_centers = residual_gaussian_owner_centers(candidate_centers, candidate_owner_indices)
    length(owner_centers) == 2 || throw(ArgumentError("residual-Gaussian tail compactness currently supports diatomics"))
    length(supplement.orbitals) == length(candidate_centers) ||
        throw(DimensionMismatch("residual-Gaussian compactness supplement count mismatch"))
    midpoint = ntuple(i -> 0.5 * (owner_centers[1][i] + owner_centers[2][i]), 3)
    distances = if metric_symbol == :midpoint_weighted_tail
        [_rg_distance(owner_centers[owner], midpoint) for owner in candidate_owner_indices]
    elseif metric_symbol == :other_center_weighted_tail
        [_rg_distance(owner_centers[owner], owner_centers[3 - owner]) for owner in candidate_owner_indices]
    else
        throw(ArgumentError("unknown residual-Gaussian compactness metric $(metric)"))
    end
    return Float64[_rg_weighted_contracted_tail(orbital, distance)
        for (orbital, distance) in zip(supplement.orbitals, distances)]
end
function residual_gaussian_compactness_keep_mask(residual_compactness, candidate_count, candidate_centers, candidate_owner_indices)
    isnothing(residual_compactness) && return trues(candidate_count)
    cutoff = Float64(getproperty(residual_compactness, :cutoff))
    values = hasproperty(residual_compactness, :values) ?
        Float64.(getproperty(residual_compactness, :values)) :
        residual_gaussian_candidate_tail_compactness_values(
            getproperty(residual_compactness, :supplement), candidate_centers,
            candidate_owner_indices, getproperty(residual_compactness, :metric))
    length(values) == candidate_count ||
        throw(DimensionMismatch("residual-Gaussian compactness value count mismatch"))
    all(isfinite, values) || throw(ArgumentError("residual-Gaussian compactness values must be finite"))
    return values .<= cutoff
end
function residual_gaussian_compactness_selector(residual_compactness)
    isnothing(residual_compactness) && return :owner_local_residual_occupation
    return hasproperty(residual_compactness, :selector) ?
        Symbol(getproperty(residual_compactness, :selector)) : :owner_local_residual_occupation
end
function residual_gaussian_order_values(residual_compactness, candidate_count, candidate_centers, candidate_owner_indices)
    isnothing(residual_compactness) &&
        throw(ArgumentError("ordered compact-first residual selection requires residual_compactness"))
    hasproperty(residual_compactness, :supplement) ||
        throw(ArgumentError("ordered compact-first residual selection requires residual_compactness.supplement"))
    supplement = getproperty(residual_compactness, :supplement)
    length(supplement.orbitals) == candidate_count ||
        throw(DimensionMismatch("residual-Gaussian compactness supplement count mismatch"))
    return (;
        width = Float64[_rg_weighted_width(orbital) for orbital in supplement.orbitals],
        midpoint_tail = residual_gaussian_candidate_tail_compactness_values(
            supplement, candidate_centers, candidate_owner_indices, :midpoint_weighted_tail),
        other_center_tail = residual_gaussian_candidate_tail_compactness_values(
            supplement, candidate_centers, candidate_owner_indices, :other_center_weighted_tail))
end
function canonicalize_residual_signs!(T_A, T_G)
    for column in axes(T_A, 2)
        row = argmax(abs.(view(T_A, :, column)))
        if T_A[row, column] < 0
            T_A[:, column] .*= -1
            T_G[:, column] .*= -1
        end
    end
    return T_A, T_G
end
residual_gaussian_overlap(T_G, T_A, X, S_AA) =
    transpose(T_G) * T_G + transpose(T_G) * X * T_A +
    transpose(T_A) * transpose(X) * T_G + transpose(T_A) * S_AA * T_A
residual_gaussian_inner(lG, lA, rG, rA, X, S_AA) =
    dot(lG, rG) + dot(lG, X * rA) + dot(lA, transpose(X) * rG) + dot(lA, S_AA * rA)
injection_complement(B, nG) = size(B, 2) == 0 ? Matrix{Float64}(I, nG, nG) :
    (qr(B).Q * I)[:, (size(B, 2) + 1):end]
injection_complement(residual) = injection_complement(residual.injected_G, residual.base_dimension)
injected_dimension(residual) = isnothing(residual.injected_G) ? 0 : size(residual.injected_G, 2)
function check_residual_gaussian_metric(values, tau_abs, tau_rel, label)
    tau = max(Float64(tau_abs), Float64(tau_rel) * max(maximum(values), 1.0))
    minimum(values) >= -tau ||
        throw(ArgumentError("$(label) has a negative eigenvalue beyond tolerance"))
    return tau
end
function owner_residual_gaussian_block(X, S_AA, owner_indices, owner, nA, cutoff,
    tau_neg_abs, tau_neg_rel; candidate_keep = trues(nA))
    indices = [index for index in findall(==(owner), owner_indices) if candidate_keep[index]]
    isempty(indices) && return nothing
    M = Matrix{Float64}(S_AA[indices, indices] - transpose(X[:, indices]) * X[:, indices])
    M = 0.5 .* (M .+ transpose(M))
    values, vectors = eigen(Symmetric(M))
    check_residual_gaussian_metric(
        values, tau_neg_abs, tau_neg_rel, "owner-local residual metric")
    keep = findall(>(cutoff), values)
    isempty(keep) && return nothing
    full_rank = length(keep) == length(values)
    natural = full_rank ? collect(eachindex(values)) :
        sort(keep; by = index -> (-values[index], index))
    transform = full_rank ? Matrix{Float64}(inv(sqrt(Symmetric(M)))) :
        Matrix{Float64}(vectors[:, natural] * Diagonal(1.0 ./ sqrt.(values[natural])))
    T_A = zeros(Float64, nA, length(natural))
    T_A[indices, :] .= transform
    labels = full_rank ?
        String["r$(owner)_$(column)_$(index)" for (column, index) in pairs(indices)] :
        String["r$(owner)_mode$(column)" for column in eachindex(natural)]
    return (; owner, T_A, occupations = Float64[values[natural]...], labels)
end
function owner_ordered_mgs_residual_gaussian_block(X, S_AA, owner_indices, owner, nA, cutoff,
    tau_neg_abs, tau_neg_rel, order_values; candidate_keep = trues(nA))
    indices = [index for index in findall(==(owner), owner_indices) if candidate_keep[index]]
    isempty(indices) && return nothing
    accepted_G = Vector{Vector{Float64}}()
    accepted_A = Vector{Vector{Float64}}()
    occupations = Float64[]
    labels = String[]
    source_indices = Int[]
    order = sort(indices; by = index ->
        (order_values.width[index], order_values.midpoint_tail[index],
            order_values.other_center_tail[index], index))
    for index in order
        G = -Vector{Float64}(X[:, index])
        A = zeros(Float64, nA); A[index] = 1.0
        for _pass in 1:2, column in eachindex(accepted_A)
            coeff = residual_gaussian_inner(
                accepted_G[column], accepted_A[column], G, A, X, S_AA)
            G .-= coeff .* accepted_G[column]
            A .-= coeff .* accepted_A[column]
        end
        residual_norm = residual_gaussian_inner(G, A, G, A, X, S_AA)
        check_residual_gaussian_metric(
            [residual_norm], tau_neg_abs, tau_neg_rel, "owner-local ordered residual metric")
        residual_norm = max(residual_norm, 0.0)
        residual_norm > cutoff || continue
        scale = inv(sqrt(residual_norm))
        G .*= scale; A .*= scale
        push!(accepted_G, G)
        push!(accepted_A, A)
        push!(occupations, residual_norm)
        push!(labels, "r$(owner)_mgs$(length(accepted_A))_$(index)")
        push!(source_indices, index)
    end
    isempty(accepted_A) && return nothing
    return (; owner, T_A = hcat(accepted_A...), occupations, labels, source_indices)
end
function injection_principal_block(X, S_AA, owner_indices, owner, nA, candidate_overlap_tol, tau_abs, tau_rel; candidate_keep = trues(nA))
    indices = [index for index in findall(==(owner), owner_indices) if candidate_keep[index]]
    isempty(indices) && return nothing
    values, vectors = eigen(Symmetric(_rg_sym(S_AA[indices, indices])))
    check_residual_gaussian_metric(values, tau_abs, tau_rel, "owner-local candidate metric")
    rank_tol = max(candidate_overlap_tol.atol, candidate_overlap_tol.rtol * maximum(values))
    keep = findall(>(rank_tol), values); isempty(keep) && return nothing
    P = zeros(Float64, nA, length(keep))
    P[indices, :] .= vectors[:, keep] * Diagonal(1.0 ./ sqrt.(values[keep]))
    C = X * P
    occupations, modes = eigen(Symmetric(_rg_sym(Matrix{Float64}(I, length(keep), length(keep)) - transpose(C) * C)))
    check_residual_gaussian_metric(occupations, tau_abs, tau_rel, "owner-local residual metric")
    return (; owner, modes = P * modes, occupations)
end
function injected_fixed_sector(X, S_AA, modes, nG, nA, atol, rtol)
    if isempty(modes)
        return (; injected_G = zeros(Float64, nG, 0), injected_A = zeros(Float64, nA, 0))
    end
    Y0 = hcat(modes...)
    values, vectors = eigen(Symmetric(_rg_sym(transpose(Y0) * S_AA * Y0)))
    check_residual_gaussian_metric(
        values, atol, rtol, "global injected-mode Gram")
    rank_tol = max(Float64(atol), Float64(rtol) * max(maximum(values), 1.0))
    keep = findall(>(rank_tol), values); isempty(keep) &&
        throw(ArgumentError("residual injection merge removed all injected modes"))
    Y = Matrix{Float64}(Y0 * vectors[:, keep] * Diagonal(1.0 ./ sqrt.(values[keep])))
    B = X * Y
    minimum(svdvals(B)) > rank_tol || throw(ArgumentError("residual injection projection is rank deficient"))
    nY = size(Y, 2); nY < nG ||
        throw(ArgumentError("residual injection subspace must be smaller than fixed sector"))
    return (; injected_G = B, injected_A = Y)
end

function _rg_supplement_orthonormal_basis(S_AA, atol, rtol)
    values, vectors = eigen(Symmetric(_rg_sym(S_AA)))
    check_residual_gaussian_metric(
        values, atol, rtol, "occupied-first supplement metric")
    rank_tol = max(Float64(atol), Float64(rtol) * max(maximum(values), 1.0))
    keep = findall(>(rank_tol), values)
    isempty(keep) && throw(ArgumentError("supplement metric has no retained rank"))
    return Matrix{Float64}(vectors[:, keep] * Diagonal(1.0 ./ sqrt.(values[keep]))),
        Float64[values...], keep
end

function _rg_capture_tolerance(S_AA, occupied_atol, occupied_rtol,
    supplement_atol, supplement_rtol)
    tolerances = Float64[occupied_atol, occupied_rtol,
        supplement_atol, supplement_rtol]
    all(isfinite, tolerances) && all(>=(0.0), tolerances) ||
        throw(ArgumentError("occupied-first overlap tolerances must be finite and nonnegative"))
    scale = max(opnorm(Matrix{Float64}(S_AA), Inf), 1.0)
    return max(max(tolerances[1], tolerances[3]),
        max(tolerances[2], tolerances[4]) * scale)
end

function _rg_validate_capture_eigenvalues(values, capture_tol, label)
    isempty(values) && return (NaN, NaN)
    raw_range = extrema(Float64.(values))
    raw_range[1] >= -capture_tol && raw_range[2] <= 1.0 + capture_tol ||
        throw(ArgumentError(
            "$(label) eigenvalues $(raw_range) fall outside physical range " *
            "[-$(capture_tol), $(1.0 + capture_tol)]"))
    return raw_range
end

function _rg_direction_metadata(labels, owners, channels, nA)
    out_labels = isnothing(labels) ? String["candidate_$(i)" for i in 1:nA] :
        String[String(label) for label in labels]
    length(out_labels) == nA ||
        throw(DimensionMismatch("occupied-first label count mismatch"))
    out_owners = isnothing(owners) ? zeros(Int, nA) : Int.(owners)
    length(out_owners) == nA ||
        throw(DimensionMismatch("occupied-first owner count mismatch"))
    out_channels = if isnothing(channels)
        String[split(replace(label, r"\d+$" => ""), "_")[end] for label in out_labels]
    else
        String[String(channel) for channel in channels]
    end
    length(out_channels) == nA ||
        throw(DimensionMismatch("occupied-first channel count mismatch"))
    return out_labels, out_owners, out_channels
end

function _rg_occupied_first_direction_rows(A, values, labels, owners, channels)
    rows = NamedTuple[]
    for column in axes(A, 2)
        weights = abs2.(view(A, :, column))
        total = sum(weights)
        index = total <= eps(Float64) ? 1 : argmax(weights)
        push!(rows, (;
            direction = column,
            value = Float64(values[column]),
            owner = owners[index],
            label = labels[index],
            channel = channels[index],
            coefficient_fraction = total <= eps(Float64) ? 0.0 :
                weights[index] / total,
        ))
    end
    return rows
end

"""
    occupied_first_injection_geometry(base_dimension, X, S_AA, Y_occ; ...)

Build an occupied-first injection geometry for supplement-space occupied
directions `Y_occ`. Occupied directions are mandatory; optional supplement
directions are selected only by capture into the represented span
`span(Y_occ) ⊕ (G ∩ Y_occ_perp)`. Weak optional directions are reported as
rejected and are not converted into residual/MWG channels.
"""
function occupied_first_injection_geometry(
    base_dimension::Integer,
    X,
    S_AA,
    Y_occ;
    optional_capture_cutoff::Real = 0.99,
    labels = nothing,
    owners = nothing,
    channels = nothing,
    provenance = (;),
    occupied_overlap_atol::Real = 1.0e-10,
    occupied_overlap_rtol::Real = 1.0e-10,
    supplement_overlap_atol::Real = 1.0e-12,
    supplement_overlap_rtol::Real = 1.0e-10,
)
    nG = Int(base_dimension)
    Xmat = Matrix{Float64}(X)
    S = Matrix{Float64}(S_AA)
    Y0 = Matrix{Float64}(Y_occ)
    nA = size(S, 1)
    size(S) == (nA, nA) ||
        throw(DimensionMismatch("occupied-first supplement metric must be square"))
    size(Xmat) == (nG, nA) ||
        throw(DimensionMismatch("occupied-first mixed overlap has wrong shape"))
    size(Y0, 1) == nA ||
        throw(DimensionMismatch("occupied-first occupied coefficients row count mismatch"))
    all(isfinite, Xmat) && all(isfinite, S) && all(isfinite, Y0) ||
        throw(ArgumentError("occupied-first inputs must be finite"))
    k = size(Y0, 2)
    k > 0 || throw(ArgumentError("occupied-first injection requires at least one occupied direction"))
    k < nG || throw(ArgumentError("occupied-first occupied subspace must be smaller than base dimension"))
    cutoff = Float64(optional_capture_cutoff)
    0.0 <= cutoff <= 1.0 + 1.0e-12 ||
        throw(ArgumentError("optional_capture_cutoff must be between 0 and 1"))
    label_data = _rg_direction_metadata(labels, owners, channels, nA)
    candidate_labels, candidate_owners, candidate_channels = label_data

    occupied_metric = _rg_sym(transpose(Y0) * S * Y0)
    occupied_orthogonality_error = maximum(abs, occupied_metric - I)
    occupied_orthogonality_error <= Float64(occupied_overlap_atol) ||
        throw(ArgumentError(
            "occupied-first Y_occ is not orthonormal in S_AA; error $(occupied_orthogonality_error)",
        ))
    capture_tol = _rg_capture_tolerance(S, occupied_overlap_atol,
        occupied_overlap_rtol, supplement_overlap_atol, supplement_overlap_rtol)
    complement_metric = _rg_sym(S - transpose(Xmat) * Xmat)
    complement_metric_values = eigvals(Symmetric(complement_metric))
    complement_metric_minimum_eigenvalue = minimum(complement_metric_values)
    complement_metric_minimum_eigenvalue >= -capture_tol ||
        throw(ArgumentError(
            "occupied-first mixed overlap has materially negative complement metric; " *
            "minimum eigenvalue $(complement_metric_minimum_eigenvalue) is below -$(capture_tol)"))
    occupied_base_capture_singular_values = Float64[svdvals(Xmat * Y0)...]
    occupied_base_capture_min =
        minimum(abs2.(occupied_base_capture_singular_values))
    fixed = injected_fixed_sector(
        Xmat,
        S,
        [Y0[:, column:column] for column in axes(Y0, 2)],
        nG,
        nA,
        occupied_overlap_atol,
        occupied_overlap_rtol,
    )
    mandatory_A = fixed.injected_A
    mandatory_G = fixed.injected_G
    Qp = injection_complement(mandatory_G, nG)
    recovered = transpose(mandatory_A) * S * Y0
    occupied_recovery_after_mandatory_inclusion_singular_values =
        Float64[svdvals(recovered)...]
    occupied_recovery_after_mandatory_inclusion_loss =
        maximum(abs.(1.0 .-
            occupied_recovery_after_mandatory_inclusion_singular_values))
    recovery_tol = max(Float64(occupied_overlap_atol), Float64(occupied_overlap_rtol))
    occupied_recovery_after_mandatory_inclusion_loss <= recovery_tol ||
        throw(ArgumentError(
            "occupied-first mandatory occupied recovery failed; loss " *
            "$(occupied_recovery_after_mandatory_inclusion_loss)",
        ))

    U, supplement_metric_values, supplement_keep = _rg_supplement_orthonormal_basis(
        S, supplement_overlap_atol, supplement_overlap_rtol)
    full_C = vcat(transpose(mandatory_A) * (S * U), transpose(Qp) * (Xmat * U))
    full_capture = _rg_sym(transpose(full_C) * full_C)
    raw_full_values = Float64[eigvals(Symmetric(full_capture))...]
    raw_full_capture_range = _rg_validate_capture_eigenvalues(
        raw_full_values, capture_tol, "occupied-first full capture")
    full_values = sort(raw_full_values; rev = true)

    U_perp0 = U - mandatory_A * (transpose(mandatory_A) * S * U)
    perp_values, perp_vectors = eigen(Symmetric(_rg_sym(transpose(U_perp0) * S * U_perp0)))
    perp_rank_tol = max(Float64(supplement_overlap_atol),
        Float64(supplement_overlap_rtol) * max(maximum(perp_values), 1.0))
    perp_keep = findall(>(perp_rank_tol), perp_values)
    Zperp = isempty(perp_keep) ? zeros(Float64, nA, 0) :
        Matrix{Float64}(U_perp0 * perp_vectors[:, perp_keep] *
            Diagonal(1.0 ./ sqrt.(perp_values[perp_keep])))
    Cperp = isempty(perp_keep) ? zeros(Float64, k + size(Qp, 2), 0) :
        vcat(transpose(mandatory_A) * (S * Zperp), transpose(Qp) * (Xmat * Zperp))
    capture_values = Float64[]
    capture_vectors = zeros(Float64, size(Zperp, 2), 0)
    if size(Zperp, 2) > 0
        values, vectors = eigen(Symmetric(_rg_sym(transpose(Cperp) * Cperp)))
        raw_complement_capture_range = _rg_validate_capture_eigenvalues(
            values, capture_tol, "occupied-first complement capture")
        order = sort(eachindex(values); by = index -> (-values[index], index))
        capture_values = Float64[values[index] for index in order]
        capture_vectors = Matrix{Float64}(vectors[:, order])
    else
        raw_complement_capture_range = (NaN, NaN)
    end
    kept = findall(>=(cutoff), capture_values)
    rejected = setdiff(eachindex(capture_values), kept)
    optional_A = isempty(kept) ? zeros(Float64, nA, 0) :
        Matrix{Float64}(Zperp * capture_vectors[:, kept])
    injected_A = hcat(mandatory_A, optional_A)
    final_rank = size(injected_A, 2)
    final_rank < nG ||
        throw(ArgumentError("occupied-first injected subspace must be smaller than base dimension"))
    injected_G = Matrix{Float64}(Xmat * injected_A)
    projection_singulars = svdvals(injected_G)
    rank_tol = max(Float64(occupied_overlap_atol),
        Float64(occupied_overlap_rtol) * max(maximum(projection_singulars), 1.0))
    minimum(projection_singulars) > rank_tol ||
        throw(ArgumentError("occupied-first injected projection is rank deficient"))
    kept_A = isempty(kept) ? zeros(Float64, nA, 0) :
        Matrix{Float64}(Zperp * capture_vectors[:, kept])
    rejected_A = isempty(rejected) ? zeros(Float64, nA, 0) :
        Matrix{Float64}(Zperp * capture_vectors[:, rejected])
    weakest_kept = isempty(kept) ? NaN : minimum(capture_values[kept])
    strongest_rejected = isempty(rejected) ? NaN : maximum(capture_values[rejected])
    diagnostics = (;
        provenance,
        occupied_orthogonality_error,
        occupied_base_capture_singular_values,
        occupied_base_capture_min,
        occupied_recovery_after_mandatory_inclusion_singular_values,
        occupied_recovery_after_mandatory_inclusion_loss,
        capture_tol,
        complement_metric_minimum_eigenvalue,
        raw_full_capture_range,
        raw_complement_capture_range,
        full_capture_eigenvalues = full_values,
        supplement_metric_eigenvalues = supplement_metric_values,
        supplement_metric_kept_indices = supplement_keep,
        complement_capture_eigenvalues = capture_values,
        optional_capture_cutoff = cutoff,
        optional_kept_count = length(kept),
        optional_rejected_count = length(rejected),
        weakest_kept_optional_capture = weakest_kept,
        strongest_rejected_weak_capture = strongest_rejected,
        final_rank,
        mandatory_count = size(mandatory_A, 2),
        rejected_weak_directions_not_mwg_residual_channels = true,
        kept_direction_summary = _rg_occupied_first_direction_rows(
            kept_A, capture_values[kept], candidate_labels, candidate_owners, candidate_channels),
        rejected_direction_summary = _rg_occupied_first_direction_rows(
            rejected_A, capture_values[rejected], candidate_labels, candidate_owners, candidate_channels),
    )
    return (;
        mandatory_A,
        optional_A,
        injected_A,
        injected_G,
        capture_eigenvalues = capture_values,
        full_capture_eigenvalues = full_values,
        optional_kept_indices = kept,
        optional_rejected_indices = rejected,
        diagnostics,
    )
end
function injected_owner_residual_block(block, fixed, X, S_AA, nA, cutoff, injection_cutoff, tau_abs, tau_rel)
    indices = findall(>(injection_cutoff), block.occupations); isempty(indices) && return nothing
    modes = block.modes[:, indices]
    Y, Qp = fixed.injected_A, injection_complement(fixed.injected_G, size(X, 1))
    C = vcat(transpose(Y) * (S_AA * modes), transpose(Qp) * (X * modes))
    values, vectors = eigen(Symmetric(_rg_sym(Matrix{Float64}(I, size(modes, 2), size(modes, 2)) - transpose(C) * C)))
    check_residual_gaussian_metric(values, tau_abs, tau_rel, "injected residual metric")
    keep = findall(>(cutoff), values); isempty(keep) && return nothing
    natural = sort(keep; by = index -> (-values[index], index))
    W = vectors[:, natural] * Diagonal(1.0 ./ sqrt.(values[natural]))
    CW = C * W
    nY = size(Y, 2)
    T_A = modes * W - Y * view(CW, 1:nY, :)
    T_G = -Qp * view(CW, (nY + 1):size(CW, 1), :)
    labels = String["r$(block.owner)_mode$(column)" for column in eachindex(natural)]
    return (; block.owner, T_G, T_A, occupations = Float64[values[natural]...], labels)
end
function residual_gaussian_block_metadata(blocks, owner_count)
    source_owners = Int[vcat((fill(block.owner, size(block.T_A, 2)) for block in blocks)...)...]
    occupations = Float64[vcat((block.occupations for block in blocks)...)...]
    labels = String[vcat((block.labels for block in blocks)...)...]
    return source_owners, occupations, [count(==(owner), source_owners) for owner in 1:owner_count], labels
end
function finalize_residual_gaussian_transform(T_G0, T_A0, X, S_AA, tau_merge_abs, tau_merge_rel, identity_atol)
    S_merge = _rg_sym(residual_gaussian_overlap(T_G0, T_A0, X, S_AA))
    merge_values = eigvals(Symmetric(S_merge))
    tau_merge = check_residual_gaussian_metric(
        merge_values, tau_merge_abs, tau_merge_rel, "residual-Gaussian final merge metric")
    minimum(merge_values) > tau_merge || throw(ArgumentError("residual-Gaussian final merge metric is near singular"))
    transform = Matrix{Float64}(inv(sqrt(Symmetric(S_merge))))
    T_A, T_G = Matrix{Float64}(T_A0 * transform), Matrix{Float64}(T_G0 * transform)
    canonicalize_residual_signs!(T_A, T_G)
    S_RR = _rg_sym(residual_gaussian_overlap(T_G, T_A, X, S_AA))
    identity_error, identity_scale = maximum(abs, S_RR - I), maximum(abs, S_RR)
    identity_error <= Float64(identity_atol) * (1.0 + max(1.0, identity_scale)) || throw(ArgumentError("residual-Gaussian R' S R validation failed"))
    return T_G, T_A
end
function protected_original_gram_clean(S_AA, indices, nA, atol, rtol)
    isempty(indices) && return (; Z = zeros(Float64, nA, 0), values = Float64[], keep = Int[])
    values, vectors = eigen(Symmetric(_rg_sym(S_AA[indices, indices])))
    rank_tol = max(Float64(atol), Float64(rtol) * maximum(values))
    keep = findall(>(rank_tol), values)
    isempty(keep) && throw(ArgumentError("protected-original candidate Gram cleanup removed all columns"))
    Z = zeros(Float64, nA, length(keep))
    Z[indices, :] .= vectors[:, keep] * Diagonal(1.0 ./ sqrt.(values[keep]))
    return (; Z, values = Float64[values...], keep)
end
function protected_original_complement(S_AA, source, fixed, atol, rtol;
    allow_empty::Bool = false)
    size(source, 2) == 0 && return (; Z = zeros(Float64, size(S_AA, 1), 0))
    Q0 = source - fixed * (transpose(fixed) * S_AA * source)
    values, vectors = eigen(Symmetric(_rg_sym(transpose(Q0) * S_AA * Q0)))
    rank_tol = max(Float64(atol), Float64(rtol) * maximum(values))
    keep = findall(>(rank_tol), values)
    isempty(keep) && !allow_empty &&
        throw(ArgumentError("protected-original complement Gram cleanup removed all columns"))
    Z = isempty(keep) ? zeros(Float64, size(S_AA, 1), 0) :
        Matrix{Float64}(Q0 * vectors[:, keep] * Diagonal(1.0 ./ sqrt.(values[keep])))
    return (; Z)
end
function protected_occupied_union(S_AA, occupied_blocks; overlap_atol::Real = 1.0e-10,
    overlap_rtol::Real = 1.0e-10)
    blocks = Matrix{Float64}[Matrix{Float64}(block) for block in occupied_blocks]
    isempty(blocks) && return nothing
    nA = size(S_AA, 1)
    all(block -> size(block, 1) == nA && size(block, 2) > 0, blocks) ||
        throw(DimensionMismatch("protected occupied blocks have incompatible dimensions"))
    tolerance = max(Float64(overlap_atol),
        Float64(overlap_rtol) * max(opnorm(Matrix{Float64}(S_AA), Inf), 1.0))
    maximum(norm(transpose(block) * S_AA * block - I, Inf) for block in blocks) <= tolerance ||
        throw(ArgumentError("protected occupied packet block is not S_AA-orthonormal"))
    all_columns = hcat(blocks...)
    gram = _rg_sym(transpose(all_columns) * S_AA * all_columns)
    values, vectors = eigen(Symmetric(gram))
    check_residual_gaussian_metric(values, overlap_atol, overlap_rtol,
        "protected occupied union Gram")
    rank_tol = max(Float64(overlap_atol),
        Float64(overlap_rtol) * max(maximum(values), 1.0))
    keep = findall(>(rank_tol), values)
    isempty(keep) && throw(ArgumentError("protected occupied union has no retained rank"))
    basis = Matrix{Float64}(all_columns * vectors[:, keep] *
        Diagonal(1.0 ./ sqrt.(values[keep])))
    recovery_singular_values = [Float64[svdvals(transpose(basis) * S_AA * block)...]
        for block in blocks]
    recovery_losses = [maximum(abs.(1.0 .- singulars))
        for singulars in recovery_singular_values]
    maximum(recovery_losses) <= tolerance ||
        throw(ArgumentError("protected occupied packet recovery failed"))
    return (; basis, gram_eigenvalues = Float64[values...], retained_rank = length(keep),
        recovery_singular_values, recovery_losses)
end
protected_original_projection_block(X, S_AA, T_G, T_A, Z) = begin
    XZ = X * Z
    vcat(XZ, transpose(T_G) * XZ + transpose(T_A) * (S_AA * Z))
end
function protected_original_fake_rdm_operator(S_AA, Z)
    W = S_AA * Z
    return _rg_sym(transpose(W) * (W ./ reshape(diag(S_AA), :, 1)))
end
function protected_original_m_errors(T_G, T_A, X, S_AA)
    GR = T_G + X * T_A
    RR = residual_gaussian_overlap(T_G, T_A, X, S_AA)
    return (; g_r_max = maximum(abs, GR), r_r_max = maximum(abs, RR - I))
end
function protected_original_q_sample_error(qrobj, nrows, start_column)
    start_column > nrows && return (; count = 0, max = 0.0, fro = 0.0)
    sample = sort(unique(vcat(collect(start_column:min(nrows, start_column + 9)),
        collect(max(start_column, nrows - 9):nrows))))
    E = zeros(Float64, nrows, length(sample))
    for (j, column) in pairs(sample)
        E[column, j] = 1.0
    end
    Q = Matrix(qrobj.Q * E)
    err = transpose(Q) * Q - I
    return (; count = length(sample), max = maximum(abs, err), fro = norm(err))
end
function protected_original_geometry_diagnostics(X, S_AA, T_G, T_A, Z_protected, Z, B)
    nZ, nM = size(Z, 2), size(B, 1)
    bsv = svdvals(B)
    zerr = transpose(Z) * S_AA * Z - I
    qrb = qr(B)
    qtb = Matrix(adjoint(qrb.Q) * B)
    qtail = qtb[(nZ + 1):end, :]
    qsample = protected_original_q_sample_error(qrb, nM, nZ + 1)
    merr = protected_original_m_errors(T_G, T_A, X, S_AA)
    nP = size(Z_protected, 2)
    protected_block = transpose(Z_protected) * S_AA * Z[:, 1:nP]
    protected_cross = nP == nZ ? zeros(Float64, nP, 0) :
        transpose(Z_protected) * S_AA * Z[:, (nP + 1):end]
    protected_svals = svdvals(protected_block)
    return (; b_singular_values = bsv, b_min = minimum(bsv), b_median = _rg_median(bsv),
        b_max = maximum(bsv), b_lt_0p999 = count(<(0.999), bsv),
        b_lt_0p99 = count(<(0.99), bsv), b_lt_0p98 = count(<(0.98), bsv),
        b_lt_0p95 = count(<(0.95), bsv), b_lt_0p9 = count(<(0.9), bsv),
        z_identity_max = maximum(abs, zerr), z_identity_fro = norm(zerr),
        z_m_qperp_max = maximum(abs, qtail), z_m_qperp_fro = norm(qtail),
        qperp_identity_sample_count = qsample.count,
        qperp_identity_sample_max = qsample.max,
        qperp_identity_sample_fro = qsample.fro,
        f_s_f_identity_block_max = maximum((maximum(abs, zerr), maximum(abs, qtail),
            qsample.max, merr.g_r_max, merr.r_r_max)),
        protected_span_min_sv = minimum(protected_svals),
        protected_span_max_sv = maximum(protected_svals),
        protected_identity_max = maximum(abs, protected_block - I),
        protected_cross_max = isempty(protected_cross) ? 0.0 : maximum(abs, protected_cross),
        m_g_r_max = merr.g_r_max, m_r_r_max = merr.r_r_max)
end
function protected_original_direction_summary(candidate_labels, candidate_owner_indices, A, values)
    rows = NamedTuple[]
    for column in axes(A, 2)
        weights = abs2.(view(A, :, column))
        index = argmax(weights)
        channel = split(replace(candidate_labels[index], r"\d+$" => ""), "_")[end]
        push!(rows, (; value = values[column], owner = candidate_owner_indices[index],
            label = candidate_labels[index], channel, coefficient_fraction =
                weights[index] / max(sum(weights), eps(Float64))))
    end
    return rows
end
function staged_protected_original_injection_geometry(base_dimension::Integer, X, S_AA,
    candidate_labels::Vector{String}, candidate_centers::Vector{NTuple{3,Float64}},
    candidate_owner_indices::Vector{Int}, residual::CartesianResidualGaussianBasis;
    occupied_blocks = Matrix{Float64}[], s_cut::Real = 0.95, occ_cut::Real = 0.003,
    candidate_overlap_atol::Real = 1.0e-12, candidate_overlap_rtol::Real = 1.0e-8,
    occupied_overlap_atol::Real = 1.0e-10, occupied_overlap_rtol::Real = 1.0e-10)
    nG, nA = Int(base_dimension), length(candidate_labels)
    size(X) == (nG, nA) || throw(DimensionMismatch("protected-original mixed overlap has wrong shape"))
    size(S_AA) == (nA, nA) || throw(DimensionMismatch("protected-original supplement overlap has wrong shape"))
    residual.base_dimension == nG && residual.candidate_count == nA ||
        throw(DimensionMismatch("protected-original compact residual dimensions do not match"))
    residual.candidate_labels == candidate_labels &&
        residual.candidate_owner_indices == candidate_owner_indices &&
        residual.candidate_centers == candidate_centers ||
        throw(ArgumentError("protected-original compact residual candidate metadata do not match"))
    protected_indices = residual.compact_source_candidate_indices
    isnothing(protected_indices) && throw(ArgumentError(
        "staged protected-original geometry requires compact residual source indices"))
    T_G, T_A = residual.T_G, residual.T_A
    Zcompact = protected_original_gram_clean(
        S_AA, protected_indices, nA, candidate_overlap_atol, candidate_overlap_rtol).Z
    occupied = protected_occupied_union(S_AA, occupied_blocks;
        overlap_atol = occupied_overlap_atol, overlap_rtol = occupied_overlap_rtol)
    Zp = if isnothing(occupied)
        Zcompact
    else
        compact_perp = protected_original_complement(S_AA, Zcompact, occupied.basis,
            candidate_overlap_atol, candidate_overlap_rtol; allow_empty = true)
        hcat(occupied.basis, compact_perp.Z)
    end
    if !isnothing(occupied)
        occupied_B = protected_original_projection_block(X, S_AA, T_G, T_A, occupied.basis)
        occupied_singulars = Float64[svdvals(occupied_B)...]
        minimum(occupied_singulars) >= Float64(s_cut) || throw(ArgumentError(
            "mandatory occupied union is insufficiently represented by compact M; minimum singular value $(minimum(occupied_singulars))"))
    else
        occupied_singulars = Float64[]
    end
    broad_indices = setdiff(collect(1:nA), protected_indices)
    E_broad = zeros(Float64, nA, length(broad_indices))
    for (column, index) in pairs(broad_indices); E_broad[index, column] = 1.0; end
    broad = protected_original_complement(
        S_AA, E_broad, Zp, candidate_overlap_atol, candidate_overlap_rtol)
    B_W = protected_original_projection_block(X, S_AA, T_G, T_A, broad.Z)
    fake_W = protected_original_fake_rdm_operator(S_AA, broad.Z)
    rep_values, rep_vectors = eigen(Symmetric(_rg_sym(transpose(B_W) * B_W)))
    rep_singulars = sqrt.(max.(rep_values, 0.0))
    keep_rep = sort(findall(>=(Float64(s_cut)), rep_singulars);
        by = index -> (-rep_singulars[index], index))
    isempty(keep_rep) && throw(ArgumentError("protected-original representability filter retained no columns"))
    Vrep = Matrix{Float64}(rep_vectors[:, keep_rep])
    Zrep = broad.Z * Vrep
    fake_rep = _rg_sym(transpose(Vrep) * fake_W * Vrep)
    fake_values, fake_vectors = eigen(Symmetric(fake_rep))
    fake_values = max.(fake_values, 0.0)
    keep_fake = sort(findall(>=(Float64(occ_cut)), fake_values);
        by = index -> (-fake_values[index], index))
    isempty(keep_fake) && throw(ArgumentError("protected-original fake-RDM filter retained no columns"))
    Zb = Zrep * Matrix{Float64}(fake_vectors[:, keep_fake])
    Z = hcat(Zp, Zb)
    B = protected_original_projection_block(X, S_AA, T_G, T_A, Z)
    diagnostics = protected_original_geometry_diagnostics(X, S_AA, T_G, T_A, Zp, Z, B)
    dropped_fake = setdiff(eachindex(fake_values), keep_fake)
    dropped_Z = Zrep * Matrix{Float64}(fake_vectors[:, dropped_fake])
    return merge(diagnostics, (;
        protected_original_count = size(Zp, 2),
        broad_gram_kept_count = size(broad.Z, 2),
        representability_kept_count = length(keep_rep),
        fake_rdm_kept_count = length(keep_fake),
        z_dimension = size(Z, 2),
        fake_rdm_trace_retained = sum(fake_values[keep_fake]; init = 0.0),
        fake_rdm_trace_dropped = sum(diag(fake_W)) - sum(fake_values[keep_fake]; init = 0.0),
        representability_singular_values = rep_singulars,
        dropped_representability_singular_values = rep_singulars[setdiff(eachindex(rep_singulars), keep_rep)],
        retained_fake_rdm_occupations = fake_values[keep_fake],
        dropped_fake_rdm_occupations = fake_values[dropped_fake],
        protected_original_source_indices = protected_indices,
        T_G, T_A,
        Z_protected = Zp,
        Z_broad = Zb,
        Z, B,
        additive_reference = isnothing(occupied) ? nothing : (;
            union_gram_eigenvalues = occupied.gram_eigenvalues,
            union_rank = occupied.retained_rank,
            packet_recovery_Z_singular_values = occupied.recovery_singular_values,
            packet_recovery_Z_losses = occupied.recovery_losses,
            occupied_representability_singular_values = occupied_singulars),
        dropped_fake_direction_summary = protected_original_direction_summary(
            candidate_labels, candidate_owner_indices, dropped_Z, fake_values[dropped_fake])))
end
function build_injected_residual_gaussian_basis(base_dimension, X, S_AA, labels, centers, owners, cutoff, injection_cutoff, candidate_overlap_atol, candidate_overlap_rtol, injected_overlap_atol, injected_overlap_rtol, tau_neg_abs, tau_neg_rel, tau_merge_abs, tau_merge_rel, orthogonality_atol, identity_atol, candidate_keep)
    injection_cutoff >= cutoff || throw(ArgumentError("residual_injection_cutoff must be at least residual_occupation_cutoff"))
    nG, nA = Int(base_dimension), length(labels)
    ctol = (; atol = Float64(candidate_overlap_atol), rtol = Float64(candidate_overlap_rtol))
    principal = filter(!isnothing, [injection_principal_block(X, S_AA, owners,
        owner, nA, ctol, tau_neg_abs, tau_neg_rel; candidate_keep) for owner in unique(owners)])
    injected = Matrix{Float64}[]
    for block in principal, index in findall(<=(injection_cutoff), block.occupations)
        push!(injected, block.modes[:, index:index])
    end
    fixed = injected_fixed_sector(X, S_AA, injected, nG, nA, injected_overlap_atol, injected_overlap_rtol)
    B, Y = fixed.injected_G, fixed.injected_A
    Qp = injection_complement(B, nG)
    YSY = transpose(Y) * S_AA * Y
    BQp = transpose(B) * Qp
    QpQp = transpose(Qp) * Qp
    y_error = isempty(YSY) ? 0.0 : maximum(abs, YSY - I)
    cross_error = isempty(BQp) ? 0.0 : maximum(abs, BQp)
    qperp_error = maximum(abs, QpQp - I)
    identity_error = max(y_error, cross_error, qperp_error)
    identity_scale = max(isempty(YSY) ? 0.0 : maximum(abs, YSY),
        isempty(BQp) ? 0.0 : maximum(abs, BQp), maximum(abs, QpQp))
    identity_error <= Float64(identity_atol) * (1.0 + max(1.0, identity_scale)) ||
        throw(ArgumentError("residual-Gaussian F' S F validation failed"))
    blocks = filter(!isnothing, [injected_owner_residual_block(block, fixed, X, S_AA, nA,
        cutoff, injection_cutoff, tau_neg_abs, tau_neg_rel) for block in principal])
    isempty(blocks) && throw(ArgumentError("residual-Gaussian candidate metric has no retained directions"))
    T_G0, T_A0 = hcat((b.T_G for b in blocks)...), hcat((b.T_A for b in blocks)...)
    T_G, T_A = finalize_residual_gaussian_transform(
        T_G0, T_A0, X, S_AA, tau_merge_abs, tau_merge_rel, identity_atol)
    FR = vcat(transpose(Y) * (transpose(X) * T_G + S_AA * T_A),
        transpose(Qp) * (T_G + X * T_A))
    norm(FR, Inf) <= orthogonality_atol ||
        throw(ArgumentError("residual-Gaussian F' S R validation failed"))
    source_owners, occupations, owner_counts, residual_labels =
        residual_gaussian_block_metadata(blocks, maximum(owners))
    return CartesianResidualGaussianBasis(nG, length(labels), size(T_A, 2), labels, owners,
        centers, nothing, source_owners, occupations, owner_counts, residual_labels, T_G, T_A,
        fixed.injected_G, fixed.injected_A, cutoff, injection_cutoff,
        Float64(tau_neg_abs), Float64(tau_neg_rel), Float64(tau_merge_abs), Float64(tau_merge_rel),
        :owner_local_residual_occupation, :injected_fixed_sector_owner_local_residual_final_merge_lowdin,
        :largest_T_A_entry_positive)
end
function build_residual_gaussian_basis(base_dimension::Integer, X, S_AA,
    candidate_labels::Vector{String}, candidate_centers::Vector{NTuple{3,Float64}},
    candidate_owner_indices::Vector{Int}; residual_occupation_cutoff::Real = 1.0e-6,
    residual_injection_cutoff::Real = 0.0, candidate_overlap_atol::Real = 1.0e-12,
    candidate_overlap_rtol::Real = 1.0e-8, injected_overlap_atol::Real = 1.0e-12,
    injected_overlap_rtol::Real = 1.0e-10,
    tau_neg_abs::Real = 1.0e-12, tau_neg_rel::Real = 1.0e-12,
    tau_merge_abs::Real = 1.0e-12, tau_merge_rel::Real = 1.0e-12,
    orthogonality_atol::Real = 1.0e-10, identity_atol::Real = 5.0e-8,
    residual_compactness = nothing)
    candidate_count = length(candidate_labels)
    length(candidate_centers) == candidate_count &&
        length(candidate_owner_indices) == candidate_count ||
        throw(DimensionMismatch("residual-Gaussian candidate metadata count mismatch"))
    size(X) == (Int(base_dimension), candidate_count) ||
        throw(DimensionMismatch("residual-Gaussian mixed overlap has wrong shape"))
    size(S_AA) == (candidate_count, candidate_count) ||
        throw(DimensionMismatch("residual-Gaussian supplement overlap has wrong shape"))
    cutoff = Float64(residual_occupation_cutoff)
    injection_cutoff = Float64(residual_injection_cutoff)
    selector = residual_gaussian_compactness_selector(residual_compactness)
    selector in (:owner_local_residual_occupation, :ordered_compact_first_mgs) ||
        throw(ArgumentError("unknown residual-Gaussian selector $(selector)"))
    candidate_keep = residual_gaussian_compactness_keep_mask(
        residual_compactness, candidate_count, candidate_centers, candidate_owner_indices)
    injection_cutoff > 0 && selector == :ordered_compact_first_mgs &&
        throw(ArgumentError("ordered compact-first residual selection does not support residual injection"))
    injection_cutoff > 0 && return build_injected_residual_gaussian_basis(
        base_dimension, X, S_AA, candidate_labels, candidate_centers,
        candidate_owner_indices, cutoff, injection_cutoff, candidate_overlap_atol,
        candidate_overlap_rtol, injected_overlap_atol, injected_overlap_rtol,
        tau_neg_abs, tau_neg_rel, tau_merge_abs, tau_merge_rel,
        orthogonality_atol, identity_atol, candidate_keep)
    owner_blocks, selection_rule, orientation = if selector == :ordered_compact_first_mgs
        order_values = residual_gaussian_order_values(
            residual_compactness, candidate_count, candidate_centers, candidate_owner_indices)
        blocks = filter(!isnothing, [owner_ordered_mgs_residual_gaussian_block(
            X, S_AA, candidate_owner_indices, owner, candidate_count, cutoff,
            tau_neg_abs, tau_neg_rel, order_values; candidate_keep)
            for owner in unique(candidate_owner_indices)])
        blocks, :owner_local_ordered_compact_first_mgs,
            :owner_local_ordered_mgs_final_merge_inverse_sqrt
    else
        blocks = filter(!isnothing, [owner_residual_gaussian_block(
            X, S_AA, candidate_owner_indices, owner, candidate_count, cutoff,
            tau_neg_abs, tau_neg_rel; candidate_keep)
            for owner in unique(candidate_owner_indices)])
        blocks, :owner_local_residual_occupation,
            :owner_local_residual_occupation_final_merge_lowdin
    end
    isempty(owner_blocks) &&
        throw(ArgumentError("residual-Gaussian candidate metric has no retained directions"))
    T_A0 = hcat((block.T_A for block in owner_blocks)...)
    T_G0 = Matrix{Float64}(-X * T_A0)
    T_G, T_A = finalize_residual_gaussian_transform(
        T_G0, T_A0, X, S_AA, tau_merge_abs, tau_merge_rel, identity_atol)
    norm(T_G + X * T_A, Inf) <= orthogonality_atol ||
        throw(ArgumentError("residual-Gaussian G' S R validation failed"))
    residual_source_owner_indices, residual_occupations, owner_counts, residual_labels =
        residual_gaussian_block_metadata(owner_blocks, maximum(candidate_owner_indices))
    compact_sources = selection_rule == :owner_local_ordered_compact_first_mgs ?
        sort(unique(vcat((block.source_indices for block in owner_blocks)...))) : nothing
    return CartesianResidualGaussianBasis(
        Int(base_dimension), candidate_count, size(T_A, 2), candidate_labels,
        candidate_owner_indices, candidate_centers, compact_sources, residual_source_owner_indices,
        residual_occupations, owner_counts, residual_labels, T_G, T_A,
        nothing, nothing, cutoff, 0.0, Float64(tau_neg_abs),
        Float64(tau_neg_rel), Float64(tau_merge_abs), Float64(tau_merge_rel),
        selection_rule, orientation, :largest_T_A_entry_positive)
end
