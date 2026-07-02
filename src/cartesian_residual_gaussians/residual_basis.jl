struct CartesianResidualGaussianBasis
    base_dimension::Int
    candidate_count::Int
    residual_dimension::Int
    candidate_labels::Vector{String}
    candidate_owner_indices::Vector{Int}
    candidate_centers::Vector{NTuple{3,Float64}}
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
function protected_original_broad_subspace(S_AA, broad_indices, protected_Z, atol, rtol)
    nA = size(S_AA, 1)
    E = zeros(Float64, nA, length(broad_indices))
    for (column, index) in pairs(broad_indices)
        E[index, column] = 1.0
    end
    Q0 = E - protected_Z * (transpose(protected_Z) * S_AA[:, broad_indices])
    values, vectors = eigen(Symmetric(_rg_sym(transpose(Q0) * S_AA * Q0)))
    rank_tol = max(Float64(atol), Float64(rtol) * maximum(values))
    keep = findall(>(rank_tol), values)
    isempty(keep) && throw(ArgumentError("protected-original broad Gram cleanup removed all columns"))
    return (; Z = Matrix{Float64}(Q0 * vectors[:, keep] * Diagonal(1.0 ./ sqrt.(values[keep]))),
        values = Float64[values...], keep)
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
    candidate_owner_indices::Vector{Int}; residual_occupation_cutoff::Real = 1.0e-6,
    residual_compactness, s_cut::Real = 0.95, occ_cut::Real = 0.003,
    candidate_overlap_atol::Real = 1.0e-12, candidate_overlap_rtol::Real = 1.0e-8,
    tau_neg_abs::Real = 1.0e-12, tau_neg_rel::Real = 1.0e-12,
    tau_merge_abs::Real = 1.0e-12, tau_merge_rel::Real = 1.0e-12,
    identity_atol::Real = 5.0e-8)
    nG, nA = Int(base_dimension), length(candidate_labels)
    size(X) == (nG, nA) || throw(DimensionMismatch("protected-original mixed overlap has wrong shape"))
    size(S_AA) == (nA, nA) || throw(DimensionMismatch("protected-original supplement overlap has wrong shape"))
    residual_gaussian_compactness_selector(residual_compactness) == :ordered_compact_first_mgs ||
        throw(ArgumentError("staged protected-original geometry requires ordered compact-first residual selection"))
    candidate_keep = residual_gaussian_compactness_keep_mask(
        residual_compactness, nA, candidate_centers, candidate_owner_indices)
    order_values = residual_gaussian_order_values(
        residual_compactness, nA, candidate_centers, candidate_owner_indices)
    blocks = filter(!isnothing, [owner_ordered_mgs_residual_gaussian_block(
        X, S_AA, candidate_owner_indices, owner, nA, Float64(residual_occupation_cutoff),
        tau_neg_abs, tau_neg_rel, order_values; candidate_keep)
        for owner in unique(candidate_owner_indices)])
    isempty(blocks) && throw(ArgumentError("protected-original compact RG selection retained no columns"))
    T_A0 = hcat((block.T_A for block in blocks)...)
    T_G0 = Matrix{Float64}(-X * T_A0)
    T_G, T_A = finalize_residual_gaussian_transform(
        T_G0, T_A0, X, S_AA, tau_merge_abs, tau_merge_rel, identity_atol)
    protected_indices = sort(unique(vcat((block.source_indices for block in blocks)...)))
    Zp = protected_original_gram_clean(
        S_AA, protected_indices, nA, candidate_overlap_atol, candidate_overlap_rtol).Z
    broad_indices = setdiff(collect(1:nA), protected_indices)
    broad = protected_original_broad_subspace(
        S_AA, broad_indices, Zp, candidate_overlap_atol, candidate_overlap_rtol)
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
    blocks = filter(!isnothing, [injected_owner_residual_block(block, fixed, X, S_AA, nA,
        cutoff, injection_cutoff, tau_neg_abs, tau_neg_rel) for block in principal])
    isempty(blocks) && throw(ArgumentError("residual-Gaussian candidate metric has no retained directions"))
    T_G0, T_A0 = hcat((b.T_G for b in blocks)...), hcat((b.T_A for b in blocks)...)
    T_G, T_A = finalize_residual_gaussian_transform(
        T_G0, T_A0, X, S_AA, tau_merge_abs, tau_merge_rel, identity_atol)
    Qp, Y = injection_complement(fixed.injected_G, nG), fixed.injected_A
    FR = vcat(transpose(Y) * (transpose(X) * T_G + S_AA * T_A),
        transpose(Qp) * (T_G + X * T_A))
    norm(FR, Inf) <= orthogonality_atol ||
        throw(ArgumentError("residual-Gaussian F' S R validation failed"))
    source_owners, occupations, owner_counts, residual_labels =
        residual_gaussian_block_metadata(blocks, maximum(owners))
    return CartesianResidualGaussianBasis(nG, length(labels), size(T_A, 2), labels, owners,
        centers, source_owners, occupations, owner_counts, residual_labels, T_G, T_A,
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
    return CartesianResidualGaussianBasis(
        Int(base_dimension), candidate_count, size(T_A, 2), candidate_labels,
        candidate_owner_indices, candidate_centers, residual_source_owner_indices,
        residual_occupations, owner_counts, residual_labels, T_G, T_A,
        nothing, nothing, cutoff, 0.0, Float64(tau_neg_abs),
        Float64(tau_neg_rel), Float64(tau_merge_abs), Float64(tau_merge_rel),
        selection_rule, orientation, :largest_T_A_entry_positive)
end
