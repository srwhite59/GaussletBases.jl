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
    injected_dimension::Int
    injected_owner_counts::Vector{Int}
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
function residual_candidate_owner(center, nuclei)
    matches = findall(==(center), nuclei)
    length(matches) == 1 ||
        throw(ArgumentError("residual-Gaussian candidate center must exactly match one nucleus"))
    return first(matches)
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
injection_complement(B, nG) = size(B, 2) == 0 ? Matrix{Float64}(I, nG, nG) :
    (qr(B).Q * I)[:, (size(B, 2) + 1):end]
injection_complement(residual) = injection_complement(residual.injected_G, residual.base_dimension)
function residual_gaussian_identity_error(T_G, T_A, X, S_AA)
    S_RR = Matrix{Float64}(residual_gaussian_overlap(T_G, T_A, X, S_AA))
    S_RR = 0.5 .* (S_RR .+ transpose(S_RR))
    return maximum(abs, S_RR - I), maximum(abs, S_RR)
end
function check_residual_gaussian_metric(values, tau_abs, tau_rel, label)
    tau = max(Float64(tau_abs), Float64(tau_rel) * max(maximum(values), 1.0))
    minimum(values) >= -tau ||
        throw(ArgumentError("$(label) has a negative eigenvalue beyond tolerance"))
    return tau
end
function owner_residual_gaussian_block(X, S_AA, owner_indices, owner, nA, cutoff,
    tau_neg_abs, tau_neg_rel)
    indices = findall(==(owner), owner_indices)
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
function injection_principal_block(X, S_AA, owner_indices, owner, nA, candidate_overlap_tol, tau_abs, tau_rel)
    indices = findall(==(owner), owner_indices)
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
        return (; injected_G = zeros(Float64, nG, 0), injected_A = zeros(Float64, nA, 0), dimension = 0)
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
    return (; injected_G = B, injected_A = Y, dimension = nY)
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
function build_injected_residual_gaussian_basis(base_dimension, X, S_AA, labels, centers, owners, cutoff, injection_cutoff, candidate_overlap_atol, candidate_overlap_rtol, injected_overlap_atol, injected_overlap_rtol, tau_neg_abs, tau_neg_rel, tau_merge_abs, tau_merge_rel, orthogonality_atol, identity_atol)
    injection_cutoff >= cutoff || throw(ArgumentError("residual_injection_cutoff must be at least residual_occupation_cutoff"))
    nG, nA = Int(base_dimension), length(labels)
    ctol = (; atol = Float64(candidate_overlap_atol), rtol = Float64(candidate_overlap_rtol))
    principal = filter(!isnothing, [injection_principal_block(X, S_AA, owners,
        owner, nA, ctol, tau_neg_abs, tau_neg_rel) for owner in unique(owners)])
    injected, injected_counts = Matrix{Float64}[], zeros(Int, maximum(owners))
    for block in principal, index in findall(<=(injection_cutoff), block.occupations)
        push!(injected, block.modes[:, index:index]); injected_counts[block.owner] += 1
    end
    fixed = injected_fixed_sector(X, S_AA, injected, nG, nA, injected_overlap_atol, injected_overlap_rtol)
    blocks = filter(!isnothing, [injected_owner_residual_block(block, fixed, X, S_AA, nA,
        cutoff, injection_cutoff, tau_neg_abs, tau_neg_rel) for block in principal])
    isempty(blocks) && throw(ArgumentError("residual-Gaussian candidate metric has no retained directions"))
    T_G0, T_A0 = hcat((b.T_G for b in blocks)...), hcat((b.T_A for b in blocks)...)
    S_merge = _rg_sym(residual_gaussian_overlap(T_G0, T_A0, X, S_AA))
    merge_values = eigvals(Symmetric(S_merge))
    tau_merge = check_residual_gaussian_metric(
        merge_values, tau_merge_abs, tau_merge_rel, "residual-Gaussian final merge metric")
    minimum(merge_values) > tau_merge ||
        throw(ArgumentError("residual-Gaussian final merge metric is near singular"))
    transform = Matrix{Float64}(inv(sqrt(Symmetric(S_merge))))
    T_A, T_G = Matrix{Float64}(T_A0 * transform), Matrix{Float64}(T_G0 * transform)
    canonicalize_residual_signs!(T_A, T_G)
    Qp, Y = injection_complement(fixed.injected_G, nG), fixed.injected_A
    FR = vcat(transpose(Y) * (transpose(X) * T_G + S_AA * T_A),
        transpose(Qp) * (T_G + X * T_A))
    norm(FR, Inf) <= orthogonality_atol ||
        throw(ArgumentError("residual-Gaussian F' S R validation failed"))
    identity_error, identity_scale = residual_gaussian_identity_error(T_G, T_A, X, S_AA)
    identity_error <= Float64(identity_atol) * (1.0 + max(1.0, identity_scale)) ||
        throw(ArgumentError("residual-Gaussian R' S R validation failed"))
    source_owners = Int[vcat((fill(b.owner, size(b.T_A, 2)) for b in blocks)...)...]
    occupations = Float64[vcat((b.occupations for b in blocks)...)...]
    residual_labels = String[vcat((b.labels for b in blocks)...)...]
    owner_counts = [count(==(owner), source_owners) for owner in 1:maximum(owners)]
    return CartesianResidualGaussianBasis(nG, length(labels), size(T_A, 2), labels, owners,
        centers, source_owners, occupations, owner_counts, residual_labels, T_G, T_A,
        fixed.injected_G, fixed.injected_A, cutoff, injection_cutoff, fixed.dimension, injected_counts,
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
    orthogonality_atol::Real = 1.0e-10, identity_atol::Real = 5.0e-8)
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
    injection_cutoff > 0 && return build_injected_residual_gaussian_basis(
        base_dimension, X, S_AA, candidate_labels, candidate_centers,
        candidate_owner_indices, cutoff, injection_cutoff, candidate_overlap_atol,
        candidate_overlap_rtol, injected_overlap_atol, injected_overlap_rtol,
        tau_neg_abs, tau_neg_rel, tau_merge_abs, tau_merge_rel,
        orthogonality_atol, identity_atol)
    owner_blocks = filter(!isnothing, [owner_residual_gaussian_block(
        X, S_AA, candidate_owner_indices, owner, candidate_count, cutoff,
        tau_neg_abs, tau_neg_rel) for owner in unique(candidate_owner_indices)])
    isempty(owner_blocks) &&
        throw(ArgumentError("residual-Gaussian candidate metric has no retained directions"))
    T_A0 = hcat((block.T_A for block in owner_blocks)...)
    T_G0 = Matrix{Float64}(-X * T_A0)
    S_merge = Matrix{Float64}(residual_gaussian_overlap(T_G0, T_A0, X, S_AA))
    S_merge = 0.5 .* (S_merge .+ transpose(S_merge))
    merge_values = eigvals(Symmetric(S_merge))
    tau_merge = check_residual_gaussian_metric(
        merge_values, tau_merge_abs, tau_merge_rel, "residual-Gaussian final merge metric")
    minimum(merge_values) > tau_merge ||
        throw(ArgumentError("residual-Gaussian final merge metric is near singular"))
    transform = Matrix{Float64}(inv(sqrt(Symmetric(S_merge))))
    T_A = Matrix{Float64}(T_A0 * transform)
    T_G = Matrix{Float64}(T_G0 * transform)
    canonicalize_residual_signs!(T_A, T_G)
    norm(T_G + X * T_A, Inf) <= orthogonality_atol ||
        throw(ArgumentError("residual-Gaussian G' S R validation failed"))
    identity_tolerance = Float64(identity_atol)
    identity_error, identity_scale = residual_gaussian_identity_error(T_G, T_A, X, S_AA)
    identity_error <= identity_tolerance * (1.0 + max(1.0, identity_scale)) ||
        throw(ArgumentError("residual-Gaussian R' S R validation failed"))
    residual_source_owner_indices = Int[vcat((fill(block.owner, size(block.T_A, 2))
        for block in owner_blocks)...)...]
    residual_occupations = Float64[vcat((block.occupations for block in owner_blocks)...)...]
    residual_labels = String[vcat((block.labels for block in owner_blocks)...)...]
    owner_counts = [count(==(owner), residual_source_owner_indices)
        for owner in 1:maximum(candidate_owner_indices)]
    return CartesianResidualGaussianBasis(
        Int(base_dimension), candidate_count, size(T_A, 2), candidate_labels,
        candidate_owner_indices, candidate_centers, residual_source_owner_indices,
        residual_occupations, owner_counts, residual_labels, T_G, T_A,
        nothing, nothing, cutoff, 0.0, 0, Int[], Float64(tau_neg_abs),
        Float64(tau_neg_rel), Float64(tau_merge_abs), Float64(tau_merge_rel),
        :owner_local_residual_occupation,
        :owner_local_residual_occupation_final_merge_lowdin, :largest_T_A_entry_positive)
end
