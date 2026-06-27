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
    occupation_cutoff::Float64
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
function build_residual_gaussian_basis(base_dimension::Integer, X, S_AA,
    candidate_labels::Vector{String}, candidate_centers::Vector{NTuple{3,Float64}},
    candidate_owner_indices::Vector{Int}; residual_occupation_cutoff::Real = 5.0e-8,
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
        residual_occupations, owner_counts, residual_labels, T_G, T_A, cutoff,
        Float64(tau_neg_abs), Float64(tau_neg_rel), Float64(tau_merge_abs),
        Float64(tau_merge_rel), :owner_local_residual_occupation,
        :owner_local_residual_occupation_final_merge_lowdin, :largest_T_A_entry_positive)
end
