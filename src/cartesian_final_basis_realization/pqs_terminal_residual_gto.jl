struct CartesianTerminalResidualGTOAugmentation
    base_dimension::Int
    candidate_count::Int
    residual_dimension::Int
    candidate_labels::Vector{String}
    candidate_owner_indices::Vector{Int}
    candidate_centers::Vector{NTuple{3,Float64}}
    retained_candidate_indices::Vector{Int}
    residual_labels::Vector{String}
    T_G::Matrix{Float64}
    T_A::Matrix{Float64}
    residual_metric_eigenvalues::Vector{Float64}
    tau_abs::Float64
    tau_rel::Float64
    tau_neg_abs::Float64
    tau_neg_rel::Float64
    rank_rule::Symbol
    orientation::Symbol
    sign_rule::Symbol
end

const _GB_PARENT = parentmodule(@__MODULE__)

function _residual_candidate_owner(center, nuclei)
    matches = findall(==(center), nuclei)
    length(matches) == 1 ||
        throw(ArgumentError("residual-GTO candidate center must exactly match one nucleus"))
    return first(matches)
end

function _residual_parent_representation(bundles)
    axis_representation = getfield(_GB_PARENT, :_cartesian_axis_representation)
    direct_product = getfield(_GB_PARENT, :_cartesian_direct_product_representation)
    return direct_product((x = axis_representation(_nested_axis_pgdg(bundles, :x).basis),
        y = axis_representation(_nested_axis_pgdg(bundles, :y).basis),
        z = axis_representation(_nested_axis_pgdg(bundles, :z).basis)))
end

function _terminal_residual_mixed_overlap(
    basis::CartesianTerminalBasisRealization,
    bundles,
    supplement,
)
    parent_representation = _residual_parent_representation(bundles)
    X = zeros(Float64, basis.final_dimension, length(supplement.orbitals))
    for block in basis.blocks
        local_overlap = gto_overlap_matrix(parent_representation, supplement;
            block_indices = block.support_indices)
        rows = block.column_range
        block.coefficients === nothing ?
            (X[rows, :] .= local_overlap) :
            mul!(view(X, rows, :), transpose(block.coefficients), local_overlap)
    end
    return X
end

function _canonicalize_residual_signs!(T_A, T_G)
    for column in axes(T_A, 2)
        row = argmax(abs.(view(T_A, :, column)))
        if T_A[row, column] < 0
            T_A[:, column] .*= -1
            T_G[:, column] .*= -1
        end
    end
    return T_A, T_G
end

function pqs_terminal_residual_gto_augmentation(
    basis::CartesianTerminalBasisRealization,
    bundles,
    supplement,
    nuclei;
    tau_abs::Real = 1.0e-10,
    tau_rel::Real = 1.0e-10,
    tau_neg_abs::Real = 1.0e-12,
    tau_neg_rel::Real = 1.0e-12,
    orthogonality_atol::Real = 1.0e-10,
    identity_atol::Real = 1.0e-10,
)
    nuclei_value = NTuple{3,Float64}[
        (Float64(nucleus[1]), Float64(nucleus[2]), Float64(nucleus[3]))
        for nucleus in nuclei
    ]
    labels = String[String(orbital.label) for orbital in supplement.orbitals]
    centers = NTuple{3,Float64}[
        (Float64(orbital.center[1]), Float64(orbital.center[2]), Float64(orbital.center[3]))
        for orbital in supplement.orbitals
    ]
    owners = Int[_residual_candidate_owner(center, nuclei_value) for center in centers]
    candidate_count = length(labels)
    X = _terminal_residual_mixed_overlap(basis, bundles, supplement)
    size(X) == (basis.final_dimension, candidate_count) ||
        throw(DimensionMismatch("residual-GTO mixed overlap has wrong shape"))
    S_AA = Matrix{Float64}(
        getfield(_GB_PARENT, :_cartesian_supplement_cross_overlap)(supplement, supplement))
    S_R = Symmetric(0.5 .* ((S_AA - transpose(X) * X) .+ transpose(S_AA - transpose(X) * X)))
    eigenvalues = eigvals(S_R)
    lambda_max = max(maximum(eigenvalues), 0.0)
    tau_keep = max(Float64(tau_abs), Float64(tau_rel) * lambda_max)
    tau_neg = max(Float64(tau_neg_abs), Float64(tau_neg_rel) * max(lambda_max, 1.0))
    minimum(eigenvalues) >= -tau_neg ||
        throw(ArgumentError("residual-GTO metric has a negative eigenvalue beyond tolerance"))
    retained = findall(>(tau_keep), eigenvalues)
    length(retained) == candidate_count ||
        throw(ArgumentError("residual-GTO rank-deficient candidate selection is not implemented in R3-A part 1"))
    T_A = Matrix{Float64}(inv(sqrt(S_R)))
    T_G = Matrix{Float64}(-X * T_A)
    _canonicalize_residual_signs!(T_A, T_G)
    norm(T_G + X * T_A, Inf) <= orthogonality_atol ||
        throw(ArgumentError("residual-GTO G' S R validation failed"))
    overlap = transpose(T_G) * T_G + transpose(T_G) * X * T_A +
              transpose(T_A) * transpose(X) * T_G + transpose(T_A) * S_AA * T_A
    norm(overlap - I, Inf) <= identity_atol ||
        throw(ArgumentError("residual-GTO R' S R validation failed"))
    retained_indices = collect(1:candidate_count)
    owner_counts = Dict{Int,Int}()
    residual_labels = String[]
    for index in retained_indices
        count = get(owner_counts, owners[index], 0) + 1
        owner_counts[owners[index]] = count
        push!(residual_labels, string("r", owners[index], "_", count, "_", labels[index]))
    end
    return CartesianTerminalResidualGTOAugmentation(
        basis.final_dimension, candidate_count, candidate_count, labels, owners, centers,
        retained_indices, residual_labels, T_G, T_A, Float64[eigenvalues...],
        Float64(tau_abs), Float64(tau_rel), Float64(tau_neg_abs), Float64(tau_neg_rel),
        :full_rank_candidate_order, :selected_candidate_order_symmetric_lowdin,
        :largest_T_A_entry_positive)
end
