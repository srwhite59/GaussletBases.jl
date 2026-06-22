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

function _r3a_dense(block, label)
    isnothing(block.dense_block) &&
        throw(ArgumentError("R3-A required CPB block $label is unavailable"))
    return Matrix{Float64}(block.dense_block)
end

function _r3a_local_index(ix, iy, iz, shape)
    return (ix - 1) * shape[2] * shape[3] + (iy - 1) * shape[3] + iz
end

function _r3a_bounding_cpb(states)
    CPB = getfield(_GB_PARENT, :CartesianCPB)
    ix = minimum(s -> s[1], states):maximum(s -> s[1], states)
    iy = minimum(s -> s[2], states):maximum(s -> s[2], states)
    iz = minimum(s -> s[3], states):maximum(s -> s[3], states)
    cpb = CPB.cpb(ix, iy, iz; role = :r3a_terminal_owned_support_bounding_cpb)
    shape = CPB.shape(cpb)
    rows = Int[
        _r3a_local_index(s[1] - first(ix) + 1, s[2] - first(iy) + 1,
            s[3] - first(iz) + 1, shape) for s in states
    ]
    return cpb, rows
end

function _r3a_mixed_block(basis, parent_basis_object, supplement, builder)
    nA = length(supplement.orbitals)
    O_GA = zeros(Float64, basis.final_dimension, nA)
    for block in basis.blocks
        cpb, rows = _r3a_bounding_cpb(block.support_states)
        local_block = _r3a_dense(builder(parent_basis_object, cpb, supplement), :mixed)[
            rows, :]
        block.coefficients === nothing ?
            (O_GA[block.column_range, :] .= local_block) :
            mul!(view(O_GA, block.column_range, :), transpose(block.coefficients), local_block)
    end
    return O_GA
end

function _r3a_augmented_operator(O_GG, O_GA, O_AA, residual)
    T_G, T_A = residual.T_G, residual.T_A
    O_GR = O_GG * T_G + O_GA * T_A
    O_RR = transpose(T_G) * O_GG * T_G +
           transpose(T_G) * O_GA * T_A +
           transpose(T_A) * transpose(O_GA) * T_G +
           transpose(T_A) * O_AA * T_A
    nG, nR = residual.base_dimension, residual.residual_dimension
    out = zeros(Float64, nG + nR, nG + nR)
    out[1:nG, 1:nG] .= O_GG
    out[1:nG, (nG + 1):(nG + nR)] .= O_GR
    out[(nG + 1):(nG + nR), 1:nG] .= transpose(O_GR)
    out[(nG + 1):(nG + nR), (nG + 1):(nG + nR)] .= O_RR
    return out
end

function _r3a_product_matrix(basis, ax, ay, az)
    matrix = zeros(Float64, basis.final_dimension, basis.final_dimension)
    assemble_terminal_product_operator!(matrix, basis, ax, ay, az)
    return matrix
end

function _r3a_centered_factor_terms(axis, expansion, center)
    center == axis.center && return axis.gaussian_factor_terms
    ops = getfield(_GB_PARENT, :mapped_ordinary_one_body_operators)(
        axis.basis; exponents = expansion.exponents, center, backend = axis.backend)
    return ops.gaussian_factors
end

function pqs_terminal_residual_gto_augmented_operators(
    basis::CartesianTerminalBasisRealization,
    bundles,
    parent_basis_object,
    supplement,
    residual::CartesianTerminalResidualGTOAugmentation,
    atom_locations,
    nuclear_charges;
    expansion = nothing,
)
    CPBP = getfield(_GB_PARENT, :CartesianCPBBlockProviders)
    expansion_value = isnothing(expansion) ?
        getfield(_GB_PARENT, :coulomb_gaussian_expansion)(doacc = false) : expansion
    pgdg = Tuple(_nested_axis_pgdg(bundles, axis) for axis in (:x, :y, :z))
    S = Tuple(axis.overlap for axis in pgdg)
    K = _r3a_product_matrix(basis, pgdg[1].kinetic, S[2], S[3]) +
        _r3a_product_matrix(basis, S[1], pgdg[2].kinetic, S[3]) +
        _r3a_product_matrix(basis, S[1], S[2], pgdg[3].kinetic)
    K_GA = _r3a_mixed_block(
        basis, parent_basis_object, supplement, CPBP.cpb_mixed_gto_kinetic_operator_block)
    K_AA = _r3a_dense(CPBP.cpb_gto_kinetic_operator_block(supplement), :gto_kinetic)
    pos_GG = (
        x = _r3a_product_matrix(basis, pgdg[1].position, S[2], S[3]),
        y = _r3a_product_matrix(basis, S[1], pgdg[2].position, S[3]),
        z = _r3a_product_matrix(basis, S[1], S[2], pgdg[3].position),
    )
    x2_GG = (
        x = _r3a_product_matrix(basis, pgdg[1].x2, S[2], S[3]),
        y = _r3a_product_matrix(basis, S[1], pgdg[2].x2, S[3]),
        z = _r3a_product_matrix(basis, S[1], S[2], pgdg[3].x2),
    )
    pos = NamedTuple{(:x, :y, :z)}(Tuple(_r3a_augmented_operator(
        pos_GG[axis],
        _r3a_mixed_block(basis, parent_basis_object, supplement,
            (parent, cpb, supp) -> CPBP.cpb_mixed_gto_position_operator_block(
                parent, cpb, supp; axis)),
        _r3a_dense(CPBP.cpb_gto_position_operator_block(supplement; axis),
            Symbol(:gto_position_, axis)),
        residual) for axis in (:x, :y, :z)))
    x2 = NamedTuple{(:x, :y, :z)}(Tuple(_r3a_augmented_operator(
        x2_GG[axis],
        _r3a_mixed_block(basis, parent_basis_object, supplement,
            (parent, cpb, supp) -> CPBP.cpb_mixed_gto_x2_operator_block(
                parent, cpb, supp; axis)),
        _r3a_dense(CPBP.cpb_gto_x2_operator_block(supplement; axis),
            Symbol(:gto_x2_, axis)),
        residual) for axis in (:x, :y, :z)))
    U = Matrix{Float64}[]
    for (center_index, location) in enumerate(atom_locations)
        center = (Float64(location[1]), Float64(location[2]), Float64(location[3]))
        U_GG = zeros(Float64, basis.final_dimension, basis.final_dimension)
        factors = ntuple(axis -> _r3a_centered_factor_terms(pgdg[axis], expansion_value,
            center[axis]), 3)
        _accumulate_terminal_gaussian_sum!(
            U_GG, basis, expansion_value.coefficients, factors[1], factors[2], factors[3])
        record = (; center_index, nuclear_charge = Float64(nuclear_charges[center_index]),
            location = center)
        U_GA = _r3a_mixed_block(basis, parent_basis_object, supplement,
            (parent, cpb, supp) -> CPBP.cpb_mixed_gto_nuclear_by_center_block(
                parent, cpb, supp, expansion_value, record))
        U_AA = _r3a_dense(CPBP.cpb_gto_nuclear_by_center_block(
            supplement, expansion_value, record), :gto_nuclear_by_center)
        push!(U, _r3a_augmented_operator(U_GG, U_GA, U_AA, residual))
    end
    return (;
        kinetic = _r3a_augmented_operator(K, K_GA, K_AA, residual),
        nuclear_attraction_unit_by_center = U,
        position = pos,
        x2,
    )
end
