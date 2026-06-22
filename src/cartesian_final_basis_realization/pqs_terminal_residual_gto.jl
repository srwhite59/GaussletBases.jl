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

_r3_center(center) = (Float64(center[1]), Float64(center[2]), Float64(center[3]))
_r3_float_centers(centers) = NTuple{3,Float64}[_r3_center(center) for center in centers]
_r3_candidate_labels(supplement) = String[String(orbital.label) for orbital in supplement.orbitals]
_r3_candidate_centers(supplement) = NTuple{3,Float64}[_r3_center(orbital.center) for orbital in supplement.orbitals]
_r3_require_size(matrix, dims, label) = size(matrix) == dims || throw(DimensionMismatch(label))
_r3_require_close(block, reference, label) = norm(block - reference, Inf) <= 1.0e-10 || throw(ArgumentError(label))
_r3_moment_ok(matrix, dims) = size(matrix) == dims && all(isfinite, matrix) && norm(matrix - transpose(matrix), Inf) <= 1.0e-10
function _r3_validate_residual_contract(
    basis,
    supplement,
    residual::CartesianTerminalResidualGTOAugmentation,
    atom_locations,
)
    nG, nR = basis.final_dimension, residual.residual_dimension
    residual.base_dimension == basis.final_dimension ||
        throw(DimensionMismatch("R3 residual base dimension must match terminal basis"))
    residual.residual_dimension == length(residual.residual_labels) ||
        throw(DimensionMismatch("R3 residual label count mismatch"))
    _r3_require_size(residual.T_G, (nG, nR), "R3 residual T_G shape mismatch")
    _r3_require_size(residual.T_A, (residual.candidate_count, nR), "R3 residual T_A shape mismatch")
    if !isnothing(supplement)
        residual.candidate_count == length(supplement.orbitals) ||
            throw(DimensionMismatch("R3 residual candidate count mismatch"))
        residual.candidate_labels == _r3_candidate_labels(supplement) ||
            throw(ArgumentError("R3 residual candidate labels do not match supplement"))
        residual.candidate_centers == _r3_candidate_centers(supplement) ||
            throw(ArgumentError("R3 residual candidate centers do not match supplement"))
    end
    locations = _r3_float_centers(atom_locations)
    for (index, owner) in pairs(residual.candidate_owner_indices)
        1 <= owner <= length(locations) ||
            throw(ArgumentError("R3 residual candidate owner index is out of range"))
        locations[owner] == residual.candidate_centers[index] ||
            throw(ArgumentError("R3 residual candidate owner center mismatch"))
    end
    return nothing
end
function _r3_validate_augmented_operator_dimensions(operators, base_hamiltonian, residual, center_count)
    nG, n = residual.base_dimension, residual.base_dimension + residual.residual_dimension
    length(operators.nuclear_attraction_unit_by_center) == center_count || throw(DimensionMismatch("R3 augmented unit nuclear center count mismatch"))
    for matrix in (operators.kinetic, operators.nuclear_attraction_unit_by_center...)
        _r3_require_size(matrix, (n, n), "R3 augmented operator dimension mismatch")
    end
    _r3_require_close(view(operators.kinetic, 1:nG, 1:nG), base_hamiltonian.kinetic, "R3 augmented kinetic G-G block mismatch")
    for (matrix, base) in zip(operators.nuclear_attraction_unit_by_center,
                              base_hamiltonian.nuclear_attraction_unit_by_center)
        _r3_require_close(view(matrix, 1:nG, 1:nG), base, "R3 augmented unit nuclear G-G block mismatch")
    end
    for matrix in (operators.position.x, operators.position.y, operators.position.z,
                   operators.x2.x, operators.x2.y, operators.x2.z)
        _r3_moment_ok(matrix, (n, n)) || throw(ArgumentError("R3 augmented moment matrix invalid"))
    end
    return nothing
end
function _r3_validate_base_hamiltonian(base_hamiltonian, residual)
    nG, center_count = residual.base_dimension, length(base_hamiltonian.nuclear_charges)
    _r3_require_size(base_hamiltonian.nuclear_positions, (center_count, 3), "R3-B base Hamiltonian center metadata mismatch")
    length(base_hamiltonian.nuclear_attraction_unit_by_center) == center_count || throw(DimensionMismatch("R3-B base Hamiltonian unit nuclear count mismatch"))
    for matrix in (base_hamiltonian.electron_electron_ida, base_hamiltonian.kinetic, base_hamiltonian.nuclear_attraction_unit_by_center...)
        _r3_require_size(matrix, (nG, nG), "R3-B base Hamiltonian matrix dimension mismatch")
    end
    return center_count
end

_r3b_default_expansion(expansion) = isnothing(expansion) ? getfield(_GB_PARENT, :coulomb_gaussian_expansion)(doacc = false) : throw(ArgumentError("R3-B augmented Hamiltonian rejects custom expansion because base V_GG has no expansion provenance"))

function _r3b_residual_mwg_descriptors(operators, residual)
    nG, nR = residual.base_dimension, residual.residual_dimension
    rows = (nG + 1):(nG + nR)
    centers = zeros(Float64, nR, 3)
    widths = zeros(Float64, nR, 3)
    for (axis_index, axis) in pairs((:x, :y, :z))
        position = operators.position[axis]
        second = operators.x2[axis]
        for (residual_index, row) in enumerate(rows)
            center = Float64(position[row, row])
            variance = Float64(second[row, row]) - center^2
            isfinite(center) && isfinite(variance) && variance > 0.0 ||
                throw(ArgumentError("R3-B residual MWG moment variance must be finite positive"))
            centers[residual_index, axis_index] = center
            widths[residual_index, axis_index] = sqrt(2.0 * variance)
            isfinite(widths[residual_index, axis_index]) &&
                widths[residual_index, axis_index] > 0.0 ||
                throw(ArgumentError("R3-B residual MWG width must be finite positive"))
        end
    end
    return centers, widths
end

_r3b_axis_bundle(bundles, axis::Int) =
    axis == 1 ? bundles.bundle_x : axis == 2 ? bundles.bundle_y : bundles.bundle_z

function _r3b_mwg_axis_pairs(bundles, expansion, residual_centers, residual_widths)
    qw = _GB_PARENT
    effective = getfield(qw, :_qwrg_effective_gaussians)(residual_centers, residual_widths)
    split = ntuple(axis -> getfield(qw, :_qwrg_split_block_matrices)(
        _r3b_axis_bundle(bundles, axis), effective[axis], expansion), 3)
    analytic = ntuple(axis -> getfield(qw, :_qwrg_gaussian_analytic_blocks)(
        effective[axis], expansion), 3)
    weights = ntuple(axis -> getfield(qw, :_qwrg_supplement_integral_weights)(
        effective[axis]), 3)
    normalize = getfield(qw, :_qwrg_density_normalized_pair_matrices)
    return (;
        ga = ntuple(axis -> normalize(split[axis].pair_ga, split[axis].weight_gg,
            weights[axis]; label = "R3-B MWG $((:x, :y, :z)[axis]) G-M"), 3),
        aa = ntuple(axis -> normalize(analytic[axis].pair_aa, weights[axis],
            weights[axis]; label = "R3-B MWG $((:x, :y, :z)[axis]) M-M"), 3),
    )
end

function _r3b_terminal_mwg_fixed_residual(basis, bundles, pair_terms, coefficients)
    nR = size(first(pair_terms.ga[1]), 2)
    V_GM = zeros(Float64, basis.final_dimension, nR)
    for block in basis.blocks
        local_values = zeros(Float64, length(block.support_states), nR)
        for residual in 1:nR
            column = view(local_values, :, residual)
            @inbounds for (row, state) in pairs(block.support_states)
                ix, iy, iz = state
                value = 0.0
                for term in eachindex(coefficients)
                    value += Float64(coefficients[term]) *
                        pair_terms.ga[1][term][ix, residual] *
                        pair_terms.ga[2][term][iy, residual] *
                        pair_terms.ga[3][term][iz, residual]
                end
                column[row] = value
            end
        end
        block.coefficients === nothing ?
            (V_GM[block.column_range, :] .= local_values) :
            begin
                support_weights = _support_weights(block.support_states, bundles)
                final_weights = vec(transpose(block.coefficients) * support_weights)
                all(weight -> isfinite(weight) && weight > 1.0e-12, final_weights) ||
                    throw(ArgumentError("R3-B final density weights must be finite positive"))
                density_coefficients =
                    block.coefficients .* reshape(support_weights, :, 1) ./
                    reshape(final_weights, 1, :)
                mul!(view(V_GM, block.column_range, :), transpose(density_coefficients),
                    local_values)
            end
    end
    return V_GM
end

function _r3b_mwg_residual_residual(pair_terms, coefficients)
    nR = size(first(pair_terms.aa[1]), 1)
    V_MM = zeros(Float64, nR, nR)
    for i in 1:nR, j in i:nR
        value = 0.0
        for term in eachindex(coefficients)
            value += Float64(coefficients[term]) *
                pair_terms.aa[1][term][i, j] *
                pair_terms.aa[2][term][i, j] *
                pair_terms.aa[3][term][i, j]
        end
        V_MM[i, j] = value
        V_MM[j, i] = value
    end
    return V_MM
end

function pqs_terminal_residual_gto_augmented_hamiltonian(
    base_hamiltonian,
    basis::CartesianTerminalBasisRealization,
    bundles,
    residual::CartesianTerminalResidualGTOAugmentation,
    augmented_operators;
    expansion = nothing,
)
    center_count = _r3_validate_base_hamiltonian(base_hamiltonian, residual)
    atom_locations = NTuple{3,Float64}[_r3_center(view(base_hamiltonian.nuclear_positions, index, :)) for index in 1:center_count]
    _r3_validate_residual_contract(basis, nothing, residual, atom_locations)
    _r3_validate_augmented_operator_dimensions(augmented_operators, base_hamiltonian, residual, center_count)
    centers, widths = _r3b_residual_mwg_descriptors(augmented_operators, residual)
    expansion_value = _r3b_default_expansion(expansion)
    pair_terms = _r3b_mwg_axis_pairs(bundles, expansion_value, centers, widths)
    V_GM = _r3b_terminal_mwg_fixed_residual(
        basis, bundles, pair_terms, expansion_value.coefficients)
    V_MM = _r3b_mwg_residual_residual(pair_terms, expansion_value.coefficients)
    nG, nR = residual.base_dimension, residual.residual_dimension
    V = zeros(Float64, nG + nR, nG + nR)
    residual_range = (nG + 1):(nG + nR)
    V[1:nG, 1:nG] .= base_hamiltonian.electron_electron_ida
    V[1:nG, residual_range] .= V_GM
    V[residual_range, 1:nG] .= transpose(V_GM)
    V[residual_range, residual_range] .= V_MM
    Hamiltonian = getfield(_GB_PARENT, :CartesianIDAHamiltonian)
    return Hamiltonian(
        augmented_operators.kinetic,
        augmented_operators.nuclear_attraction_unit_by_center,
        V,
        base_hamiltonian.nup,
        base_hamiltonian.ndn;
        nuclear_charges = base_hamiltonian.nuclear_charges,
        nuclear_positions = base_hamiltonian.nuclear_positions,
    )
end

function _residual_candidate_owner(center, nuclei)
    matches = findall(==(center), nuclei)
    length(matches) == 1 ||
        throw(ArgumentError("residual-GTO candidate center must exactly match one nucleus"))
    return first(matches)
end

function _r3a_qw_orbital(orbital)
    orbital.primitive_normalization === :axiswise_normalized_cartesian_gaussian ||
        throw(ArgumentError("R3-A QW donor requires axiswise-normalized Cartesian Gaussian primitives"))
    ctor = getfield(_GB_PARENT, :_AtomicCartesianShellOrbital3D)
    lx, ly, lz = orbital.angular_powers
    return ctor(orbital.label, lx, ly, lz, orbital.exponents, orbital.coefficients,
        orbital.center)
end

_r3a_qw_supplement(supplement) =
    (; orbitals = [_r3a_qw_orbital(orbital) for orbital in supplement.orbitals])

function _r3a_qw_proxy_layers(bundles)
    pgdg = (; x = _nested_axis_pgdg(bundles, :x),
        y = _nested_axis_pgdg(bundles, :y),
        z = _nested_axis_pgdg(bundles, :z))
    proxy = getfield(_GB_PARENT, :_qwrg_diatomic_supplement_proxy_layer)
    return (;
        x = proxy(pgdg.x.basis, bundles.bundle_x, :x),
        y = proxy(pgdg.y.basis, bundles.bundle_y, :y),
        z = proxy(pgdg.z.basis, bundles.bundle_z, :z),
        ncart = size(pgdg.x.overlap, 1) * size(pgdg.y.overlap, 1) *
                size(pgdg.z.overlap, 1),
    )
end

_r3a_sym(matrix) = Matrix{Float64}(0.5 .* (matrix .+ transpose(matrix)))

function _r3a_project_parent_ga(basis, parent_ga)
    out = zeros(Float64, basis.final_dimension, size(parent_ga, 2))
    for block in basis.blocks
        rows = view(parent_ga, block.support_indices, :)
        target = view(out, block.column_range, :)
        block.coefficients === nothing ?
            (target .= rows) :
            mul!(target, transpose(block.coefficients), rows)
    end
    return out
end

function _r3a_qw_nuclear_blocks(proxy, supplement, expansion, atom_locations)
    qw = _GB_PARENT
    ncart, norbital = proxy.ncart, length(supplement.orbitals)
    ga = [zeros(Float64, ncart, norbital) for _ in atom_locations]
    aa = [zeros(Float64, norbital, norbital) for _ in atom_locations]
    cross_cache = Dict()
    scratch = zeros(Float64, ncart)
    fill_product! = getfield(qw, :_qwrg_fill_product_column!)
    axis_cross = getfield(qw, :_qwrg_atomic_axis_factor_cross_data)
    for (orbital_index, orbital) in pairs(supplement.orbitals)
        for (center_index, center) in pairs(atom_locations)
            factors = (
                get!(cross_cache, (orbital_index, :x, center[1])) do
                    axis_cross(proxy.x, orbital, :x, expansion, center[1])
                end,
                get!(cross_cache, (orbital_index, :y, center[2])) do
                    axis_cross(proxy.y, orbital, :y, expansion, center[2])
                end,
                get!(cross_cache, (orbital_index, :z, center[3])) do
                    axis_cross(proxy.z, orbital, :z, expansion, center[3])
                end,
            )
            for term in eachindex(expansion.coefficients),
                    primitive in eachindex(orbital.coefficients)
                fill_product!(scratch,
                    view(factors[1][term], :, primitive),
                    view(factors[2][term], :, primitive),
                    view(factors[3][term], :, primitive))
                ga[center_index][:, orbital_index] .-=
                    expansion.coefficients[term] *
                    Float64(orbital.coefficients[primitive]) .* scratch
            end
        end
    end
    axis_aa = getfield(qw, :_qwrg_atomic_axis_factor_aa_data)
    weighted = getfield(qw, :_qwrg_atomic_weighted_hadamard)
    for (left_index, left) in pairs(supplement.orbitals),
            (right_index, right) in pairs(supplement.orbitals)
        local_cache = Dict()
        for (center_index, center) in pairs(atom_locations)
            fx = get!(local_cache, (:x, center[1])) do
                axis_aa(left, right, :x, expansion, center[1])
            end
            fy = get!(local_cache, (:y, center[2])) do
                axis_aa(left, right, :y, expansion, center[2])
            end
            fz = get!(local_cache, (:z, center[3])) do
                axis_aa(left, right, :z, expansion, center[3])
            end
            value = 0.0
            for term in eachindex(expansion.coefficients)
                value -= expansion.coefficients[term] *
                    weighted(left.coefficients, fx[term], fy[term], fz[term],
                        right.coefficients)
            end
            aa[center_index][left_index, right_index] = value
        end
    end
    return (; ga, aa = [_r3a_sym(matrix) for matrix in aa])
end

function _r3a_qw_blocks(basis, bundles, supplement, atom_locations, expansion)
    donor = _r3a_qw_supplement(supplement)
    proxy = _r3a_qw_proxy_layers(bundles)
    cross = getfield(_GB_PARENT, :_qwrg_cartesian_shell_cross_moment_blocks_3d)(
        (x = proxy.x, y = proxy.y, z = proxy.z), donor, expansion, proxy.ncart;
        include_factor_terms = false)
    self = getfield(_GB_PARENT, :_qwrg_cartesian_shell_self_moment_blocks_3d)(
        donor, expansion; include_factor_terms = false)
    nuclear = _r3a_qw_nuclear_blocks(proxy, donor, expansion, _r3_float_centers(atom_locations))
    return (;
        mixed = (;
            overlap = _r3a_project_parent_ga(basis, cross.overlap_ga),
            kinetic = _r3a_project_parent_ga(basis, cross.kinetic_ga),
            position = (x = _r3a_project_parent_ga(basis, cross.position_x_ga),
                y = _r3a_project_parent_ga(basis, cross.position_y_ga),
                z = _r3a_project_parent_ga(basis, cross.position_z_ga)),
            x2 = (x = _r3a_project_parent_ga(basis, cross.x2_x_ga),
                y = _r3a_project_parent_ga(basis, cross.x2_y_ga),
                z = _r3a_project_parent_ga(basis, cross.x2_z_ga)),
            nuclear = [_r3a_project_parent_ga(basis, matrix) for matrix in nuclear.ga],
        ),
        self = (;
            overlap = _r3a_sym(self.overlap_aa),
            kinetic = _r3a_sym(self.kinetic_aa),
            position = (x = _r3a_sym(self.position_x_aa),
                y = _r3a_sym(self.position_y_aa), z = _r3a_sym(self.position_z_aa)),
            x2 = (x = _r3a_sym(self.x2_x_aa),
                y = _r3a_sym(self.x2_y_aa), z = _r3a_sym(self.x2_z_aa)),
            nuclear = nuclear.aa,
        ),
    )
end

function _terminal_residual_mixed_overlap(
    basis::CartesianTerminalBasisRealization,
    bundles,
    supplement,
)
    expansion = getfield(_GB_PARENT, :coulomb_gaussian_expansion)(doacc = false)
    return _r3a_qw_blocks(basis, bundles, supplement, (), expansion).mixed.overlap
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
    nuclei_value = _r3_float_centers(nuclei)
    labels = _r3_candidate_labels(supplement)
    centers = _r3_candidate_centers(supplement)
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
    center == axis.center && Float64.(axis.exponents) == Float64.(expansion.exponents) &&
        return axis.gaussian_factor_terms
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
    length(atom_locations) == length(nuclear_charges) || throw(DimensionMismatch("R3-A atom location count must match nuclear charges"))
    _r3_validate_residual_contract(basis, supplement, residual, atom_locations)
    expansion_value = isnothing(expansion) ?
        getfield(_GB_PARENT, :coulomb_gaussian_expansion)(doacc = false) : expansion
    supplement_blocks = _r3a_qw_blocks(basis, bundles, supplement, atom_locations,
        expansion_value)
    pgdg = Tuple(_nested_axis_pgdg(bundles, axis) for axis in (:x, :y, :z))
    S = Tuple(axis.overlap for axis in pgdg)
    K = _r3a_product_matrix(basis, pgdg[1].kinetic, S[2], S[3]) +
        _r3a_product_matrix(basis, S[1], pgdg[2].kinetic, S[3]) +
        _r3a_product_matrix(basis, S[1], S[2], pgdg[3].kinetic)
    K_GA = supplement_blocks.mixed.kinetic
    K_AA = supplement_blocks.self.kinetic
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
        supplement_blocks.mixed.position[axis],
        supplement_blocks.self.position[axis],
        residual) for axis in (:x, :y, :z)))
    x2 = NamedTuple{(:x, :y, :z)}(Tuple(_r3a_augmented_operator(
        x2_GG[axis],
        supplement_blocks.mixed.x2[axis],
        supplement_blocks.self.x2[axis],
        residual) for axis in (:x, :y, :z)))
    U = Matrix{Float64}[]
    for (center_index, center) in enumerate(_r3_float_centers(atom_locations))
        U_GG = zeros(Float64, basis.final_dimension, basis.final_dimension)
        factors = ntuple(axis -> _r3a_centered_factor_terms(pgdg[axis], expansion_value,
            center[axis]), 3)
        _accumulate_terminal_gaussian_sum!(
            U_GG, basis, expansion_value.coefficients, factors[1], factors[2], factors[3])
        U_GA = supplement_blocks.mixed.nuclear[center_index]
        U_AA = supplement_blocks.self.nuclear[center_index]
        push!(U, _r3a_augmented_operator(U_GG, U_GA, U_AA, residual))
    end
    return (;
        kinetic = _r3a_augmented_operator(K, K_GA, K_AA, residual),
        nuclear_attraction_unit_by_center = U,
        position = pos,
        x2,
    )
end
