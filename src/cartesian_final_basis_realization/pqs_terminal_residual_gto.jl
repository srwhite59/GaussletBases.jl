const _GB_PARENT = parentmodule(@__MODULE__)
const CartesianTerminalResidualGTOAugmentation = CRG.CartesianResidualGaussianBasis

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
    length(residual.candidate_labels) == residual.candidate_count &&
        length(residual.candidate_owner_indices) == residual.candidate_count &&
        length(residual.candidate_centers) == residual.candidate_count ||
        throw(DimensionMismatch("R3 residual candidate metadata count mismatch"))
    residual.residual_dimension == length(residual.residual_labels) ||
        throw(DimensionMismatch("R3 residual label count mismatch"))
    length(residual.residual_source_owner_indices) == nR &&
        length(residual.residual_occupations) == nR &&
        sum(residual.owner_retained_counts) == nR ||
        throw(DimensionMismatch("R3 residual owner-local metadata count mismatch"))
    all(>(residual.occupation_cutoff), residual.residual_occupations) ||
        throw(ArgumentError("R3 retained residual occupations must exceed the cutoff"))
    _r3_require_size(residual.T_G, (nG, nR), "R3 residual T_G shape mismatch")
    _r3_require_size(residual.T_A, (residual.candidate_count, nR), "R3 residual T_A shape mismatch")
    if !isnothing(supplement)
        residual.candidate_count == length(supplement.orbitals) ||
            throw(DimensionMismatch("R3 residual candidate count mismatch"))
        residual.candidate_labels == CRG.residual_gaussian_candidate_labels(supplement) ||
            throw(ArgumentError("R3 residual candidate labels do not match supplement"))
        residual.candidate_centers == CRG.residual_gaussian_candidate_centers(supplement) ||
            throw(ArgumentError("R3 residual candidate centers do not match supplement"))
    end
    locations = CRG.residual_gaussian_float_centers(atom_locations)
    for (index, owner) in pairs(residual.candidate_owner_indices)
        1 <= owner <= length(locations) ||
            throw(ArgumentError("R3 residual candidate owner index is out of range"))
        locations[owner] == residual.candidate_centers[index] ||
            throw(ArgumentError("R3 residual candidate owner center mismatch"))
    end
    for owner in residual.residual_source_owner_indices
        1 <= owner <= length(locations) ||
            throw(ArgumentError("R3 residual source owner index is out of range"))
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

_r3b_residual_mwg_descriptors(operators, residual) =
    CRG.moment_matched_gaussians(operators, residual)
_r3b_mwg_axis_pairs(bundles, expansion, centers, widths) =
    CRG._mwg_axis_pairs(bundles, expansion, centers, widths)

function pqs_terminal_residual_gto_augmented_hamiltonian(
    base_hamiltonian,
    basis::CartesianTerminalBasisRealization,
    bundles,
    residual::CartesianTerminalResidualGTOAugmentation,
    augmented_operators;
    expansion = nothing,
)
    center_count = _r3_validate_base_hamiltonian(base_hamiltonian, residual)
    atom_locations = NTuple{3,Float64}[CRG.residual_gaussian_center(
        view(base_hamiltonian.nuclear_positions, index, :)) for index in 1:center_count]
    _r3_validate_residual_contract(basis, nothing, residual, atom_locations)
    _r3_validate_augmented_operator_dimensions(augmented_operators, base_hamiltonian, residual, center_count)
    V = CRG.assemble_residual_ida_interaction(
        base_hamiltonian.electron_electron_ida, basis, bundles, residual,
        augmented_operators; expansion)
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

function pqs_terminal_residual_gto_augmented_hamiltonian(
    base_hamiltonian,
    basis::CartesianTerminalBasisRealization,
    bundles,
    supplement,
    atom_locations,
    nuclear_charges;
    expansion = nothing,
)
    residual = pqs_terminal_residual_gto_augmentation(
        basis, bundles, supplement, atom_locations)
    augmented_operators = pqs_terminal_residual_gto_augmented_operators(
        basis, bundles, nothing, supplement, residual, atom_locations, nuclear_charges;
        expansion)
    return pqs_terminal_residual_gto_augmented_hamiltonian(
        base_hamiltonian, basis, bundles, residual, augmented_operators; expansion)
end

function _r3_supplement_owner_counts(residual, center_count)
    counts = zeros(Int, center_count)
    length(residual.candidate_owner_indices) == residual.candidate_count ||
        throw(DimensionMismatch("R3 supplement owner count mismatch"))
    for owner in residual.candidate_owner_indices
        1 <= owner <= center_count ||
            throw(ArgumentError("R3 supplement owner index out of range"))
        counts[owner] += 1
    end
    return counts
end

function write_pqs_terminal_residual_gto_augmented_hamiltonian(
    path,
    hamiltonian,
    residual::CartesianTerminalResidualGTOAugmentation;
    basis_by_center,
    lmax::Integer,
    uncontracted::Bool,
    width_filtering = nothing,
    validation_check_labels = (),
    h2_validation_self_coulomb = nothing,
)
    Hamiltonian = getfield(_GB_PARENT, :CartesianIDAHamiltonian)
    hamiltonian isa Hamiltonian{Float64} ||
        throw(ArgumentError("R3 augmented artifact writer requires CartesianIDAHamiltonian{Float64}"))
    augmented_dimension = size(hamiltonian.kinetic, 1)
    augmented_dimension == residual.base_dimension + residual.residual_dimension ||
        throw(DimensionMismatch("R3 augmented artifact dimension mismatch"))
    basis_labels = String[String(label) for label in basis_by_center]
    owner_counts = _r3_supplement_owner_counts(residual, length(basis_labels))
    values = (;
        provenance_version = 1,
        producer = :cartesian_residual_gto_mwg_augmentation,
        supplement_policy = :mwg_residual_gto,
        basis_by_center = basis_labels,
        lmax = Int(lmax),
        uncontracted = Bool(uncontracted),
        width_filtering,
        candidate_count = residual.candidate_count,
        owner_counts,
        base_dimension = residual.base_dimension,
        residual_dimension = residual.residual_dimension,
        augmented_dimension,
        augmented_basis_order = :base_then_residual,
        residual_basis_convention = residual.orientation,
        rank_rule = residual.selection_rule,
        occupation_cutoff = residual.occupation_cutoff,
        tau_neg_abs = residual.tau_neg_abs,
        tau_neg_rel = residual.tau_neg_rel,
        tau_merge_abs = residual.tau_merge_abs,
        tau_merge_rel = residual.tau_merge_rel,
        mwg_convention_version = 1,
        mwg_convention = :separable_moment_matched_density_normalized,
        one_body_source = :exact_transformed_raw_blocks,
        interaction_source = :weight_aware_residual_mwg_ida_blocks,
        validation_check_labels = Symbol[Symbol(label) for label in validation_check_labels],
        h2_self_coulomb_reference = isnothing(h2_validation_self_coulomb) ?
            nothing : Float64(h2_validation_self_coulomb),
    )
    getfield(_GB_PARENT, :write_cartesian_ida_hamiltonian)(String(path), hamiltonian)
    getfield(_GB_PARENT, :jldopen)(String(path), "r+") do file
        for key in keys(values)
            file["supplement_provenance/$(key)"] = getproperty(values, key)
        end
    end
    return path
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
    return (; ga, aa = [CRG.symmetrize_operator(matrix) for matrix in aa])
end

function _r3a_qw_blocks(basis, bundles, supplement, atom_locations, expansion)
    donor = _r3a_qw_supplement(supplement)
    proxy = _r3a_qw_proxy_layers(bundles)
    cross = getfield(_GB_PARENT, :_qwrg_cartesian_shell_cross_moment_blocks_3d)(
        (x = proxy.x, y = proxy.y, z = proxy.z), donor, expansion, proxy.ncart;
        include_factor_terms = false)
    self = getfield(_GB_PARENT, :_qwrg_cartesian_shell_self_moment_blocks_3d)(
        donor, expansion; include_factor_terms = false)
    nuclear = _r3a_qw_nuclear_blocks(proxy, donor, expansion,
        CRG.residual_gaussian_float_centers(atom_locations))
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
            overlap = CRG.symmetrize_operator(self.overlap_aa),
            kinetic = CRG.symmetrize_operator(self.kinetic_aa),
            position = (x = CRG.symmetrize_operator(self.position_x_aa),
                y = CRG.symmetrize_operator(self.position_y_aa),
                z = CRG.symmetrize_operator(self.position_z_aa)),
            x2 = (x = CRG.symmetrize_operator(self.x2_x_aa),
                y = CRG.symmetrize_operator(self.x2_y_aa),
                z = CRG.symmetrize_operator(self.x2_z_aa)),
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
    return CRG.terminal_residual_mixed_overlap(basis, bundles, supplement,
        _r3a_qw_blocks, expansion)
end

function pqs_terminal_residual_gto_augmentation(
    basis::CartesianTerminalBasisRealization,
    bundles,
    supplement,
    nuclei;
    residual_occupation_cutoff::Real = 1.0e-8,
    tau_neg_abs::Real = 1.0e-12,
    tau_neg_rel::Real = 1.0e-12,
    tau_merge_abs::Real = 1.0e-12,
    tau_merge_rel::Real = 1.0e-12,
    orthogonality_atol::Real = 1.0e-10,
    identity_atol::Real = 1.0e-10,
)
    nuclei_value = CRG.residual_gaussian_float_centers(nuclei)
    labels = CRG.residual_gaussian_candidate_labels(supplement)
    centers = CRG.residual_gaussian_candidate_centers(supplement)
    owners = Int[CRG.residual_candidate_owner(center, nuclei_value) for center in centers]
    X = _terminal_residual_mixed_overlap(basis, bundles, supplement)
    S_AA = Matrix{Float64}(
        getfield(_GB_PARENT, :_cartesian_supplement_cross_overlap)(supplement, supplement))
    return CRG.build_residual_gaussian_basis(basis.final_dimension, X, S_AA,
        labels, centers, owners; residual_occupation_cutoff, tau_neg_abs,
        tau_neg_rel, tau_merge_abs, tau_merge_rel, orthogonality_atol, identity_atol)
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
    pos = NamedTuple{(:x, :y, :z)}(Tuple(CRG.transform_augmented_operator(
        pos_GG[axis],
        supplement_blocks.mixed.position[axis],
        supplement_blocks.self.position[axis],
        residual) for axis in (:x, :y, :z)))
    x2 = NamedTuple{(:x, :y, :z)}(Tuple(CRG.transform_augmented_operator(
        x2_GG[axis],
        supplement_blocks.mixed.x2[axis],
        supplement_blocks.self.x2[axis],
        residual) for axis in (:x, :y, :z)))
    U = Matrix{Float64}[]
    for (center_index, center) in enumerate(CRG.residual_gaussian_float_centers(atom_locations))
        U_GG = zeros(Float64, basis.final_dimension, basis.final_dimension)
        factors = ntuple(axis -> _r3a_centered_factor_terms(pgdg[axis], expansion_value,
            center[axis]), 3)
        _accumulate_terminal_gaussian_sum!(
            U_GG, basis, expansion_value.coefficients, factors[1], factors[2], factors[3])
        U_GA = supplement_blocks.mixed.nuclear[center_index]
        U_AA = supplement_blocks.self.nuclear[center_index]
        push!(U, CRG.transform_augmented_operator(U_GG, U_GA, U_AA, residual))
    end
    return (;
        kinetic = CRG.transform_augmented_operator(K, K_GA, K_AA, residual),
        nuclear_attraction_unit_by_center = U,
        position = pos,
        x2,
    )
end
