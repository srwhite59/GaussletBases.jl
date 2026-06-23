const _GB_PARENT = parentmodule(@__MODULE__)
const CGRB = getfield(_GB_PARENT, :CartesianGaussianRawBlocks)
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

function pqs_terminal_residual_gto_augmented_hamiltonian(
    base_hamiltonian,
    basis::CartesianTerminalBasisRealization,
    bundles,
    residual::CartesianTerminalResidualGTOAugmentation,
    augmented_operators;
    expansion = nothing,
)
    V = pqs_terminal_residual_gto_augmented_vee(
        base_hamiltonian, basis, bundles, residual, augmented_operators;
        expansion)
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

function pqs_terminal_residual_gto_augmented_vee(base_hamiltonian, basis::CartesianTerminalBasisRealization, bundles, residual::CartesianTerminalResidualGTOAugmentation, augmented_operators; expansion = nothing)
    center_count = _r3_validate_base_hamiltonian(base_hamiltonian, residual)
    atom_locations = NTuple{3,Float64}[CRG.residual_gaussian_center(
        view(base_hamiltonian.nuclear_positions, index, :)) for index in 1:center_count]
    _r3_validate_residual_contract(basis, nothing, residual, atom_locations)
    _r3_validate_augmented_operator_dimensions(augmented_operators, base_hamiltonian, residual, center_count)
    return CRG.assemble_residual_ida_interaction(
        base_hamiltonian.electron_electron_ida, basis, bundles, residual,
        augmented_operators; expansion)
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
    expansion_value = isnothing(expansion) ?
        getfield(_GB_PARENT, :coulomb_gaussian_expansion)(doacc = false) : expansion
    supplement_blocks = _r3a_qw_blocks(basis, bundles, supplement, atom_locations,
        expansion_value)
    nuclei_value = CRG.residual_gaussian_float_centers(atom_locations)
    labels = CRG.residual_gaussian_candidate_labels(supplement)
    centers = CRG.residual_gaussian_candidate_centers(supplement)
    owners = Int[CRG.residual_candidate_owner(center, nuclei_value) for center in centers]
    S_AA = Matrix{Float64}(
        getfield(_GB_PARENT, :_cartesian_supplement_cross_overlap)(supplement, supplement))
    residual = CRG.build_residual_gaussian_basis(
        basis.final_dimension, supplement_blocks.mixed.overlap, S_AA, labels, centers, owners)
    augmented_operators = pqs_terminal_residual_gto_augmented_operators(
        basis, bundles, nothing, supplement, residual, atom_locations, nuclear_charges;
        expansion = expansion_value, supplement_blocks,
        base_kinetic = base_hamiltonian.kinetic,
        base_unit_nuclear = base_hamiltonian.nuclear_attraction_unit_by_center)
    return pqs_terminal_residual_gto_augmented_hamiltonian(
        base_hamiltonian, basis, bundles, residual, augmented_operators;
        expansion)
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

function _r3a_qw_blocks(basis, bundles, supplement, atom_locations, expansion)
    donor = _r3a_qw_supplement(supplement)
    proxy = _r3a_qw_proxy_layers(bundles)
    non_nuclear = CGRB.gaussian_non_nuclear_raw_blocks(proxy, donor, expansion)
    nuclear = CGRB.gaussian_nuclear_raw_blocks_by_center(proxy, donor, expansion,
        CRG.residual_gaussian_float_centers(atom_locations))
    return (;
        mixed = (;
            overlap = _r3a_project_parent_ga(basis, non_nuclear.ga.overlap),
            kinetic = _r3a_project_parent_ga(basis, non_nuclear.ga.kinetic),
            position = (x = _r3a_project_parent_ga(basis, non_nuclear.ga.position.x),
                y = _r3a_project_parent_ga(basis, non_nuclear.ga.position.y),
                z = _r3a_project_parent_ga(basis, non_nuclear.ga.position.z)),
            x2 = (x = _r3a_project_parent_ga(basis, non_nuclear.ga.x2.x),
                y = _r3a_project_parent_ga(basis, non_nuclear.ga.x2.y),
                z = _r3a_project_parent_ga(basis, non_nuclear.ga.x2.z)),
            nuclear = [_r3a_project_parent_ga(basis, matrix) for matrix in nuclear.ga],
        ),
        self = (;
            overlap = non_nuclear.aa.overlap,
            kinetic = non_nuclear.aa.kinetic,
            position = non_nuclear.aa.position,
            x2 = non_nuclear.aa.x2,
            nuclear = nuclear.aa,
        ),
    )
end

function _terminal_residual_mixed_overlap(
    basis::CartesianTerminalBasisRealization,
    bundles,
    supplement,
)
    donor = _r3a_qw_supplement(supplement)
    proxy = _r3a_qw_proxy_layers(bundles)
    blocks = CGRB.gaussian_non_nuclear_overlap_blocks(proxy, donor)
    return _r3a_project_parent_ga(basis, blocks.ga.overlap)
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

function _r3a_centered_factor_terms(axis, expansion, center)
    center == axis.center && Float64.(axis.exponents) == Float64.(expansion.exponents) &&
        return axis.gaussian_factor_terms
    ops = getfield(_GB_PARENT, :mapped_ordinary_one_body_operators)(
        axis.basis; exponents = expansion.exponents, center, backend = axis.backend)
    return ops.gaussian_factors
end

function pqs_terminal_residual_gto_augmented_products(basis::CartesianTerminalBasisRealization, bundles, parent_basis_object, supplement, residual::CartesianTerminalResidualGTOAugmentation, atom_locations, nuclear_charges; expansion = nothing, supplement_blocks = nothing, base_kinetic = nothing)
    length(atom_locations) == length(nuclear_charges) || throw(DimensionMismatch("R3-A atom location count must match nuclear charges"))
    _r3_validate_residual_contract(basis, supplement, residual, atom_locations)
    expansion_value = isnothing(expansion) ?
        getfield(_GB_PARENT, :coulomb_gaussian_expansion)(doacc = false) : expansion
    supplement_blocks_value = isnothing(supplement_blocks) ?
        _r3a_qw_blocks(basis, bundles, supplement, atom_locations, expansion_value) :
        supplement_blocks
    pgdg = Tuple(_nested_axis_pgdg(bundles, axis) for axis in (:x, :y, :z))
    S = Tuple(axis.overlap for axis in pgdg)
    scratch_GG = zeros(Float64, basis.final_dimension, basis.final_dimension)
    product_action_buffer = Ref(Matrix{Float64}(undef, 0, 0))
    product_tile_buffer = Ref(Matrix{Float64}(undef, 0, 0))
    product_block_buffer = Ref(Matrix{Float64}(undef, 0, 0))
    product!(ax, ay, az) = _assemble_terminal_product_operator!(
        scratch_GG, basis, ax, ay, az,
        product_action_buffer, product_tile_buffer, product_block_buffer)
    kinetic_GG = if isnothing(base_kinetic)
        product!(pgdg[1].kinetic, S[2], S[3])
        product!(S[1], pgdg[2].kinetic, S[3])
        product!(S[1], S[2], pgdg[3].kinetic)
        scratch_GG
    else
        _r3_require_size(base_kinetic, (basis.final_dimension, basis.final_dimension),
            "R3 trusted base kinetic dimension mismatch")
        base_kinetic
    end
    kinetic = CRG.transform_augmented_operator(kinetic_GG,
        supplement_blocks_value.mixed.kinetic, supplement_blocks_value.self.kinetic, residual)
    fill!(scratch_GG, 0.0)
    product!(pgdg[1].position, S[2], S[3])
    pos_x = CRG.transform_augmented_operator(scratch_GG,
        supplement_blocks_value.mixed.position.x, supplement_blocks_value.self.position.x, residual)
    fill!(scratch_GG, 0.0)
    product!(S[1], pgdg[2].position, S[3])
    pos_y = CRG.transform_augmented_operator(scratch_GG,
        supplement_blocks_value.mixed.position.y, supplement_blocks_value.self.position.y, residual)
    fill!(scratch_GG, 0.0)
    product!(S[1], S[2], pgdg[3].position)
    pos_z = CRG.transform_augmented_operator(scratch_GG,
        supplement_blocks_value.mixed.position.z, supplement_blocks_value.self.position.z, residual)
    fill!(scratch_GG, 0.0)
    product!(pgdg[1].x2, S[2], S[3])
    x2_x = CRG.transform_augmented_operator(scratch_GG,
        supplement_blocks_value.mixed.x2.x, supplement_blocks_value.self.x2.x, residual)
    fill!(scratch_GG, 0.0)
    product!(S[1], pgdg[2].x2, S[3])
    x2_y = CRG.transform_augmented_operator(scratch_GG,
        supplement_blocks_value.mixed.x2.y, supplement_blocks_value.self.x2.y, residual)
    fill!(scratch_GG, 0.0)
    product!(S[1], S[2], pgdg[3].x2)
    x2_z = CRG.transform_augmented_operator(scratch_GG,
        supplement_blocks_value.mixed.x2.z, supplement_blocks_value.self.x2.z, residual)
    pos = (x = pos_x, y = pos_y, z = pos_z)
    x2 = (x = x2_x, y = x2_y, z = x2_z)
    return (; kinetic, position = pos, x2, supplement_blocks = supplement_blocks_value)
end

function pqs_terminal_residual_gto_augmented_unit_nuclear(basis::CartesianTerminalBasisRealization, bundles, residual::CartesianTerminalResidualGTOAugmentation, atom_locations, nuclear_charges, augmented_products; expansion = nothing, base_unit_nuclear = nothing)
    length(atom_locations) == length(nuclear_charges) || throw(DimensionMismatch("R3-A atom location count must match nuclear charges"))
    _r3_validate_residual_contract(basis, nothing, residual, atom_locations)
    recompute_unit_GG = isnothing(base_unit_nuclear)
    if !recompute_unit_GG
        length(base_unit_nuclear) == length(atom_locations) ||
            throw(DimensionMismatch("R3 trusted base unit nuclear center count mismatch"))
    end
    expansion_value = recompute_unit_GG ? (isnothing(expansion) ?
        getfield(_GB_PARENT, :coulomb_gaussian_expansion)(doacc = false) : expansion) : nothing
    supplement_blocks_value = augmented_products.supplement_blocks
    pgdg = recompute_unit_GG ? Tuple(_nested_axis_pgdg(bundles, axis) for axis in (:x, :y, :z)) : nothing
    U = Matrix{Float64}[]
    gaussian_sum_action_buffer = recompute_unit_GG ? Ref(Matrix{Float64}(undef, 0, 0)) : nothing
    gaussian_sum_tile_buffer = recompute_unit_GG ? Ref(Matrix{Float64}(undef, 0, 0)) : nothing
    gaussian_sum_block_buffer = recompute_unit_GG ? Ref(Matrix{Float64}(undef, 0, 0)) : nothing
    for (center_index, center) in enumerate(CRG.residual_gaussian_float_centers(atom_locations))
        U_GG = if recompute_unit_GG
            matrix = zeros(Float64, basis.final_dimension, basis.final_dimension)
            factors = ntuple(axis -> _r3a_centered_factor_terms(pgdg[axis], expansion_value,
                center[axis]), 3)
            _accumulate_terminal_gaussian_sum!(
                matrix, basis, expansion_value.coefficients, factors[1], factors[2], factors[3],
                gaussian_sum_action_buffer, gaussian_sum_tile_buffer, gaussian_sum_block_buffer)
            matrix
        else
            matrix = base_unit_nuclear[center_index]
            _r3_require_size(matrix, (basis.final_dimension, basis.final_dimension),
                "R3 trusted base unit nuclear dimension mismatch")
            matrix
        end
        U_GA = supplement_blocks_value.mixed.nuclear[center_index]
        U_AA = supplement_blocks_value.self.nuclear[center_index]
        push!(U, CRG.transform_augmented_operator(U_GG, U_GA, U_AA, residual))
    end
    return U
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
    supplement_blocks = nothing,
    base_kinetic = nothing,
    base_unit_nuclear = nothing,
)
    products = pqs_terminal_residual_gto_augmented_products(
        basis, bundles, parent_basis_object, supplement, residual, atom_locations,
        nuclear_charges; expansion, supplement_blocks, base_kinetic)
    U = pqs_terminal_residual_gto_augmented_unit_nuclear(
        basis, bundles, residual, atom_locations, nuclear_charges, products;
        expansion, base_unit_nuclear)
    return (;
        kinetic = products.kinetic,
        nuclear_attraction_unit_by_center = U,
        position = products.position,
        x2 = products.x2,
    )
end
