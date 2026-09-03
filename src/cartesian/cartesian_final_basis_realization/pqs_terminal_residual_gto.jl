const _GB_PARENT = parentmodule(@__MODULE__)
const CGRB = getfield(_GB_PARENT, :CartesianGaussianRawBlocks)
const CartesianTerminalResidualGTOAugmentation = CRG.CartesianResidualGaussianBasis

_r3_require_size(matrix, dims, label) = size(matrix) == dims || throw(DimensionMismatch(label))
_r3_require_close(block, reference, label) = norm(block - reference, Inf) <= 1.0e-10 || throw(ArgumentError(label))
_r3_moment_ok(matrix, dims) = size(matrix) == dims && all(isfinite, matrix) && norm(matrix - transpose(matrix), Inf) <= 1.0e-10
function _r3_validate_pgdg_expansion(bundles, expansion)
    expansion isa getfield(_GB_PARENT, :CoulombGaussianExpansion) ||
        throw(ArgumentError("R3 requires a producer-owned Coulomb expansion"))
    pgdg = Tuple(_nested_axis_pgdg(bundles, axis) for axis in (:x, :y, :z))
    for (axis_name, axis) in zip((:x, :y, :z), pgdg)
        length(axis.exponents) == length(expansion) &&
            axis.exponents == expansion.exponents ||
            throw(ArgumentError("R3 $(axis_name)-axis Coulomb exponent sequence mismatch"))
    end
    return pgdg
end
function _r3_validate_residual_contract(
    base_dimension::Integer,
    supplement,
    residual::CartesianTerminalResidualGTOAugmentation,
    atom_locations,
)
    nG, nR = Int(base_dimension), residual.residual_dimension
    residual.base_dimension == nG ||
        throw(DimensionMismatch("R3 residual base dimension mismatch"))
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
_r3_validate_residual_contract(basis::CartesianTerminalBasisRealization,
    supplement, residual, atom_locations) = _r3_validate_residual_contract(
        basis.final_dimension, supplement, residual, atom_locations)
function _r3_validate_augmented_operator_dimensions(operators, base_hamiltonian, residual, center_count)
    nG, n = residual.base_dimension, residual.base_dimension + residual.residual_dimension
    length(operators.nuclear_attraction_unit_by_center) == center_count || throw(DimensionMismatch("R3 augmented unit nuclear center count mismatch"))
    for matrix in (operators.kinetic, operators.nuclear_attraction_unit_by_center...)
        _r3_require_size(matrix, (n, n), "R3 augmented operator dimension mismatch")
    end
    if CRG.injected_dimension(residual) == 0
        _r3_require_close(view(operators.kinetic, 1:nG, 1:nG), base_hamiltonian.kinetic, "R3 augmented kinetic G-G block mismatch")
        for (matrix, base) in zip(operators.nuclear_attraction_unit_by_center,
                                  base_hamiltonian.nuclear_attraction_unit_by_center)
            _r3_require_close(view(matrix, 1:nG, 1:nG), base, "R3 augmented unit nuclear G-G block mismatch")
        end
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

function _r3_validate_same_construction_nuclear_charges(nuclear_charges, base_hamiltonian)
    charges = Float64.(collect(nuclear_charges))
    base_charges = Float64.(collect(base_hamiltonian.nuclear_charges))
    length(charges) == length(base_charges) ||
        throw(DimensionMismatch("R3 same-construction nuclear charge count must match base Hamiltonian"))
    charges == base_charges ||
        throw(ArgumentError("R3 same-construction nuclear charges must match base Hamiltonian"))
    return nothing
end

function pqs_terminal_residual_gto_augmented_hamiltonian(
    base_hamiltonian,
    basis::CartesianTerminalBasisRealization,
    bundles,
    residual::CartesianTerminalResidualGTOAugmentation,
    augmented_operators;
    expansion,
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

function pqs_terminal_residual_gto_augmented_vee(base_hamiltonian, basis::CartesianTerminalBasisRealization, bundles, residual::CartesianTerminalResidualGTOAugmentation, augmented_operators; expansion)
    center_count = _r3_validate_base_hamiltonian(base_hamiltonian, residual)
    atom_locations = NTuple{3,Float64}[CRG.residual_gaussian_center(
        view(base_hamiltonian.nuclear_positions, index, :)) for index in 1:center_count]
    _r3_validate_residual_contract(basis, nothing, residual, atom_locations)
    _r3_validate_augmented_operator_dimensions(augmented_operators, base_hamiltonian, residual, center_count)
    _r3_validate_pgdg_expansion(bundles, expansion)
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
    expansion,
)
    _r3_validate_same_construction_nuclear_charges(nuclear_charges, base_hamiltonian)
    _r3_validate_pgdg_expansion(bundles, expansion)
    supplement_blocks = _r3a_qw_blocks(basis, bundles, supplement, atom_locations,
        expansion)
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
        expansion, supplement_blocks,
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
    residual.residual_injection_cutoff <= 0 || throw(ArgumentError("R3 augmented artifact writer does not support injected residual sectors"))
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

function _r3a_project_parent_ga(composition::CartesianParentBackedInjectedComposition,
    parent_ga)
    terminal = _r3a_project_parent_ga(composition.terminal_basis, parent_ga)
    residual = zeros(Float64,
        composition.parent_backed_dimension - size(terminal, 1), size(parent_ga, 2))
    for (prf, range) in zip(composition.parent_residual_blocks,
            composition.parent_residual_column_ranges)
        residual[range, :] .= transpose(prf.coefficients) *
            @view(parent_ga[prf.support_indices, :])
    end
    return vcat(terminal, residual)
end

_r3a_validate_representation(::CartesianTerminalBasisRealization, bundles) = nothing
_r3a_validate_representation(
    composition::CartesianParentBackedInjectedComposition, bundles) =
    _validate_parent_backed_injected_composition(composition, bundles)

function _r3a_qw_blocks(representation, bundles, supplement, atom_locations, expansion)
    _r3a_validate_representation(representation, bundles)
    donor = _r3a_qw_supplement(supplement)
    proxy = _r3a_qw_proxy_layers(bundles)
    non_nuclear = CGRB.gaussian_non_nuclear_raw_blocks(proxy, donor, expansion)
    nuclear = CGRB.gaussian_nuclear_raw_blocks_by_center(proxy, donor, expansion,
        CRG.residual_gaussian_float_centers(atom_locations))
    return (;
        mixed = (;
            overlap = _r3a_project_parent_ga(representation, non_nuclear.ga.overlap),
            kinetic = _r3a_project_parent_ga(representation, non_nuclear.ga.kinetic),
            position = (x = _r3a_project_parent_ga(representation, non_nuclear.ga.position.x),
                y = _r3a_project_parent_ga(representation, non_nuclear.ga.position.y),
                z = _r3a_project_parent_ga(representation, non_nuclear.ga.position.z)),
            x2 = (x = _r3a_project_parent_ga(representation, non_nuclear.ga.x2.x),
                y = _r3a_project_parent_ga(representation, non_nuclear.ga.x2.y),
                z = _r3a_project_parent_ga(representation, non_nuclear.ga.x2.z)),
            nuclear = [_r3a_project_parent_ga(representation, matrix) for matrix in nuclear.ga],
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


function parent_backed_injected_residual_gto_augmentation(
    composition::CartesianParentBackedInjectedComposition,
    bundles,
    supplement,
    nuclei;
    expansion,
)
    _validate_parent_backed_injected_composition(composition, bundles)
    nuclei_value = CRG.residual_gaussian_float_centers(nuclei)
    labels = CRG.residual_gaussian_candidate_labels(supplement)
    centers = CRG.residual_gaussian_candidate_centers(supplement)
    owners = Int[CRG.residual_candidate_owner(center, nuclei_value) for center in centers]
    blocks = _r3a_qw_blocks(composition, bundles, supplement, nuclei, expansion)
    X = blocks.mixed.overlap
    S_AA = blocks.self.overlap
    residual = CRG.build_residual_gaussian_basis(
        composition.parent_backed_dimension, X, S_AA, labels, centers, owners;
        residual_occupation_cutoff = 1.0e-10,
        residual_injection_cutoff = 0.0,
        residual_compactness = nothing)
    return (; residual, mixed_overlap = X, supplement_overlap = S_AA,
        supplement_blocks = blocks)
end


function _validate_parent_backed_residual_augmentation(
    composition::CartesianParentBackedInjectedComposition,
    bundles,
    augmentation,
)
    _validate_parent_backed_injected_composition(composition, bundles)
    hasproperty(augmentation, :residual) &&
        hasproperty(augmentation, :mixed_overlap) &&
        hasproperty(augmentation, :supplement_overlap) &&
        hasproperty(augmentation, :supplement_blocks) || throw(ArgumentError(
        "parent-backed augmentation is structurally incomplete"))
    blocks = augmentation.supplement_blocks
    augmentation.mixed_overlap == blocks.mixed.overlap || throw(ArgumentError(
        "parent-backed augmentation mixed overlap does not match its raw blocks"))
    augmentation.supplement_overlap == blocks.self.overlap || throw(ArgumentError(
        "parent-backed augmentation supplement overlap does not match its raw blocks"))
    residual = augmentation.residual
    size(augmentation.mixed_overlap) ==
        (composition.parent_backed_dimension, residual.candidate_count) ||
        throw(DimensionMismatch("parent-backed augmentation mixed-overlap dimensions differ"))
    size(augmentation.supplement_overlap) ==
        (residual.candidate_count, residual.candidate_count) ||
        throw(DimensionMismatch(
            "parent-backed augmentation supplement-overlap dimensions differ"))
    norm(residual.T_G + augmentation.mixed_overlap * residual.T_A, Inf) <= 1.0e-10 ||
        throw(ArgumentError("parent-backed augmentation residual projection is stale"))
    residual_metric = CRG.residual_gaussian_overlap(residual.T_G, residual.T_A,
        augmentation.mixed_overlap, augmentation.supplement_overlap)
    norm(residual_metric - I, Inf) <= 5.0e-8 || throw(ArgumentError(
        "parent-backed augmentation residual metric is not identity"))
    return residual, blocks
end

function parent_backed_injected_residual_gto_augmented_operators(
    composition::CartesianParentBackedInjectedComposition,
    bundles,
    supplement,
    augmentation,
    atom_locations,
    nuclear_charges;
    expansion,
)
    residual, blocks = _validate_parent_backed_residual_augmentation(
        composition, bundles, augmentation)
    locations = CRG.residual_gaussian_float_centers(atom_locations)
    charges = Float64.(nuclear_charges)
    length(locations) == length(charges) || throw(DimensionMismatch(
        "parent-backed atom location and nuclear charge counts differ"))
    all(isfinite, charges) || throw(ArgumentError(
        "parent-backed nuclear charges must be finite"))
    center_count = length(locations)
    length(blocks.mixed.nuclear) == center_count || throw(DimensionMismatch(
        "parent-backed mixed nuclear-block count differs from atom count"))
    length(blocks.self.nuclear) == center_count || throw(DimensionMismatch(
        "parent-backed self nuclear-block count differs from atom count"))
    _r3_validate_residual_contract(composition.parent_backed_dimension,
        supplement, residual, locations)
    _r3_validate_pgdg_expansion(bundles, expansion)
    one_body = parent_backed_injected_one_body_operators(composition, bundles,
        locations, charges; expansion)
    length(one_body.nuclear_attraction_unit_by_center) == center_count ||
        throw(DimensionMismatch(
            "parent-backed one-body unit-nuclear count differs from atom count"))
    kinetic = CRG.transform_augmented_operator(one_body.kinetic,
        blocks.mixed.kinetic, blocks.self.kinetic, residual)
    unit_nuclear = Matrix{Float64}[CRG.transform_augmented_operator(
        one_body.nuclear_attraction_unit_by_center[index],
        blocks.mixed.nuclear[index], blocks.self.nuclear[index], residual)
        for index in 1:center_count]
    position = (;
        x = CRG.transform_augmented_operator(one_body.position.x,
            blocks.mixed.position.x, blocks.self.position.x, residual),
        y = CRG.transform_augmented_operator(one_body.position.y,
            blocks.mixed.position.y, blocks.self.position.y, residual),
        z = CRG.transform_augmented_operator(one_body.position.z,
            blocks.mixed.position.z, blocks.self.position.z, residual))
    x2 = (;
        x = CRG.transform_augmented_operator(one_body.x2.x,
            blocks.mixed.x2.x, blocks.self.x2.x, residual),
        y = CRG.transform_augmented_operator(one_body.x2.y,
            blocks.mixed.x2.y, blocks.self.x2.y, residual),
        z = CRG.transform_augmented_operator(one_body.x2.z,
            blocks.mixed.x2.z, blocks.self.x2.z, residual))
    H1 = copy(kinetic)
    for index in eachindex(charges)
        H1 .+= charges[index] .* unit_nuclear[index]
    end
    matrices = Any[kinetic, unit_nuclear..., H1,
        position.x, position.y, position.z, x2.x, x2.y, x2.z]
    all(matrix -> all(isfinite, matrix) &&
        norm(matrix - transpose(matrix), Inf) <= 1.0e-10, matrices) ||
        throw(ArgumentError(
            "parent-backed augmented one-body matrices must be finite and symmetric"))
    return (; kinetic, nuclear_attraction_unit_by_center = unit_nuclear,
        one_body_hamiltonian = H1, position, x2,
        parent_one_body = one_body, supplement_blocks = blocks)
end

function _validate_parent_backed_interaction_operators(
    composition::CartesianParentBackedInjectedComposition,
    bundles,
    augmentation,
    operators,
)
    residual, blocks = _validate_parent_backed_residual_augmentation(
        composition, bundles, augmentation)
    hasproperty(operators, :position) && hasproperty(operators, :x2) &&
        hasproperty(operators, :parent_one_body) || throw(ArgumentError(
        "parent-backed interaction operators are structurally incomplete"))
    nB = composition.parent_backed_dimension
    n = nB + residual.residual_dimension
    for family in (:position, :x2), axis in (:x, :y, :z)
        supplied = getproperty(getproperty(operators, family), axis)
        parent = getproperty(getproperty(operators.parent_one_body, family), axis)
        size(parent) == (nB, nB) || throw(DimensionMismatch(
            "parent-backed interaction parent moment dimensions differ"))
        size(supplied) == (n, n) || throw(DimensionMismatch(
            "parent-backed interaction moment dimensions differ"))
        all(isfinite, parent) && all(isfinite, supplied) || throw(ArgumentError(
            "parent-backed interaction moments must be finite"))
        expected = CRG.transform_augmented_operator(parent,
            getproperty(getproperty(blocks.mixed, family), axis),
            getproperty(getproperty(blocks.self, family), axis), residual)
        expected .-= supplied
        norm(expected, Inf) <= 1.0e-10 ||
            throw(ArgumentError(
                "parent-backed interaction moments do not match their augmentation"))
    end
    return residual, blocks
end

function parent_backed_injected_residual_gto_interaction(
    composition::CartesianParentBackedInjectedComposition,
    bundles,
    augmentation,
    operators;
    expansion,
)
    residual, _ = _validate_parent_backed_interaction_operators(
        composition, bundles, augmentation, operators)
    _r3_validate_pgdg_expansion(bundles, expansion)
    base = parent_backed_injected_interaction_base_blocks(
        composition, bundles; expansion)
    mwg = CRG.parent_backed_injected_mwg_blocks(
        composition.terminal_basis, bundles, residual, operators; expansion)
    nG = composition.terminal_basis.final_dimension
    nB = composition.parent_backed_dimension
    nR = nB - nG
    nE = residual.residual_dimension
    size(mwg.terminal_external) == (nG, nE) &&
        size(mwg.parent_residual_external) == (nR, nE) &&
        size(mwg.external_external) == (nE, nE) || throw(DimensionMismatch(
        "parent-backed separated MWG block dimensions differ"))
    g_range = 1:nG
    r_range = (nG + 1):nB
    e_range = (nB + 1):(nB + nE)
    V = zeros(Float64, nB + nE, nB + nE)
    V[g_range, g_range] .= base.terminal
    V[g_range, r_range] .= base.terminal_parent_residual
    V[r_range, g_range] .= transpose(base.terminal_parent_residual)
    V[r_range, r_range] .= base.parent_residual
    V[g_range, e_range] .= mwg.terminal_external
    V[e_range, g_range] .= transpose(mwg.terminal_external)
    V[r_range, e_range] .= mwg.parent_residual_external
    V[e_range, r_range] .= transpose(mwg.parent_residual_external)
    V[e_range, e_range] .= mwg.external_external
    all(isfinite, V) && norm(V - transpose(V), Inf) <= 1.0e-10 ||
        throw(ArgumentError(
            "parent-backed injected interaction must be finite and symmetric"))
    diagnostics = (; base = base.diagnostics,
        parent_residual_centers = mwg.parent_residual_centers,
        parent_residual_widths = mwg.parent_residual_widths,
        external_centers = mwg.external_centers,
        external_widths = mwg.external_widths)
    return (; electron_electron_ida = V, diagnostics)
end

function parent_backed_injected_gaussian_potential_raw_blocks(
    composition::CartesianParentBackedInjectedComposition,
    bundles,
    supplement,
    augmentation,
    potential_expansion,
    center,
)
    _validate_parent_backed_residual_augmentation(composition, bundles, augmentation)
    placement = Tuple(Float64.(center))
    length(placement) == 3 && all(isfinite, placement) || throw(ArgumentError(
        "parent-backed fitted-potential center must be a finite 3-vector"))
    proxy = _r3a_qw_proxy_layers(bundles)
    raw = CGRB.placed_spherical_gaussian_potential_raw_blocks(
        composition.terminal_basis, bundles, proxy, _r3a_qw_supplement(supplement),
        potential_expansion, placement)
    pgdg = Tuple(_nested_axis_pgdg(bundles, axis) for axis in (:x, :y, :z))
    factors = ntuple(axis -> _r3a_centered_factor_terms(
        pgdg[axis], potential_expansion, placement[axis]), 3)
    ranges = composition.parent_residual_column_ranges
    prf = _parent_residual_gaussian_sum_blocks(composition.terminal_basis,
        composition.parent_residual_blocks, ranges, potential_expansion.coefficients,
        factors...; scale = 1.0)
    GG = _parent_backed_operator_matrix(raw.GG, prf, composition)
    GA = _r3a_project_parent_ga(composition, raw.GA)
    nB = composition.parent_backed_dimension
    nA = augmentation.residual.candidate_count
    size(GG) == (nB, nB) && size(GA) == (nB, nA) && size(raw.AA) == (nA, nA) ||
        throw(DimensionMismatch(
            "parent-backed fitted-potential raw block dimensions differ"))
    all(isfinite, GG) && all(isfinite, GA) && all(isfinite, raw.AA) ||
        throw(ArgumentError("parent-backed fitted-potential raw blocks must be finite"))
    return (; GG, GA, AA = raw.AA)
end

function pqs_terminal_residual_gto_augmentation(
    basis::CartesianTerminalBasisRealization,
    bundles,
    supplement,
    nuclei;
    residual_occupation_cutoff::Real = 1.0e-6,
    tau_neg_abs::Real = 1.0e-12,
    tau_neg_rel::Real = 1.0e-12,
    tau_merge_abs::Real = 1.0e-12,
    tau_merge_rel::Real = 1.0e-12,
    orthogonality_atol::Real = 1.0e-10, identity_atol::Real = 5.0e-8,
    residual_injection_cutoff::Real = 0.0,
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
        tau_neg_rel, tau_merge_abs, tau_merge_rel, orthogonality_atol, identity_atol,
        residual_injection_cutoff)
end

function _r3a_centered_factor_terms(axis, expansion, center)
    center == axis.center && Float64.(axis.exponents) == Float64.(expansion.exponents) &&
        return axis.gaussian_factor_terms
    ops = getfield(_GB_PARENT, :mapped_ordinary_one_body_operators)(
        axis.basis; exponents = expansion.exponents, center, backend = axis.backend)
    return ops.gaussian_factors
end

function pqs_terminal_residual_gto_augmented_products(basis::CartesianTerminalBasisRealization, bundles, parent_basis_object, supplement, residual::CartesianTerminalResidualGTOAugmentation, atom_locations, nuclear_charges; expansion, supplement_blocks = nothing, base_kinetic = nothing)
    length(atom_locations) == length(nuclear_charges) || throw(DimensionMismatch("R3-A atom location count must match nuclear charges"))
    _r3_validate_residual_contract(basis, supplement, residual, atom_locations)
    _r3_validate_pgdg_expansion(bundles, expansion)
    supplement_blocks_value = isnothing(supplement_blocks) ?
        _r3a_qw_blocks(basis, bundles, supplement, atom_locations, expansion) :
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

function pqs_terminal_residual_gto_augmented_unit_nuclear(basis::CartesianTerminalBasisRealization, bundles, residual::CartesianTerminalResidualGTOAugmentation, atom_locations, nuclear_charges, augmented_products; expansion, base_unit_nuclear = nothing)
    length(atom_locations) == length(nuclear_charges) || throw(DimensionMismatch("R3-A atom location count must match nuclear charges"))
    _r3_validate_residual_contract(basis, nothing, residual, atom_locations)
    recompute_unit_GG = isnothing(base_unit_nuclear)
    if !recompute_unit_GG
        length(base_unit_nuclear) == length(atom_locations) ||
            throw(DimensionMismatch("R3 trusted base unit nuclear center count mismatch"))
    end
    pgdg = _r3_validate_pgdg_expansion(bundles, expansion)
    supplement_blocks_value = augmented_products.supplement_blocks
    U = Matrix{Float64}[]
    gaussian_sum_action_buffer = recompute_unit_GG ? Ref(Matrix{Float64}(undef, 0, 0)) : nothing
    gaussian_sum_tile_buffer = recompute_unit_GG ? Ref(Matrix{Float64}(undef, 0, 0)) : nothing
    gaussian_sum_block_buffer = recompute_unit_GG ? Ref(Matrix{Float64}(undef, 0, 0)) : nothing
    for (center_index, center) in enumerate(CRG.residual_gaussian_float_centers(atom_locations))
        U_GG = if recompute_unit_GG
            matrix = zeros(Float64, basis.final_dimension, basis.final_dimension)
            factors = ntuple(axis -> _r3a_centered_factor_terms(pgdg[axis], expansion,
                center[axis]), 3)
            _accumulate_terminal_gaussian_sum!(
                matrix, basis, expansion.coefficients, factors[1], factors[2], factors[3],
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
    expansion,
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
