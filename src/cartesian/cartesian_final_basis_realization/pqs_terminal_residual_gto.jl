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
    residual = CRG._build_residual_gaussian_basis(
        _terminal_residual_finalizer(basis, bundles, supplement),
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
    residual_occupation_cutoff::Real = 1.0e-8,
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
    return CRG._build_residual_gaussian_basis(
        _terminal_residual_finalizer(basis, bundles, supplement), basis.final_dimension, X, S_AA,
        labels, centers, owners; residual_occupation_cutoff, tau_neg_abs,
        tau_neg_rel, tau_merge_abs, tau_merge_rel, orthogonality_atol, identity_atol,
        residual_injection_cutoff)
end

_terminal_residual_product(xs) = foldl(Base.checked_mul, xs; init = 1)
function _terminal_residual_transform(c, rs)
    nx, ny, nz = map(r -> size(r, 2), rs); n = size(c, 2)
    y = rs[3] * reshape(c, nz, :)
    y = permutedims(reshape(y, size(rs[3], 1), ny, nx, n), (2, 1, 3, 4))
    y = rs[2] * reshape(y, ny, :)
    y = permutedims(reshape(y, size(rs[2], 1), size(rs[3], 1), nx, n), (3, 1, 2, 4))
    y = rs[1] * reshape(y, nx, :)
    return reshape(permutedims(reshape(y, size(rs[1], 1), size(rs[2], 1),
        size(rs[3], 1), n), (3, 2, 1, 4)), :, n)
end
function _terminal_residual_anchors(gs, keys, tail)
    pts = Float64[g.center_value + t*g.width for g in gs for t in (-tail,-6.,-2.,0.,2.,6.,tail)]
    append!(pts, [c+t/sqrt(2a) for (a,c,l) in keys for t in (-tail,-6.,-2.,0.,2.,6.,tail)])
    all(isfinite, pts) || throw(ArgumentError("residual quadrature anchors must be finite"))
    return sort!(unique!(pts))
end
function _terminal_residual_axis(layer, gs, keys, pts, order)
    count = Base.checked_mul(order, length(pts)-1)
    rule = eigen(SymTridiagonal(zeros(order), [i/sqrt(4i^2-1) for i in 1:order-1]))
    x, sw = Vector{Float64}(undef,count), Vector{Float64}(undef,count)
    for j in 1:length(pts)-1
        h = (pts[j+1]-pts[j])/2; mid = pts[j]+h
        for i in 1:order
            row = (j-1)*order+i; x[row] = mid+h*rule.values[i]
            sw[row] = sqrt(2h)*abs(rule.vectors[1,i])
        end
    end
    stencil = getfield(_GB_PARENT, :stencil_matrix)(layer); np = size(stencil,2)
    centers = [g.center_value for g in gs]
    widths = [g.width for g in gs]
    ev = zeros(count,np+length(keys))
    for first in 1:8:length(gs)
        cols = first:min(first+7,length(gs))
        prim = [sw[i]*exp(-((x[i]-centers[j])/widths[j])^2/2) for i in eachindex(x), j in cols]
        mul!(view(ev,:,1:np),prim,view(stencil,cols,:),1.,1.)
    end
    pref = getfield(_GB_PARENT, :_qwrg_atomic_shell_prefactor)
    for (j,(a,c,l)) in enumerate(keys)
        scale = pref(a,l)
        for i in eachindex(x)
            d = x[i]-c; ev[i,np+j] = sw[i]*scale*d^l*exp(-a*d*d)
        end
    end
    all(isfinite, ev) || throw(ArgumentError("nonfinite residual axis evaluation"))
    return Matrix(qr!(ev).R)
end
function _terminal_residual_finalizer(basis, bundles, supplement)
    return function (G0,A0,X,S,ta,tr,it,ot)
        entry_rss = Sys.maxrss(); nR = size(A0,2); batch = min(8,nR)
        ax = ntuple(k -> _nested_axis_pgdg(bundles,(:x,:y,:z)[k]),3)
        proxy = _r3a_qw_proxy_layers(bundles)
        layers = ntuple(k -> getproperty(proxy,(:x,:y,:z)[k]),3)
        all(layers[k] === ax[k].auxiliary_layer for k in 1:3) ||
            throw(ArgumentError("residual evaluator does not represent the parent axes"))
        gs = map(layer -> getfield(_GB_PARENT,:primitives)(getfield(_GB_PARENT,:primitive_set)(layer)),layers)
        orbs = supplement.orbitals
        all(o.primitive_normalization === :axiswise_normalized_cartesian_gaussian for o in orbs) ||
            throw(ArgumentError("residual evaluator requires normalized Cartesian primitives"))
        keys = ntuple(k -> sort!(unique([(a,o.center[k],o.angular_powers[k]) for o in orbs for a in o.exponents])),3)
        all(a>0 && isfinite(a) && isfinite(c) && l>=0 for ks in keys for (a,c,l) in ks) ||
            throw(ArgumentError("invalid residual Gaussian axis input"))
        p = map(a -> size(a.overlap,1),ax); kd = map(length,keys); cd = p.+kd
        tails = map(ks -> max(12.,sqrt(2maximum(t[3] for t in ks)+1)+10),keys)
        anchors = ntuple(k -> _terminal_residual_anchors(gs[k],keys[k],tails[k]),3)
        outer = ntuple(k -> _terminal_residual_anchors(gs[k],keys[k],tails[k]+4),3)
        nodes = ntuple(k -> Base.checked_mul(16,max(length(anchors[k]),length(outer[k]))-1),3)
        N,nP,nK = map(_terminal_residual_product,(cd,p,kd)); nG = size(G0,1)
        m = maximum(_terminal_residual_product(i>=stage ? dst[i] : src[i] for i in 1:3)
            for (src,dst) in ((p,cd),(kd,cd),(cd,p)) for stage in 1:4)
        E = big(512)*1024^2+8*(4big(N)*nR+12big(batch)*m+4sum(big.(nodes).*cd)+
            6big(nG)*nR+12big(nR)^2+2big(nP)*batch+2big(nK)*batch)
        E<=12big(1024)^3 && entry_rss+E<=16big(1024)^3 ||
            throw(ArgumentError("stable residual memory admission failed: estimated bytes=$E, entry RSS=$entry_rss"))
        rows = [[(s[1]-1)*p[2]*p[3]+(s[2]-1)*p[3]+s[3] for s in b.support_states] for b in basis.blocks]
        covered = zeros(Int,nP)
        for row in rows; covered[row] .+= 1; end
        all(==(1),covered) || throw(ArgumentError("terminal residual support is not a partition"))
        km = map(ks -> Dict(key=>i for (i,key) in enumerate(ks)),keys)
        function lift(g)
            v = zeros(nP,size(g,2))
            for (b,row) in zip(basis.blocks,rows)
                v[row,:] = isnothing(b.coefficients) ? g[b.column_range,:] : b.coefficients*view(g,b.column_range,:)
            end
            return v
        end
        function restrict(v)
            g = zeros(nG,size(v,2))
            for (b,row) in zip(basis.blocks,rows)
                g[b.column_range,:] = isnothing(b.coefficients) ? v[row,:] : b.coefficients'*view(v,row,:)
            end
            return g
        end
        function vectors(g,a,rs)
            v = _terminal_residual_transform(lift(g),ntuple(k -> view(rs[k],:,1:p[k]),3))
            ct = zeros(nK,size(a,2))
            for (j,o) in enumerate(orbs),(ex,co) in zip(o.exponents,o.coefficients)
                inds = ntuple(k -> km[k][(ex,o.center[k],o.angular_powers[k])],3)
                row = (inds[1]-1)*kd[2]*kd[3]+(inds[2]-1)*kd[3]+inds[3]
                for col in axes(a,2); ct[row,col] += co*a[j,col]; end
            end
            v .+= _terminal_residual_transform(ct,ntuple(k -> view(rs[k],:,p[k]+1:cd[k]),3))
            return v
        end
        localerr = maximum(isnothing(b.coefficients) ? 0. : norm(b.coefficients'*b.coefficients-I) for b in basis.blocks)
        G,A = copy(G0),copy(A0); previous = nothing
        for (order,tailcheck) in ((8,false),(12,false),(16,false),(16,true))
            rs = ntuple(k -> _terminal_residual_axis(layers[k],gs[k],keys[k],tailcheck ? outer[k] : anchors[k],order),3)
            rg = ntuple(k -> view(rs[k],:,1:p[k]),3)
            axiserr = [opnorm(rg[k]'*rg[k]-I,2) for k in 1:3]
            localerr+(prod(1 .+ axiserr)-1)*(1+localerr)<=it ||
                throw(ArgumentError("terminal base is not orthonormal on the residual grid"))
            cross = zeros(nG,nR)
            if order==8
                V = Matrix{Float64}(undef,N,nR)
                for first in 1:batch:nR
                    cols = first:min(first+batch-1,nR); v = vectors(view(G,:,cols),view(A,:,cols),rs)
                    for pass in 1:2
                        delta = restrict(_terminal_residual_transform(v,map(adjoint,rg)))
                        G[:,cols] .-= delta
                        v .-= _terminal_residual_transform(lift(delta),rg)
                    end
                    V[:,cols] = v
                end
                factor = svd(Matrix(qr!(V).R)); vals = factor.S.^2
                threshold = CRG.check_residual_gaussian_metric(vals,ta,tr,"physical residual final merge metric")
                minimum(vals)>threshold || throw(ArgumentError("physical residual final merge metric is near singular"))
                U = factor.V*Diagonal(inv.(factor.S))*factor.V'
                G,A = G*U,A*U; CRG.canonicalize_residual_signs!(A,G)
                V = nothing; GC.gc()
            end
            V = Matrix{Float64}(undef,N,nR)
            for first in 1:batch:nR
                cols = first:min(first+batch-1,nR); v = vectors(view(G,:,cols),view(A,:,cols),rs)
                V[:,cols] = v; cross[:,cols] = restrict(_terminal_residual_transform(v,map(adjoint,rg)))
            end
            gram = V'*V; V = nothing; GC.gc()
            maximum(abs,gram-I)<=it*(1+max(1.,maximum(abs,gram))) && maximum(abs,cross)<=ot ||
                throw(ArgumentError("physical residual identity/cross validation failed"))
            max(opnorm(gram-I,Inf),opnorm(cross,Inf),opnorm(cross',Inf))<=1e-5 ||
                throw(ArgumentError("physical residual row-sum validation failed"))
            if !isnothing(previous)
                dg,dc = gram-previous[1],cross-previous[2]
                maximum(abs,dg)<=it/10 && maximum(abs,dc)<=ot/10 &&
                    max(opnorm(dg,Inf),opnorm(dc,Inf),opnorm(dc',Inf))<=1e-6 ||
                    throw(ArgumentError("physical residual quadrature/tail stabilization failed"))
            end
            previous = (gram,cross)
        end
        return G,A
    end
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
