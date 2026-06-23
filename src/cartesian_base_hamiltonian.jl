const _CARTESIAN_BASE_SYSTEM_KEYS = Set((:atom_symbols, :nuclear_charges, :atom_locations, :nup, :ndn))
const _CARTESIAN_BASE_H_BASIS_REQUIRED_KEYS = Set((:q, :core_spacing, :radius))
const _CARTESIAN_BASE_H2_BASIS_REQUIRED_KEYS = Set((:q, :core_spacing, :xmax_parallel, :xmax_transverse))
const _CARTESIAN_BASE_OPTIONAL_BASIS_KEYS = Set((:parent_axis_family, :reference_spacing, :tail_spacing))
const _CARTESIAN_BASE_H_OPTIONAL_BASIS_KEYS = union(_CARTESIAN_BASE_OPTIONAL_BASIS_KEYS, Set((:d,)))
_cartesian_base_check_keys(input, expected, label) =
    Set(Symbol.(keys(input))) == expected ||
        throw(ArgumentError("$(label) has unsupported keys"))

function _cartesian_base_check_basis_keys(basis, required, optional = _CARTESIAN_BASE_OPTIONAL_BASIS_KEYS)
    supplied = Set(Symbol.(keys(basis)))
    required ⊆ supplied && supplied ⊆ union(required, optional) ||
        throw(ArgumentError("basis has missing or unsupported keys"))
end

_cartesian_base_get_positive(basis, key, default) =
    haskey(basis, key) ? _cartesian_base_positive(getproperty(basis, key), "basis.$(key)") : default

function _cartesian_base_parent_axis_family(basis)
    family = haskey(basis, :parent_axis_family) ? basis.parent_axis_family : :G10
    family === :G10 || throw(ArgumentError("only parent_axis_family=:G10 is supported"))
    return family
end

function _cartesian_base_positive(value, label)
    x = Float64(value)
    isfinite(x) && x > 0.0 || throw(ArgumentError("$(label) must be finite and positive"))
    return x
end

function _cartesian_base_atom_mapping_d(basis, core_spacing)
    haskey(basis, :d) || return core_spacing
    d = _cartesian_base_positive(basis.d, "basis.d")
    d == core_spacing || throw(ArgumentError("legacy basis.d must equal basis.core_spacing"))
    return core_spacing
end

function _cartesian_base_q(value)
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("basis.q must be a positive integer"))
    q = Int(value)
    q > 0 || throw(ArgumentError("basis.q must be positive"))
    return q
end

function _cartesian_base_location(value)
    value isa Tuple || throw(ArgumentError("nuclear center must be a fixed 3-tuple"))
    length(value) == 3 || throw(ArgumentError("nuclear center must have three coordinates"))
    location = ntuple(axis -> Float64(value[axis]), 3)
    all(isfinite, location) || throw(ArgumentError("nuclear coordinates must be finite"))
    return location
end

function _cartesian_base_inputs(system::NamedTuple, basis::NamedTuple)
    _cartesian_base_check_keys(system, _CARTESIAN_BASE_SYSTEM_KEYS, "system")
    for field in (:atom_symbols, :nuclear_charges, :atom_locations)
        getproperty(system, field) isa AbstractVector ||
            throw(ArgumentError("system.$(field) must be an AbstractVector"))
    end
    symbols = String.(system.atom_symbols)
    charges = Float64.(system.nuclear_charges)
    locations = _cartesian_base_location.(system.atom_locations)
    length(symbols) == length(charges) == length(locations) ||
        throw(ArgumentError("center collection lengths must match"))
    all(isfinite, charges) && all(>(0.0), charges) ||
        throw(ArgumentError("nuclear charges must be finite and positive"))
    system.nup isa Integer && system.ndn isa Integer &&
        !(system.nup isa Bool) && !(system.ndn isa Bool) ||
        throw(ArgumentError("nup and ndn must be nonnegative integers"))
    nup = Int(system.nup)
    ndn = Int(system.ndn)
    nup >= 0 && ndn >= 0 ||
        throw(ArgumentError("nup and ndn must be nonnegative integers"))
    base = (; symbols, charges, locations, nup, ndn)
    if length(symbols) == 1
        _cartesian_base_check_basis_keys(
            basis, _CARTESIAN_BASE_H_BASIS_REQUIRED_KEYS, _CARTESIAN_BASE_H_OPTIONAL_BASIS_KEYS)
        symbols == ["H"] && charges == [1.0] ||
            throw(ArgumentError("only origin-centered H is supported"))
        locations[1] == (0.0, 0.0, 0.0) ||
            throw(ArgumentError("H must be centered at (0,0,0)"))
        nup + ndn == 1 || throw(ArgumentError("H requires one electron"))
        core_spacing = _cartesian_base_positive(basis.core_spacing, "basis.core_spacing")
        return merge(base, (; kind = :h, q = _cartesian_base_q(basis.q),
            core_spacing,
            radius = _cartesian_base_positive(basis.radius, "basis.radius"),
            d = _cartesian_base_atom_mapping_d(basis, core_spacing), xmax_parallel = nothing,
            xmax_transverse = nothing, parent_axis_family = _cartesian_base_parent_axis_family(basis),
            reference_spacing = _cartesian_base_get_positive(basis, :reference_spacing, 1.0),
            tail_spacing = _cartesian_base_get_positive(basis, :tail_spacing, 10.0)))
    elseif length(symbols) == 2
        _cartesian_base_check_basis_keys(basis, _CARTESIAN_BASE_H2_BASIS_REQUIRED_KEYS)
        symbols == ["H", "H"] && charges == [1.0, 1.0] ||
            throw(ArgumentError("only z-axis H2 is supported"))
        all(location -> location[1] == 0.0 && location[2] == 0.0, locations) ||
            throw(ArgumentError("H2 centers must lie on the Cartesian z axis"))
        locations[1][3] != locations[2][3] ||
            throw(ArgumentError("H2 centers must be distinct"))
        nup == 1 && ndn == 1 || throw(ArgumentError("H2 requires nup=1 and ndn=1"))
        return merge(base, (; kind = :h2, q = _cartesian_base_q(basis.q),
            core_spacing = _cartesian_base_positive(basis.core_spacing, "basis.core_spacing"),
            radius = nothing, d = nothing,
            xmax_parallel = _cartesian_base_positive(basis.xmax_parallel, "basis.xmax_parallel"),
            xmax_transverse = _cartesian_base_positive(basis.xmax_transverse, "basis.xmax_transverse"),
            parent_axis_family = _cartesian_base_parent_axis_family(basis),
            reference_spacing = _cartesian_base_get_positive(basis, :reference_spacing, 1.0),
            tail_spacing = _cartesian_base_get_positive(basis, :tail_spacing, 10.0)))
    end
    throw(ArgumentError("only H and H2 are supported"))
end

function _cartesian_base_route(kind)
    route_shape = kind === :h ? (:pqs_left, :product, :pqs_right) :
        (:atom_contact_core, :shared_shell_1, :shared_shell_2)
    route_kind = kind === :h ? :one_center_fixed_q_complete_core_shell :
        :bond_aligned_diatomic_independent_pqs_source_box_core_shell
    return (;
        route_family = :pqs_source_box, route_kind, route_shape, product_body_rule = :centered_single_z_slab,
        pqs_retained_rule = :boundary_comx_product_mode_selection,
        product_retained_rule = :product_doside_retained_unit, terms = (:overlap,),
        pair_factor_normalization = :density_normalized,
        white_lindsey_route_shape = (:standard_cartesian_units, :low_order_comx_coarsening),
        white_lindsey_mapping_rule = :standard_unit_backbone_mapping_family,
        white_lindsey_nesting_rule = :unit_box_low_order_comx_coarsening,
        white_lindsey_retained_rule = :low_order_unit_comx_retained_basis,
        white_lindsey_operator_rule = :low_order_unit_operator_blocks,
        supplement_policy = nothing, run_final_basis = false, run_h1 = false, run_h1_j = false)
end

function _cartesian_base_stages(input)
    atom_symbols = Tuple(Symbol.(input.symbols))
    nuclear_charges = Tuple(input.charges)
    atom_locations = Tuple(input.locations)
    spacing_inputs = (;
        q = input.q, n_s = input.q, reference_spacing = input.reference_spacing, tail_spacing = input.tail_spacing,
        q_to_core_spacing_rule = :standard_pqs_ns_equals_q, core_spacing = input.core_spacing,
        xmax_parallel = input.xmax_parallel, xmax_transverse = input.xmax_transverse)
    parent_inputs = (;
        parent_axis_bundle_backend = :pgdg_localized_experimental, parent_axis_family = input.parent_axis_family,
        parent_mapping_tail_spacing = input.tail_spacing,
        parent_mapping_rule = input.kind === :h ? :white_lindsey_atomic_mapping : :identity_mapping,
        parent_mapping_Z = input.kind === :h ? only(input.charges) : nothing,
        parent_mapping_d = input.kind === :h ? input.core_spacing : nothing)
    if input.kind === :h
        side = 2 * input.q + 1
        system_inputs = (; atom_symbols, nuclear_charges, atom_locations,
            nup = input.nup, ndn = input.ndn, bond_axis = nothing, bond_length = nothing,
            radius = input.radius,
            parent_axis_counts = (x = side, y = side, z = side),
            map_backend = :pgdg_localized_experimental)
    else
        z = (input.locations[1][3], input.locations[2][3])
        system_inputs = (; atom_symbols, nuclear_charges, atom_locations,
            nup = input.nup, ndn = input.ndn, bond_axis = :z, bond_length = abs(z[2] - z[1]),
            radius = max(input.xmax_parallel, input.xmax_transverse),
            parent_axis_counts = nothing, map_backend = :pgdg_localized_experimental)
    end
    system = cartesian_system(system_inputs)
    recipe = cartesian_recipe(_cartesian_base_route(input.kind))
    parent = cartesian_parent(system, spacing_inputs, parent_inputs, recipe)
    shells = cartesian_shells(parent, spacing_inputs, recipe)
    units = cartesian_units(parent, shells, recipe)
    transforms = cartesian_transforms(units, recipe)
    return parent, transforms
end

_cartesian_base_atom_location_matrix(input) =
    [input.locations[row][axis] for row in eachindex(input.locations), axis in 1:3]

_cartesian_base_axis_counts_tuple(axis_counts) =
    (Int(axis_counts.x), Int(axis_counts.y), Int(axis_counts.z))

function _cartesian_base_write_hamiltonian(path, ham, input, parent)
    write_cartesian_ida_hamiltonian(path, ham)
    route = input.kind === :h ? :one_center_pqs_base : :z_axis_diatomic_pqs_base
    mapping_kind = input.kind === :h ? :white_lindsey_atomic_mapping : :multicenter_pqs_mapping
    values = (; provenance_version = 1, producer = :cartesian_base_hamiltonian,
        route, q = input.q, core_spacing = input.core_spacing,
        reference_spacing = input.reference_spacing, tail_spacing = input.tail_spacing,
        parent_axis_family = input.parent_axis_family,
        parent_axis_counts = _cartesian_base_axis_counts_tuple(parent.axis_counts),
        mapping_kind, mapping_d = input.d, radius = input.radius,
        xmax_parallel = input.xmax_parallel, xmax_transverse = input.xmax_transverse,
        atom_symbols = copy(input.symbols), nuclear_charges = copy(input.charges),
        atom_locations = _cartesian_base_atom_location_matrix(input),
        nup = input.nup, ndn = input.ndn, final_dimension = size(ham.kinetic, 1))
    jldopen(String(path), "r+") do file
        for key in keys(values)
            file["producer_provenance/$(key)"] = getproperty(values, key)
        end
    end
    return path
end

function cartesian_base_hamiltonian(
    system::NamedTuple;
    basis::NamedTuple,
    hamfile::Union{Nothing,AbstractString} = nothing,
)::CartesianIDAHamiltonian{Float64}
    !(!isnothing(hamfile) && isempty(String(hamfile))) ||
        throw(ArgumentError("hamfile must not be empty"))
    base = cartesian_base_working_basis(system; basis)
    return cartesian_base_hamiltonian_assembly(base; hamfile)
end

_cartesian_base_expansion() = coulomb_gaussian_expansion(doacc = false)

_cartesian_base_pgdg(base) =
    Tuple(_nested_axis_pgdg(base.parent.parent_axis_bundle_object, axis) for axis in (:x, :y, :z))

function cartesian_base_working_basis(system::NamedTuple; basis::NamedTuple, supplemented::Bool = false)
    input = supplemented ? _cartesian_r3_diatomic_inputs(system, basis) :
        _cartesian_base_inputs(system, basis)
    parent, transforms = _cartesian_base_stages(input)
    return (; input, parent, terminal_basis = transforms.terminal_basis_realization)
end

cartesian_base_products(base) =
    _pqs_source_box_route_driver_terminal_products(base.terminal_basis, _cartesian_base_pgdg(base))

cartesian_base_unit_nuclear(base) =
    _pqs_source_box_route_driver_terminal_unit_nuclear(base.terminal_basis,
        _cartesian_base_expansion(), _cartesian_base_pgdg(base), base.input.locations)

cartesian_base_vee(base) =
    _pqs_source_box_route_driver_terminal_vee(base.terminal_basis,
        _cartesian_base_expansion(), _cartesian_base_pgdg(base))

function cartesian_base_hamiltonian_assembly(base; hamfile::Union{Nothing,AbstractString} = nothing)::CartesianIDAHamiltonian{Float64}
    products = cartesian_base_products(base)
    unit_nuclear = cartesian_base_unit_nuclear(base)
    vee = cartesian_base_vee(base)
    return cartesian_base_hamiltonian_assembly(
        base, products, unit_nuclear, vee; hamfile)
end

function cartesian_base_hamiltonian_assembly(base, products, unit_nuclear, vee; hamfile::Union{Nothing,AbstractString} = nothing)::CartesianIDAHamiltonian{Float64}
    !(!isnothing(hamfile) && isempty(String(hamfile))) ||
        throw(ArgumentError("hamfile must not be empty"))
    input = base.input
    parent = base.parent
    ham = CartesianIDAHamiltonian(products.kinetic, unit_nuclear, vee, input.nup, input.ndn;
        nuclear_charges = input.charges, nuclear_positions = input.locations)
    isnothing(hamfile) ||
        _cartesian_base_write_hamiltonian(String(hamfile), ham, input, parent)
    return ham
end

function _cartesian_r3_diatomic_inputs(system::NamedTuple, basis::NamedTuple)
    _cartesian_base_check_keys(system, _CARTESIAN_BASE_SYSTEM_KEYS, "system")
    _cartesian_base_check_basis_keys(basis, _CARTESIAN_BASE_H2_BASIS_REQUIRED_KEYS)
    for field in (:atom_symbols, :nuclear_charges, :atom_locations)
        getproperty(system, field) isa AbstractVector ||
            throw(ArgumentError("system.$(field) must be an AbstractVector"))
    end
    symbols = String.(system.atom_symbols)
    charges = Float64.(system.nuclear_charges)
    locations = _cartesian_base_location.(system.atom_locations)
    length(symbols) == length(charges) == length(locations) == 2 ||
        throw(ArgumentError("R3 usability facade supports two-center systems only"))
    symbols[1] == symbols[2] ||
        throw(ArgumentError("heteronuclear supplements are unsupported"))
    all(charge -> isfinite(charge) && charge > 0.0, charges) ||
        throw(ArgumentError("nuclear charges must be finite and positive"))
    charges[1] == charges[2] ||
        throw(ArgumentError("homonuclear supplemented diatomics require equal nuclear charges"))
    system.nup isa Integer && system.ndn isa Integer &&
        !(system.nup isa Bool) && !(system.ndn isa Bool) ||
        throw(ArgumentError("nup and ndn must be nonnegative integers"))
    nup, ndn = Int(system.nup), Int(system.ndn)
    nup >= 0 && ndn >= 0 || throw(ArgumentError("nup and ndn must be nonnegative integers"))
    total_charge = sum(charges)
    electron_count = round(Int, total_charge)
    isapprox(total_charge, electron_count; atol = 1.0e-12, rtol = 0.0) ||
        throw(ArgumentError("neutral all-electron count requires integer total nuclear charge"))
    nup + ndn == electron_count ||
        throw(ArgumentError("neutral all-electron count requires nup + ndn == total nuclear charge"))
    all(location -> all(isfinite, location) && location[1] == 0.0 && location[2] == 0.0, locations) ||
        throw(ArgumentError("R3 usability facade supports Cartesian z-axis diatomics only"))
    locations[1][3] != locations[2][3] ||
        throw(ArgumentError("diatomic centers must be distinct"))
    return (; symbols, charges, locations, nup, ndn, kind = :h2,
        q = _cartesian_base_q(basis.q),
        core_spacing = _cartesian_base_positive(basis.core_spacing, "basis.core_spacing"),
        radius = nothing, d = nothing,
        xmax_parallel = _cartesian_base_positive(basis.xmax_parallel, "basis.xmax_parallel"),
        xmax_transverse = _cartesian_base_positive(basis.xmax_transverse, "basis.xmax_transverse"),
        parent_axis_family = _cartesian_base_parent_axis_family(basis),
        reference_spacing = _cartesian_base_get_positive(basis, :reference_spacing, 1.0),
        tail_spacing = _cartesian_base_get_positive(basis, :tail_spacing, 10.0))
end

function _cartesian_r3_supplement_inputs(input, supplement::NamedTuple)
    required = Set((:basis_by_center, :lmax))
    optional = Set((:uncontracted, :width_filtering, :basisfile))
    supplied = Set(Symbol.(keys(supplement)))
    required ⊆ supplied && supplied ⊆ union(required, optional) ||
        throw(ArgumentError("supplement has missing or unsupported keys"))
    supplement.basis_by_center isa AbstractVector ||
        throw(ArgumentError("supplement.basis_by_center must be an AbstractVector"))
    all(label -> label isa AbstractString, supplement.basis_by_center) ||
        throw(ArgumentError("supplement basis labels must be strings"))
    basis_by_center = String.(supplement.basis_by_center)
    length(basis_by_center) == length(input.symbols) ||
        throw(ArgumentError("supplement basis count must match centers"))
    all(==(first(basis_by_center)), basis_by_center) ||
        throw(ArgumentError("heteronuclear basis labels are unsupported"))
    supplement.lmax isa Integer && !(supplement.lmax isa Bool) ||
        throw(ArgumentError("supplement.lmax must be an integer"))
    lmax = Int(supplement.lmax)
    0 <= lmax <= 6 || throw(ArgumentError("supplement.lmax must be between 0 and 6"))
    uncontracted = haskey(supplement, :uncontracted) ? supplement.uncontracted : false
    uncontracted isa Bool || throw(ArgumentError("supplement.uncontracted must be Bool"))
    width_filtering = haskey(supplement, :width_filtering) ? supplement.width_filtering : nothing
    max_width = nothing
    if !isnothing(width_filtering)
        width_filtering isa NamedTuple &&
            Set(Symbol.(keys(width_filtering))) == Set((:max_width,)) ||
            throw(ArgumentError("supplement.width_filtering must be nothing or (; max_width)"))
        max_width = _cartesian_base_positive(
            width_filtering.max_width, "supplement.width_filtering.max_width")
    end
    basisfile = haskey(supplement, :basisfile) ? supplement.basisfile : nothing
    isnothing(basisfile) || basisfile isa AbstractString ||
        throw(ArgumentError("supplement.basisfile must be nothing or an AbstractString"))
    return (; basis_by_center, lmax, uncontracted, width_filtering, max_width,
        basisfile = isnothing(basisfile) ? nothing : String(basisfile))
end

function _cartesian_r3_h2_validation_fixture(input, supplement_input)
    return input.symbols == ["H", "H"] && input.locations == [(0.0, 0.0, -2.0), (0.0, 0.0, 2.0)] &&
           input.q == 5 && input.core_spacing == 0.5 &&
           input.xmax_parallel == 6.0 && input.xmax_transverse == 4.0 &&
           supplement_input.basis_by_center == ["cc-pVTZ", "cc-pVTZ"] &&
           supplement_input.lmax == 1 && !supplement_input.uncontracted &&
           isnothing(supplement_input.width_filtering)
end

function cartesian_residual_gto_mwg_hamiltonian(
    system::NamedTuple;
    basis::NamedTuple,
    supplement::NamedTuple,
    hamfile::Union{Nothing,AbstractString} = nothing,
)::CartesianIDAHamiltonian{Float64}
    !(!isnothing(hamfile) && isempty(String(hamfile))) ||
        throw(ArgumentError("hamfile must not be empty"))
    base = cartesian_base_working_basis(system; basis, supplemented = true)
    base_products = cartesian_base_products(base)
    base_unit_nuclear = cartesian_base_unit_nuclear(base)
    base_vee = cartesian_base_vee(base)
    base_ham = cartesian_base_hamiltonian_assembly(base, base_products, base_unit_nuclear, base_vee)
    supplement_basis = cartesian_residual_gto_supplement_basis(base, supplement)
    residual = cartesian_residual_gto_augmentation(base, supplement_basis)
    augmented_products = cartesian_residual_gto_augmented_products(base, supplement_basis, residual)
    augmented_unit_nuclear = cartesian_residual_gto_augmented_unit_nuclear(base, residual, augmented_products)
    augmented_vee = cartesian_residual_gto_augmented_vee(base, base_ham, residual, augmented_products, augmented_unit_nuclear)
    return cartesian_residual_gto_mwg_hamiltonian_assembly(base, base_ham, supplement_basis,
        residual, augmented_products, augmented_unit_nuclear, augmented_vee; hamfile)
end

function cartesian_residual_gto_supplement_basis(base, supplement::NamedTuple)
    input = base.input
    supplement_input = _cartesian_r3_supplement_inputs(input, supplement)
    raw_supplement = legacy_bond_aligned_diatomic_gaussian_supplement(
        first(input.symbols), first(supplement_input.basis_by_center), input.locations;
        lmax = supplement_input.lmax,
        basisfile = supplement_input.basisfile,
        uncontracted = supplement_input.uncontracted,
        max_width = supplement_input.max_width)
    return (; input = supplement_input, basis = basis_representation(raw_supplement))
end

function cartesian_residual_gto_augmentation(base, supplement_basis)
    C = CartesianFinalBasisRealization
    input = base.input
    return C.pqs_terminal_residual_gto_augmentation(
        base.terminal_basis, base.parent.parent_axis_bundle_object,
        supplement_basis.basis, input.locations)
end

function cartesian_residual_gto_augmented_operators(base, supplement_basis, residual)
    products = cartesian_residual_gto_augmented_products(base, supplement_basis, residual)
    unit_nuclear = cartesian_residual_gto_augmented_unit_nuclear(base, residual, products)
    return _cartesian_residual_augmented_operator_pack(products, unit_nuclear)
end

function cartesian_residual_gto_augmented_products(base, supplement_basis, residual)
    C = CartesianFinalBasisRealization
    input = base.input
    parent = base.parent
    return C.pqs_terminal_residual_gto_augmented_products(
        base.terminal_basis, parent.parent_axis_bundle_object,
        parent.parent_basis_object, supplement_basis.basis, residual,
        input.locations, input.charges)
end

function cartesian_residual_gto_augmented_unit_nuclear(base, residual, augmented_products)
    C = CartesianFinalBasisRealization
    input = base.input
    parent = base.parent
    return C.pqs_terminal_residual_gto_augmented_unit_nuclear(
        base.terminal_basis, parent.parent_axis_bundle_object,
        residual, input.locations, input.charges, augmented_products)
end

_cartesian_residual_augmented_operator_pack(products, unit_nuclear) =
    (; kinetic = products.kinetic, nuclear_attraction_unit_by_center = unit_nuclear,
        position = products.position, x2 = products.x2)

function cartesian_residual_gto_augmented_vee(base, base_ham::CartesianIDAHamiltonian{Float64}, residual, augmented_products, augmented_unit_nuclear)
    C = CartesianFinalBasisRealization
    return C.pqs_terminal_residual_gto_augmented_vee(
        base_ham, base.terminal_basis, base.parent.parent_axis_bundle_object,
        residual, _cartesian_residual_augmented_operator_pack(augmented_products, augmented_unit_nuclear))
end

function cartesian_residual_gto_mwg_hamiltonian_assembly(base, base_ham::CartesianIDAHamiltonian{Float64}, supplement_basis, residual, augmented_products, augmented_unit_nuclear, augmented_vee; hamfile::Union{Nothing,AbstractString} = nothing)::CartesianIDAHamiltonian{Float64}
    !(!isnothing(hamfile) && isempty(String(hamfile))) ||
        throw(ArgumentError("hamfile must not be empty"))
    C = CartesianFinalBasisRealization
    input = base.input
    supplement_input = supplement_basis.input
    operators = _cartesian_residual_augmented_operator_pack(augmented_products, augmented_unit_nuclear)
    ham = CartesianIDAHamiltonian(operators.kinetic,
        operators.nuclear_attraction_unit_by_center, augmented_vee, base_ham.nup,
        base_ham.ndn; nuclear_charges = base_ham.nuclear_charges,
        nuclear_positions = base_ham.nuclear_positions)
    h2_fixture = _cartesian_r3_h2_validation_fixture(input, supplement_input)
    isnothing(hamfile) || C.write_pqs_terminal_residual_gto_augmented_hamiltonian(
        String(hamfile), ham, residual;
        basis_by_center = supplement_input.basis_by_center,
        lmax = supplement_input.lmax,
        uncontracted = supplement_input.uncontracted,
        width_filtering = supplement_input.width_filtering,
        validation_check_labels = h2_fixture ?
            (:h2_lowest_augmented_one_body_orbital_ida_self_coulomb,) : (),
        h2_validation_self_coulomb = h2_fixture ? 0.4574265214362075 : nothing)
    return ham
end
