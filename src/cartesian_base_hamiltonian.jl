const _CARTESIAN_BASE_SYSTEM_KEYS = Set((:atom_symbols, :nuclear_charges, :atom_locations, :nup, :ndn))
const _CARTESIAN_BASE_H_BASIS_REQUIRED_KEYS = Set((:q, :core_spacing, :radius, :d))
const _CARTESIAN_BASE_H2_BASIS_REQUIRED_KEYS = Set((:q, :core_spacing, :xmax_parallel, :xmax_transverse))
const _CARTESIAN_BASE_OPTIONAL_BASIS_KEYS = Set((:parent_axis_family, :reference_spacing, :tail_spacing))
_cartesian_base_check_keys(input, expected, label) =
    Set(Symbol.(keys(input))) == expected ||
        throw(ArgumentError("$(label) has unsupported keys"))

function _cartesian_base_check_basis_keys(basis, required)
    supplied = Set(Symbol.(keys(basis)))
    required ⊆ supplied && supplied ⊆ union(required, _CARTESIAN_BASE_OPTIONAL_BASIS_KEYS) ||
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
        _cartesian_base_check_basis_keys(basis, _CARTESIAN_BASE_H_BASIS_REQUIRED_KEYS)
        symbols == ["H"] && charges == [1.0] ||
            throw(ArgumentError("only origin-centered H is supported"))
        locations[1] == (0.0, 0.0, 0.0) ||
            throw(ArgumentError("H must be centered at (0,0,0)"))
        nup + ndn == 1 || throw(ArgumentError("H requires one electron"))
        return merge(base, (; kind = :h, q = _cartesian_base_q(basis.q),
            core_spacing = _cartesian_base_positive(basis.core_spacing, "basis.core_spacing"),
            radius = _cartesian_base_positive(basis.radius, "basis.radius"),
            d = _cartesian_base_positive(basis.d, "basis.d"), xmax_parallel = nothing,
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
        parent_mapping_d = input.kind === :h ? input.d : nothing)
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
    input = _cartesian_base_inputs(system, basis)
    parent, transforms = _cartesian_base_stages(input)
    ham = _cartesian_base_ida_hamiltonian(
        transforms.terminal_basis_realization, parent.parent_axis_bundle_object,
        input.locations, input.charges, input.nup, input.ndn)
    isnothing(hamfile) ||
        _cartesian_base_write_hamiltonian(String(hamfile), ham, input, parent)
    return ham
end
