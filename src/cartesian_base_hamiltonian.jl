const _CARTESIAN_BASE_SYSTEM_KEYS = Set((:atom_symbols, :nuclear_charges, :atom_locations, :nup, :ndn))
const _CARTESIAN_BASE_SIZE_KEYS = Set((:ns, :q))
const _CARTESIAN_BASE_H_BASIS_REQUIRED_KEYS = Set((:core_spacing, :radius))
const _CARTESIAN_BASE_H2_BASIS_REQUIRED_KEYS = Set((:core_spacing, :xmax_parallel, :xmax_transverse))
const _CARTESIAN_BASE_OPTIONAL_BASIS_KEYS = Set((:parent_axis_family, :reference_spacing, :tail_spacing, :nesting, :source_span))
const _CARTESIAN_BASE_H_OPTIONAL_BASIS_KEYS = union(_CARTESIAN_BASE_OPTIONAL_BASIS_KEYS, Set((:d,)))
_cartesian_base_check_keys(input, expected, label) =
    Set(Symbol.(keys(input))) == expected ||
        throw(ArgumentError("$(label) has unsupported keys"))

function _cartesian_base_check_basis_keys(basis, required, optional = _CARTESIAN_BASE_OPTIONAL_BASIS_KEYS)
    supplied = Set(Symbol.(keys(basis)))
    required ⊆ supplied && supplied ⊆ union(required, optional, _CARTESIAN_BASE_SIZE_KEYS) &&
        !isempty(intersect(supplied, _CARTESIAN_BASE_SIZE_KEYS)) ||
        throw(ArgumentError("basis has missing or unsupported keys"))
end

_cartesian_base_get_positive(basis, key, default) =
    haskey(basis, key) ? _cartesian_base_positive(getproperty(basis, key), "basis.$(key)") : default

function _cartesian_base_parent_axis_family(basis)
    family = haskey(basis, :parent_axis_family) ? basis.parent_axis_family : :G10
    family === :G10 || throw(ArgumentError("only parent_axis_family=:G10 is supported"))
    return family
end

function _cartesian_base_nesting(basis)
    value = haskey(basis, :nesting) ? basis.nesting : :pqs
    (value isa Symbol || value isa AbstractString) ||
        throw(ArgumentError("basis.nesting must be :pqs or :wl"))
    nesting = Symbol(value)
    nesting in (:pqs, :wl) || throw(ArgumentError("basis.nesting must be :pqs or :wl"))
    return nesting
end

function _cartesian_base_source_span(basis, nesting)
    value = haskey(basis, :source_span) ? basis.source_span : :ordinary
    (value isa Symbol || value isa AbstractString) ||
        throw(ArgumentError("basis.source_span must be :ordinary or :mapped_comx"))
    source_span = Symbol(value)
    source_span in (:ordinary, :mapped_comx) ||
        throw(ArgumentError("basis.source_span must be :ordinary or :mapped_comx"))
    source_span === :mapped_comx && nesting === :wl &&
        throw(ArgumentError("basis.source_span=:mapped_comx is supported only with nesting=:pqs"))
    return source_span
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

function _cartesian_base_ns(value)
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("basis.ns must be a positive integer"))
    ns = Int(value)
    ns > 0 || throw(ArgumentError("basis.ns must be positive"))
    return ns
end

function _cartesian_base_size_parts(basis, nesting)
    has_ns = haskey(basis, :ns)
    has_q = haskey(basis, :q)
    has_ns || has_q || throw(ArgumentError("basis has missing or unsupported keys"))
    q_rule = nesting === :pqs ? :pqs_ns_equals_q : :wl_ns_minus_2
    if has_ns
        ns = _cartesian_base_ns(basis.ns)
        nesting === :wl && ns < 3 &&
            throw(ArgumentError("basis.ns must be at least 3 for nesting=:wl"))
        q = nesting === :pqs ? ns : ns - 2
        has_q && _cartesian_base_q(basis.q) != q &&
            throw(ArgumentError("basis.q is inconsistent with basis.ns and basis.nesting"))
        return (; ns, q, q_rule, ns_source = :public_ns)
    end
    q = _cartesian_base_q(basis.q)
    ns = nesting === :pqs ? q : q + 2
    return (; ns, q, q_rule, ns_source = :legacy_q_compatibility)
end

function _cartesian_base_integer_charge(charge, label)
    electron_count = round(Int, charge)
    isapprox(charge, electron_count; atol = 1.0e-12, rtol = 0.0) ||
        throw(ArgumentError("$(label) requires integer-valued nuclear charge"))
    return electron_count
end

function _cartesian_base_diatomic_basis_parts(basis)
    _cartesian_base_check_basis_keys(basis, _CARTESIAN_BASE_H2_BASIS_REQUIRED_KEYS)
    nesting = _cartesian_base_nesting(basis)
    size_parts = _cartesian_base_size_parts(basis, nesting)
    nesting === :wl && size_parts.ns < 4 &&
        throw(ArgumentError("basis.ns must be at least 4 for z-axis diatomic nesting=:wl"))
    source_span = _cartesian_base_source_span(basis, nesting)
    return (; size_parts...,
        core_spacing = _cartesian_base_positive(basis.core_spacing, "basis.core_spacing"),
        radius = nothing, d = nothing,
        xmax_parallel = _cartesian_base_positive(basis.xmax_parallel, "basis.xmax_parallel"),
        xmax_transverse = _cartesian_base_positive(basis.xmax_transverse, "basis.xmax_transverse"),
        nesting, source_span,
        parent_axis_family = _cartesian_base_parent_axis_family(basis),
        reference_spacing = _cartesian_base_get_positive(basis, :reference_spacing, 1.0),
        tail_spacing = _cartesian_base_get_positive(basis, :tail_spacing, 10.0))
end

function _cartesian_base_location(value)
    value isa Tuple || throw(ArgumentError("nuclear center must be a fixed 3-tuple"))
    length(value) == 3 || throw(ArgumentError("nuclear center must have three coordinates"))
    location = ntuple(axis -> Float64(value[axis]), 3)
    all(isfinite, location) || throw(ArgumentError("nuclear coordinates must be finite"))
    return location
end

function _cartesian_base_system_parts(system::NamedTuple)
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
    return (; symbols, charges, locations, nup, ndn)
end

function _cartesian_base_inputs(system::NamedTuple, basis::NamedTuple)
    base = _cartesian_base_system_parts(system)
    symbols, charges, locations, nup, ndn =
        base.symbols, base.charges, base.locations, base.nup, base.ndn
    if length(symbols) == 1
        _cartesian_base_check_basis_keys(
            basis, _CARTESIAN_BASE_H_BASIS_REQUIRED_KEYS, _CARTESIAN_BASE_H_OPTIONAL_BASIS_KEYS)
        locations[1] == (0.0, 0.0, 0.0) ||
            throw(ArgumentError("one-center atom must be centered at (0,0,0)"))
        electron_count = _cartesian_base_integer_charge(only(charges), "one-center atom")
        nup + ndn == electron_count ||
            throw(ArgumentError("one-center atom requires neutral all-electron count"))
        nesting = _cartesian_base_nesting(basis)
        size_parts = _cartesian_base_size_parts(basis, nesting)
        source_span = _cartesian_base_source_span(basis, nesting)
        core_spacing = _cartesian_base_positive(basis.core_spacing, "basis.core_spacing")
        return merge(base, (; kind = :h, size_parts...,
            core_spacing,
            radius = _cartesian_base_positive(basis.radius, "basis.radius"),
            d = _cartesian_base_atom_mapping_d(basis, core_spacing), xmax_parallel = nothing,
            xmax_transverse = nothing, nesting, source_span,
            parent_axis_family = _cartesian_base_parent_axis_family(basis),
            reference_spacing = _cartesian_base_get_positive(basis, :reference_spacing, 1.0),
            tail_spacing = _cartesian_base_get_positive(basis, :tail_spacing, 10.0)))
    elseif length(symbols) == 2
        basis_parts = _cartesian_base_diatomic_basis_parts(basis)
        symbols[1] == symbols[2] ||
            throw(ArgumentError("base diatomic requires equal atom-symbol labels"))
        charges[1] == charges[2] ||
            throw(ArgumentError("homonuclear base diatomic requires equal nuclear charges"))
        electron_count = 2 * _cartesian_base_integer_charge(charges[1], "base diatomic")
        all(location -> location[1] == 0.0 && location[2] == 0.0, locations) ||
            throw(ArgumentError("base diatomic centers must lie on the Cartesian z axis"))
        locations[1][3] != locations[2][3] ||
            throw(ArgumentError("base diatomic centers must be distinct"))
        nup + ndn == electron_count ||
            throw(ArgumentError("base diatomic requires neutral all-electron count"))
        return merge(base, (; kind = :h2), basis_parts)
    end
    throw(ArgumentError("only one-center atoms and z-axis homonuclear diatomics are supported"))
end

function _cartesian_base_route(kind, nesting, source_span)
    route_shape = kind === :h ? (:pqs_left, :product, :pqs_right) :
        (:atom_contact_core, :shared_shell_1, :shared_shell_2)
    route_kind = kind === :h ? :one_center_fixed_q_complete_core_shell :
        :bond_aligned_diatomic_independent_pqs_source_box_core_shell
    common = (; route_kind, terms = (:overlap,),
        pair_factor_normalization = :density_normalized, supplement_policy = nothing,
        run_final_basis = false, run_h1 = false, run_h1_j = false)
    nesting === :pqs && return merge(common, (;
        route_family = :pqs_source_box, route_shape,
        source_span,
        product_body_rule = :centered_single_z_slab,
        pqs_retained_rule = :boundary_comx_product_mode_selection,
        product_retained_rule = :product_doside_retained_unit))
    nesting === :wl && return merge(common, (;
        route_family = :white_lindsey_low_order,
        white_lindsey_route_shape = (:standard_cartesian_units, :low_order_comx_coarsening),
        white_lindsey_mapping_rule = :standard_unit_backbone_mapping_family,
        white_lindsey_nesting_rule = :unit_box_low_order_comx_coarsening,
        white_lindsey_retained_rule = :low_order_unit_comx_retained_basis,
        white_lindsey_operator_rule = :low_order_unit_operator_blocks))
    throw(ArgumentError("basis.nesting must be :pqs or :wl"))
end

function _cartesian_base_direct_core_side(ns)
    side = isodd(ns) ? ns : ns + 1
    side > 0 || throw(ArgumentError("derived public ns must be positive"))
    return side
end

function _cartesian_base_atom_parent_axis_counts(input)
    mapping = white_lindsey_atomic_mapping(;
        Z = only(input.charges),
        d = input.core_spacing,
        tail_spacing = input.tail_spacing)
    count = _qwrg_mapped_odd_count_for_extent(
        mapping,
        input.radius;
        reference_spacing = input.reference_spacing)
    side = max(count, _cartesian_base_direct_core_side(input.ns))
    return (x = side, y = side, z = side)
end

function _cartesian_base_stages(input)
    atom_symbols = Tuple(Symbol.(input.symbols))
    nuclear_charges = Tuple(input.charges)
    atom_locations = Tuple(input.locations)
    spacing_inputs = (;
        q = input.q, n_s = input.ns, reference_spacing = input.reference_spacing, tail_spacing = input.tail_spacing,
        q_to_core_spacing_rule = input.q_rule, core_spacing = input.core_spacing,
        xmax_parallel = input.xmax_parallel, xmax_transverse = input.xmax_transverse)
    parent_inputs = (;
        parent_axis_bundle_backend = :pgdg_localized_experimental, parent_axis_family = input.parent_axis_family,
        parent_mapping_tail_spacing = input.tail_spacing,
        parent_mapping_rule = input.kind === :h ? :white_lindsey_atomic_mapping : :identity_mapping,
        parent_mapping_Z = input.kind === :h ? only(input.charges) : nothing,
        parent_mapping_d = input.kind === :h ? input.core_spacing : nothing)
    if input.kind === :h
        system_inputs = (; atom_symbols, nuclear_charges, atom_locations,
            nup = input.nup, ndn = input.ndn, bond_axis = nothing, bond_length = nothing,
            radius = input.radius,
            parent_axis_counts = _cartesian_base_atom_parent_axis_counts(input),
            map_backend = :pgdg_localized_experimental)
    else
        z = (input.locations[1][3], input.locations[2][3])
        system_inputs = (; atom_symbols, nuclear_charges, atom_locations,
            nup = input.nup, ndn = input.ndn, bond_axis = :z, bond_length = abs(z[2] - z[1]),
            radius = max(input.xmax_parallel, input.xmax_transverse),
            parent_axis_counts = nothing, map_backend = :pgdg_localized_experimental)
    end
    system = cartesian_system(system_inputs)
    recipe = cartesian_recipe(_cartesian_base_route(input.kind, input.nesting, input.source_span))
    parent = cartesian_parent(system, spacing_inputs, parent_inputs, recipe)
    shells = cartesian_shells(parent, spacing_inputs, recipe)
    units = cartesian_units(parent, shells, recipe)
    transforms = cartesian_transforms(units, recipe)
    return parent, transforms
end

function _cartesian_record_source_mode_facts(record, contract)
    support = record.support_record
    if record.transform_kind === :direct_identity_transform_contract
        dims = Tuple(Int.(length.(support.outer_box)))::NTuple{3,Int}
        return (; dims, modes = nothing,
            ordering = :x_major_y_major_z_fast, lowdin = false,
            construction_kind = support.lowering_contract_kind,
            axis_intervals = support.outer_box, seed_modes = nothing)
    elseif record.transform_kind === :pqs_source_modes_boundary_selection_shell_realization_contract
        isnothing(contract) && return nothing
        raw_plan = get(contract.metadata, :raw_product_source_plan, nothing)
        dims = isnothing(raw_plan) ?
            Tuple(Int.(support.source_mode_shape))::NTuple{3,Int} :
            raw_plan.source_mode_dims
        ordering = isnothing(raw_plan) ? :x_major_y_major_z_fast :
            raw_plan.source_mode_ordering
        modes = isnothing(raw_plan) ? nothing : raw_plan.source_mode_indices
        retained_rule = get(contract.metadata, :raw_product_source_retained_rule, nothing)
        seed_modes = isnothing(retained_rule) ? nothing : retained_rule.retained_mode_indices
        return (; dims, modes, ordering, lowdin = true,
            construction_kind = :pqs_boundary_comx_source_box,
            axis_intervals = support.outer_box, seed_modes)
    end
    return nothing
end

function _cartesian_append_source_modes!(mode_shells, mode_indices, mode_units, mode_labels,
    local_axes, mode_lowdin, id, unit_key, source_label, facts)
    function append_mode(mode_index, mode)
        push!(mode_shells, id); push!(mode_indices, mode_index); push!(mode_units, unit_key)
        push!(mode_labels, string(source_label, ":", mode[1], "_", mode[2], "_", mode[3]))
        push!(local_axes, mode); push!(mode_lowdin, facts.lowdin)
    end
    if isnothing(facts.modes)
        mode_index = 0
        for ix in 1:facts.dims[1], iy in 1:facts.dims[2], iz in 1:facts.dims[3]
            mode_index += 1
            append_mode(mode_index, (ix, iy, iz))
        end
    else
        for (mode_index, mode) in pairs(facts.modes)
            append_mode(mode_index, mode)
        end
    end
    return nothing
end

function _cartesian_append_seed_relations!(rel_cols, rel_shells, rel_labels, rel_axes,
    block, id, source_label, facts)
    isnothing(facts.seed_modes) && return nothing
    length(facts.seed_modes) == length(block.column_range) ||
        throw(DimensionMismatch("retained seed mode count does not match terminal columns"))
    for (local_col, col) in enumerate(block.column_range)
        mode = facts.seed_modes[local_col]
        push!(rel_cols, col); push!(rel_shells, id)
        push!(rel_labels, string(source_label, ":", mode[1], "_", mode[2], "_", mode[3]))
        push!(rel_axes, mode)
    end
    return nothing
end

function _cartesian_source_mode_provenance(transforms)
    low = get(transforms, :low_order_transforms, nothing)
    basis = get(transforms, :terminal_basis_realization, nothing)
    (isnothing(low) || isnothing(basis)) && return nothing
    plan = low.terminal_retained_rule_plan
    plan.status === :available_terminal_retained_rule_plan || return nothing
    contract_plan = low.retained_unit_transform_contract_plan
    contracts = Dict(contract.unit_key => contract for contract in contract_plan.contracts)
    blocks = Dict(block.unit_key => block for block in basis.blocks)
    shell_id = Int[]; unit_labels = Symbol[]; unit_kinds = Symbol[]
    starts = Int[]; stops = Int[]; labels = Symbol[]; kinds = Symbol[]
    axis_starts = Vector{NTuple{3,Int}}(); axis_stops = Vector{NTuple{3,Int}}()
    dims_rows = Vector{NTuple{3,Int}}(); mode_counts = Int[]; orderings = Symbol[]
    lowdin = Bool[]; shell_by_unit = Dict{Symbol,Tuple{Int,Symbol}}()
    mode_shells = Int[]; mode_indices = Int[]; mode_units = Symbol[]
    mode_labels = String[]; local_axes = Vector{NTuple{3,Int}}(); mode_lowdin = Bool[]
    rel_cols = Int[]; rel_shells = Int[]; rel_labels = String[]
    rel_axes = Vector{NTuple{3,Int}}()
    for record in plan.records
        block = get(blocks, record.support_record.unit_key, nothing)
        isnothing(block) && continue
        contract = get(contracts, record.transform_contract_unit_key, nothing)
        facts = _cartesian_record_source_mode_facts(record, contract)
        isnothing(facts) && continue
        id = length(shell_id) + 1
        source_label = Symbol(string(record.support_record.unit_key), "_source_shell")
        push!(shell_id, id); push!(unit_labels, record.support_record.unit_key)
        push!(unit_kinds, record.lowering_contract_kind)
        push!(starts, first(block.column_range)); push!(stops, last(block.column_range))
        push!(labels, source_label); push!(kinds, facts.construction_kind)
        push!(axis_starts, ntuple(axis -> first(facts.axis_intervals[axis]), 3))
        push!(axis_stops, ntuple(axis -> last(facts.axis_intervals[axis]), 3))
        push!(dims_rows, facts.dims)
        push!(mode_counts, isnothing(facts.modes) ? prod(facts.dims) : length(facts.modes))
        push!(orderings, facts.ordering); push!(lowdin, facts.lowdin)
        shell_by_unit[record.support_record.unit_key] = (id, source_label)
        _cartesian_append_source_modes!(mode_shells, mode_indices, mode_units, mode_labels,
            local_axes, mode_lowdin, id, record.support_record.unit_key, source_label, facts)
        _cartesian_append_seed_relations!(
            rel_cols, rel_shells, rel_labels, rel_axes, block, id, source_label, facts)
    end
    isempty(shell_id) && return nothing
    n_shells = length(shell_id); n_modes = length(mode_shells)
    shell_matrix(rows) = [rows[row][axis] for row in eachindex(rows), axis in 1:3]
    n_relations = length(rel_cols)
    return (;
        shell_by_unit,
        source_shells = (;
            status = :native_terminal_source_shells, schema_version = 1,
            row_count = n_shells, source_shell_id = shell_id,
            unit_label = unit_labels, unit_kind = unit_kinds,
            final_basis_start = starts, final_basis_stop = stops,
            source_shell_label = labels, construction_kind = kinds,
            axis_start = shell_matrix(axis_starts), axis_stop = shell_matrix(axis_stops),
            contracted_dims = shell_matrix(dims_rows), source_mode_count = mode_counts,
            source_mode_ordering = orderings, center_definition = fill(:unavailable, n_shells),
            center_status = fill(:unavailable, n_shells), lowdin_correction_applied = lowdin,
            shell_label_status = fill(:native, n_shells),
            ray_label_status = fill(:unavailable, n_shells),
            radial_order_status = fill(:unavailable, n_shells),
            inferred_from_centers = falses(n_shells), inferred_from_nearest_grid = falses(n_shells),
            inferred_from_support_order = falses(n_shells), inferred_from_support_indices = falses(n_shells),
            inferred_from_raw_to_final_support = falses(n_shells)),
        source_modes = (;
            status = :native_terminal_source_modes, schema_version = 1,
            row_count = n_modes, source_shell_id = mode_shells, mode_index = mode_indices,
            unit_label = mode_units, native_source_id_label = mode_labels,
            local_axis_x = [axis[1] for axis in local_axes],
            local_axis_y = [axis[2] for axis in local_axes],
            local_axis_z = [axis[3] for axis in local_axes],
            center_x = fill(NaN, n_modes), center_y = fill(NaN, n_modes),
            center_z = fill(NaN, n_modes), center_definition = fill(:unavailable, n_modes),
            center_status = fill(:unavailable, n_modes), lowdin_correction_applied = mode_lowdin,
            source_mode_status = fill(:native, n_modes), shell_label_status = fill(:native, n_modes),
            ray_label_status = fill(:unavailable, n_modes),
            radial_order_status = fill(:unavailable, n_modes),
            inferred_from_centers = falses(n_modes), inferred_from_nearest_grid = falses(n_modes),
            inferred_from_support_order = falses(n_modes), inferred_from_support_indices = falses(n_modes),
            inferred_from_raw_to_final_support = falses(n_modes)),
        final_basis_source_relations = (;
            final_basis_col = rel_cols, relation_index = ones(Int, n_relations),
            relation_kind = fill(:boundary_mode, n_relations),
            source_shell_id = rel_shells, source_mode_label = rel_labels,
            local_axis_x = [axis[1] for axis in rel_axes],
            local_axis_y = [axis[2] for axis in rel_axes],
            local_axis_z = [axis[3] for axis in rel_axes],
            relation_status = fill(:native_retained_boundary_seed, n_relations),
            shell_label_status = fill(:native, n_relations),
            ray_label_status = fill(:unavailable, n_relations),
            radial_order_status = fill(:unavailable, n_relations),
            coefficient_status = fill(:not_serialized, n_relations),
            weight_status = fill(:not_serialized, n_relations),
            span_status = fill(:unavailable, n_relations),
            inferred_from_centers = falses(n_relations),
            inferred_from_nearest_grid = falses(n_relations),
            inferred_from_support_order = falses(n_relations),
            inferred_from_support_indices = falses(n_relations),
            inferred_from_raw_to_final_support = falses(n_relations)))
end

function _cartesian_inventory_axis_edge(values, index, side)
    n = length(values)
    n == 1 && return Float64(only(values))
    spacing = index == 1 ? abs(values[2] - values[1]) :
        index == n ? abs(values[n] - values[n - 1]) :
        (abs(values[index] - values[index - 1]) +
            abs(values[index + 1] - values[index])) / 2
    return Float64(values[index] + (side === :low ? -spacing / 2 : spacing / 2))
end

function _cartesian_inventory_geometry(block, centers)
    box = ntuple(axis -> begin
        vals = (state[axis] for state in block.support_states)
        minimum(vals):maximum(state[axis] for state in block.support_states)
    end, 3)
    axis_values = (centers.x, centers.y, centers.z)
    physical = ntuple(axis -> (
        _cartesian_inventory_axis_edge(axis_values[axis], first(box[axis]), :low),
        _cartesian_inventory_axis_edge(axis_values[axis], last(box[axis]), :high),
    ), 3)
    return (; index_ranges = (; x = box[1], y = box[2], z = box[3]),
        physical_ranges = (; x = physical[1], y = physical[2], z = physical[3]))
end

function _cartesian_terminal_inventory_rows(transforms, bundles)
    low = get(transforms, :low_order_transforms, nothing)
    basis = get(transforms, :terminal_basis_realization, nothing)
    (isnothing(low) || isnothing(basis)) && return NamedTuple[]
    contract_plan = low.retained_unit_transform_contract_plan
    isnothing(contract_plan) && return NamedTuple[]
    units = CartesianRetainedUnits.units(contract_plan.retained_unit_plan)
    contracts = Dict(contract.unit_key => contract for contract in
        CartesianRetainedUnitTransformContracts.transform_contracts(contract_plan))
    blocks = Dict(block.unit_key => block for block in basis.blocks)
    block_key_by_retained = Dict{Symbol,Symbol}()
    plan = low.terminal_retained_rule_plan
    if !isnothing(plan) && hasproperty(plan, :records)
        for record in plan.records
            isnothing(record.retained_unit_key) && continue
            block_key_by_retained[record.retained_unit_key] = record.support_record.unit_key
        end
    end
    rows = NamedTuple[]
    key_index = Dict{Tuple,Int}()
    centers = _cartesian_manifest_axis_centers(bundles)
    for unit in units
        contract = get(contracts, unit.unit_key, nothing)
        block = get(blocks, get(block_key_by_retained, unit.unit_key, unit.unit_key), nothing)
        (isnothing(contract) || isnothing(block)) && continue
        meta = contract.metadata
        geometry = _cartesian_inventory_geometry(block, centers)
        shell_index = get(meta, :terminal_region_shell_index, :unavailable)
        key = (
            unit.terminal_region_key, unit.terminal_region_kind, unit.lowering_kind,
            unit.unit_kind, contract.transform_path,
            isnothing(block.coefficients) ? :identity : :compact_product,
            shell_index, geometry.index_ranges, geometry.physical_ranges,
            get(meta, :slab_normal_axis, :unavailable),
            get(meta, :slab_side, :unavailable),
            Int(get(meta, :slab_thickness, 0)),
            Int(get(meta, :slab_stack_index, 0)),
            Int(get(meta, :slab_stack_count, 0)),
        )
        support_rows = length(block.support_indices)
        final_cols = length(block.column_range)
        if haskey(key_index, key)
            old = rows[key_index[key]]
            rows[key_index[key]] = merge(old, (;
                support_rows = old.support_rows + support_rows,
                final_cols = old.final_cols + final_cols,
                compression_ratio = (old.support_rows + support_rows) /
                    (old.final_cols + final_cols),
            ))
        else
            key_index[key] = length(rows) + 1
            push!(rows, (;
                region_key = key[1], region_kind = key[2], lowering_kind = key[3],
                retained_unit_kind = key[4], realization_kind = key[5],
                realization_class = key[6], shell_index = key[7],
                index_ranges = key[8], physical_ranges = key[9], support_rows, final_cols,
                compression_ratio = support_rows / final_cols,
                slab_axis = key[10], slab_side = key[11], slab_thickness = key[12],
                slab_stack_index = key[13], slab_stack_count = key[14],
            ))
        end
    end
    return rows
end

_cartesian_base_atom_location_matrix(input) =
    [input.locations[row][axis] for row in eachindex(input.locations), axis in 1:3]

_cartesian_base_axis_counts_tuple(axis_counts) =
    (Int(axis_counts.x), Int(axis_counts.y), Int(axis_counts.z))

function _cartesian_write_sidecar_values!(file, group, values)
    for key in keys(values); file["$(group)/$(key)"] = getproperty(values, key); end
end

function _cartesian_manifest_label_table(n)
    return (;
        final_basis_col = collect(1:n), sector = fill(:base, n),
        unit_label = fill(:unavailable, n), unit_kind = fill(:unavailable, n),
        source_region_label = fill(:unavailable, n), source_region_label_status = fill(:unavailable, n),
        source_box_label = fill(:unavailable, n), source_box_label_status = fill(:unavailable, n),
        owner_nucleus_index = zeros(Int, n), owner_label_status = fill(:unavailable, n),
        shell_label_status = fill(:unavailable, n), shell_index = zeros(Int, n),
        ray_label_status = fill(:unavailable, n), ray_id = zeros(Int, n),
        ray_family_label = fill(:unavailable, n), radial_order_status = fill(:unavailable, n),
        radial_order = zeros(Int, n),
        center_x = zeros(Float64, n), center_y = zeros(Float64, n), center_z = zeros(Float64, n),
        center_definition = fill(:unavailable, n), center_status = fill(:unavailable, n),
        lowdin_correction_applied = falses(n), supplement_label = fill("unavailable", n),
        angular_power_x = fill(-1, n), angular_power_y = fill(-1, n), angular_power_z = fill(-1, n),
        inferred_from_centers = falses(n), inferred_from_nearest_grid = falses(n),
        inferred_from_support_order = falses(n), inferred_from_support_indices = falses(n),
        inferred_from_raw_to_final_support = falses(n),
    )
end

_cartesian_manifest_axis_centers(bundles) =
    (; x = Float64.(_nested_axis_pgdg(bundles, :x).centers),
        y = Float64.(_nested_axis_pgdg(bundles, :y).centers),
        z = Float64.(_nested_axis_pgdg(bundles, :z).centers))

_cartesian_manifest_state_center(centers, state) =
    (centers.x[state[1]], centers.y[state[2]], centers.z[state[3]])

function _cartesian_manifest_block_center(centers, block, local_col)
    isnothing(block.coefficients) &&
        return _cartesian_manifest_state_center(centers, block.support_states[local_col])
    total = cx = cy = cz = 0.0
    for (row, state) in pairs(block.support_states)
        x, y, z = _cartesian_manifest_state_center(centers, state)
        weight = abs2(block.coefficients[row, local_col])
        total += weight
        cx += weight * x
        cy += weight * y
        cz += weight * z
    end
    total > 0.0 || return _cartesian_manifest_state_center(centers, first(block.support_states))
    return (cx / total, cy / total, cz / total)
end

function _cartesian_manifest_fill_base_labels!(labels, base)
    basis = base.terminal_basis
    bundles = base.parent.parent_axis_bundle_object
    centers = _cartesian_manifest_axis_centers(bundles)
    source_mode_provenance =
        hasproperty(base, :source_mode_provenance) ? base.source_mode_provenance : nothing
    shell_by_unit = isnothing(source_mode_provenance) ?
        Dict{Symbol,Tuple{Int,Symbol}}() : source_mode_provenance.shell_by_unit
    for block in basis.blocks, (local_col, col) in enumerate(block.column_range)
        direct = isnothing(block.coefficients)
        x, y, z = _cartesian_manifest_block_center(centers, block, local_col)
        labels.unit_label[col] = block.unit_key
        labels.unit_kind[col] = direct ? :direct_parent_support : :lowdin_support_mixture
        labels.source_region_label[col], labels.source_region_label_status[col] = block.unit_key, :available
        labels.center_x[col], labels.center_y[col], labels.center_z[col] = x, y, z
        labels.center_definition[col] =
            direct ? :parent_support_row_center : :coefficient_abs2_support_centroid
        labels.center_status[col], labels.lowdin_correction_applied[col] =
            direct ? :available : :representative, !direct
        if base.input.kind === :h
            labels.owner_nucleus_index[col], labels.owner_label_status[col] = 1, :available
        end
        shell = get(shell_by_unit, block.unit_key, nothing)
        if !isnothing(shell)
            labels.source_box_label[col], labels.source_box_label_status[col] =
                shell[2], :available
            labels.shell_index[col], labels.shell_label_status[col] =
                shell[1], :native
        end
    end
    return labels
end

function _cartesian_manifest_fill_residual_labels!(labels, residual, augmented_products)
    offset = residual.base_dimension
    position = augmented_products.position
    for row in 1:residual.residual_dimension
        col = offset + row
        owner = residual.residual_source_owner_indices[row]
        labels.sector[col], labels.unit_label[col] = :residual, :residual_gaussian
        labels.unit_kind[col] = :owner_local_residual_gaussian
        labels.source_region_label[col], labels.source_region_label_status[col] =
            Symbol("owner_", string(owner)), :available
        labels.owner_nucleus_index[col], labels.owner_label_status[col] = owner, :available
        labels.center_x[col] = position.x[col, col]
        labels.center_y[col] = position.y[col, col]
        labels.center_z[col] = position.z[col, col]
        labels.center_definition[col], labels.center_status[col] =
            :exact_position_expectation, :representative
        labels.lowdin_correction_applied[col], labels.supplement_label[col] =
            true, residual.residual_labels[row]
    end
    return labels
end

function _cartesian_manifest_final_basis_labels(base; residual = nothing, augmented_products = nothing)
    base_dimension = base.terminal_basis.final_dimension
    residual_dimension = isnothing(residual) ? 0 : residual.residual_dimension
    labels = _cartesian_manifest_label_table(base_dimension + residual_dimension)
    _cartesian_manifest_fill_base_labels!(labels, base)
    isnothing(residual) ||
        _cartesian_manifest_fill_residual_labels!(labels, residual, augmented_products)
    return labels
end

function _cartesian_recipe_provenance(input, parent, base_dimension;
    residual_dimension = 0, supplement_input = nothing, producer, route)
    return (;
        provenance_version = 1, producer, nesting = input.nesting, route,
        ns = input.ns, q = input.q, q_rule = input.q_rule, ns_source = input.ns_source,
        core_spacing = input.core_spacing, padding = input.kind === :h ? input.radius : input.xmax_transverse,
        radius = input.radius, xmax_parallel = input.xmax_parallel,
        xmax_transverse = input.xmax_transverse,
        extent_source = input.kind === :h ? :one_center_radius : :z_axis_diatomic_padding_extents,
        parent_axis_counts = _cartesian_base_axis_counts_tuple(parent.axis_counts),
        atom_symbols = copy(input.symbols), nuclear_charges = copy(input.charges),
        atom_locations = _cartesian_base_atom_location_matrix(input), nup = input.nup, ndn = input.ndn,
        basisname = isnothing(supplement_input) ? nothing : first(supplement_input.basis_by_center),
        basisfile = isnothing(supplement_input) ? nothing : supplement_input.basisfile,
        lmax = isnothing(supplement_input) ? nothing : supplement_input.lmax,
        uncontracted = isnothing(supplement_input) ? nothing : supplement_input.uncontracted,
        width_filtering = isnothing(supplement_input) ? nothing : supplement_input.width_filtering,
        base_dimension, residual_dimension, augmented_dimension = base_dimension + residual_dimension,
    )
end

function _cartesian_write_hamiltonian_manifest(path, base, ham; residual = nothing,
    augmented_products = nothing, supplement_input = nothing, producer, route)
    labels = _cartesian_manifest_final_basis_labels(base; residual, augmented_products)
    size(ham.kinetic, 1) == length(labels.final_basis_col) ||
        throw(DimensionMismatch("Hamiltonian manifest dimension mismatch"))
    recipe = _cartesian_recipe_provenance(base.input, base.parent,
        base.terminal_basis.final_dimension; residual_dimension =
        isnothing(residual) ? 0 : residual.residual_dimension,
        supplement_input, producer, route)
    jldopen(String(path), "r+") do file
        file["hamiltonian_manifest/manifest_version"] = 1
        _cartesian_write_sidecar_values!(file, "hamiltonian_manifest/final_basis_labels", labels)
        if hasproperty(base, :source_mode_provenance) && !isnothing(base.source_mode_provenance)
            _cartesian_write_sidecar_values!(
                file, "hamiltonian_manifest/source_shells",
                base.source_mode_provenance.source_shells)
            _cartesian_write_sidecar_values!(
                file, "hamiltonian_manifest/source_modes",
                base.source_mode_provenance.source_modes)
            relations = base.source_mode_provenance.final_basis_source_relations
            isempty(relations.final_basis_col) || _cartesian_write_sidecar_values!(
                file, "hamiltonian_manifest/final_basis_source_relations", relations)
        end
        _cartesian_write_sidecar_values!(file, "recipe_provenance", recipe)
    end
    return path
end

function _cartesian_base_route_label(input)
    input.kind === :h && input.nesting === :pqs && return :one_center_pqs_base
    input.kind === :h && input.nesting === :wl && return :one_center_wl_base
    input.kind === :h2 && input.nesting === :pqs && return :z_axis_diatomic_pqs_base
    input.kind === :h2 && input.nesting === :wl && return :z_axis_diatomic_wl_base
    throw(ArgumentError("unsupported base route label for kind=$(input.kind), nesting=$(input.nesting)"))
end

function _cartesian_supplemented_route_label(input)
    input.kind === :h && input.nesting === :pqs && return :one_center_pqs_residual_gto_mwg
    input.kind === :h && input.nesting === :wl && return :one_center_wl_residual_gto_mwg
    input.kind === :h2 && input.nesting === :pqs && return :z_axis_diatomic_pqs_residual_gto_mwg
    input.kind === :h2 && input.nesting === :wl && return :z_axis_diatomic_wl_residual_gto_mwg
    throw(ArgumentError("unsupported supplemented route label for kind=$(input.kind), nesting=$(input.nesting)"))
end

function _cartesian_base_write_hamiltonian(path, ham, base)
    input = base.input
    parent = base.parent
    write_cartesian_ida_hamiltonian(path, ham)
    route = _cartesian_base_route_label(input)
    mapping_kind = input.kind === :h ? :white_lindsey_atomic_mapping : :multicenter_pqs_mapping
    values = (; provenance_version = 1, producer = :cartesian_base_hamiltonian,
        nesting = input.nesting, route, ns = input.ns, q = input.q,
        q_rule = input.q_rule, ns_source = input.ns_source,
        core_spacing = input.core_spacing,
        reference_spacing = input.reference_spacing, tail_spacing = input.tail_spacing,
        parent_axis_family = input.parent_axis_family,
        parent_axis_counts = _cartesian_base_axis_counts_tuple(parent.axis_counts),
        mapping_kind, mapping_d = input.d, radius = input.radius,
        xmax_parallel = input.xmax_parallel, xmax_transverse = input.xmax_transverse,
        atom_symbols = copy(input.symbols), nuclear_charges = copy(input.charges),
        atom_locations = _cartesian_base_atom_location_matrix(input),
        nup = input.nup, ndn = input.ndn, final_dimension = size(ham.kinetic, 1))
    jldopen(String(path), "r+") do file
        _cartesian_write_sidecar_values!(file, "producer_provenance", values)
    end
    _cartesian_write_hamiltonian_manifest(path, base, ham; producer = :cartesian_base_hamiltonian,
        route)
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

function _cartesian_base_terminal_basis(base)
    basis = base.terminal_basis
    !isnothing(basis) && return basis
    nesting = hasproperty(base.input, :nesting) ? base.input.nesting : :pqs
    throw(ArgumentError("nesting=$(repr(nesting)) base Hamiltonian path is not yet wired to a terminal basis"))
end

function cartesian_base_working_basis(system::NamedTuple; basis::NamedTuple, supplemented::Bool = false)
    input = supplemented ? _cartesian_supplemented_inputs(system, basis) :
        _cartesian_base_inputs(system, basis)
    parent, transforms = _cartesian_base_stages(input)
    source_mode_provenance = _cartesian_source_mode_provenance(transforms)
    basis_realization = transforms.terminal_basis_realization
    terminal_inventory = isnothing(basis_realization) ? nothing :
        (; final_dimension = basis_realization.final_dimension,
            rows = _cartesian_terminal_inventory_rows(transforms, parent.parent_axis_bundle_object))
    return (; input, parent, terminal_basis = transforms.terminal_basis_realization,
        source_mode_provenance, terminal_inventory)
end

cartesian_base_products(base) =
    _pqs_source_box_route_driver_terminal_products(_cartesian_base_terminal_basis(base), _cartesian_base_pgdg(base))

cartesian_base_unit_nuclear(base) =
    _pqs_source_box_route_driver_terminal_unit_nuclear(_cartesian_base_terminal_basis(base),
        _cartesian_base_expansion(), _cartesian_base_pgdg(base), base.input.locations)

cartesian_base_vee(base) =
    _pqs_source_box_route_driver_terminal_vee(_cartesian_base_terminal_basis(base),
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
        _cartesian_base_write_hamiltonian(String(hamfile), ham, base)
    return ham
end

function _cartesian_r3_diatomic_inputs(system::NamedTuple, basis::NamedTuple)
    basis_parts = _cartesian_base_diatomic_basis_parts(basis)
    base = _cartesian_base_system_parts(system)
    symbols, charges, locations, nup, ndn =
        base.symbols, base.charges, base.locations, base.nup, base.ndn
    length(symbols) == length(charges) == length(locations) == 2 ||
        throw(ArgumentError("R3 usability facade supports two-center systems only"))
    symbols[1] == symbols[2] ||
        throw(ArgumentError("heteronuclear supplements are unsupported"))
    charges[1] == charges[2] ||
        throw(ArgumentError("homonuclear supplemented diatomics require equal nuclear charges"))
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
    return merge(base, (; kind = :h2), basis_parts)
end

function _cartesian_supplemented_inputs(system::NamedTuple, basis::NamedTuple)
    base = _cartesian_base_system_parts(system)
    length(base.symbols) == 1 && return _cartesian_base_inputs(system, basis)
    return _cartesian_r3_diatomic_inputs(system, basis)
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
    augmented_products = cartesian_residual_gto_augmented_products(
        base, supplement_basis, residual; base_kinetic = base_ham.kinetic)
    augmented_unit_nuclear = cartesian_residual_gto_augmented_unit_nuclear(
        base, residual, augmented_products;
        base_unit_nuclear = base_ham.nuclear_attraction_unit_by_center)
    augmented_vee = cartesian_residual_gto_augmented_vee(base, base_ham, residual, augmented_products, augmented_unit_nuclear)
    return cartesian_residual_gto_mwg_hamiltonian_assembly(base, base_ham, supplement_basis,
        residual, augmented_products, augmented_unit_nuclear, augmented_vee; hamfile)
end

function cartesian_residual_gto_supplement_basis(base, supplement::NamedTuple)
    input = base.input
    supplement_input = _cartesian_r3_supplement_inputs(input, supplement)
    raw_supplement = if length(input.symbols) == 1
        legacy_atomic_gaussian_supplement(
            first(input.symbols), first(supplement_input.basis_by_center);
            lmax = supplement_input.lmax,
            basisfile = supplement_input.basisfile,
            uncontracted = supplement_input.uncontracted,
            max_width = supplement_input.max_width)
    else
        legacy_bond_aligned_diatomic_gaussian_supplement(
            first(input.symbols), first(supplement_input.basis_by_center), input.locations;
            lmax = supplement_input.lmax,
            basisfile = supplement_input.basisfile,
            uncontracted = supplement_input.uncontracted,
            max_width = supplement_input.max_width)
    end
    return (; input = supplement_input, basis = basis_representation(raw_supplement))
end

function cartesian_residual_gto_augmentation(base, supplement_basis)
    C = CartesianFinalBasisRealization
    input = base.input
    return C.pqs_terminal_residual_gto_augmentation(
        base.terminal_basis, base.parent.parent_axis_bundle_object,
        supplement_basis.basis, input.locations)
end

function cartesian_residual_gto_augmented_products(base, supplement_basis, residual; base_kinetic = nothing)
    C = CartesianFinalBasisRealization
    input = base.input
    parent = base.parent
    return C.pqs_terminal_residual_gto_augmented_products(
        base.terminal_basis, parent.parent_axis_bundle_object,
        parent.parent_basis_object, supplement_basis.basis, residual,
        input.locations, input.charges; base_kinetic)
end

function cartesian_residual_gto_augmented_unit_nuclear(base, residual, augmented_products; base_unit_nuclear = nothing)
    C = CartesianFinalBasisRealization
    input = base.input
    parent = base.parent
    return C.pqs_terminal_residual_gto_augmented_unit_nuclear(
        base.terminal_basis, parent.parent_axis_bundle_object,
        residual, input.locations, input.charges, augmented_products;
        base_unit_nuclear)
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
    isnothing(hamfile) || _cartesian_write_hamiltonian_manifest(
        String(hamfile), base, ham; residual, augmented_products, supplement_input,
        producer = :cartesian_residual_gto_mwg_hamiltonian,
        route = _cartesian_supplemented_route_label(input))
    return ham
end
