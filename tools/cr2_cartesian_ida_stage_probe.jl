#!/usr/bin/env julia

# Developer stage probe for the Cr2-facing Cartesian IDA producer route.
# This intentionally uses only the public cartesian_* stage sequence.

using GaussletBases

const ANGSTROM_TO_BOHR = 1.8897261246257702

R = 1.68 * ANGSTROM_TO_BOHR

atom_symbols = ("Cr", "Cr")
nuclear_charges = (24, 24)
atom_locations = ((0.0, 0.0, -R / 2), (0.0, 0.0, R / 2))
nup = 24
ndn = 24
bond_axis = :z
bond_length = R
radius = 4.0
parent_axis_counts = nothing
map_backend = :pgdg_localized_experimental

q = 5
n_s = 5
reference_spacing = 1.0
tail_spacing = 10.0
q_to_core_spacing_rule = :standard_pqs_ns_equals_q
core_spacing = 1.2 / (24 * (q - 3))
xmax_parallel = 20.0
xmax_transverse = 20.0

parent_axis_family = :G10
parent_axis_bundle_backend = :pgdg_localized_experimental
parent_mapping_rule = :identity_mapping
parent_mapping_Z = nothing
parent_mapping_d = nothing
parent_mapping_tail_spacing = tail_spacing

route_family = :pqs_source_box
route_kind = :bond_aligned_diatomic_independent_pqs_source_box_core_shell
route_shape = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
product_body_rule = :centered_single_z_slab
pqs_retained_rule = :boundary_comx_product_mode_selection
product_retained_rule = :product_doside_retained_unit
terms = (:overlap, :position_x, :position_y, :position_z,
    :x2_x, :x2_y, :x2_z, :kinetic)
pair_factor_normalization = :density_normalized

white_lindsey_route_shape = (:standard_cartesian_units, :low_order_comx_coarsening)
white_lindsey_mapping_rule = :standard_unit_backbone_mapping_family
white_lindsey_nesting_rule = :unit_box_low_order_comx_coarsening
white_lindsey_retained_rule = :low_order_unit_comx_retained_basis
white_lindsey_operator_rule = :low_order_unit_operator_blocks

supplement_policy = :mwg_residual_gto
supplement_basis = "cc-pV5Z"
supplement_lmax = 1
supplement_uncontracted = false
materialize_route = true
save_ida_hamiltonian = true
hamiltonian_output = :cartesian_ida_hamiltonian
hamfile = joinpath(tempdir(), "cr2_all_electron_ida_probe.jld2")
run_final_basis = true
run_h1 = true

system_inputs = (; atom_symbols, nuclear_charges, atom_locations,
    nup, ndn, bond_axis, bond_length, radius, parent_axis_counts, map_backend)
spacing_inputs = (; q, n_s, reference_spacing, tail_spacing,
    q_to_core_spacing_rule, core_spacing, xmax_parallel, xmax_transverse)
parent_inputs = (; parent_axis_bundle_backend,
    parent_axis_family, parent_mapping_rule,
    parent_mapping_Z, parent_mapping_d, parent_mapping_tail_spacing)
route_inputs = (; route_family, route_kind, route_shape, product_body_rule,
    pqs_retained_rule, product_retained_rule, terms, pair_factor_normalization,
    white_lindsey_route_shape, white_lindsey_mapping_rule,
    white_lindsey_nesting_rule, white_lindsey_retained_rule,
    white_lindsey_operator_rule,
    supplement_policy, supplement_basis, supplement_lmax,
    supplement_uncontracted, run_final_basis, run_h1)
materialization_inputs = (; materialize_route,
    save_ham_artifact = save_ida_hamiltonian, hamfile,
    hamiltonian_output)

last_successful_public_stage = nothing
first_blocker = nothing

function _probe_get(obj, key, default = nothing)
    isnothing(obj) && return default
    hasproperty(obj, key) && return getproperty(obj, key)
    return default
end

function _probe_dig(obj, keys...)
    value = obj
    for key in keys
        value = _probe_get(value, key, nothing)
        isnothing(value) && return nothing
    end
    return value
end

function _probe_first_notnothing(values...)
    for value in values
        isnothing(value) || return value
    end
    return nothing
end

function _probe_print(label, value)
    isnothing(value) && return nothing
    println("  ", label, " = ", value)
    return nothing
end

function _probe_count_summary(counts)
    isnothing(counts) && return nothing
    counts isa NamedTuple && return join(
        ("$(key)=$(getproperty(counts, key))" for key in keys(counts)),
        ", ",
    )
    return counts
end

function _probe_dense_memory(dim, center_count)
    isnothing(dim) && return nothing
    dim isa Integer || return nothing
    ncenter = isnothing(center_count) ? 2 : Int(center_count)
    matrix_mb = 8 * dim * dim / 1024^2
    return (;
        K_MB = round(matrix_mb; digits = 3),
        U_MB = round(ncenter * matrix_mb; digits = 3),
        Vee_MB = round(matrix_mb; digits = 3),
        total_MB = round((ncenter + 2) * matrix_mb; digits = 3),
    )
end

_probe_axis_index(axis::Symbol) =
    axis === :x ? 1 :
    axis === :y ? 2 :
    axis === :z ? 3 :
    throw(ArgumentError("unsupported bond axis $(repr(axis))"))

_probe_axis_symbol(axis::Int) = (:x, :y, :z)[axis]

function _probe_axis_centers(parent)
    counts = _probe_get(parent, :axis_counts)
    fallback =
        isnothing(counts) ?
        nothing :
        (collect(1:counts.x), collect(1:counts.y), collect(1:counts.z))
    basis = _probe_get(parent, :parent_basis_object)
    isnothing(basis) && return fallback
    axes = GaussletBases.CartesianParentGaussletBases.parent_axes(basis)
    return (
        collect(axes.x.center_data),
        collect(axes.y.center_data),
        collect(axes.z.center_data),
    )
end

function _probe_nearest_index(axis_values, coordinate)
    distances = abs.(Float64.(axis_values) .- Float64(coordinate))
    return findmin(distances)[2]
end

function _probe_axis_window(axis_values, index::Int; radius::Int = 2)
    lo = max(1, index - radius)
    hi = min(length(axis_values), index + radius)
    return Tuple(
        (index = i, center = Float64(axis_values[i]), spacing_left =
            i == 1 ? nothing : Float64(axis_values[i] - axis_values[i - 1]),
            spacing_right =
                i == length(axis_values) ?
                nothing :
                Float64(axis_values[i + 1] - axis_values[i]))
        for i in lo:hi
    )
end

function _probe_core_box(center::NTuple{3,Int}, side::Int)
    radius = div(side - 1, 2)
    return ntuple(axis -> (center[axis] - radius):(center[axis] + radius), 3)
end

_probe_ranges_overlap(a, b) = first(a) <= last(b) && first(b) <= last(a)

function _probe_bond_length(locations)
    length(locations) == 2 || return nothing
    return sqrt(sum(
        (Float64(locations[2][axis]) - Float64(locations[1][axis]))^2
        for axis in 1:3
    ))
end

function _probe_geometry_facts(parent)
    setup = _probe_get(parent, :standard_setup)
    axes = _probe_axis_centers(parent)
    isnothing(setup) && return nothing
    isnothing(axes) && return nothing
    locations = Tuple(_probe_get(parent, :atom_locations, ()))
    length(locations) == 2 || return nothing
    core_side = setup.core_cube_side
    bond_axis = _probe_get(parent, :bond_axis, :z)
    bond_axis_index = _probe_axis_index(bond_axis)
    nuclear_indices = Tuple(
        ntuple(axis -> _probe_nearest_index(axes[axis], location[axis]), 3)
        for location in locations
    )
    atom_axis_indices = Tuple(index[bond_axis_index] for index in nuclear_indices)
    atom_order = Tuple(sortperm(collect(atom_axis_indices)))
    left_atom, right_atom = atom_order
    core_boxes = Tuple(_probe_core_box(index, core_side) for index in nuclear_indices)
    left_core = core_boxes[left_atom]
    right_core = core_boxes[right_atom]
    index_separation =
        abs(nuclear_indices[right_atom][bond_axis_index] -
            nuclear_indices[left_atom][bond_axis_index])
    bond_axis_gap =
        first(right_core[bond_axis_index]) -
        last(left_core[bond_axis_index]) - 1
    transverse_ranges_overlap = Tuple(
        (
            axis = _probe_axis_symbol(axis),
            overlap = _probe_ranges_overlap(left_core[axis], right_core[axis]),
            left_range = left_core[axis],
            right_range = right_core[axis],
        ) for axis in 1:3 if axis != bond_axis_index
    )
    proposed_contact_core = ntuple(
        axis -> begin
            lo = min(first(left_core[axis]), first(right_core[axis]))
            hi = max(last(left_core[axis]), last(right_core[axis]))
            lo:hi
        end,
        3,
    )
    proposed_contact_core_dimensions =
        Tuple(length(proposed_contact_core[axis]) for axis in 1:3)
    coordinate_errors = Tuple(
        ntuple(
            axis ->
                Float64(axes[axis][nuclear_indices[atom][axis]]) -
                Float64(locations[atom][axis]),
            3,
        ) for atom in eachindex(locations)
    )
    bond_axis_windows = Tuple(
        _probe_axis_window(axes[bond_axis_index], index)
        for index in atom_axis_indices
    )
    classification =
        bond_axis_gap < q ?
        :atom_seed_boxes_use_combined_contact_core :
        index_separation >= core_side && bond_axis_gap >= 0 ?
        :passes_disjoint_atom_core_check :
        :atom_seed_boxes_need_combined_contact_core
    return (;
        q = setup.q,
        n_s = setup.n_s,
        standard_setup_core_cube_side = core_side,
        direct_core_side = core_side,
        terminal_shellification_core_side = core_side,
        source_plan_geometry_core_side = core_side,
        core_spacing = setup.core_spacing,
        bond_length_bohr = _probe_bond_length(locations),
        bond_axis,
        parent_axis_counts = _probe_get(parent, :axis_counts),
        nuclear_continuous_coordinates = locations,
        nuclear_grid_indices = nuclear_indices,
        nearest_index_coordinate_errors = coordinate_errors,
        mapped_bond_axis_center_values_near_nuclei = bond_axis_windows,
        bond_axis_index_separation = index_separation,
        left_core_box_bond_axis_range = left_core[bond_axis_index],
        right_core_box_bond_axis_range = right_core[bond_axis_index],
        minimum_disjoint_index_separation = core_side,
        bond_axis_gap,
        seed_overlap_touch_or_gap =
            bond_axis_gap < 0 ? :overlap :
            bond_axis_gap == 0 ? :touch :
            :gap,
        proposed_combined_core_bond_axis_range =
            proposed_contact_core[bond_axis_index],
        proposed_combined_core_dimensions = proposed_contact_core_dimensions,
        proposed_combined_core_support_count =
            prod(proposed_contact_core_dimensions),
        transverse_ranges_overlap,
        core_boxes,
        classification,
        call_path =
            :pqs_independent_diatomic_support_plan_to_raw_terminal_geometry,
    )
end

function _probe_print_geometry(label, parent)
    facts = _probe_geometry_facts(parent)
    isnothing(facts) && return nothing
    println(label)
    for key in keys(facts)
        println("  ", key, " = ", getproperty(facts, key))
    end
    return facts
end

_probe_count_equal(values, target) =
    isnothing(values) ? nothing : count(==(target), values)

function _probe_outer_mismatch_roles(scaffold)
    regions = _probe_get(scaffold, :regions, ())
    isempty(regions) && return ()
    return Tuple(
        _probe_get(region, :role)
        for region in regions
        if _probe_get(region, :region_kind) === :outer_mismatch_slab
    )
end

function _probe_central_gap_type(scaffold)
    isnothing(scaffold) && return nothing
    _probe_get(scaffold, :central_midpoint_slab_count, 0) > 0 &&
        return :midpoint_slab
    _probe_get(scaffold, :central_distorted_product_box_count, 0) > 0 &&
        return :distorted_product_box
    return :none
end

function _probe_topology_facts(stage)
    low_order_shellification = _probe_get(stage, :low_order_shellification)
    low_order_units = _probe_get(stage, :low_order_units)
    low_order_transforms = _probe_get(stage, :low_order_transforms)
    scaffold =
        _probe_first_notnothing(
            _probe_get(stage, :shellification_scaffold, nothing),
            _probe_get(low_order_shellification, :shellification_scaffold, nothing),
            _probe_get(low_order_shellification, :scaffold, nothing),
            _probe_get(low_order_units, :shellification_scaffold, nothing),
            _probe_get(low_order_transforms, :shellification_scaffold, nothing),
            nothing,
        )
    terminal =
        _probe_first_notnothing(
            _probe_get(low_order_shellification, :terminal_shellification, nothing),
            low_order_shellification,
            nothing,
        )
    unit_inventory =
        _probe_first_notnothing(
            _probe_get(low_order_units, :unit_inventory, nothing),
            _probe_get(low_order_transforms, :unit_inventory, nothing),
            nothing,
        )
    lowering_inventory =
        _probe_first_notnothing(
            _probe_get(low_order_units, :lowering_contract_inventory, nothing),
            _probe_get(low_order_transforms, :lowering_contract_inventory, nothing),
            nothing,
        )
    ordered_region_roles =
        _probe_first_notnothing(
            _probe_get(scaffold, :ordered_region_roles, nothing),
            _probe_get(unit_inventory, :terminal_region_roles, nothing),
            (),
        )
    support_counts =
        _probe_first_notnothing(
            _probe_get(unit_inventory, :support_counts, nothing),
            !isnothing(scaffold) ?
            Tuple(_probe_get(region, :support_count) for region in scaffold.regions) :
            nothing,
            nothing,
        )
    atom_contact_core_count =
        _probe_count_equal(ordered_region_roles, :atom_contact_core)
    atom_contact_core_support_count =
        isnothing(support_counts) ?
        nothing :
        sum(
            index -> ordered_region_roles[index] === :atom_contact_core ?
                     support_counts[index] :
                     0,
            eachindex(support_counts);
            init = 0,
        )
    coverage = _probe_get(scaffold, :coverage)
    return (;
        shellification_status =
            _probe_first_notnothing(
                _probe_get(stage, :shellification_status, nothing),
                _probe_get(low_order_shellification, :status, nothing),
                _probe_get(terminal, :status, nothing),
                isnothing(scaffold) ? :not_available : :available,
            ),
        shellification_blocker =
            _probe_first_notnothing(
                _probe_get(stage, :shellification_blocker, nothing),
                _probe_get(low_order_shellification, :blocker, nothing),
                _probe_get(terminal, :blocker, nothing),
                nothing,
            ),
        shellification_blocker_message =
            _probe_first_notnothing(
                _probe_get(stage, :shellification_blocker_message, nothing),
                _probe_get(low_order_shellification, :blocker_message, nothing),
                _probe_get(terminal, :blocker_message, nothing),
                nothing,
            ),
        core_side =
            _probe_first_notnothing(
                _probe_get(scaffold, :core_side, nothing),
                _probe_get(terminal, :core_side, nothing),
                nothing,
            ),
        q =
            _probe_first_notnothing(
                _probe_get(scaffold, :q, nothing),
                _probe_get(terminal, :q, nothing),
                nothing,
            ),
        parent_axis_counts =
            _probe_first_notnothing(
                _probe_dig(stage, :route_skeleton, :parent_axis_counts),
                _probe_dig(stage, :parent, :axis_counts),
                _probe_get(terminal, :parent_axis_counts, nothing),
                nothing,
            ),
        ordered_region_roles,
        region_count = _probe_get(scaffold, :region_count),
        support_counts_by_region = support_counts,
        atom_contact_core_count,
        atom_contact_core_support_count,
        unit_count = _probe_get(unit_inventory, :unit_count),
        unit_roles = _probe_get(unit_inventory, :unit_roles),
        unit_kinds = _probe_get(unit_inventory, :unit_kinds),
        lowering_contract_count =
            _probe_get(lowering_inventory, :lowering_contract_count),
        lowering_contract_kinds =
            _probe_get(lowering_inventory, :lowering_contract_kinds),
        lowering_contract_kind_counts =
            _probe_get(lowering_inventory, :lowering_contract_kind_counts),
        central_gap_type = _probe_central_gap_type(scaffold),
        central_midpoint_slab_count =
            _probe_get(scaffold, :central_midpoint_slab_count),
        central_distorted_product_box_count =
            _probe_get(scaffold, :central_distorted_product_box_count),
        shared_molecular_shell_count =
            _probe_count_equal(ordered_region_roles, :shared_molecular_shell),
        outer_mismatch_slab_count = length(_probe_outer_mismatch_roles(scaffold)),
        outer_mismatch_slab_roles = _probe_outer_mismatch_roles(scaffold),
        coverage_duplicate_count = _probe_get(coverage, :duplicate_count),
        coverage_missing_count = _probe_get(coverage, :missing_count),
        coverage_outside_count = _probe_get(coverage, :outside_count),
        coverage_complete = _probe_get(coverage, :coverage_complete),
    )
end

function _probe_print_topology(label, stage)
    facts = _probe_topology_facts(stage)
    println(label)
    for key in keys(facts)
        _probe_print(String(key), getproperty(facts, key))
    end
    return facts
end

function _probe_h2_geometry_comparison()
    h2_spacing_inputs = (; q, n_s, reference_spacing, tail_spacing,
        q_to_core_spacing_rule, core_spacing = 0.5, xmax_parallel = 6.0,
        xmax_transverse = 4.0)
    h2_system_inputs = (;
        atom_symbols = ("H", "H"),
        nuclear_charges = (1, 1),
        atom_locations = ((0.0, 0.0, -2.0), (0.0, 0.0, 2.0)),
        nup = 1,
        ndn = 1,
        bond_axis = :z,
        bond_length = 4.0,
        radius = 4.0,
        parent_axis_counts = nothing,
        map_backend,
    )
    system = GaussletBases.cartesian_system(h2_system_inputs)
    recipe = GaussletBases.cartesian_recipe(route_inputs)
    parent = GaussletBases.cartesian_parent(
        system,
        h2_spacing_inputs,
        parent_inputs,
        recipe,
    )
    _probe_print_geometry("h2_materialized_fixture_geometry", parent)
    shells = GaussletBases.cartesian_shells(parent, h2_spacing_inputs, recipe)
    _probe_print_topology("h2_shells_topology", shells)
    units = GaussletBases.cartesian_units(parent, shells, recipe)
    _probe_print_topology("h2_units_topology", units)
    transforms = GaussletBases.cartesian_transforms(units, recipe)
    _probe_print_topology("h2_transforms_topology", transforms)
    return nothing
end

function _probe_stage_summary(stage, value)
    if stage === :cartesian_system
        _probe_print("atom_symbols", _probe_get(value, :atom_symbols))
        _probe_print("nuclear_charges", _probe_get(value, :nuclear_charges))
        _probe_print("bond_length_bohr", _probe_get(value, :bond_length))
        _probe_print("spin_counts", (; nup = _probe_get(value, :nup),
            ndn = _probe_get(value, :ndn)))
    elseif stage === :cartesian_recipe
        _probe_print("route_kind", _probe_get(value, :route_kind))
        _probe_print("supplement_policy", _probe_get(value, :supplement_policy))
        _probe_print("supplement_basis", _probe_get(value, :supplement_basis))
        _probe_print("supplement_lmax", _probe_get(value, :supplement_lmax))
    elseif stage === :cartesian_parent
        _probe_print("system_classification",
            _probe_get(value, :system_classification))
        _probe_print("parent_axis_dimensions", _probe_get(value, :axis_counts))
        _probe_print("bond_axis", _probe_get(value, :bond_axis))
        _probe_print_geometry("cr2_probe_geometry", value)
    elseif stage === :cartesian_shells
        _probe_print_topology("cr2_shells_topology", value)
        scaffold = _probe_get(value, :shellification_scaffold)
        _probe_print("terminal_region_roles",
            _probe_dig(scaffold, :ordered_region_roles))
        _probe_print("shared_shell_count",
            _probe_dig(scaffold, :shared_molecular_shell_count))
        _probe_print("route_shape", _probe_get(value, :route_shape))
    elseif stage in (:cartesian_units, :cartesian_transforms)
        _probe_print_topology("cr2_$(stage)_topology", value)
        _probe_print("retained_counts",
            _probe_count_summary(_probe_get(value, :retained_counts)))
        _probe_print("estimated_final_dimension",
            _probe_get(value, :retained_dimension))
    elseif stage === :cartesian_pair_terms
        _probe_print("pair_family_counts", _probe_get(value, :pair_family_counts))
        _probe_print("estimated_final_dimension",
            _probe_get(value, :retained_dimension))
    elseif stage === :cartesian_assembly
        target = _probe_dig(value, :diatomic_physical_gausslet_target_payload,
            :summary)
        supplement = _probe_dig(
            value,
            :diatomic_physical_gausslet_supplement_representation_payload,
            :summary,
        )
        final_basis = _probe_dig(
            value,
            :diatomic_physical_gausslet_final_basis_payload,
            :summary,
        )
        _probe_print("target_status", _probe_get(target, :status))
        _probe_print("target_blocker", _probe_get(target, :blocker))
        _probe_print("support_counts",
            _probe_count_summary(_probe_get(target, :support_counts)))
        _probe_print("retained_counts",
            _probe_count_summary(_probe_get(target, :retained_counts)))
        _probe_print("estimated_final_dimension",
            _probe_get(target, :expected_final_dimension))
        _probe_print("supplement_status", _probe_get(supplement, :status))
        _probe_print("supplement_blocker", _probe_get(supplement, :blocker))
        _probe_print("supplement_orbital_count",
            _probe_get(supplement, :orbital_count))
        _probe_print("final_basis_status", _probe_get(final_basis, :status))
        _probe_print("projected_shell_overlap_errors",
            _probe_get(final_basis, :projected_shell_overlap_errors))
    elseif stage === :cartesian_report
        _probe_print("target_status",
            _probe_get(value, :physical_gausslet_target_status))
        _probe_print("target_blocker",
            _probe_get(value, :physical_gausslet_target_blocker))
        dim = _probe_get(value, :retained_dimension)
        _probe_print("estimated_final_dimension", dim)
        _probe_print("estimated_dense_K_U_V_memory_MB",
            _probe_dense_memory(dim, length(atom_symbols)))
    elseif stage === :cartesian_materialization
        _probe_print("result_kind", _probe_get(value, :result_kind))
        _probe_print("materialized", _probe_get(value, :materialized))
        _probe_print("residual_rank", _probe_get(value, :residual_rank))
        _probe_print("ida_orbital_dimension",
            _probe_get(value, :ida_orbital_dimension))
        _probe_print("hamfile", _probe_get(value, :hamfile))
    end
    return nothing
end

function _probe_run_stage(thunk, name)
    global last_successful_public_stage, first_blocker
    started = time()
    try
        value = thunk()
        elapsed = time() - started
        println(
            "[cr2 probe stage] ",
            name,
            " elapsed_s=",
            round(elapsed; digits = 3),
            " status=ok",
        )
        _probe_stage_summary(name, value)
        last_successful_public_stage = name
        return value, true
    catch error
        elapsed = time() - started
        message = sprint(showerror, error)
        println(
            "[cr2 probe stage] ",
            name,
            " elapsed_s=",
            round(elapsed; digits = 3),
            " status=blocked",
        )
        println("  blocker_type = ", typeof(error))
        println("  blocker_message = ", message)
        first_blocker = (; stage = name, type = typeof(error), message)
        return nothing, false
    end
end

function _probe_finish()
    println("cr2_probe_summary")
    println("  last_successful_public_stage = ", last_successful_public_stage)
    if isnothing(first_blocker)
        println("  first_blocker_stage = none")
        println("  first_blocker_message = none")
    else
        println("  first_blocker_stage = ", first_blocker.stage)
        println("  first_blocker_message = ", first_blocker.message)
    end
    return nothing
end

_probe_h2_geometry_comparison()

system, ok = _probe_run_stage(:cartesian_system) do
    GaussletBases.cartesian_system(system_inputs)
end
ok || (_probe_finish(); exit(0))

recipe, ok = _probe_run_stage(:cartesian_recipe) do
    GaussletBases.cartesian_recipe(route_inputs)
end
ok || (_probe_finish(); exit(0))

parent, ok = _probe_run_stage(:cartesian_parent) do
    GaussletBases.cartesian_parent(system, spacing_inputs, parent_inputs, recipe)
end
ok || (_probe_finish(); exit(0))

shells, ok = _probe_run_stage(:cartesian_shells) do
    GaussletBases.cartesian_shells(parent, spacing_inputs, recipe)
end
ok || (_probe_finish(); exit(0))

units, ok = _probe_run_stage(:cartesian_units) do
    GaussletBases.cartesian_units(parent, shells, recipe)
end
ok || (_probe_finish(); exit(0))

transforms, ok = _probe_run_stage(:cartesian_transforms) do
    GaussletBases.cartesian_transforms(units, recipe)
end
ok || (_probe_finish(); exit(0))

pairs, ok = _probe_run_stage(:cartesian_pair_terms) do
    GaussletBases.cartesian_pair_terms(units, transforms, recipe)
end
ok || (_probe_finish(); exit(0))

assembly, ok = _probe_run_stage(:cartesian_assembly) do
    GaussletBases.cartesian_assembly(parent, shells, units, transforms, pairs, recipe)
end
ok || (_probe_finish(); exit(0))

report, ok = _probe_run_stage(:cartesian_report) do
    GaussletBases.cartesian_report(system, parent, assembly, recipe)
end
ok || (_probe_finish(); exit(0))

materialization, ok = _probe_run_stage(:cartesian_materialization) do
    GaussletBases.cartesian_materialization(report, materialization_inputs)
end
ok || (_probe_finish(); exit(0))

_probe_finish()
