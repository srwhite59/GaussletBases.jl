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
core_spacing = 0.5
xmax_parallel = 6.0
xmax_transverse = 4.0

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
supplement_basis = "cc-pVTZ"
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
    elseif stage === :cartesian_shells
        scaffold = _probe_get(value, :shellification_scaffold)
        _probe_print("terminal_region_roles",
            _probe_dig(scaffold, :ordered_region_roles))
        _probe_print("shared_shell_count",
            _probe_dig(scaffold, :shared_molecular_shell_count))
        _probe_print("route_shape", _probe_get(value, :route_shape))
    elseif stage in (:cartesian_units, :cartesian_transforms)
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
