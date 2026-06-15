#!/usr/bin/env julia

using GaussletBases

# Initialization of parameters

route_family = :pqs_source_box
route_kind = :be2_cartesian_nesting_route_driver_spine
atom_symbols = ("Be", "Be")
nuclear_charges = (4, 4)
atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0))
bond_axis = nothing
bond_length = nothing
radius = 15.0
parent_axis_counts = (x = 9, y = 7, z = 9)
map_backend = :pgdg_localized_experimental

q = 5
n_s = q
reference_spacing = 1.0
tail_spacing = 10.0
q_to_core_spacing_rule = :standard_pqs_ns_equals_q
core_spacing = nothing
xmax_parallel = nothing
xmax_transverse = nothing

parent_axis_family = :G10
parent_axis_bundle_backend = :pgdg_localized_experimental
parent_mapping_rule = :identity_mapping
parent_mapping_Z = nothing
parent_mapping_d = nothing
parent_mapping_tail_spacing = tail_spacing

route_shape = (:pqs_left, :product, :pqs_right)
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
white_lindsey_Z = 2.0
white_lindsey_expansion = coulomb_gaussian_expansion(doacc = false)
supplement_policy = nothing
run_final_basis = nothing
run_h1 = true
run_h1_j = true
run_private_rhf = false
private_rhf_electron_count = nothing
private_rhf_fixture_role = :route_smoke
private_rhf_mixing_kind = :fock_diis
private_rhf_max_iterations = 25
private_rhf_density_atol = 1.0e-8
private_rhf_energy_atol = 1.0e-10
private_rhf_residual_atol = 1.0e-8
private_rhf_trace_atol = private_rhf_density_atol
private_rhf_idempotency_atol = private_rhf_density_atol
private_rhf_max_history = nothing
private_rhf_diis_start_iteration = 2
private_rhf_diis_regularization = 1.0e-12
private_rhf_diis_coefficient_max_abs = 25.0

save_artifact = false
save_tsv = false
materialize_route = false
probe_route_configured_one_center_materializer = false
save_basis_artifact = false
save_ham_artifact = false
stop_after_stage = nothing
# `nothing` means the helper derives these from route metadata:
# backend from map/probe backend, nside from n_s.
materializer_backend = nothing
materializer_nside = nothing
route_configured_diatomic_ham_interaction_treatment = :ggt_nearest
outfile = "cartesian_ham_builder_report.jld2"
tsvfile = "cartesian_ham_builder_report.tsv"
basisfile = "cartesian_nesting_route_driver_basis_bundle.jld2"
hamfile = "cartesian_nesting_route_driver_ham_bundle.jld2"
driver_input_path = nothing

inputs = String[]
if length(ARGS) > 0
    if occursin("=", ARGS[1])
        inputs = ARGS
    else
        fullpath = isabspath(ARGS[1]) ? ARGS[1] : joinpath(pwd(), ARGS[1])
        driver_input_path = fullpath
        println("including ", fullpath)
        include(fullpath)
        length(ARGS) > 1 && (inputs = ARGS[2:end])
    end
end
eval.(Meta.parse.(inputs))

system_inputs = (; atom_symbols, nuclear_charges, atom_locations,
    bond_axis, bond_length, radius, parent_axis_counts, map_backend)
spacing_inputs = (; q, n_s, reference_spacing, tail_spacing,
    q_to_core_spacing_rule, core_spacing, xmax_parallel, xmax_transverse)
parent_inputs = (; parent_axis_bundle_backend,
    parent_axis_family, parent_mapping_rule,
    parent_mapping_Z, parent_mapping_d, parent_mapping_tail_spacing)
private_rhf_inputs = (; run_private_rhf, private_rhf_electron_count,
    private_rhf_fixture_role, private_rhf_mixing_kind, private_rhf_max_iterations,
    private_rhf_density_atol, private_rhf_energy_atol, private_rhf_residual_atol,
    private_rhf_trace_atol, private_rhf_idempotency_atol, private_rhf_max_history,
    private_rhf_diis_start_iteration, private_rhf_diis_regularization,
    private_rhf_diis_coefficient_max_abs)
route_inputs = (; route_family, route_kind, route_shape, product_body_rule,
    pqs_retained_rule, product_retained_rule, terms, pair_factor_normalization,
    white_lindsey_route_shape, white_lindsey_mapping_rule,
    white_lindsey_nesting_rule, white_lindsey_retained_rule,
    white_lindsey_operator_rule,
    supplement_policy, run_final_basis, run_h1,
    run_h1_j, private_rhf_inputs)
materialization_inputs = (; materialize_route, probe_route_configured_one_center_materializer,
    save_basis_artifact, save_ham_artifact, basisfile, hamfile,
    materializer_backend, materializer_nside,
    route_configured_diatomic_ham_interaction_treatment,
    white_lindsey_expansion, white_lindsey_Z)
save_inputs = (; save_artifact, save_tsv, outfile, tsvfile,
    input_path = driver_input_path)


# Begin actual construction

function _cartesian_driver_stage(name)
    println("[driver stage] ", name)
end

function _cartesian_driver_stop_after(name)
    stop_after_stage === name || return nothing
    println("[driver stop] after ", name)
    exit(0)
end

_cartesian_driver_stage(:cartesian_system)
system = GaussletBases.cartesian_system(system_inputs)
_cartesian_driver_stop_after(:cartesian_system)
_cartesian_driver_stage(:cartesian_recipe)
recipe = GaussletBases.cartesian_recipe(route_inputs)
_cartesian_driver_stop_after(:cartesian_recipe)
_cartesian_driver_stage(:cartesian_parent)
parent = GaussletBases.cartesian_parent(system, spacing_inputs, parent_inputs, recipe)
_cartesian_driver_stop_after(:cartesian_parent)
_cartesian_driver_stage(:cartesian_shells)
shells = GaussletBases.cartesian_shells(
    parent,
    spacing_inputs,
    recipe,
)
_cartesian_driver_stop_after(:cartesian_shells)
_cartesian_driver_stage(:cartesian_units)
units = GaussletBases.cartesian_units(parent, shells, recipe)
_cartesian_driver_stop_after(:cartesian_units)
_cartesian_driver_stage(:cartesian_transforms)
@time "Transforming: " transforms = GaussletBases.cartesian_transforms(units, recipe)
_cartesian_driver_stop_after(:cartesian_transforms)
_cartesian_driver_stage(:cartesian_pair_terms)
@time "Pair terms: " pairs = GaussletBases.cartesian_pair_terms(units, transforms, recipe)
_cartesian_driver_stop_after(:cartesian_pair_terms)
_cartesian_driver_stage(:cartesian_assembly)
@time "Assembly: " assembly = GaussletBases.cartesian_assembly(parent, shells, units, transforms, pairs, recipe)
_cartesian_driver_stop_after(:cartesian_assembly)
_cartesian_driver_stage(:cartesian_report)
report = GaussletBases.cartesian_report(system, parent, assembly, recipe)
_cartesian_driver_stop_after(:cartesian_report)
_cartesian_driver_stage(:cartesian_materialization)
materialization = GaussletBases.cartesian_materialization(report, materialization_inputs)
_cartesian_driver_stop_after(:cartesian_materialization)

_cartesian_driver_stage(Symbol("cartesian_print/save"))
GaussletBases.cartesian_print_summary(report, materialization)
GaussletBases.cartesian_print_details(report, materialization)
GaussletBases.cartesian_save(report, save_inputs, materialization)
_cartesian_driver_stop_after(Symbol("cartesian_print/save"))

println("driver complete")
