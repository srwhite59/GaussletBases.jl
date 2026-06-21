#!/usr/bin/env julia

# Canonical human-facing Cartesian producer template.
# Keep the public stage sequence visible here. Test instrumentation, private
# solver controls, fixture values, and route-internal provider switches belong
# in the driver harness or route implementation, not this copyable script.

using GaussletBases

# Physical system

atom_symbols = ("Be", "Be")
nuclear_charges = (4, 4)
atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0))
nup = nothing
ndn = nothing
bond_axis = nothing
bond_length = nothing
radius = 15.0
parent_axis_counts = (x = 9, y = 7, z = 9)
map_backend = :pgdg_localized_experimental

# Basis and spacing

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

# Route recipe

route_family = :pqs_source_box
route_kind = :be2_cartesian_nesting_route_driver_spine
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
supplement_policy = nothing
run_final_basis = nothing
run_h1 = true
run_h1_j = false

# Output request

save_artifact = false
save_tsv = false
materialize_route = false
save_ida_hamiltonian = false
hamiltonian_output =
    save_ida_hamiltonian ? :cartesian_ida_hamiltonian : nothing
outfile = "cartesian_ham_builder_report.jld2"
tsvfile = "cartesian_ham_builder_report.tsv"
hamfile = "cartesian_ida_hamiltonian.jld2"

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
    supplement_policy, run_final_basis, run_h1, run_h1_j)
materialization_inputs = (; materialize_route,
    save_ham_artifact = save_ida_hamiltonian, hamfile,
    hamiltonian_output)
save_inputs = (; save_artifact, save_tsv, outfile, tsvfile,
    input_path = nothing)

# Begin actual construction

system = GaussletBases.cartesian_system(system_inputs)
recipe = GaussletBases.cartesian_recipe(route_inputs)
parent = GaussletBases.cartesian_parent(system, spacing_inputs, parent_inputs, recipe)
shells = GaussletBases.cartesian_shells(parent, spacing_inputs, recipe)
units = GaussletBases.cartesian_units(parent, shells, recipe)
@time "Transforming: " transforms = GaussletBases.cartesian_transforms(units, recipe)
@time "Pair terms: " pairs = GaussletBases.cartesian_pair_terms(units, transforms, recipe)
@time "Assembly: " assembly =
    GaussletBases.cartesian_assembly(parent, shells, units, transforms, pairs, recipe)
report = GaussletBases.cartesian_report(system, parent, assembly, recipe)
materialization = GaussletBases.cartesian_materialization(
    report, transforms.terminal_basis_realization, materialization_inputs)

GaussletBases.cartesian_print_summary(report, materialization)
GaussletBases.cartesian_print_details(report, materialization)
GaussletBases.cartesian_save(report, save_inputs, materialization)

println("driver complete")
