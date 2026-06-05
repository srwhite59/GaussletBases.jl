#!/usr/bin/env julia

using GaussletBases

# Initialization of parameters

route_family = :pqs_source_box
route_kind = :be2_cartesian_nesting_route_driver_spine
atom_symbols = ("Be", "Be")
nuclear_charges = (4, 4)
atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0))
radius = 15.0
parent_axis_counts = (x = 9, y = 7, z = 9)
map_backend = :pgdg_localized_experimental

q = 5
n_s = q
reference_spacing = 1.0
tail_spacing = 10.0
q_to_core_spacing_rule = :standard_pqs_ns_equals_q
core_spacing = nothing

probe_parent_axis_construction = :auto
parent_axis_probe_backend = :pgdg_localized_experimental
parent_axis_probe_family = :G10
probe_raw_product_box_plans = :auto
raw_product_box_probe_backend = :pgdg_localized_experimental

route_shape = (:pqs_left, :product, :pqs_right)
product_body_rule = :centered_single_z_slab
pqs_retained_rule = :boundary_comx_product_mode_selection
product_retained_rule = :product_doside_retained_unit
terms = (:overlap, :position_x, :position_y, :position_z,
    :x2_x, :x2_y, :x2_z, :kinetic)
pair_factor_normalization = :density_normalized

support_dense_direct_allowed = false
reference_only_authorities = (:support_row_oracle, :dense_parent_projection)

white_lindsey_route_shape = (:standard_cartesian_units, :low_order_comx_coarsening)
white_lindsey_mapping_rule = :standard_unit_backbone_mapping_family
white_lindsey_nesting_rule = :unit_box_low_order_comx_coarsening
white_lindsey_retained_rule = :low_order_unit_comx_retained_basis
white_lindsey_operator_rule = :low_order_unit_operator_blocks
white_lindsey_benchmark_role = :published_cartesian_baseline_for_pqs_comparison
white_lindsey_Z = 2.0
white_lindsey_expansion = coulomb_gaussian_expansion(doacc = false)

save_artifact = false
save_tsv = false
materialize_route = false
probe_route_configured_one_center_materializer = false
probe_route_configured_diatomic_atom_growth_materializer = false
low_order_shellization_policy = nothing
save_basis_artifact = false
save_ham_artifact = false
# `nothing` means the helper derives these from route metadata:
# backend from map/probe backend, nside from n_s.
materializer_backend = nothing
materializer_nside = nothing
route_configured_diatomic_ham_interaction_treatment = :ggt_nearest
outfile = "cartesian_ham_builder_report.jld2"
tsvfile = "cartesian_ham_builder_report.tsv"
basisfile = "cartesian_nesting_route_driver_basis_bundle.jld2"
hamfile = "cartesian_nesting_route_driver_ham_bundle.jld2"

inputs = String[]
if length(ARGS) > 0
    if occursin("=", ARGS[1])
        inputs = ARGS
    else
        fullpath = isabspath(ARGS[1]) ? ARGS[1] : joinpath(pwd(), ARGS[1])
        println("including ", fullpath)
        include(fullpath)
        length(ARGS) > 1 && (inputs = ARGS[2:end])
    end
end
eval.(Meta.parse.(inputs))

system_inputs = (; atom_symbols, nuclear_charges, atom_locations,
    radius, parent_axis_counts, map_backend)
spacing_inputs = (; q, n_s, reference_spacing, tail_spacing,
    q_to_core_spacing_rule, core_spacing)
parent_inputs = (; probe_parent_axis_construction, parent_axis_probe_backend,
    parent_axis_probe_family)
route_probe_inputs = (; probe_raw_product_box_plans, raw_product_box_probe_backend)
route_inputs = (; route_family, route_kind, route_shape, product_body_rule,
    pqs_retained_rule, product_retained_rule, terms, pair_factor_normalization,
    support_dense_direct_allowed, reference_only_authorities,
    white_lindsey_route_shape, white_lindsey_mapping_rule,
    white_lindsey_nesting_rule, white_lindsey_retained_rule,
    white_lindsey_operator_rule, white_lindsey_benchmark_role)
materialization_inputs = (; materialize_route, probe_route_configured_one_center_materializer,
    probe_route_configured_diatomic_atom_growth_materializer,
    low_order_shellization_policy,
    save_basis_artifact, save_ham_artifact, basisfile, hamfile,
    materializer_backend, materializer_nside,
    route_configured_diatomic_ham_interaction_treatment,
    white_lindsey_expansion, white_lindsey_Z)
save_inputs = (; save_artifact, save_tsv, outfile, tsvfile)


# Begin actual construction

system = GaussletBases.cartesian_system(system_inputs)
recipe = GaussletBases.cartesian_recipe(route_inputs)
parent = GaussletBases.cartesian_parent(system, spacing_inputs, parent_inputs, recipe)
shells = GaussletBases.cartesian_shells(parent, spacing_inputs, recipe)
units = GaussletBases.cartesian_units(parent, shells, route_probe_inputs, recipe)
@time "Transforming: " transforms = GaussletBases.cartesian_transforms(units, recipe)
@time "Pair terms: " pairs = GaussletBases.cartesian_pair_terms(units, transforms, recipe)
@time "Assembly: " assembly = GaussletBases.cartesian_assembly(parent, shells, units, transforms, pairs, recipe)
report = GaussletBases.cartesian_report(system, parent, assembly, recipe)
materialization = GaussletBases.cartesian_materialization(report, materialization_inputs)

GaussletBases.cartesian_print_summary(report, materialization)
GaussletBases.cartesian_print_details(report, materialization)
GaussletBases.cartesian_save(report, save_inputs, materialization)

println("driver complete")
