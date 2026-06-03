#!/usr/bin/env julia

using GaussletBases

route_kind = :be2_pqs_source_box_development_spine
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
terms = (
    :overlap,
    :position_x,
    :position_y,
    :position_z,
    :x2_x,
    :x2_y,
    :x2_z,
    :kinetic,
)
pair_factor_normalization = :density_normalized
support_dense_direct_allowed = false
reference_only_authorities = (
    :support_row_oracle,
    :dense_parent_projection,
)
save_artifact = false
save_tsv = false
outfile = "pqs_source_box_route_driver_report.jld2"
tsvfile = "pqs_source_box_route_driver_report.tsv"

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

report = GaussletBases._pqs_source_box_route_driver_dry_run(
    ;
    route_kind,
    atom_symbols,
    nuclear_charges,
    atom_locations,
    radius,
    parent_axis_counts,
    map_backend,
    q,
    n_s,
    reference_spacing,
    tail_spacing,
    q_to_core_spacing_rule,
    core_spacing,
    probe_parent_axis_construction,
    parent_axis_probe_backend,
    parent_axis_probe_family,
    probe_raw_product_box_plans,
    raw_product_box_probe_backend,
    route_shape,
    product_body_rule,
    pqs_retained_rule,
    product_retained_rule,
    terms,
    pair_factor_normalization,
    support_dense_direct_allowed,
    reference_only_authorities,
)

standard_setup = report.standard_setup
parent_axis_readiness = report.parent_axis_readiness
route_axis_counts = report.route_axis_counts
diagnostics = report.diagnostics
retained_counts = report.retained_counts
retained_dimension = report.retained_dimension

println("PQS source-box route driver")
@show route_kind
@show q
@show radius
@show route_shape
@show product_body_rule
@show pair_factor_normalization
@show standard_setup.n_s
@show standard_setup.core_cube_side
@show standard_setup.core_spacing
@show standard_setup.spacing.q_to_core_spacing_rule_status
@show parent_axis_readiness.status
@show parent_axis_readiness.parent_axis_counts_status
@show route_axis_counts.parent_axis_counts_source
@show route_axis_counts.parent_axis_counts
@show diagnostics.parent_axis_probe_requested
@show diagnostics.parent_axis_probe_status
@show diagnostics.raw_product_box_probe_requested
@show diagnostics.raw_product_box_probe_status
@show retained_counts
@show retained_dimension

GaussletBases._pqs_source_box_route_driver_print_details(report)
GaussletBases._pqs_source_box_route_driver_save(
    report;
    save_artifact,
    save_tsv,
    outfile,
    tsvfile,
)

println("driver complete")
