#!/usr/bin/env julia

using GaussletBases

# Editable route defaults. Override these with a config file or `name=value`
# arguments after the script path.
    route_family = :pqs_source_box
    route_kind = :be2_cartesian_nesting_route_driver_spine
    atom_symbols = ("Be", "Be")
    nuclear_charges = (4, 4)
    atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0))
    radius = 15.0
    parent_axis_counts = (x = 9, y = 7, z = 9)
    map_backend = :pgdg_localized_experimental

# Route sizing and spacing defaults.
    q = 5
    n_s = q
    reference_spacing = 1.0
    tail_spacing = 10.0
    q_to_core_spacing_rule = :standard_pqs_ns_equals_q
    core_spacing = nothing

# Optional construction probes. `:auto` runs the probe when the prerequisite
# spacing or parent-axis metadata is available.
    probe_parent_axis_construction = :auto
    parent_axis_probe_backend = :pgdg_localized_experimental
    parent_axis_probe_family = :G10
    probe_raw_product_box_plans = :auto
    raw_product_box_probe_backend = :pgdg_localized_experimental

# Source-box route recipe. The driver reports this route shape; it does not
# build a Hamiltonian or adopt a public/default route.
    route_shape = (:pqs_left, :product, :pqs_right)
    product_body_rule = :centered_single_z_slab
    pqs_retained_rule = :boundary_comx_product_mode_selection
    product_retained_rule = :product_doside_retained_unit
    terms = ( :overlap, :position_x, :position_y, :position_z,
        :x2_x, :x2_y, :x2_z, :kinetic,)
    pair_factor_normalization = :density_normalized

# These names document reference/debug paths only. They are not algorithmic
# source-box rules and should remain reported as no-go/advisory metadata.
    support_dense_direct_allowed = false
    reference_only_authorities = (:support_row_oracle, :dense_parent_projection,)

# White-Lindsey low-order benchmark recipe. This uses the same standard
# unit/box organization as the modern routes; it does not preserve the old
# code's atom-count split heuristic as a route contract.
    white_lindsey_route_shape = (:standard_cartesian_units, :low_order_comx_coarsening,)
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
    save_basis_artifact = false
    save_ham_artifact = false
    outfile = "pqs_source_box_route_driver_report.jld2"
    tsvfile = "pqs_source_box_route_driver_report.tsv"
    basisfile = "cartesian_nesting_route_driver_basis_bundle.jld2"
    hamfile = "cartesian_nesting_route_driver_ham_bundle.jld2"

# Optional config include, followed by simple `name=value` overrides.
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
        radius, parent_axis_counts, map_backend,)
    spacing_inputs = (; q, n_s, reference_spacing, tail_spacing,
        q_to_core_spacing_rule, core_spacing,)
    probe_inputs = (; probe_parent_axis_construction, parent_axis_probe_backend,
        parent_axis_probe_family, probe_raw_product_box_plans, raw_product_box_probe_backend,)

    source_box_recipe = (; route_shape, product_body_rule,
        pqs_retained_rule, product_retained_rule,
        support_dense_direct_allowed, reference_only_authorities,)
    white_lindsey_recipe = (; route_shape = white_lindsey_route_shape,
        mapping_rule = white_lindsey_mapping_rule, nesting_rule = white_lindsey_nesting_rule,
        retained_rule = white_lindsey_retained_rule, operator_rule = white_lindsey_operator_rule,
        benchmark_role = white_lindsey_benchmark_role,)
    route_recipe = (; route_family, route_kind, terms, pair_factor_normalization,
        source_box = source_box_recipe, white_lindsey = white_lindsey_recipe,)

# Build the metadata-only route report in visible stages. The helpers below
# hide field-packing boilerplate, not the route order.
    GaussletBases._pqs_source_box_route_driver_check_inputs(route_recipe)

    standard_setup = GaussletBases._pqs_source_box_route_driver_standard_setup(
        system_inputs, spacing_inputs)

    parent_axis = GaussletBases._pqs_source_box_route_driver_parent_axis(
        standard_setup, system_inputs, probe_inputs)
    parent_axis_readiness = parent_axis.parent_axis_readiness

    route_axis_counts = GaussletBases._pqs_source_box_route_driver_route_axis_counts(
        standard_setup, parent_axis, system_inputs, route_recipe)

    system_metadata = GaussletBases._pqs_source_box_route_driver_system_metadata(
        standard_setup, route_axis_counts, system_inputs)

    route_skeleton = GaussletBases._pqs_source_box_route_driver_route_skeleton(
        route_axis_counts, spacing_inputs, route_recipe)

    raw_box = GaussletBases._pqs_source_box_route_driver_raw_box_probe(
        standard_setup, route_skeleton, parent_axis, route_axis_counts,
        probe_inputs, route_recipe)

    recipe_metadata = GaussletBases._pqs_source_box_route_driver_recipe_metadata(
        standard_setup, route_axis_counts, parent_axis, raw_box,
        spacing_inputs, probe_inputs, route_recipe)

    parent_description = GaussletBases._pqs_source_box_route_driver_parent_description(
        standard_setup, parent_axis, route_axis_counts, route_skeleton, raw_box)
    route_facts = GaussletBases._pqs_source_box_route_driver_route_facts(route_skeleton)
    contract = GaussletBases._pqs_source_box_route_driver_contract_metadata(route_recipe)
    diagnostics = GaussletBases._pqs_source_box_route_driver_diagnostics(
        standard_setup, parent_axis, route_axis_counts, route_skeleton, raw_box, contract)

    report = GaussletBases._pqs_source_box_route_driver_report(
        standard_setup, parent_axis, route_axis_counts, raw_box,
        system_metadata, recipe_metadata, parent_description,
        route_skeleton, route_facts, contract, diagnostics)

    materialization = GaussletBases._pqs_source_box_route_driver_materialization(
        report;
        materialize_route, save_basis_artifact, save_ham_artifact, basisfile, hamfile,
        white_lindsey_expansion, white_lindsey_Z,)

# Short screen summary first; detailed sections follow below.
    retained_counts = route_facts.retained_counts
    retained_dimension = route_facts.retained_dimension

    println("Cartesian nesting route driver")
    @show route_family route_kind q radius
    if route_family == :pqs_source_box
        @show route_shape product_body_rule pair_factor_normalization
    else
        @show white_lindsey_route_shape white_lindsey_mapping_rule white_lindsey_nesting_rule
        @show white_lindsey_Z length(white_lindsey_expansion.exponents)
    end
    @show standard_setup.n_s standard_setup.core_cube_side standard_setup.core_spacing
    @show standard_setup.spacing.q_to_core_spacing_rule_status
    @show parent_axis_readiness.status parent_axis_readiness.parent_axis_counts_status
    @show route_axis_counts.parent_axis_counts_source route_axis_counts.parent_axis_counts
    @show diagnostics.parent_axis_probe_requested diagnostics.parent_axis_probe_status
    @show diagnostics.raw_product_box_probe_requested diagnostics.raw_product_box_probe_status
    @show retained_counts retained_dimension
    @show materialization.basis_artifact_status materialization.basis_artifact_written
    @show materialization.status materialization.ham_bundle_export_status
    @show materialization.route_configured_system_classification
    @show materialization.route_configured_shellization_request_status
    @show materialization.route_configured_shellization_planning_family
    @show materialization.route_configured_midpoint_slab_status
    @show materialization.shellization_source materialization.route_configured_shellization_consumed
    @show materialization.ham_artifact_status materialization.ham_artifact_written

# Detailed sections and optional artifacts stay behind helpers so the driver
# remains an editable run recipe rather than a report-schema implementation.
    GaussletBases._pqs_source_box_route_driver_print_details(report)
    if materialize_route || save_basis_artifact || save_ham_artifact
        GaussletBases._pqs_source_box_route_driver_print_materialization(materialization)
    end
    GaussletBases._pqs_source_box_route_driver_save( report;
        save_artifact, save_tsv, outfile, tsvfile, materialization,)

println("driver complete")
