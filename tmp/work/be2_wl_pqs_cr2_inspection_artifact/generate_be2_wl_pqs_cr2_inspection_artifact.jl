using Dates
using GaussletBases

const OUTPUT_DIR = @__DIR__
const BUNDLE_PATH = joinpath(
    OUTPUT_DIR,
    "be2_wl_pqs_handoff_inspection_bundle.jld2",
)
const FINGERPRINT_PATH = joinpath(
    OUTPUT_DIR,
    "be2_wl_pqs_handoff_fingerprint.tsv",
)

function _git_text(args...)
    try
        return strip(read(Cmd(["git", string.(args)...]), String))
    catch
        return "unavailable"
    end
end

function _git_dirty()
    try
        return !isempty(strip(read(`git status --porcelain --untracked-files=no`, String)))
    catch
        return true
    end
end

function _be2_pqs_cr2_inspection_assembly()
    system_inputs = (;
        atom_symbols = ("Be", "Be"), nuclear_charges = (4, 4),
        atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        radius = 15.0, parent_axis_counts = (x = 9, y = 7, z = 9),
        map_backend = :pgdg_localized_experimental,
    )
    spacing_inputs = (;
        q = 5, n_s = 5, reference_spacing = 1.0, tail_spacing = 10.0,
        q_to_core_spacing_rule = :standard_pqs_ns_equals_q,
        core_spacing = nothing,
    )
    probe_inputs = (;
        probe_parent_axis_construction = :auto,
        parent_axis_probe_backend = :pgdg_localized_experimental,
        parent_axis_probe_family = :G10,
        probe_raw_product_box_plans = false,
        raw_product_box_probe_backend = :pgdg_localized_experimental,
    )
    route_inputs = (;
        route_family = :pqs_source_box,
        route_kind = :be2_cartesian_nesting_route_driver_spine,
        route_shape = (:pqs_left, :product, :pqs_right),
        product_body_rule = :centered_single_z_slab,
        pqs_retained_rule = :boundary_comx_product_mode_selection,
        product_retained_rule = :product_doside_retained_unit,
        terms = (
            :overlap, :position_x, :position_y, :position_z,
            :x2_x, :x2_y, :x2_z, :kinetic,
        ),
        pair_factor_normalization = :density_normalized,
        support_dense_direct_allowed = false,
        reference_only_authorities = (:support_row_oracle, :dense_parent_projection),
        white_lindsey_route_shape =
            (:standard_cartesian_units, :low_order_comx_coarsening),
        white_lindsey_mapping_rule = :standard_unit_backbone_mapping_family,
        white_lindsey_nesting_rule = :unit_box_low_order_comx_coarsening,
        white_lindsey_retained_rule = :low_order_unit_comx_retained_basis,
        white_lindsey_operator_rule = :low_order_unit_operator_blocks,
        white_lindsey_benchmark_role =
            :published_cartesian_baseline_for_pqs_comparison,
    )

    system = GaussletBases.cartesian_system(system_inputs)
    recipe = GaussletBases.cartesian_recipe(route_inputs)
    parent = GaussletBases.cartesian_parent(
        system,
        spacing_inputs,
        probe_inputs,
        recipe,
    )
    shells = GaussletBases.cartesian_shells(parent, spacing_inputs, recipe)
    units = GaussletBases.cartesian_units(parent, shells, probe_inputs, recipe)
    transforms = GaussletBases.cartesian_transforms(units, recipe)
    pairs = GaussletBases.cartesian_pair_terms(units, transforms, recipe)
    assembly = GaussletBases.cartesian_assembly(
        parent,
        shells,
        units,
        transforms,
        pairs,
        recipe,
    )
    return (; assembly, system_inputs, spacing_inputs, probe_inputs, route_inputs)
end

function _add_generator_provenance!(payload, inputs)
    payload.jld2_values["producer"] = (;
        package = "GaussletBases",
        repo_commit = _git_text("rev-parse", "HEAD"),
        dirty = _git_dirty(),
        generated_at = string(now(UTC)),
        generator_entrypoint = @__FILE__,
    )
    payload.jld2_values["fixture"] = (;
        atom_symbols = collect(String.(inputs.system_inputs.atom_symbols)),
        q = inputs.spacing_inputs.q,
        n_s = inputs.spacing_inputs.n_s,
        reference_spacing = inputs.spacing_inputs.reference_spacing,
        tail_spacing = inputs.spacing_inputs.tail_spacing,
        parent_axis_counts = collect(values(inputs.system_inputs.parent_axis_counts)),
        map_backend = inputs.system_inputs.map_backend,
        route_kind = inputs.route_inputs.route_kind,
    )
    return payload
end

function main()
    mkpath(OUTPUT_DIR)
    inputs = _be2_pqs_cr2_inspection_assembly()
    payload =
        GaussletBases._pqs_source_box_route_driver_be2_cr2_inspection_bundle_payload(
            inputs.assembly,
        )
    _add_generator_provenance!(payload, inputs)
    result = GaussletBases._pqs_source_box_route_driver_write_be2_cr2_inspection_bundle(
        BUNDLE_PATH,
        FINGERPRINT_PATH,
        payload,
    )
    pqs = payload.rows[1]
    println("bundle_path=", result.jld2_path)
    println("fingerprint_path=", result.tsv_path)
    println("pqs_status=", pqs.status)
    println("cr2_read_only_inspector_ready=", pqs.cr2_read_only_inspector_ready)
    println("cr2_solver_ready=", pqs.cr2_solver_ready)
    println("white_lindsey_status=", payload.rows[2].status)
    return nothing
end

main()
