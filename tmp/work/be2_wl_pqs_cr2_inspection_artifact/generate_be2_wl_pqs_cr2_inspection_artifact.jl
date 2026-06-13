using Dates
using GaussletBases
using LinearAlgebra

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

function _be2_white_lindsey_atom_growth_route()
    system_inputs = (;
        atom_symbols = ("Be", "Be"), nuclear_charges = (4, 4),
        atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        radius = 15.0, parent_axis_counts = (x = 9, y = 7, z = 9),
        map_backend = :pgdg_localized_experimental,
    )
    spacing_inputs = (;
        q = 5, n_s = 5, reference_spacing = 1.0, tail_spacing = 10.0,
        q_to_core_spacing_rule = :standard_pqs_ns_equals_q,
        core_spacing = 0.15,
    )
    parent_inputs = (;
        probe_parent_axis_construction = :auto,
        parent_axis_probe_backend = :pgdg_localized_experimental,
        parent_axis_probe_family = :G10,
    )
    route_probe_inputs = (;
        probe_raw_product_box_plans = :auto,
        raw_product_box_probe_backend = :pgdg_localized_experimental,
    )
    route_inputs = (;
        route_family = :white_lindsey_low_order,
        route_kind = :be2_cartesian_diatomic_driver_config_smoke,
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
    materialization_inputs = (;
        materialize_route = true,
        probe_route_configured_one_center_materializer = false,
        private_global_overlap_requested = false,
        private_global_overlap_global_dimension = nothing,
        private_global_overlap_inputs = (;),
        probe_route_configured_diatomic_atom_growth_materializer = true,
        low_order_shellization_policy = :atom_growth_complete_rectangular,
        save_basis_artifact = false,
        save_ham_artifact = false,
        basisfile = joinpath(OUTPUT_DIR, "unused_wl_basis.jld2"),
        hamfile = joinpath(OUTPUT_DIR, "unused_wl_ham.jld2"),
        materializer_backend = :pgdg_localized_experimental,
        materializer_nside = 5,
        route_configured_diatomic_ham_interaction_treatment = :ggt_nearest,
        white_lindsey_expansion = GaussletBases.coulomb_gaussian_expansion(doacc = false),
        white_lindsey_Z = 2.0,
    )

    system = GaussletBases.cartesian_system(system_inputs)
    recipe = GaussletBases.cartesian_recipe(route_inputs)
    parent = GaussletBases.cartesian_parent(system, spacing_inputs, parent_inputs, recipe)
    shells = GaussletBases.cartesian_shells(
        parent,
        spacing_inputs,
        recipe;
        low_order_shellization_policy = :atom_growth_complete_rectangular,
        probe_route_configured_diatomic_atom_growth_materializer = true,
    )
    units = GaussletBases.cartesian_units(parent, shells, route_probe_inputs, recipe)
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
    report = GaussletBases.cartesian_report(system, parent, assembly, recipe)
    materialization =
        GaussletBases.cartesian_materialization(report, materialization_inputs)
    return (;
        materialization, system_inputs, spacing_inputs, parent_inputs,
        route_probe_inputs, route_inputs, materialization_inputs,
    )
end

function _be2_wl_unavailable_payload(payload, blocker)
    pqs_system = payload.jld2_values["routes/pqs_source_box/system"]
    wl_fingerprint = (route_label = :white_lindsey, status = :unavailable,
        blocker, final_dimension = nothing, pre_final_dimension = nothing,
        support_weight_count = nothing, one_body_shape = nothing,
        two_body_shape = nothing, h1_lowest = nothing,
        h1_symmetry_defect = nothing, one_body_finite = false,
        two_body_finite = false, density_gauge = :not_applicable,
        raw_pair_factor_convention = :not_applicable,
        cr2_read_only_inspector_ready = false, cr2_solver_ready = false,
        cr2_export_ready = false, cr2_handoff_blocker = blocker,
        nuclear_repulsion = pqs_system.nuclear_repulsion,
        electron_count = pqs_system.electron_count,
        spin_sector = pqs_system.spin_sector)
    payload.jld2_values["routes/white_lindsey/route"] = (;
        label = :white_lindsey,
        family = :white_lindsey_low_order,
        kind = :route_configured_diatomic_atom_growth,
        status = :unavailable,
        blocker,
    )
    payload.jld2_values["routes/white_lindsey/readiness"] = (;
        cr2_read_only_inspector_ready = false,
        cr2_solver_ready = false,
        cr2_export_ready = false,
        cr2_handoff_blocker = blocker,
    )
    return (;
        jld2_values = payload.jld2_values,
        rows = (payload.rows[1], wl_fingerprint),
        fingerprint_columns = payload.fingerprint_columns,
    )
end

function _be2_wl_symmetry_defect(matrix)
    return norm(matrix - matrix')
end

function _be2_populate_white_lindsey_route(payload, wl_route)
    materialization = wl_route.materialization
    probe = materialization.route_configured_diatomic_atom_growth_materializer_probe
    basis_adapter = probe.basis_adapter
    ham_adapter = probe.ham_adapter
    if isnothing(basis_adapter) || isnothing(ham_adapter)
        blocker = something(
            materialization.route_configured_diatomic_atom_growth_ham_adapter_blocker,
            :missing_route_configured_diatomic_atom_growth_ham_payload,
        )
        return _be2_wl_unavailable_payload(payload, blocker)
    elseif basis_adapter.status != :available_route_configured_diatomic_atom_growth_basis_adapter
        return _be2_wl_unavailable_payload(payload, basis_adapter.blocker)
    elseif ham_adapter.status != :available_route_configured_diatomic_ham_adapter
        return _be2_wl_unavailable_payload(payload, ham_adapter.blocker)
    end

    pqs_system = payload.jld2_values["routes/pqs_source_box/system"]
    materialized_report = materialization.materialized_report
    representation = basis_adapter.representation
    operators = ham_adapter.operators
    overlap = Matrix{Float64}(operators.overlap)
    one_body = Matrix{Float64}(operators.one_body_hamiltonian)
    interaction = Matrix{Float64}(operators.interaction_matrix)
    final_weights = Vector{Float64}(basis_adapter.final_integral_weights)
    final_dimension = basis_adapter.retained_dimension
    expected_shape = (final_dimension, final_dimension)
    one_body_finite = all(isfinite, overlap) && all(isfinite, one_body)
    two_body_finite = all(isfinite, interaction)
    one_body_symmetry_defect = _be2_wl_symmetry_defect(one_body)
    two_body_symmetry_defect = _be2_wl_symmetry_defect(interaction)
    hf_convention_blocker =
        :missing_reviewed_density_density_hf_fock_energy_convention
    overlap_identity_defect = norm(overlap - I, Inf)
    route_available =
        size(overlap) == expected_shape &&
        size(one_body) == expected_shape &&
        size(interaction) == expected_shape &&
        length(final_weights) == final_dimension &&
        all(isfinite, final_weights) &&
        all(>(0.0), final_weights) &&
        one_body_finite &&
        two_body_finite
    status =
        route_available ?
        :available_route_configured_diatomic_atom_growth_ham_payload :
        :blocked_route_configured_diatomic_atom_growth_ham_payload
    blocker =
        route_available ?
        nothing :
        :route_configured_diatomic_atom_growth_validation_failed

    wl_fingerprint = (route_label = :white_lindsey, status,
        blocker, final_dimension, pre_final_dimension = nothing,
        support_weight_count = length(final_weights),
        one_body_shape = size(one_body), two_body_shape = size(interaction),
        h1_lowest = nothing, h1_symmetry_defect = one_body_symmetry_defect,
        one_body_finite, two_body_finite, density_gauge = :not_applicable,
        raw_pair_factor_convention = :not_applicable,
        cr2_read_only_inspector_ready = route_available,
        cr2_solver_ready = false, cr2_export_ready = false,
        cr2_handoff_blocker =
            route_available ? :missing_hfdmrg_density_density_contract : blocker,
        nuclear_repulsion = pqs_system.nuclear_repulsion,
        electron_count = pqs_system.electron_count,
        spin_sector = pqs_system.spin_sector)

    route_blocker =
        route_available ? nothing : :route_configured_diatomic_atom_growth_validation_failed
    payload.jld2_values["routes/white_lindsey/route"] = (;
        label = :white_lindsey,
        family = :white_lindsey_low_order,
        kind = :route_configured_diatomic_atom_growth,
        status,
        blocker = route_blocker,
        data_authority = :final_basis_ordinary_cartesian_qiu_white,
        source_path = :route_configured_diatomic_atom_growth_materialization,
        route_default_behavior_changed = false,
        density_density_hf_convention_status = hf_convention_blocker,
        density_density_hf_convention_blocker = hf_convention_blocker,
    )
    payload.jld2_values["routes/white_lindsey/readiness"] = (;
        cr2_read_only_inspector_ready = route_available,
        cr2_solver_ready = false,
        cr2_export_ready = false,
        cr2_handoff_blocker =
            route_available ? :missing_hfdmrg_density_density_contract : route_blocker,
        solver_ready = false,
        hfdmrg_ready = false,
        rhf_ready = false,
        dmrg_ready = false,
        hamv6_export_ready = false,
    )
    payload.jld2_values["routes/white_lindsey/system"] = (;
        status,
        blocker = route_blocker,
        atom_symbols = collect(String.(wl_route.system_inputs.atom_symbols)),
        nuclear_charges = collect(Float64.(wl_route.system_inputs.nuclear_charges)),
        nuclear_coordinates = reduce(
            vcat,
            (permutedims(collect(coord)) for coord in wl_route.system_inputs.atom_locations),
        ),
        bond_axis = materialization.route_configured_bond_axis,
        nuclear_repulsion = pqs_system.nuclear_repulsion,
        electron_count = pqs_system.electron_count,
        spin_sector = pqs_system.spin_sector,
        route_family = :white_lindsey_low_order,
        q = wl_route.spacing_inputs.q,
        n_s = wl_route.spacing_inputs.n_s,
        reference_spacing = wl_route.spacing_inputs.reference_spacing,
        tail_spacing = wl_route.spacing_inputs.tail_spacing,
        core_spacing = wl_route.spacing_inputs.core_spacing,
        parent_axis_counts = collect(values(wl_route.system_inputs.parent_axis_counts)),
        map_backend = wl_route.system_inputs.map_backend,
    )
    payload.jld2_values["routes/white_lindsey/final_basis"] = (;
        status,
        blocker = route_blocker,
        representation_kind = :final_basis_ordinary_cartesian_qiu_white,
        basis_kind = representation.metadata.basis_kind,
        parent_kind = representation.metadata.parent_kind,
        final_dimension,
        parent_dimension = representation.metadata.parent_dimension,
        final_integral_weights = final_weights,
        basis_labels = String[String(label) for label in representation.metadata.basis_labels],
        basis_centers = Matrix{Float64}(representation.metadata.basis_centers),
        materialized_report_kind = materialization.materialized_report_kind,
        shellization_source = materialization.shellization_source,
        shellization_authority = materialized_report.shellization_authority,
        shellification_materialization_kind =
            materialization.route_configured_diatomic_atom_growth_materializer_probe.materialization.object_kind,
        low_order_shellization_policy =
            materialization.low_order_shellization_policy_resolved,
        overlap_convention = :orthonormal_identity_diagnostic_checked_not_stored,
        overlap_matrix_stored = false,
        overlap_identity_defect,
    )
    payload.jld2_values["routes/white_lindsey/one_body"] = (;
        status,
        blocker = route_blocker,
        representation_kind = :final_basis_one_body_hamiltonian,
        hamiltonian = one_body,
        kinetic_one_body =
            isnothing(operators.kinetic_one_body) ?
            Float64[] :
            Matrix{Float64}(operators.kinetic_one_body),
        nuclear_term_storage =
            isnothing(operators.nuclear_term_storage) ?
            :unavailable :
            operators.nuclear_term_storage,
        nuclear_one_body_by_center_count =
            isnothing(operators.nuclear_one_body_by_center) ?
            0 :
            length(operators.nuclear_one_body_by_center),
        nuclear_one_body_by_center =
            isnothing(operators.nuclear_one_body_by_center) ?
            Matrix{Float64}[] :
            Matrix{Float64}[Matrix{Float64}(matrix) for matrix in operators.nuclear_one_body_by_center],
    )
    payload.jld2_values["routes/white_lindsey/two_body"] = (;
        status,
        blocker = route_blocker,
        representation_kind = :final_basis_density_density_matrix,
        interaction_matrix_representation_kind =
            :final_basis_density_density_matrix,
        interaction_model = :density_density,
        interaction_treatment = ham_adapter.interaction_treatment,
        interaction_matrix = interaction,
        pre_final_pair_matrix = Float64[],
        final_to_pre_final_coefficients = Float64[],
        pre_final_weights = Float64[],
        support_weights = Float64[],
        support_raw_pair_numerator = Float64[],
        density_gauge = :not_applicable,
        raw_pair_factor_convention = :not_applicable,
    )
    payload.jld2_values["routes/white_lindsey/validation"] = (;
        status,
        blocker = route_blocker,
        overlap_shape = size(overlap),
        one_body_shape = size(one_body),
        two_body_shape = size(interaction),
        expected_shape,
        overlap_finite = all(isfinite, overlap),
        one_body_finite,
        two_body_finite,
        final_integral_weights_positive = all(>(0.0), final_weights),
        one_body_symmetry_defect,
        two_body_symmetry_defect,
        nuclear_metadata_status = ham_adapter.nuclear_metadata_status,
        operator_payload_status = ham_adapter.operator_payload_status,
        interaction_status = ham_adapter.interaction_status,
        basis_adapter_status = basis_adapter.status,
        ham_adapter_status = ham_adapter.status,
    )
    payload.jld2_values["routes/white_lindsey/metadata"] = (;
        status,
        blocker = route_blocker,
        route_family = :white_lindsey_low_order,
        private_development_only = true,
        materialization_status = materialization.status,
        materialized_report_kind = materialization.materialized_report_kind,
        atom_growth_probe_status =
            materialization.route_configured_diatomic_atom_growth_materializer_probe_status,
        route_configured_diatomic_atom_growth_probe_consumed =
            materialization.route_configured_diatomic_atom_growth_materializer_probe_consumed,
        route_configured_diatomic_atom_growth_shellification_consumed =
            materialization.route_configured_diatomic_atom_growth_shellification_consumed,
        route_configured_legacy_diatomic_source_consumed =
            materialization.route_configured_legacy_diatomic_source_consumed,
        materializer_backend = wl_route.materialization_inputs.materializer_backend,
        materializer_nside = wl_route.materialization_inputs.materializer_nside,
        low_order_shellization_policy_requested =
            materialization.low_order_shellization_policy_requested,
        low_order_shellization_policy_resolved =
            materialization.low_order_shellization_policy_resolved,
        low_order_shellization_policy_status =
            materialization.low_order_shellization_policy_status,
        interaction_treatment_requested =
            wl_route.materialization_inputs.route_configured_diatomic_ham_interaction_treatment,
        interaction_treatment_consumed = ham_adapter.interaction_treatment,
        supplement_residual_gto_status = :unavailable,
        qiu_white_atom_local_hf_inputs_status = :unavailable,
        correction_egoi_stationary_cusp_status = :unavailable,
        mwg_ida_route_configured_diatomic_ham_status =
            :pending_route_configured_diatomic_mwg_operator_support,
        old_seed_one_center_promoted = false,
    )
    payload.jld2_values["routes/white_lindsey/hf_convention"] = (;
        density_density_hf_convention_status = hf_convention_blocker,
        density_density_hf_convention_blocker = hf_convention_blocker,
    )
    return (;
        jld2_values = payload.jld2_values,
        rows = (payload.rows[1], wl_fingerprint),
        fingerprint_columns = payload.fingerprint_columns,
    )
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
    wl_route = _be2_white_lindsey_atom_growth_route()
    payload = _be2_populate_white_lindsey_route(payload, wl_route)
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
