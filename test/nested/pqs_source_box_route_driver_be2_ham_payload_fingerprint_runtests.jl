using Test
using LinearAlgebra
using JLD2
using GaussletBases

function _pqs_be2_ham_payload_fingerprint_assembly(;
    probe_parent_axis_construction = false,
)
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
        probe_parent_axis_construction,
        parent_axis_probe_backend = :pgdg_localized_experimental,
        parent_axis_probe_family = :G10, probe_raw_product_box_plans = false,
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
        white_lindsey_route_shape = (:standard_cartesian_units, :low_order_comx_coarsening),
        white_lindsey_mapping_rule = :standard_unit_backbone_mapping_family,
        white_lindsey_nesting_rule = :unit_box_low_order_comx_coarsening,
        white_lindsey_retained_rule = :low_order_unit_comx_retained_basis,
        white_lindsey_operator_rule = :low_order_unit_operator_blocks,
        white_lindsey_benchmark_role = :published_cartesian_baseline_for_pqs_comparison,
    )

    system = GaussletBases.cartesian_system(system_inputs)
    recipe = GaussletBases.cartesian_recipe(route_inputs)
    parent = GaussletBases.cartesian_parent(system, spacing_inputs, probe_inputs, recipe)
    shells = GaussletBases.cartesian_shells(parent, spacing_inputs, recipe)
    units = GaussletBases.cartesian_units(parent, shells, probe_inputs, recipe)
    transforms = GaussletBases.cartesian_transforms(units, recipe)
    pairs = GaussletBases.cartesian_pair_terms(units, transforms, recipe)
    return GaussletBases.cartesian_assembly(parent, shells, units, transforms, pairs, recipe)
end

@testset "Be2 PQS probe-enabled Ham readiness fingerprint" begin
    assembly = _pqs_be2_ham_payload_fingerprint_assembly(
        probe_parent_axis_construction = :auto)

    source_plan = assembly.diatomic_complete_core_shell_source_plan_payload.source_plan
    final_basis = assembly.diatomic_complete_core_shell_final_basis_payload.final_basis
    @test source_plan.object_kind == :pqs_diatomic_complete_core_shell_source_plan
    @test source_plan.support_order == (:product, :pqs_left, :pqs_right)
    @test source_plan.route_retained_order == (:pqs_left, :pqs_right, :product)
    @test final_basis.final_retained_count == 221
    @test final_basis.support_row_order == :core_then_shell

    h1_payload = assembly.diatomic_complete_core_shell_h1_payload

    h1_matrix = h1_payload.final_hamiltonian.hamiltonian_matrix
    @test size(h1_matrix) == (221, 221)
    @test all(isfinite, h1_matrix)
    @test norm(h1_matrix - h1_matrix') <= 1.0e-10
    @test isapprox(h1_payload.summary.lowest_energy, -0.27746109235228694;
        atol = 1.0e-12, rtol = 0.0)

    payload =
        GaussletBases._pqs_source_box_route_driver_be2_cr2_inspection_bundle_payload(
            assembly)
    mktempdir() do dir
        jld2_path = joinpath(dir, "be2_cr2_inspection.jld2")
        tsv_path = joinpath(dir, "be2_cr2_inspection.tsv")
        GaussletBases._pqs_source_box_route_driver_write_be2_cr2_inspection_bundle(
            jld2_path, tsv_path, payload)
        @test isfile(jld2_path)
        @test isfile(tsv_path)
        jldopen(jld2_path, "r") do file
            @test String(file["schema/name"]) ==
                  "be2_wl_pqs_handoff_inspection_bundle"
            @test file["schema/version"] == 1
            @test all(route -> route in keys(file["routes"]),
                ("pqs_source_box", "white_lindsey"))
            @test size(file["routes/pqs_source_box/one_body/hamiltonian"]) == (221, 221)
            @test size(file["routes/pqs_source_box/two_body/pre_final_pair_matrix"]) ==
                  (221, 221)
            @test length(file["routes/pqs_source_box/two_body/support_weights"]) == 275
            @test Bool(file["routes/pqs_source_box/readiness/cr2_read_only_inspector_ready"])
            @test !Bool(file["routes/pqs_source_box/readiness/cr2_solver_ready"])
            @test !Bool(file["routes/pqs_source_box/readiness/cr2_export_ready"])
            @test String(file["routes/white_lindsey/route/status"]) ==
                  "unavailable"
            @test all(
                family -> family in keys(file["routes/white_lindsey"]),
                ("system", "final_basis", "one_body", "two_body", "validation"),
            )
        end
        tsv_lines = readlines(tsv_path)
        @test String.(split(tsv_lines[1], '\t')) == collect(String.(payload.fingerprint_columns))
        pqs_row = split(tsv_lines[2], '\t')
        wl_row = split(tsv_lines[3], '\t')
        pqs = Dict(Symbol(key) => value for (key, value) in
                   zip(payload.fingerprint_columns, pqs_row))
        wl = Dict(Symbol(key) => value for (key, value) in
                  zip(payload.fingerprint_columns, wl_row))
        expected = (;
            route_label = "pqs_source_box", final_dimension = "221",
            support_weight_count = "275",
            density_gauge = "pre_final_localized_positive_weight",
            raw_pair_factor_convention = "raw_numerator",
            cr2_read_only_inspector_ready = "true", cr2_solver_ready = "false",
            cr2_export_ready = "false", nuclear_repulsion = "4.0",
            electron_count = "8", spin_sector = "closed_shell_singlet",
        )
        @test all(pqs[key] == value for (key, value) in pairs(expected))
        @test (wl[:route_label], wl[:status]) == ("white_lindsey", "unavailable")
    end
end
