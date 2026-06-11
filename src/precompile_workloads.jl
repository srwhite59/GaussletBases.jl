# Narrow precompile workloads for production acceptance routes.
#
# These workloads compile existing route code only. They must not define route
# behavior, introduce fallbacks, or become scientific acceptance checks.

function _precompile_decomposed_wl_one_body_context()
    expansion = coulomb_gaussian_expansion(doacc = false)
    seed_report = _white_lindsey_low_order_materialized_seed_report(
        parent_side_count = 7,
        Z = 2.0,
        d = 0.2,
        tail_spacing = 10.0,
    )
    doside_source_1d = _mapped_ordinary_gausslet_1d_bundle(
        seed_report.fixture.basis;
        exponents = expansion.exponents,
        center = 0.0,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    parent_axis_bundle_object = (;
        x = doside_source_1d,
        y = doside_source_1d,
        z = doside_source_1d,
    )
    pgdg = doside_source_1d.pgdg_intermediate
    parent_axis_counts = ntuple(_ -> seed_report.fixture.parent_side_count, 3)
    inventory = CartesianPairBlockMaterialization.white_lindsey_decomposed_unit_pair_inventory(
        seed_report;
        metadata = (;
            parent_axis_counts,
            parent_axis_bundle_object,
        ),
    )
    center_record = (;
        center_key = :precompile_helium_nucleus,
        center_index = 1,
        nuclear_charge = 2.0,
        location = (0.0, 0.0, 0.0),
    )
    return (;
        expansion,
        parent_axis_counts,
        parent_axis_bundle_object,
        overlap_1d = (; x = pgdg.overlap, y = pgdg.overlap, z = pgdg.overlap),
        kinetic_1d = (; x = pgdg.kinetic, y = pgdg.kinetic, z = pgdg.kinetic),
        inventory,
        center_record,
    )
end

function _precompile_decomposed_wl_one_body_route!(ctx)
    CartesianPairBlockMaterialization.route_global_decomposed_wl_overlap_matrix(
        ctx.inventory;
        parent_axis_counts = ctx.parent_axis_counts,
        parent_axis_bundle_object = ctx.parent_axis_bundle_object,
        overlap_1d = ctx.overlap_1d,
    )
    CartesianPairBlockMaterialization.route_global_decomposed_wl_kinetic_matrix(
        ctx.inventory;
        parent_axis_counts = ctx.parent_axis_counts,
        parent_axis_bundle_object = ctx.parent_axis_bundle_object,
        overlap_1d = ctx.overlap_1d,
        kinetic_1d = ctx.kinetic_1d,
    )
    CartesianPairBlockMaterialization.route_global_electron_nuclear_by_center_matrices(
        ctx.inventory;
        parent_axis_counts = ctx.parent_axis_counts,
        parent_axis_bundle_object = ctx.parent_axis_bundle_object,
        coulomb_expansion = ctx.expansion,
        center_records = (ctx.center_record,),
    )
    return nothing
end

@setup_workload begin
    old_live_timing = timing_live_enabled()
    set_timing_live!(false)
    try
        ctx = _precompile_decomposed_wl_one_body_context()
        @compile_workload begin
            _precompile_decomposed_wl_one_body_route!(ctx)
        end
    finally
        set_timing_live!(old_live_timing)
    end
end
