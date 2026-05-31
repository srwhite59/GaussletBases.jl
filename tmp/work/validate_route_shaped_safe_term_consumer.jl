using LinearAlgebra
using GaussletBases

const CCPM = GaussletBases.CartesianContractedParentMetrics

function pqs_test_bundle(count::Int)
    xmax = 8.0
    tail = 10.0
    endpoint = (count - 1) / 2
    basis = build_basis(MappedUniformBasisSpec(:G10;
        count,
        mapping = AsinhMapping(
            a = 0.25,
            s = asinh(xmax / 0.25) / (endpoint - xmax / tail),
            tail_spacing = tail,
        ),
        reference_spacing = 1.0,
    ))
    expansion = coulomb_gaussian_expansion(doacc = false)
    return GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :pgdg_localized_experimental,
        refinement_levels = 0,
    )
end

function pqs_axis_metrics(bundles)
    pgdg_x = GaussletBases._nested_axis_pgdg(bundles, :x)
    pgdg_y = GaussletBases._nested_axis_pgdg(bundles, :y)
    pgdg_z = GaussletBases._nested_axis_pgdg(bundles, :z)
    return (
        x = (
            overlap = pgdg_x.overlap,
            position = pgdg_x.position,
            x2 = pgdg_x.x2,
            weights = pgdg_x.weights,
            centers = pgdg_x.centers,
            kinetic = pgdg_x.kinetic,
            source = :nested_pgdg_axis,
        ),
        y = (
            overlap = pgdg_y.overlap,
            position = pgdg_y.position,
            x2 = pgdg_y.x2,
            weights = pgdg_y.weights,
            centers = pgdg_y.centers,
            kinetic = pgdg_y.kinetic,
            source = :nested_pgdg_axis,
        ),
        z = (
            overlap = pgdg_z.overlap,
            position = pgdg_z.position,
            x2 = pgdg_z.x2,
            weights = pgdg_z.weights,
            centers = pgdg_z.centers,
            kinetic = pgdg_z.kinetic,
            source = :nested_pgdg_axis,
        ),
    )
end

function route_shaped_fixture()
    bundle5 = pqs_test_bundle(5)
    bundle7 = pqs_test_bundle(7)
    bundles = GaussletBases._CartesianNestedAxisBundles3D(bundle5, bundle5, bundle7)
    expansion = coulomb_gaussian_expansion(doacc = false)
    term_coefficients = Float64.(expansion.coefficients)
    left_layer = GaussletBases._nested_projected_q_shell_layer(
        bundles,
        (1:5, 1:5, 1:5),
        (2:4, 2:4, 2:4);
        bond_axis = :z,
        q = 5,
        L = 5,
        term_coefficients,
    )
    right_layer = GaussletBases._nested_projected_q_shell_layer(
        bundles,
        (1:5, 1:5, 3:7),
        (2:4, 2:4, 4:6);
        bond_axis = :z,
        q = 5,
        L = 5,
        term_coefficients,
    )
    left_descriptor =
        GaussletBases._nested_projected_q_shell_staged_unit_descriptor(left_layer)
    right_descriptor =
        GaussletBases._nested_projected_q_shell_staged_unit_descriptor(right_layer)
    left_shared = GaussletBases._cartesian_raw_product_box_plan(
        bundles,
        left_descriptor.axis_intervals,
        (5, 5, 5);
        enforce_symmetric_odd = false,
    )
    right_shared = GaussletBases._cartesian_raw_product_box_plan(
        bundles,
        right_descriptor.axis_intervals,
        (5, 5, 5);
        enforce_symmetric_odd = false,
    )
    metrics = pqs_axis_metrics(bundles)
    left_plan = CCPM._pqs_raw_product_box_plan(
        left_descriptor,
        left_shared,
        metrics,
    )
    right_plan = CCPM._pqs_raw_product_box_plan(
        right_descriptor,
        right_shared,
        metrics,
    )

    dims = (5, 5, 7)
    product_states = NTuple{3,Int}[
        (ix, iy, 4) for ix in 1:5 for iy in 1:5
    ]
    product_indices = [
        GaussletBases._cartesian_flat_index(state..., dims) for
        state in product_states
    ]
    identity_axis = Matrix{Float64}(I, 5, 5)
    product_axes = (
        GaussletBases._nested_product_staged_active_axis(1:5, identity_axis),
        GaussletBases._nested_product_staged_active_axis(1:5, identity_axis),
        GaussletBases._nested_product_staged_fixed_axis(4),
    )
    product_axis_indices =
        GaussletBases._nested_product_axis_function_indices(3, 1, 5, 2, 5)
    product_unit = GaussletBases._CartesianNestedProductStagedByCenterUnit3D(
        :middle_body_product_slab,
        :product_doside,
        1:25,
        product_indices,
        product_states,
        Matrix{Float64}(I, 25, 25),
        product_axes,
        product_axis_indices,
        (source = :route_shaped_safe_term_consumer_probe,),
        (support_count = 25, retained_count = 25),
    )
    route_units = (
        route_kind = :homonuclear_pqs_product_source_box_safe_term_fixture,
        units = (
            pqs_left = left_plan,
            pqs_right = right_plan,
            product = product_unit,
        ),
        roles = (:pqs_left, :pqs_right, :product),
        metadata = (
            parent_dims = dims,
            bond_axis = :z,
            pqs_left_box = (1:5, 1:5, 1:5),
            pqs_right_box = (1:5, 1:5, 3:7),
            product_slab_fixed_index = 4,
            pqs_source_mode_dims = (5, 5, 5),
        ),
        provenance = (source = :route_shaped_safe_term_consumer_probe,),
    )
    return (route_units = route_units, metrics = metrics)
end

function main()
    fixture = route_shaped_fixture()
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
    consumer = CCPM._pqs_pqs_product_route_shaped_safe_term_consumer(
        fixture.route_units,
        fixture.metrics;
        terms,
    )
    shadow = CCPM._pqs_pqs_product_source_box_shadow_blocks(
        fixture.route_units.units.pqs_left,
        fixture.route_units.units.pqs_right,
        fixture.route_units.units.product,
        fixture.metrics;
        terms,
    )

    max_full_error = 0.0
    max_component_error = 0.0
    for term in terms
        max_full_error = max(
            max_full_error,
            norm(consumer.blocks[term] - shadow.blocks[term], Inf),
        )
        for key in keys(consumer.component_blocks[term])
            max_component_error = max(
                max_component_error,
                norm(
                    getproperty(consumer.component_blocks[term], key) -
                    getproperty(shadow.component_blocks[term], key),
                    Inf,
                ),
            )
        end
        @assert size(consumer.blocks[term]) == (221, 221)
        @assert all(isfinite, consumer.blocks[term])
    end
    @assert consumer.retained_dimension == 221
    @assert consumer.pair_count == 6
    @assert consumer.term_count == 8
    @assert max_full_error == 0.0
    @assert max_component_error == 0.0
    try
        CCPM._pqs_pqs_product_route_shaped_safe_term_consumer(
            fixture.route_units,
            fixture.metrics;
            terms = (:weights,),
        )
        error("unsupported :weights term did not throw")
    catch err
        err isa ArgumentError || rethrow()
    end

    println((
        retained_dimension = consumer.retained_dimension,
        pair_count = consumer.pair_count,
        term_count = consumer.term_count,
        elapsed_seconds = consumer.performance.elapsed_seconds,
        allocated_bytes = consumer.performance.allocated_bytes,
        gc_time_seconds = consumer.performance.gc_time_seconds,
        max_full_error = max_full_error,
        max_component_error = max_component_error,
        private_shadow_only = consumer.diagnostics.private_shadow_only,
        packet_adoption = consumer.diagnostics.packet_adoption,
        qwhamiltonian_consumes = consumer.diagnostics.qwhamiltonian_consumes,
    ))
end

main()
