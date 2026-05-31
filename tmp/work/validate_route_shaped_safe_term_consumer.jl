using LinearAlgebra
using Printf
using GaussletBases

const CCPM = GaussletBases.CartesianContractedParentMetrics

const SAFE_TERMS = (
    :overlap,
    :position_x,
    :position_y,
    :position_z,
    :x2_x,
    :x2_y,
    :x2_z,
    :kinetic,
)

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

function route_fixture(; route_name::Symbol, q::Int, L::Int)
    parent_dims = (q, q, L + 2)
    bundles = GaussletBases._CartesianNestedAxisBundles3D(
        pqs_test_bundle(q),
        pqs_test_bundle(q),
        pqs_test_bundle(L + 2),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    term_coefficients = Float64.(expansion.coefficients)
    left_current = (1:q, 1:q, 1:L)
    left_inner = (2:(q - 1), 2:(q - 1), 2:(L - 1))
    right_current = (1:q, 1:q, 3:(L + 2))
    right_inner = (2:(q - 1), 2:(q - 1), 4:(L + 1))
    left_layer = GaussletBases._nested_projected_q_shell_layer(
        bundles,
        left_current,
        left_inner;
        bond_axis = :z,
        q,
        L,
        term_coefficients,
    )
    right_layer = GaussletBases._nested_projected_q_shell_layer(
        bundles,
        right_current,
        right_inner;
        bond_axis = :z,
        q,
        L,
        term_coefficients,
    )
    left_descriptor =
        GaussletBases._nested_projected_q_shell_staged_unit_descriptor(left_layer)
    right_descriptor =
        GaussletBases._nested_projected_q_shell_staged_unit_descriptor(right_layer)
    source_mode_dims = (q, q, L)
    left_shared = GaussletBases._cartesian_raw_product_box_plan(
        bundles,
        left_descriptor.axis_intervals,
        source_mode_dims;
        enforce_symmetric_odd = false,
    )
    right_shared = GaussletBases._cartesian_raw_product_box_plan(
        bundles,
        right_descriptor.axis_intervals,
        source_mode_dims;
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

    slab_z = (L + 3) ÷ 2
    product_states = NTuple{3,Int}[
        (ix, iy, slab_z) for ix in 1:q for iy in 1:q
    ]
    product_indices = [
        GaussletBases._cartesian_flat_index(state..., parent_dims) for
        state in product_states
    ]
    identity_axis = Matrix{Float64}(I, q, q)
    product_axes = (
        GaussletBases._nested_product_staged_active_axis(1:q, identity_axis),
        GaussletBases._nested_product_staged_active_axis(1:q, identity_axis),
        GaussletBases._nested_product_staged_fixed_axis(slab_z),
    )
    product_axis_indices =
        GaussletBases._nested_product_axis_function_indices(3, 1, q, 2, q)
    product_unit = GaussletBases._CartesianNestedProductStagedByCenterUnit3D(
        :middle_body_product_slab,
        :product_doside,
        1:(q * q),
        product_indices,
        product_states,
        Matrix{Float64}(I, q * q, q * q),
        product_axes,
        product_axis_indices,
        (source = :route_shaped_safe_term_consumer_scaling_probe,),
        (support_count = q * q, retained_count = q * q),
    )
    route_descriptor = CCPM._pqs_pqs_product_safe_term_route_descriptor(
        left_plan,
        right_plan,
        product_unit;
        route_name,
        parent_dims,
        bond_axis = :z,
        metadata = (
            q = q,
            L = L,
            pqs_left_box = left_current,
            pqs_right_box = right_current,
            product_slab_fixed_index = slab_z,
            pqs_source_mode_dims = source_mode_dims,
        ),
        provenance = (source = :route_shaped_safe_term_consumer_scaling_probe,),
    )
    return (
        route_descriptor = route_descriptor,
        metrics = metrics,
        left_plan = left_plan,
        right_plan = right_plan,
        product_unit = product_unit,
    )
end

function unsupported_weights_rejected(route_descriptor, metrics)
    try
        CCPM._pqs_pqs_product_route_shaped_safe_term_consumer(
            route_descriptor,
            metrics;
            terms = (:weights,),
        )
        return false
    catch err
        err isa ArgumentError || rethrow()
        return true
    end
end

function compare_route(; route_name::Symbol, q::Int, L::Int)
    fixture = route_fixture(; route_name, q, L)
    consumer = CCPM._pqs_pqs_product_route_shaped_safe_term_consumer(
        fixture.route_descriptor,
        fixture.metrics;
        terms = SAFE_TERMS,
    )
    shadow = CCPM._pqs_pqs_product_source_box_shadow_blocks(
        fixture.left_plan,
        fixture.right_plan,
        fixture.product_unit,
        fixture.metrics;
        terms = SAFE_TERMS,
    )

    max_full_error = 0.0
    max_component_error = 0.0
    for term in SAFE_TERMS
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
        @assert size(consumer.blocks[term]) ==
                (consumer.retained_dimension, consumer.retained_dimension)
        @assert all(isfinite, consumer.blocks[term])
    end
    @assert consumer.ranges == shadow.ranges
    @assert consumer.retained_dimension == shadow.retained_dimension
    @assert consumer.pair_count == 6
    @assert consumer.term_count == length(SAFE_TERMS)
    @assert length(consumer.retained_units) == 3
    @assert max_full_error == 0.0
    @assert max_component_error == 0.0
    @assert !consumer.performance.dense_raw_source_box_pair_matrix_materialized
    @assert consumer.performance.dense_raw_pair_storage_avoided
    @assert !consumer.diagnostics.packet_adoption
    @assert !consumer.diagnostics.qwhamiltonian_consumes
    @assert consumer.diagnostics.private_shadow_only

    return (
        route_name = route_name,
        q = q,
        L = L,
        parent_dims = fixture.route_descriptor.metadata.parent_dims,
        pqs_source_mode_dims = fixture.route_descriptor.metadata.pqs_source_mode_dims,
        pqs_retained_count = fixture.route_descriptor.unit_summaries[1].retained_count,
        product_retained_count = fixture.route_descriptor.unit_summaries[3].retained_count,
        retained_dimension = consumer.retained_dimension,
        retained_unit_count = length(consumer.retained_units),
        pair_count = consumer.pair_count,
        term_count = consumer.term_count,
        elapsed_seconds = consumer.performance.elapsed_seconds,
        allocated_bytes = consumer.performance.allocated_bytes,
        gc_time_seconds = consumer.performance.gc_time_seconds,
        dense_raw_source_box_pair_matrix_materialized =
            consumer.performance.dense_raw_source_box_pair_matrix_materialized,
        dense_raw_pair_storage_avoided =
            consumer.performance.dense_raw_pair_storage_avoided,
        max_full_error = max_full_error,
        max_component_error = max_component_error,
        unsupported_weights_rejected =
            unsupported_weights_rejected(fixture.route_descriptor, fixture.metrics),
    )
end

function tsv_value(value)
    return replace(string(value), '\t' => ' ', '\n' => ' ')
end

function write_tsv(path::AbstractString, rows)
    headers = (
        :route_name,
        :q,
        :L,
        :parent_dims,
        :pqs_source_mode_dims,
        :pqs_retained_count,
        :product_retained_count,
        :retained_dimension,
        :retained_unit_count,
        :pair_count,
        :term_count,
        :elapsed_seconds,
        :allocated_bytes,
        :gc_time_seconds,
        :dense_raw_source_box_pair_matrix_materialized,
        :dense_raw_pair_storage_avoided,
        :max_full_error,
        :max_component_error,
        :unsupported_weights_rejected,
    )
    open(path, "w") do io
        println(io, join(headers, '\t'))
        for row in rows
            println(io, join((tsv_value(getproperty(row, header)) for header in headers), '\t'))
        end
    end
    return path
end

function main()
    specs = (
        (route_name = :q5_L5_slab5, q = 5, L = 5),
        (route_name = :q5_L7_slab5, q = 5, L = 7),
        (route_name = :q5_L9_slab5, q = 5, L = 9),
        (route_name = :q7_L7_slab7, q = 7, L = 7),
    )
    rows = map(spec -> compare_route(; spec...), specs)
    for row in rows
        @printf(
            "%s dim=%d pairs=%d terms=%d time=%.3fs alloc=%d maxerr=%.3e unsupported_weights_rejected=%s\n",
            row.route_name,
            row.retained_dimension,
            row.pair_count,
            row.term_count,
            row.elapsed_seconds,
            row.allocated_bytes,
            row.max_full_error,
            row.unsupported_weights_rejected,
        )
    end
    tsv_path = joinpath(@__DIR__, "route_shaped_safe_term_scaling.tsv")
    write_tsv(tsv_path, rows)
    println("wrote $(tsv_path)")
end

main()
