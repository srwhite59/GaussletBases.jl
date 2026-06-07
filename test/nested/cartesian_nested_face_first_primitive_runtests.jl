# Integration/slow test. Do not include in default nested runner.

@testset "Cartesian nested face first primitive" begin
    function _fixed_a_nested_test_basis(count::Int; a::Float64 = 0.25, xmax::Float64 = 10.0, tail_spacing::Float64 = 10.0)
        endpoint = (count - 1) / 2
        s = asinh(xmax / a) / (endpoint - xmax / tail_spacing)
        basis = build_basis(MappedUniformBasisSpec(:G10;
            count,
            mapping = AsinhMapping(; a, s, tail_spacing),
            reference_spacing = 1.0,
        ))
        return basis, s
    end

    expansion = coulomb_gaussian_expansion(doacc = false)
    basis, s = _fixed_a_nested_test_basis(13)
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    pgdg = bundle.pgdg_intermediate
    interval = 2:(length(basis) - 1)
    side = GaussletBases._nested_doside_1d(bundle, interval, 4)
    pgdg_exact_bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :pgdg_localized_experimental,
        refinement_levels = 0,
    )
    source_axis = GaussletBases._cartesian_source_box_axis_transform(
        pgdg_exact_bundle,
        interval,
        4;
        axis = :x,
        enforce_symmetric_odd = false,
    )
    source_axis_plan = GaussletBases._cartesian_source_box_axis_transform_plan(
        GaussletBases._CartesianNestedAxisBundles3D(
            pgdg_exact_bundle,
            pgdg_exact_bundle,
            pgdg_exact_bundle,
        ),
        (interval, interval, interval),
        (4, 4, 5);
        enforce_symmetric_odd = false,
    )
    pgdg_exact_bundles = GaussletBases._CartesianNestedAxisBundles3D(
        pgdg_exact_bundle,
        pgdg_exact_bundle,
        pgdg_exact_bundle,
    )
    cubic_raw_box_plan = GaussletBases._cartesian_raw_product_box_plan(
        pgdg_exact_bundles,
        (interval, interval, interval),
        (5, 5, 5);
        enforce_symmetric_odd = false,
    )
    rectangular_direct_axis_plan =
        GaussletBases._cartesian_source_box_axis_transform_plan(
            pgdg_exact_bundles,
            (interval, interval, interval),
            (5, 5, 7);
            enforce_symmetric_odd = false,
        )
    rectangular_raw_box_plan = GaussletBases._cartesian_raw_product_box_plan(
        pgdg_exact_bundles,
        (interval, interval, interval),
        (5, 5, 7);
        enforce_symmetric_odd = false,
    )

    @test s > 0.0
    @test side isa GaussletBases._CartesianNestedDoSide1D
    @test side.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test side.interval == interval
    @test side.retained_count == 3
    @test size(side.local_coefficients) == (length(interval), 3)
    @test size(side.coefficient_matrix) == (length(basis), 3)
    @test maximum(abs.(side.coefficient_matrix[1:(first(interval) - 1), :])) == 0.0
    @test maximum(abs.(side.coefficient_matrix[(last(interval) + 1):end, :])) == 0.0
    @test norm(transpose(side.local_coefficients) * side.local_overlap * side.local_coefficients - I, Inf) < 1.0e-10
    @test norm(transpose(side.coefficient_matrix) * pgdg.overlap * side.coefficient_matrix - I, Inf) < 1.0e-10
    @test issorted(side.localized_centers)
    @test length(side.localized_weights) == 3
    @test any(abs.(side.localized_centers) .< 1.0e-10)
    @test source_axis.object_kind == :cartesian_source_box_axis_transform_1d
    @test source_axis.axis == :x
    @test source_axis.interval == interval
    @test source_axis.source_mode_dim_requested == 4
    @test source_axis.source_mode_dim == 4
    @test !source_axis.source_mode_dim_adjusted
    @test source_axis.integration_contract == :pgdg_exact
    @test source_axis.integration_contract_label == "pgdg-exact"
    @test source_axis.diagnostics.pgdg_backend == :pgdg_localized_experimental
    @test source_axis.diagnostics.exact_with_respect_to_pgdg_proxy_basis
    @test !source_axis.diagnostics.numerical_reference_fallback
    @test source_axis.diagnostics.coefficient_overlap_error < 1.0e-10
    @test source_axis_plan.object_kind ==
          :cartesian_source_box_axis_transform_plan_3d
    @test source_axis_plan.source_box == (interval, interval, interval)
    @test source_axis_plan.source_mode_dims == (4, 4, 5)
    @test source_axis_plan.source_mode_count == 80
    @test !source_axis_plan.diagnostics.source_mode_dims_adjusted
    @test source_axis_plan.diagnostics.integration_contract == :pgdg_exact
    @test source_axis_plan.diagnostics.max_axis_overlap_error < 1.0e-10
    @test cubic_raw_box_plan.object_kind == :cartesian_raw_product_box_plan_3d
    @test cubic_raw_box_plan.source_box == (interval, interval, interval)
    @test cubic_raw_box_plan.axis_intervals == (interval, interval, interval)
    @test cubic_raw_box_plan.source_mode_dims == (5, 5, 5)
    @test cubic_raw_box_plan.source_mode_count == 125
    @test length(cubic_raw_box_plan.source_mode_indices) == 125
    @test first(cubic_raw_box_plan.source_mode_indices) == (1, 1, 1)
    @test cubic_raw_box_plan.source_mode_indices[2] == (1, 1, 2)
    @test cubic_raw_box_plan.source_mode_indices[6] == (1, 2, 1)
    @test last(cubic_raw_box_plan.source_mode_indices) == (5, 5, 5)
    @test cubic_raw_box_plan.source_mode_column_indices == collect(1:125)
    @test cubic_raw_box_plan.source_mode_ordering == :x_major_y_major_z_fast
    @test all(
        axis -> size(cubic_raw_box_plan.axis_local_coefficients[axis]) ==
                (length(interval), 5),
        1:3,
    )
    @test cubic_raw_box_plan.axis_transform_plan.source_mode_dims == (5, 5, 5)
    @test cubic_raw_box_plan.diagnostics.source_mode_dims_are_total_lengths
    @test cubic_raw_box_plan.diagnostics.deterministic_given_box_and_dims
    @test cubic_raw_box_plan.diagnostics.integration_contract == :pgdg_exact
    @test !cubic_raw_box_plan.diagnostics.numerical_reference_fallback
    @test !cubic_raw_box_plan.diagnostics.retained_rule_attached
    @test !cubic_raw_box_plan.diagnostics.packet_adoption
    @test cubic_raw_box_plan.diagnostics.max_axis_overlap_error < 1.0e-10
    @test cubic_raw_box_plan.diagnostics.source_product_modes_orthogonal
    @test rectangular_raw_box_plan.source_mode_dims == (5, 5, 7)
    @test rectangular_raw_box_plan.source_mode_count == 175
    @test length(rectangular_raw_box_plan.source_mode_indices) == 175
    @test rectangular_raw_box_plan.source_mode_indices[7] == (1, 1, 7)
    @test rectangular_raw_box_plan.source_mode_indices[8] == (1, 2, 1)
    @test last(rectangular_raw_box_plan.source_mode_indices) == (5, 5, 7)
    @test all(
        axis -> rectangular_raw_box_plan.axis_transform_plan.axes[axis].source_mode_dim ==
                rectangular_direct_axis_plan.axes[axis].source_mode_dim,
        1:3,
    )
    @test all(
        axis -> rectangular_raw_box_plan.axis_local_coefficients[axis] ≈
                rectangular_direct_axis_plan.axes[axis].local_coefficients,
        1:3,
    )
    @test rectangular_raw_box_plan.diagnostics.max_axis_overlap_error < 1.0e-10
    @test rectangular_raw_box_plan.diagnostics.source_product_modes_orthogonal

    face_lo = GaussletBases._nested_xy_face_product(
        pgdg,
        interval,
        interval,
        1;
        retain_x = 4,
        retain_y = 3,
    )
    face_hi = GaussletBases._nested_xy_face_product(
        pgdg,
        interval,
        interval,
        length(basis);
        retain_x = 4,
        retain_y = 3,
    )
    face_overlap = GaussletBases._nested_xy_face_overlap(face_lo, pgdg.overlap)
    face_cross = GaussletBases._nested_xy_face_cross_overlap(face_lo, face_hi, pgdg.overlap)

    @test face_lo isa GaussletBases._CartesianNestedXYFace3D
    @test face_lo.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test size(face_lo.coefficient_matrix) == (length(basis)^3, 9)
    @test length(face_lo.support_indices) == length(interval)^2
    @test isempty(intersect(face_lo.support_indices, face_hi.support_indices))
    @test norm(face_overlap - I, Inf) < 1.0e-10
    @test norm(face_cross, Inf) < 1.0e-10
end
