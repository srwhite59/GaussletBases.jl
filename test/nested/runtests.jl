include("pqs_source_metadata_real_artifact_acceptance_runtests.jl")
include("pqs_component_route_report_adapter_runtests.jl")
include("pqs_standard_source_box_route_setup_runtests.jl")
include("pqs_standard_parent_axis_readiness_runtests.jl")
include("pqs_explicit_core_spacing_parent_axis_probe_runtests.jl")
include("pqs_route_axis_count_selection_runtests.jl")
include("pqs_raw_product_box_plan_probe_runtests.jl")
include("pqs_source_box_route_skeleton_runtests.jl")
include("pqs_source_box_route_driver_report_runtests.jl")
include("pqs_source_box_route_driver_crc_print_line_runtests.jl")
include("cartesian_terminal_shellification_geometry_runtests.jl")
include("cartesian_shellification_plan_runtests.jl")
include("cartesian_ham_builder_one_center_config_smoke_runtests.jl")
include("cartesian_ham_builder_diatomic_config_smoke_runtests.jl")
include("cartesian_route_diatomic_materializer_probe_runtests.jl")
include("white_lindsey_materialized_seed_runtests.jl")

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

@testset "Cartesian nested owned-unit coverage audit" begin
    dense_unit = GaussletBases._CartesianNestedOwnedUnit3D(
        :endcap_a,
        [1, 2],
        [1.0 0.0; 0.0 1.0];
        metadata = (side = :left,),
    )
    sparse_unit = GaussletBases._CartesianNestedOwnedUnit3D(
        :panel_b,
        [3, 4],
        sparse([1, 2], [1, 1], [0.5, 0.5], 2, 1);
        metadata = (side = :right,),
    )
    exact = GaussletBases._nested_owned_unit_coverage_audit(
        [dense_unit, sparse_unit],
        [1, 2, 3, 4],
    )

    @test dense_unit.coefficient_matrix isa Matrix{Float64}
    @test sparse_unit.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test dense_unit.role == :endcap_a
    @test dense_unit.metadata.side == :left
    @test exact.expected_support_count == 4
    @test exact.owned_support_count == 4
    @test exact.duplicate_count == 0
    @test exact.missing_count == 0
    @test exact.outside_count == 0
    @test exact.retained_count == 3
    @test exact.coverage_ok

    duplicate_unit = GaussletBases._CartesianNestedOwnedUnit3D(:duplicate_panel, [2, 3], ones(2, 1))
    duplicate = GaussletBases._nested_owned_unit_coverage_audit(
        [dense_unit, duplicate_unit],
        [1, 2, 3],
    )
    @test duplicate.duplicate_count == 1
    @test duplicate.missing_count == 0
    @test duplicate.outside_count == 0
    @test !duplicate.coverage_ok

    missing = GaussletBases._nested_owned_unit_coverage_audit([dense_unit], [1, 2, 3])
    @test missing.missing_count == 1
    @test !missing.coverage_ok

    outside = GaussletBases._nested_owned_unit_coverage_audit([dense_unit], [1])
    @test outside.outside_count == 1
    @test !outside.coverage_ok

    zero_retained = GaussletBases._CartesianNestedOwnedUnit3D(:empty_panel, [1, 2], zeros(2, 0))
    nonfinite = GaussletBases._CartesianNestedOwnedUnit3D(:bad_panel, [1], [Inf;;])
    @test_throws DimensionMismatch GaussletBases._CartesianNestedOwnedUnit3D(:bad_rows, [1, 2], ones(1, 1))
    @test_throws ArgumentError GaussletBases._nested_owned_unit_coverage_audit([zero_retained], [1, 2])
    @test_throws ArgumentError GaussletBases._nested_owned_unit_coverage_audit([nonfinite], [1])
    @test_throws ArgumentError GaussletBases._nested_owned_unit_coverage_audit([dense_unit], [1, 1, 2])
end

@testset "Cartesian nested projected q-shell local layer" begin
    function _pqs_test_bundle(count::Int)
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
            backend = :numerical_reference,
            refinement_levels = 0,
        )
    end

    function _pqs_one_hot_selector_columns(matrix::AbstractMatrix{<:Real})
        dense = Matrix{Float64}(matrix)
        for column in axes(dense, 2)
            rows = findall(!iszero, dense[:, column])
            length(rows) == 1 || return false
            dense[only(rows), column] == 1.0 || return false
        end
        return true
    end

    function _pqs_boundary_mode_rule_ok(mode::NTuple{3,Int}, sides::NTuple{3,Int})
        return any(axis -> mode[axis] == 1 || mode[axis] == sides[axis], 1:3)
    end

    function _pqs_axis_metrics(bundles)
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

    function _pqs_cross_product_matrix(left_states, right_states, mx, my, mz)
        result = Matrix{Float64}(undef, length(left_states), length(right_states))
        for right_index in eachindex(right_states), left_index in eachindex(left_states)
            ix, iy, iz = left_states[left_index]
            jx, jy, jz = right_states[right_index]
            result[left_index, right_index] =
                mx[ix, jx] * my[iy, jy] * mz[iz, jz]
        end
        return result
    end

    function _pqs_product_box_column_selection_axis_matrices(metrics, term::Symbol)
        term == :overlap && return ((:overlap, :overlap, :overlap),)
        term == :position_x && return ((:position, :overlap, :overlap),)
        term == :position_y && return ((:overlap, :position, :overlap),)
        term == :position_z && return ((:overlap, :overlap, :position),)
        term == :x2_x && return ((:x2, :overlap, :overlap),)
        term == :x2_y && return ((:overlap, :x2, :overlap),)
        term == :x2_z && return ((:overlap, :overlap, :x2),)
        term == :kinetic && return (
            (:kinetic, :overlap, :overlap),
            (:overlap, :kinetic, :overlap),
            (:overlap, :overlap, :kinetic),
        )
        throw(ArgumentError("unsupported PQS product-box column-selection reference term $(term)"))
    end

    function _pqs_product_box_column_selection_coefficients(descriptor)
        source_mode_dims = ntuple(
            axis -> size(descriptor.axis_local_coefficients[axis], 2),
            3,
        )
        states = NTuple{3,Int}[]
        for ix in descriptor.current_box[1],
            iy in descriptor.current_box[2],
            iz in descriptor.current_box[3]
            push!(states, (ix, iy, iz))
        end
        coefficients = zeros(Float64, length(states), prod(source_mode_dims))
        for (row, state) in pairs(states)
            local_rows = ntuple(
                axis -> state[axis] - first(descriptor.axis_intervals[axis]) + 1,
                3,
            )
            column = 0
            for mode_x in 1:source_mode_dims[1],
                mode_y in 1:source_mode_dims[2],
                mode_z in 1:source_mode_dims[3]
                column += 1
                coefficients[row, column] =
                    descriptor.axis_local_coefficients[1][local_rows[1], mode_x] *
                    descriptor.axis_local_coefficients[2][local_rows[2], mode_y] *
                    descriptor.axis_local_coefficients[3][local_rows[3], mode_z]
            end
        end
        return (states = states, coefficients = coefficients)
    end

    function _pqs_source_box_gto_dense_reference(
        descriptor,
        parent_representation,
        probes,
    )
        raw = _pqs_product_box_column_selection_coefficients(descriptor)
        probe_representation = GaussletBases._gto_probe_representation(probes)
        parent_cross = GaussletBases._cartesian_basis_supplement_cross(
            parent_representation,
            probe_representation,
        )
        parent_states = GaussletBases._cartesian_parent_state_basis(parent_representation).states
        parent_row_by_state = Dict{NTuple{3,Int},Int}()
        for (row, state) in pairs(parent_states)
            parent_row_by_state[state] = row
        end
        parent_rows = Int[parent_row_by_state[state] for state in raw.states]
        source_cross =
            transpose(raw.coefficients) * Matrix{Float64}(parent_cross[parent_rows, :])
        return Matrix{Float64}(source_cross[descriptor.boundary_column_indices, :])
    end

    function _pqs_product_box_column_selection_reference(descriptor, metrics; term::Symbol)
        raw = _pqs_product_box_column_selection_coefficients(descriptor)
        selected = descriptor.boundary_column_indices
        result = zeros(Float64, length(selected), length(selected))
        for axis_kinds in _pqs_product_box_column_selection_axis_matrices(
            metrics,
            term,
        )
            axis_matrices = ntuple(
                axis -> getproperty(
                    getproperty(metrics, (:x, :y, :z)[axis]),
                    axis_kinds[axis],
                ),
                3,
            )
            parent_block = _pqs_cross_product_matrix(
                raw.states,
                raw.states,
                axis_matrices[1],
                axis_matrices[2],
                axis_matrices[3],
            )
            # Oracle: O_boundary = P_boundary' * O_product_box * P_boundary.
            source_block =
                transpose(raw.coefficients) * parent_block * raw.coefficients
            result .+= source_block[selected, selected]
        end
        return result
    end

    function _pqs_pqs_source_box_explicit_reference(
        left_descriptor,
        right_descriptor,
        metrics;
        term::Symbol,
    )
        left_raw = _pqs_product_box_column_selection_coefficients(left_descriptor)
        right_raw = _pqs_product_box_column_selection_coefficients(right_descriptor)
        left_selected = left_descriptor.boundary_column_indices
        right_selected = right_descriptor.boundary_column_indices
        result = zeros(Float64, length(left_selected), length(right_selected))
        for axis_kinds in _pqs_product_box_column_selection_axis_matrices(
            metrics,
            term,
        )
            axis_matrices = ntuple(
                axis -> getproperty(
                    getproperty(metrics, (:x, :y, :z)[axis]),
                    axis_kinds[axis],
                ),
                3,
            )
            parent_block = _pqs_cross_product_matrix(
                left_raw.states,
                right_raw.states,
                axis_matrices[1],
                axis_matrices[2],
                axis_matrices[3],
            )
            source_block =
                transpose(left_raw.coefficients) * parent_block * right_raw.coefficients
            result .+= source_block[left_selected, right_selected]
        end
        return result
    end

    function _pqs_product_source_box_explicit_reference(
        descriptor,
        product_unit,
        metrics;
        term::Symbol,
    )
        raw = _pqs_product_box_column_selection_coefficients(descriptor)
        selected = descriptor.boundary_column_indices
        product_coefficients = Matrix{Float64}(product_unit.coefficient_matrix)
        result = zeros(Float64, length(selected), size(product_coefficients, 2))
        for axis_kinds in _pqs_product_box_column_selection_axis_matrices(
            metrics,
            term,
        )
            axis_matrices = ntuple(
                axis -> getproperty(
                    getproperty(metrics, (:x, :y, :z)[axis]),
                    axis_kinds[axis],
                ),
                3,
            )
            parent_block = _pqs_cross_product_matrix(
                raw.states,
                product_unit.support_states,
                axis_matrices[1],
                axis_matrices[2],
                axis_matrices[3],
            )
            source_block =
                transpose(raw.coefficients) * parent_block * product_coefficients
            result .+= source_block[selected, :]
        end
        return result
    end

    function _pqs_product_density_density_explicit_reference(
        descriptor,
        product_unit,
        term_coefficients,
        axis_pair_factor_terms,
    )
        raw = _pqs_product_box_column_selection_coefficients(descriptor)
        selected = descriptor.boundary_column_indices
        product_coefficients = Matrix{Float64}(product_unit.coefficient_matrix)
        result = zeros(Float64, length(selected), size(product_coefficients, 2))
        for term in eachindex(term_coefficients)
            parent_block = _pqs_cross_product_matrix(
                raw.states,
                product_unit.support_states,
                @view(axis_pair_factor_terms.x[term, :, :]),
                @view(axis_pair_factor_terms.y[term, :, :]),
                @view(axis_pair_factor_terms.z[term, :, :]),
            )
            source_block =
                transpose(raw.coefficients) * parent_block * product_coefficients
            result .+= term_coefficients[term] .* source_block[selected, :]
        end
        return result
    end

    function _pqs_pqs_density_density_explicit_reference(
        left_descriptor,
        right_descriptor,
        term_coefficients,
        axis_pair_factor_terms,
    )
        left_raw = _pqs_product_box_column_selection_coefficients(left_descriptor)
        right_raw = _pqs_product_box_column_selection_coefficients(right_descriptor)
        left_selected = left_descriptor.boundary_column_indices
        right_selected = right_descriptor.boundary_column_indices
        result = zeros(Float64, length(left_selected), length(right_selected))
        for term in eachindex(term_coefficients)
            parent_block = _pqs_cross_product_matrix(
                left_raw.states,
                right_raw.states,
                @view(axis_pair_factor_terms.x[term, :, :]),
                @view(axis_pair_factor_terms.y[term, :, :]),
                @view(axis_pair_factor_terms.z[term, :, :]),
            )
            source_block =
                transpose(left_raw.coefficients) * parent_block * right_raw.coefficients
            result .+= term_coefficients[term] .* source_block[left_selected, right_selected]
        end
        return result
    end

    function _product_source_box_explicit_reference(product_unit, metrics; term::Symbol)
        product_coefficients = Matrix{Float64}(product_unit.coefficient_matrix)
        result = zeros(Float64, size(product_coefficients, 2), size(product_coefficients, 2))
        for axis_kinds in _pqs_product_box_column_selection_axis_matrices(
            metrics,
            term,
        )
            axis_matrices = ntuple(
                axis -> getproperty(
                    getproperty(metrics, (:x, :y, :z)[axis]),
                    axis_kinds[axis],
                ),
                3,
            )
            parent_block = _pqs_cross_product_matrix(
                product_unit.support_states,
                product_unit.support_states,
                axis_matrices[1],
                axis_matrices[2],
                axis_matrices[3],
            )
            result .+=
                transpose(product_coefficients) * parent_block * product_coefficients
        end
        return result
    end

    function _check_pqs_product_source_box_shadow_blocks(
        metrics_module,
        descriptor,
        pqs_plan,
        product_unit,
        metrics;
        expected_source_mode_dims,
        expected_retained_count,
        shared_raw_product_box_plan = nothing,
    )
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
        shadow = metrics_module._pqs_product_source_box_shadow_blocks(
            pqs_plan,
            product_unit,
            metrics;
            terms,
        )
        descriptor_shadow = metrics_module._pqs_product_source_box_shadow_blocks(
            descriptor,
            pqs_plan,
            product_unit,
            metrics;
            terms,
        )
        product_count = length(product_unit.column_range)
        pqs_range = 1:expected_retained_count
        product_range = (expected_retained_count + 1):(expected_retained_count + product_count)
        @test shadow.path == :pqs_product_source_box_shadow_blocks
        @test descriptor_shadow.path == shadow.path
        @test shadow.terms == terms
        @test descriptor_shadow.terms == terms
        @test shadow.ranges.pqs == pqs_range
        @test shadow.ranges.product == product_range
        @test descriptor_shadow.ranges == shadow.ranges
        @test shadow.retained_dimension == expected_retained_count + product_count
        @test descriptor_shadow.retained_dimension == shadow.retained_dimension
        @test shadow.diagnostics.raw_plan_first_path
        @test !shadow.diagnostics.descriptor_wrapper
        @test shadow.diagnostics.source_box_shadow_only
        @test !shadow.diagnostics.packet_adoption
        @test !shadow.diagnostics.fixed_block_routing
        @test shadow.diagnostics.pqs_representation == :mode_selected_raw_product_box
        if !isnothing(shared_raw_product_box_plan)
            @test pqs_plan.shared_raw_product_box_plan === shared_raw_product_box_plan
            @test pqs_plan.shared_raw_product_box_plan_used
            @test pqs_plan.source_mode_indices ==
                  shared_raw_product_box_plan.source_mode_indices
            @test pqs_plan.source_mode_ordering ==
                  shared_raw_product_box_plan.source_mode_ordering
            @test shadow.diagnostics.shared_raw_product_box_plan_available
            @test shadow.diagnostics.shared_raw_product_box_plan_used
            @test shadow.diagnostics.source_mode_ordering ==
                  shared_raw_product_box_plan.source_mode_ordering
        else
            @test !pqs_plan.shared_raw_product_box_plan_used
            @test !shadow.diagnostics.shared_raw_product_box_plan_used
        end
        @test shadow.diagnostics.product_doside_retained_transform_used
        @test shadow.diagnostics.reverse_pqs_product_transpose_only
        @test shadow.diagnostics.all_pairs_inventory_private
        @test shadow.diagnostics.pair_inventory_complete_for_units ==
              (:pqs, :product)
        @test shadow.diagnostics.all_pairs_inventory_pair_count == 3
        @test shadow.all_pairs_inventory.object_kind ==
              :pqs_product_source_box_all_pairs_inventory
        @test shadow.all_pairs_inventory.diagnostics.all_pairs_inventory_private
        @test shadow.all_pairs_inventory.diagnostics.pair_inventory_complete_for_units ==
              (:pqs, :product)
        @test length(shadow.all_pairs_inventory.retained_units) == 2
        @test length(shadow.all_pairs_inventory.pair_entries) == 3
        @test shadow.all_pairs_inventory.retained_units[1].unit_key == :pqs
        @test shadow.all_pairs_inventory.retained_units[1].source_family ==
              :mode_selected_raw_product_box
        @test shadow.all_pairs_inventory.retained_units[1].retained_range ==
              pqs_range
        @test shadow.all_pairs_inventory.retained_units[1].source_dimensions ==
              expected_source_mode_dims
        @test shadow.all_pairs_inventory.retained_units[1].retained_count ==
              expected_retained_count
        @test shadow.all_pairs_inventory.retained_units[2].unit_key == :product
        @test shadow.all_pairs_inventory.retained_units[2].source_family ==
              :product_doside
        @test shadow.all_pairs_inventory.retained_units[2].retained_range ==
              product_range
        @test shadow.all_pairs_inventory.retained_units[2].source_dimensions ==
              (2, 2, 1)
        @test shadow.all_pairs_inventory.retained_units[2].retained_count ==
              product_count
        @test map(entry -> entry.pair_key, shadow.all_pairs_inventory.pair_entries) ==
              ((:pqs, :pqs), (:pqs, :product), (:product, :product))
        @test map(entry -> entry.pair_kind, shadow.all_pairs_inventory.pair_entries) ==
              (:pqs_pqs_source_box, :pqs_product_source_box, :product_doside_source_box_pair)
        @test map(entry -> entry.block_helper, shadow.all_pairs_inventory.pair_entries) ==
              (
                  :_pqs_pqs_source_box_reference_blocks,
                  :_pqs_product_source_box_reference_blocks,
                  :_product_doside_source_box_reference_block,
              )
        @test map(entry -> entry.pair_policy, shadow.all_pairs_inventory.pair_entries) ==
              (
                  :source_box_algorithm_available,
                  :source_box_algorithm_available,
                  :source_box_algorithm_available,
              )
        @test all(
            entry -> entry.source_box_algorithmic,
            shadow.all_pairs_inventory.pair_entries,
        )
        @test shadow.all_pairs_inventory.diagnostics.every_pair_uses_source_box_algorithmic_policy
        @test shadow.all_pairs_inventory.diagnostics.source_box_algorithmic_pair_count == 3
        @test !shadow.all_pairs_inventory.diagnostics.packet_adoption
        @test !shadow.all_pairs_inventory.diagnostics.fixed_block_routing
        @test !shadow.all_pairs_inventory.diagnostics.qwhamiltonian_consumes
        @test !shadow.all_pairs_inventory.diagnostics.public_default_consumes
        @test !shadow.all_pairs_inventory.diagnostics.shell_projection_used
        @test !shadow.all_pairs_inventory.diagnostics.lowdin_cleanup_used
        @test !shadow.all_pairs_inventory.diagnostics.support_local_pqs_oracle_used
        @test !shadow.all_pairs_inventory.diagnostics.retained_pqs_weights_used
        @test !shadow.all_pairs_inventory.diagnostics.ida_weight_division_allowed
        @test !shadow.all_pairs_inventory.diagnostics.dense_raw_source_box_pair_matrix_materialized
        @test shadow.all_pairs_inventory.diagnostics.dense_raw_pair_storage_avoided
        @test shadow.diagnostics.pqs_pqs_block_source ==
              :pqs_pqs_source_box_reference_blocks
        @test shadow.diagnostics.pqs_pqs_raw_box_self_reference_compared
        @test shadow.pqs_pqs_reference_blocks.diagnostics.raw_box_self_reference_compared
        @test shadow.pqs_pqs_reference_blocks.diagnostics.pair_plan_reused_for_terms
        @test shadow.diagnostics.pqs_product_block_source ==
              :pqs_product_source_box_reference_blocks
        @test shadow.diagnostics.pqs_product_pair_plan_reused_for_terms
        @test shadow.diagnostics.pair_plan_reused_for_terms
        @test shadow.diagnostics.product_product_block_source ==
              :_product_doside_source_box_reference_block
        @test shadow.diagnostics.every_pair_uses_source_box_algorithmic_policy
        @test shadow.diagnostics.source_box_algorithmic_pair_count == 3
        @test !shadow.diagnostics.dense_raw_source_box_pair_matrix_materialized
        @test shadow.diagnostics.dense_raw_pair_storage_avoided
        @test shadow.diagnostics.retained_block_assembled_directly_from_1d_factors
        @test shadow.diagnostics.source_box_pair_storage_scaling ==
              :one_dimensional_factors_plus_retained_block
        @test shadow.pqs_product_reference_blocks.diagnostics.pair_plan_reused_for_terms
        @test shadow.pqs_product_reference_blocks.diagnostics.pair_plan_reuse_term_count ==
              length(terms)
        @test !shadow.diagnostics.shell_projection_used
        @test !shadow.diagnostics.lowdin_cleanup_used
        @test !shadow.diagnostics.support_coefficient_matrix_used
        @test !shadow.diagnostics.support_local_pqs_oracle_used
        @test shadow.diagnostics.retained_weight_semantics ==
              :not_positive_quadrature_weights
        @test !shadow.diagnostics.retained_pqs_weights_used
        @test !shadow.diagnostics.retained_pqs_weights_positive_checked
        @test !shadow.diagnostics.ida_weight_division_allowed
        @test !shadow.diagnostics.qwhamiltonian_consumes
        @test !shadow.diagnostics.public_default_consumes
        @test !shadow.diagnostics.cr2_science_status_changed
        @test !shadow.diagnostics.local_ecp_gaussian_mwg_implemented
        @test !shadow.diagnostics.generic_retained_unit_framework
        @test pqs_plan.source_mode_dims == expected_source_mode_dims
        for term in terms
            block = shadow.blocks[term]
            descriptor_block = descriptor_shadow.blocks[term]
            components = shadow.component_blocks[term]
            descriptor_components = descriptor_shadow.component_blocks[term]
            expected_pqs = _pqs_product_box_column_selection_reference(
                descriptor,
                metrics;
                term,
            )
            expected_pqs_product = _pqs_product_source_box_explicit_reference(
                descriptor,
                product_unit,
                metrics;
                term,
            )
            expected_product = _product_source_box_explicit_reference(
                product_unit,
                metrics;
                term,
            )
            @test size(block) ==
                  (
                      expected_retained_count + product_count,
                      expected_retained_count + product_count,
                  )
            @test descriptor_block ≈ block atol = 0.0 rtol = 0.0
            @test block[pqs_range, pqs_range] ≈ expected_pqs atol = 1.0e-10 rtol = 1.0e-10
            @test block[pqs_range, product_range] ≈ expected_pqs_product atol = 1.0e-10 rtol = 1.0e-10
            @test block[product_range, pqs_range] ≈ transpose(expected_pqs_product) atol = 1.0e-10 rtol = 1.0e-10
            @test block[product_range, product_range] ≈ expected_product atol = 1.0e-10 rtol = 1.0e-10
            @test components.pqs_pqs ≈ expected_pqs atol = 1.0e-10 rtol = 1.0e-10
            @test components.pqs_product ≈ expected_pqs_product atol = 1.0e-10 rtol = 1.0e-10
            @test components.product_pqs ≈ transpose(expected_pqs_product) atol = 1.0e-10 rtol = 1.0e-10
            @test components.product_product ≈ expected_product atol = 1.0e-10 rtol = 1.0e-10
            @test descriptor_components.pqs_pqs ≈ components.pqs_pqs atol = 0.0 rtol = 0.0
            @test descriptor_components.pqs_product ≈ components.pqs_product atol = 0.0 rtol = 0.0
            @test descriptor_components.product_pqs ≈ components.product_pqs atol = 0.0 rtol = 0.0
            @test descriptor_components.product_product ≈ components.product_product atol = 0.0 rtol = 0.0
        end
        @test_throws ArgumentError metrics_module._pqs_product_source_box_shadow_blocks(
            pqs_plan,
            product_unit,
            metrics;
            terms = (:weights,),
        )
    end

    function _check_pqs_product_source_box_mixed_block(
        metrics_module,
        descriptor,
        pqs_plan,
        product_unit,
        metrics;
        expected_source_mode_dims,
        expected_retained_count,
        shared_raw_product_box_plan = nothing,
    )
        pair_plan = metrics_module._pqs_product_source_box_pair_plan(
            pqs_plan,
            product_unit,
            metrics,
        )
        @test pair_plan.pair_kind == :pqs_product_source_box
        @test pair_plan.left_source_family == :mode_selected_raw_product_box
        @test pair_plan.right_source_family == :product_doside
        @test pair_plan.left_source_dimensions == expected_source_mode_dims
        @test pair_plan.left_source_dimension == prod(expected_source_mode_dims)
        @test pair_plan.right_source_dimensions == (2, 2, 1)
        @test pair_plan.right_source_dimension == 4
        @test pair_plan.left_retained_count == expected_retained_count
        @test pair_plan.right_retained_count == length(product_unit.column_range)
        @test pair_plan.axis_intervals.pqs == pqs_plan.axis_intervals
        @test pair_plan.axis_intervals.product == (1:2, 1:2, 1:1)
        @test length(pair_plan.axis_centers.pqs) == 3
        @test length(pair_plan.axis_centers.product) == 3
        @test pair_plan.pqs_boundary_mode_selector.mode_indices ==
              descriptor.boundary_mode_indices
        @test pair_plan.pqs_boundary_mode_selector.column_indices ==
              descriptor.boundary_column_indices
        @test pair_plan.product_retained_transform.kind == :product_doside
        @test pair_plan.product_retained_transform.object_kind ==
              :product_doside_retained_unit_plan
        @test pair_plan.product_retained_transform.retained_rule_kind ==
              :product_doside
        @test pair_plan.product_retained_unit_plan ===
              pair_plan.product_retained_transform
        @test pair_plan.product_retained_unit_plan.source_axis_intervals ==
              (1:2, 1:2, 1:1)
        @test pair_plan.product_retained_unit_plan.source_axis_lengths ==
              (2, 2, 1)
        @test pair_plan.product_retained_unit_plan.source_dimension == 4
        @test pair_plan.product_retained_unit_plan.retained_axis_counts ==
              (2, 2, 1)
        @test pair_plan.product_retained_unit_plan.retained_count ==
              length(product_unit.column_range)
        @test pair_plan.product_retained_unit_plan.column_range ==
              product_unit.column_range
        @test pair_plan.product_retained_unit_plan.axis_coefficient_matrices[1] ==
              product_unit.axes[1].coefficient_matrix
        @test pair_plan.product_retained_transform.axis_function_indices ==
              product_unit.axis_function_indices
        @test size(pair_plan.one_dimensional_cross_factors.x.overlap) ==
              (expected_source_mode_dims[1], 2)
        @test size(pair_plan.one_dimensional_cross_factors.y.position) ==
              (expected_source_mode_dims[2], 2)
        @test size(pair_plan.one_dimensional_cross_factors.z.kinetic) ==
              (expected_source_mode_dims[3], 1)
        @test pair_plan.diagnostics.private_shadow_only
        @test pair_plan.diagnostics.pqs_representation ==
              :mode_selected_raw_product_box
        @test pair_plan.diagnostics.raw_product_box_plan_used
        @test pair_plan.diagnostics.pqs_raw_product_box_plan_used
        if !isnothing(shared_raw_product_box_plan)
            @test pqs_plan.shared_raw_product_box_plan === shared_raw_product_box_plan
            @test pqs_plan.shared_raw_product_box_plan_used
            @test pqs_plan.source_mode_indices ==
                  shared_raw_product_box_plan.source_mode_indices
            @test pqs_plan.source_mode_ordering ==
                  shared_raw_product_box_plan.source_mode_ordering
            @test pair_plan.diagnostics.shared_raw_product_box_plan_available
            @test pair_plan.diagnostics.shared_raw_product_box_plan_used
            @test pair_plan.diagnostics.source_mode_ordering ==
                  shared_raw_product_box_plan.source_mode_ordering
        else
            @test !pqs_plan.shared_raw_product_box_plan_used
            @test !pair_plan.diagnostics.shared_raw_product_box_plan_used
        end
        @test pair_plan.diagnostics.pqs_boundary_mode_selection_used
        @test pair_plan.diagnostics.product_doside_retained_transform_used
        @test pair_plan.diagnostics.product_doside_retained_unit_plan_used
        @test pair_plan.product_retained_unit_plan.diagnostics.private_adapter
        @test pair_plan.product_retained_unit_plan.diagnostics.metadata_only
        @test !pair_plan.product_retained_unit_plan.diagnostics.coefficients_rebuilt
        @test !pair_plan.product_retained_unit_plan.diagnostics.block_math_changed
        @test !pair_plan.product_retained_unit_plan.diagnostics.packet_adoption
        @test !pair_plan.product_retained_unit_plan.diagnostics.ida_weight_semantics_changed
        @test !pair_plan.product_retained_unit_plan.diagnostics.generic_retained_unit_framework
        @test pair_plan.diagnostics.raw_product_box_operators_use_1d_factors
        @test !pair_plan.diagnostics.shell_projection_used
        @test !pair_plan.diagnostics.lowdin_cleanup_used
        @test !pair_plan.diagnostics.support_coefficient_matrix_used
        @test !pair_plan.diagnostics.support_local_pqs_oracle_used
        @test !pair_plan.diagnostics.retained_pqs_weights_used
        @test !pair_plan.diagnostics.retained_pqs_weights_positive_checked
        @test !pair_plan.diagnostics.ida_weight_division_allowed
        @test !pair_plan.diagnostics.packet_adoption
        @test !pair_plan.diagnostics.fixed_block_routing
        @test !pair_plan.diagnostics.qwhamiltonian_consumes
        @test !pair_plan.diagnostics.public_default_consumes
        @test !pair_plan.diagnostics.cr2_science_status_changed
        @test !pair_plan.diagnostics.local_ecp_gaussian_mwg_implemented
        @test !pair_plan.diagnostics.generic_retained_unit_framework

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
        multi_term_blocks =
            metrics_module._pqs_product_source_box_reference_blocks_from_pair_plan(
                pair_plan;
                terms,
            )
        wrapper_multi_term_blocks =
            metrics_module._pqs_product_source_box_reference_blocks(
                pqs_plan,
                product_unit,
                metrics;
                terms,
            )
        @test multi_term_blocks.path == :pqs_product_source_box_reference_blocks
        @test multi_term_blocks.terms == terms
        @test wrapper_multi_term_blocks.terms == terms
        @test multi_term_blocks.pair_plan === pair_plan
        @test wrapper_multi_term_blocks.pair_plan.pair_kind ==
              :pqs_product_source_box
        @test multi_term_blocks.diagnostics.pair_plan_reused_for_terms
        @test multi_term_blocks.diagnostics.pair_plan_reuse_term_count ==
              length(terms)
        @test !multi_term_blocks.diagnostics.dense_raw_source_box_pair_matrix_materialized
        @test multi_term_blocks.diagnostics.dense_raw_pair_storage_avoided
        @test multi_term_blocks.diagnostics.retained_block_assembled_directly_from_1d_factors
        @test multi_term_blocks.diagnostics.source_box_pair_storage_scaling ==
              :one_dimensional_factors_plus_retained_block
        @test !multi_term_blocks.diagnostics.shell_projection_used
        @test !multi_term_blocks.diagnostics.lowdin_cleanup_used
        @test !multi_term_blocks.diagnostics.support_coefficient_matrix_used
        @test !multi_term_blocks.diagnostics.support_local_pqs_oracle_used
        @test !multi_term_blocks.diagnostics.retained_pqs_weights_used
        @test !multi_term_blocks.diagnostics.ida_weight_division_allowed
        @test !multi_term_blocks.diagnostics.packet_adoption
        @test !multi_term_blocks.diagnostics.fixed_block_routing
        @test !multi_term_blocks.diagnostics.qwhamiltonian_consumes
        @test !multi_term_blocks.diagnostics.public_default_consumes
        @test :weights in multi_term_blocks.diagnostics.unsupported_terms
        for term in terms
            source_box_block =
                metrics_module._pqs_product_source_box_reference_block(
                    pqs_plan,
                    product_unit,
                    metrics;
                    term,
                )
            expected = _pqs_product_source_box_explicit_reference(
                descriptor,
                product_unit,
                metrics;
                term,
            )
            @test source_box_block.path == :pqs_product_source_box_reference
            @test source_box_block.term == term
            @test size(source_box_block.block) ==
                  (expected_retained_count, length(product_unit.column_range))
            @test source_box_block.block ≈ expected atol = 1.0e-10 rtol = 1.0e-10
            @test multi_term_blocks.blocks[term] ≈ source_box_block.block atol = 1.0e-14 rtol = 1.0e-14
            @test wrapper_multi_term_blocks.blocks[term] ≈ source_box_block.block atol = 1.0e-14 rtol = 1.0e-14
            @test source_box_block.diagnostics.raw_product_box_plan_used
            @test source_box_block.diagnostics.pqs_raw_product_box_plan_used
            if !isnothing(shared_raw_product_box_plan)
                @test source_box_block.diagnostics.shared_raw_product_box_plan_used
                @test source_box_block.diagnostics.source_mode_ordering ==
                      shared_raw_product_box_plan.source_mode_ordering
            else
                @test !source_box_block.diagnostics.shared_raw_product_box_plan_used
            end
            @test !source_box_block.diagnostics.shell_projection_used
            @test !source_box_block.diagnostics.lowdin_cleanup_used
            @test !source_box_block.diagnostics.support_coefficient_matrix_used
            @test !source_box_block.diagnostics.support_local_pqs_oracle_used
            @test !source_box_block.diagnostics.retained_pqs_weights_used
            @test !source_box_block.diagnostics.ida_weight_division_allowed
            @test !source_box_block.diagnostics.packet_adoption
            @test !source_box_block.diagnostics.fixed_block_routing
            @test !source_box_block.diagnostics.qwhamiltonian_consumes
            @test !source_box_block.diagnostics.public_default_consumes
            @test !source_box_block.diagnostics.cr2_science_status_changed
            @test !source_box_block.diagnostics.pair_plan_reused_for_terms
            @test !source_box_block.diagnostics.dense_raw_source_box_pair_matrix_materialized
            @test source_box_block.diagnostics.dense_raw_pair_storage_avoided
            @test source_box_block.diagnostics.retained_block_assembled_directly_from_1d_factors
            @test source_box_block.diagnostics.source_box_pair_storage_scaling ==
                  :one_dimensional_factors_plus_retained_block
            @test :weights in source_box_block.diagnostics.unsupported_terms
        end
        @test_throws ArgumentError metrics_module._pqs_product_source_box_reference_blocks_from_pair_plan(
            pair_plan;
            terms = (:weights,),
        )
        @test_throws ArgumentError metrics_module._pqs_product_source_box_reference_blocks(
            pqs_plan,
            product_unit,
            metrics;
            terms = (:weights,),
        )
        @test_throws ArgumentError metrics_module._pqs_product_source_box_reference_block(
            pqs_plan,
            product_unit,
            metrics;
            term = :weights,
        )
    end

    function _check_pqs_source_box_gto_cross_overlap_shadow(
        descriptor,
        parent_representation,
        probes;
        expected_source_mode_dims,
        expected_retained_count,
        shared_raw_product_box_plan = nothing,
    )
        descriptor_shadow = GaussletBases._pqs_source_box_gto_cross_overlap_shadow(
            descriptor,
            parent_representation,
            probes;
            provenance = :pqs_source_box_gto_test_descriptor_adapter,
        )
        raw_plan =
            GaussletBases.CartesianContractedParentMetrics._pqs_raw_product_box_plan(
                descriptor,
            )
        fallback_raw_plan_shadow = GaussletBases._pqs_source_box_gto_cross_overlap_shadow(
            raw_plan,
            parent_representation,
            probes;
            provenance = :pqs_source_box_gto_test_raw_plan,
        )
        shared_raw_plan = isnothing(shared_raw_product_box_plan) ?
                          nothing :
                          GaussletBases.CartesianContractedParentMetrics._pqs_raw_product_box_plan(
                              descriptor,
                              shared_raw_product_box_plan,
                          )
        shared_raw_plan_shadow = isnothing(shared_raw_plan) ?
                                 nothing :
                                 GaussletBases._pqs_source_box_gto_cross_overlap_shadow(
                                     shared_raw_plan,
                                     parent_representation,
                                     probes;
                                     provenance = :pqs_source_box_gto_test_shared_raw_plan,
                                 )
        shadow = isnothing(shared_raw_plan_shadow) ?
                 fallback_raw_plan_shadow :
                 shared_raw_plan_shadow
        expected = _pqs_source_box_gto_dense_reference(
            descriptor,
            parent_representation,
            probes,
        )
        probe_representation = basis_representation(probes)
        @test shadow.path == :pqs_source_box_gto_cross_overlap_shadow
        @test descriptor_shadow.path == shadow.path
        @test size(shadow.cross_overlap) ==
              (expected_retained_count, length(probe_representation.orbitals))
        @test shadow.cross_overlap ≈ expected atol = 1.0e-10 rtol = 1.0e-10
        @test descriptor_shadow.cross_overlap ≈ shadow.cross_overlap atol = 1.0e-10 rtol = 1.0e-10
        @test fallback_raw_plan_shadow.cross_overlap ≈ shadow.cross_overlap atol = 1.0e-10 rtol = 1.0e-10
        if !isnothing(shared_raw_plan_shadow)
            @test shared_raw_plan.shared_raw_product_box_plan ===
                  shared_raw_product_box_plan
            @test shared_raw_plan.shared_raw_product_box_plan_used
            @test shared_raw_plan.source_mode_indices ==
                  shared_raw_product_box_plan.source_mode_indices
            @test shared_raw_plan.source_mode_ordering ==
                  shared_raw_product_box_plan.source_mode_ordering
            @test shared_raw_plan_shadow.cross_overlap ≈
                  shadow.cross_overlap atol = 1.0e-10 rtol = 1.0e-10
            @test shared_raw_plan_shadow.diagnostics.shared_raw_product_box_plan_available
            @test shared_raw_plan_shadow.diagnostics.shared_raw_product_box_plan_used
            @test shared_raw_plan_shadow.diagnostics.source_mode_ordering ==
                  shared_raw_product_box_plan.source_mode_ordering
        end
        @test shadow.probe_representation === probe_representation
        @test descriptor_shadow.probe_representation === probe_representation
        @test shadow.diagnostics.path == :pqs_source_box_gto_cross_overlap_shadow
        @test shadow.diagnostics.raw_plan_first_path
        @test !shadow.diagnostics.descriptor_adapter
        @test shadow.diagnostics.private_shadow_only
        @test shadow.diagnostics.pqs_representation ==
              :mode_selected_raw_product_box
        @test shadow.diagnostics.raw_product_box_plan_used
        @test descriptor_shadow.diagnostics.raw_plan_first_path
        @test fallback_raw_plan_shadow.diagnostics.raw_plan_first_path
        @test fallback_raw_plan_shadow.diagnostics.raw_product_box_plan_used
        if isnothing(shared_raw_product_box_plan)
            @test !shadow.diagnostics.shared_raw_product_box_plan_used
        else
            @test shadow.diagnostics.shared_raw_product_box_plan_used
        end
        @test !descriptor_shadow.diagnostics.shared_raw_product_box_plan_used
        @test !fallback_raw_plan_shadow.diagnostics.shared_raw_product_box_plan_used
        @test shadow.diagnostics.source_mode_dims == expected_source_mode_dims
        @test shadow.diagnostics.boundary_mode_count == expected_retained_count
        @test shadow.diagnostics.gto_count == length(probe_representation.orbitals)
        @test shadow.diagnostics.primitive_axis_overlap_source ==
              :_cartesian_basis_supplement_axis_primitive_cross
        @test shadow.diagnostics.projected_1d_axis_tables_used
        @test shadow.diagnostics.product_box_column_selection_reference
        @test !shadow.diagnostics.final_basis_handoff_adoption
        @test !shadow.diagnostics.shell_projection_used
        @test !shadow.diagnostics.lowdin_cleanup_used
        @test !shadow.diagnostics.support_coefficient_matrix_used
        @test !shadow.diagnostics.support_local_pqs_oracle_used
        @test shadow.diagnostics.retained_weight_semantics ==
              :not_positive_quadrature_weights
        @test !shadow.diagnostics.retained_pqs_weights_used
        @test !shadow.diagnostics.retained_pqs_weights_positive_checked
        @test !shadow.diagnostics.ida_weight_division_allowed
        @test !shadow.diagnostics.packet_adoption
        @test !shadow.diagnostics.fixed_block_routing
        @test !shadow.diagnostics.qwhamiltonian_consumes
        @test !shadow.diagnostics.public_default_consumes
        @test !shadow.diagnostics.cr2_science_status_changed
        @test !shadow.diagnostics.local_ecp_gaussian_mwg_implemented
        @test !shadow.diagnostics.generic_retained_unit_framework
        @test :lowdin_cleanup in shadow.diagnostics.unsupported_terms
        @test shadow.diagnostics.output_finite
        @test isfinite(shadow.diagnostics.max_abs_cross_overlap)
        bad_parent = (parent_kind = :not_a_cartesian_representation,)
        @test_throws MethodError GaussletBases._pqs_source_box_gto_cross_overlap_shadow(
            descriptor,
            bad_parent,
            probes,
        )
    end

    function _check_pqs_raw_product_box_reference(
        metrics_module,
        descriptor,
        metrics;
        expected_source_mode_dims,
        expected_retained_count,
        shared_raw_product_box_plan = nothing,
    )
        plan = metrics_module._pqs_product_box_realization_plan(
            descriptor,
            metrics;
            shared_raw_product_box_plan = shared_raw_product_box_plan,
        )
        raw_plan = isnothing(shared_raw_product_box_plan) ?
                   metrics_module._pqs_raw_product_box_plan(
                       descriptor,
                       metrics,
                   ) :
                   metrics_module._pqs_raw_product_box_plan(
                       descriptor,
                       shared_raw_product_box_plan,
                       metrics,
                   )
        shell_plan = metrics_module._pqs_shell_realization_plan(
            descriptor,
            metrics,
        )
        structural_raw_plan = metrics_module._pqs_raw_product_box_plan(descriptor)
        shared_structural_raw_plan = isnothing(shared_raw_product_box_plan) ?
                                     nothing :
                                     metrics_module._pqs_raw_product_box_plan(
                                         descriptor,
                                         shared_raw_product_box_plan,
                                     )
        @test raw_plan.path == :pqs_raw_product_box_plan
        @test raw_plan.representation == :orthogonal_raw_product_box
        @test raw_plan.source_mode_dims == expected_source_mode_dims
        @test raw_plan.source_mode_count == prod(expected_source_mode_dims)
        @test raw_plan.source_mode_ordering == :x_major_y_major_z_fast
        @test raw_plan.operator_factors_available
        @test raw_plan.source_product_modes_orthogonal
        @test raw_plan.max_1d_source_overlap_error < 1.0e-10
        @test raw_plan.max_product_overlap_error < 1.0e-10
        @test raw_plan.selected_overlap_error < 1.0e-10
        @test !raw_plan.row_projected_shell_support
        @test !raw_plan.lowdin_cleanup_used
        @test raw_plan.boundary_selector.mode_indices ==
              descriptor.boundary_mode_indices
        @test raw_plan.boundary_selector.column_indices ==
              descriptor.boundary_column_indices
        @test structural_raw_plan.path == :pqs_raw_product_box_plan
        @test !structural_raw_plan.operator_factors_available
        @test structural_raw_plan.source_mode_dims == expected_source_mode_dims
        @test structural_raw_plan.source_mode_ordering == :x_major_y_major_z_fast
        @test structural_raw_plan.boundary_selector.selected_count ==
              expected_retained_count
        @test !structural_raw_plan.shared_raw_product_box_plan_available
        @test !structural_raw_plan.shared_raw_product_box_plan_used
        @test structural_raw_plan.diagnostics.shared_raw_product_box_plan_status ==
              :unavailable_descriptor_only_path_has_no_axis_bundles
        @test structural_raw_plan.diagnostics.shared_raw_product_box_plan_unavailable_reason ==
              :descriptor_only_path_has_no_axis_bundles
        if !isnothing(shared_raw_product_box_plan)
            @test shared_raw_product_box_plan.object_kind ==
                  :cartesian_raw_product_box_plan_3d
            @test raw_plan.shared_raw_product_box_plan === shared_raw_product_box_plan
            @test raw_plan.shared_raw_product_box_plan_available
            @test raw_plan.shared_raw_product_box_plan_used
            @test raw_plan.source_mode_indices ==
                  shared_raw_product_box_plan.source_mode_indices
            @test raw_plan.source_mode_ordering ==
                  shared_raw_product_box_plan.source_mode_ordering
            @test raw_plan.diagnostics.shared_raw_product_box_plan_status == :available
            @test isnothing(raw_plan.diagnostics.shared_raw_product_box_plan_unavailable_reason)
            @test raw_plan.diagnostics.shared_raw_product_box_plan_available
            @test raw_plan.diagnostics.shared_raw_product_box_plan_used
            @test !shared_structural_raw_plan.operator_factors_available
            @test shared_structural_raw_plan.shared_raw_product_box_plan ===
                  shared_raw_product_box_plan
            @test shared_structural_raw_plan.shared_raw_product_box_plan_available
            @test shared_structural_raw_plan.shared_raw_product_box_plan_used
            @test shared_structural_raw_plan.source_mode_indices ==
                  shared_raw_product_box_plan.source_mode_indices
            for axis in 1:3
                @test raw_plan.axis_local_coefficients[axis] ≈
                      shared_raw_product_box_plan.axis_local_coefficients[axis]
                @test shared_structural_raw_plan.axis_local_coefficients[axis] ≈
                      shared_raw_product_box_plan.axis_local_coefficients[axis]
            end
        else
            @test !raw_plan.shared_raw_product_box_plan_available
            @test !raw_plan.shared_raw_product_box_plan_used
            @test raw_plan.diagnostics.shared_raw_product_box_plan_status ==
                  :unavailable_descriptor_only_path_has_no_axis_bundles
            @test raw_plan.diagnostics.shared_raw_product_box_plan_unavailable_reason ==
                  :descriptor_only_path_has_no_axis_bundles
        end
        @test shell_plan.path == :pqs_shell_realization_plan
        @test shell_plan.representation == :shell_projection_lowdin_isometry
        @test shell_plan.shell_projection_used
        @test shell_plan.lowdin_cleanup_used
        @test shell_plan.isometry_error < 1.0e-10
        @test shell_plan.isometric
        @test plan.raw_product_box_plan.source_mode_dims ==
              raw_plan.source_mode_dims
        @test plan.raw_product_box_plan.shared_raw_product_box_plan_used ==
              raw_plan.shared_raw_product_box_plan_used
        @test plan.shell_realization_plan.isometry_error ==
              shell_plan.isometry_error
        @test plan.source_box_plan.representation == :orthogonal_raw_product_box
        @test plan.source_box_plan.source_mode_dims == expected_source_mode_dims
        @test plan.source_box_plan.source_mode_count == prod(expected_source_mode_dims)
        @test plan.source_box_plan.source_product_modes_orthogonal
        @test plan.source_box_plan.max_1d_source_overlap_error < 1.0e-10
        @test plan.source_box_plan.max_product_overlap_error < 1.0e-10
        @test plan.source_box_plan.selected_overlap_error < 1.0e-10
        @test !plan.source_box_plan.row_projected_shell_support
        @test !plan.source_box_plan.lowdin_cleanup_used
        @test plan.boundary_selector.mode_indices == descriptor.boundary_mode_indices
        @test plan.boundary_selector.column_indices == descriptor.boundary_column_indices
        @test plan.boundary_selector.selected_count == expected_retained_count
        @test plan.boundary_selector.preserves_orthogonality
        @test size(plan.one_dimensional_operator_factors.x.overlap) ==
              (expected_source_mode_dims[1], expected_source_mode_dims[1])
        @test size(plan.one_dimensional_operator_factors.y.position) ==
              (expected_source_mode_dims[2], expected_source_mode_dims[2])
        @test size(plan.one_dimensional_operator_factors.z.kinetic) ==
              (expected_source_mode_dims[3], expected_source_mode_dims[3])
        @test plan.shell_realization_plan.representation ==
              :shell_projection_lowdin_isometry
        @test plan.shell_realization_plan.shell_projection_used
        @test plan.shell_realization_plan.lowdin_cleanup_used
        @test size(plan.shell_projection_matrix) ==
              (descriptor.support_count, expected_retained_count)
        @test size(plan.lowdin_cleanup) ==
              (expected_retained_count, expected_retained_count)
        @test plan.shell_realization_plan.isometry_error < 1.0e-10
        @test plan.shell_realization_plan.isometric
        @test plan.diagnostics.private_shadow_only
        @test !plan.diagnostics.raw_product_box_stage_lowdin_cleanup_used
        @test plan.diagnostics.raw_product_box_operators_use_1d_factors
        @test !plan.diagnostics.shell_projection_used_for_raw_box_operators
        @test plan.diagnostics.shell_projection_realization_available
        @test plan.diagnostics.shell_projection_realization_requires_lowdin
        @test plan.diagnostics.shell_projection_realization_isometric
        @test plan.diagnostics.retained_functions_live_in_product_box_mode_span
        @test plan.diagnostics.shell_realized_functions_live_in_shell_row_support_subspace
        @test plan.diagnostics.retained_weight_semantics ==
              :not_positive_quadrature_weights
        @test !plan.diagnostics.ida_weight_division_allowed
        @test !plan.diagnostics.packet_adoption
        @test !plan.diagnostics.fixed_block_sidecar_installation
        @test !plan.diagnostics.qwhamiltonian_consumes
        @test !plan.diagnostics.public_default_consumes
        @test !plan.diagnostics.cr2_science_status_changed
        @test !plan.diagnostics.generic_retained_unit_framework

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
        for term in terms
            raw_box_block = metrics_module._pqs_raw_product_box_reference_block(
                raw_plan;
                term,
            )
            descriptor_raw_box_block =
                metrics_module._pqs_raw_product_box_reference_block(
                    descriptor,
                    metrics;
                    term,
                )
            expected = _pqs_product_box_column_selection_reference(
                descriptor,
                metrics;
                term,
            )
            @test raw_box_block.path == :pqs_mode_selected_raw_product_box_reference
            @test raw_box_block.term == term
            @test size(raw_box_block.block) ==
                  (expected_retained_count, expected_retained_count)
            @test raw_box_block.block ≈ expected atol = 1.0e-10 rtol = 1.0e-10
            @test descriptor_raw_box_block.block ≈ raw_box_block.block atol = 1.0e-10 rtol = 1.0e-10
            if term == :overlap
                @test raw_box_block.block ≈ Matrix{Float64}(
                    I,
                    expected_retained_count,
                    expected_retained_count,
                ) atol = 1.0e-10 rtol = 1.0e-10
            end
            @test raw_box_block.diagnostics.pqs_representation ==
                  :mode_selected_raw_product_box
            @test raw_box_block.raw_product_box_plan.source_box_plan_contract ==
                  :RawProductBoxPlan
            @test raw_box_block.raw_product_box_plan.retained_rule_contract ==
                  :RetainedRule
            @test raw_box_block.raw_product_box_plan.retained_rule_kind ==
                  :boundary_comx_product_mode_selection
            @test raw_box_block.raw_product_box_plan.retained_rule_algorithmic
            @test raw_box_block.diagnostics.private_shadow_only
            @test !raw_box_block.diagnostics.production_supported
            @test !raw_box_block.diagnostics.row_projected_shell_support
            @test !raw_box_block.diagnostics.shell_row_projection_used
            @test !raw_box_block.diagnostics.lowdin_cleanup_used
            @test !raw_box_block.diagnostics.raw_product_box_stage_lowdin_cleanup_used
            @test raw_box_block.diagnostics.shell_projection_realization_stage ==
                  :postponed
            @test !raw_box_block.diagnostics.shell_projection_realization_applied
            @test raw_box_block.diagnostics.shell_projection_realization_requires_lowdin
            @test raw_box_block.diagnostics.lowdin_cleanup_scope ==
                  :shell_projection_realization_stage_only
            @test raw_box_block.diagnostics.descriptor_cleanup_transform_ignored
            @test raw_box_block.diagnostics.source_product_modes_orthogonal
            @test raw_box_block.diagnostics.max_1d_source_overlap_error < 1.0e-10
            @test raw_box_block.diagnostics.max_product_overlap_error < 1.0e-10
            @test raw_box_block.diagnostics.selected_overlap_error < 1.0e-10
            @test raw_box_block.diagnostics.overlap_identity_error < 1.0e-10
            @test raw_box_block.diagnostics.boundary_selection_preserves_orthogonality
            @test raw_box_block.diagnostics.product_box_column_selection_reference
            @test raw_box_block.diagnostics.boundary_column_selection_only
            @test raw_box_block.diagnostics.source_mode_dims == expected_source_mode_dims
            @test raw_box_block.diagnostics.boundary_mode_count == expected_retained_count
            @test raw_box_block.diagnostics.retained_count == expected_retained_count
            @test raw_box_block.diagnostics.retained_functions_live_in_product_box_mode_span
            @test !raw_box_block.diagnostics.retained_functions_live_in_shell_row_support_subspace
            @test raw_box_block.diagnostics.retained_columns_have_full_product_box_support
            @test !raw_box_block.diagnostics.support_coefficient_matrix_oracle_used
            @test !raw_box_block.diagnostics.shell_row_support_oracle_used
            @test raw_box_block.diagnostics.retained_weight_semantics ==
                  :not_positive_quadrature_weights
            @test !raw_box_block.diagnostics.retained_pqs_weights_used
            @test !raw_box_block.diagnostics.retained_pqs_weights_positive_checked
            @test !raw_box_block.diagnostics.ida_weight_division_allowed
            @test !raw_box_block.diagnostics.packet_adoption
            @test !raw_box_block.diagnostics.fixed_block_sidecar_installation
            @test !raw_box_block.diagnostics.qwhamiltonian_consumes
            @test !raw_box_block.diagnostics.public_default_consumes
            @test !raw_box_block.diagnostics.cr2_science_status_changed
            @test !raw_box_block.diagnostics.local_ecp_gaussian_mwg_implemented
            @test !raw_box_block.diagnostics.generic_retained_unit_framework
        end
        @test_throws ArgumentError metrics_module._pqs_raw_product_box_reference_block(
            descriptor,
            metrics;
            term = :weights,
        )
    end

    function _check_pqs_pqs_source_box_self_blocks(
        metrics_module,
        pqs_plan,
        metrics;
        expected_source_mode_dims,
        expected_retained_count,
    )
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
        pair_plan = metrics_module._pqs_pqs_source_box_pair_plan(
            pqs_plan,
            pqs_plan,
            metrics,
        )
        @test pair_plan.pair_kind == :pqs_pqs_source_box
        @test pair_plan.object_contract == :SourceBoxPairOperatorPlan
        @test pair_plan.pair_policy == :source_box_algorithm_available
        @test pair_plan.left_source_family == :mode_selected_raw_product_box
        @test pair_plan.right_source_family == :mode_selected_raw_product_box
        @test pair_plan.left_raw_product_box_plan_contract == :RawProductBoxPlan
        @test pair_plan.right_raw_product_box_plan_contract == :RawProductBoxPlan
        @test pair_plan.left_retained_rule_contract == :RetainedRule
        @test pair_plan.right_retained_rule_contract == :RetainedRule
        @test pair_plan.left_retained_rule_kind ==
              :boundary_comx_product_mode_selection
        @test pair_plan.right_retained_rule_kind ==
              :boundary_comx_product_mode_selection
        @test pair_plan.left_source_dimensions == expected_source_mode_dims
        @test pair_plan.right_source_dimensions == expected_source_mode_dims
        @test pair_plan.left_source_dimension == prod(expected_source_mode_dims)
        @test pair_plan.right_source_dimension == prod(expected_source_mode_dims)
        @test pair_plan.left_retained_count == expected_retained_count
        @test pair_plan.right_retained_count == expected_retained_count
        @test pair_plan.axis_intervals.left == pqs_plan.axis_intervals
        @test pair_plan.axis_intervals.right == pqs_plan.axis_intervals
        @test pair_plan.left_raw_product_box_plan.source_mode_dims ==
              pqs_plan.source_mode_dims
        @test pair_plan.right_raw_product_box_plan.source_mode_dims ==
              pqs_plan.source_mode_dims
        @test pair_plan.left_boundary_mode_selector.mode_indices ==
              pqs_plan.boundary_selector.mode_indices
        @test pair_plan.right_boundary_mode_selector.mode_indices ==
              pqs_plan.boundary_selector.mode_indices
        @test pair_plan.supported_terms == terms
        @test size(pair_plan.one_dimensional_cross_factors.x.overlap) ==
              (expected_source_mode_dims[1], expected_source_mode_dims[1])
        @test size(pair_plan.one_dimensional_cross_factors.y.position) ==
              (expected_source_mode_dims[2], expected_source_mode_dims[2])
        @test size(pair_plan.one_dimensional_cross_factors.z.kinetic) ==
              (expected_source_mode_dims[3], expected_source_mode_dims[3])
        @test pair_plan.diagnostics.source == :pqs_pqs_source_box_pair_plan
        @test pair_plan.diagnostics.pair_kind == :pqs_pqs_source_box
        @test pair_plan.diagnostics.source_box_pair_operator_plan_contract ==
              :SourceBoxPairOperatorPlan
        @test pair_plan.diagnostics.pair_policy ==
              :source_box_algorithm_available
        @test pair_plan.diagnostics.algorithmic_pair_policy ==
              :source_box_algorithm_available
        @test pair_plan.diagnostics.left_raw_product_box_plan_contract ==
              :RawProductBoxPlan
        @test pair_plan.diagnostics.right_raw_product_box_plan_contract ==
              :RawProductBoxPlan
        @test pair_plan.diagnostics.left_retained_rule_contract ==
              :RetainedRule
        @test pair_plan.diagnostics.right_retained_rule_contract ==
              :RetainedRule
        @test pair_plan.diagnostics.left_retained_rule_kind ==
              :boundary_comx_product_mode_selection
        @test pair_plan.diagnostics.right_retained_rule_kind ==
              :boundary_comx_product_mode_selection
        @test pair_plan.diagnostics.private_shadow_only
        @test pair_plan.diagnostics.self_same_plan_only
        @test !pair_plan.diagnostics.cross_pqs_inputs_supported
        @test pair_plan.diagnostics.same_raw_product_box_plan
        @test pair_plan.diagnostics.equal_source_mode_dims_required
        @test pair_plan.diagnostics.pqs_representation ==
              :mode_selected_raw_product_box
        @test pair_plan.diagnostics.left_boundary_mode_selection_used
        @test pair_plan.diagnostics.right_boundary_mode_selection_used
        @test pair_plan.diagnostics.raw_product_box_operators_use_1d_factors
        @test !pair_plan.diagnostics.shell_realization_adapter_used
        @test !pair_plan.diagnostics.support_row_adapter_used
        @test !pair_plan.diagnostics.shell_projection_used
        @test !pair_plan.diagnostics.lowdin_cleanup_used
        @test !pair_plan.diagnostics.support_coefficient_matrix_used
        @test !pair_plan.diagnostics.support_local_pqs_oracle_used
        @test !pair_plan.diagnostics.retained_pqs_weights_used
        @test !pair_plan.diagnostics.retained_pqs_weights_positive_checked
        @test pair_plan.diagnostics.retained_weight_semantics ==
              :not_positive_quadrature_weights
        @test !pair_plan.diagnostics.ida_weight_division_allowed
        @test !pair_plan.diagnostics.packet_adoption
        @test !pair_plan.diagnostics.fixed_block_routing
        @test !pair_plan.diagnostics.qwhamiltonian_consumes
        @test !pair_plan.diagnostics.public_default_consumes
        @test !pair_plan.diagnostics.cr2_science_status_changed
        @test !pair_plan.diagnostics.local_ecp_gaussian_mwg_implemented
        @test !pair_plan.diagnostics.generic_retained_unit_framework
        @test !pair_plan.diagnostics.dense_raw_source_box_pair_matrix_materialized
        @test pair_plan.diagnostics.dense_raw_pair_storage_avoided
        @test pair_plan.diagnostics.retained_block_assembled_directly_from_1d_factors
        @test pair_plan.diagnostics.source_box_pair_storage_scaling ==
              :one_dimensional_factors_plus_retained_block

        blocks = metrics_module._pqs_pqs_source_box_reference_blocks_from_pair_plan(
            pair_plan;
            terms,
        )
        wrapper_blocks = metrics_module._pqs_pqs_source_box_reference_blocks(
            pqs_plan,
            pqs_plan,
            metrics;
            terms,
        )
        @test blocks.path == :pqs_pqs_source_box_reference_blocks
        @test blocks.terms == terms
        @test wrapper_blocks.terms == terms
        @test blocks.diagnostics.raw_box_self_reference_compared
        @test blocks.diagnostics.raw_box_self_reference_helper ==
              :_pqs_raw_product_box_reference_block
        @test blocks.diagnostics.validation_reference_contract ==
              :explicit_raw_product_box_boundary_column_selection
        @test blocks.diagnostics.internal_validation_reference_compared
        @test blocks.diagnostics.explicit_raw_product_box_boundary_column_selection_reference_compared
        @test blocks.diagnostics.explicit_raw_product_box_boundary_column_selection_reference_helper ==
              :_pqs_pqs_source_box_explicit_boundary_selection_reference
        @test blocks.diagnostics.explicit_source_box_oracle_tested
        @test !blocks.diagnostics.cross_box_external_raw_product_oracle_required
        @test !blocks.diagnostics.cross_box_external_raw_product_oracle_compared_by_helper
        @test blocks.diagnostics.dense_raw_source_box_pair_matrix_materialized_for_validation
        @test blocks.diagnostics.pair_plan_reused_for_terms
        @test blocks.diagnostics.pair_plan_reuse_term_count == length(terms)
        @test blocks.diagnostics.max_block_error < 1.0e-10
        @test !blocks.diagnostics.dense_raw_source_box_pair_matrix_materialized
        @test blocks.diagnostics.dense_raw_pair_storage_avoided
        @test blocks.diagnostics.retained_block_assembled_directly_from_1d_factors
        @test !blocks.diagnostics.shell_projection_used
        @test !blocks.diagnostics.lowdin_cleanup_used
        @test !blocks.diagnostics.support_local_pqs_oracle_used
        @test !blocks.diagnostics.retained_pqs_weights_used
        @test !blocks.diagnostics.ida_weight_division_allowed
        for term in terms
            expected = metrics_module._pqs_raw_product_box_reference_block(
                pqs_plan;
                term,
            ).block
            single = metrics_module._pqs_pqs_source_box_reference_block(
                pqs_plan,
                pqs_plan,
                metrics;
                term,
            )
            @test size(blocks.blocks[term]) ==
                  (expected_retained_count, expected_retained_count)
            @test blocks.blocks[term] ≈ expected atol = 1.0e-10 rtol = 1.0e-10
            @test blocks.raw_box_reference_blocks[term] ≈ expected atol = 0.0 rtol = 0.0
            @test wrapper_blocks.blocks[term] ≈ blocks.blocks[term] atol = 0.0 rtol = 0.0
            @test blocks.block_errors[term] < 1.0e-10
            @test single.path == :pqs_pqs_source_box_reference
            @test single.block ≈ expected atol = 1.0e-10 rtol = 1.0e-10
            @test single.raw_box_reference_block ≈ expected atol = 0.0 rtol = 0.0
            @test single.block_error < 1.0e-10
            @test !single.diagnostics.pair_plan_reused_for_terms
            if term == :overlap
                @test single.block ≈ Matrix{Float64}(
                    I,
                    expected_retained_count,
                    expected_retained_count,
                ) atol = 1.0e-10 rtol = 1.0e-10
            end
        end
        @test_throws ArgumentError metrics_module._pqs_pqs_source_box_reference_blocks(
            pqs_plan,
            pqs_plan,
            metrics;
            terms = (:weights,),
        )
    end

    function _check_pqs_pqs_source_box_cross_blocks(
        metrics_module,
        left_descriptor,
        right_descriptor,
        left_pqs_plan,
        right_pqs_plan,
        metrics;
        expected_source_mode_dims,
        expected_retained_count,
    )
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
        pair_plan = metrics_module._pqs_pqs_source_box_pair_plan(
            left_pqs_plan,
            right_pqs_plan,
            metrics,
        )
        @test pair_plan.pair_kind == :pqs_pqs_source_box
        @test pair_plan.object_contract == :SourceBoxPairOperatorPlan
        @test pair_plan.pair_policy == :source_box_algorithm_available
        @test pair_plan.left_source_family == :mode_selected_raw_product_box
        @test pair_plan.right_source_family == :mode_selected_raw_product_box
        @test pair_plan.left_raw_product_box_plan_contract == :RawProductBoxPlan
        @test pair_plan.right_raw_product_box_plan_contract == :RawProductBoxPlan
        @test pair_plan.left_retained_rule_contract == :RetainedRule
        @test pair_plan.right_retained_rule_contract == :RetainedRule
        @test pair_plan.left_retained_rule_kind ==
              :boundary_comx_product_mode_selection
        @test pair_plan.right_retained_rule_kind ==
              :boundary_comx_product_mode_selection
        @test pair_plan.left_source_dimensions == expected_source_mode_dims
        @test pair_plan.right_source_dimensions == expected_source_mode_dims
        @test pair_plan.left_source_dimension == prod(expected_source_mode_dims)
        @test pair_plan.right_source_dimension == prod(expected_source_mode_dims)
        @test pair_plan.left_retained_count == expected_retained_count
        @test pair_plan.right_retained_count == expected_retained_count
        @test pair_plan.axis_intervals.left == left_pqs_plan.axis_intervals
        @test pair_plan.axis_intervals.right == right_pqs_plan.axis_intervals
        @test pair_plan.left_raw_product_box_plan.source_mode_dims ==
              left_pqs_plan.source_mode_dims
        @test pair_plan.right_raw_product_box_plan.source_mode_dims ==
              right_pqs_plan.source_mode_dims
        @test pair_plan.left_boundary_mode_selector.mode_indices ==
              left_pqs_plan.boundary_selector.mode_indices
        @test pair_plan.right_boundary_mode_selector.mode_indices ==
              right_pqs_plan.boundary_selector.mode_indices
        @test pair_plan.supported_terms == terms
        @test size(pair_plan.one_dimensional_cross_factors.x.overlap) ==
              (expected_source_mode_dims[1], expected_source_mode_dims[1])
        @test size(pair_plan.one_dimensional_cross_factors.y.position) ==
              (expected_source_mode_dims[2], expected_source_mode_dims[2])
        @test size(pair_plan.one_dimensional_cross_factors.z.kinetic) ==
              (expected_source_mode_dims[3], expected_source_mode_dims[3])
        @test pair_plan.diagnostics.source == :pqs_pqs_source_box_pair_plan
        @test pair_plan.diagnostics.pair_kind == :pqs_pqs_source_box
        @test pair_plan.diagnostics.source_box_pair_operator_plan_contract ==
              :SourceBoxPairOperatorPlan
        @test pair_plan.diagnostics.pair_policy ==
              :source_box_algorithm_available
        @test pair_plan.diagnostics.algorithmic_pair_policy ==
              :source_box_algorithm_available
        @test pair_plan.diagnostics.left_raw_product_box_plan_contract ==
              :RawProductBoxPlan
        @test pair_plan.diagnostics.right_raw_product_box_plan_contract ==
              :RawProductBoxPlan
        @test pair_plan.diagnostics.left_retained_rule_contract ==
              :RetainedRule
        @test pair_plan.diagnostics.right_retained_rule_contract ==
              :RetainedRule
        @test pair_plan.diagnostics.left_retained_rule_kind ==
              :boundary_comx_product_mode_selection
        @test pair_plan.diagnostics.right_retained_rule_kind ==
              :boundary_comx_product_mode_selection
        @test pair_plan.diagnostics.private_shadow_only
        @test !pair_plan.diagnostics.self_same_plan_only
        @test pair_plan.diagnostics.cross_pqs_inputs_supported
        @test !pair_plan.diagnostics.same_raw_product_box_plan
        @test pair_plan.diagnostics.equal_source_mode_dims_required
        @test pair_plan.diagnostics.pqs_representation ==
              :mode_selected_raw_product_box
        @test pair_plan.diagnostics.left_boundary_mode_selection_used
        @test pair_plan.diagnostics.right_boundary_mode_selection_used
        @test pair_plan.diagnostics.raw_product_box_operators_use_1d_factors
        @test !pair_plan.diagnostics.shell_realization_adapter_used
        @test !pair_plan.diagnostics.support_row_adapter_used
        @test !pair_plan.diagnostics.shell_projection_used
        @test !pair_plan.diagnostics.lowdin_cleanup_used
        @test !pair_plan.diagnostics.support_coefficient_matrix_used
        @test !pair_plan.diagnostics.support_local_pqs_oracle_used
        @test !pair_plan.diagnostics.retained_pqs_weights_used
        @test !pair_plan.diagnostics.retained_pqs_weights_positive_checked
        @test pair_plan.diagnostics.retained_weight_semantics ==
              :not_positive_quadrature_weights
        @test !pair_plan.diagnostics.ida_weight_division_allowed
        @test !pair_plan.diagnostics.packet_adoption
        @test !pair_plan.diagnostics.fixed_block_routing
        @test !pair_plan.diagnostics.qwhamiltonian_consumes
        @test !pair_plan.diagnostics.public_default_consumes
        @test !pair_plan.diagnostics.cr2_science_status_changed
        @test !pair_plan.diagnostics.local_ecp_gaussian_mwg_implemented
        @test !pair_plan.diagnostics.generic_retained_unit_framework
        @test !pair_plan.diagnostics.dense_raw_source_box_pair_matrix_materialized
        @test pair_plan.diagnostics.dense_raw_pair_storage_avoided
        @test pair_plan.diagnostics.retained_block_assembled_directly_from_1d_factors
        @test pair_plan.diagnostics.source_box_pair_storage_scaling ==
              :one_dimensional_factors_plus_retained_block

        blocks = metrics_module._pqs_pqs_source_box_reference_blocks_from_pair_plan(
            pair_plan;
            terms,
        )
        wrapper_blocks = metrics_module._pqs_pqs_source_box_reference_blocks(
            left_pqs_plan,
            right_pqs_plan,
            metrics;
            terms,
        )
        reverse_blocks = metrics_module._pqs_pqs_source_box_reference_blocks(
            right_pqs_plan,
            left_pqs_plan,
            metrics;
            terms,
        )
        @test blocks.path == :pqs_pqs_source_box_reference_blocks
        @test blocks.terms == terms
        @test wrapper_blocks.terms == terms
        @test reverse_blocks.terms == terms
        @test !blocks.diagnostics.raw_box_self_reference_compared
        @test isnothing(blocks.diagnostics.raw_box_self_reference_helper)
        @test blocks.diagnostics.source_box_algorithm_formula_available
        @test blocks.diagnostics.validation_reference_contract ==
              :explicit_raw_product_box_boundary_column_selection
        @test blocks.diagnostics.internal_validation_reference_compared
        @test blocks.diagnostics.explicit_raw_product_box_boundary_column_selection_reference_compared
        @test blocks.diagnostics.explicit_raw_product_box_boundary_column_selection_reference_helper ==
              :_pqs_pqs_source_box_explicit_boundary_selection_reference
        @test !blocks.diagnostics.cross_box_external_raw_product_oracle_required
        @test blocks.diagnostics.cross_box_external_raw_product_oracle_compared_by_helper
        @test blocks.diagnostics.pair_plan_reused_for_terms
        @test blocks.diagnostics.pair_plan_reuse_term_count == length(terms)
        @test blocks.diagnostics.max_block_error < 1.0e-10
        @test blocks.diagnostics.explicit_source_box_oracle_tested
        @test blocks.diagnostics.dense_raw_source_box_pair_matrix_materialized_for_validation
        @test !blocks.diagnostics.dense_raw_source_box_pair_matrix_materialized
        @test blocks.diagnostics.dense_raw_pair_storage_avoided
        @test blocks.diagnostics.retained_block_assembled_directly_from_1d_factors
        @test !blocks.diagnostics.shell_projection_used
        @test !blocks.diagnostics.lowdin_cleanup_used
        @test !blocks.diagnostics.support_local_pqs_oracle_used
        @test !blocks.diagnostics.retained_pqs_weights_used
        @test !blocks.diagnostics.ida_weight_division_allowed

        max_oracle_error = 0.0
        max_transpose_error = 0.0
        for term in terms
            expected = _pqs_pqs_source_box_explicit_reference(
                left_descriptor,
                right_descriptor,
                metrics;
                term,
            )
            reverse_expected = _pqs_pqs_source_box_explicit_reference(
                right_descriptor,
                left_descriptor,
                metrics;
                term,
            )
            single = metrics_module._pqs_pqs_source_box_reference_block(
                left_pqs_plan,
                right_pqs_plan,
                metrics;
                term,
            )
            @test size(blocks.blocks[term]) ==
                  (expected_retained_count, expected_retained_count)
            @test blocks.blocks[term] ≈ expected atol = 1.0e-10 rtol = 1.0e-10
            @test wrapper_blocks.blocks[term] ≈ blocks.blocks[term] atol = 0.0 rtol = 0.0
            @test reverse_blocks.blocks[term] ≈ reverse_expected atol = 1.0e-10 rtol = 1.0e-10
            @test reverse_blocks.blocks[term] ≈ transpose(blocks.blocks[term]) atol = 1.0e-10 rtol = 1.0e-10
            @test single.path == :pqs_pqs_source_box_reference
            @test single.block ≈ expected atol = 1.0e-10 rtol = 1.0e-10
            @test single.raw_box_reference_block ≈ expected atol = 1.0e-10 rtol = 1.0e-10
            @test single.block_error < 1.0e-10
            @test !single.diagnostics.pair_plan_reused_for_terms
            max_oracle_error = max(
                max_oracle_error,
                LinearAlgebra.norm(blocks.blocks[term] - expected, Inf),
            )
            max_transpose_error = max(
                max_transpose_error,
                LinearAlgebra.norm(
                    reverse_blocks.blocks[term] - transpose(blocks.blocks[term]),
                    Inf,
                ),
            )
        end
        @test max_oracle_error < 1.0e-10
        @test max_transpose_error < 1.0e-10
        @test_throws ArgumentError metrics_module._pqs_pqs_source_box_reference_blocks(
            left_pqs_plan,
            right_pqs_plan,
            metrics;
            terms = (:weights,),
        )
    end

    function _check_pqs_pqs_product_source_box_shadow_blocks(
        metrics_module,
        left_descriptor,
        right_descriptor,
        left_pqs_plan,
        right_pqs_plan,
        product_unit,
        metrics;
        expected_source_mode_dims,
        expected_retained_count,
    )
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
        shadow = metrics_module._pqs_pqs_product_source_box_shadow_blocks(
            left_pqs_plan,
            right_pqs_plan,
            product_unit,
            metrics;
            terms,
        )
        product_count = length(product_unit.column_range)
        left_range = 1:expected_retained_count
        right_range =
            (last(left_range) + 1):(last(left_range) + expected_retained_count)
        product_range =
            (last(right_range) + 1):(last(right_range) + product_count)
        retained_dimension = expected_retained_count + expected_retained_count + product_count

        @test shadow.path == :pqs_pqs_product_source_box_shadow_blocks
        @test shadow.terms == terms
        @test shadow.ranges.pqs_left == left_range
        @test shadow.ranges.pqs_right == right_range
        @test shadow.ranges.product == product_range
        @test shadow.retained_dimension == retained_dimension
        @test shadow.diagnostics.source_box_shadow_only
        @test shadow.diagnostics.private_shadow_only
        @test shadow.diagnostics.all_pairs_inventory_private
        @test shadow.diagnostics.pair_inventory_complete_for_units ==
              (:pqs_left, :pqs_right, :product)
        @test shadow.diagnostics.all_pairs_inventory_pair_count == 6
        @test shadow.diagnostics.retained_unit_count == 3
        @test !shadow.diagnostics.packet_adoption
        @test !shadow.diagnostics.fixed_block_routing
        @test shadow.diagnostics.pqs_representation == :mode_selected_raw_product_box
        @test shadow.diagnostics.product_doside_retained_transform_used
        @test shadow.diagnostics.pqs_self_block_source ==
              :pqs_pqs_source_box_reference_blocks
        @test shadow.diagnostics.cross_pqs_block_source ==
              :pqs_pqs_source_box_reference_blocks
        @test shadow.diagnostics.pqs_product_block_source ==
              :pqs_product_source_box_reference_blocks
        @test shadow.diagnostics.product_product_block_source ==
              :_product_doside_source_box_reference_block
        @test !shadow.diagnostics.pqs_left_right_raw_box_self_reference_compared
        @test shadow.diagnostics.pqs_left_left_raw_box_self_reference_compared
        @test shadow.diagnostics.pqs_right_right_raw_box_self_reference_compared
        @test shadow.diagnostics.pqs_left_left_explicit_source_box_oracle_tested
        @test shadow.diagnostics.pqs_left_right_explicit_source_box_oracle_tested
        @test shadow.diagnostics.pqs_right_right_explicit_source_box_oracle_tested
        @test shadow.diagnostics.lower_triangular_cross_blocks_transpose_only
        @test shadow.diagnostics.pair_plan_reused_for_terms
        @test shadow.diagnostics.explicit_source_box_oracle_tested
        @test !shadow.diagnostics.pqs_cross_box_external_raw_product_oracle_required
        @test shadow.diagnostics.pqs_cross_box_internal_raw_product_oracle_compared
        @test shadow.diagnostics.every_pair_uses_source_box_algorithmic_policy
        @test shadow.diagnostics.source_box_algorithmic_pair_count == 6
        @test !shadow.diagnostics.shell_projection_used
        @test !shadow.diagnostics.lowdin_cleanup_used
        @test !shadow.diagnostics.support_coefficient_matrix_used
        @test !shadow.diagnostics.support_local_pqs_oracle_used
        @test shadow.diagnostics.retained_weight_semantics ==
              :not_positive_quadrature_weights
        @test !shadow.diagnostics.retained_pqs_weights_used
        @test !shadow.diagnostics.retained_pqs_weights_positive_checked
        @test !shadow.diagnostics.ida_weight_division_allowed
        @test !shadow.diagnostics.qwhamiltonian_consumes
        @test !shadow.diagnostics.public_default_consumes
        @test !shadow.diagnostics.cr2_science_status_changed
        @test !shadow.diagnostics.local_ecp_gaussian_mwg_implemented
        @test !shadow.diagnostics.generic_retained_unit_framework
        @test !shadow.diagnostics.dense_raw_source_box_pair_matrix_materialized
        @test shadow.diagnostics.dense_raw_source_box_pair_matrix_materialized_for_validation
        @test shadow.diagnostics.dense_raw_pair_storage_avoided
        @test shadow.diagnostics.retained_block_assembled_directly_from_1d_factors
        @test shadow.diagnostics.source_box_pair_storage_scaling ==
              :one_dimensional_factors_plus_retained_block

        inventory = shadow.all_pairs_inventory
        @test inventory.object_kind ==
              :pqs_pqs_product_source_box_all_pairs_inventory
        @test inventory.diagnostics.all_pairs_inventory_private
        @test inventory.diagnostics.pair_inventory_complete_for_units ==
              (:pqs_left, :pqs_right, :product)
        @test inventory.diagnostics.retained_unit_count == 3
        @test inventory.diagnostics.upper_triangular_pair_count == 6
        @test inventory.diagnostics.expected_upper_triangular_pair_count == 6
        @test length(inventory.retained_units) == 3
        @test length(inventory.pair_entries) == 6
        @test map(unit -> unit.unit_key, inventory.retained_units) ==
              (:pqs_left, :pqs_right, :product)
        @test map(unit -> unit.retained_range, inventory.retained_units) ==
              (left_range, right_range, product_range)
        @test inventory.retained_units[1].source_dimensions ==
              expected_source_mode_dims
        @test inventory.retained_units[2].source_dimensions ==
              expected_source_mode_dims
        @test inventory.retained_units[1].retained_count ==
              expected_retained_count
        @test inventory.retained_units[2].retained_count ==
              expected_retained_count
        @test inventory.retained_units[3].source_dimensions == (2, 2, 1)
        @test inventory.retained_units[3].retained_count == product_count
        @test map(entry -> entry.pair_key, inventory.pair_entries) ==
              (
                  (:pqs_left, :pqs_left),
                  (:pqs_left, :pqs_right),
                  (:pqs_left, :product),
                  (:pqs_right, :pqs_right),
                  (:pqs_right, :product),
                  (:product, :product),
              )
        @test map(entry -> entry.pair_kind, inventory.pair_entries) ==
              (
                  :pqs_pqs_source_box,
                  :pqs_pqs_source_box,
                  :pqs_product_source_box,
                  :pqs_pqs_source_box,
                  :pqs_product_source_box,
                  :product_doside_source_box_pair,
              )
        @test map(entry -> entry.block_helper, inventory.pair_entries) ==
              (
                  :_pqs_pqs_source_box_reference_blocks,
                  :_pqs_pqs_source_box_reference_blocks,
                  :_pqs_product_source_box_reference_blocks,
                  :_pqs_pqs_source_box_reference_blocks,
                  :_pqs_product_source_box_reference_blocks,
                  :_product_doside_source_box_reference_block,
              )
        @test map(entry -> entry.pair_policy, inventory.pair_entries) ==
              (
                  :source_box_algorithm_available,
                  :source_box_algorithm_available,
                  :source_box_algorithm_available,
                  :source_box_algorithm_available,
                  :source_box_algorithm_available,
                  :source_box_algorithm_available,
              )
        @test all(entry -> entry.source_box_algorithmic, inventory.pair_entries)
        @test inventory.diagnostics.every_pair_uses_source_box_algorithmic_policy
        @test inventory.diagnostics.source_box_algorithmic_pair_count == 6
        @test !inventory.diagnostics.packet_adoption
        @test !inventory.diagnostics.fixed_block_routing
        @test !inventory.diagnostics.qwhamiltonian_consumes
        @test !inventory.diagnostics.public_default_consumes
        @test !inventory.diagnostics.shell_projection_used
        @test !inventory.diagnostics.lowdin_cleanup_used
        @test !inventory.diagnostics.support_local_pqs_oracle_used
        @test !inventory.diagnostics.retained_pqs_weights_used
        @test !inventory.diagnostics.ida_weight_division_allowed
        @test !inventory.diagnostics.dense_raw_source_box_pair_matrix_materialized
        @test inventory.diagnostics.dense_raw_pair_storage_avoided

        @test shadow.pqs_left_left_reference_blocks.diagnostics.raw_box_self_reference_compared
        @test !shadow.pqs_left_right_reference_blocks.diagnostics.raw_box_self_reference_compared
        @test shadow.pqs_right_right_reference_blocks.diagnostics.raw_box_self_reference_compared
        @test shadow.pqs_left_left_reference_blocks.diagnostics.explicit_source_box_oracle_tested
        @test shadow.pqs_left_right_reference_blocks.diagnostics.explicit_source_box_oracle_tested
        @test shadow.pqs_right_right_reference_blocks.diagnostics.explicit_source_box_oracle_tested
        @test shadow.pqs_left_product_reference_blocks.diagnostics.pair_plan_reused_for_terms
        @test shadow.pqs_right_product_reference_blocks.diagnostics.pair_plan_reused_for_terms

        max_component_error = 0.0
        for term in terms
            block = shadow.blocks[term]
            components = shadow.component_blocks[term]
            expected_left_left = _pqs_product_box_column_selection_reference(
                left_descriptor,
                metrics;
                term,
            )
            expected_left_right = _pqs_pqs_source_box_explicit_reference(
                left_descriptor,
                right_descriptor,
                metrics;
                term,
            )
            expected_right_right = _pqs_product_box_column_selection_reference(
                right_descriptor,
                metrics;
                term,
            )
            expected_left_product = _pqs_product_source_box_explicit_reference(
                left_descriptor,
                product_unit,
                metrics;
                term,
            )
            expected_right_product = _pqs_product_source_box_explicit_reference(
                right_descriptor,
                product_unit,
                metrics;
                term,
            )
            expected_product = _product_source_box_explicit_reference(
                product_unit,
                metrics;
                term,
            )
            @test size(block) == (retained_dimension, retained_dimension)
            @test all(isfinite, block)
            @test block[left_range, left_range] ≈ expected_left_left atol = 1.0e-10 rtol = 1.0e-10
            @test block[left_range, right_range] ≈ expected_left_right atol = 1.0e-10 rtol = 1.0e-10
            @test block[right_range, left_range] ≈ transpose(expected_left_right) atol = 1.0e-10 rtol = 1.0e-10
            @test block[right_range, right_range] ≈ expected_right_right atol = 1.0e-10 rtol = 1.0e-10
            @test block[left_range, product_range] ≈ expected_left_product atol = 1.0e-10 rtol = 1.0e-10
            @test block[product_range, left_range] ≈ transpose(expected_left_product) atol = 1.0e-10 rtol = 1.0e-10
            @test block[right_range, product_range] ≈ expected_right_product atol = 1.0e-10 rtol = 1.0e-10
            @test block[product_range, right_range] ≈ transpose(expected_right_product) atol = 1.0e-10 rtol = 1.0e-10
            @test block[product_range, product_range] ≈ expected_product atol = 1.0e-10 rtol = 1.0e-10
            @test components.pqs_left_pqs_left ≈ expected_left_left atol = 1.0e-10 rtol = 1.0e-10
            @test components.pqs_left_pqs_right ≈ expected_left_right atol = 1.0e-10 rtol = 1.0e-10
            @test components.pqs_right_pqs_left ≈ transpose(expected_left_right) atol = 1.0e-10 rtol = 1.0e-10
            @test components.pqs_left_product ≈ expected_left_product atol = 1.0e-10 rtol = 1.0e-10
            @test components.product_pqs_left ≈ transpose(expected_left_product) atol = 1.0e-10 rtol = 1.0e-10
            @test components.pqs_right_pqs_right ≈ expected_right_right atol = 1.0e-10 rtol = 1.0e-10
            @test components.pqs_right_product ≈ expected_right_product atol = 1.0e-10 rtol = 1.0e-10
            @test components.product_pqs_right ≈ transpose(expected_right_product) atol = 1.0e-10 rtol = 1.0e-10
            @test components.product_product ≈ expected_product atol = 1.0e-10 rtol = 1.0e-10
            max_component_error = max(
                max_component_error,
                LinearAlgebra.norm(
                    block[left_range, right_range] - expected_left_right,
                    Inf,
                ),
                LinearAlgebra.norm(
                    block[left_range, product_range] - expected_left_product,
                    Inf,
                ),
                LinearAlgebra.norm(
                    block[right_range, product_range] - expected_right_product,
                    Inf,
                ),
                LinearAlgebra.norm(
                    block[product_range, product_range] - expected_product,
                    Inf,
                ),
            )
        end
        @test max_component_error < 1.0e-10
        @test_throws ArgumentError metrics_module._pqs_pqs_product_source_box_shadow_blocks(
            left_pqs_plan,
            right_pqs_plan,
            product_unit,
            metrics;
            terms = (:weights,),
        )
    end

    function _check_pqs_pqs_product_route_shaped_safe_term_consumer(
        metrics_module,
        left_descriptor,
        right_descriptor,
        route_units,
        metrics,
    )
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
        units = route_units.units
        shadow = metrics_module._pqs_pqs_product_source_box_shadow_blocks(
            units.pqs_left,
            units.pqs_right,
            units.product,
            metrics;
            terms,
        )
        consumer = metrics_module._pqs_pqs_product_route_shaped_safe_term_consumer(
            route_units,
            metrics;
            terms,
        )

        @test consumer.path == :pqs_pqs_product_route_shaped_safe_term_consumer
        @test consumer.route_kind ==
              :pqs_pqs_product_source_box_safe_term_route
        @test consumer.route_name == :q5_L5_slab5_test_route
        @test consumer.route_units === route_units
        @test route_units.object_kind ==
              :pqs_pqs_product_safe_term_route_descriptor
        @test route_units.expected_ranges == shadow.ranges
        @test route_units.retained_dimension == 221
        @test route_units.retained_unit_count == 3
        @test route_units.expected_pair_count == 6
        @test route_units.supported_terms == terms
        @test length(route_units.unit_summaries) == 3
        @test route_units.unit_summaries[1].unit_key == :pqs_left
        @test route_units.unit_summaries[1].source_family ==
              :mode_selected_raw_product_box
        @test route_units.unit_summaries[1].source_dimensions == (5, 5, 5)
        @test route_units.unit_summaries[1].retained_count == 98
        @test route_units.unit_summaries[1].retained_rule_kind ==
              :boundary_comx_product_mode_selection
        @test route_units.unit_summaries[2].unit_key == :pqs_right
        @test route_units.unit_summaries[2].source_dimensions == (5, 5, 5)
        @test route_units.unit_summaries[2].retained_count == 98
        @test route_units.unit_summaries[3].unit_key == :product
        @test route_units.unit_summaries[3].source_family == :product_doside
        @test route_units.unit_summaries[3].source_dimensions == (5, 5, 1)
        @test route_units.unit_summaries[3].retained_count == 25
        @test consumer.terms == terms
        @test consumer.ranges == shadow.ranges
        @test consumer.retained_dimension == 221
        @test consumer.retained_dimension == shadow.retained_dimension
        @test consumer.pair_count == 6
        @test consumer.term_count == length(terms)
        @test length(consumer.retained_units) == 3
        @test length(consumer.all_pairs_inventory.pair_entries) == 6
        @test consumer.all_pairs_inventory.object_kind ==
              :pqs_pqs_product_source_box_all_pairs_inventory
        @test consumer.safe_term_matrices === consumer.blocks
        @test consumer.complete_retained_space_matrices === consumer.blocks
        @test consumer.shadow.path == :pqs_pqs_product_source_box_shadow_blocks
        @test consumer.metadata.parent_dims == (5, 5, 7)
        @test consumer.metadata.bond_axis == :z
        @test consumer.metadata.product_slab_fixed_index == 4
        @test consumer.metadata.pqs_source_mode_dims == (5, 5, 5)
        @test consumer.provenance.source ==
              :route_shaped_safe_term_consumer_test_fixture

        @test consumer.performance.elapsed_seconds >= 0.0
        @test consumer.performance.allocated_bytes >= 0
        @test consumer.performance.gc_time_seconds >= 0.0
        @test consumer.performance.retained_dimension == 221
        @test consumer.performance.pair_count == 6
        @test consumer.performance.term_count == length(terms)
        @test !consumer.performance.dense_raw_source_box_pair_matrix_materialized
        @test consumer.performance.dense_raw_source_box_pair_matrix_materialized_for_validation
        @test consumer.performance.dense_raw_pair_storage_avoided

        @test consumer.diagnostics.source ==
              :pqs_pqs_product_route_shaped_safe_term_consumer
        @test consumer.diagnostics.route_shaped_consumer
        @test consumer.diagnostics.route_kind ==
              :pqs_pqs_product_source_box_safe_term_route
        @test consumer.diagnostics.route_descriptor_object_kind ==
              :pqs_pqs_product_safe_term_route_descriptor
        @test consumer.diagnostics.descriptor_expected_ranges_checked
        @test consumer.diagnostics.descriptor_retained_dimension_checked
        @test consumer.diagnostics.descriptor_pair_count_checked
        @test consumer.diagnostics.descriptor_supported_terms_checked
        @test consumer.diagnostics.route_roles ==
              (:pqs_left, :pqs_right, :product)
        @test consumer.diagnostics.private_shadow_only
        @test !consumer.diagnostics.packet_adoption
        @test !consumer.diagnostics.fixed_block_routing
        @test !consumer.diagnostics.qwhamiltonian_consumes
        @test !consumer.diagnostics.public_default_consumes
        @test !consumer.diagnostics.cr2_science_status_changed
        @test !consumer.diagnostics.shell_projection_used
        @test !consumer.diagnostics.lowdin_cleanup_used
        @test !consumer.diagnostics.support_local_pqs_oracle_used
        @test consumer.diagnostics.retained_weight_semantics ==
              :not_positive_quadrature_weights
        @test !consumer.diagnostics.ida_weight_division_allowed
        @test !consumer.diagnostics.dense_raw_source_box_pair_matrix_materialized
        @test consumer.diagnostics.dense_raw_source_box_pair_matrix_materialized_for_validation
        @test consumer.diagnostics.dense_raw_pair_storage_avoided
        @test consumer.diagnostics.explicit_source_box_oracle_tested
        @test !consumer.diagnostics.pqs_cross_box_external_raw_product_oracle_required
        @test consumer.diagnostics.pqs_cross_box_internal_raw_product_oracle_compared
        @test consumer.diagnostics.every_pair_uses_source_box_algorithmic_policy
        @test consumer.diagnostics.source_box_algorithmic_pair_count == 6
        @test consumer.diagnostics.pair_policies ==
              (
                  :source_box_algorithm_available,
                  :source_box_algorithm_available,
                  :source_box_algorithm_available,
                  :source_box_algorithm_available,
                  :source_box_algorithm_available,
                  :source_box_algorithm_available,
              )
        @test consumer.diagnostics.performance_recorded
        @test consumer.diagnostics.retained_dimension == 221
        @test consumer.diagnostics.retained_unit_count == 3
        @test consumer.diagnostics.pair_count == 6
        @test consumer.diagnostics.term_count == length(terms)
        @test consumer.diagnostics.complete_retained_space_matrices_built
        @test consumer.diagnostics.source_box_shadow_helper ==
              :_pqs_pqs_product_source_box_shadow_blocks

        max_full_error = 0.0
        max_cross_pqs_error = 0.0
        for term in terms
            @test haskey(consumer.blocks, term)
            @test haskey(consumer.component_blocks, term)
            block = consumer.blocks[term]
            @test size(block) == (221, 221)
            @test all(isfinite, block)
            @test block ≈ shadow.blocks[term] atol = 0.0 rtol = 0.0
            expected_left_right = _pqs_pqs_source_box_explicit_reference(
                left_descriptor,
                right_descriptor,
                metrics;
                term,
            )
            @test consumer.component_blocks[term].pqs_left_pqs_right ≈
                  expected_left_right atol = 1.0e-10 rtol = 1.0e-10
            @test consumer.component_blocks[term].pqs_right_pqs_left ≈
                  transpose(expected_left_right) atol = 1.0e-10 rtol = 1.0e-10
            @test block[consumer.ranges.pqs_left, consumer.ranges.pqs_right] ≈
                  expected_left_right atol = 1.0e-10 rtol = 1.0e-10
            @test block[consumer.ranges.pqs_right, consumer.ranges.pqs_left] ≈
                  transpose(expected_left_right) atol = 1.0e-10 rtol = 1.0e-10
            max_full_error = max(
                max_full_error,
                LinearAlgebra.norm(block - shadow.blocks[term], Inf),
            )
            max_cross_pqs_error = max(
                max_cross_pqs_error,
                LinearAlgebra.norm(
                    block[consumer.ranges.pqs_left, consumer.ranges.pqs_right] -
                    expected_left_right,
                    Inf,
                ),
            )
        end
        @test max_full_error == 0.0
        @test max_cross_pqs_error < 1.0e-10
        @test_throws ArgumentError metrics_module._pqs_pqs_product_route_shaped_safe_term_consumer(
            route_units,
            metrics;
            terms = (:weights,),
        )
    end

    function _density_density_route_factor_fixture(parent_dims::NTuple{3,Int})
        term_coefficients = [0.47, 0.19]
        function _symmetric_pair_term_tensor(source_count::Int, axis_scale::Float64)
            terms = Array{Float64,3}(undef, 2, source_count, source_count)
            for term in 1:2, col in 1:source_count, row in 1:source_count
                terms[term, row, col] =
                    axis_scale / (1.0 + abs(row - col)) +
                    0.025 * term +
                    0.004 * (row + col)
            end
            return terms
        end
        function _raw_terms_from_density_normalized(normalized_terms, weights)
            raw_terms = similar(normalized_terms)
            weight_outer = weights * transpose(weights)
            for term in axes(normalized_terms, 1)
                raw_terms[term, :, :] .= normalized_terms[term, :, :] .* weight_outer
            end
            return raw_terms
        end
        density_pair_terms = (
            x = _symmetric_pair_term_tensor(parent_dims[1], 0.41),
            y = _symmetric_pair_term_tensor(parent_dims[2], 0.37),
            z = _symmetric_pair_term_tensor(parent_dims[3], 0.31),
        )
        density_source_weights = (
            x = [1.0 + 0.03 * index for index in 1:parent_dims[1]],
            y = [0.9 + 0.04 * index for index in 1:parent_dims[2]],
            z = [1.2 + 0.02 * index for index in 1:parent_dims[3]],
        )
        raw_pair_terms = (
            x = _raw_terms_from_density_normalized(
                density_pair_terms.x,
                density_source_weights.x,
            ),
            y = _raw_terms_from_density_normalized(
                density_pair_terms.y,
                density_source_weights.y,
            ),
            z = _raw_terms_from_density_normalized(
                density_pair_terms.z,
                density_source_weights.z,
            ),
        )
        return (
            term_coefficients = term_coefficients,
            axis_pair_factor_terms = density_pair_terms,
            axis_weights = density_source_weights,
            raw_axis_pair_factor_terms = raw_pair_terms,
        )
    end

    function _check_pqs_pqs_product_route_shaped_density_density_consumer(
        metrics_module,
        route_units,
    )
        density_fixture =
            _density_density_route_factor_fixture(route_units.metadata.parent_dims)
        term_coefficients = density_fixture.term_coefficients
        density_pair_terms = density_fixture.axis_pair_factor_terms
        density_source_weights = density_fixture.axis_weights
        consumer =
            metrics_module._pqs_pqs_product_route_shaped_density_density_consumer(
                route_units;
                term_coefficients,
                axis_pair_factor_terms = density_pair_terms,
                axis_weights = density_source_weights,
            )
        units = route_units.units
        ranges = consumer.ranges
        left_left =
            metrics_module._pqs_pqs_source_box_density_density_interaction_block(
                units.pqs_left,
                units.pqs_left;
                term_coefficients,
                axis_pair_factor_terms = density_pair_terms,
                axis_weights = density_source_weights,
            )
        left_right =
            metrics_module._pqs_pqs_source_box_density_density_interaction_block(
                units.pqs_left,
                units.pqs_right;
                term_coefficients,
                axis_pair_factor_terms = density_pair_terms,
                axis_weights = density_source_weights,
            )
        right_right =
            metrics_module._pqs_pqs_source_box_density_density_interaction_block(
                units.pqs_right,
                units.pqs_right;
                term_coefficients,
                axis_pair_factor_terms = density_pair_terms,
                axis_weights = density_source_weights,
            )
        left_product =
            metrics_module._pqs_product_source_box_density_density_interaction_block(
                units.pqs_left,
                units.product;
                term_coefficients,
                axis_pair_factor_terms = density_pair_terms,
                axis_weights = density_source_weights,
            )
        right_product =
            metrics_module._pqs_product_source_box_density_density_interaction_block(
                units.pqs_right,
                units.product;
                term_coefficients,
                axis_pair_factor_terms = density_pair_terms,
                axis_weights = density_source_weights,
            )
        product_product =
            metrics_module._product_doside_source_box_density_density_interaction_block(
                units.product,
                units.product;
                term_coefficients,
                axis_pair_factor_terms = density_pair_terms,
                axis_weights = density_source_weights,
            )
        explicit = zeros(Float64, consumer.retained_dimension, consumer.retained_dimension)
        explicit[ranges.pqs_left, ranges.pqs_left] .= left_left.block
        explicit[ranges.pqs_left, ranges.pqs_right] .= left_right.block
        explicit[ranges.pqs_right, ranges.pqs_left] .= transpose(left_right.block)
        explicit[ranges.pqs_left, ranges.product] .= left_product.block
        explicit[ranges.product, ranges.pqs_left] .= transpose(left_product.block)
        explicit[ranges.pqs_right, ranges.pqs_right] .= right_right.block
        explicit[ranges.pqs_right, ranges.product] .= right_product.block
        explicit[ranges.product, ranges.pqs_right] .= transpose(right_product.block)
        explicit[ranges.product, ranges.product] .= product_product.block

        @test consumer.path ==
              :pqs_pqs_product_route_shaped_density_density_consumer
        @test consumer.route_kind ==
              :pqs_pqs_product_source_box_safe_term_route
        @test consumer.route_units === route_units
        @test consumer.retained_dimension == route_units.retained_dimension
        @test consumer.retained_dimension == 221
        @test consumer.pair_count == 6
        @test consumer.pair_family_counts ==
              (pqs_pqs = 3, pqs_product = 2, product_product = 1)
        @test consumer.term_count == length(term_coefficients)
        @test consumer.pair_factor_normalization == :density_normalized
        @test consumer.output_finite
        @test all(isfinite, consumer.block)
        @test consumer.symmetry_error <= 1.0e-10
        @test consumer.block ≈ explicit atol = 1.0e-12 rtol = 1.0e-12
        @test consumer.component_blocks.pqs_left_pqs_left ≈
              left_left.block atol = 0.0 rtol = 0.0
        @test consumer.component_blocks.pqs_left_pqs_right ≈
              left_right.block atol = 0.0 rtol = 0.0
        @test consumer.component_blocks.pqs_right_pqs_left ≈
              transpose(left_right.block) atol = 0.0 rtol = 0.0
        @test consumer.component_blocks.pqs_left_product ≈
              left_product.block atol = 0.0 rtol = 0.0
        @test consumer.component_blocks.product_pqs_left ≈
              transpose(left_product.block) atol = 0.0 rtol = 0.0
        @test consumer.component_blocks.product_product ≈
              product_product.block atol = 0.0 rtol = 0.0
        @test length(consumer.retained_units) == 3
        @test map(unit -> unit.unit_key, consumer.retained_units) ==
              (:pqs_left, :pqs_right, :product)
        @test map(unit -> unit.retained_range, consumer.retained_units) ==
              (ranges.pqs_left, ranges.pqs_right, ranges.product)
        @test consumer.all_pairs_inventory.object_kind ==
              :pqs_pqs_product_density_density_all_pairs_inventory
        @test length(consumer.all_pairs_inventory.pair_entries) == 6
        @test consumer.all_pairs_inventory.diagnostics.pair_family_counts ==
              consumer.pair_family_counts
        @test consumer.all_pairs_inventory.diagnostics.block_helper_by_family ==
              (
                  pqs_pqs =
                      :_pqs_pqs_source_box_density_density_interaction_block,
                  pqs_product =
                      :_pqs_product_source_box_density_density_interaction_block,
                  product_product =
                      :_product_doside_source_box_density_density_interaction_block,
              )

        @test consumer.diagnostics.source ==
              :pqs_pqs_product_route_shaped_density_density_consumer
        @test consumer.diagnostics.route_shaped_consumer
        @test consumer.diagnostics.route_shape ==
              (:pqs_left, :pqs_right, :product)
        @test consumer.diagnostics.retained_dimension == 221
        @test consumer.diagnostics.retained_unit_count == 3
        @test consumer.diagnostics.pair_count == 6
        @test consumer.diagnostics.pair_family_counts ==
              (pqs_pqs = 3, pqs_product = 2, product_product = 1)
        @test consumer.diagnostics.pair_factor_normalization ==
              :density_normalized
        @test consumer.diagnostics.density_normalized_pair_factors
        @test !consumer.diagnostics.raw_weighted_pair_factors
        @test !consumer.diagnostics.source_weight_division_applied_by_helper
        @test consumer.diagnostics.source_weight_division_owner ==
              :caller_supplied_density_normalized_pair_factors
        @test consumer.diagnostics.output_representation ==
              :two_index_density_density
        @test !consumer.diagnostics.four_index_galerkin_tensor
        @test consumer.diagnostics.interaction_operator ==
              :electron_electron_density_density
        @test consumer.diagnostics.source_box_first
        @test consumer.diagnostics.source_box_algorithmic_path_true_for_every_pair
        @test consumer.diagnostics.every_pair_uses_source_box_algorithmic_policy
        @test consumer.diagnostics.source_box_algorithmic_pair_count == 6
        @test consumer.diagnostics.helper_used_for_pair_families ==
              consumer.all_pairs_inventory.diagnostics.block_helper_by_family
        @test consumer.diagnostics.product_pqs_blocks_transpose_only
        @test consumer.diagnostics.lower_triangular_cross_blocks_transpose_only
        @test consumer.diagnostics.pair_factor_terms_symmetric
        @test consumer.diagnostics.symmetric_same_route_input
        @test consumer.diagnostics.output_finite
        @test consumer.diagnostics.symmetry_error <= 1.0e-10
        @test consumer.diagnostics.descriptor_expected_ranges_checked
        @test consumer.diagnostics.descriptor_retained_dimension_checked
        @test consumer.diagnostics.descriptor_pair_count_checked
        @test consumer.diagnostics.private_shadow_only
        @test consumer.diagnostics.input_pair_factor_data ==
              :caller_supplied_explicit_data
        @test !consumer.diagnostics.input_pair_factor_data_pgdg_checked
        @test !consumer.diagnostics.real_mwg_ida_pair_factor_provenance_adapted
        @test consumer.diagnostics.source_weights_are_raw_source_weights
        @test !consumer.diagnostics.retained_pqs_weights_used
        @test !consumer.diagnostics.retained_pqs_weights_positive_checked
        @test !consumer.diagnostics.retained_weight_division_allowed
        @test !consumer.diagnostics.retained_pqs_weight_division_allowed
        @test !consumer.diagnostics.ida_weight_division_allowed
        @test consumer.diagnostics.retained_weight_semantics ==
              :not_positive_quadrature_weights
        @test !consumer.diagnostics.shell_projection_used
        @test !consumer.diagnostics.lowdin_cleanup_used
        @test !consumer.diagnostics.support_local_oracle_used
        @test !consumer.diagnostics.support_local_pqs_oracle_used
        @test !consumer.diagnostics.support_local_shell_row_algorithm
        @test !consumer.diagnostics.support_coefficient_matrix_used
        @test !consumer.diagnostics.shell_row_algorithm
        @test !consumer.diagnostics.packet_adoption
        @test !consumer.diagnostics.fixed_block_routing
        @test !consumer.diagnostics.qwhamiltonian_consumes
        @test !consumer.diagnostics.public_default_consumes
        @test !consumer.diagnostics.ecp_terms_implemented
        @test !consumer.diagnostics.cr2_science_status_changed
        @test !consumer.diagnostics.ida_mwg_semantics_changed
        @test !consumer.diagnostics.mwg_ida_semantics_changed
        @test !consumer.diagnostics.mwg_interaction_implemented
        @test consumer.diagnostics.complete_retained_space_matrix_built
        @test !consumer.diagnostics.dense_raw_source_box_pair_matrix_materialized
        @test consumer.diagnostics.dense_raw_pair_storage_avoided

        raw_pair_terms = density_fixture.raw_axis_pair_factor_terms
        raw_consumer =
            metrics_module._pqs_pqs_product_route_shaped_density_density_consumer(
                route_units;
                term_coefficients,
                raw_axis_pair_factor_terms = raw_pair_terms,
                axis_weights = density_source_weights,
                pair_factor_normalization = :raw_weighted,
            )
        @test raw_consumer.block ≈ consumer.block atol = 1.0e-12 rtol = 1.0e-12
        @test raw_consumer.diagnostics.pair_factor_normalization ==
              :raw_weighted
        @test raw_consumer.diagnostics.raw_weighted_pair_factors
        @test !raw_consumer.diagnostics.density_normalized_pair_factors
        @test raw_consumer.diagnostics.density_normalized_pair_factors_generated
        @test raw_consumer.diagnostics.source_weight_division_owner ==
              :source_box_raw_weights
        @test raw_consumer.diagnostics.source_weight_division_applied_by_helper
        @test raw_consumer.diagnostics.source_weight_division_shape ==
              :axis_pair_weight_outer
        @test raw_consumer.all_pairs_inventory.diagnostics.block_helper_by_family ==
              (
                  pqs_pqs =
                      :_pqs_pqs_source_box_raw_weighted_density_density_interaction_block,
                  pqs_product =
                      :_pqs_product_source_box_raw_weighted_density_density_interaction_block,
                  product_product =
                      :_product_doside_source_box_raw_weighted_density_density_interaction_block,
              )
        @test !raw_consumer.diagnostics.retained_pqs_weights_used
        @test !raw_consumer.diagnostics.retained_weight_division_allowed
        @test !raw_consumer.diagnostics.retained_pqs_weight_division_allowed
        @test !raw_consumer.diagnostics.ida_weight_division_allowed
        @test !raw_consumer.diagnostics.packet_adoption
        @test !raw_consumer.diagnostics.fixed_block_routing
        @test !raw_consumer.diagnostics.qwhamiltonian_consumes
        @test !raw_consumer.diagnostics.public_default_consumes
        @test !raw_consumer.diagnostics.ida_mwg_semantics_changed
        @test !raw_consumer.diagnostics.mwg_ida_semantics_changed
        @test !raw_consumer.diagnostics.real_mwg_ida_pair_factor_provenance_adapted

        nonsymmetric_terms = (
            x = copy(density_pair_terms.x),
            y = density_pair_terms.y,
            z = density_pair_terms.z,
        )
        nonsymmetric_terms.x[1, 1, 2] += 0.125
        @test_throws ArgumentError metrics_module._pqs_pqs_product_route_shaped_density_density_consumer(
            route_units;
            term_coefficients,
            axis_pair_factor_terms = nonsymmetric_terms,
            axis_weights = density_source_weights,
        )
        @test_throws ArgumentError metrics_module._pqs_pqs_product_route_shaped_density_density_consumer(
            route_units;
            term_coefficients,
            axis_pair_factor_terms = density_pair_terms,
            axis_weights = (
                x = density_source_weights.x,
                y = density_source_weights.y,
                z = [1.22, 1.24, 1.26, 0.0, 1.3, 1.32, 1.34],
            ),
        )
    end

    function _check_pqs_pqs_product_raw_box_density_density_route_producer(
        metrics_module,
        bundles,
        route_units,
        metrics,
        left_source_box::NTuple{3,UnitRange{Int}},
        right_source_box::NTuple{3,UnitRange{Int}},
        product_source_box::NTuple{3,UnitRange{Int}};
        source_mode_dims::NTuple{3,Int},
        parent_dims::NTuple{3,Int},
        bond_axis::Symbol,
        ida_term_coefficients::AbstractVector{<:Real},
        ida_dense_parent_matrix::AbstractMatrix{<:Real},
    )
        density_fixture = _density_density_route_factor_fixture(parent_dims)
        direct_consumer =
            metrics_module._pqs_pqs_product_route_shaped_density_density_consumer(
                route_units;
                term_coefficients = density_fixture.term_coefficients,
                axis_pair_factor_terms = density_fixture.axis_pair_factor_terms,
                axis_weights = density_fixture.axis_weights,
            )
        producer =
            metrics_module._pqs_pqs_product_raw_box_density_density_route_producer(
                bundles,
                left_source_box,
                right_source_box,
                product_source_box,
                metrics;
                source_mode_dims,
                route_name = route_units.route_name,
                parent_dims,
                bond_axis,
                metadata = (
                    pqs_left_box = left_source_box,
                    pqs_right_box = right_source_box,
                    product_slab_fixed_index =
                        only(product_source_box[findfirst(
                            axis -> length(product_source_box[axis]) == 1,
                            1:3,
                        )]),
                    pqs_source_mode_dims = source_mode_dims,
                    density_density_route_producer_test = true,
                ),
                provenance = (
                    source =
                        :route_shaped_density_density_route_producer_test_fixture,
                ),
                term_coefficients = density_fixture.term_coefficients,
                axis_pair_factor_terms = density_fixture.axis_pair_factor_terms,
                axis_weights = density_fixture.axis_weights,
            )

        @test producer.object_kind ==
              :pqs_pqs_product_raw_box_density_density_route_producer
        @test producer.status == :private_density_density_reference_only
        @test producer.raw_box_route_producer.object_kind ==
              :pqs_pqs_product_raw_box_route_producer
        @test producer.descriptor === producer.raw_box_route_producer.descriptor
        @test producer.route_descriptor === producer.descriptor
        @test producer.consumer === producer.density_density_consumer
        @test producer.descriptor.object_kind ==
              :pqs_pqs_product_safe_term_route_descriptor
        @test producer.descriptor.expected_ranges == route_units.expected_ranges
        @test producer.descriptor.retained_dimension == route_units.retained_dimension
        @test producer.descriptor.retained_dimension == 221
        @test producer.descriptor.expected_pair_count == 6
        @test producer.retained_dimension == 221
        @test producer.pair_count == 6
        @test producer.pair_family_counts ==
              (pqs_pqs = 3, pqs_product = 2, product_product = 1)
        @test producer.term_count == length(density_fixture.term_coefficients)
        @test producer.pair_factor_normalization == :density_normalized
        @test producer.output_finite
        @test all(isfinite, producer.block)
        @test producer.symmetry_error <= 1.0e-10
        @test producer.block ≈ direct_consumer.block atol = 1.0e-12 rtol = 1.0e-12
        @test producer.density_density_matrix ===
              producer.density_density_consumer.density_density_matrix
        @test producer.complete_retained_space_matrix ===
              producer.density_density_consumer.complete_retained_space_matrix
        @test producer.all_pairs_inventory.object_kind ==
              :pqs_pqs_product_density_density_all_pairs_inventory
        @test producer.all_pairs_inventory.diagnostics.pair_family_counts ==
              producer.pair_family_counts

        @test producer.diagnostics.source ==
              :pqs_pqs_product_raw_box_density_density_route_producer
        @test producer.diagnostics.private_density_density_reference_only
        @test producer.diagnostics.private_shadow_only
        @test producer.diagnostics.raw_box_route_producer_called
        @test producer.diagnostics.route_descriptor_emitted
        @test producer.diagnostics.route_descriptor_source ==
              :pqs_pqs_product_raw_box_route_producer
        @test producer.diagnostics.route_descriptor_object_kind ==
              :pqs_pqs_product_safe_term_route_descriptor
        @test producer.diagnostics.route_descriptor_provenance.source ==
              :route_shaped_density_density_route_producer_test_fixture
        @test producer.diagnostics.route_descriptor_built_from_explicit_fixture_facts
        @test producer.diagnostics.density_density_consumer_called
        @test producer.diagnostics.density_density_consumer_path ==
              :pqs_pqs_product_route_shaped_density_density_consumer
        @test producer.diagnostics.returns_descriptor_and_density_density_consumer_result
        @test producer.diagnostics.retained_dimension == 221
        @test producer.diagnostics.pair_count == 6
        @test producer.diagnostics.pair_family_counts ==
              (pqs_pqs = 3, pqs_product = 2, product_product = 1)
        @test producer.diagnostics.pair_factor_normalization ==
              :density_normalized
        @test producer.diagnostics.input_pair_factor_data ==
              :caller_supplied_explicit_data
        @test producer.diagnostics.synthetic_or_caller_supplied_pair_factors
        @test !producer.diagnostics.real_mwg_ida_pair_factor_provenance_adapted
        @test producer.diagnostics.source_box_first
        @test producer.diagnostics.source_box_algorithmic_path_true_for_every_pair
        @test producer.diagnostics.every_pair_uses_source_box_algorithmic_policy
        @test producer.diagnostics.source_box_algorithmic_pair_count == 6
        @test !producer.diagnostics.shell_projection_used
        @test !producer.diagnostics.lowdin_cleanup_used
        @test !producer.diagnostics.support_local_oracle_used
        @test !producer.diagnostics.support_local_pqs_oracle_used
        @test !producer.diagnostics.support_local_shell_row_algorithm
        @test !producer.diagnostics.support_coefficient_matrix_used
        @test !producer.diagnostics.shell_row_algorithm
        @test !producer.diagnostics.retained_pqs_weights_used
        @test !producer.diagnostics.retained_pqs_weights_positive_checked
        @test !producer.diagnostics.retained_weight_division_allowed
        @test !producer.diagnostics.retained_pqs_weight_division_allowed
        @test !producer.diagnostics.ida_weight_division_allowed
        @test producer.diagnostics.retained_weight_semantics ==
              :not_positive_quadrature_weights
        @test !producer.diagnostics.packet_adoption
        @test !producer.diagnostics.fixed_block_routing
        @test !producer.diagnostics.qwhamiltonian_consumes
        @test !producer.diagnostics.public_default_consumes
        @test !producer.diagnostics.ecp_terms_implemented
        @test !producer.diagnostics.cr2_science_status_changed
        @test !producer.diagnostics.ida_mwg_semantics_changed
        @test !producer.diagnostics.mwg_ida_semantics_changed
        @test !producer.diagnostics.mwg_interaction_implemented

        raw_producer =
            metrics_module._pqs_pqs_product_raw_box_density_density_route_producer(
                bundles,
                left_source_box,
                right_source_box,
                product_source_box,
                metrics;
                source_mode_dims,
                route_name = route_units.route_name,
                parent_dims,
                bond_axis,
                metadata = (
                    density_density_raw_weighted_route_producer_test = true,
                ),
                provenance = (
                    source =
                        :route_shaped_density_density_route_producer_test_fixture,
                ),
                term_coefficients = density_fixture.term_coefficients,
                raw_axis_pair_factor_terms =
                    density_fixture.raw_axis_pair_factor_terms,
                axis_weights = density_fixture.axis_weights,
                pair_factor_normalization = :raw_weighted,
            )
        @test raw_producer.block ≈ producer.block atol = 1.0e-12 rtol = 1.0e-12
        @test raw_producer.diagnostics.pair_factor_normalization ==
              :raw_weighted
        @test raw_producer.diagnostics.raw_weighted_pair_factors
        @test !raw_producer.diagnostics.density_normalized_pair_factors
        @test raw_producer.diagnostics.density_normalized_pair_factors_generated
        @test raw_producer.diagnostics.source_weight_division_owner ==
              :source_box_raw_weights
        @test raw_producer.diagnostics.source_weight_division_applied_by_helper
        @test raw_producer.diagnostics.source_weight_division_shape ==
              :axis_pair_weight_outer
        @test raw_producer.diagnostics.synthetic_or_caller_supplied_pair_factors
        @test !raw_producer.diagnostics.real_mwg_ida_pair_factor_provenance_adapted
        @test !raw_producer.diagnostics.retained_pqs_weights_used
        @test !raw_producer.diagnostics.retained_weight_division_allowed
        @test !raw_producer.diagnostics.retained_pqs_weight_division_allowed
        @test !raw_producer.diagnostics.ida_weight_division_allowed
        @test !raw_producer.diagnostics.packet_adoption
        @test !raw_producer.diagnostics.fixed_block_routing
        @test !raw_producer.diagnostics.qwhamiltonian_consumes
        @test !raw_producer.diagnostics.public_default_consumes

        ida_provenance = metrics_module._pqs_source_box_ida_factor_provenance(
            bundles;
            expected_term_count = length(ida_term_coefficients),
        )
        explicit_ida_producer =
            metrics_module._pqs_pqs_product_raw_box_density_density_route_producer(
                bundles,
                left_source_box,
                right_source_box,
                product_source_box,
                metrics;
                source_mode_dims,
                route_name = route_units.route_name,
                parent_dims,
                bond_axis,
                metadata = (
                    density_density_ida_explicit_route_producer_test = true,
                ),
                provenance = (
                    source =
                        :route_shaped_density_density_ida_provenance_test_fixture,
                ),
                term_coefficients = ida_term_coefficients,
                axis_pair_factor_terms = ida_provenance.axis_pair_factor_terms,
                axis_weights = ida_provenance.axis_weights,
            )
        ida_adapter =
            metrics_module._pqs_pqs_product_raw_box_density_density_route_producer_from_ida_provenance(
                bundles,
                left_source_box,
                right_source_box,
                product_source_box,
                metrics,
                ida_provenance;
                source_mode_dims,
                route_name = route_units.route_name,
                parent_dims,
                bond_axis,
                metadata = (
                    density_density_ida_provenance_adapter_test = true,
                ),
                provenance = (
                    source =
                        :route_shaped_density_density_ida_provenance_test_fixture,
                ),
                term_coefficients = ida_term_coefficients,
            )
        @test ida_adapter.object_kind ==
              :pqs_pqs_product_raw_box_density_density_route_producer_from_ida_provenance
        @test ida_adapter.status ==
              :private_ida_provenance_density_density_reference_only
        @test ida_adapter.explicit_input_route_producer.block ≈
              explicit_ida_producer.block atol = 0.0 rtol = 0.0
        @test ida_adapter.block ≈ explicit_ida_producer.block atol = 1.0e-12 rtol = 1.0e-12
        @test ida_adapter.pair_factor_normalization == :density_normalized
        @test ida_adapter.term_count == length(ida_term_coefficients)
        @test ida_adapter.diagnostics.input_pair_factor_data ==
              :ida_gausslet_source_box_provenance
        @test ida_adapter.diagnostics.input_pair_factor_data_pgdg_checked
        @test ida_adapter.diagnostics.interaction_path ==
              :ida_gausslet_source_box
        @test !ida_adapter.diagnostics.mwg_supplement_residual_path
        @test !ida_adapter.diagnostics.mwg_supplement_residual_provenance_adapted
        @test ida_adapter.diagnostics.real_ida_gausslet_source_box_provenance_adapted
        @test !ida_adapter.diagnostics.real_mwg_ida_pair_factor_provenance_adapted
        @test !ida_adapter.diagnostics.synthetic_or_caller_supplied_pair_factors
        @test !ida_adapter.diagnostics.retained_pqs_weights_used
        @test !ida_adapter.diagnostics.retained_weight_division_allowed
        @test !ida_adapter.diagnostics.retained_pqs_weight_division_allowed
        @test !ida_adapter.diagnostics.ida_weight_division_allowed
        @test !ida_adapter.diagnostics.shell_projection_used
        @test !ida_adapter.diagnostics.lowdin_cleanup_used
        @test !ida_adapter.diagnostics.support_local_oracle_used
        @test !ida_adapter.diagnostics.support_local_pqs_oracle_used
        @test !ida_adapter.diagnostics.support_local_shell_row_algorithm
        @test !ida_adapter.diagnostics.support_coefficient_matrix_used
        @test !ida_adapter.diagnostics.shell_row_algorithm
        @test !ida_adapter.diagnostics.packet_adoption
        @test !ida_adapter.diagnostics.fixed_block_routing
        @test !ida_adapter.diagnostics.qwhamiltonian_consumes
        @test !ida_adapter.diagnostics.public_default_consumes
        @test !ida_adapter.diagnostics.ecp_terms_implemented
        @test !ida_adapter.diagnostics.cr2_science_status_changed

        ida_parent_authority =
            metrics_module._pqs_pqs_product_dense_parent_ida_authority_comparison(
                ida_adapter,
                ida_dense_parent_matrix;
                parent_dims,
                dense_parent_matrix_source = :qwrg_diatomic_interaction_matrix,
                comparison_atol = 1.0e-10,
            )
        @test ida_parent_authority.object_kind ==
              :pqs_pqs_product_dense_parent_ida_authority_comparison
        @test ida_parent_authority.status == :private_validation_only
        @test ida_parent_authority.authority_kind ==
              :dense_parent_ida_projection
        @test size(ida_parent_authority.coefficient_matrix) ==
              (prod(parent_dims), ida_adapter.retained_dimension)
        @test size(ida_parent_authority.projected_block) ==
              size(ida_adapter.block)
        @test ida_parent_authority.output_finite
        @test ida_parent_authority.within_tolerance
        @test ida_parent_authority.max_error <= 1.0e-10
        @test ida_parent_authority.diagnostics.authority_kind ==
              :dense_parent_ida_projection
        @test ida_parent_authority.diagnostics.dense_parent_matrix_source ==
              :qwrg_diatomic_interaction_matrix
        @test ida_parent_authority.diagnostics.dense_parent_matrix_used_for_validation
        @test !ida_parent_authority.diagnostics.dense_parent_matrix_algorithmic
        @test !ida_parent_authority.diagnostics.source_box_algorithm_changed
        @test !ida_parent_authority.diagnostics.support_local_algorithm_used
        @test !ida_parent_authority.diagnostics.support_local_oracle_used
        @test !ida_parent_authority.diagnostics.support_coefficient_matrix_used
        @test !ida_parent_authority.diagnostics.shell_projection_used
        @test !ida_parent_authority.diagnostics.lowdin_cleanup_used
        @test !ida_parent_authority.diagnostics.retained_pqs_weights_used
        @test !ida_parent_authority.diagnostics.retained_pqs_weights_positive_checked
        @test !ida_parent_authority.diagnostics.retained_weight_division_allowed
        @test !ida_parent_authority.diagnostics.retained_pqs_weight_division_allowed
        @test !ida_parent_authority.diagnostics.ida_weight_division_allowed
        @test !ida_parent_authority.diagnostics.packet_adoption
        @test !ida_parent_authority.diagnostics.fixed_block_routing
        @test !ida_parent_authority.diagnostics.qwhamiltonian_consumes
        @test !ida_parent_authority.diagnostics.public_default_consumes
        @test !ida_parent_authority.diagnostics.ecp_terms_implemented
        @test !ida_parent_authority.diagnostics.cr2_science_status_changed
        @test !ida_parent_authority.diagnostics.mwg_supplement_residual_path
        @test !ida_parent_authority.diagnostics.mwg_supplement_residual_provenance_adapted
        @test ida_parent_authority.diagnostics.output_representation ==
              :two_index_density_density
        @test !ida_parent_authority.diagnostics.four_index_galerkin_tensor

        explicit_raw_ida_producer =
            metrics_module._pqs_pqs_product_raw_box_density_density_route_producer(
                bundles,
                left_source_box,
                right_source_box,
                product_source_box,
                metrics;
                source_mode_dims,
                route_name = route_units.route_name,
                parent_dims,
                bond_axis,
                metadata = (
                    density_density_ida_raw_explicit_route_producer_test = true,
                ),
                provenance = (
                    source =
                        :route_shaped_density_density_ida_provenance_test_fixture,
                ),
                term_coefficients = ida_term_coefficients,
                raw_axis_pair_factor_terms =
                    ida_provenance.raw_axis_pair_factor_terms,
                axis_weights = ida_provenance.axis_weights,
                pair_factor_normalization = :raw_weighted,
            )
        raw_ida_adapter =
            metrics_module._pqs_pqs_product_raw_box_density_density_route_producer_from_ida_provenance(
                bundles,
                left_source_box,
                right_source_box,
                product_source_box,
                metrics,
                ida_provenance;
                source_mode_dims,
                route_name = route_units.route_name,
                parent_dims,
                bond_axis,
                metadata = (
                    density_density_ida_raw_provenance_adapter_test = true,
                ),
                provenance = (
                    source =
                        :route_shaped_density_density_ida_provenance_test_fixture,
                ),
                term_coefficients = ida_term_coefficients,
                pair_factor_normalization = :raw_weighted,
            )
        @test raw_ida_adapter.block ≈ explicit_raw_ida_producer.block atol = 1.0e-12 rtol = 1.0e-12
        @test raw_ida_adapter.block ≈ ida_adapter.block atol = 1.0e-12 rtol = 1.0e-12
        @test raw_ida_adapter.pair_factor_normalization == :raw_weighted
        @test raw_ida_adapter.diagnostics.input_pair_factor_data ==
              :ida_gausslet_source_box_provenance
        @test raw_ida_adapter.diagnostics.interaction_path ==
              :ida_gausslet_source_box
        @test !raw_ida_adapter.diagnostics.mwg_supplement_residual_path
        @test raw_ida_adapter.diagnostics.source_weight_division_owner ==
              :pgdg_auxiliary_source_weights
        @test raw_ida_adapter.diagnostics.source_weight_division_shape ==
              :axis_pair_weight_outer
        @test raw_ida_adapter.diagnostics.real_ida_gausslet_source_box_provenance_adapted
        @test !raw_ida_adapter.diagnostics.real_mwg_ida_pair_factor_provenance_adapted
        @test !raw_ida_adapter.diagnostics.synthetic_or_caller_supplied_pair_factors
        @test !raw_ida_adapter.diagnostics.retained_pqs_weights_used
        @test !raw_ida_adapter.diagnostics.retained_weight_division_allowed
        @test !raw_ida_adapter.diagnostics.retained_pqs_weight_division_allowed
        @test !raw_ida_adapter.diagnostics.ida_weight_division_allowed
        @test !raw_ida_adapter.diagnostics.packet_adoption
        @test !raw_ida_adapter.diagnostics.fixed_block_routing
        @test !raw_ida_adapter.diagnostics.qwhamiltonian_consumes
        @test !raw_ida_adapter.diagnostics.public_default_consumes
    end

    function _check_pqs_pqs_product_source_box_component_route_smoke(
        metrics_module,
        bundles,
        metrics,
        left_source_box::NTuple{3,UnitRange{Int}},
        right_source_box::NTuple{3,UnitRange{Int}},
        product_source_box::NTuple{3,UnitRange{Int}};
        source_mode_dims::NTuple{3,Int},
        parent_dims::NTuple{3,Int},
        bond_axis::Symbol,
        term_coefficients::AbstractVector{<:Real},
        ida_dense_parent_matrix::AbstractMatrix{<:Real},
        nuclear_expansion::CoulombGaussianExpansion,
    )
        ida_provenance = metrics_module._pqs_source_box_ida_factor_provenance(
            bundles;
            expected_term_count = length(term_coefficients),
        )
        nuclear_axis_layers = (
            x = build_basis(UniformBasisSpec(
                :G10;
                xmin = -2.0,
                xmax = 2.0,
                spacing = 1.0,
            )),
            y = build_basis(UniformBasisSpec(
                :G10;
                xmin = -2.0,
                xmax = 2.0,
                spacing = 1.0,
            )),
            z = build_basis(UniformBasisSpec(
                :G10;
                xmin = -3.0,
                xmax = 3.0,
                spacing = 1.0,
            )),
        )
        nuclear_centers = (
            (0.0, 0.0, -0.45),
            (0.0, 0.0, 0.55),
        )
        nuclear_charges = (2.0, 2.0)
        center_labels = (:left_center, :right_center)
        component =
            metrics_module._pqs_pqs_product_source_box_component_route_smoke(
                bundles,
                left_source_box,
                right_source_box,
                product_source_box,
                metrics,
                ida_provenance,
                nuclear_axis_layers,
                nuclear_expansion;
                source_mode_dims,
                centers = nuclear_centers,
                nuclear_charges,
                center_labels,
                term_coefficients,
                route_name = :q5_L5_slab5_component_route_smoke_test,
                parent_dims,
                bond_axis,
                metadata = (component_route_smoke_test = true,),
                provenance = (source = :component_route_smoke_test_fixture,),
                dense_parent_ida_matrix = ida_dense_parent_matrix,
                dense_parent_matrix_source = :qwrg_diatomic_interaction_matrix,
            )
        raw_component =
            metrics_module._pqs_pqs_product_source_box_component_route_smoke(
                bundles,
                left_source_box,
                right_source_box,
                product_source_box,
                metrics,
                ida_provenance,
                nuclear_axis_layers,
                nuclear_expansion;
                source_mode_dims,
                centers = nuclear_centers,
                nuclear_charges,
                center_labels,
                term_coefficients,
                pair_factor_normalization = :raw_weighted,
                route_name =
                    :q5_L5_slab5_component_route_smoke_raw_weighted_test,
                parent_dims,
                bond_axis,
                metadata = (component_route_smoke_raw_weighted_test = true,),
                provenance = (source = :component_route_smoke_test_fixture,),
                dense_parent_ida_matrix = ida_dense_parent_matrix,
                dense_parent_matrix_source = :qwrg_diatomic_interaction_matrix,
            )
        direct_ida =
            metrics_module._pqs_pqs_product_raw_box_density_density_route_producer_from_ida_provenance(
                bundles,
                left_source_box,
                right_source_box,
                product_source_box,
                metrics,
                ida_provenance;
                source_mode_dims,
                term_coefficients,
                route_name = :q5_L5_slab5_component_route_smoke_test,
                parent_dims,
                bond_axis,
                metadata = (component_route_smoke_direct_ida_test = true,),
                provenance = (source = :component_route_smoke_test_fixture,),
            )
        raw_direct_ida =
            metrics_module._pqs_pqs_product_raw_box_density_density_route_producer_from_ida_provenance(
                bundles,
                left_source_box,
                right_source_box,
                product_source_box,
                metrics,
                ida_provenance;
                source_mode_dims,
                term_coefficients,
                pair_factor_normalization = :raw_weighted,
                route_name =
                    :q5_L5_slab5_component_route_smoke_raw_weighted_test,
                parent_dims,
                bond_axis,
                metadata =
                    (component_route_smoke_direct_raw_weighted_ida_test = true,),
                provenance = (source = :component_route_smoke_test_fixture,),
            )
        component_summary =
            metrics_module._pqs_pqs_product_component_route_smoke_summary(
                component,
            )
        raw_component_summary =
            metrics_module._pqs_pqs_product_component_route_smoke_summary(
                raw_component,
            )

        units = component.route_descriptor.units
        ranges = component.ranges
        pqs_left_left =
            metrics_module._pqs_pqs_source_box_nuclear_attraction_by_center(
                units.pqs_left,
                units.pqs_left,
                nuclear_axis_layers,
                nuclear_expansion;
                centers = nuclear_centers,
                nuclear_charges,
            )
        pqs_left_right =
            metrics_module._pqs_pqs_source_box_nuclear_attraction_by_center(
                units.pqs_left,
                units.pqs_right,
                nuclear_axis_layers,
                nuclear_expansion;
                centers = nuclear_centers,
                nuclear_charges,
            )
        pqs_right_right =
            metrics_module._pqs_pqs_source_box_nuclear_attraction_by_center(
                units.pqs_right,
                units.pqs_right,
                nuclear_axis_layers,
                nuclear_expansion;
                centers = nuclear_centers,
                nuclear_charges,
            )
        pqs_left_product =
            metrics_module._pqs_product_source_box_nuclear_attraction_by_center(
                units.pqs_left,
                units.product,
                nuclear_axis_layers,
                nuclear_expansion;
                centers = nuclear_centers,
                nuclear_charges,
            )
        pqs_right_product =
            metrics_module._pqs_product_source_box_nuclear_attraction_by_center(
                units.pqs_right,
                units.product,
                nuclear_axis_layers,
                nuclear_expansion;
                centers = nuclear_centers,
                nuclear_charges,
            )
        product_product =
            metrics_module._product_doside_source_box_nuclear_attraction_by_center(
                units.product,
                units.product,
                nuclear_axis_layers,
                nuclear_expansion;
                centers = nuclear_centers,
                nuclear_charges,
            )
        expected_total =
            zeros(Float64, component.retained_dimension, component.retained_dimension)
        for center_index in eachindex(nuclear_centers)
            explicit =
                zeros(Float64, component.retained_dimension, component.retained_dimension)
            explicit[ranges.pqs_left, ranges.pqs_left] .=
                pqs_left_left.blocks_by_center[center_index].block
            explicit[ranges.pqs_left, ranges.pqs_right] .=
                pqs_left_right.blocks_by_center[center_index].block
            explicit[ranges.pqs_right, ranges.pqs_left] .=
                transpose(pqs_left_right.blocks_by_center[center_index].block)
            explicit[ranges.pqs_left, ranges.product] .=
                pqs_left_product.blocks_by_center[center_index].block
            explicit[ranges.product, ranges.pqs_left] .=
                transpose(pqs_left_product.blocks_by_center[center_index].block)
            explicit[ranges.pqs_right, ranges.pqs_right] .=
                pqs_right_right.blocks_by_center[center_index].block
            explicit[ranges.pqs_right, ranges.product] .=
                pqs_right_product.blocks_by_center[center_index].block
            explicit[ranges.product, ranges.pqs_right] .=
                transpose(pqs_right_product.blocks_by_center[center_index].block)
            explicit[ranges.product, ranges.product] .=
                product_product.blocks_by_center[center_index].block
            expected_total .+= explicit

            @test component.per_center_nuclear_matrices[center_index].center_label ==
                  center_labels[center_index]
            @test component.per_center_nuclear_matrices[center_index].center ==
                  nuclear_centers[center_index]
            @test component.per_center_nuclear_matrices[center_index].nuclear_charge ==
                  nuclear_charges[center_index]
            @test component.per_center_nuclear_matrices[center_index].block ≈
                  explicit atol = 1.0e-12 rtol = 1.0e-12
            @test component.per_center_nuclear_matrices[center_index].symmetry_error <=
                  1.0e-10
            @test raw_component.per_center_nuclear_matrices[center_index].center_label ==
                  center_labels[center_index]
            @test raw_component.per_center_nuclear_matrices[center_index].center ==
                  nuclear_centers[center_index]
            @test raw_component.per_center_nuclear_matrices[center_index].nuclear_charge ==
                  nuclear_charges[center_index]
            @test raw_component.per_center_nuclear_matrices[center_index].block ≈
                  explicit atol = 1.0e-12 rtol = 1.0e-12
            @test raw_component.per_center_nuclear_matrices[center_index].symmetry_error <=
                  1.0e-10
        end

        @test component.object_kind ==
              :pqs_pqs_product_source_box_component_route_smoke
        @test raw_component.object_kind ==
              :pqs_pqs_product_source_box_component_route_smoke
        @test component_summary.object_kind ==
              :pqs_pqs_product_component_route_smoke_summary
        @test raw_component_summary.object_kind ==
              :pqs_pqs_product_component_route_smoke_summary
        @test component.status == :private_component_route_smoke
        @test raw_component.status == :private_component_route_smoke
        @test component_summary.status ==
              :private_component_route_smoke_summary
        @test raw_component_summary.status ==
              :private_component_route_smoke_summary
        @test component_summary.component_object_kind == component.object_kind
        @test raw_component_summary.component_object_kind ==
              raw_component.object_kind
        @test component_summary.component_status == component.status
        @test raw_component_summary.component_status == raw_component.status
        @test component.route_shape == (:pqs_left, :pqs_right, :product)
        @test raw_component.route_shape == (:pqs_left, :pqs_right, :product)
        @test component_summary.route_shape == component.route_shape
        @test raw_component_summary.route_shape == raw_component.route_shape
        @test component.retained_dimension == 221
        @test raw_component.retained_dimension == 221
        @test component_summary.retained_dimension == 221
        @test raw_component_summary.retained_dimension == 221
        @test component.route_descriptor.retained_dimension == 221
        @test raw_component.route_descriptor.retained_dimension == 221
        @test component.route_descriptor.expected_pair_count == 6
        @test raw_component.route_descriptor.expected_pair_count == 6
        @test map(unit -> unit.unit_key, component.route_descriptor.unit_summaries) ==
              (:pqs_left, :pqs_right, :product)
        @test map(unit -> unit.retained_range, component.route_descriptor.unit_summaries) ==
              (ranges.pqs_left, ranges.pqs_right, ranges.product)
        @test map(unit -> unit.unit_key, raw_component.route_descriptor.unit_summaries) ==
              (:pqs_left, :pqs_right, :product)
        @test map(unit -> unit.retained_range, raw_component.route_descriptor.unit_summaries) ==
              (ranges.pqs_left, ranges.pqs_right, ranges.product)
        @test length(component.per_center_nuclear_matrices) == length(nuclear_centers)
        @test length(raw_component.per_center_nuclear_matrices) ==
              length(nuclear_centers)
        @test component.center_labels == center_labels
        @test raw_component.center_labels == center_labels
        @test component_summary.center_labels == center_labels
        @test raw_component_summary.center_labels == center_labels
        @test component.centers == nuclear_centers
        @test raw_component.centers == nuclear_centers
        @test component.nuclear_charges == nuclear_charges
        @test raw_component.nuclear_charges == nuclear_charges
        @test component_summary.nuclear_charges == nuclear_charges
        @test raw_component_summary.nuclear_charges == nuclear_charges
        @test component_summary.center_count == length(nuclear_centers)
        @test raw_component_summary.center_count == length(nuclear_centers)
        @test component.ida_term_count == length(term_coefficients)
        @test raw_component.ida_term_count == length(term_coefficients)
        @test component_summary.ida_term_count == length(term_coefficients)
        @test raw_component_summary.ida_term_count == length(term_coefficients)
        @test component.pair_factor_normalization == :density_normalized
        @test raw_component.pair_factor_normalization == :raw_weighted
        @test component_summary.pair_factor_normalization ==
              :density_normalized
        @test raw_component_summary.pair_factor_normalization ==
              :raw_weighted
        @test component.output_finite
        @test raw_component.output_finite
        @test component_summary.finite_checks.output_finite
        @test component_summary.finite_checks.nuclear_output_finite
        @test component_summary.finite_checks.electron_electron_output_finite
        @test raw_component_summary.finite_checks.output_finite
        @test raw_component_summary.finite_checks.nuclear_output_finite
        @test raw_component_summary.finite_checks.electron_electron_output_finite
        @test all(isfinite, component.nuclear_total_matrix)
        @test all(isfinite, raw_component.nuclear_total_matrix)
        @test all(isfinite, component.electron_electron_matrix)
        @test all(isfinite, raw_component.electron_electron_matrix)
        @test component.nuclear_total_matrix ≈ expected_total atol = 1.0e-12 rtol = 1.0e-12
        @test raw_component.nuclear_total_matrix ≈ expected_total atol = 1.0e-12 rtol = 1.0e-12
        @test raw_component.nuclear_total_matrix ≈ component.nuclear_total_matrix atol = 1.0e-12 rtol = 1.0e-12
        @test component.nuclear_symmetry_error <= 1.0e-10
        @test raw_component.nuclear_symmetry_error <= 1.0e-10
        @test component_summary.symmetry_errors.nuclear ==
              component.nuclear_symmetry_error
        @test raw_component_summary.symmetry_errors.nuclear ==
              raw_component.nuclear_symmetry_error
        @test component.electron_electron_matrix ≈ direct_ida.block atol = 1.0e-12 rtol = 1.0e-12
        @test raw_component.electron_electron_matrix ≈ raw_direct_ida.block atol = 1.0e-12 rtol = 1.0e-12
        @test raw_component.electron_electron_matrix ≈ component.electron_electron_matrix atol = 1.0e-12 rtol = 1.0e-12
        @test component.electron_electron_symmetry_error <= 1.0e-10
        @test raw_component.electron_electron_symmetry_error <= 1.0e-10
        @test component_summary.symmetry_errors.electron_electron ==
              component.electron_electron_symmetry_error
        @test raw_component_summary.symmetry_errors.electron_electron ==
              raw_component.electron_electron_symmetry_error
        @test component.dense_parent_ida_authority.object_kind ==
              :pqs_pqs_product_dense_parent_ida_authority_comparison
        @test component.dense_parent_ida_authority.within_tolerance
        @test component.dense_parent_ida_authority.max_error <= 1.0e-10
        @test raw_component.dense_parent_ida_authority === nothing
        @test component_summary.dense_parent_ida_authority.available
        @test component_summary.dense_parent_ida_authority.max_error <=
              1.0e-10
        @test component_summary.dense_parent_ida_authority.within_tolerance
        @test component_summary.dense_parent_ida_authority.skip_reason ===
              nothing
        @test component_summary.dense_parent_ida_authority.validation_only
        @test !raw_component_summary.dense_parent_ida_authority.available
        @test raw_component_summary.dense_parent_ida_authority.max_error ===
              nothing
        @test raw_component_summary.dense_parent_ida_authority.within_tolerance ===
              nothing
        @test raw_component_summary.dense_parent_ida_authority.skip_reason ==
              :density_normalized_authority_only
        @test !raw_component_summary.dense_parent_ida_authority.validation_only
        @test component.authority_max_errors.electron_electron_dense_parent <=
              1.0e-10
        @test raw_component.authority_max_errors.electron_electron_dense_parent ===
              nothing
        @test component.authority_max_errors.nuclear_total_from_center_sum <=
              1.0e-14
        @test raw_component.authority_max_errors.nuclear_total_from_center_sum <=
              1.0e-14
        @test component_summary.nuclear_total_from_center_error <= 1.0e-14
        @test raw_component_summary.nuclear_total_from_center_error <= 1.0e-14

        @test component.nuclear_attraction_by_center.pair_count == 6
        @test raw_component.nuclear_attraction_by_center.pair_count == 6
        @test component_summary.nuclear_pair_count == 6
        @test raw_component_summary.nuclear_pair_count == 6
        @test component.nuclear_attraction_by_center.pair_family_counts ==
              (pqs_pqs = 3, pqs_product = 2, product_product = 1)
        @test raw_component.nuclear_attraction_by_center.pair_family_counts ==
              (pqs_pqs = 3, pqs_product = 2, product_product = 1)
        @test component_summary.nuclear_pair_family_counts ==
              (pqs_pqs = 3, pqs_product = 2, product_product = 1)
        @test raw_component_summary.nuclear_pair_family_counts ==
              (pqs_pqs = 3, pqs_product = 2, product_product = 1)
        @test component.nuclear_attraction_by_center.diagnostics.helper_used_for_pair_families ==
              (
                  pqs_pqs = :_pqs_pqs_source_box_nuclear_attraction_by_center,
                  pqs_product = :_pqs_product_source_box_nuclear_attraction_by_center,
                  product_product =
                      :_product_doside_source_box_nuclear_attraction_by_center,
              )
        @test raw_component.nuclear_attraction_by_center.diagnostics.helper_used_for_pair_families ==
              component.nuclear_attraction_by_center.diagnostics.helper_used_for_pair_families
        @test component.electron_electron_density_density.pair_count == 6
        @test raw_component.electron_electron_density_density.pair_count == 6
        @test component_summary.electron_electron_pair_count == 6
        @test raw_component_summary.electron_electron_pair_count == 6
        @test component.electron_electron_density_density.pair_family_counts ==
              (pqs_pqs = 3, pqs_product = 2, product_product = 1)
        @test raw_component.electron_electron_density_density.pair_family_counts ==
              (pqs_pqs = 3, pqs_product = 2, product_product = 1)
        @test component_summary.electron_electron_pair_family_counts ==
              (pqs_pqs = 3, pqs_product = 2, product_product = 1)
        @test raw_component_summary.electron_electron_pair_family_counts ==
              (pqs_pqs = 3, pqs_product = 2, product_product = 1)
        @test component.electron_electron_density_density.diagnostics.helper_used_for_pair_families ==
              (
                  pqs_pqs =
                      :_pqs_pqs_source_box_density_density_interaction_block,
                  pqs_product =
                      :_pqs_product_source_box_density_density_interaction_block,
                  product_product =
                      :_product_doside_source_box_density_density_interaction_block,
              )
        @test raw_component.electron_electron_density_density.diagnostics.helper_used_for_pair_families ==
              component.electron_electron_density_density.diagnostics.helper_used_for_pair_families
        @test raw_component.electron_electron_density_density.pair_factor_normalization ==
              :raw_weighted
        @test raw_component.electron_electron_density_density.diagnostics.source_weight_division_owner ==
              :pgdg_auxiliary_source_weights
        @test raw_component.electron_electron_density_density.diagnostics.source_weight_division_applied_by_helper
        @test raw_component.electron_electron_density_density.diagnostics.source_weight_division_shape ==
              :axis_pair_weight_outer
        @test component_summary.source_weight_division_owner ==
              :caller_supplied_density_normalized_pair_factors
        @test !component_summary.source_weight_division_applied_by_helper
        @test component_summary.source_weight_division_shape === nothing
        @test raw_component_summary.source_weight_division_owner ==
              :pgdg_auxiliary_source_weights
        @test raw_component_summary.source_weight_division_applied_by_helper
        @test raw_component_summary.source_weight_division_shape ==
              :axis_pair_weight_outer

        @test component.diagnostics.private_component_route_smoke
        @test raw_component.diagnostics.private_component_route_smoke
        @test component.diagnostics.route_shape == (:pqs_left, :pqs_right, :product)
        @test raw_component.diagnostics.route_shape ==
              (:pqs_left, :pqs_right, :product)
        @test component.diagnostics.retained_dimension == 221
        @test raw_component.diagnostics.retained_dimension == 221
        @test component.diagnostics.unit_roles == (:pqs_left, :pqs_right, :product)
        @test raw_component.diagnostics.unit_roles ==
              (:pqs_left, :pqs_right, :product)
        @test component.diagnostics.center_labels == center_labels
        @test raw_component.diagnostics.center_labels == center_labels
        @test component.diagnostics.nuclear_charges == nuclear_charges
        @test raw_component.diagnostics.nuclear_charges == nuclear_charges
        @test component.diagnostics.by_center_nuclear_terms_preserved
        @test raw_component.diagnostics.by_center_nuclear_terms_preserved
        @test component.diagnostics.nuclear_total_is_explicit_sum_of_center_pieces
        @test raw_component.diagnostics.nuclear_total_is_explicit_sum_of_center_pieces
        @test component.diagnostics.ida_term_count == length(term_coefficients)
        @test raw_component.diagnostics.ida_term_count == length(term_coefficients)
        @test component.diagnostics.pair_factor_normalization ==
              :density_normalized
        @test component.diagnostics.density_normalized_pair_factors
        @test !component.diagnostics.raw_weighted_pair_factors
        @test !component.diagnostics.density_normalized_pair_factors_generated
        @test component.diagnostics.source_weight_division_owner ==
              :caller_supplied_density_normalized_pair_factors
        @test !component.diagnostics.source_weight_division_applied_by_helper
        @test component.diagnostics.source_weight_division_shape === nothing
        @test component.diagnostics.source_weights_are_raw_source_weights
        @test raw_component.diagnostics.pair_factor_normalization ==
              :raw_weighted
        @test !raw_component.diagnostics.density_normalized_pair_factors
        @test raw_component.diagnostics.raw_weighted_pair_factors
        @test raw_component.diagnostics.density_normalized_pair_factors_generated
        @test raw_component.diagnostics.source_weight_division_owner ==
              :pgdg_auxiliary_source_weights
        @test raw_component.diagnostics.source_weight_division_applied_by_helper
        @test raw_component.diagnostics.source_weight_division_shape ==
              :axis_pair_weight_outer
        @test raw_component.diagnostics.source_weights_are_raw_source_weights
        @test component.diagnostics.input_pair_factor_data ==
              :ida_gausslet_source_box_provenance
        @test raw_component.diagnostics.input_pair_factor_data ==
              :ida_gausslet_source_box_provenance
        @test component.diagnostics.interaction_path == :ida_gausslet_source_box
        @test raw_component.diagnostics.interaction_path == :ida_gausslet_source_box
        @test !component.diagnostics.mwg_supplement_residual_path
        @test !raw_component.diagnostics.mwg_supplement_residual_path
        @test component.diagnostics.real_ida_gausslet_source_box_provenance_adapted
        @test raw_component.diagnostics.real_ida_gausslet_source_box_provenance_adapted
        @test !component.diagnostics.real_mwg_ida_pair_factor_provenance_adapted
        @test !raw_component.diagnostics.real_mwg_ida_pair_factor_provenance_adapted
        @test component.diagnostics.source_box_first
        @test raw_component.diagnostics.source_box_first
        @test component.diagnostics.source_box_algorithmic_path_true_for_every_pair
        @test raw_component.diagnostics.source_box_algorithmic_path_true_for_every_pair
        @test component.diagnostics.nuclear_source_box_algorithmic_pair_count == 6
        @test raw_component.diagnostics.nuclear_source_box_algorithmic_pair_count == 6
        @test component.diagnostics.electron_electron_source_box_algorithmic_pair_count == 6
        @test raw_component.diagnostics.electron_electron_source_box_algorithmic_pair_count == 6
        @test component.diagnostics.dense_parent_projection_validation_only
        @test !raw_component.diagnostics.dense_parent_projection_validation_only
        @test component.diagnostics.dense_parent_ida_authority_available
        @test component.diagnostics.dense_parent_ida_authority_supported_for_pair_factor_normalization
        @test component.diagnostics.dense_parent_ida_authority_skipped_reason ===
              nothing
        @test !raw_component.diagnostics.dense_parent_ida_authority_available
        @test !raw_component.diagnostics.dense_parent_ida_authority_supported_for_pair_factor_normalization
        @test raw_component.diagnostics.dense_parent_ida_authority_skipped_reason ==
              :density_normalized_authority_only
        @test !component.diagnostics.dense_parent_projection_algorithmic
        @test !raw_component.diagnostics.dense_parent_projection_algorithmic
        @test component.diagnostics.electron_electron_output_representation ==
              :two_index_density_density
        @test raw_component.diagnostics.electron_electron_output_representation ==
              :two_index_density_density
        @test !component.diagnostics.electron_electron_four_index_galerkin_tensor
        @test !raw_component.diagnostics.electron_electron_four_index_galerkin_tensor
        @test !component.diagnostics.hamiltonian_matrix_built
        @test !raw_component.diagnostics.hamiltonian_matrix_built
        @test !component.diagnostics.packet_adoption
        @test !raw_component.diagnostics.packet_adoption
        @test !component.diagnostics.fixed_block_routing
        @test !raw_component.diagnostics.fixed_block_routing
        @test !component.diagnostics.qwhamiltonian_consumes
        @test !raw_component.diagnostics.qwhamiltonian_consumes
        @test !component.diagnostics.public_default_consumes
        @test !raw_component.diagnostics.public_default_consumes
        @test !component.diagnostics.ecp_terms_implemented
        @test !raw_component.diagnostics.ecp_terms_implemented
        @test !component.diagnostics.cr2_science_status_changed
        @test !raw_component.diagnostics.cr2_science_status_changed
        @test !component.diagnostics.mwg_interaction_implemented
        @test !raw_component.diagnostics.mwg_interaction_implemented
        @test !component.diagnostics.mwg_supplement_residual_provenance_adapted
        @test !raw_component.diagnostics.mwg_supplement_residual_provenance_adapted
        @test !component.diagnostics.ida_mwg_semantics_changed
        @test !raw_component.diagnostics.ida_mwg_semantics_changed
        @test !component.diagnostics.mwg_ida_semantics_changed
        @test !raw_component.diagnostics.mwg_ida_semantics_changed
        @test !component.diagnostics.retained_pqs_weights_used
        @test !raw_component.diagnostics.retained_pqs_weights_used
        @test !component.diagnostics.retained_pqs_weights_positive_checked
        @test !raw_component.diagnostics.retained_pqs_weights_positive_checked
        @test !component.diagnostics.retained_weight_division_allowed
        @test !raw_component.diagnostics.retained_weight_division_allowed
        @test !component.diagnostics.retained_pqs_weight_division_allowed
        @test !raw_component.diagnostics.retained_pqs_weight_division_allowed
        @test !component.diagnostics.ida_weight_division_allowed
        @test !raw_component.diagnostics.ida_weight_division_allowed
        @test component.diagnostics.retained_weight_semantics ==
              :not_positive_quadrature_weights
        @test raw_component.diagnostics.retained_weight_semantics ==
              :not_positive_quadrature_weights
        @test !component.diagnostics.shell_projection_used
        @test !raw_component.diagnostics.shell_projection_used
        @test !component.diagnostics.lowdin_cleanup_used
        @test !raw_component.diagnostics.lowdin_cleanup_used
        @test !component.diagnostics.support_local_oracle_used
        @test !raw_component.diagnostics.support_local_oracle_used
        @test !component.diagnostics.support_local_pqs_oracle_used
        @test !raw_component.diagnostics.support_local_pqs_oracle_used
        @test !component.diagnostics.support_local_shell_row_algorithm
        @test !raw_component.diagnostics.support_local_shell_row_algorithm
        @test !component.diagnostics.support_coefficient_matrix_used
        @test !raw_component.diagnostics.support_coefficient_matrix_used
        @test !component.diagnostics.shell_row_algorithm
        @test !raw_component.diagnostics.shell_row_algorithm
        for summary in (component_summary, raw_component_summary)
            @test summary.no_go_diagnostics.source_box_first
            @test summary.no_go_diagnostics.source_box_algorithmic_path_true_for_every_pair
            @test !summary.no_go_diagnostics.shell_projection_used
            @test !summary.no_go_diagnostics.lowdin_cleanup_used
            @test !summary.no_go_diagnostics.support_local_oracle_used
            @test !summary.no_go_diagnostics.support_local_pqs_oracle_used
            @test !summary.no_go_diagnostics.support_local_shell_row_algorithm
            @test !summary.no_go_diagnostics.support_coefficient_matrix_used
            @test !summary.no_go_diagnostics.shell_row_algorithm
            @test !summary.no_go_diagnostics.retained_pqs_weights_used
            @test !summary.no_go_diagnostics.retained_pqs_weights_positive_checked
            @test !summary.no_go_diagnostics.retained_weight_division_allowed
            @test !summary.no_go_diagnostics.retained_pqs_weight_division_allowed
            @test !summary.no_go_diagnostics.ida_weight_division_allowed
            @test !summary.no_go_diagnostics.packet_adoption
            @test !summary.no_go_diagnostics.fixed_block_routing
            @test !summary.no_go_diagnostics.qwhamiltonian_consumes
            @test !summary.no_go_diagnostics.hamiltonian_matrix_built
            @test !summary.no_go_diagnostics.public_default_consumes
            @test !summary.no_go_diagnostics.mwg_supplement_residual_path
            @test !summary.no_go_diagnostics.mwg_supplement_residual_provenance_adapted
            @test !summary.no_go_diagnostics.ecp_terms_implemented
            @test !summary.no_go_diagnostics.cr2_science_status_changed
            @test !summary.no_go_diagnostics.dense_parent_projection_algorithmic
            @test summary.performance.nuclear.elapsed_seconds >= 0.0
            @test summary.performance.nuclear.allocated_bytes >= 0
            @test summary.performance.nuclear.gc_time_seconds >= 0.0
            @test summary.performance.electron_electron.elapsed_seconds >= 0.0
            @test summary.performance.electron_electron.allocated_bytes >= 0
            @test summary.performance.electron_electron.gc_time_seconds >= 0.0
            @test summary.diagnostics.private_component_route_smoke_summary
            @test summary.diagnostics.source_box_first
            @test summary.diagnostics.source_box_algorithmic_path_true_for_every_pair
            @test !summary.diagnostics.retained_pqs_weights_used
            @test !summary.diagnostics.retained_weight_division_allowed
            @test !summary.diagnostics.ida_weight_division_allowed
            @test !summary.diagnostics.packet_adoption
            @test !summary.diagnostics.fixed_block_routing
            @test !summary.diagnostics.qwhamiltonian_consumes
            @test !summary.diagnostics.public_default_consumes
            @test !summary.diagnostics.mwg_supplement_residual_provenance_adapted
            @test !summary.diagnostics.ecp_terms_implemented
            @test !summary.diagnostics.cr2_science_status_changed
        end
    end

    function _product_staged_comparison_axis_row(axis, state_index::Int)
        axis.kind == :fixed && return 1
        axis.kind == :active && return state_index - first(axis.interval) + 1
        throw(ArgumentError("unsupported staged axis kind $(axis.kind)"))
    end

    function _product_staged_comparison_coefficient_matrix(
        support_states::AbstractVector{<:NTuple{3,Int}},
        axes::NTuple{3,Any},
        axis_function_indices::AbstractVector{<:NTuple{3,Int}},
    )
        coefficients = zeros(Float64, length(support_states), length(axis_function_indices))
        for (column, axis_function_index) in pairs(axis_function_indices)
            for (row, state) in pairs(support_states)
                value = 1.0
                for axis_index in 1:3
                    axis = axes[axis_index]
                    axis_row =
                        _product_staged_comparison_axis_row(axis, state[axis_index])
                    value *= axis.coefficient_matrix[axis_row, axis_function_index[axis_index]]
                end
                coefficients[row, column] = value
            end
        end
        return coefficients
    end

    function _product_staged_comparison_unit_with_range(unit, column_range::UnitRange{Int})
        length(column_range) == length(unit.column_range) || throw(
            ArgumentError("test product-staged column range must preserve retained count"),
        )
        return GaussletBases._CartesianNestedProductStagedByCenterUnit3D(
            unit.role,
            unit.kind,
            column_range,
            copy(unit.support_indices),
            copy(unit.support_states),
            unit.coefficient_matrix,
            unit.axes,
            copy(unit.axis_function_indices),
            (
                source = :product_staged_metric_comparison_test_packaging,
                original_provenance = unit.provenance,
                original_column_range = unit.column_range,
            ),
            (
                source = :product_staged_metric_comparison_test_packaging,
                original_diagnostics = unit.diagnostics,
            ),
        )
    end

    function _product_staged_comparison_retained_blocks(left_unit, right_unit, axis_metrics)
        left_range = 1:length(left_unit.column_range)
        same_unit = left_unit === right_unit
        if same_unit
            left_packaged = _product_staged_comparison_unit_with_range(left_unit, left_range)
            right_packaged = left_packaged
            right_range = left_range
            dimension = length(left_range)
        else
            right_range =
                (last(left_range) + 1):(last(left_range) + length(right_unit.column_range))
            left_packaged = _product_staged_comparison_unit_with_range(left_unit, left_range)
            right_packaged = _product_staged_comparison_unit_with_range(right_unit, right_range)
            dimension = last(right_range)
        end
        overlap = zeros(Float64, dimension, dimension)
        position_x = zeros(Float64, dimension, dimension)
        position_y = zeros(Float64, dimension, dimension)
        position_z = zeros(Float64, dimension, dimension)
        GaussletBases.CartesianContractedParentMetrics._fill_product_staged_metric_blocks!(
            overlap,
            position_x,
            position_y,
            position_z,
            left_packaged,
            right_packaged,
            axis_metrics,
        )
        return (
            helper_path = :fill_product_staged_metric_blocks,
            fixture_scope = :private_test_only,
            overlap = Matrix{Float64}(overlap[left_range, right_range]),
            position_x = Matrix{Float64}(position_x[left_range, right_range]),
            position_y = Matrix{Float64}(position_y[left_range, right_range]),
            position_z = Matrix{Float64}(position_z[left_range, right_range]),
        )
    end

    function _product_staged_comparison_block_for_term(blocks, term::Symbol)
        term == :overlap && return blocks.overlap
        term == :position_x && return blocks.position_x
        term == :position_y && return blocks.position_y
        term == :position_z && return blocks.position_z
        throw(ArgumentError("unsupported product-staged comparison term $term"))
    end

    expansion = coulomb_gaussian_expansion(doacc = false)
    term_coefficients = Float64.(expansion.coefficients)
    bundle5 = _pqs_test_bundle(5)
    bundle7 = _pqs_test_bundle(7)
    CCP = GaussletBases.CartesianContractedParents
    CCPM = GaussletBases.CartesianContractedParentMetrics

    cubic_bundles = GaussletBases._CartesianNestedAxisBundles3D(bundle5, bundle5, bundle5)
    cubic_current = (1:5, 1:5, 1:5)
    cubic_inner = (2:4, 2:4, 2:4)
    cubic_expected = setdiff(
        GaussletBases._nested_box_support_indices(cubic_current..., (5, 5, 5)),
        GaussletBases._nested_box_support_indices(cubic_inner..., (5, 5, 5)),
    )
    sort!(cubic_expected)
    cubic = GaussletBases._nested_projected_q_shell_layer(
        cubic_bundles,
        cubic_current,
        cubic_inner;
        bond_axis = :z,
        q = 5,
        L = 5,
        term_coefficients,
    )

    @test cubic isa GaussletBases._CartesianNestedProjectedQShellLayer3D
    @test cubic.support_indices == cubic_expected
    @test cubic.diagnostics.support_contract == :projected_q_shell_raw_boundary
    @test cubic.diagnostics.coefficient_contract ==
          :full_block_boundary_comx_product_mode_projection
    @test cubic.diagnostics.primitive_family == :projected_q_shell
    @test cubic.diagnostics.boundary_support_count == 98
    @test cubic.diagnostics.full_block_column_count == 125
    @test cubic.diagnostics.boundary_comx_product_mode_count == 98
    @test cubic.diagnostics.boundary_comx_product_mode_selection_rule ==
          :any_axis_mode_index_first_or_last
    @test cubic.diagnostics.retained_count == 98
    @test cubic.diagnostics.rank_count == 98
    @test cubic.diagnostics.rank_drop_count == 0
    @test cubic.diagnostics.duplicate_count == 0
    @test cubic.diagnostics.missing_count == 0
    @test cubic.diagnostics.outside_count == 0
    @test cubic.diagnostics.coverage_ok
    @test cubic.diagnostics.cleanup_method == :projected_boundary_symmetric_lowdin
    @test cubic.diagnostics.overlap_error < 1.0e-8
    @test cubic.diagnostics.packet_overlap_error < 1.0e-8
    @test !cubic.diagnostics.dense_parent_matrix_used
    @test !cubic.diagnostics.endcap_panel_stitching
    @test !cubic.diagnostics.active_builder_consumes
    @test cubic.provenance.construction_contract ==
          :raw_boundary_projection_of_boundary_comx_product_modes_from_full_local_block_transform
    @test :contracted_inner_cube_subtraction in cubic.provenance.not_contract
    @test :locked_prior_span_projection in cubic.provenance.not_contract
    @test size(cubic.coefficient_matrix) == (125, 98)
    @test all(isfinite, Matrix{Float64}(cubic.coefficient_matrix[cubic.support_indices, :]))
    @test !_pqs_one_hot_selector_columns(cubic.coefficient_matrix[cubic.support_indices, :])
    @test norm(cubic.packet.overlap - I, Inf) < 1.0e-8
    @test all(isfinite, cubic.packet.weights)
    @test minimum(cubic.packet.weights) > -1.0e-10
    cubic_descriptor = GaussletBases._nested_projected_q_shell_staged_unit_descriptor(cubic)
    @test cubic_descriptor.kind == :projected_q_shell
    @test isnothing(cubic_descriptor.role)
    @test cubic_descriptor.support_indices == cubic.support_indices
    @test cubic_descriptor.support_states == cubic.support_states
    @test cubic_descriptor.current_box == cubic_current
    @test cubic_descriptor.inner_box == cubic_inner
    @test cubic_descriptor.bond_axis == :z
    @test cubic_descriptor.q == 5
    @test cubic_descriptor.L == 5
    @test cubic_descriptor.support_count == 98
    @test cubic_descriptor.mode_count == 98
    @test cubic_descriptor.retained_count == 98
    @test length(cubic_descriptor.boundary_mode_indices) == 98
    @test all(
        mode -> _pqs_boundary_mode_rule_ok(mode, (5, 5, 5)),
        cubic_descriptor.boundary_mode_indices,
    )
    @test cubic_descriptor.selection_rule == :any_axis_mode_index_first_or_last
    @test cubic_descriptor.cleanup_method == :projected_boundary_symmetric_lowdin
    @test cubic_descriptor.cleanup_matrix_size == (98, 98)
    @test cubic_descriptor.cleanup_rank_count == 98
    @test cubic_descriptor.cleanup_rank_drop_count == 0
    @test length(cubic_descriptor.cleanup_eigenvalues) == 98
    @test all(>(cubic_descriptor.cleanup_cutoff), cubic_descriptor.cleanup_eigenvalues)
    @test cubic_descriptor.support_local_coefficient_shape == (98, 98)
    @test :product_doside_unit in cubic_descriptor.non_contracts
    @test :dense_full_parent_fallback in cubic_descriptor.non_contracts
    @test cubic_descriptor.diagnostics.metadata_only
    @test !cubic_descriptor.diagnostics.product_doside_unit
    @test !cubic_descriptor.active_consumption.fixed_block_sidecar_installed
    @test !cubic_descriptor.active_consumption.metric_packet_consumes
    @test !cubic_descriptor.active_consumption.by_center_consumes
    cubic_shared_raw_product_box_plan = GaussletBases._cartesian_raw_product_box_plan(
        cubic_bundles,
        cubic_descriptor.axis_intervals,
        (5, 5, 5);
        enforce_symmetric_odd = false,
    )
    @test_throws DimensionMismatch CCPM._pqs_raw_product_box_plan(
        cubic_descriptor,
        merge(cubic_shared_raw_product_box_plan, (source_mode_dims = (5, 5, 4),)),
    )
    @test_throws DimensionMismatch CCPM._pqs_raw_product_box_plan(
        cubic_descriptor,
        merge(cubic_shared_raw_product_box_plan, (source_mode_count = 124,)),
    )
    @test_throws DimensionMismatch CCPM._pqs_raw_product_box_plan(
        cubic_descriptor,
        merge(cubic_shared_raw_product_box_plan, (axis_intervals = (1:5, 1:5, 1:4),)),
    )
    bad_cubic_shared_coefficients = let
        bad_x = copy(cubic_shared_raw_product_box_plan.axis_local_coefficients[1])
        bad_x[1, 1] += 1.0e-3
        (
            bad_x,
            cubic_shared_raw_product_box_plan.axis_local_coefficients[2],
            cubic_shared_raw_product_box_plan.axis_local_coefficients[3],
        )
    end
    @test_throws ArgumentError CCPM._pqs_raw_product_box_plan(
        cubic_descriptor,
        merge(
            cubic_shared_raw_product_box_plan,
            (axis_local_coefficients = bad_cubic_shared_coefficients,),
        ),
    )
    cubic_region = CCP.cartesian_shell_region(
        cubic_descriptor;
        parent_dimension = 5 * 5 * 5,
    )
    @test cubic_region isa CCP.CartesianShellRegion3D
    @test cubic_region.region_family == :projected_q_shell_boundary_modes
    @test cubic_region.status == :prototype
    @test cubic_region.box == cubic_current
    @test cubic_region.inner_exclusion_box == cubic_inner
    @test cubic_region.support_indices == cubic.support_indices
    @test cubic_region.support_summary.entry_count == 98
    @test cubic_region.support_summary.missing_count == 27
    @test cubic_region.ownership_coverage_contract == :boundary_only
    @test cubic_region.retention.retention_rule ==
          :boundary_comx_product_mode_selection
    @test cubic_region.retention.cleanup_rule == :shell_projection_lowdin
    @test cubic_region.retention.preferred_contraction_rule ==
          :shell_realized_support_local_fixture
    @test cubic_region.retention.expected_unit_family == :projected_q_shell_descriptor
    @test cubic_region.retention.metric_capability == :pqs_low_order_product_metric_prototype
    @test cubic_region.retention.diagnostics.raw_box_retained_rule ==
          :boundary_comx_product_mode_selection
    @test cubic_region.retention.diagnostics.shell_realization_rule ==
          :shell_projection_lowdin
    @test cubic_region.retention.diagnostics.representation_stage ==
          :shell_realized_support_local_fixture
    @test !cubic_region.retention.diagnostics.raw_product_box_operator_contract
    @test cubic_region.retention.diagnostics.shell_projection_lowdin_realization
    @test cubic_region.retention.diagnostics.raw_box_boundary_selector_transform
    @test :fixed_block_sidecar_payload in cubic_region.retention.missing_payload_fields
    @test !cubic_region.current_route_consumes
    @test !cubic_region.descriptor_drives_builder
    @test cubic_region.descriptor_only
    @test cubic_region.diagnostics.prototype_only
    @test !cubic_region.diagnostics.fixed_block_sidecar_installed
    @test cubic_region.geometry.boundary_mode_count == 98
    @test cubic_region.geometry.selection_rule == :any_axis_mode_index_first_or_last
    cubic_resolved_payload = CCP._cartesian_resolved_contraction_payload(
        cubic_descriptor;
        parent_dimension = 5 * 5 * 5,
    )
    @test cubic_resolved_payload isa CCP._CartesianResolvedContractionPayload3D
    @test cubic_resolved_payload.metric_path == :unsupported_prototype
    @test !cubic_resolved_payload.ready_for_metric_execution
    @test cubic_resolved_payload.payload_kind == :projected_q_shell
    @test isnothing(cubic_resolved_payload.column_range)
    @test cubic_resolved_payload.payload === cubic_descriptor
    @test :fixed_block_sidecar_payload in cubic_resolved_payload.missing_fields
    @test :product_staged_metric_payload in cubic_resolved_payload.missing_fields
    @test cubic_resolved_payload.diagnostics.prototype
    @test cubic_resolved_payload.diagnostics.metric_capability ==
          :pqs_low_order_product_metric_prototype
    cubic_rule = CCP.cartesian_contraction_rule(
        cubic_descriptor;
        parent_dimension = 5 * 5 * 5,
    )
    @test cubic_rule isa CCP.CartesianContractionRule3D
    @test CCP.contraction_rule_family(cubic_rule) == :projected_q_shell_boundary_modes
    @test CCP.contraction_rule_kind(cubic_rule) == :projected_q_shell
    @test isnothing(cubic_rule.role)
    @test cubic_rule.support_indices == cubic.support_indices
    @test cubic_rule.support_summary.entry_count == 98
    @test cubic_rule.support_summary.unique_count == 98
    @test cubic_rule.support_summary.duplicate_count == 0
    @test cubic_rule.support_summary.outside_count == 0
    @test cubic_rule.support_summary.missing_count == 27
    @test !cubic_rule.support_summary.support_complete
    @test cubic_rule.local_geometry.current_box == cubic_current
    @test cubic_rule.local_geometry.inner_box == cubic_inner
    @test cubic_rule.local_geometry.axis_local_coefficient_shapes == ((5, 5), (5, 5), (5, 5))
    @test CCP.contraction_rule_column_range(cubic_rule) === nothing
    @test CCP.contraction_rule_source_dimension(cubic_rule) == 125
    @test CCP.contraction_rule_retained_dimension(cubic_rule) == 98
    @test CCP.contraction_rule_transform_rule(cubic_rule) ==
          :boundary_comx_product_mode_selection
    @test CCP.contraction_rule_cleanup_rule(cubic_rule) == :shell_projection_lowdin
    @test CCP.contraction_rule_metric_capability(cubic_rule) ==
          :pqs_low_order_product_metric_prototype
    @test cubic_rule.diagnostics.raw_box_retained_rule ==
          :boundary_comx_product_mode_selection
    @test cubic_rule.diagnostics.shell_realization_rule == :shell_projection_lowdin
    @test cubic_rule.diagnostics.representation_stage ==
          :shell_realized_support_local_fixture
    @test !cubic_rule.diagnostics.raw_product_box_operator_contract
    @test cubic_rule.diagnostics.shell_projection_lowdin_realization
    @test cubic_rule.diagnostics.boundary_mode_count == 98
    @test cubic_rule.diagnostics.cleanup_rank_drop_count == 0
    cubic_rule_inventory = CCP.cartesian_contraction_rule_inventory(
        [cubic_rule];
        parent_dimension = 5 * 5 * 5,
        contracted_dimension = nothing,
        unit_count = 0,
        every_unit_has_rule_metadata = false,
        every_unit_rule_derivable = false,
        provenance = (source = :pqs_descriptor_rule_path,),
    )
    @test cubic_rule_inventory.rule_count == 1
    @test cubic_rule_inventory.unit_count == 0
    @test cubic_rule_inventory.contracted_dimension === nothing
    @test Dict(cubic_rule_inventory.rule_family_counts) ==
          Dict(:projected_q_shell_boundary_modes => 1)
    @test cubic_rule_inventory.metric_capabilities ==
          [:pqs_low_order_product_metric_prototype]
    @test cubic_rule_inventory.total_source_dimension == 125
    @test cubic_rule_inventory.total_retained_dimension == 98
    @test cubic_rule_inventory.support_summary.entry_count == 98
    @test cubic_rule_inventory.support_summary.missing_count == 27
    @test !cubic_rule_inventory.every_unit_has_rule_metadata
    @test !cubic_rule_inventory.every_unit_rule_derivable
    @test cubic_rule_inventory.any_metadata_only_rule
    @test cubic_rule_inventory.any_prototype_rule
    @test !cubic_rule_inventory.diagnostics.parent_level_unit_inventory
    @test !cubic_rule_inventory.diagnostics.all_rules_have_column_ranges
    @test cubic_rule_inventory.diagnostics.q_shell_rule_present
    @test !cubic_rule_inventory.diagnostics.q_shell_installed_as_contracted_parent_unit
    cubic_rule_dispatch =
        CCPM._contracted_parent_metric_dispatch_plan_from_rules([cubic_rule])
    @test cubic_rule_dispatch.unit_count == 1
    @test cubic_rule_dispatch.unsupported_unit_count == 1
    @test cubic_rule_dispatch.prototype_rule_count == 1
    @test cubic_rule_dispatch.unsupported_block_count == 1
    @test !cubic_rule_dispatch.plan_supported
    @test only(cubic_rule_dispatch.unit_paths).block_role == :unsupported_prototype
    @test only(cubic_rule_dispatch.block_paths).path == :unsupported
    @test cubic.diagnostics.pqs_staged_unit_descriptor_available
    @test cubic.diagnostics.pqs_staged_unit_kind == :projected_q_shell
    @test cubic.provenance.pqs_staged_unit_descriptor === cubic_descriptor
    cubic_metric_prototype =
        GaussletBases._nested_projected_q_shell_descriptor_metric_prototype(
            cubic,
            cubic_bundles,
        )
    @test cubic_metric_prototype.overlap ≈ cubic.packet.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test cubic_metric_prototype.weights ≈ cubic.packet.weights atol = 1.0e-10 rtol = 1.0e-10
    @test cubic_metric_prototype.diagnostics.descriptor_kind == :projected_q_shell
    @test cubic_metric_prototype.diagnostics.boundary_comx_product_modes_used
    @test cubic_metric_prototype.diagnostics.raw_boundary_projection_used
    @test cubic_metric_prototype.diagnostics.lowdin_cleanup_applied
    @test cubic_metric_prototype.diagnostics.support_local_boundary_matrix_used
    @test !cubic_metric_prototype.diagnostics.slab_decomposed_product_contraction
    @test !cubic_metric_prototype.diagnostics.overlap_invariant_applied
    @test cubic_metric_prototype.diagnostics.overlap_invariant_debug_check
    @test cubic_metric_prototype.diagnostics.overlap_invariant_error < 1.0e-10
    @test !cubic_metric_prototype.diagnostics.overlap_is_operator_target
    @test !cubic_metric_prototype.diagnostics.dense_full_parent_matrix_used
    @test !cubic_metric_prototype.diagnostics.product_doside_unit
    @test !cubic_metric_prototype.diagnostics.fixed_block_sidecar_installed
    @test !cubic_metric_prototype.diagnostics.optimized_sidecar_installed
    cubic_product_metric =
        GaussletBases._nested_projected_q_shell_descriptor_metric_product_contraction(
            cubic,
            cubic_bundles,
        )
    @test cubic_product_metric.overlap == Matrix{Float64}(I, 98, 98)
    @test cubic_product_metric.weights ≈ cubic_metric_prototype.weights atol = 1.0e-10 rtol = 1.0e-10
    @test cubic_product_metric.weights ≈ cubic.packet.weights atol = 1.0e-10 rtol = 1.0e-10
    @test cubic_product_metric.position_x ≈ cubic.packet.position_x atol = 1.0e-10 rtol = 1.0e-10
    @test cubic_product_metric.position_y ≈ cubic.packet.position_y atol = 1.0e-10 rtol = 1.0e-10
    @test cubic_product_metric.position_z ≈ cubic.packet.position_z atol = 1.0e-10 rtol = 1.0e-10
    @test cubic_product_metric.coverage.piece_count == 6
    @test cubic_product_metric.coverage.piece_roles == (:xlo, :xhi, :ylo, :yhi, :zlo, :zhi)
    @test cubic_product_metric.coverage.support_count == 98
    @test cubic_product_metric.coverage.unique_support_count == 98
    @test cubic_product_metric.coverage.duplicate_count == 0
    @test cubic_product_metric.coverage.missing_count == 0
    @test cubic_product_metric.coverage.outside_count == 0
    @test cubic_product_metric.coverage.coverage_ok
    @test cubic_product_metric.diagnostics.descriptor_kind == :projected_q_shell
    @test cubic_product_metric.diagnostics.slab_decomposed_product_contraction
    @test !cubic_product_metric.diagnostics.support_local_boundary_matrix_used
    @test !cubic_product_metric.diagnostics.dense_full_parent_matrix_used
    @test cubic_product_metric.diagnostics.boundary_comx_product_modes_used
    @test cubic_product_metric.diagnostics.raw_boundary_projection_used
    @test cubic_product_metric.diagnostics.lowdin_cleanup_applied
    @test cubic_product_metric.diagnostics.overlap_invariant_applied
    @test !cubic_product_metric.diagnostics.overlap_invariant_debug_check
    @test cubic_product_metric.diagnostics.overlap_invariant_error == 0.0
    @test !cubic_product_metric.diagnostics.overlap_is_operator_target
    @test cubic_product_metric.diagnostics.nontrivial_product_contracted_terms ==
          (:weights, :position_x, :position_y, :position_z)
    @test cubic_product_metric.diagnostics.mode_axis_indices_cached
    @test cubic_product_metric.diagnostics.symmetric_mode_matrix_assembly
    @test cubic_product_metric.diagnostics.axis_piece_blocks_use_views
    @test !cubic_product_metric.diagnostics.product_doside_unit
    @test !cubic_product_metric.diagnostics.fixed_block_sidecar_installed
    @test !cubic_product_metric.diagnostics.optimized_sidecar_installed
    @test cubic_product_metric.diagnostics.prototype_only
    cubic_pqs_payload = CCP._cartesian_executable_projected_q_shell_payload_fixture(
        cubic_descriptor;
        column_range = 1:98,
        parent_dimension = 5 * 5 * 5,
    )
    @test cubic_pqs_payload.kind == :projected_q_shell
    @test cubic_pqs_payload.column_range == 1:98
    @test cubic_pqs_payload.support_indices == cubic.support_indices
    @test cubic_pqs_payload.support_states == cubic.support_states
    @test !hasproperty(cubic_pqs_payload, :coefficient_matrix)
    @test hasproperty(cubic_pqs_payload, :support_coefficient_matrix)
    @test size(cubic_pqs_payload.support_coefficient_matrix) == (98, 98)
    @test size(cubic_pqs_payload.support_coefficient_matrix, 1) ==
          length(cubic_pqs_payload.support_indices)
    @test size(cubic_pqs_payload.support_coefficient_matrix, 1) != 5 * 5 * 5
    cubic_pqs_entries = CCPM._staged_unit_entries(cubic_pqs_payload)
    @test length(cubic_pqs_entries) == length(cubic_pqs_payload.column_range)
    @test sum(length, cubic_pqs_entries) ==
          count(!iszero, Matrix{Float64}(cubic_pqs_payload.support_coefficient_matrix))
    @test all(
        entry -> (entry.ix, entry.iy, entry.iz) in cubic_pqs_payload.support_states,
        Iterators.flatten(cubic_pqs_entries),
    )
    @test cubic_pqs_payload.diagnostics.fixture_only
    @test !cubic_pqs_payload.diagnostics.production_supported
    @test cubic_pqs_payload.diagnostics.coefficient_scope == :support_local_boundary_rows
    @test !cubic_pqs_payload.diagnostics.parent_dimension_coefficient_map
    @test cubic_pqs_payload.diagnostics.support_coefficient_shape == (98, 98)
    @test cubic_pqs_payload.diagnostics.retained_column_weight_role ==
          :debug_reference_only
    @test cubic_pqs_payload.diagnostics.retained_weight_semantics ==
          :debug_reference_only
    @test !cubic_pqs_payload.diagnostics.retained_weights_used_for_ida_division
    @test !cubic_pqs_payload.diagnostics.ida_weight_division_allowed
    @test !cubic_pqs_payload.diagnostics.quadrature_weight_semantics_claimed
    @test cubic_pqs_payload.diagnostics.active_interaction_path == :none_fixture_only
    @test !cubic_pqs_payload.diagnostics.fixed_block_sidecar_installed
    @test !cubic_pqs_payload.diagnostics.default_builder_consumes
    @test !cubic_pqs_payload.diagnostics.pqs_product_optimized_path_ready
    cubic_pqs_sidecar = CCP._cartesian_projected_q_shell_sidecar_fixture(
        cubic_descriptor;
        column_range = 1:98,
        dims = (5, 5, 5),
    )
    @test cubic_pqs_sidecar.dims == (5, 5, 5)
    @test length(cubic_pqs_sidecar.payloads) == 1
    @test only(cubic_pqs_sidecar.payloads).support_coefficient_matrix ==
          cubic_pqs_payload.support_coefficient_matrix
    @test cubic_pqs_sidecar.diagnostics.fixture_only
    @test !cubic_pqs_sidecar.diagnostics.production_supported
    @test !cubic_pqs_sidecar.diagnostics.fixed_block_sidecar_installed
    @test !cubic_pqs_sidecar.diagnostics.default_builder_consumes
    @test !cubic_pqs_sidecar.diagnostics.qw_consumes
    @test !cubic_pqs_sidecar.diagnostics.hamiltonian_consumes
    @test cubic_pqs_sidecar.diagnostics.support_local_reference_only
    @test !cubic_pqs_sidecar.diagnostics.pqs_product_optimized_path_ready
    @test !cubic_pqs_sidecar.provenance.installed_in_fixed_block
    cubic_pqs_resolved = CCP._cartesian_resolved_contraction_payload(
        cubic_pqs_payload;
        parent_dimension = 5 * 5 * 5,
    )
    cubic_pqs_sidecar_resolved = CCP._cartesian_resolved_contraction_payloads(
        cubic_pqs_sidecar,
    )
    @test length(cubic_pqs_sidecar_resolved) == 1
    @test only(cubic_pqs_sidecar_resolved).ready_for_metric_execution
    @test only(cubic_pqs_sidecar_resolved).payload === only(cubic_pqs_sidecar.payloads)
    @test only(cubic_pqs_sidecar_resolved).metric_path ==
          :pqs_low_order_support_local_reference
    @test cubic_pqs_resolved.metric_path == :pqs_low_order_support_local_reference
    @test cubic_pqs_resolved.ready_for_metric_execution
    @test cubic_pqs_resolved.payload_kind == :projected_q_shell
    @test cubic_pqs_resolved.column_range == 1:98
    @test isempty(cubic_pqs_resolved.missing_fields)
    @test cubic_pqs_resolved.diagnostics.block_role == :pqs
    @test cubic_pqs_resolved.diagnostics.metric_capability ==
          :pqs_low_order_support_local_reference
    @test cubic_pqs_resolved.diagnostics.coefficient_scope == :support_local_boundary_rows
    @test !cubic_pqs_resolved.diagnostics.parent_dimension_coefficient_map
    @test cubic_pqs_resolved.diagnostics.retained_column_weight_role ==
          :debug_reference_only
    @test cubic_pqs_resolved.diagnostics.retained_weight_semantics ==
          :debug_reference_only
    @test !cubic_pqs_resolved.diagnostics.retained_weights_used_for_ida_division
    @test !cubic_pqs_resolved.diagnostics.ida_weight_division_allowed
    @test !cubic_pqs_resolved.diagnostics.quadrature_weight_semantics_claimed
    @test cubic_pqs_resolved.diagnostics.active_interaction_path == :none_fixture_only
    @test !cubic_pqs_resolved.diagnostics.production_supported
    @test !cubic_pqs_resolved.diagnostics.fixed_block_sidecar_installed
    pqs_self_dispatch = CCPM._metric_dispatch_plan_from_resolved_payloads(
        [cubic_pqs_resolved],
    )
    @test pqs_self_dispatch.plan_supported
    @test pqs_self_dispatch.pqs_unit_count == 1
    @test pqs_self_dispatch.pqs_pqs_block_count == 1
    @test only(pqs_self_dispatch.block_paths).path == :pqs_pqs_low_order_reference
    cubic_metrics = _pqs_axis_metrics(cubic_bundles)
    _check_pqs_raw_product_box_reference(
        CCPM,
        cubic_descriptor,
        cubic_metrics;
        expected_source_mode_dims = (5, 5, 5),
        expected_retained_count = 98,
        shared_raw_product_box_plan = cubic_shared_raw_product_box_plan,
    )
    cubic_pqs_pqs_source_box_plan = CCPM._pqs_raw_product_box_plan(
        cubic_descriptor,
        cubic_shared_raw_product_box_plan,
        cubic_metrics,
    )
    _check_pqs_pqs_source_box_self_blocks(
        CCPM,
        cubic_pqs_pqs_source_box_plan,
        cubic_metrics;
        expected_source_mode_dims = (5, 5, 5),
        expected_retained_count = 98,
    )
    shifted_cubic_bundles =
        GaussletBases._CartesianNestedAxisBundles3D(bundle7, bundle5, bundle5)
    shifted_left_current = (1:5, 1:5, 1:5)
    shifted_left_inner = (2:4, 2:4, 2:4)
    shifted_right_current = (3:7, 1:5, 1:5)
    shifted_right_inner = (4:6, 2:4, 2:4)
    shifted_left = GaussletBases._nested_projected_q_shell_layer(
        shifted_cubic_bundles,
        shifted_left_current,
        shifted_left_inner;
        bond_axis = :z,
        q = 5,
        L = 5,
        term_coefficients,
    )
    shifted_right = GaussletBases._nested_projected_q_shell_layer(
        shifted_cubic_bundles,
        shifted_right_current,
        shifted_right_inner;
        bond_axis = :z,
        q = 5,
        L = 5,
        term_coefficients,
    )
    shifted_left_descriptor =
        GaussletBases._nested_projected_q_shell_staged_unit_descriptor(
            shifted_left,
        )
    shifted_right_descriptor =
        GaussletBases._nested_projected_q_shell_staged_unit_descriptor(
            shifted_right,
        )
    shifted_left_shared_raw_product_box_plan =
        GaussletBases._cartesian_raw_product_box_plan(
            shifted_cubic_bundles,
            shifted_left_descriptor.axis_intervals,
            (5, 5, 5);
            enforce_symmetric_odd = false,
        )
    shifted_right_shared_raw_product_box_plan =
        GaussletBases._cartesian_raw_product_box_plan(
            shifted_cubic_bundles,
            shifted_right_descriptor.axis_intervals,
            (5, 5, 5);
            enforce_symmetric_odd = false,
        )
    shifted_metrics = _pqs_axis_metrics(shifted_cubic_bundles)
    shifted_left_pqs_pqs_source_box_plan =
        CCPM._pqs_raw_product_box_plan(
            shifted_left_descriptor,
            shifted_left_shared_raw_product_box_plan,
            shifted_metrics,
        )
    shifted_right_pqs_pqs_source_box_plan =
        CCPM._pqs_raw_product_box_plan(
            shifted_right_descriptor,
            shifted_right_shared_raw_product_box_plan,
            shifted_metrics,
        )
    _check_pqs_pqs_source_box_cross_blocks(
        CCPM,
        shifted_left_descriptor,
        shifted_right_descriptor,
        shifted_left_pqs_pqs_source_box_plan,
        shifted_right_pqs_pqs_source_box_plan,
        shifted_metrics;
        expected_source_mode_dims = (5, 5, 5),
        expected_retained_count = 98,
    )
    pqs_self_block = CCPM._resolved_payload_low_order_metric_block(
        cubic_pqs_resolved,
        cubic_pqs_resolved,
        cubic_metrics,
    )
    @test pqs_self_block.path == :pqs_pqs_low_order_reference
    @test pqs_self_block.overlap == Matrix{Float64}(I, 98, 98)
    @test pqs_self_block.weights ≈ cubic_product_metric.weights atol = 1.0e-10 rtol = 1.0e-10
    @test pqs_self_block.position_x ≈ cubic.packet.position_x atol = 1.0e-10 rtol = 1.0e-10
    @test pqs_self_block.position_y ≈ cubic.packet.position_y atol = 1.0e-10 rtol = 1.0e-10
    @test pqs_self_block.position_z ≈ cubic.packet.position_z atol = 1.0e-10 rtol = 1.0e-10
    support_coefficients = Matrix{Float64}(cubic_pqs_payload.support_coefficient_matrix)
    expected_first_moments = hcat(
        vec(
            transpose(support_coefficients) * Float64[
                cubic_metrics.x.centers[state[1]] *
                cubic_metrics.x.weights[state[1]] *
                cubic_metrics.y.weights[state[2]] *
                cubic_metrics.z.weights[state[3]] for state in cubic.support_states
            ],
        ),
        vec(
            transpose(support_coefficients) * Float64[
                cubic_metrics.x.weights[state[1]] *
                cubic_metrics.y.centers[state[2]] *
                cubic_metrics.y.weights[state[2]] *
                cubic_metrics.z.weights[state[3]] for state in cubic.support_states
            ],
        ),
        vec(
            transpose(support_coefficients) * Float64[
                cubic_metrics.x.weights[state[1]] *
                cubic_metrics.y.weights[state[2]] *
                cubic_metrics.z.centers[state[3]] *
                cubic_metrics.z.weights[state[3]] for state in cubic.support_states
            ],
        ),
    )
    @test pqs_self_block.first_moments ≈ expected_first_moments atol = 1.0e-10 rtol = 1.0e-10
    @test pqs_self_block.diagnostics.pqs_self_overlap_invariant_applied
    @test pqs_self_block.diagnostics.support_overlap_debug_error < 1.0e-10
    @test !pqs_self_block.diagnostics.production_optimized_pqs_product_path
    support_dense_coefficients = zeros(Float64, 5 * 5 * 5, 2)
    support_dense_coefficients[cubic.support_indices[1], 1] = 1.0
    support_dense_coefficients[cubic.support_indices[end], 2] = 1.0
    support_dense_unit = GaussletBases._nested_product_staged_generic_unit(
        :pqs_mixed_reference,
        support_dense_coefficients,
        99:100,
        (5, 5, 5),
    )
    support_dense_resolved = CCP._cartesian_resolved_contraction_payload(
        support_dense_unit;
        parent_dimension = 5 * 5 * 5,
    )
    pqs_mixed_dispatch = CCPM._metric_dispatch_plan_from_resolved_payloads(
        [cubic_pqs_resolved, support_dense_resolved],
    )
    @test pqs_mixed_dispatch.plan_supported
    @test pqs_mixed_dispatch.pqs_unit_count == 1
    @test pqs_mixed_dispatch.support_fallback_unit_count == 1
    @test pqs_mixed_dispatch.pqs_pqs_block_count == 1
    @test pqs_mixed_dispatch.pqs_support_block_count == 1
    @test :pqs_support_local_reference in [path.path for path in pqs_mixed_dispatch.block_paths]
    @test CCPM._metric_dispatch_block_path(
        (unsupported = false, block_role = :pqs),
        (unsupported = false, block_role = :product),
    ) == :unsupported_pqs_product_optimized
    pqs_mixed_block = CCPM._resolved_payload_low_order_metric_block(
        cubic_pqs_resolved,
        support_dense_resolved,
        cubic_metrics,
    )
    dense_coefficients = Matrix{Float64}(support_dense_unit.coefficient_matrix)
    expected_mixed_overlap =
        transpose(support_coefficients) *
        _pqs_cross_product_matrix(
            cubic.support_states,
            support_dense_unit.support_states,
            cubic_metrics.x.overlap,
            cubic_metrics.y.overlap,
            cubic_metrics.z.overlap,
        ) *
        dense_coefficients
    expected_mixed_position_x =
        transpose(support_coefficients) *
        _pqs_cross_product_matrix(
            cubic.support_states,
            support_dense_unit.support_states,
            cubic_metrics.x.position,
            cubic_metrics.y.overlap,
            cubic_metrics.z.overlap,
        ) *
        dense_coefficients
    @test pqs_mixed_block.path == :pqs_support_local_reference
    @test pqs_mixed_block.overlap ≈ expected_mixed_overlap atol = 1.0e-12 rtol = 1.0e-12
    @test pqs_mixed_block.position_x ≈ expected_mixed_position_x atol = 1.0e-12 rtol = 1.0e-12
    @test pqs_mixed_block.weights === nothing
    @test pqs_mixed_block.first_moments === nothing
    @test pqs_mixed_block.diagnostics.support_local_reference_used
    @test !pqs_mixed_block.diagnostics.pqs_self_overlap_invariant_applied
    sidecar_metric_reference =
        CCPM._projected_q_shell_sidecar_low_order_metric_reference(
            cubic_pqs_sidecar,
            cubic_metrics;
            mixed_payloads = (support_dense_resolved,),
        )
    @test sidecar_metric_reference.diagnostics.fixture_only
    @test sidecar_metric_reference.diagnostics.reference_scoped
    @test !sidecar_metric_reference.diagnostics.production_supported
    @test sidecar_metric_reference.diagnostics.support_local_reference_only
    @test !sidecar_metric_reference.diagnostics.fixed_block_sidecar_installed
    @test !sidecar_metric_reference.diagnostics.default_builder_consumes
    @test !sidecar_metric_reference.diagnostics.qw_consumes
    @test !sidecar_metric_reference.diagnostics.hamiltonian_consumes
    @test sidecar_metric_reference.diagnostics.payload_count == 1
    @test sidecar_metric_reference.diagnostics.self_block_count == 1
    @test sidecar_metric_reference.diagnostics.mixed_block_count == 1
    @test !sidecar_metric_reference.diagnostics.pqs_product_optimized_path_ready
    @test only(sidecar_metric_reference.resolved_payloads).ready_for_metric_execution
    @test only(sidecar_metric_reference.resolved_payloads).payload ===
          only(cubic_pqs_sidecar.payloads)
    sidecar_self_block = only(sidecar_metric_reference.self_blocks)
    sidecar_mixed_block = only(sidecar_metric_reference.mixed_blocks)
    @test sidecar_self_block.path == :pqs_pqs_low_order_reference
    @test sidecar_mixed_block.path == :pqs_support_local_reference
    @test sidecar_self_block.overlap == pqs_self_block.overlap
    @test sidecar_self_block.weights ≈ pqs_self_block.weights atol = 1.0e-12 rtol = 1.0e-12
    @test sidecar_self_block.position_x ≈ pqs_self_block.position_x atol = 1.0e-12 rtol = 1.0e-12
    @test sidecar_self_block.position_y ≈ pqs_self_block.position_y atol = 1.0e-12 rtol = 1.0e-12
    @test sidecar_self_block.position_z ≈ pqs_self_block.position_z atol = 1.0e-12 rtol = 1.0e-12
    @test sidecar_mixed_block.overlap ≈ pqs_mixed_block.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test sidecar_mixed_block.position_x ≈ pqs_mixed_block.position_x atol = 1.0e-12 rtol = 1.0e-12
    pqs_fixed_sidecar_block =
        CCP._cartesian_projected_q_shell_fixed_block_sidecar_fixture(
            bundle5.basis,
            cubic;
            gausslet_backend = :pgdg_localized_experimental,
        )
    @test pqs_fixed_sidecar_block isa GaussletBases._NestedFixedBlock3D
    @test pqs_fixed_sidecar_block.shell === cubic
    @test pqs_fixed_sidecar_block.gausslet_backend == :pgdg_localized_experimental
    @test size(pqs_fixed_sidecar_block.coefficient_matrix) == (5 * 5 * 5, 98)
    @test pqs_fixed_sidecar_block.support_indices == cubic.support_indices
    installed_pqs_sidecar =
        CCP._nested_projected_q_shell_sidecar_fixture(pqs_fixed_sidecar_block)
    @test installed_pqs_sidecar === pqs_fixed_sidecar_block.staged_by_center_sidecar[]
    @test installed_pqs_sidecar.dims == (5, 5, 5)
    @test installed_pqs_sidecar.provenance.source ==
          :projected_q_shell_fixed_block_sidecar_fixture
    @test installed_pqs_sidecar.provenance.installed_in_fixed_block
    @test installed_pqs_sidecar.provenance.fixed_block_sidecar_slot ==
          :staged_by_center_sidecar_fixture_only
    @test installed_pqs_sidecar.diagnostics.source ==
          :projected_q_shell_fixed_block_sidecar_fixture
    @test installed_pqs_sidecar.diagnostics.fixture_only
    @test installed_pqs_sidecar.diagnostics.fixed_block_sidecar_installed
    @test !installed_pqs_sidecar.diagnostics.by_center_consumes
    @test !installed_pqs_sidecar.diagnostics.default_builder_consumes
    @test !installed_pqs_sidecar.diagnostics.qw_consumes
    @test !installed_pqs_sidecar.diagnostics.hamiltonian_consumes
    @test !installed_pqs_sidecar.diagnostics.production_supported
    @test installed_pqs_sidecar.diagnostics.metric_capability ==
          :pqs_low_order_support_local_reference
    @test_throws ArgumentError GaussletBases._nested_staged_by_center_sidecar(
        pqs_fixed_sidecar_block,
    )
    @test GaussletBases._nested_by_center_sidecar_path(pqs_fixed_sidecar_block) ==
          :unknown_staged_sidecar
    installed_pqs_resolved = CCP._cartesian_resolved_contraction_payloads(
        installed_pqs_sidecar,
    )
    @test length(installed_pqs_resolved) == 1
    @test only(installed_pqs_resolved).ready_for_metric_execution
    @test only(installed_pqs_resolved).diagnostics.fixture_only
    @test !only(installed_pqs_resolved).diagnostics.production_supported
    @test only(installed_pqs_resolved).diagnostics.fixed_block_sidecar_installed
    @test only(installed_pqs_resolved).diagnostics.block_role == :pqs
    installed_metric_reference =
        CCPM._projected_q_shell_sidecar_low_order_metric_reference(
            installed_pqs_sidecar,
            cubic_metrics;
            mixed_payloads = (support_dense_resolved,),
        )
    @test installed_metric_reference.diagnostics.fixed_block_sidecar_installed
    @test !installed_metric_reference.diagnostics.default_builder_consumes
    @test !installed_metric_reference.diagnostics.qw_consumes
    @test !installed_metric_reference.diagnostics.hamiltonian_consumes
    installed_self_block = only(installed_metric_reference.self_blocks)
    installed_mixed_block = only(installed_metric_reference.mixed_blocks)
    @test installed_self_block.overlap == pqs_self_block.overlap
    @test installed_self_block.weights ≈ pqs_self_block.weights atol = 1.0e-12 rtol = 1.0e-12
    @test installed_self_block.position_x ≈ pqs_self_block.position_x atol = 1.0e-12 rtol = 1.0e-12
    @test installed_self_block.position_y ≈ pqs_self_block.position_y atol = 1.0e-12 rtol = 1.0e-12
    @test installed_self_block.position_z ≈ pqs_self_block.position_z atol = 1.0e-12 rtol = 1.0e-12
    @test installed_mixed_block.overlap ≈ pqs_mixed_block.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test installed_mixed_block.position_x ≈ pqs_mixed_block.position_x atol = 1.0e-12 rtol = 1.0e-12
    installed_source_transform =
        only(CCP._cartesian_raw_product_source_retained_transform(pqs_fixed_sidecar_block))
    pqs_raw_source = installed_source_transform.raw_source
    pqs_retained_transform = installed_source_transform.retained_transform
    @test pqs_raw_source.source_dimension == 5 * 5 * 5
    @test pqs_raw_source.local_box == cubic_current
    @test pqs_raw_source.axis_intervals == cubic_current
    @test pqs_raw_source.raw_source_weight_role == :raw_source_positive
    @test pqs_raw_source.support_indices === nothing
    @test pqs_raw_source.support_summary.support_scope == :full_local_product_box
    @test pqs_raw_source.support_summary.raw_product_support_count == 5 * 5 * 5
    @test pqs_raw_source.support_summary.boundary_support_count == 98
    @test pqs_raw_source.diagnostics.raw_source_contract == :full_local_product_block
    @test pqs_raw_source.diagnostics.retained_transform_expected_stage ==
          :shell_realized_support_local_fixture
    @test !pqs_raw_source.diagnostics.raw_product_box_operator_contract
    @test pqs_raw_source.diagnostics.shell_projection_lowdin_realization
    @test !pqs_raw_source.diagnostics.raw_box_boundary_selector_transform
    @test !pqs_raw_source.diagnostics.operator_matrices_built
    @test pqs_retained_transform.source_id == pqs_raw_source.source_id
    @test pqs_retained_transform.retained_dimension == 98
    @test pqs_retained_transform.transform_kind == :boundary_projection_lowdin
    @test pqs_retained_transform.transform_matrix === nothing
    @test pqs_retained_transform.transform_stages == (
        :raw_product_modes,
        :raw_boundary_projection,
        :full_rank_symmetric_lowdin_cleanup,
        :retained_columns,
    )
    @test pqs_retained_transform.retained_column_weight_role == :debug_reference_only
    @test pqs_retained_transform.diagnostics.transform_contract ==
          :factored_raw_product_to_retained
    @test pqs_retained_transform.diagnostics.representation_stage ==
          :shell_realized_support_local_fixture
    @test pqs_retained_transform.diagnostics.support_local_fixture_adapter
    @test !pqs_retained_transform.diagnostics.raw_product_box_operator_contract
    @test pqs_retained_transform.diagnostics.shell_projection_lowdin_realization
    @test !pqs_retained_transform.diagnostics.raw_box_boundary_selector_transform
    @test pqs_retained_transform.diagnostics.transform_matrix_scope ==
          :factored_raw_product_to_retained
    @test !pqs_retained_transform.diagnostics.full_raw_to_retained_matrix_materialized
    @test !pqs_retained_transform.diagnostics.full_raw_to_retained_matrix_available
    @test pqs_retained_transform.diagnostics.cleanup_stage_matrix_scope ==
          :boundary_mode_to_retained
    @test pqs_retained_transform.diagnostics.boundary_projection_stage ==
          :raw_product_modes_to_boundary_rows
    @test pqs_retained_transform.diagnostics.cleanup_stage ==
          :full_rank_symmetric_lowdin_boundary_cleanup
    @test pqs_retained_transform.diagnostics.raw_source_dimension == 5 * 5 * 5
    @test pqs_retained_transform.diagnostics.boundary_support_indices ==
          cubic.support_indices
    @test pqs_retained_transform.diagnostics.boundary_support_count == 98
    @test pqs_retained_transform.diagnostics.boundary_mode_indices ==
          cubic_descriptor.boundary_mode_indices
    @test pqs_retained_transform.diagnostics.boundary_column_indices ==
          cubic_descriptor.boundary_column_indices
    @test pqs_retained_transform.diagnostics.boundary_mode_count == 98
    @test pqs_retained_transform.diagnostics.retained_count == 98
    @test pqs_retained_transform.provenance.cleanup_stage_matrix ==
          cubic_pqs_payload.cleanup_transform
    @test pqs_retained_transform.diagnostics.cleanup_method ==
          :projected_boundary_symmetric_lowdin
    @test pqs_retained_transform.diagnostics.retained_weight_positive_checked == false
    @test pqs_retained_transform.diagnostics.ida_weight_division_allowed == false
    raw_box_pqs_transform = CCP._cartesian_raw_box_pqs_retained_rule_transform(
        cubic_descriptor;
        source_id = pqs_raw_source.source_id,
    )
    @test raw_box_pqs_transform.source_id == pqs_raw_source.source_id
    @test raw_box_pqs_transform.retained_dimension == 98
    @test raw_box_pqs_transform.transform_kind ==
          :boundary_comx_product_mode_selection
    @test raw_box_pqs_transform.transform_matrix === nothing
    @test raw_box_pqs_transform.transform_stages == (
        :raw_product_box_modes,
        :boundary_comx_product_mode_selection,
        :retained_columns,
    )
    @test raw_box_pqs_transform.retained_column_weight_role ==
          :not_positive_quadrature_weights
    @test raw_box_pqs_transform.diagnostics.representation_stage ==
          :mode_selected_raw_product_box
    @test raw_box_pqs_transform.diagnostics.raw_product_box_operator_contract
    @test raw_box_pqs_transform.diagnostics.raw_box_boundary_selector_transform
    @test !raw_box_pqs_transform.diagnostics.shell_projection_used
    @test !raw_box_pqs_transform.diagnostics.lowdin_cleanup_used
    @test !raw_box_pqs_transform.diagnostics.shell_projection_lowdin_realization
    @test !raw_box_pqs_transform.diagnostics.support_coefficient_matrix_used
    @test !raw_box_pqs_transform.diagnostics.support_local_fixture_adapter
    @test !raw_box_pqs_transform.diagnostics.ida_weight_division_allowed
    pqs_pair_packet = CCP._cartesian_raw_product_source_pair_operator_packet(
        pqs_raw_source,
        pqs_raw_source;
        operator_kind = :low_order_metric,
        supported_terms = (:overlap, :weights, :position_x, :position_y, :position_z),
    )
    @test pqs_pair_packet.left_source_id == pqs_raw_source.source_id
    @test pqs_pair_packet.right_source_id == pqs_raw_source.source_id
    @test pqs_pair_packet.operator_kind == :low_order_metric
    @test pqs_pair_packet.supported_terms ==
          (:overlap, :weights, :position_x, :position_y, :position_z)
    @test pqs_pair_packet.symmetry_status == :symmetric_upper_triangle_placeholder
    @test pqs_pair_packet.backend == :metadata_only
    @test pqs_pair_packet.operator_matrices === nothing
    @test pqs_pair_packet.diagnostics.placeholder_only
    @test !pqs_pair_packet.diagnostics.operator_matrices_built
    @test pqs_pair_packet.diagnostics.ids_terms_provenance_only
    @test !pqs_pair_packet.diagnostics.all_pairs_inventory_built
    @test !pqs_pair_packet.diagnostics.left_right_sources_embedded
    @test !pqs_pair_packet.diagnostics.left_right_retained_transforms_embedded
    @test pqs_pair_packet.diagnostics.future_inventory_must_resolve_sources_and_transforms
    @test !pqs_pair_packet.diagnostics.raw_operator_block_ready
    pqs_pair_plan = CCP._cartesian_raw_product_source_pair_plan(
        pqs_fixed_sidecar_block;
        operator_kind = :low_order_metric,
        supported_terms = (:overlap, :weights, :position_x, :position_y, :position_z),
    )
    @test pqs_pair_plan.operator_kind == :low_order_metric
    @test pqs_pair_plan.supported_terms ==
          (:overlap, :weights, :position_x, :position_y, :position_z)
    @test pqs_pair_plan.symmetry_status == :symmetric_upper_triangle_placeholder
    @test pqs_pair_plan.source_ids == [pqs_raw_source.source_id]
    @test length(pqs_pair_plan.raw_sources) == 1
    @test length(pqs_pair_plan.retained_transforms) == 1
    @test length(pqs_pair_plan.pair_packets) == 1
    @test pqs_pair_plan.pair_keys == [(pqs_raw_source.source_id, pqs_raw_source.source_id)]
    @test pqs_pair_plan.raw_sources[pqs_raw_source.source_id].source_id ==
          pqs_raw_source.source_id
    @test pqs_pair_plan.retained_transforms[pqs_raw_source.source_id].source_id ==
          pqs_raw_source.source_id
    @test pqs_pair_plan.retained_transforms[pqs_raw_source.source_id].transform_stages ==
          pqs_retained_transform.transform_stages
    @test pqs_pair_plan.retained_transforms[pqs_raw_source.source_id].transform_matrix ===
          nothing
    @test only(pqs_pair_plan.pair_packets).left_source_id == pqs_raw_source.source_id
    @test only(pqs_pair_plan.pair_packets).right_source_id == pqs_raw_source.source_id
    @test only(pqs_pair_plan.pair_packets).operator_matrices === nothing
    @test pqs_pair_plan.diagnostics.source == :projected_q_shell_fixed_block_pair_plan
    @test pqs_pair_plan.diagnostics.source_count == 1
    @test pqs_pair_plan.diagnostics.retained_transform_count == 1
    @test pqs_pair_plan.diagnostics.pair_count == 1
    @test pqs_pair_plan.diagnostics.expected_upper_triangle_pair_count == 1
    @test pqs_pair_plan.diagnostics.upper_triangle_only
    @test pqs_pair_plan.diagnostics.all_pairs_resolve_sources
    @test pqs_pair_plan.diagnostics.all_pairs_resolve_retained_transforms
    @test pqs_pair_plan.diagnostics.pair_packets_placeholder_only
    @test !pqs_pair_plan.diagnostics.raw_operator_matrices_built
    @test !pqs_pair_plan.diagnostics.retained_operator_blocks_built
    @test !pqs_pair_plan.diagnostics.metric_execution_changed
    @test !pqs_pair_plan.diagnostics.qwhamiltonian_consumes
    @test !pqs_pair_plan.diagnostics.public_default_consumes
    @test !pqs_pair_plan.diagnostics.backend_policy_changed
    @test !pqs_pair_plan.diagnostics.quadrature_policy_changed
    @test !pqs_pair_plan.diagnostics.cr2_science_status_changed
    pqs_resolved_pair = CCP._cartesian_resolve_raw_product_source_pair(pqs_pair_plan, 1)
    @test pqs_resolved_pair.pair_key ==
          (pqs_raw_source.source_id, pqs_raw_source.source_id)
    @test pqs_resolved_pair.left_raw_source.source_dimension == 5 * 5 * 5
    @test pqs_resolved_pair.right_raw_source.source_dimension == 5 * 5 * 5
    @test pqs_resolved_pair.left_retained_transform.retained_dimension == 98
    @test pqs_resolved_pair.right_retained_transform.retained_dimension == 98
    @test pqs_resolved_pair.left_retained_transform.transform_kind ==
          :boundary_projection_lowdin
    @test pqs_resolved_pair.right_retained_transform.transform_kind ==
          :boundary_projection_lowdin
    @test pqs_resolved_pair.diagnostics.upper_triangular
    @test pqs_resolved_pair.diagnostics.placeholder_only
    @test !pqs_resolved_pair.diagnostics.operator_matrices_built
    @test pqs_resolved_pair.diagnostics.left_raw_source_weight_role ==
          :raw_source_positive
    @test pqs_resolved_pair.diagnostics.right_raw_source_weight_role ==
          :raw_source_positive
    @test pqs_resolved_pair.diagnostics.left_retained_column_weight_role ==
          :debug_reference_only
    @test pqs_resolved_pair.diagnostics.right_retained_column_weight_role ==
          :debug_reference_only
    @test pqs_resolved_pair.diagnostics.raw_weight_roles_explicit
    @test pqs_resolved_pair.diagnostics.retained_weight_roles_explicit
    @test !pqs_resolved_pair.diagnostics.retained_ida_weight_division_allowed
    @test pqs_resolved_pair.diagnostics.pqs_factored_transform_present
    pqs_raw_overlap_packet = CCP._cartesian_raw_low_order_operator_packet(
        pqs_resolved_pair;
        term = :overlap,
    )
    @test pqs_raw_overlap_packet.left_source_id == pqs_raw_source.source_id
    @test pqs_raw_overlap_packet.right_source_id == pqs_raw_source.source_id
    @test pqs_raw_overlap_packet.operator_kind == :low_order_metric
    @test pqs_raw_overlap_packet.term == :overlap
    @test pqs_raw_overlap_packet.source_dimensions == (5 * 5 * 5, 5 * 5 * 5)
    @test pqs_raw_overlap_packet.symmetry_status == :symmetric_upper_triangle_placeholder
    @test pqs_raw_overlap_packet.backend == :private_raw_product_reference
    @test pqs_raw_overlap_packet.raw_operator_matrix ==
          Matrix{Float64}(I, 5 * 5 * 5, 5 * 5 * 5)
    @test pqs_raw_overlap_packet.diagnostics.raw_reference ==
          :orthonormal_raw_product_mode_overlap_identity
    @test pqs_raw_overlap_packet.diagnostics.raw_reference_error == 0.0
    @test pqs_raw_overlap_packet.diagnostics.raw_operator_matrix_built
    @test !pqs_raw_overlap_packet.diagnostics.retained_operator_block_built
    @test !pqs_raw_overlap_packet.diagnostics.retained_transform_applied
    @test !pqs_raw_overlap_packet.diagnostics.all_pair_matrices_built
    @test !pqs_raw_overlap_packet.diagnostics.metric_execution_changed
    @test !pqs_raw_overlap_packet.diagnostics.qwhamiltonian_consumes
    @test !pqs_raw_overlap_packet.diagnostics.public_default_consumes
    @test !pqs_raw_overlap_packet.diagnostics.backend_policy_changed
    @test !pqs_raw_overlap_packet.diagnostics.quadrature_policy_changed
    @test !pqs_raw_overlap_packet.diagnostics.cr2_science_status_changed
    @test only(pqs_pair_plan.pair_packets).operator_matrices === nothing
    pqs_plan_audit = CCP._cartesian_raw_product_source_pair_plan_audit(pqs_pair_plan)
    @test length(pqs_plan_audit.resolved_pairs) == 1
    @test pqs_plan_audit.diagnostics.every_pair_resolves_raw_sources
    @test pqs_plan_audit.diagnostics.every_pair_resolves_retained_transforms
    @test pqs_plan_audit.diagnostics.every_pair_upper_triangular
    @test pqs_plan_audit.diagnostics.every_pair_placeholder_only
    @test pqs_plan_audit.diagnostics.raw_weight_roles_explicit
    @test pqs_plan_audit.diagnostics.retained_weight_roles_explicit
    @test !pqs_plan_audit.diagnostics.retained_ida_weight_division_allowed
    @test !pqs_plan_audit.diagnostics.raw_operator_matrices_built
    @test !pqs_plan_audit.diagnostics.retained_operator_blocks_built
    @test !pqs_plan_audit.diagnostics.metric_execution_changed
    @test !pqs_plan_audit.diagnostics.qwhamiltonian_consumes
    @test !pqs_plan_audit.diagnostics.public_default_consumes
    @test !pqs_plan_audit.diagnostics.backend_policy_changed
    @test !pqs_plan_audit.diagnostics.quadrature_policy_changed
    @test !pqs_plan_audit.diagnostics.cr2_science_status_changed
    identity_axis = Matrix{Float64}(I, 2, 2)
    product_unit = GaussletBases._CartesianNestedProductStagedByCenterUnit3D(
        :identity_product_slab,
        :product_doside,
        1:4,
        collect(1:4),
        NTuple{3,Int}[(1, 1, 1), (1, 2, 1), (2, 1, 1), (2, 2, 1)],
        Matrix{Float64}(I, 4, 4),
        (
            GaussletBases._nested_product_staged_active_axis(1:2, identity_axis),
            GaussletBases._nested_product_staged_active_axis(1:2, identity_axis),
            GaussletBases._nested_product_staged_fixed_axis(1),
        ),
        GaussletBases._nested_product_axis_function_indices(3, 1, 2, 2, 2),
        (source = :raw_product_source_test_fixture,),
        (support_count = 4, retained_count = 4),
    )
    cubic_pqs_product_source_box_plan = CCPM._pqs_raw_product_box_plan(
        cubic_descriptor,
        cubic_metrics;
        shared_raw_product_box_plan = cubic_shared_raw_product_box_plan,
    )
    _check_pqs_product_source_box_mixed_block(
        CCPM,
        cubic_descriptor,
        cubic_pqs_product_source_box_plan,
        product_unit,
        cubic_metrics;
        expected_source_mode_dims = (5, 5, 5),
        expected_retained_count = 98,
        shared_raw_product_box_plan = cubic_shared_raw_product_box_plan,
    )
    _check_pqs_pqs_product_source_box_shadow_blocks(
        CCPM,
        shifted_left_descriptor,
        shifted_right_descriptor,
        shifted_left_pqs_pqs_source_box_plan,
        shifted_right_pqs_pqs_source_box_plan,
        product_unit,
        shifted_metrics;
        expected_source_mode_dims = (5, 5, 5),
        expected_retained_count = 98,
    )
    route_bundles =
        GaussletBases._CartesianNestedAxisBundles3D(bundle5, bundle5, bundle7)
    route_left_current = (1:5, 1:5, 1:5)
    route_left_inner = (2:4, 2:4, 2:4)
    route_right_current = (1:5, 1:5, 3:7)
    route_right_inner = (2:4, 2:4, 4:6)
    route_left = GaussletBases._nested_projected_q_shell_layer(
        route_bundles,
        route_left_current,
        route_left_inner;
        bond_axis = :z,
        q = 5,
        L = 5,
        term_coefficients,
    )
    route_right = GaussletBases._nested_projected_q_shell_layer(
        route_bundles,
        route_right_current,
        route_right_inner;
        bond_axis = :z,
        q = 5,
        L = 5,
        term_coefficients,
    )
    route_left_descriptor =
        GaussletBases._nested_projected_q_shell_staged_unit_descriptor(route_left)
    route_right_descriptor =
        GaussletBases._nested_projected_q_shell_staged_unit_descriptor(route_right)
    route_left_shared_raw_product_box_plan =
        GaussletBases._cartesian_raw_product_box_plan(
            route_bundles,
            route_left_descriptor.axis_intervals,
            (5, 5, 5);
            enforce_symmetric_odd = false,
        )
    route_right_shared_raw_product_box_plan =
        GaussletBases._cartesian_raw_product_box_plan(
            route_bundles,
            route_right_descriptor.axis_intervals,
            (5, 5, 5);
            enforce_symmetric_odd = false,
        )
    route_metrics = _pqs_axis_metrics(route_bundles)
    route_left_pqs_plan = CCPM._pqs_raw_product_box_plan(
        route_left_descriptor,
        route_left_shared_raw_product_box_plan,
        route_metrics,
    )
    route_right_pqs_plan = CCPM._pqs_raw_product_box_plan(
        route_right_descriptor,
        route_right_shared_raw_product_box_plan,
        route_metrics,
    )
    route_dims = (5, 5, 7)
    route_product_states = NTuple{3,Int}[
        (ix, iy, 4) for ix in 1:5 for iy in 1:5
    ]
    route_product_indices = [
        GaussletBases._cartesian_flat_index(state..., route_dims) for
        state in route_product_states
    ]
    route_identity_axis = Matrix{Float64}(I, 5, 5)
    route_product_axes = (
        GaussletBases._nested_product_staged_active_axis(
            1:5,
            route_identity_axis,
        ),
        GaussletBases._nested_product_staged_active_axis(
            1:5,
            route_identity_axis,
        ),
        GaussletBases._nested_product_staged_fixed_axis(4),
    )
    route_product_axis_indices =
        GaussletBases._nested_product_axis_function_indices(3, 1, 5, 2, 5)
    route_product_unit =
        GaussletBases._CartesianNestedProductStagedByCenterUnit3D(
            :middle_body_product_slab,
            :product_doside,
            1:25,
            route_product_indices,
            route_product_states,
            Matrix{Float64}(I, 25, 25),
            route_product_axes,
            route_product_axis_indices,
            (source = :route_shaped_safe_term_consumer_test_fixture,),
            (support_count = 25, retained_count = 25),
        )
    route_units = CCPM._pqs_pqs_product_safe_term_route_descriptor(
        route_left_pqs_plan,
        route_right_pqs_plan,
        route_product_unit;
        route_name = :q5_L5_slab5_test_route,
        parent_dims = route_dims,
        bond_axis = :z,
        metadata = (
            pqs_left_box = route_left_current,
            pqs_right_box = route_right_current,
            product_slab_fixed_index = 4,
            pqs_source_mode_dims = (5, 5, 5),
        ),
        provenance = (source = :route_shaped_safe_term_consumer_test_fixture,),
    )
    produced_route = CCPM._pqs_pqs_product_raw_box_route_producer(
        route_bundles,
        route_left_current,
        route_right_current,
        (1:5, 1:5, 4:4),
        route_metrics;
        source_mode_dims = (5, 5, 5),
        route_name = :q5_L5_slab5_test_route,
        parent_dims = route_dims,
        bond_axis = :z,
        metadata = (
            pqs_left_box = route_left_current,
            pqs_right_box = route_right_current,
            product_slab_fixed_index = 4,
            pqs_source_mode_dims = (5, 5, 5),
        ),
        provenance = (source = :route_shaped_safe_term_consumer_test_fixture,),
    )
    @test produced_route.object_kind == :pqs_pqs_product_raw_box_route_producer
    @test produced_route.status == :private_shadow_only
    @test produced_route.descriptor.object_kind ==
          :pqs_pqs_product_safe_term_route_descriptor
    @test produced_route.descriptor.expected_ranges == route_units.expected_ranges
    @test produced_route.descriptor.retained_dimension == route_units.retained_dimension
    @test produced_route.descriptor.expected_pair_count == route_units.expected_pair_count
    @test produced_route.descriptor.supported_terms == route_units.supported_terms
    @test map(summary -> summary.source_dimensions,
              produced_route.descriptor.unit_summaries) ==
          map(summary -> summary.source_dimensions, route_units.unit_summaries)
    @test map(summary -> summary.retained_count,
              produced_route.descriptor.unit_summaries) ==
          map(summary -> summary.retained_count, route_units.unit_summaries)
    @test produced_route.all_pairs_inventory.diagnostics.every_pair_uses_source_box_algorithmic_policy
    @test produced_route.all_pairs_inventory.diagnostics.source_box_algorithmic_pair_count == 6
    @test produced_route.diagnostics.raw_product_box_plan_built
    @test produced_route.diagnostics.retained_rule_built
    @test produced_route.diagnostics.route_descriptor_emitted
    @test produced_route.diagnostics.every_pair_uses_source_box_algorithmic_policy
    @test produced_route.diagnostics.source_box_algorithmic_pair_count == 6
    @test !produced_route.diagnostics.dense_raw_source_box_pair_matrix_materialized
    @test !produced_route.diagnostics.dense_raw_source_box_pair_matrix_materialized_by_producer
    @test produced_route.diagnostics.dense_raw_source_box_pair_matrices_validation_only
    @test !produced_route.diagnostics.shell_projection_used
    @test !produced_route.diagnostics.lowdin_cleanup_used
    @test !produced_route.diagnostics.support_local_pqs_oracle_used
    @test !produced_route.diagnostics.support_coefficient_matrix_used
    @test !produced_route.diagnostics.retained_pqs_weights_used
    @test produced_route.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !produced_route.diagnostics.ida_weight_division_allowed
    @test !produced_route.diagnostics.packet_adoption
    @test !produced_route.diagnostics.qwhamiltonian_consumes
    route_geometry_facts =
        CCPM._pqs_pqs_product_raw_box_homonuclear_geometry_facts(
            parent_dims = route_dims,
            bond_axis = :z,
            q = 5,
            L = 5,
            left_start = (1, 1, 1),
            right_shift = (0, 0, 2),
            product_slab_fixed_index = 4,
            route_name = :q5_L5_slab5_test_route,
            metadata = (
                pqs_left_box = route_left_current,
                pqs_right_box = route_right_current,
                product_slab_fixed_index = 4,
                pqs_source_mode_dims = (5, 5, 5),
            ),
            provenance = (source = :route_shaped_safe_term_consumer_test_fixture,),
        )
    @test route_geometry_facts.object_kind ==
          :pqs_pqs_product_raw_box_homonuclear_geometry_facts
    @test route_geometry_facts.status == :private_fixture_geometry_producer
    @test route_geometry_facts.left_source_box == route_left_current
    @test route_geometry_facts.right_source_box == route_right_current
    @test route_geometry_facts.product_source_box == (1:5, 1:5, 4:4)
    @test route_geometry_facts.source_mode_dims == (5, 5, 5)
    @test route_geometry_facts.diagnostics.emits_explicit_source_box_inputs
    @test route_geometry_facts.diagnostics.source_boxes_inside_parent_dims
    @test route_geometry_facts.diagnostics.raw_box_route_producer_consumes
    @test !route_geometry_facts.diagnostics.shell_projection_used
    @test !route_geometry_facts.diagnostics.lowdin_cleanup_used
    @test !route_geometry_facts.diagnostics.support_local_pqs_oracle_used
    @test !route_geometry_facts.diagnostics.support_coefficient_matrix_used
    @test !route_geometry_facts.diagnostics.retained_pqs_weights_used
    @test !route_geometry_facts.diagnostics.ida_weight_division_allowed
    @test !route_geometry_facts.diagnostics.packet_adoption
    geometry_route = CCPM._pqs_pqs_product_raw_box_route_from_geometry_facts(
        route_bundles,
        route_geometry_facts,
        route_metrics,
    )
    @test geometry_route.object_kind ==
          :pqs_pqs_product_raw_box_geometry_route_producer
    @test geometry_route.geometry_facts === route_geometry_facts
    @test geometry_route.produced_route.object_kind ==
          :pqs_pqs_product_raw_box_route_producer
    @test geometry_route.descriptor.expected_ranges ==
          produced_route.descriptor.expected_ranges
    @test geometry_route.descriptor.retained_dimension ==
          produced_route.descriptor.retained_dimension
    @test geometry_route.descriptor.expected_pair_count ==
          produced_route.descriptor.expected_pair_count
    @test geometry_route.descriptor.supported_terms ==
          produced_route.descriptor.supported_terms
    @test map(summary -> summary.source_dimensions,
              geometry_route.descriptor.unit_summaries) ==
          map(summary -> summary.source_dimensions,
              produced_route.descriptor.unit_summaries)
    @test map(summary -> summary.retained_count,
              geometry_route.descriptor.unit_summaries) ==
          map(summary -> summary.retained_count,
              produced_route.descriptor.unit_summaries)
    @test geometry_route.all_pairs_inventory.diagnostics.source_box_algorithmic_pair_count == 6
    @test geometry_route.diagnostics.geometry_facts_consumed
    @test geometry_route.diagnostics.private_fixture_geometry_producer
    @test !geometry_route.diagnostics.shell_projection_used
    @test !geometry_route.diagnostics.lowdin_cleanup_used
    @test !geometry_route.diagnostics.support_local_pqs_oracle_used
    @test !geometry_route.diagnostics.support_coefficient_matrix_used
    @test !geometry_route.diagnostics.retained_pqs_weights_used
    @test !geometry_route.diagnostics.ida_weight_division_allowed
    @test !geometry_route.diagnostics.packet_adoption
    geometry_route_consumer =
        CCPM._pqs_pqs_product_route_shaped_safe_term_consumer(
            geometry_route.descriptor,
            route_metrics,
        )
    produced_route_consumer_for_geometry =
        CCPM._pqs_pqs_product_route_shaped_safe_term_consumer(
            produced_route.descriptor,
            route_metrics,
        )
    max_geometry_route_error = maximum(
        norm(
            geometry_route_consumer.blocks[term] -
            produced_route_consumer_for_geometry.blocks[term],
            Inf,
        ) for term in geometry_route_consumer.terms
    )
    @test max_geometry_route_error < 1.0e-10
    @test_throws ArgumentError CCPM._pqs_pqs_product_raw_box_homonuclear_geometry_facts(
        parent_dims = route_dims,
        bond_axis = :z,
        q = 5,
        L = 5,
        left_start = (1, 1, 1),
        right_shift = (0, 0, 2),
        product_slab_fixed_index = 8,
    )
    @test_throws ArgumentError CCPM._pqs_pqs_product_raw_box_route_from_geometry_facts(
        route_bundles,
        merge(route_geometry_facts, (parent_dims = (5, 5, 8),)),
        route_metrics,
    )
    @test_throws ArgumentError CCPM._pqs_pqs_product_raw_box_homonuclear_geometry_facts(
        parent_dims = route_dims,
        bond_axis = :w,
        q = 5,
        L = 5,
        left_start = (1, 1, 1),
        right_shift = (0, 0, 2),
        product_slab_fixed_index = 4,
    )
    @test_throws ArgumentError CCPM._pqs_pqs_product_raw_box_homonuclear_geometry_facts(
        parent_dims = route_dims,
        bond_axis = :z,
        left_start = (1, 1, 1),
        right_shift = (0, 0, 2),
        product_slab_fixed_index = 4,
    )
    @test_throws ArgumentError CCPM._pqs_pqs_product_raw_box_homonuclear_geometry_facts(
        parent_dims = route_dims,
        bond_axis = :z,
        source_mode_dims = (5, 5),
        left_start = (1, 1, 1),
        right_shift = (0, 0, 2),
        product_slab_fixed_index = 4,
    )
    @test_throws ArgumentError CCPM._pqs_pqs_product_raw_box_homonuclear_geometry_facts(
        parent_dims = route_dims,
        bond_axis = :z,
        source_mode_dims = (5, 1, 5),
        left_start = (1, 1, 1),
        right_shift = (0, 0, 2),
        product_slab_fixed_index = 4,
    )
    @test_throws ArgumentError CCPM._pqs_pqs_product_raw_box_route_producer(
        route_bundles,
        (0:4, 1:5, 1:5),
        route_right_current,
        (1:5, 1:5, 4:4),
        route_metrics;
        source_mode_dims = (5, 5, 5),
        parent_dims = route_dims,
    )
    @test_throws ArgumentError CCPM._pqs_pqs_product_raw_box_route_producer(
        route_bundles,
        route_left_current,
        (1:5, 1:5, 4:8),
        (1:5, 1:5, 4:4),
        route_metrics;
        source_mode_dims = (5, 5, 5),
        parent_dims = route_dims,
    )
    @test_throws ArgumentError CCPM._pqs_pqs_product_raw_box_route_producer(
        route_bundles,
        route_left_current,
        route_right_current,
        (1:5, 1:5, 8:8),
        route_metrics;
        source_mode_dims = (5, 5, 5),
        parent_dims = route_dims,
    )
    @test_throws ArgumentError CCPM._pqs_pqs_product_raw_box_route_producer(
        route_bundles,
        route_left_current,
        route_right_current,
        (1:5, 1:5, 4:4),
        route_metrics;
        source_mode_dims = (1, 5, 5),
        parent_dims = route_dims,
    )
    @test_throws ArgumentError CCPM._pqs_pqs_product_raw_box_route_producer(
        route_bundles,
        route_left_current,
        route_right_current,
        (1:5, 1:5, 1:5),
        route_metrics;
        source_mode_dims = (5, 5, 5),
        parent_dims = route_dims,
    )
    @test_throws ArgumentError CCPM._pqs_pqs_product_raw_box_route_producer(
        route_bundles,
        route_left_current,
        route_right_current,
        (1:5, 4:4, 4:4),
        route_metrics;
        source_mode_dims = (5, 5, 5),
        parent_dims = route_dims,
    )
    try
        CCPM._pqs_pqs_product_raw_box_route_producer(
            route_bundles,
            (0:4, 1:5, 1:5),
            route_right_current,
            (1:5, 1:5, 4:4),
            route_metrics;
            source_mode_dims = (5, 5, 5),
            parent_dims = route_dims,
            supported_terms = (:weights,),
        )
        error("unsupported route-producer term was accepted")
    catch err
        @test err isa ArgumentError
        @test occursin("unsupported term", sprint(showerror, err))
    end
    _check_pqs_pqs_product_route_shaped_safe_term_consumer(
        CCPM,
        route_left_descriptor,
        route_right_descriptor,
        route_units,
        route_metrics,
    )
    _check_pqs_pqs_product_route_shaped_density_density_consumer(
        CCPM,
        route_units,
    )
    _check_pqs_pqs_product_raw_box_density_density_route_producer(
        CCPM,
        route_bundles,
        route_units,
        route_metrics,
        route_left_current,
        route_right_current,
        (1:5, 1:5, 4:4);
        source_mode_dims = (5, 5, 5),
        parent_dims = route_dims,
        bond_axis = :z,
        ida_term_coefficients = term_coefficients,
        ida_dense_parent_matrix = GaussletBases._qwrg_diatomic_interaction_matrix(
            route_bundles.bundle_x,
            route_bundles.bundle_y,
            route_bundles.bundle_z,
            expansion,
        ),
    )
    _check_pqs_pqs_product_source_box_component_route_smoke(
        CCPM,
        route_bundles,
        route_metrics,
        route_left_current,
        route_right_current,
        (1:5, 1:5, 4:4);
        source_mode_dims = (5, 5, 5),
        parent_dims = route_dims,
        bond_axis = :z,
        term_coefficients,
        ida_dense_parent_matrix = GaussletBases._qwrg_diatomic_interaction_matrix(
            route_bundles.bundle_x,
            route_bundles.bundle_y,
            route_bundles.bundle_z,
            expansion,
        ),
        nuclear_expansion = expansion,
    )
    produced_route_consumer =
        CCPM._pqs_pqs_product_route_shaped_safe_term_consumer(
            produced_route.descriptor,
            route_metrics,
        )
    hand_built_route_consumer =
        CCPM._pqs_pqs_product_route_shaped_safe_term_consumer(
            route_units,
            route_metrics,
        )
    max_produced_route_error = maximum(
        norm(
            produced_route_consumer.blocks[term] -
            hand_built_route_consumer.blocks[term],
            Inf,
        ) for term in produced_route_consumer.terms
    )
    @test max_produced_route_error < 1.0e-10
    @test produced_route_consumer.diagnostics.every_pair_uses_source_box_algorithmic_policy
    @test produced_route_consumer.diagnostics.source_box_algorithmic_pair_count == 6
    route_producer_sample_specs = (
        (
            name = :shifted_cubic_q5_L5,
            left_box = route_left_current,
            right_box = route_right_current,
            product_box = (1:5, 1:5, 4:4),
            source_mode_dims = (5, 5, 5),
        ),
        (
            name = :same_box_rectangular_q5_L7,
            left_box = (1:5, 1:5, 1:7),
            right_box = (1:5, 1:5, 1:7),
            product_box = (1:5, 1:5, 4:4),
            source_mode_dims = (5, 5, 7),
        ),
    )
    route_producer_sample_summaries = Any[]
    for sample in route_producer_sample_specs
        sample_producer_timed = @timed CCPM._pqs_pqs_product_raw_box_route_producer(
            route_bundles,
            sample.left_box,
            sample.right_box,
            sample.product_box,
            route_metrics;
            source_mode_dims = sample.source_mode_dims,
            route_name = sample.name,
            parent_dims = route_dims,
            bond_axis = :z,
            metadata = (
                sample_name = sample.name,
                sampled_route_producer_validation = true,
            ),
            provenance = (source = :sampled_route_producer_validation,),
        )
        sample_producer = sample_producer_timed.value
        sample_consumer_timed =
            @timed CCPM._pqs_pqs_product_route_shaped_safe_term_consumer(
                sample_producer.descriptor,
                route_metrics,
            )
        sample_consumer = sample_consumer_timed.value
        sample_shadow = CCPM._pqs_pqs_product_source_box_shadow_blocks(
            sample_producer.raw_pqs_plans.pqs_left,
            sample_producer.raw_pqs_plans.pqs_right,
            sample_producer.product_unit,
            route_metrics,
        )
        sample_error = maximum(
            norm(
                sample_consumer.blocks[term] - sample_shadow.blocks[term],
                Inf,
            ) for term in sample_consumer.terms
        )
        @test sample_error < 1.0e-10
        @test sample_producer.diagnostics.every_pair_uses_source_box_algorithmic_policy
        @test sample_consumer.diagnostics.every_pair_uses_source_box_algorithmic_policy
        @test sample_consumer.pair_count == 6
        @test !sample_producer.diagnostics.dense_raw_source_box_pair_matrix_materialized
        @test !sample_producer.diagnostics.shell_projection_used
        @test !sample_producer.diagnostics.lowdin_cleanup_used
        @test !sample_producer.diagnostics.support_local_pqs_oracle_used
        @test !sample_producer.diagnostics.support_coefficient_matrix_used
        @test !sample_producer.diagnostics.retained_pqs_weights_used
        @test !sample_producer.diagnostics.ida_weight_division_allowed
        @test !sample_consumer.diagnostics.packet_adoption
        @test !sample_consumer.diagnostics.fixed_block_routing
        @test !sample_consumer.diagnostics.qwhamiltonian_consumes
        sample_left_start = ntuple(axis -> first(sample.left_box[axis]), 3)
        sample_right_shift = ntuple(
            axis -> first(sample.right_box[axis]) - first(sample.left_box[axis]),
            3,
        )
        sample_geometry_facts =
            CCPM._pqs_pqs_product_raw_box_homonuclear_geometry_facts(
                parent_dims = route_dims,
                bond_axis = :z,
                source_mode_dims = sample.source_mode_dims,
                left_start = sample_left_start,
                right_shift = sample_right_shift,
                product_slab_fixed_index = first(sample.product_box[3]),
                route_name = sample.name,
                metadata = (
                    sample_name = sample.name,
                    sampled_geometry_producer_validation = true,
                ),
                provenance = (source = :sampled_geometry_producer_validation,),
            )
        @test sample_geometry_facts.left_source_box == sample.left_box
        @test sample_geometry_facts.right_source_box == sample.right_box
        @test sample_geometry_facts.product_source_box == sample.product_box
        @test sample_geometry_facts.source_mode_dims == sample.source_mode_dims
        @test !sample_geometry_facts.diagnostics.shell_projection_used
        @test !sample_geometry_facts.diagnostics.lowdin_cleanup_used
        @test !sample_geometry_facts.diagnostics.support_local_pqs_oracle_used
        @test !sample_geometry_facts.diagnostics.retained_pqs_weights_used
        @test !sample_geometry_facts.diagnostics.ida_weight_division_allowed
        sample_geometry_route =
            CCPM._pqs_pqs_product_raw_box_route_from_geometry_facts(
                route_bundles,
                sample_geometry_facts,
                route_metrics,
            )
        @test sample_geometry_route.descriptor.expected_ranges ==
              sample_producer.descriptor.expected_ranges
        @test sample_geometry_route.descriptor.retained_dimension ==
              sample_producer.descriptor.retained_dimension
        @test sample_geometry_route.descriptor.expected_pair_count ==
              sample_producer.descriptor.expected_pair_count
        @test sample_geometry_route.descriptor.supported_terms ==
              sample_producer.descriptor.supported_terms
        @test map(summary -> summary.source_dimensions,
                  sample_geometry_route.descriptor.unit_summaries) ==
              map(summary -> summary.source_dimensions,
                  sample_producer.descriptor.unit_summaries)
        @test map(summary -> summary.retained_count,
                  sample_geometry_route.descriptor.unit_summaries) ==
              map(summary -> summary.retained_count,
                  sample_producer.descriptor.unit_summaries)
        sample_geometry_consumer =
            CCPM._pqs_pqs_product_route_shaped_safe_term_consumer(
                sample_geometry_route.descriptor,
                route_metrics,
            )
        sample_geometry_error = maximum(
            norm(
                sample_geometry_consumer.blocks[term] -
                sample_consumer.blocks[term],
                Inf,
            ) for term in sample_geometry_consumer.terms
        )
        @test sample_geometry_error < 1.0e-10
        @test sample_geometry_route.diagnostics.geometry_facts_consumed
        @test !sample_geometry_route.diagnostics.shell_projection_used
        @test !sample_geometry_route.diagnostics.lowdin_cleanup_used
        @test !sample_geometry_route.diagnostics.support_local_pqs_oracle_used
        @test !sample_geometry_route.diagnostics.retained_pqs_weights_used
        @test !sample_geometry_route.diagnostics.ida_weight_division_allowed
        push!(
            route_producer_sample_summaries,
            (
                name = sample.name,
                source_mode_dims = sample.source_mode_dims,
                left_box_lengths = ntuple(axis -> length(sample.left_box[axis]), 3),
                right_box_lengths = ntuple(axis -> length(sample.right_box[axis]), 3),
                product_box_lengths = ntuple(axis -> length(sample.product_box[axis]), 3),
                retained_dimension = sample_consumer.retained_dimension,
                pair_count = sample_consumer.pair_count,
                producer_elapsed_seconds = Float64(sample_producer_timed.time),
                producer_allocated_bytes = Int(sample_producer_timed.bytes),
                consumer_elapsed_seconds = Float64(sample_consumer_timed.time),
                consumer_allocated_bytes = Int(sample_consumer_timed.bytes),
                dense_raw_source_box_pair_matrix_validation_only =
                    sample_consumer.diagnostics.dense_raw_source_box_pair_matrix_materialized_for_validation,
                max_consumer_shadow_error = sample_error,
                max_geometry_consumer_error = sample_geometry_error,
            ),
        )
    end
    @test length(route_producer_sample_summaries) == 2
    @test any(
        summary -> summary.source_mode_dims[3] != summary.source_mode_dims[1],
        route_producer_sample_summaries,
    )
    @test all(summary -> summary.pair_count == 6, route_producer_sample_summaries)
    @test all(
        summary -> summary.producer_allocated_bytes > 0 &&
                   summary.consumer_allocated_bytes > 0,
        route_producer_sample_summaries,
    )
    @test all(
        summary -> summary.max_consumer_shadow_error < 1.0e-10,
        route_producer_sample_summaries,
    )
    @test all(
        summary -> summary.max_geometry_consumer_error < 1.0e-10,
        route_producer_sample_summaries,
    )
    x_axis_route_bundles =
        GaussletBases._CartesianNestedAxisBundles3D(bundle7, bundle5, bundle5)
    x_axis_route_metrics = _pqs_axis_metrics(x_axis_route_bundles)
    x_axis_route_dims = (7, 5, 5)
    x_axis_left_box = (1:5, 1:5, 1:5)
    x_axis_right_box = (3:7, 1:5, 1:5)
    x_axis_product_box = (4:4, 1:5, 1:5)
    x_axis_explicit_route =
        CCPM._pqs_pqs_product_raw_box_route_producer(
            x_axis_route_bundles,
            x_axis_left_box,
            x_axis_right_box,
            x_axis_product_box,
            x_axis_route_metrics;
            source_mode_dims = (5, 5, 5),
            route_name = :x_axis_shifted_cubic_q5_L5,
            parent_dims = x_axis_route_dims,
            bond_axis = :x,
            metadata = (
                sample_name = :x_axis_shifted_cubic_q5_L5,
                axis_general_geometry_validation = true,
            ),
            provenance = (source = :axis_general_geometry_validation,),
        )
    x_axis_geometry_facts =
        CCPM._pqs_pqs_product_raw_box_homonuclear_geometry_facts(
            parent_dims = x_axis_route_dims,
            bond_axis = :x,
            q = 5,
            L = 5,
            left_start = (1, 1, 1),
            right_shift = (2, 0, 0),
            product_slab_fixed_index = 4,
            route_name = :x_axis_shifted_cubic_q5_L5,
            metadata = (
                sample_name = :x_axis_shifted_cubic_q5_L5,
                axis_general_geometry_validation = true,
            ),
            provenance = (source = :axis_general_geometry_validation,),
        )
    @test x_axis_geometry_facts.bond_axis == :x
    @test x_axis_geometry_facts.bond_axis_index == 1
    @test x_axis_geometry_facts.left_source_box == x_axis_left_box
    @test x_axis_geometry_facts.right_source_box == x_axis_right_box
    @test x_axis_geometry_facts.product_source_box == x_axis_product_box
    @test x_axis_geometry_facts.source_mode_dims == (5, 5, 5)
    @test x_axis_geometry_facts.diagnostics.product_slab_fixed_axis == 1
    @test x_axis_geometry_facts.diagnostics.source_boxes_inside_parent_dims
    @test !x_axis_geometry_facts.diagnostics.shell_projection_used
    @test !x_axis_geometry_facts.diagnostics.lowdin_cleanup_used
    @test !x_axis_geometry_facts.diagnostics.support_local_pqs_oracle_used
    @test !x_axis_geometry_facts.diagnostics.support_coefficient_matrix_used
    @test !x_axis_geometry_facts.diagnostics.retained_pqs_weights_used
    @test !x_axis_geometry_facts.diagnostics.ida_weight_division_allowed
    x_axis_geometry_route =
        CCPM._pqs_pqs_product_raw_box_route_from_geometry_facts(
            x_axis_route_bundles,
            x_axis_geometry_facts,
            x_axis_route_metrics,
        )
    @test x_axis_geometry_route.descriptor.expected_ranges ==
          x_axis_explicit_route.descriptor.expected_ranges
    @test x_axis_geometry_route.descriptor.retained_dimension ==
          x_axis_explicit_route.descriptor.retained_dimension
    @test x_axis_geometry_route.descriptor.expected_pair_count ==
          x_axis_explicit_route.descriptor.expected_pair_count
    @test x_axis_geometry_route.descriptor.supported_terms ==
          x_axis_explicit_route.descriptor.supported_terms
    @test map(summary -> summary.source_dimensions,
              x_axis_geometry_route.descriptor.unit_summaries) ==
          map(summary -> summary.source_dimensions,
              x_axis_explicit_route.descriptor.unit_summaries)
    @test map(summary -> summary.retained_count,
              x_axis_geometry_route.descriptor.unit_summaries) ==
          map(summary -> summary.retained_count,
              x_axis_explicit_route.descriptor.unit_summaries)
    @test x_axis_geometry_route.all_pairs_inventory.diagnostics.every_pair_uses_source_box_algorithmic_policy
    @test x_axis_geometry_route.all_pairs_inventory.diagnostics.source_box_algorithmic_pair_count == 6
    @test x_axis_geometry_route.diagnostics.geometry_facts_consumed
    @test !x_axis_geometry_route.diagnostics.shell_projection_used
    @test !x_axis_geometry_route.diagnostics.lowdin_cleanup_used
    @test !x_axis_geometry_route.diagnostics.support_local_pqs_oracle_used
    @test !x_axis_geometry_route.diagnostics.support_coefficient_matrix_used
    @test !x_axis_geometry_route.diagnostics.retained_pqs_weights_used
    @test !x_axis_geometry_route.diagnostics.ida_weight_division_allowed
    x_axis_explicit_consumer =
        CCPM._pqs_pqs_product_route_shaped_safe_term_consumer(
            x_axis_explicit_route.descriptor,
            x_axis_route_metrics,
        )
    x_axis_geometry_consumer =
        CCPM._pqs_pqs_product_route_shaped_safe_term_consumer(
            x_axis_geometry_route.descriptor,
            x_axis_route_metrics,
        )
    x_axis_geometry_error = maximum(
        norm(
            x_axis_geometry_consumer.blocks[term] -
            x_axis_explicit_consumer.blocks[term],
            Inf,
        ) for term in x_axis_geometry_consumer.terms
    )
    @test x_axis_geometry_error < 1.0e-10
    @test x_axis_geometry_consumer.pair_count == 6
    @test !x_axis_geometry_consumer.diagnostics.packet_adoption
    @test !x_axis_geometry_consumer.diagnostics.fixed_block_routing
    @test !x_axis_geometry_consumer.diagnostics.qwhamiltonian_consumes
    @test x_axis_geometry_consumer.diagnostics.dense_raw_source_box_pair_matrix_materialized_for_validation
    route_fact_diagnostic =
        CCPM._pqs_pqs_product_route_descriptor_diagnostic(
            route_units,
            route_metrics,
        )
    @test route_fact_diagnostic.status == :descriptor_available
    @test route_fact_diagnostic.descriptor === route_units
    @test route_fact_diagnostic.diagnostics.pqs_descriptor_count == 2
    @test route_fact_diagnostic.diagnostics.pqs_raw_plan_convertible_count == 2
    @test route_fact_diagnostic.diagnostics.product_doside_unit_count == 1
    @test route_fact_diagnostic.diagnostics.direct_or_support_body_piece_count == 0
    @test route_fact_diagnostic.diagnostics.descriptor_emitted
    @test !route_fact_diagnostic.diagnostics.packet_adoption
    @test !route_fact_diagnostic.diagnostics.fixed_block_construction_changed
    @test !route_fact_diagnostic.diagnostics.qwhamiltonian_changed
    @test !route_fact_diagnostic.diagnostics.shell_projection_used
    @test !route_fact_diagnostic.diagnostics.lowdin_cleanup_used
    @test !route_fact_diagnostic.diagnostics.support_local_pqs_oracle_used
    @test route_fact_diagnostic.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !route_fact_diagnostic.diagnostics.ida_weight_division_allowed
    @test !route_fact_diagnostic.diagnostics.direct_support_reinterpreted_as_product_doside
    route_fact_consumer =
        CCPM._pqs_pqs_product_route_shaped_safe_term_consumer(
            route_fact_diagnostic.descriptor,
            route_metrics,
        )
    @test route_fact_consumer.retained_dimension == route_units.retained_dimension
    @test route_fact_consumer.pair_count == route_units.expected_pair_count
    nonidentity_product_axes = (
        GaussletBases._nested_product_staged_active_axis(
            1:2,
            [
                inv(sqrt(2.0)) inv(sqrt(2.0))
                inv(sqrt(2.0)) -inv(sqrt(2.0))
            ],
        ),
        GaussletBases._nested_product_staged_active_axis(1:2, identity_axis),
        GaussletBases._nested_product_staged_fixed_axis(1),
    )
    nonidentity_product_axis_indices =
        GaussletBases._nested_product_axis_function_indices(3, 1, 2, 2, 2)
    nonidentity_product_coefficients =
        _product_staged_comparison_coefficient_matrix(
            NTuple{3,Int}[(1, 1, 1), (1, 2, 1), (2, 1, 1), (2, 2, 1)],
            nonidentity_product_axes,
            nonidentity_product_axis_indices,
        )
    @test !(nonidentity_product_coefficients ≈ Matrix{Float64}(I, 4, 4))
    nonidentity_axis_product_unit =
        GaussletBases._CartesianNestedProductStagedByCenterUnit3D(
            :nonidentity_axis_product_slab,
            :product_doside,
            1:4,
            collect(1:4),
            NTuple{3,Int}[(1, 1, 1), (1, 2, 1), (2, 1, 1), (2, 2, 1)],
            nonidentity_product_coefficients,
            nonidentity_product_axes,
            nonidentity_product_axis_indices,
            (source = :nonidentity_axis_source_box_test_fixture,),
            (support_count = 4, retained_count = 4),
        )
    for term in (:overlap, :position_x, :x2_x, :kinetic)
        nonidentity_block = CCPM._pqs_product_source_box_reference_block(
            cubic_pqs_product_source_box_plan,
            nonidentity_axis_product_unit,
            cubic_metrics;
            term,
        )
        nonidentity_expected = _pqs_product_source_box_explicit_reference(
            cubic_descriptor,
            nonidentity_axis_product_unit,
            cubic_metrics;
            term,
        )
        @test nonidentity_block.block ≈ nonidentity_expected atol = 1.0e-10 rtol = 1.0e-10
        @test !nonidentity_block.diagnostics.shell_projection_used
        @test !nonidentity_block.diagnostics.lowdin_cleanup_used
        @test nonidentity_block.diagnostics.shared_raw_product_box_plan_used
        @test !nonidentity_block.diagnostics.support_local_pqs_oracle_used
        @test !nonidentity_block.diagnostics.retained_pqs_weights_used
        @test !nonidentity_block.diagnostics.ida_weight_division_allowed
        @test !nonidentity_block.diagnostics.packet_adoption
        @test !nonidentity_block.diagnostics.qwhamiltonian_consumes
        @test !nonidentity_block.diagnostics.public_default_consumes
        @test !nonidentity_block.diagnostics.cr2_science_status_changed
    end
    product_source_transform = CCP._cartesian_raw_product_source_retained_transform(
        product_unit;
        source_id = :identity_product_slab_source,
        parent_dims = (2, 2, 1),
    )
    @test product_source_transform.raw_source.source_id == :identity_product_slab_source
    @test product_source_transform.raw_source.source_dimension == 4
    @test product_source_transform.raw_source.support_indices == collect(1:4)
    @test product_source_transform.raw_source.raw_source_weight_role == :raw_source_positive
    @test product_source_transform.raw_source.local_box == (1:2, 1:2, 1:1)
    @test product_source_transform.retained_transform.source_id ==
          :identity_product_slab_source
    @test product_source_transform.retained_transform.retained_dimension == 4
    @test product_source_transform.retained_transform.transform_kind ==
          :product_axis_transform
    @test product_source_transform.retained_transform.transform_stages == (
        :raw_product_source,
        :separable_axis_transform_metadata,
        :retained_columns,
    )
    @test product_source_transform.retained_transform.retained_column_weight_role ==
          :debug_reference_only
    @test !product_source_transform.retained_transform.diagnostics.ida_weight_division_allowed
    @test product_source_transform.retained_transform.diagnostics.full_raw_to_retained_matrix_materialized
    @test product_source_transform.retained_transform.diagnostics.fast_product_path_requires_separable_axis_transforms
    @test product_source_transform.retained_transform.diagnostics.separable_axis_transforms_available
    pqs_product_support_coefficients =
        Matrix{Float64}(cubic_pqs_payload.support_coefficient_matrix)
    function _pqs_product_expected_low_order_block(term::Symbol)
        axis_matrices = term == :overlap ? (
            cubic_metrics.x.overlap,
            cubic_metrics.y.overlap,
            cubic_metrics.z.overlap,
        ) : term == :position_x ? (
            cubic_metrics.x.position,
            cubic_metrics.y.overlap,
            cubic_metrics.z.overlap,
        ) : term == :position_y ? (
            cubic_metrics.x.overlap,
            cubic_metrics.y.position,
            cubic_metrics.z.overlap,
        ) : term == :position_z ? (
            cubic_metrics.x.overlap,
            cubic_metrics.y.overlap,
            cubic_metrics.z.position,
        ) : throw(ArgumentError("unsupported PQS/product test term"))
        return transpose(pqs_product_support_coefficients) *
               _pqs_cross_product_matrix(
                   cubic_pqs_payload.support_states,
                   product_unit.support_states,
                   axis_matrices[1],
                   axis_matrices[2],
                   axis_matrices[3],
               ) *
               Matrix{Float64}(product_unit.coefficient_matrix)
    end
    cubic_pgdg_x = GaussletBases._nested_axis_pgdg(cubic_bundles, :x)
    cubic_pgdg_y = GaussletBases._nested_axis_pgdg(cubic_bundles, :y)
    cubic_pgdg_z = GaussletBases._nested_axis_pgdg(cubic_bundles, :z)
    cubic_kinetic_axis_ops = (
        x = (overlap = cubic_pgdg_x.overlap, kinetic = cubic_pgdg_x.kinetic),
        y = (overlap = cubic_pgdg_y.overlap, kinetic = cubic_pgdg_y.kinetic),
        z = (overlap = cubic_pgdg_z.overlap, kinetic = cubic_pgdg_z.kinetic),
    )
    function _pqs_product_expected_kinetic_block()
        x_term = _pqs_cross_product_matrix(
            cubic_pqs_payload.support_states,
            product_unit.support_states,
            cubic_kinetic_axis_ops.x.kinetic,
            cubic_kinetic_axis_ops.y.overlap,
            cubic_kinetic_axis_ops.z.overlap,
        )
        y_term = _pqs_cross_product_matrix(
            cubic_pqs_payload.support_states,
            product_unit.support_states,
            cubic_kinetic_axis_ops.x.overlap,
            cubic_kinetic_axis_ops.y.kinetic,
            cubic_kinetic_axis_ops.z.overlap,
        )
        z_term = _pqs_cross_product_matrix(
            cubic_pqs_payload.support_states,
            product_unit.support_states,
            cubic_kinetic_axis_ops.x.overlap,
            cubic_kinetic_axis_ops.y.overlap,
            cubic_kinetic_axis_ops.z.kinetic,
        )
        return transpose(pqs_product_support_coefficients) *
               (x_term + y_term + z_term) *
               Matrix{Float64}(product_unit.coefficient_matrix)
    end
    pqs_product_reference_errors = Float64[]
    for term in (:overlap, :position_x, :position_y, :position_z)
        pqs_product_block = CCPM._pqs_product_low_order_reference_block(
            cubic_pqs_payload,
            product_unit,
            cubic_metrics;
            term,
        )
        product_pqs_block = CCPM._pqs_product_low_order_reference_block(
            product_unit,
            cubic_pqs_payload,
            cubic_metrics;
            term,
        )
        expected_pqs_product_block = _pqs_product_expected_low_order_block(term)
        @test pqs_product_block.path == :pqs_product_low_order_reference
        @test pqs_product_block.term == term
        @test size(pqs_product_block.block) == (98, 4)
        @test pqs_product_block.block ≈ pqs_product_block.oracle_block atol = 1.0e-12 rtol = 1.0e-12
        @test pqs_product_block.block ≈ expected_pqs_product_block atol = 1.0e-12 rtol = 1.0e-12
        @test pqs_product_block.coefficient_error < 1.0e-12
        @test pqs_product_block.block_error < 1.0e-12
        @test pqs_product_block.diagnostics.fixture_reference_only
        @test !pqs_product_block.diagnostics.production_supported
        @test !pqs_product_block.diagnostics.packet_adoption
        @test !pqs_product_block.diagnostics.fixed_block_sidecar_installation
        @test !pqs_product_block.diagnostics.qwhamiltonian_consumes
        @test !pqs_product_block.diagnostics.public_default_consumes
        @test !pqs_product_block.diagnostics.cr2_science_status_changed
        @test !pqs_product_block.diagnostics.ida_weight_division_allowed
        @test pqs_product_block.diagnostics.retained_pqs_weights_role ==
              :debug_reference_only
        @test !pqs_product_block.diagnostics.retained_pqs_weights_used
        @test !pqs_product_block.diagnostics.retained_pqs_weights_positive_checked
        @test !pqs_product_block.diagnostics.pqs_self_overlap_identity_shortcut_used
        @test pqs_product_block.diagnostics.factored_pqs_transform_used
        @test pqs_product_block.diagnostics.seed_reconstructed_from_descriptor
        @test pqs_product_block.diagnostics.cleanup_transform_stage_applied
        @test pqs_product_block.diagnostics.support_coefficient_matrix_compared
        @test pqs_product_block.diagnostics.support_local_oracle_used
        @test !pqs_product_block.diagnostics.optimized_pqs_product_path
        @test :kinetic in pqs_product_block.diagnostics.unsupported_terms
        @test product_pqs_block.path == :product_pqs_low_order_reference
        @test product_pqs_block.block ≈ transpose(pqs_product_block.block) atol = 1.0e-12 rtol = 1.0e-12
        @test product_pqs_block.oracle_block ≈ transpose(pqs_product_block.oracle_block) atol = 1.0e-12 rtol = 1.0e-12
        @test product_pqs_block.diagnostics.transposed_from_pqs_product_reference
        @test !product_pqs_block.diagnostics.pqs_self_overlap_identity_shortcut_used
        push!(pqs_product_reference_errors, pqs_product_block.block_error)
    end
    @test maximum(pqs_product_reference_errors) < 1.0e-12
    @test_throws ArgumentError CCPM._pqs_product_low_order_reference_block(
        cubic_pqs_payload,
        product_unit,
        cubic_metrics;
        term = :kinetic,
    )
    pqs_product_kinetic_block = CCPM._pqs_product_kinetic_reference_block(
        cubic_pqs_payload,
        product_unit,
        cubic_kinetic_axis_ops,
    )
    product_pqs_kinetic_block = CCPM._pqs_product_kinetic_reference_block(
        product_unit,
        cubic_pqs_payload,
        cubic_kinetic_axis_ops,
    )
    expected_pqs_product_kinetic = _pqs_product_expected_kinetic_block()
    @test pqs_product_kinetic_block.path == :pqs_product_kinetic_reference
    @test pqs_product_kinetic_block.term == :kinetic
    @test size(pqs_product_kinetic_block.block) == (98, 4)
    @test pqs_product_kinetic_block.block ≈ pqs_product_kinetic_block.oracle_block atol = 1.0e-12 rtol = 1.0e-12
    @test pqs_product_kinetic_block.block ≈ expected_pqs_product_kinetic atol = 1.0e-12 rtol = 1.0e-12
    @test pqs_product_kinetic_block.coefficient_error < 1.0e-12
    @test pqs_product_kinetic_block.block_error < 1.0e-12
    @test pqs_product_kinetic_block.diagnostics.fixture_reference_only
    @test !pqs_product_kinetic_block.diagnostics.production_supported
    @test pqs_product_kinetic_block.diagnostics.signed_operator_reference
    @test pqs_product_kinetic_block.diagnostics.retained_weight_semantics == :not_used
    @test pqs_product_kinetic_block.diagnostics.retained_pqs_weights_role ==
          :debug_reference_only
    @test !pqs_product_kinetic_block.diagnostics.retained_pqs_weights_used
    @test !pqs_product_kinetic_block.diagnostics.retained_pqs_weights_positive_checked
    @test !pqs_product_kinetic_block.diagnostics.ida_weight_division_allowed
    @test !pqs_product_kinetic_block.diagnostics.quadrature_weight_semantics_claimed
    @test !pqs_product_kinetic_block.diagnostics.packet_adoption
    @test !pqs_product_kinetic_block.diagnostics.fixed_block_sidecar_installation
    @test !pqs_product_kinetic_block.diagnostics.qwhamiltonian_consumes
    @test !pqs_product_kinetic_block.diagnostics.public_default_consumes
    @test !pqs_product_kinetic_block.diagnostics.cr2_science_status_changed
    @test pqs_product_kinetic_block.diagnostics.factored_pqs_transform_used
    @test pqs_product_kinetic_block.diagnostics.seed_reconstructed_from_descriptor
    @test pqs_product_kinetic_block.diagnostics.cleanup_transform_stage_applied
    @test pqs_product_kinetic_block.diagnostics.support_coefficient_matrix_compared
    @test pqs_product_kinetic_block.diagnostics.support_local_oracle_used
    @test !pqs_product_kinetic_block.diagnostics.optimized_pqs_product_path
    @test pqs_product_kinetic_block.diagnostics.kinetic_factor_form == (
        (:kinetic, :overlap, :overlap),
        (:overlap, :kinetic, :overlap),
        (:overlap, :overlap, :kinetic),
    )
    @test product_pqs_kinetic_block.path == :product_pqs_kinetic_reference
    @test product_pqs_kinetic_block.block ≈ transpose(pqs_product_kinetic_block.block) atol = 1.0e-12 rtol = 1.0e-12
    @test product_pqs_kinetic_block.oracle_block ≈ transpose(pqs_product_kinetic_block.oracle_block) atol = 1.0e-12 rtol = 1.0e-12
    @test product_pqs_kinetic_block.diagnostics.transposed_from_pqs_product_reference
    @test_throws ArgumentError CCPM._pqs_product_kinetic_reference_block(
        product_unit,
        product_unit,
        cubic_kinetic_axis_ops,
    )
    @test_throws ArgumentError CCPM._pqs_product_kinetic_reference_block(
        cubic_pqs_payload,
        support_dense_unit,
        cubic_kinetic_axis_ops,
    )
    generic_product_source_transform = CCP._cartesian_raw_product_source_retained_transform(
        product_unit;
        source_id = :generic_product_slab_source,
        parent_dims = (2, 2, 1),
    )
    @test generic_product_source_transform.raw_source.source_id ==
          :generic_product_slab_source
    @test generic_product_source_transform.raw_source.source_dimension == 4
    @test generic_product_source_transform.raw_source.support_indices == collect(1:4)
    @test generic_product_source_transform.raw_source.raw_source_weight_role ==
          :raw_source_positive
    @test generic_product_source_transform.retained_transform.source_id ==
          :generic_product_slab_source
    @test generic_product_source_transform.retained_transform.transform_kind ==
          :product_axis_transform
    nonidentity_transform = [
        inv(sqrt(2.0)) 0.0
        0.0 inv(sqrt(2.0))
        inv(sqrt(2.0)) 0.0
        0.0 -inv(sqrt(2.0))
    ]
    nonidentity_product_unit = GaussletBases._CartesianNestedProductStagedByCenterUnit3D(
        :nonidentity_product_slab,
        :product_doside,
        1:2,
        collect(1:4),
        NTuple{3,Int}[(1, 1, 1), (1, 2, 1), (2, 1, 1), (2, 2, 1)],
        nonidentity_transform,
        (
            GaussletBases._nested_product_staged_active_axis(1:2, identity_axis),
            GaussletBases._nested_product_staged_active_axis(1:2, identity_axis),
            GaussletBases._nested_product_staged_fixed_axis(1),
        ),
        NTuple{3,Int}[(1, 1, 1), (1, 2, 1)],
        (source = :nonidentity_raw_product_source_test_fixture,),
        (support_count = 4, retained_count = 2),
    )
    nonidentity_product_source_transform = CCP._cartesian_raw_product_source_retained_transform(
        nonidentity_product_unit;
        source_id = :nonidentity_product_slab_source,
        parent_dims = (2, 2, 1),
    )
    @test nonidentity_product_source_transform.raw_source.source_id ==
          :nonidentity_product_slab_source
    @test nonidentity_product_source_transform.raw_source.source_dimension == 4
    @test nonidentity_product_source_transform.raw_source.support_indices == collect(1:4)
    @test nonidentity_product_source_transform.retained_transform.source_id ==
          :nonidentity_product_slab_source
    @test nonidentity_product_source_transform.retained_transform.retained_dimension == 2
    @test nonidentity_product_source_transform.retained_transform.transform_kind ==
          :product_axis_transform
    @test nonidentity_product_source_transform.retained_transform.diagnostics.full_raw_to_retained_matrix_materialized
    mixed_pair_plan = CCP._cartesian_raw_product_source_pair_plan(
        (installed_source_transform, product_source_transform);
        operator_kind = :low_order_metric,
        supported_terms = (:overlap, :axis_index_x, :position_x, :position_y, :position_z),
        source = :pqs_product_raw_source_pair_plan_test,
    )
    mixed_source_ids = sort!(
        [pqs_raw_source.source_id, :identity_product_slab_source];
        by = string,
    )
    expected_mixed_pair_keys = NTuple{2,Symbol}[]
    for right_index in eachindex(mixed_source_ids)
        right_id = mixed_source_ids[right_index]
        for left_index in 1:right_index
            push!(expected_mixed_pair_keys, (mixed_source_ids[left_index], right_id))
        end
    end
    @test length(mixed_pair_plan.raw_sources) == 2
    @test length(mixed_pair_plan.retained_transforms) == 2
    @test length(mixed_pair_plan.pair_packets) == 3
    @test mixed_pair_plan.source_ids == mixed_source_ids
    @test all(id -> haskey(mixed_pair_plan.raw_sources, id), mixed_source_ids)
    @test all(id -> haskey(mixed_pair_plan.retained_transforms, id), mixed_source_ids)
    @test mixed_pair_plan.pair_keys == expected_mixed_pair_keys
    @test mixed_pair_plan.diagnostics.source == :pqs_product_raw_source_pair_plan_test
    @test mixed_pair_plan.diagnostics.expected_upper_triangle_pair_count == 3
    @test mixed_pair_plan.diagnostics.upper_triangle_only
    @test mixed_pair_plan.diagnostics.all_pairs_resolve_sources
    @test mixed_pair_plan.diagnostics.all_pairs_resolve_retained_transforms
    @test mixed_pair_plan.diagnostics.pair_packets_placeholder_only
    @test all(packet -> packet.operator_matrices === nothing, mixed_pair_plan.pair_packets)
    @test all(
        packet -> packet.symmetry_status == :symmetric_upper_triangle_placeholder,
        mixed_pair_plan.pair_packets,
    )
    mixed_resolved_pairs =
        CCP._cartesian_resolved_raw_product_source_pairs(mixed_pair_plan)
    @test length(mixed_resolved_pairs) == 3
    @test [pair.pair_key for pair in mixed_resolved_pairs] == expected_mixed_pair_keys
    @test all(pair -> pair.diagnostics.upper_triangular, mixed_resolved_pairs)
    @test all(pair -> pair.diagnostics.placeholder_only, mixed_resolved_pairs)
    @test all(pair -> pair.diagnostics.raw_weight_roles_explicit, mixed_resolved_pairs)
    @test all(
        pair -> pair.diagnostics.retained_weight_roles_explicit,
        mixed_resolved_pairs,
    )
    @test all(
        pair -> !pair.diagnostics.retained_ida_weight_division_allowed,
        mixed_resolved_pairs,
    )
    @test all(
        pair -> pair.left_raw_source.source_id == pair.pair_key[1] &&
                pair.right_raw_source.source_id == pair.pair_key[2],
        mixed_resolved_pairs,
    )
    @test all(
        pair -> pair.left_retained_transform.source_id == pair.pair_key[1] &&
                pair.right_retained_transform.source_id == pair.pair_key[2],
        mixed_resolved_pairs,
    )
    product_self_pair = only(
        pair for pair in mixed_resolved_pairs
        if pair.pair_key == (:identity_product_slab_source, :identity_product_slab_source)
    )
    generic_product_pair_plan = CCP._cartesian_raw_product_source_pair_plan(
        (generic_product_source_transform,);
        operator_kind = :low_order_metric,
        supported_terms = (:overlap, :axis_index_x, :position_x, :position_y, :position_z),
        source = :generic_product_raw_source_pair_plan_test,
    )
    generic_product_self_pair = only(
        CCP._cartesian_resolved_raw_product_source_pairs(generic_product_pair_plan)
    )
    @test generic_product_self_pair.pair_key ==
          (:generic_product_slab_source, :generic_product_slab_source)
    nonidentity_product_pair_plan = CCP._cartesian_raw_product_source_pair_plan(
        (nonidentity_product_source_transform,);
        operator_kind = :low_order_metric,
        supported_terms = (:overlap, :axis_index_x, :position_x, :position_y, :position_z),
        source = :nonidentity_product_raw_source_pair_plan_test,
    )
    nonidentity_product_self_pair = only(
        CCP._cartesian_resolved_raw_product_source_pairs(nonidentity_product_pair_plan)
    )
    @test nonidentity_product_self_pair.pair_key ==
          (:nonidentity_product_slab_source, :nonidentity_product_slab_source)
    cross_product_pair_plan = CCP._cartesian_raw_product_source_pair_plan(
        (product_source_transform, nonidentity_product_source_transform);
        operator_kind = :low_order_metric,
        supported_terms = (:position_x, :position_y, :position_z),
        source = :identity_nonidentity_product_cross_pair_plan_test,
    )
    cross_product_pairs =
        CCP._cartesian_resolved_raw_product_source_pairs(cross_product_pair_plan)
    cross_product_pair = only(
        pair for pair in cross_product_pairs
        if pair.pair_key ==
           (:identity_product_slab_source, :nonidentity_product_slab_source)
    )
    @test length(cross_product_pairs) == 3
    @test cross_product_pair.left_raw_source.source_id ==
          :identity_product_slab_source
    @test cross_product_pair.right_raw_source.source_id ==
          :nonidentity_product_slab_source
    @test cross_product_pair.left_raw_source.source_dimension == 4
    @test cross_product_pair.right_raw_source.source_dimension == 4
    @test cross_product_pair.left_retained_transform.retained_dimension == 4
    @test cross_product_pair.right_retained_transform.retained_dimension == 2
    @test !cross_product_pair.diagnostics.pqs_factored_transform_present
    cross_overlap_pair_plan = CCP._cartesian_raw_product_source_pair_plan(
        (product_source_transform, nonidentity_product_source_transform);
        operator_kind = :low_order_metric,
        supported_terms = (:overlap,),
        source = :identity_nonidentity_product_cross_overlap_pair_plan_test,
    )
    cross_overlap_pair = only(
        pair for pair in CCP._cartesian_resolved_raw_product_source_pairs(
            cross_overlap_pair_plan,
        )
        if pair.pair_key ==
           (:identity_product_slab_source, :nonidentity_product_slab_source)
    )
    position_omitting_pair_plan = CCP._cartesian_raw_product_source_pair_plan(
        (product_source_transform,);
        operator_kind = :low_order_metric,
        supported_terms = (:overlap, :axis_index_x),
        source = :position_omitting_product_raw_source_pair_plan_test,
    )
    position_omitting_pair = only(
        CCP._cartesian_resolved_raw_product_source_pairs(position_omitting_pair_plan)
    )
    axis_index_omitting_pair_plan = CCP._cartesian_raw_product_source_pair_plan(
        (product_source_transform,);
        operator_kind = :low_order_metric,
        supported_terms = (:overlap, :position_x, :position_y, :position_z),
        source = :axis_index_omitting_product_raw_source_pair_plan_test,
    )
    axis_index_omitting_pair = only(
        CCP._cartesian_resolved_raw_product_source_pairs(axis_index_omitting_pair_plan)
    )
    overlap_omitting_pair_plan = CCP._cartesian_raw_product_source_pair_plan(
        (product_source_transform,);
        operator_kind = :low_order_metric,
        supported_terms = (:position_x, :position_y, :position_z),
        source = :overlap_omitting_product_raw_source_pair_plan_test,
    )
    overlap_omitting_pair = only(
        CCP._cartesian_resolved_raw_product_source_pairs(overlap_omitting_pair_plan)
    )
    cross_position_omitting_pair_plan = CCP._cartesian_raw_product_source_pair_plan(
        (product_source_transform, nonidentity_product_source_transform);
        operator_kind = :low_order_metric,
        supported_terms = (:position_y, :position_z),
        source = :position_omitting_product_cross_pair_plan_test,
    )
    cross_position_omitting_pair = only(
        pair for pair in CCP._cartesian_resolved_raw_product_source_pairs(
            cross_position_omitting_pair_plan,
        )
        if pair.pair_key ==
           (:identity_product_slab_source, :nonidentity_product_slab_source)
    )
    product_raw_overlap_packet = CCP._cartesian_raw_low_order_operator_packet(
        product_self_pair;
        term = :overlap,
    )
    @test product_raw_overlap_packet.source_dimensions == (4, 4)
    @test product_raw_overlap_packet.raw_operator_matrix == Matrix{Float64}(I, 4, 4)
    expected_product_axis_index_x = Float64[
        1 0 0 0
        0 1 0 0
        0 0 2 0
        0 0 0 2
    ]
    product_raw_axis_index_x_packet = CCP._cartesian_raw_low_order_operator_packet(
        product_self_pair;
        term = :axis_index_x,
    )
    @test product_raw_axis_index_x_packet.left_source_id == :identity_product_slab_source
    @test product_raw_axis_index_x_packet.right_source_id == :identity_product_slab_source
    @test product_raw_axis_index_x_packet.operator_kind == :low_order_metric
    @test product_raw_axis_index_x_packet.term == :axis_index_x
    @test product_raw_axis_index_x_packet.source_dimensions == (4, 4)
    @test product_raw_axis_index_x_packet.raw_operator_matrix == expected_product_axis_index_x
    @test product_raw_axis_index_x_packet.diagnostics.raw_reference ==
          :identity_product_slab_axis_index_x
    @test product_raw_axis_index_x_packet.diagnostics.raw_reference_error == 0.0
    @test product_raw_axis_index_x_packet.diagnostics.separable_axis_metadata_used
    @test product_raw_axis_index_x_packet.diagnostics.axis_index_diagnostic
    @test !product_raw_axis_index_x_packet.diagnostics.physical_position_operator
    @test !product_raw_axis_index_x_packet.diagnostics.dense_parent_matrix_used
    @test product_raw_axis_index_x_packet.diagnostics.raw_operator_matrix_built
    @test !product_raw_axis_index_x_packet.diagnostics.retained_operator_block_built
    @test !product_raw_axis_index_x_packet.diagnostics.retained_transform_applied
    @test_throws ArgumentError CCP._cartesian_raw_low_order_operator_packet(
        product_self_pair;
        term = :position_x,
    )
    @test_throws ArgumentError CCP._cartesian_raw_low_order_operator_packet(
        axis_index_omitting_pair;
        term = :axis_index_x,
    )
    @test_throws ArgumentError CCP._cartesian_raw_low_order_operator_packet(
        generic_product_self_pair;
        term = :axis_index_x,
    )
    @test_throws ArgumentError CCP._cartesian_raw_low_order_operator_packet(
        nonidentity_product_self_pair;
        term = :axis_index_x,
    )
    @test_throws ArgumentError CCP._cartesian_raw_low_order_operator_packet(
        pqs_resolved_pair;
        term = :axis_index_x,
    )
    physical_axis_metrics = (
        x = (
            overlap = Matrix{Float64}(I, 2, 2),
            position = Matrix(Diagonal([0.25, 1.75])),
            x2 = Matrix(Diagonal([0.25^2, 1.75^2])),
            kinetic = [
                2.0 -0.25
                -0.25 2.5
            ],
            weights = [1.0, 1.0],
            centers = [0.25, 1.75],
            source = :synthetic_noninteger_axis_metric_fixture,
        ),
        y = (
            overlap = Matrix{Float64}(I, 2, 2),
            position = Matrix(Diagonal([-0.5, 0.5])),
            x2 = Matrix(Diagonal([0.25, 0.25])),
            kinetic = [
                1.5 0.125
                0.125 1.25
            ],
            weights = [1.0, 1.0],
            centers = [-0.5, 0.5],
            source = :synthetic_noninteger_axis_metric_fixture,
        ),
        z = (
            overlap = Matrix{Float64}(I, 1, 1),
            position = Matrix(Diagonal([3.25])),
            x2 = Matrix(Diagonal([3.25^2])),
            kinetic = Matrix(Diagonal([0.75])),
            weights = [1.0],
            centers = [3.25],
            source = :synthetic_noninteger_axis_metric_fixture,
        ),
    )
    product_support_states = product_self_pair.left_raw_source.provenance.unit.support_states
    function support_local_physical_position_reference(
        term::Symbol,
        left_support_states = product_support_states,
        right_support_states = product_support_states,
    )
        axis_matrices = term == :position_x ? (
            physical_axis_metrics.x.position,
            physical_axis_metrics.y.overlap,
            physical_axis_metrics.z.overlap,
        ) : term == :position_y ? (
            physical_axis_metrics.x.overlap,
            physical_axis_metrics.y.position,
            physical_axis_metrics.z.overlap,
        ) : term == :position_z ? (
            physical_axis_metrics.x.overlap,
            physical_axis_metrics.y.overlap,
            physical_axis_metrics.z.position,
        ) : throw(ArgumentError("unsupported test physical position term"))
        reference =
            zeros(Float64, length(left_support_states), length(right_support_states))
        for col in eachindex(right_support_states)
            jx, jy, jz = right_support_states[col]
            for row in eachindex(left_support_states)
                ix, iy, iz = left_support_states[row]
                reference[row, col] =
                    axis_matrices[1][ix, jx] *
                    axis_matrices[2][iy, jy] *
                    axis_matrices[3][iz, jz]
            end
        end
        return reference
    end
    function support_local_overlap_reference(
        axis_metrics,
        left_support_states = product_support_states,
        right_support_states = product_support_states,
    )
        axis_matrices = (
            axis_metrics.x.overlap,
            axis_metrics.y.overlap,
            axis_metrics.z.overlap,
        )
        reference =
            zeros(Float64, length(left_support_states), length(right_support_states))
        for col in eachindex(right_support_states)
            jx, jy, jz = right_support_states[col]
            for row in eachindex(left_support_states)
                ix, iy, iz = left_support_states[row]
                reference[row, col] =
                    axis_matrices[1][ix, jx] *
                    axis_matrices[2][iy, jy] *
                    axis_matrices[3][iz, jz]
            end
        end
        return reference
    end
    expected_physical_positions = (
        position_x = Diagonal([0.25, 0.25, 1.75, 1.75]),
        position_y = Diagonal([-0.5, 0.5, -0.5, 0.5]),
        position_z = Diagonal([3.25, 3.25, 3.25, 3.25]),
    )
    factorized_product_overlap_packet =
        CCP._cartesian_factorized_product_doside_raw_low_order_operator_packet(
            product_self_pair;
            term = :overlap,
            axis_metrics = physical_axis_metrics,
        )
    factorized_cross_overlap_packet =
        CCP._cartesian_factorized_product_doside_raw_low_order_operator_packet(
            cross_overlap_pair;
            term = :overlap,
            axis_metrics = physical_axis_metrics,
        )
    product_overlap_reference = support_local_overlap_reference(physical_axis_metrics)
    cross_overlap_reference = support_local_overlap_reference(
        physical_axis_metrics,
        cross_overlap_pair.left_raw_source.provenance.unit.support_states,
        cross_overlap_pair.right_raw_source.provenance.unit.support_states,
    )
    @test factorized_product_overlap_packet.left_source_id ==
          :identity_product_slab_source
    @test factorized_product_overlap_packet.right_source_id ==
          :identity_product_slab_source
    @test factorized_product_overlap_packet.term == :overlap
    @test factorized_product_overlap_packet.source_dimensions == (4, 4)
    @test factorized_product_overlap_packet.raw_operator_matrix ≈
          product_overlap_reference atol = 1.0e-14 rtol = 1.0e-14
    @test factorized_product_overlap_packet.diagnostics.raw_reference ==
          :product_doside_factorized_overlap
    @test factorized_product_overlap_packet.diagnostics.factorized_axis_path_used
    @test factorized_product_overlap_packet.diagnostics.axis_factor_scope ==
          :staged_local_axis_intervals
    @test !factorized_product_overlap_packet.diagnostics.support_row_reference_used
    @test factorized_product_overlap_packet.diagnostics.raw_basis_scope ==
          :raw_product_source_rows
    @test factorized_product_overlap_packet.diagnostics.fixture_only
    @test !factorized_product_overlap_packet.diagnostics.production_supported
    @test !factorized_product_overlap_packet.diagnostics.dense_parent_matrix_used
    @test !factorized_product_overlap_packet.diagnostics.metric_execution_changed
    @test !factorized_product_overlap_packet.diagnostics.retained_positive_weight_claim
    @test !factorized_product_overlap_packet.diagnostics.ida_weight_division_allowed
    @test factorized_cross_overlap_packet.left_source_id ==
          :identity_product_slab_source
    @test factorized_cross_overlap_packet.right_source_id ==
          :nonidentity_product_slab_source
    @test factorized_cross_overlap_packet.source_dimensions == (4, 4)
    @test factorized_cross_overlap_packet.raw_operator_matrix ≈
          cross_overlap_reference atol = 1.0e-14 rtol = 1.0e-14
    @test factorized_cross_overlap_packet.diagnostics.cross_pair
    @test factorized_cross_overlap_packet.diagnostics.left_source_dimension == 4
    @test factorized_cross_overlap_packet.diagnostics.right_source_dimension == 4
    @test factorized_cross_overlap_packet.diagnostics.factorized_axis_path_used
    factorized_product_retained_overlap =
        CCP._cartesian_retained_low_order_operator_block(
            factorized_product_overlap_packet,
            product_source_transform.retained_transform,
        )
    factorized_cross_retained_overlap =
        CCP._cartesian_retained_low_order_operator_block(
            factorized_cross_overlap_packet,
            product_source_transform.retained_transform,
            nonidentity_product_source_transform.retained_transform,
        )
    @test factorized_product_retained_overlap.retained_operator_matrix ≈
          product_overlap_reference atol = 1.0e-14 rtol = 1.0e-14
    @test factorized_cross_retained_overlap.retained_dimensions == (4, 2)
    @test factorized_cross_retained_overlap.retained_operator_matrix ≈
          cross_overlap_reference * nonidentity_transform atol = 1.0e-14 rtol = 1.0e-14
    product_direct_self_overlap = CCPM._product_doside_retained_low_order_block(
        product_unit,
        product_unit,
        physical_axis_metrics;
        term = :overlap,
    )
    product_source_box_pair_plan =
        CCPM._product_doside_source_box_pair_plan(
            product_unit,
            product_unit,
            physical_axis_metrics,
        )
    @test product_source_box_pair_plan.pair_kind ==
          :product_doside_source_box_pair
    @test product_source_box_pair_plan.left_source_family == :product_doside
    @test product_source_box_pair_plan.right_source_family == :product_doside
    @test product_source_box_pair_plan.left_source_dimensions == (2, 2, 1)
    @test product_source_box_pair_plan.right_source_dimensions == (2, 2, 1)
    @test product_source_box_pair_plan.left_source_dimension == 4
    @test product_source_box_pair_plan.right_source_dimension == 4
    @test product_source_box_pair_plan.left_retained_axis_counts == (2, 2, 1)
    @test product_source_box_pair_plan.right_retained_axis_counts == (2, 2, 1)
    @test product_source_box_pair_plan.left_retained_count == 4
    @test product_source_box_pair_plan.right_retained_count == 4
    @test product_source_box_pair_plan.left_column_range == product_unit.column_range
    @test product_source_box_pair_plan.right_column_range == product_unit.column_range
    @test product_source_box_pair_plan.axis_intervals.left == (1:2, 1:2, 1:1)
    @test product_source_box_pair_plan.axis_intervals.right == (1:2, 1:2, 1:1)
    @test product_source_box_pair_plan.axis_centers.left[1] == [0.25, 1.75]
    @test product_source_box_pair_plan.axis_centers.right[2] == [-0.5, 0.5]
    @test product_source_box_pair_plan.left_retained_unit_plan.object_kind ==
          :product_doside_retained_unit_plan
    @test product_source_box_pair_plan.right_retained_unit_plan.object_kind ==
          :product_doside_retained_unit_plan
    @test product_source_box_pair_plan.left_retained_unit_plan.retained_rule_kind ==
          :product_doside
    @test product_source_box_pair_plan.right_retained_unit_plan.retained_rule_kind ==
          :product_doside
    @test product_source_box_pair_plan.left_retained_transform ===
          product_source_box_pair_plan.left_retained_unit_plan
    @test product_source_box_pair_plan.right_retained_transform ===
          product_source_box_pair_plan.right_retained_unit_plan
    @test size(product_source_box_pair_plan.one_dimensional_cross_factors.x.overlap) ==
          (2, 2)
    @test size(product_source_box_pair_plan.one_dimensional_cross_factors.y.position) ==
          (2, 2)
    @test size(product_source_box_pair_plan.one_dimensional_cross_factors.x.x2) ==
          (2, 2)
    @test size(product_source_box_pair_plan.one_dimensional_cross_factors.x.kinetic) ==
          (2, 2)
    @test size(product_source_box_pair_plan.one_dimensional_cross_factors.z.overlap) ==
          (1, 1)
    @test product_source_box_pair_plan.diagnostics.private_shadow_only
    @test product_source_box_pair_plan.diagnostics.source_box_pair_plan
    @test product_source_box_pair_plan.diagnostics.product_doside_retained_unit_plan_used
    @test product_source_box_pair_plan.diagnostics.existing_product_staged_retained_helpers_authoritative
    @test product_source_box_pair_plan.diagnostics.operator_factor_source ==
          :explicit_metric_operator_data
    @test product_source_box_pair_plan.diagnostics.operator_metric_sources ==
          (
              :synthetic_noninteger_axis_metric_fixture,
              :synthetic_noninteger_axis_metric_fixture,
              :synthetic_noninteger_axis_metric_fixture,
          )
    @test !product_source_box_pair_plan.diagnostics.numerical_reference_fallback
    @test !product_source_box_pair_plan.diagnostics.product_staged_metric_execution_changed
    @test !product_source_box_pair_plan.diagnostics.product_doside_retained_block_math_changed
    @test product_source_box_pair_plan.diagnostics.raw_product_box_operators_use_1d_factors
    @test !product_source_box_pair_plan.diagnostics.packet_adoption
    @test !product_source_box_pair_plan.diagnostics.fixed_block_routing
    @test !product_source_box_pair_plan.diagnostics.qwhamiltonian_consumes
    @test !product_source_box_pair_plan.diagnostics.public_default_consumes
    @test !product_source_box_pair_plan.diagnostics.ida_mwg_semantics_changed
    @test !product_source_box_pair_plan.diagnostics.retained_weight_semantics_changed
    @test !product_source_box_pair_plan.diagnostics.generic_retained_unit_framework
    for term in (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
        source_box_block = CCPM._product_doside_source_box_reference_block(
            product_unit,
            product_unit,
            physical_axis_metrics;
            term,
        )
        direct_block = term == :kinetic ?
                       CCPM._product_doside_retained_kinetic_block(
                           product_unit,
                           product_unit,
                           physical_axis_metrics,
                       ) :
                       CCPM._product_doside_retained_low_order_block(
                           product_unit,
                           product_unit,
                           physical_axis_metrics;
                           term,
                       )
        explicit_block = _product_source_box_explicit_reference(
            product_unit,
            physical_axis_metrics;
            term,
        )
        @test source_box_block.path == :product_doside_source_box_reference
        @test source_box_block.term == term
        @test source_box_block.block ≈ direct_block atol = 1.0e-14 rtol = 1.0e-14
        @test source_box_block.block ≈ explicit_block atol = 1.0e-14 rtol = 1.0e-14
        @test source_box_block.authoritative_block ≈ direct_block atol = 0.0 rtol = 0.0
        @test source_box_block.block_error <= 1.0e-14
        @test source_box_block.diagnostics.authoritative_helper ==
              (
                  term == :kinetic ?
                  :_product_doside_retained_kinetic_block :
                  :_product_doside_retained_low_order_block
              )
        @test source_box_block.diagnostics.authoritative_block_compared
        @test source_box_block.diagnostics.operator_factor_source ==
              :explicit_metric_operator_data
        @test !source_box_block.diagnostics.numerical_reference_fallback
        @test source_box_block.diagnostics.kinetic_factor_form ==
              (
                  (:kinetic, :overlap, :overlap),
                  (:overlap, :kinetic, :overlap),
                  (:overlap, :overlap, :kinetic),
              )
        @test source_box_block.diagnostics.private_shadow_only
        @test !source_box_block.diagnostics.product_staged_metric_execution_changed
        @test !source_box_block.diagnostics.packet_adoption
        @test !source_box_block.diagnostics.fixed_block_routing
        @test !source_box_block.diagnostics.qwhamiltonian_consumes
        @test !source_box_block.diagnostics.public_default_consumes
        @test !source_box_block.diagnostics.ida_mwg_semantics_changed
        @test !source_box_block.diagnostics.retained_weight_semantics_changed
        @test !source_box_block.diagnostics.generic_retained_unit_framework
    end
    nonidentity_product_source_box_pair_plan =
        CCPM._product_doside_source_box_pair_plan(
            nonidentity_axis_product_unit,
            product_unit,
            physical_axis_metrics,
        )
    @test nonidentity_product_source_box_pair_plan.left_retained_unit_plan.axis_coefficient_matrices[1] ==
          nonidentity_axis_product_unit.axes[1].coefficient_matrix
    @test nonidentity_product_source_box_pair_plan.right_retained_unit_plan.axis_coefficient_matrices[1] ==
          product_unit.axes[1].coefficient_matrix
    nonidentity_product_source_box_block =
        CCPM._product_doside_source_box_reference_block(
            nonidentity_axis_product_unit,
            product_unit,
            physical_axis_metrics;
            term = :kinetic,
        )
    nonidentity_product_direct_block =
        CCPM._product_doside_retained_kinetic_block(
            nonidentity_axis_product_unit,
            product_unit,
            physical_axis_metrics,
        )
    @test nonidentity_product_source_box_block.block ≈
          nonidentity_product_direct_block atol = 1.0e-14 rtol = 1.0e-14
    @test nonidentity_product_source_box_block.diagnostics.authoritative_block_compared
    product_source_box_shadow = CCPM._product_doside_source_box_shadow_blocks(
        nonidentity_axis_product_unit,
        product_unit,
        physical_axis_metrics;
        terms = (:overlap, :position_x, :x2_x, :kinetic),
    )
    @test product_source_box_shadow.path ==
          :product_doside_source_box_shadow_blocks
    @test product_source_box_shadow.terms ==
          (:overlap, :position_x, :x2_x, :kinetic)
    @test product_source_box_shadow.retained_dimension == 8
    @test product_source_box_shadow.ranges.left == 1:4
    @test product_source_box_shadow.ranges.right == 5:8
    @test product_source_box_shadow.diagnostics.source_box_shadow_only
    @test product_source_box_shadow.diagnostics.private_shadow_only
    @test product_source_box_shadow.diagnostics.existing_product_staged_retained_helpers_authoritative
    @test product_source_box_shadow.diagnostics.operator_factor_source ==
          :explicit_metric_operator_data
    @test !product_source_box_shadow.diagnostics.numerical_reference_fallback
    @test product_source_box_shadow.diagnostics.component_reference_helper ==
          :_product_doside_source_box_reference_block
    @test product_source_box_shadow.diagnostics.right_left_block_source ==
          :_product_doside_source_box_reference_block
    @test product_source_box_shadow.diagnostics.symmetric_real_transpose_checked
    @test product_source_box_shadow.diagnostics.max_right_left_transpose_error <= 1.0e-14
    @test !product_source_box_shadow.diagnostics.product_staged_metric_execution_changed
    @test !product_source_box_shadow.diagnostics.product_doside_retained_block_math_changed
    @test !product_source_box_shadow.diagnostics.packet_adoption
    @test !product_source_box_shadow.diagnostics.fixed_block_routing
    @test !product_source_box_shadow.diagnostics.qwhamiltonian_consumes
    @test !product_source_box_shadow.diagnostics.public_default_consumes
    @test !product_source_box_shadow.diagnostics.ida_mwg_semantics_changed
    @test !product_source_box_shadow.diagnostics.retained_weight_semantics_changed
    @test !product_source_box_shadow.diagnostics.generic_retained_unit_framework
    @test product_source_box_shadow.diagnostics.output_finite
    for term in product_source_box_shadow.terms
        block = product_source_box_shadow.blocks[term]
        components = product_source_box_shadow.component_blocks[term]
        expected_left_left = CCPM._product_doside_source_box_reference_block(
            nonidentity_axis_product_unit,
            nonidentity_axis_product_unit,
            physical_axis_metrics;
            term,
        ).block
        expected_left_right = CCPM._product_doside_source_box_reference_block(
            nonidentity_axis_product_unit,
            product_unit,
            physical_axis_metrics;
            term,
        ).block
        expected_right_left = CCPM._product_doside_source_box_reference_block(
            product_unit,
            nonidentity_axis_product_unit,
            physical_axis_metrics;
            term,
        ).block
        expected_right_right = CCPM._product_doside_source_box_reference_block(
            product_unit,
            product_unit,
            physical_axis_metrics;
            term,
        ).block
        @test all(isfinite, block)
        @test size(block) == (8, 8)
        @test components.left_left ≈ expected_left_left atol = 1.0e-14 rtol = 1.0e-14
        @test components.left_right ≈ expected_left_right atol = 1.0e-14 rtol = 1.0e-14
        @test components.right_left ≈ expected_right_left atol = 1.0e-14 rtol = 1.0e-14
        @test components.right_right ≈ expected_right_right atol = 1.0e-14 rtol = 1.0e-14
        @test block[1:4, 1:4] ≈ components.left_left atol = 1.0e-14 rtol = 1.0e-14
        @test block[1:4, 5:8] ≈ components.left_right atol = 1.0e-14 rtol = 1.0e-14
        @test block[5:8, 1:4] ≈ components.right_left atol = 1.0e-14 rtol = 1.0e-14
        @test block[5:8, 5:8] ≈ components.right_right atol = 1.0e-14 rtol = 1.0e-14
        @test components.right_left_transpose_error <= 1.0e-14
        @test product_source_box_shadow.transpose_errors[term] <= 1.0e-14
    end
    density_density_term_coefficients = [0.6, 0.25]
    function _test_pair_term_tensor(first_term, second_term)
        first_matrix = Matrix{Float64}(first_term)
        second_matrix = Matrix{Float64}(second_term)
        terms = Array{Float64,3}(undef, 2, size(first_matrix, 1), size(first_matrix, 2))
        terms[1, :, :] .= first_matrix
        terms[2, :, :] .= second_matrix
        return terms
    end
    density_normalized_pair_terms = (
        x = _test_pair_term_tensor(
            [1.0 0.2; 0.2 1.5],
            [0.7 -0.1; -0.1 0.9],
        ),
        y = _test_pair_term_tensor(
            [0.8 0.05; 0.05 1.1],
            [1.2 0.15; 0.15 0.6],
        ),
        z = reshape([1.3, 0.4], 2, 1, 1),
    )
    density_source_weights = (
        x = [1.25, 0.75],
        y = [0.9, 1.1],
        z = [2.0],
    )
    function _product_source_box_density_density_explicit_reference(
        left_unit,
        right_unit,
        term_coefficients,
        axis_pair_factor_terms,
    )
        projected_terms = ntuple(axis -> begin
            axis_name = (:x, :y, :z)[axis]
            terms = getproperty(axis_pair_factor_terms, axis_name)
            left_axis = left_unit.axes[axis]
            right_axis = right_unit.axes[axis]
            left_interval = CCPM._staged_axis_interval(left_axis)
            right_interval = CCPM._staged_axis_interval(right_axis)
            projected = Array{Float64,3}(
                undef,
                size(terms, 1),
                size(left_axis.coefficient_matrix, 2),
                size(right_axis.coefficient_matrix, 2),
            )
            for term_index in axes(terms, 1)
                local_terms = Matrix{Float64}(
                    @view terms[term_index, left_interval, right_interval]
                )
                projected[term_index, :, :] .=
                    transpose(left_axis.coefficient_matrix) *
                    local_terms *
                    right_axis.coefficient_matrix
            end
            projected
        end, 3)
        projected_x, projected_y, projected_z = projected_terms
        left_modes = left_unit.axis_function_indices
        right_modes = right_unit.axis_function_indices
        reference = zeros(Float64, length(left_modes), length(right_modes))
        for col in eachindex(right_modes)
            xj, yj, zj = right_modes[col]
            for row in eachindex(left_modes)
                xi, yi, zi = left_modes[row]
                reference[row, col] = sum(
                    term_coefficients[term] *
                    projected_x[term, xi, xj] *
                    projected_y[term, yi, yj] *
                    projected_z[term, zi, zj]
                    for term in eachindex(term_coefficients)
                )
            end
        end
        return reference
    end
    product_density_density =
        CCPM._product_doside_source_box_density_density_interaction_block(
            product_unit,
            product_unit;
            term_coefficients = density_density_term_coefficients,
            axis_pair_factor_terms = density_normalized_pair_terms,
            axis_weights = density_source_weights,
        )
    product_density_reference =
        _product_source_box_density_density_explicit_reference(
            product_unit,
            product_unit,
            density_density_term_coefficients,
            density_normalized_pair_terms,
        )
    @test product_density_density.path ==
          :product_doside_source_box_density_density_interaction
    @test product_density_density.interaction_operator ==
          :electron_electron_density_density
    @test product_density_density.block ≈ product_density_reference atol = 1.0e-14 rtol = 1.0e-14
    @test product_density_density.pair_plan.pair_kind ==
          :product_doside_source_box_density_density_pair
    @test product_density_density.pair_plan.supported_terms == (:pair_sum,)
    @test product_density_density.pair_plan.source_weights.left[1] == [1.25, 0.75]
    @test product_density_density.pair_plan.source_weights.right[2] == [0.9, 1.1]
    @test product_density_density.diagnostics.interaction_operator ==
          :electron_electron_density_density
    @test product_density_density.diagnostics.output_representation ==
          :two_index_density_density
    @test !product_density_density.diagnostics.four_index_galerkin_tensor
    @test product_density_density.diagnostics.raw_source_weights_available
    @test product_density_density.diagnostics.density_normalized_pair_factors
    @test !product_density_density.diagnostics.raw_weighted_pair_factors
    @test product_density_density.diagnostics.source_weight_division_owner ==
          :caller_supplied_density_normalized_pair_factors
    @test !product_density_density.diagnostics.source_weight_division_applied_by_helper
    @test !product_density_density.diagnostics.retained_weight_division_allowed
    @test !product_density_density.diagnostics.retained_pqs_weight_division_allowed
    @test !product_density_density.diagnostics.ida_weight_division_allowed
    @test !product_density_density.diagnostics.packet_adoption
    @test !product_density_density.diagnostics.fixed_block_routing
    @test !product_density_density.diagnostics.qwhamiltonian_consumes
    @test !product_density_density.diagnostics.public_default_consumes
    @test !product_density_density.diagnostics.mwg_ida_semantics_changed
    @test !product_density_density.diagnostics.numerical_reference_fallback
    @test product_density_density.diagnostics.electron_electron_terms_implemented
    nonidentity_density_density =
        CCPM._product_doside_source_box_density_density_interaction_block(
            nonidentity_axis_product_unit,
            product_unit;
            term_coefficients = density_density_term_coefficients,
            axis_pair_factor_terms = density_normalized_pair_terms,
            axis_weights = density_source_weights,
        )
    nonidentity_density_reference =
        _product_source_box_density_density_explicit_reference(
            nonidentity_axis_product_unit,
            product_unit,
            density_density_term_coefficients,
            density_normalized_pair_terms,
        )
    @test nonidentity_density_density.block ≈
          nonidentity_density_reference atol = 1.0e-14 rtol = 1.0e-14
    @test !nonidentity_density_density.diagnostics.source_weight_division_applied_by_helper
    function _raw_pair_terms_from_density_normalized(normalized_terms, weights)
        raw_terms = similar(normalized_terms)
        weight_outer = weights * transpose(weights)
        for term in axes(normalized_terms, 1)
            raw_terms[term, :, :] .= normalized_terms[term, :, :] .* weight_outer
        end
        return raw_terms
    end
    raw_weighted_pair_terms = (
        x = _raw_pair_terms_from_density_normalized(
            density_normalized_pair_terms.x,
            density_source_weights.x,
        ),
        y = _raw_pair_terms_from_density_normalized(
            density_normalized_pair_terms.y,
            density_source_weights.y,
        ),
        z = _raw_pair_terms_from_density_normalized(
            density_normalized_pair_terms.z,
            density_source_weights.z,
        ),
    )
    raw_weighted_density_density =
        CCPM._product_doside_source_box_raw_weighted_density_density_interaction_block(
            product_unit,
            product_unit;
            term_coefficients = density_density_term_coefficients,
            raw_axis_pair_factor_terms = raw_weighted_pair_terms,
            axis_weights = density_source_weights,
        )
    @test raw_weighted_density_density.path ==
          :product_doside_source_box_raw_weighted_density_density_interaction
    @test raw_weighted_density_density.block ≈
          product_density_density.block atol = 1.0e-14 rtol = 1.0e-14
    @test raw_weighted_density_density.normalized_axis_pair_factor_terms.x ≈
          density_normalized_pair_terms.x atol = 1.0e-14 rtol = 1.0e-14
    @test raw_weighted_density_density.diagnostics.pair_factor_normalization ==
          :raw_weighted
    @test raw_weighted_density_density.diagnostics.raw_weighted_pair_factors
    @test !raw_weighted_density_density.diagnostics.density_normalized_pair_factors
    @test raw_weighted_density_density.diagnostics.density_normalized_pair_factors_generated
    @test raw_weighted_density_density.diagnostics.source_weight_division_owner ==
          :source_box_raw_weights
    @test raw_weighted_density_density.diagnostics.source_weight_division_applied_by_helper
    @test raw_weighted_density_density.diagnostics.source_weight_division_shape ==
          :axis_pair_weight_outer
    @test !raw_weighted_density_density.diagnostics.retained_weight_division_allowed
    @test !raw_weighted_density_density.diagnostics.retained_pqs_weight_division_allowed
    @test !raw_weighted_density_density.diagnostics.ida_weight_division_allowed
    @test !raw_weighted_density_density.diagnostics.mwg_ida_semantics_changed
    @test !raw_weighted_density_density.diagnostics.packet_adoption
    @test !raw_weighted_density_density.diagnostics.qwhamiltonian_consumes
    @test !raw_weighted_density_density.diagnostics.numerical_reference_fallback
    @test_throws ArgumentError CCPM._product_doside_source_box_raw_weighted_density_density_interaction_block(
        product_unit,
        product_unit;
        term_coefficients = density_density_term_coefficients,
        raw_axis_pair_factor_terms = raw_weighted_pair_terms,
        axis_weights = (x = [1.25, 0.0], y = [0.9, 1.1], z = [2.0]),
    )
    @test_throws ArgumentError CCPM._product_doside_source_box_density_density_interaction_block(
        product_unit,
        product_unit;
        term_coefficients = density_density_term_coefficients,
        axis_pair_factor_terms = density_normalized_pair_terms,
        axis_weights = density_source_weights,
        pair_factor_normalization = :raw_weighted,
    )
    @test_throws ArgumentError CCPM._product_doside_source_box_density_density_interaction_block(
        product_unit,
        product_unit;
        term_coefficients = density_density_term_coefficients,
        axis_pair_factor_terms = density_normalized_pair_terms,
        axis_weights = (x = [1.25, -0.75], y = [0.9, 1.1], z = [2.0]),
    )
    @test_throws ArgumentError CCPM._product_doside_source_box_density_density_interaction_block(
        product_unit,
        product_unit;
        term_coefficients = density_density_term_coefficients[1:1],
        axis_pair_factor_terms = density_normalized_pair_terms,
        axis_weights = density_source_weights,
    )
    @test_throws ArgumentError CCPM._product_doside_source_box_shadow_blocks(
        product_unit,
        product_unit,
        physical_axis_metrics;
        terms = (:weights,),
    )
    @test_throws ArgumentError CCPM._product_doside_source_box_reference_block(
        product_unit,
        product_unit,
        physical_axis_metrics;
        term = :weights,
    )
    product_staged_self_blocks = _product_staged_comparison_retained_blocks(
        product_unit,
        product_unit,
        physical_axis_metrics,
    )
    @test product_staged_self_blocks.helper_path == :fill_product_staged_metric_blocks
    @test product_staged_self_blocks.fixture_scope == :private_test_only
    @test size(product_direct_self_overlap) == (4, 4)
    @test product_direct_self_overlap ≈ product_staged_self_blocks.overlap atol = 1.0e-14 rtol = 1.0e-14
    @test factorized_product_retained_overlap.retained_operator_matrix ≈
          product_staged_self_blocks.overlap atol = 1.0e-14 rtol = 1.0e-14
    @test_throws ArgumentError CCPM._product_doside_retained_low_order_block(
        product_unit,
        product_unit,
        physical_axis_metrics;
        term = :axis_index_x,
    )
    @test_throws ArgumentError CCPM._product_doside_retained_low_order_block(
        product_unit,
        product_unit,
        physical_axis_metrics;
        term = :kinetic,
    )
    support_dense_retained_block_guard_unit =
        GaussletBases._CartesianNestedProductStagedByCenterUnit3D(
            :support_dense_retained_block_guard,
            :support_dense,
            1:4,
            copy(product_unit.support_indices),
            copy(product_unit.support_states),
            product_unit.coefficient_matrix,
            product_unit.axes,
            copy(product_unit.axis_function_indices),
            (source = :product_retained_block_guard_test,),
            (support_count = 4, retained_count = 4),
        )
    @test_throws ArgumentError CCPM._product_doside_retained_low_order_block(
        support_dense_retained_block_guard_unit,
        product_unit,
        physical_axis_metrics;
        term = :overlap,
    )
    mismatched_axis_metadata_retained_block_guard_unit =
        GaussletBases._CartesianNestedProductStagedByCenterUnit3D(
            :mismatched_axis_metadata_retained_block_guard,
            :product_doside,
            1:4,
            copy(product_unit.support_indices),
            copy(product_unit.support_states),
            product_unit.coefficient_matrix,
            product_unit.axes,
            product_unit.axis_function_indices[1:3],
            (source = :product_retained_block_guard_test,),
            (support_count = 4, retained_count = 4),
        )
    @test_throws ArgumentError CCPM._product_doside_retained_low_order_block(
        mismatched_axis_metadata_retained_block_guard_unit,
        product_unit,
        physical_axis_metrics;
        term = :overlap,
    )
    bad_retained_block_axis_metrics = (
        x = (
            overlap = Matrix{Float64}(I, 1, 1),
            position = Matrix{Float64}(I, 1, 1),
            weights = [1.0],
            centers = [0.0],
            source = :bad_product_retained_block_axis_metric_fixture,
        ),
        y = physical_axis_metrics.y,
        z = physical_axis_metrics.z,
    )
    @test_throws ArgumentError CCPM._product_doside_retained_low_order_block(
        product_unit,
        product_unit,
        bad_retained_block_axis_metrics;
        term = :overlap,
    )
    @test_throws ArgumentError CCP._cartesian_factorized_product_doside_raw_low_order_operator_packet(
        overlap_omitting_pair;
        term = :overlap,
        axis_metrics = physical_axis_metrics,
    )
    @test_throws ArgumentError CCP._cartesian_factorized_product_doside_raw_low_order_operator_packet(
        cross_product_pair;
        term = :overlap,
        axis_metrics = physical_axis_metrics,
    )
    @test_throws ArgumentError CCP._cartesian_factorized_product_doside_raw_low_order_operator_packet(
        pqs_resolved_pair;
        term = :overlap,
        axis_metrics = physical_axis_metrics,
    )
    physical_position_packets = Dict{Symbol,Any}()
    generic_physical_position_packets = Dict{Symbol,Any}()
    nonidentity_physical_position_packets = Dict{Symbol,Any}()
    cross_physical_position_packets = Dict{Symbol,Any}()
    factorized_position_packets = Dict{Symbol,Any}()
    factorized_cross_position_packets = Dict{Symbol,Any}()
    for term in (:position_x, :position_y, :position_z)
        packet = CCP._cartesian_physical_raw_low_order_operator_packet(
            product_self_pair;
            term,
            axis_metrics = physical_axis_metrics,
        )
        generic_packet = CCP._cartesian_physical_raw_low_order_operator_packet(
            generic_product_self_pair;
            term,
            axis_metrics = physical_axis_metrics,
        )
        nonidentity_packet = CCP._cartesian_physical_raw_low_order_operator_packet(
            nonidentity_product_self_pair;
            term,
            axis_metrics = physical_axis_metrics,
        )
        cross_packet = CCP._cartesian_physical_raw_low_order_operator_packet(
            cross_product_pair;
            term,
            axis_metrics = physical_axis_metrics,
        )
        factorized_packet =
            CCP._cartesian_factorized_product_doside_raw_low_order_operator_packet(
                product_self_pair;
                term,
                axis_metrics = physical_axis_metrics,
            )
        factorized_cross_packet =
            CCP._cartesian_factorized_product_doside_raw_low_order_operator_packet(
                cross_product_pair;
                term,
                axis_metrics = physical_axis_metrics,
            )
        physical_position_packets[term] = packet
        generic_physical_position_packets[term] = generic_packet
        nonidentity_physical_position_packets[term] = nonidentity_packet
        cross_physical_position_packets[term] = cross_packet
        factorized_position_packets[term] = factorized_packet
        factorized_cross_position_packets[term] = factorized_cross_packet
        expected = Matrix{Float64}(getproperty(expected_physical_positions, term))
        reference = support_local_physical_position_reference(term)
        cross_reference = support_local_physical_position_reference(
            term,
            cross_product_pair.left_raw_source.provenance.unit.support_states,
            cross_product_pair.right_raw_source.provenance.unit.support_states,
        )
        @test packet.left_source_id == :identity_product_slab_source
        @test packet.right_source_id == :identity_product_slab_source
        @test packet.operator_kind == :low_order_metric
        @test packet.term == term
        @test packet.source_dimensions == (4, 4)
        @test packet.raw_operator_matrix ≈ reference atol = 1.0e-14 rtol = 1.0e-14
        @test packet.raw_operator_matrix ≈ expected atol = 1.0e-14 rtol = 1.0e-14
        @test generic_packet.left_source_id == :generic_product_slab_source
        @test generic_packet.right_source_id == :generic_product_slab_source
        @test generic_packet.operator_kind == :low_order_metric
        @test generic_packet.term == term
        @test generic_packet.source_dimensions == (4, 4)
        @test generic_packet.raw_operator_matrix ≈ reference atol = 1.0e-14 rtol = 1.0e-14
        @test generic_packet.raw_operator_matrix ≈ expected atol = 1.0e-14 rtol = 1.0e-14
        @test nonidentity_packet.left_source_id == :nonidentity_product_slab_source
        @test nonidentity_packet.right_source_id == :nonidentity_product_slab_source
        @test nonidentity_packet.operator_kind == :low_order_metric
        @test nonidentity_packet.term == term
        @test nonidentity_packet.source_dimensions == (4, 4)
        @test nonidentity_packet.raw_operator_matrix ≈ reference atol = 1.0e-14 rtol = 1.0e-14
        @test nonidentity_packet.raw_operator_matrix ≈ expected atol = 1.0e-14 rtol = 1.0e-14
        @test cross_packet.left_source_id == :identity_product_slab_source
        @test cross_packet.right_source_id == :nonidentity_product_slab_source
        @test cross_packet.operator_kind == :low_order_metric
        @test cross_packet.term == term
        @test cross_packet.source_dimensions == (4, 4)
        @test cross_packet.raw_operator_matrix ≈ cross_reference atol = 1.0e-14 rtol = 1.0e-14
        @test cross_packet.raw_operator_matrix ≈ expected atol = 1.0e-14 rtol = 1.0e-14
        @test factorized_packet.left_source_id == :identity_product_slab_source
        @test factorized_packet.right_source_id == :identity_product_slab_source
        @test factorized_packet.operator_kind == :low_order_metric
        @test factorized_packet.term == term
        @test factorized_packet.source_dimensions == (4, 4)
        @test factorized_packet.raw_operator_matrix ≈ reference atol = 1.0e-14 rtol = 1.0e-14
        @test factorized_packet.raw_operator_matrix ≈ packet.raw_operator_matrix atol = 1.0e-14 rtol = 1.0e-14
        @test factorized_cross_packet.left_source_id == :identity_product_slab_source
        @test factorized_cross_packet.right_source_id == :nonidentity_product_slab_source
        @test factorized_cross_packet.operator_kind == :low_order_metric
        @test factorized_cross_packet.term == term
        @test factorized_cross_packet.source_dimensions == (4, 4)
        @test factorized_cross_packet.raw_operator_matrix ≈ cross_reference atol = 1.0e-14 rtol = 1.0e-14
        @test factorized_cross_packet.raw_operator_matrix ≈ cross_packet.raw_operator_matrix atol = 1.0e-14 rtol = 1.0e-14
        @test packet.raw_operator_matrix != expected_product_axis_index_x
        @test generic_packet.raw_operator_matrix != expected_product_axis_index_x
        @test nonidentity_packet.raw_operator_matrix != expected_product_axis_index_x
        @test cross_packet.raw_operator_matrix != expected_product_axis_index_x
        @test factorized_packet.raw_operator_matrix != expected_product_axis_index_x
        @test factorized_cross_packet.raw_operator_matrix != expected_product_axis_index_x
        @test packet.diagnostics.raw_reference ==
              Symbol(:product_doside_physical_, term)
        @test generic_packet.diagnostics.raw_reference ==
              Symbol(:product_doside_physical_, term)
        @test nonidentity_packet.diagnostics.raw_reference ==
              Symbol(:product_doside_physical_, term)
        @test cross_packet.diagnostics.raw_reference ==
              Symbol(:product_doside_physical_, term)
        @test factorized_packet.diagnostics.raw_reference ==
              Symbol(:product_doside_factorized_, term)
        @test factorized_cross_packet.diagnostics.raw_reference ==
              Symbol(:product_doside_factorized_, term)
        @test packet.diagnostics.raw_reference_error == 0.0
        @test generic_packet.diagnostics.raw_reference_error == 0.0
        @test nonidentity_packet.diagnostics.raw_reference_error == 0.0
        @test cross_packet.diagnostics.raw_reference_error == 0.0
        @test factorized_packet.diagnostics.raw_reference_error == 0.0
        @test factorized_cross_packet.diagnostics.raw_reference_error == 0.0
        @test packet.diagnostics.physical_position_operator
        @test generic_packet.diagnostics.physical_position_operator
        @test nonidentity_packet.diagnostics.physical_position_operator
        @test cross_packet.diagnostics.physical_position_operator
        @test factorized_packet.diagnostics.physical_position_operator
        @test factorized_cross_packet.diagnostics.physical_position_operator
        @test !packet.diagnostics.axis_index_diagnostic
        @test !generic_packet.diagnostics.axis_index_diagnostic
        @test !nonidentity_packet.diagnostics.axis_index_diagnostic
        @test !cross_packet.diagnostics.axis_index_diagnostic
        @test packet.diagnostics.position_axis ==
              (term == :position_x ? :x : term == :position_y ? :y : :z)
        @test generic_packet.diagnostics.position_axis ==
              (term == :position_x ? :x : term == :position_y ? :y : :z)
        @test nonidentity_packet.diagnostics.position_axis ==
              (term == :position_x ? :x : term == :position_y ? :y : :z)
        @test cross_packet.diagnostics.position_axis ==
              (term == :position_x ? :x : term == :position_y ? :y : :z)
        @test factorized_packet.diagnostics.position_axis ==
              (term == :position_x ? :x : term == :position_y ? :y : :z)
        @test factorized_cross_packet.diagnostics.position_axis ==
              (term == :position_x ? :x : term == :position_y ? :y : :z)
        @test packet.diagnostics.axis_metric_sources.x ==
              :synthetic_noninteger_axis_metric_fixture
        @test generic_packet.diagnostics.axis_metric_sources.x ==
              :synthetic_noninteger_axis_metric_fixture
        @test nonidentity_packet.diagnostics.axis_metric_sources.x ==
              :synthetic_noninteger_axis_metric_fixture
        @test cross_packet.diagnostics.axis_metric_sources.x ==
              :synthetic_noninteger_axis_metric_fixture
        @test factorized_packet.diagnostics.axis_metric_sources.x ==
              :synthetic_noninteger_axis_metric_fixture
        @test factorized_cross_packet.diagnostics.axis_metric_sources.x ==
              :synthetic_noninteger_axis_metric_fixture
        @test packet.diagnostics.raw_basis_scope == :raw_product_source_rows
        @test generic_packet.diagnostics.raw_basis_scope == :raw_product_source_rows
        @test nonidentity_packet.diagnostics.raw_basis_scope == :raw_product_source_rows
        @test cross_packet.diagnostics.raw_basis_scope == :raw_product_source_rows
        @test factorized_packet.diagnostics.raw_basis_scope == :raw_product_source_rows
        @test factorized_cross_packet.diagnostics.raw_basis_scope ==
              :raw_product_source_rows
        @test factorized_packet.diagnostics.factorized_axis_path_used
        @test factorized_cross_packet.diagnostics.factorized_axis_path_used
        @test factorized_packet.diagnostics.axis_factor_scope ==
              :staged_local_axis_intervals
        @test factorized_cross_packet.diagnostics.axis_factor_scope ==
              :staged_local_axis_intervals
        @test !factorized_packet.diagnostics.support_row_reference_used
        @test !factorized_cross_packet.diagnostics.support_row_reference_used
        @test !packet.diagnostics.cross_pair
        @test !generic_packet.diagnostics.cross_pair
        @test !nonidentity_packet.diagnostics.cross_pair
        @test cross_packet.diagnostics.cross_pair
        @test !factorized_packet.diagnostics.cross_pair
        @test factorized_cross_packet.diagnostics.cross_pair
        @test cross_packet.diagnostics.left_source_dimension == 4
        @test cross_packet.diagnostics.right_source_dimension == 4
        @test factorized_cross_packet.diagnostics.left_source_dimension == 4
        @test factorized_cross_packet.diagnostics.right_source_dimension == 4
        @test !packet.diagnostics.dense_parent_matrix_used
        @test !generic_packet.diagnostics.dense_parent_matrix_used
        @test !nonidentity_packet.diagnostics.dense_parent_matrix_used
        @test !cross_packet.diagnostics.dense_parent_matrix_used
        @test !factorized_packet.diagnostics.dense_parent_matrix_used
        @test !factorized_cross_packet.diagnostics.dense_parent_matrix_used
        @test packet.diagnostics.fixture_only
        @test generic_packet.diagnostics.fixture_only
        @test nonidentity_packet.diagnostics.fixture_only
        @test cross_packet.diagnostics.fixture_only
        @test factorized_packet.diagnostics.fixture_only
        @test factorized_cross_packet.diagnostics.fixture_only
        @test !packet.diagnostics.production_supported
        @test !generic_packet.diagnostics.production_supported
        @test !nonidentity_packet.diagnostics.production_supported
        @test !cross_packet.diagnostics.production_supported
        @test !factorized_packet.diagnostics.production_supported
        @test !factorized_cross_packet.diagnostics.production_supported
        @test !packet.diagnostics.metric_execution_changed
        @test !generic_packet.diagnostics.metric_execution_changed
        @test !nonidentity_packet.diagnostics.metric_execution_changed
        @test !cross_packet.diagnostics.metric_execution_changed
        @test !factorized_packet.diagnostics.metric_execution_changed
        @test !factorized_cross_packet.diagnostics.metric_execution_changed
        @test !packet.diagnostics.qwhamiltonian_consumes
        @test !generic_packet.diagnostics.qwhamiltonian_consumes
        @test !nonidentity_packet.diagnostics.qwhamiltonian_consumes
        @test !cross_packet.diagnostics.qwhamiltonian_consumes
        @test !factorized_packet.diagnostics.qwhamiltonian_consumes
        @test !factorized_cross_packet.diagnostics.qwhamiltonian_consumes
        @test !packet.diagnostics.public_default_consumes
        @test !generic_packet.diagnostics.public_default_consumes
        @test !nonidentity_packet.diagnostics.public_default_consumes
        @test !cross_packet.diagnostics.public_default_consumes
        @test !factorized_packet.diagnostics.public_default_consumes
        @test !factorized_cross_packet.diagnostics.public_default_consumes
        @test !packet.diagnostics.backend_policy_changed
        @test !generic_packet.diagnostics.backend_policy_changed
        @test !nonidentity_packet.diagnostics.backend_policy_changed
        @test !cross_packet.diagnostics.backend_policy_changed
        @test !factorized_packet.diagnostics.backend_policy_changed
        @test !factorized_cross_packet.diagnostics.backend_policy_changed
        @test !packet.diagnostics.quadrature_policy_changed
        @test !generic_packet.diagnostics.quadrature_policy_changed
        @test !nonidentity_packet.diagnostics.quadrature_policy_changed
        @test !cross_packet.diagnostics.quadrature_policy_changed
        @test !factorized_packet.diagnostics.quadrature_policy_changed
        @test !factorized_cross_packet.diagnostics.quadrature_policy_changed
        @test !packet.diagnostics.cr2_science_status_changed
        @test !generic_packet.diagnostics.cr2_science_status_changed
        @test !nonidentity_packet.diagnostics.cr2_science_status_changed
        @test !cross_packet.diagnostics.cr2_science_status_changed
        @test !factorized_packet.diagnostics.cr2_science_status_changed
        @test !factorized_cross_packet.diagnostics.cr2_science_status_changed
    end
    @test_throws ArgumentError CCP._cartesian_physical_raw_low_order_operator_packet(
        pqs_resolved_pair;
        term = :position_x,
        axis_metrics = physical_axis_metrics,
    )
    @test_throws ArgumentError CCP._cartesian_physical_raw_low_order_operator_packet(
        position_omitting_pair;
        term = :position_x,
        axis_metrics = physical_axis_metrics,
    )
    @test_throws ArgumentError CCP._cartesian_factorized_product_doside_raw_low_order_operator_packet(
        position_omitting_pair;
        term = :position_x,
        axis_metrics = physical_axis_metrics,
    )
    @test_throws ArgumentError CCP._cartesian_physical_raw_low_order_operator_packet(
        cross_position_omitting_pair;
        term = :position_x,
        axis_metrics = physical_axis_metrics,
    )
    @test_throws ArgumentError CCP._cartesian_factorized_product_doside_raw_low_order_operator_packet(
        cross_position_omitting_pair;
        term = :position_x,
        axis_metrics = physical_axis_metrics,
    )
    @test_throws ArgumentError CCP._cartesian_physical_raw_low_order_operator_packet(
        product_self_pair;
        term = :axis_index_x,
        axis_metrics = physical_axis_metrics,
    )
    @test_throws ArgumentError CCP._cartesian_factorized_product_doside_raw_low_order_operator_packet(
        product_self_pair;
        term = :axis_index_x,
        axis_metrics = physical_axis_metrics,
    )
    @test_throws ArgumentError CCP._cartesian_physical_raw_low_order_operator_packet(
        cross_product_pair;
        term = :axis_index_x,
        axis_metrics = physical_axis_metrics,
    )
    @test_throws ArgumentError CCP._cartesian_raw_low_order_operator_packet(
        cross_product_pair;
        term = :axis_index_x,
    )
    product_retained_overlap = CCP._cartesian_retained_low_order_operator_block(
        product_raw_overlap_packet,
        product_source_transform.retained_transform,
    )
    @test product_retained_overlap.left_source_id == :identity_product_slab_source
    @test product_retained_overlap.right_source_id == :identity_product_slab_source
    @test product_retained_overlap.operator_kind == :low_order_metric
    @test product_retained_overlap.term == :overlap
    @test product_retained_overlap.retained_dimensions == (4, 4)
    @test product_retained_overlap.retained_operator_matrix ≈
          Matrix{Float64}(I, 4, 4) atol = 1.0e-14 rtol = 1.0e-14
    @test product_retained_overlap.diagnostics.left_transform_kind ==
          :product_axis_transform
    @test product_retained_overlap.diagnostics.right_transform_kind ==
          :product_axis_transform
    @test product_retained_overlap.diagnostics.left_transform_materialized
    @test product_retained_overlap.diagnostics.right_transform_materialized
    @test product_retained_overlap.diagnostics.retained_transform_applied
    @test product_retained_overlap.diagnostics.retained_operator_block_built
    @test !product_retained_overlap.diagnostics.all_pair_matrices_built
    @test !product_retained_overlap.diagnostics.metric_execution_changed
    @test !product_retained_overlap.diagnostics.qwhamiltonian_consumes
    @test !product_retained_overlap.diagnostics.public_default_consumes
    @test !product_retained_overlap.diagnostics.backend_policy_changed
    @test !product_retained_overlap.diagnostics.quadrature_policy_changed
    @test !product_retained_overlap.diagnostics.cr2_science_status_changed
    @test !product_retained_overlap.diagnostics.ida_weight_division_allowed
    product_retained_axis_index_x = CCP._cartesian_retained_low_order_operator_block(
        product_raw_axis_index_x_packet,
        product_source_transform.retained_transform,
    )
    @test product_retained_axis_index_x.left_source_id == :identity_product_slab_source
    @test product_retained_axis_index_x.right_source_id == :identity_product_slab_source
    @test product_retained_axis_index_x.operator_kind == :low_order_metric
    @test product_retained_axis_index_x.term == :axis_index_x
    @test product_retained_axis_index_x.retained_dimensions == (4, 4)
    @test product_retained_axis_index_x.retained_operator_matrix ≈
          expected_product_axis_index_x atol = 1.0e-14 rtol = 1.0e-14
    @test product_retained_axis_index_x.diagnostics.left_transform_kind ==
          :product_axis_transform
    @test product_retained_axis_index_x.diagnostics.right_transform_kind ==
          :product_axis_transform
    @test product_retained_axis_index_x.diagnostics.retained_transform_applied
    @test product_retained_axis_index_x.diagnostics.retained_operator_block_built
    @test product_retained_axis_index_x.diagnostics.retained_reference ==
          :explicit_materialized_transform
    @test !product_retained_axis_index_x.diagnostics.all_pair_matrices_built
    @test !product_retained_axis_index_x.diagnostics.metric_execution_changed
    @test !product_retained_axis_index_x.diagnostics.qwhamiltonian_consumes
    @test !product_retained_axis_index_x.diagnostics.public_default_consumes
    @test !product_retained_axis_index_x.diagnostics.backend_policy_changed
    @test !product_retained_axis_index_x.diagnostics.quadrature_policy_changed
    @test !product_retained_axis_index_x.diagnostics.cr2_science_status_changed
    @test !product_retained_axis_index_x.diagnostics.ida_weight_division_allowed
    for term in (:position_x, :position_y, :position_z)
        retained = CCP._cartesian_retained_low_order_operator_block(
            physical_position_packets[term],
            product_source_transform.retained_transform,
        )
        generic_retained = CCP._cartesian_retained_low_order_operator_block(
            generic_physical_position_packets[term],
            generic_product_source_transform.retained_transform,
        )
        nonidentity_retained = CCP._cartesian_retained_low_order_operator_block(
            nonidentity_physical_position_packets[term],
            nonidentity_product_source_transform.retained_transform,
        )
        cross_retained = CCP._cartesian_retained_low_order_operator_block(
            cross_physical_position_packets[term],
            product_source_transform.retained_transform,
            nonidentity_product_source_transform.retained_transform,
        )
        factorized_retained = CCP._cartesian_retained_low_order_operator_block(
            factorized_position_packets[term],
            product_source_transform.retained_transform,
        )
        factorized_cross_retained = CCP._cartesian_retained_low_order_operator_block(
            factorized_cross_position_packets[term],
            product_source_transform.retained_transform,
            nonidentity_product_source_transform.retained_transform,
        )
        expected = Matrix{Float64}(getproperty(expected_physical_positions, term))
        reference = support_local_physical_position_reference(term)
        nonidentity_reference = transpose(nonidentity_transform) * reference * nonidentity_transform
        cross_reference = support_local_physical_position_reference(
            term,
            cross_product_pair.left_raw_source.provenance.unit.support_states,
            cross_product_pair.right_raw_source.provenance.unit.support_states,
        )
        cross_retained_reference = cross_reference * nonidentity_transform
        @test retained.left_source_id == :identity_product_slab_source
        @test retained.right_source_id == :identity_product_slab_source
        @test retained.operator_kind == :low_order_metric
        @test retained.term == term
        @test retained.retained_dimensions == (4, 4)
        @test retained.retained_operator_matrix ≈ reference atol = 1.0e-14 rtol = 1.0e-14
        @test retained.retained_operator_matrix ≈ expected atol = 1.0e-14 rtol = 1.0e-14
        @test generic_retained.left_source_id == :generic_product_slab_source
        @test generic_retained.right_source_id == :generic_product_slab_source
        @test generic_retained.operator_kind == :low_order_metric
        @test generic_retained.term == term
        @test generic_retained.retained_dimensions == (4, 4)
        @test generic_retained.retained_operator_matrix ≈ reference atol = 1.0e-14 rtol = 1.0e-14
        @test generic_retained.retained_operator_matrix ≈ expected atol = 1.0e-14 rtol = 1.0e-14
        @test nonidentity_retained.left_source_id == :nonidentity_product_slab_source
        @test nonidentity_retained.right_source_id == :nonidentity_product_slab_source
        @test nonidentity_retained.operator_kind == :low_order_metric
        @test nonidentity_retained.term == term
        @test nonidentity_retained.retained_dimensions == (2, 2)
        @test nonidentity_retained.retained_operator_matrix ≈ nonidentity_reference atol = 1.0e-14 rtol = 1.0e-14
        @test cross_retained.left_source_id == :identity_product_slab_source
        @test cross_retained.right_source_id == :nonidentity_product_slab_source
        @test cross_retained.operator_kind == :low_order_metric
        @test cross_retained.term == term
        @test cross_retained.retained_dimensions == (4, 2)
        @test cross_retained.retained_operator_matrix ≈ cross_retained_reference atol = 1.0e-14 rtol = 1.0e-14
        @test factorized_retained.left_source_id == :identity_product_slab_source
        @test factorized_retained.right_source_id == :identity_product_slab_source
        @test factorized_retained.operator_kind == :low_order_metric
        @test factorized_retained.term == term
        @test factorized_retained.retained_dimensions == (4, 4)
        @test factorized_retained.retained_operator_matrix ≈ reference atol = 1.0e-14 rtol = 1.0e-14
        product_direct_self_block = CCPM._product_doside_retained_low_order_block(
            product_unit,
            product_unit,
            physical_axis_metrics;
            term,
        )
        @test size(product_direct_self_block) == (4, 4)
        @test product_direct_self_block ≈
              _product_staged_comparison_block_for_term(product_staged_self_blocks, term) atol = 1.0e-14 rtol = 1.0e-14
        @test factorized_retained.retained_operator_matrix ≈
              _product_staged_comparison_block_for_term(product_staged_self_blocks, term) atol = 1.0e-14 rtol = 1.0e-14
        @test factorized_cross_retained.left_source_id ==
              :identity_product_slab_source
        @test factorized_cross_retained.right_source_id ==
              :nonidentity_product_slab_source
        @test factorized_cross_retained.operator_kind == :low_order_metric
        @test factorized_cross_retained.term == term
        @test factorized_cross_retained.retained_dimensions == (4, 2)
        @test factorized_cross_retained.retained_operator_matrix ≈ cross_retained_reference atol = 1.0e-14 rtol = 1.0e-14
        @test retained.diagnostics.left_transform_kind == :product_axis_transform
        @test generic_retained.diagnostics.left_transform_kind == :product_axis_transform
        @test nonidentity_retained.diagnostics.left_transform_kind == :product_axis_transform
        @test cross_retained.diagnostics.left_transform_kind == :product_axis_transform
        @test factorized_retained.diagnostics.left_transform_kind ==
              :product_axis_transform
        @test factorized_cross_retained.diagnostics.left_transform_kind ==
              :product_axis_transform
        @test retained.diagnostics.right_transform_kind == :product_axis_transform
        @test generic_retained.diagnostics.right_transform_kind == :product_axis_transform
        @test nonidentity_retained.diagnostics.right_transform_kind == :product_axis_transform
        @test cross_retained.diagnostics.right_transform_kind == :product_axis_transform
        @test factorized_retained.diagnostics.right_transform_kind ==
              :product_axis_transform
        @test factorized_cross_retained.diagnostics.right_transform_kind ==
              :product_axis_transform
        @test retained.diagnostics.left_transform_materialized
        @test generic_retained.diagnostics.left_transform_materialized
        @test nonidentity_retained.diagnostics.left_transform_materialized
        @test cross_retained.diagnostics.left_transform_materialized
        @test factorized_retained.diagnostics.left_transform_materialized
        @test factorized_cross_retained.diagnostics.left_transform_materialized
        @test retained.diagnostics.right_transform_materialized
        @test generic_retained.diagnostics.right_transform_materialized
        @test nonidentity_retained.diagnostics.right_transform_materialized
        @test cross_retained.diagnostics.right_transform_materialized
        @test factorized_retained.diagnostics.right_transform_materialized
        @test factorized_cross_retained.diagnostics.right_transform_materialized
        @test retained.diagnostics.retained_transform_applied
        @test generic_retained.diagnostics.retained_transform_applied
        @test nonidentity_retained.diagnostics.retained_transform_applied
        @test cross_retained.diagnostics.retained_transform_applied
        @test factorized_retained.diagnostics.retained_transform_applied
        @test factorized_cross_retained.diagnostics.retained_transform_applied
        @test retained.diagnostics.retained_operator_block_built
        @test generic_retained.diagnostics.retained_operator_block_built
        @test nonidentity_retained.diagnostics.retained_operator_block_built
        @test cross_retained.diagnostics.retained_operator_block_built
        @test factorized_retained.diagnostics.retained_operator_block_built
        @test factorized_cross_retained.diagnostics.retained_operator_block_built
        @test retained.diagnostics.retained_reference == :explicit_materialized_transform
        @test generic_retained.diagnostics.retained_reference == :explicit_materialized_transform
        @test nonidentity_retained.diagnostics.retained_reference == :explicit_materialized_transform
        @test cross_retained.diagnostics.retained_reference == :explicit_materialized_transform
        @test factorized_retained.diagnostics.retained_reference ==
              :explicit_materialized_transform
        @test factorized_cross_retained.diagnostics.retained_reference ==
              :explicit_materialized_transform
        @test retained.diagnostics.raw_operator_matrix_source ==
              :private_physical_raw_low_order_operator_packet
        @test generic_retained.diagnostics.raw_operator_matrix_source ==
              :private_physical_raw_low_order_operator_packet
        @test nonidentity_retained.diagnostics.raw_operator_matrix_source ==
              :private_physical_raw_low_order_operator_packet
        @test cross_retained.diagnostics.raw_operator_matrix_source ==
              :private_physical_raw_low_order_operator_packet
        @test factorized_retained.diagnostics.raw_operator_matrix_source ==
              :private_factorized_product_doside_raw_low_order_operator_packet
        @test factorized_cross_retained.diagnostics.raw_operator_matrix_source ==
              :private_factorized_product_doside_raw_low_order_operator_packet
        @test !retained.diagnostics.all_pair_matrices_built
        @test !generic_retained.diagnostics.all_pair_matrices_built
        @test !nonidentity_retained.diagnostics.all_pair_matrices_built
        @test !cross_retained.diagnostics.all_pair_matrices_built
        @test !factorized_retained.diagnostics.all_pair_matrices_built
        @test !factorized_cross_retained.diagnostics.all_pair_matrices_built
        @test !retained.diagnostics.metric_execution_changed
        @test !generic_retained.diagnostics.metric_execution_changed
        @test !nonidentity_retained.diagnostics.metric_execution_changed
        @test !cross_retained.diagnostics.metric_execution_changed
        @test !factorized_retained.diagnostics.metric_execution_changed
        @test !factorized_cross_retained.diagnostics.metric_execution_changed
        @test !retained.diagnostics.qwhamiltonian_consumes
        @test !generic_retained.diagnostics.qwhamiltonian_consumes
        @test !nonidentity_retained.diagnostics.qwhamiltonian_consumes
        @test !cross_retained.diagnostics.qwhamiltonian_consumes
        @test !factorized_retained.diagnostics.qwhamiltonian_consumes
        @test !factorized_cross_retained.diagnostics.qwhamiltonian_consumes
        @test !retained.diagnostics.public_default_consumes
        @test !generic_retained.diagnostics.public_default_consumes
        @test !nonidentity_retained.diagnostics.public_default_consumes
        @test !cross_retained.diagnostics.public_default_consumes
        @test !factorized_retained.diagnostics.public_default_consumes
        @test !factorized_cross_retained.diagnostics.public_default_consumes
        @test !retained.diagnostics.backend_policy_changed
        @test !generic_retained.diagnostics.backend_policy_changed
        @test !nonidentity_retained.diagnostics.backend_policy_changed
        @test !cross_retained.diagnostics.backend_policy_changed
        @test !factorized_retained.diagnostics.backend_policy_changed
        @test !factorized_cross_retained.diagnostics.backend_policy_changed
        @test !retained.diagnostics.quadrature_policy_changed
        @test !generic_retained.diagnostics.quadrature_policy_changed
        @test !nonidentity_retained.diagnostics.quadrature_policy_changed
        @test !cross_retained.diagnostics.quadrature_policy_changed
        @test !factorized_retained.diagnostics.quadrature_policy_changed
        @test !factorized_cross_retained.diagnostics.quadrature_policy_changed
        @test !retained.diagnostics.cr2_science_status_changed
        @test !generic_retained.diagnostics.cr2_science_status_changed
        @test !nonidentity_retained.diagnostics.cr2_science_status_changed
        @test !cross_retained.diagnostics.cr2_science_status_changed
        @test !factorized_retained.diagnostics.cr2_science_status_changed
        @test !factorized_cross_retained.diagnostics.cr2_science_status_changed
        @test !retained.diagnostics.ida_weight_division_allowed
        @test !generic_retained.diagnostics.ida_weight_division_allowed
        @test !nonidentity_retained.diagnostics.ida_weight_division_allowed
        @test !cross_retained.diagnostics.ida_weight_division_allowed
        @test !factorized_retained.diagnostics.ida_weight_division_allowed
        @test !factorized_cross_retained.diagnostics.ida_weight_division_allowed
    end
    distinct_left_transform = [
        1.0 0.0 0.0
        0.0 1.0 0.0
        0.0 0.0 1.0
        0.5 -0.25 0.75
    ]
    distinct_right_transform = [
        0.6 0.0
        0.0 0.7
        0.8 0.1
        0.0 -0.6
    ]
    distinct_left_unit = GaussletBases._CartesianNestedProductStagedByCenterUnit3D(
        :distinct_left_product_slab,
        :product_doside,
        1:3,
        collect(1:4),
        NTuple{3,Int}[(1, 1, 1), (1, 2, 1), (2, 1, 1), (2, 2, 1)],
        distinct_left_transform,
        (
            GaussletBases._nested_product_staged_active_axis(1:2, identity_axis),
            GaussletBases._nested_product_staged_active_axis(1:2, identity_axis),
            GaussletBases._nested_product_staged_fixed_axis(1),
        ),
        NTuple{3,Int}[(1, 1, 1), (2, 1, 1), (2, 2, 1)],
        (source = :distinct_left_raw_product_source_test_fixture,),
        (support_count = 4, retained_count = 3),
    )
    distinct_right_unit = GaussletBases._CartesianNestedProductStagedByCenterUnit3D(
        :distinct_right_product_slab,
        :product_doside,
        1:2,
        collect(3:6),
        NTuple{3,Int}[(2, 1, 1), (2, 2, 1), (3, 1, 1), (3, 2, 1)],
        distinct_right_transform,
        (
            GaussletBases._nested_product_staged_active_axis(2:3, identity_axis),
            GaussletBases._nested_product_staged_active_axis(1:2, identity_axis),
            GaussletBases._nested_product_staged_fixed_axis(1),
        ),
        NTuple{3,Int}[(1, 1, 1), (2, 2, 1)],
        (source = :distinct_right_raw_product_source_test_fixture,),
        (support_count = 4, retained_count = 2),
    )
    distinct_left_source_transform = CCP._cartesian_raw_product_source_retained_transform(
        distinct_left_unit;
        source_id = :distinct_left_product_slab_source,
        parent_dims = (3, 2, 1),
    )
    distinct_right_source_transform = CCP._cartesian_raw_product_source_retained_transform(
        distinct_right_unit;
        source_id = :distinct_right_product_slab_source,
        parent_dims = (3, 2, 1),
    )
    @test distinct_left_source_transform.raw_source.support_indices == collect(1:4)
    @test distinct_right_source_transform.raw_source.support_indices == collect(3:6)
    @test distinct_left_source_transform.raw_source.support_indices !=
          distinct_right_source_transform.raw_source.support_indices
    @test distinct_left_source_transform.retained_transform.retained_dimension == 3
    @test distinct_right_source_transform.retained_transform.retained_dimension == 2
    distinct_pair_plan = CCP._cartesian_raw_product_source_pair_plan(
        (distinct_left_source_transform, distinct_right_source_transform);
        operator_kind = :low_order_metric,
        supported_terms = (:position_x, :position_y, :position_z),
        source = :distinct_support_product_cross_pair_plan_test,
    )
    distinct_cross_pair = only(
        pair for pair in CCP._cartesian_resolved_raw_product_source_pairs(
            distinct_pair_plan,
        )
        if pair.pair_key ==
           (:distinct_left_product_slab_source, :distinct_right_product_slab_source)
    )
    @test distinct_cross_pair.left_raw_source.parent_dims == (3, 2, 1)
    @test distinct_cross_pair.right_raw_source.parent_dims == (3, 2, 1)
    @test distinct_cross_pair.left_raw_source.source_dimension == 4
    @test distinct_cross_pair.right_raw_source.source_dimension == 4
    @test distinct_cross_pair.left_retained_transform.retained_dimension == 3
    @test distinct_cross_pair.right_retained_transform.retained_dimension == 2
    @test !distinct_cross_pair.diagnostics.pqs_factored_transform_present
    distinct_overlap_pair_plan = CCP._cartesian_raw_product_source_pair_plan(
        (distinct_left_source_transform, distinct_right_source_transform);
        operator_kind = :low_order_metric,
        supported_terms = (:overlap,),
        source = :distinct_support_product_cross_overlap_pair_plan_test,
    )
    distinct_overlap_pair = only(
        pair for pair in CCP._cartesian_resolved_raw_product_source_pairs(
            distinct_overlap_pair_plan,
        )
        if pair.pair_key ==
           (:distinct_left_product_slab_source, :distinct_right_product_slab_source)
    )
    distinct_position_omitting_pair_plan = CCP._cartesian_raw_product_source_pair_plan(
        (distinct_left_source_transform, distinct_right_source_transform);
        operator_kind = :low_order_metric,
        supported_terms = (:position_y, :position_z),
        source = :distinct_support_position_omitting_pair_plan_test,
    )
    distinct_position_omitting_pair = only(
        pair for pair in CCP._cartesian_resolved_raw_product_source_pairs(
            distinct_position_omitting_pair_plan,
        )
        if pair.pair_key ==
           (:distinct_left_product_slab_source, :distinct_right_product_slab_source)
    )
    distinct_axis_metrics = (
        x = (
            overlap = [
                1.0 0.1 0.0
                0.1 1.0 0.2
                0.0 0.2 1.0
            ],
            position = [
                0.0 0.05 0.0
                0.05 1.0 0.3
                0.0 0.3 2.0
            ],
            weights = [1.0, 1.0, 1.0],
            centers = [0.0, 1.0, 2.0],
            source = :distinct_support_axis_metric_fixture,
        ),
        y = (
            overlap = [
                1.0 0.05
                0.05 1.0
            ],
            position = [
                -0.5 0.02
                0.02 0.5
            ],
            weights = [1.0, 1.0],
            centers = [-0.5, 0.5],
            source = :distinct_support_axis_metric_fixture,
        ),
        z = (
            overlap = Matrix{Float64}(I, 1, 1),
            position = Matrix(Diagonal([2.75])),
            weights = [1.0],
            centers = [2.75],
            source = :distinct_support_axis_metric_fixture,
        ),
    )
    function distinct_support_physical_position_reference(term::Symbol)
        axis_matrices = term == :position_x ? (
            distinct_axis_metrics.x.position,
            distinct_axis_metrics.y.overlap,
            distinct_axis_metrics.z.overlap,
        ) : term == :position_y ? (
            distinct_axis_metrics.x.overlap,
            distinct_axis_metrics.y.position,
            distinct_axis_metrics.z.overlap,
        ) : term == :position_z ? (
            distinct_axis_metrics.x.overlap,
            distinct_axis_metrics.y.overlap,
            distinct_axis_metrics.z.position,
        ) : throw(ArgumentError("unsupported distinct-support physical term"))
        left_states = distinct_cross_pair.left_raw_source.provenance.unit.support_states
        right_states = distinct_cross_pair.right_raw_source.provenance.unit.support_states
        reference = zeros(Float64, length(left_states), length(right_states))
        for col in eachindex(right_states)
            jx, jy, jz = right_states[col]
            for row in eachindex(left_states)
                ix, iy, iz = left_states[row]
                reference[row, col] =
                    axis_matrices[1][ix, jx] *
                    axis_matrices[2][iy, jy] *
                    axis_matrices[3][iz, jz]
            end
        end
        return reference
    end
    function distinct_support_overlap_reference()
        left_states = distinct_overlap_pair.left_raw_source.provenance.unit.support_states
        right_states = distinct_overlap_pair.right_raw_source.provenance.unit.support_states
        axis_matrices = (
            distinct_axis_metrics.x.overlap,
            distinct_axis_metrics.y.overlap,
            distinct_axis_metrics.z.overlap,
        )
        reference = zeros(Float64, length(left_states), length(right_states))
        for col in eachindex(right_states)
            jx, jy, jz = right_states[col]
            for row in eachindex(left_states)
                ix, iy, iz = left_states[row]
                reference[row, col] =
                    axis_matrices[1][ix, jx] *
                    axis_matrices[2][iy, jy] *
                    axis_matrices[3][iz, jz]
            end
        end
        return reference
    end
    distinct_overlap_packet =
        CCP._cartesian_factorized_product_doside_raw_low_order_operator_packet(
            distinct_overlap_pair;
            term = :overlap,
            axis_metrics = distinct_axis_metrics,
        )
    distinct_overlap_reference = distinct_support_overlap_reference()
    distinct_overlap_retained = CCP._cartesian_retained_low_order_operator_block(
        distinct_overlap_packet,
        distinct_left_source_transform.retained_transform,
        distinct_right_source_transform.retained_transform,
    )
    distinct_overlap_retained_reference =
        transpose(distinct_left_transform) *
        distinct_overlap_reference *
        distinct_right_transform
    @test distinct_overlap_packet.left_source_id == :distinct_left_product_slab_source
    @test distinct_overlap_packet.right_source_id == :distinct_right_product_slab_source
    @test distinct_overlap_packet.source_dimensions == (4, 4)
    @test distinct_overlap_packet.raw_operator_matrix ≈ distinct_overlap_reference atol = 1.0e-14 rtol = 1.0e-14
    @test distinct_overlap_packet.raw_operator_matrix != Matrix{Float64}(I, 4, 4)
    @test distinct_overlap_packet.diagnostics.cross_pair
    @test distinct_overlap_packet.diagnostics.factorized_axis_path_used
    @test distinct_overlap_packet.diagnostics.axis_factor_scope ==
          :staged_local_axis_intervals
    @test distinct_overlap_packet.diagnostics.raw_basis_scope ==
          :raw_product_source_rows
    @test !distinct_overlap_packet.diagnostics.support_row_reference_used
    @test distinct_overlap_packet.diagnostics.fixture_only
    @test !distinct_overlap_packet.diagnostics.production_supported
    @test !distinct_overlap_packet.diagnostics.dense_parent_matrix_used
    @test !distinct_overlap_packet.diagnostics.metric_execution_changed
    @test !distinct_overlap_packet.diagnostics.qwhamiltonian_consumes
    @test !distinct_overlap_packet.diagnostics.backend_policy_changed
    @test !distinct_overlap_packet.diagnostics.quadrature_policy_changed
    @test !distinct_overlap_packet.diagnostics.cr2_science_status_changed
    @test distinct_overlap_retained.retained_dimensions == (3, 2)
    @test distinct_overlap_retained.retained_operator_matrix ≈ distinct_overlap_retained_reference atol = 1.0e-14 rtol = 1.0e-14
    bad_distinct_axis_metrics = (
        x = (
            overlap = Matrix{Float64}(I, 2, 2),
            source = :bad_axis_metric_dimension_fixture,
        ),
        y = distinct_axis_metrics.y,
        z = distinct_axis_metrics.z,
    )
    @test_throws DimensionMismatch CCP._cartesian_factorized_product_doside_raw_low_order_operator_packet(
        distinct_overlap_pair;
        term = :overlap,
        axis_metrics = bad_distinct_axis_metrics,
    )
    bad_distinct_position_axis_metrics = (
        x = (
            overlap = Matrix{Float64}(I, 3, 3),
            position = Matrix{Float64}(I, 2, 2),
            source = :bad_position_axis_metric_dimension_fixture,
        ),
        y = distinct_axis_metrics.y,
        z = distinct_axis_metrics.z,
    )
    @test_throws DimensionMismatch CCP._cartesian_factorized_product_doside_raw_low_order_operator_packet(
        distinct_cross_pair;
        term = :position_x,
        axis_metrics = bad_distinct_position_axis_metrics,
    )
    outside_interval_unit = GaussletBases._CartesianNestedProductStagedByCenterUnit3D(
        :outside_interval_product_slab,
        :product_doside,
        1:4,
        [1, 2, 5, 4],
        NTuple{3,Int}[(1, 1, 1), (1, 2, 1), (3, 1, 1), (2, 2, 1)],
        Matrix{Float64}(I, 4, 4),
        (
            GaussletBases._nested_product_staged_active_axis(1:2, identity_axis),
            GaussletBases._nested_product_staged_active_axis(1:2, identity_axis),
            GaussletBases._nested_product_staged_fixed_axis(1),
        ),
        GaussletBases._nested_product_axis_function_indices(3, 1, 2, 2, 2),
        (source = :outside_interval_raw_product_source_test_fixture,),
        (support_count = 4, retained_count = 4),
    )
    outside_interval_source_transform =
        CCP._cartesian_raw_product_source_retained_transform(
            outside_interval_unit;
            source_id = :outside_interval_product_slab_source,
            parent_dims = (3, 2, 1),
        )
    outside_interval_pair_plan = CCP._cartesian_raw_product_source_pair_plan(
        (outside_interval_source_transform, distinct_right_source_transform);
        operator_kind = :low_order_metric,
        supported_terms = (:overlap,),
        source = :outside_interval_product_overlap_pair_plan_test,
    )
    outside_interval_pair = only(
        pair for pair in CCP._cartesian_resolved_raw_product_source_pairs(
            outside_interval_pair_plan,
        )
        if pair.pair_key ==
           (:distinct_right_product_slab_source, :outside_interval_product_slab_source)
    )
    @test_throws ArgumentError CCP._cartesian_factorized_product_doside_raw_low_order_operator_packet(
        outside_interval_pair;
        term = :overlap,
        axis_metrics = distinct_axis_metrics,
    )
    outside_interval_position_pair_plan = CCP._cartesian_raw_product_source_pair_plan(
        (outside_interval_source_transform, distinct_right_source_transform);
        operator_kind = :low_order_metric,
        supported_terms = (:position_x,),
        source = :outside_interval_product_position_pair_plan_test,
    )
    outside_interval_position_pair = only(
        pair for pair in CCP._cartesian_resolved_raw_product_source_pairs(
            outside_interval_position_pair_plan,
        )
        if pair.pair_key ==
           (:distinct_right_product_slab_source, :outside_interval_product_slab_source)
    )
    @test_throws ArgumentError CCP._cartesian_factorized_product_doside_raw_low_order_operator_packet(
        outside_interval_position_pair;
        term = :position_x,
        axis_metrics = distinct_axis_metrics,
    )
    for term in (:position_x, :position_y, :position_z)
        distinct_packet = CCP._cartesian_physical_raw_low_order_operator_packet(
            distinct_cross_pair;
            term,
            axis_metrics = distinct_axis_metrics,
        )
        distinct_factorized_packet =
            CCP._cartesian_factorized_product_doside_raw_low_order_operator_packet(
                distinct_cross_pair;
                term,
                axis_metrics = distinct_axis_metrics,
            )
        distinct_reference = distinct_support_physical_position_reference(term)
        distinct_retained = CCP._cartesian_retained_low_order_operator_block(
            distinct_packet,
            distinct_left_source_transform.retained_transform,
            distinct_right_source_transform.retained_transform,
        )
        distinct_factorized_retained = CCP._cartesian_retained_low_order_operator_block(
            distinct_factorized_packet,
            distinct_left_source_transform.retained_transform,
            distinct_right_source_transform.retained_transform,
        )
        distinct_retained_reference =
            transpose(distinct_left_transform) *
            distinct_reference *
            distinct_right_transform
        @test distinct_packet.left_source_id == :distinct_left_product_slab_source
        @test distinct_packet.right_source_id == :distinct_right_product_slab_source
        @test distinct_packet.source_dimensions == (4, 4)
        @test distinct_packet.raw_operator_matrix ≈ distinct_reference atol = 1.0e-14 rtol = 1.0e-14
        @test distinct_packet.raw_operator_matrix !=
              Matrix{Float64}(getproperty(expected_physical_positions, term))
        @test distinct_factorized_packet.left_source_id ==
              :distinct_left_product_slab_source
        @test distinct_factorized_packet.right_source_id ==
              :distinct_right_product_slab_source
        @test distinct_factorized_packet.source_dimensions == (4, 4)
        @test distinct_factorized_packet.raw_operator_matrix ≈ distinct_reference atol = 1.0e-14 rtol = 1.0e-14
        @test distinct_factorized_packet.raw_operator_matrix ≈ distinct_packet.raw_operator_matrix atol = 1.0e-14 rtol = 1.0e-14
        @test distinct_factorized_packet.diagnostics.factorized_axis_path_used
        @test distinct_factorized_packet.diagnostics.physical_position_operator
        @test distinct_factorized_packet.diagnostics.position_axis ==
              (term == :position_x ? :x : term == :position_y ? :y : :z)
        @test !distinct_factorized_packet.diagnostics.support_row_reference_used
        @test distinct_factorized_packet.diagnostics.raw_basis_scope ==
              :raw_product_source_rows
        @test distinct_factorized_packet.diagnostics.fixture_only
        @test !distinct_factorized_packet.diagnostics.production_supported
        @test !distinct_factorized_packet.diagnostics.dense_parent_matrix_used
        @test !distinct_factorized_packet.diagnostics.metric_execution_changed
        @test distinct_packet.diagnostics.cross_pair
        @test distinct_packet.diagnostics.left_source_dimension == 4
        @test distinct_packet.diagnostics.right_source_dimension == 4
        @test distinct_packet.diagnostics.raw_basis_scope == :raw_product_source_rows
        @test distinct_packet.diagnostics.fixture_only
        @test !distinct_packet.diagnostics.production_supported
        @test !distinct_packet.diagnostics.dense_parent_matrix_used
        @test !distinct_packet.diagnostics.metric_execution_changed
        @test !distinct_packet.diagnostics.qwhamiltonian_consumes
        @test !distinct_packet.diagnostics.backend_policy_changed
        @test !distinct_packet.diagnostics.quadrature_policy_changed
        @test !distinct_packet.diagnostics.cr2_science_status_changed
        @test distinct_retained.retained_dimensions == (3, 2)
        @test distinct_retained.retained_operator_matrix ≈ distinct_retained_reference atol = 1.0e-14 rtol = 1.0e-14
        @test distinct_factorized_retained.retained_dimensions == (3, 2)
        @test distinct_factorized_retained.retained_operator_matrix ≈ distinct_retained_reference atol = 1.0e-14 rtol = 1.0e-14
        @test distinct_factorized_retained.diagnostics.raw_operator_matrix_source ==
              :private_factorized_product_doside_raw_low_order_operator_packet
        @test distinct_retained.diagnostics.retained_transform_applied
        @test distinct_retained.diagnostics.retained_operator_block_built
        @test !distinct_retained.diagnostics.all_pair_matrices_built
        @test !distinct_retained.diagnostics.metric_execution_changed
        @test !distinct_retained.diagnostics.qwhamiltonian_consumes
        @test !distinct_retained.diagnostics.public_default_consumes
        @test !distinct_retained.diagnostics.backend_policy_changed
        @test !distinct_retained.diagnostics.quadrature_policy_changed
        @test !distinct_retained.diagnostics.cr2_science_status_changed
    end
    consistent_left_axes = (
        GaussletBases._nested_product_staged_active_axis(
            1:2,
            [
                1.0 0.2
                0.1 0.9
            ],
        ),
        GaussletBases._nested_product_staged_active_axis(
            1:2,
            [
                0.8 0.3
                -0.2 1.1
            ],
        ),
        GaussletBases._nested_product_staged_fixed_axis(1),
    )
    consistent_left_support_states =
        NTuple{3,Int}[(1, 1, 1), (1, 2, 1), (2, 1, 1), (2, 2, 1)]
    consistent_left_axis_function_indices =
        GaussletBases._nested_product_axis_function_indices(3, 1, 2, 2, 2)
    consistent_left_unit = GaussletBases._CartesianNestedProductStagedByCenterUnit3D(
        :consistent_left_product_slab,
        :product_doside,
        1:4,
        collect(1:4),
        consistent_left_support_states,
        _product_staged_comparison_coefficient_matrix(
            consistent_left_support_states,
            consistent_left_axes,
            consistent_left_axis_function_indices,
        ),
        consistent_left_axes,
        consistent_left_axis_function_indices,
        (source = :consistent_left_product_staged_metric_comparison_fixture,),
        (support_count = 4, retained_count = 4),
    )
    consistent_right_axes = (
        GaussletBases._nested_product_staged_active_axis(
            2:3,
            [
                0.6 -0.1
                0.4 0.9
            ],
        ),
        GaussletBases._nested_product_staged_active_axis(
            1:2,
            reshape([0.7, -0.6], 2, 1),
        ),
        GaussletBases._nested_product_staged_fixed_axis(1),
    )
    consistent_right_support_states =
        NTuple{3,Int}[(2, 1, 1), (2, 2, 1), (3, 1, 1), (3, 2, 1)]
    consistent_right_axis_function_indices =
        GaussletBases._nested_product_axis_function_indices(3, 1, 2, 2, 1)
    consistent_right_unit = GaussletBases._CartesianNestedProductStagedByCenterUnit3D(
        :consistent_right_product_slab,
        :product_doside,
        1:2,
        collect(3:6),
        consistent_right_support_states,
        _product_staged_comparison_coefficient_matrix(
            consistent_right_support_states,
            consistent_right_axes,
            consistent_right_axis_function_indices,
        ),
        consistent_right_axes,
        consistent_right_axis_function_indices,
        (source = :consistent_right_product_staged_metric_comparison_fixture,),
        (support_count = 4, retained_count = 2),
    )
    consistent_left_source_transform = CCP._cartesian_raw_product_source_retained_transform(
        consistent_left_unit;
        source_id = :consistent_left_product_slab_source,
        parent_dims = (3, 2, 1),
    )
    consistent_right_source_transform = CCP._cartesian_raw_product_source_retained_transform(
        consistent_right_unit;
        source_id = :consistent_right_product_slab_source,
        parent_dims = (3, 2, 1),
    )
    consistent_pair_plan = CCP._cartesian_raw_product_source_pair_plan(
        (consistent_left_source_transform, consistent_right_source_transform);
        operator_kind = :low_order_metric,
        supported_terms = (:overlap, :position_x, :position_y, :position_z),
        source = :consistent_product_staged_metric_comparison_pair_plan_test,
    )
    consistent_cross_pair = only(
        pair for pair in CCP._cartesian_resolved_raw_product_source_pairs(
            consistent_pair_plan,
        )
        if pair.pair_key ==
           (:consistent_left_product_slab_source, :consistent_right_product_slab_source)
    )
    consistent_product_staged_blocks = _product_staged_comparison_retained_blocks(
        consistent_left_unit,
        consistent_right_unit,
        distinct_axis_metrics,
    )
    @test consistent_product_staged_blocks.helper_path ==
          :fill_product_staged_metric_blocks
    @test consistent_product_staged_blocks.fixture_scope == :private_test_only
    consistent_left_helper_weights, consistent_left_helper_first_moments =
        CCPM._product_doside_retained_linear_vectors(
            consistent_left_unit,
            distinct_axis_metrics,
        )
    consistent_left_staged_weights, consistent_left_staged_first_moments =
        CCPM._staged_unit_linear_vectors(
            consistent_left_unit,
            distinct_axis_metrics,
        )
    consistent_left_support_dense_linear_reference_unit =
        GaussletBases._CartesianNestedProductStagedByCenterUnit3D(
            :consistent_left_support_dense_linear_vector_reference,
            :support_dense,
            consistent_left_unit.column_range,
            copy(consistent_left_unit.support_indices),
            copy(consistent_left_unit.support_states),
            consistent_left_unit.coefficient_matrix,
            consistent_left_unit.axes,
            copy(consistent_left_unit.axis_function_indices),
            (source = :consistent_left_product_staged_linear_vector_reference_fixture,),
            (support_count = 4, retained_count = 4),
        )
    consistent_left_reference_weights, consistent_left_reference_first_moments =
        CCPM._staged_unit_linear_vectors(
            consistent_left_support_dense_linear_reference_unit,
            distinct_axis_metrics,
        )
    @test consistent_left_helper_weights ≈ consistent_left_staged_weights atol = 1.0e-14 rtol = 1.0e-14
    @test consistent_left_helper_first_moments ≈ consistent_left_staged_first_moments atol = 1.0e-14 rtol = 1.0e-14
    @test consistent_left_helper_weights ≈ consistent_left_reference_weights atol = 1.0e-14 rtol = 1.0e-14
    @test consistent_left_helper_first_moments ≈ consistent_left_reference_first_moments atol = 1.0e-14 rtol = 1.0e-14
    @test length(consistent_left_helper_weights) == length(consistent_left_unit.column_range)
    @test size(consistent_left_helper_first_moments) ==
          (length(consistent_left_unit.column_range), 3)
    kinetic_axis_ops = (
        x = (
            overlap = distinct_axis_metrics.x.overlap,
            kinetic = [
                1.2 -0.4 0.1
                -0.4 1.6 -0.35
                0.1 -0.35 2.1
            ],
        ),
        y = (
            overlap = distinct_axis_metrics.y.overlap,
            kinetic = [
                0.9 -0.25
                -0.25 1.3
            ],
        ),
        z = (
            overlap = distinct_axis_metrics.z.overlap,
            kinetic = [0.7;;],
        ),
    )
    kinetic_axis_factor_terms = (
        (kinetic_axis_ops.x.kinetic, kinetic_axis_ops.y.overlap, kinetic_axis_ops.z.overlap),
        (kinetic_axis_ops.x.overlap, kinetic_axis_ops.y.kinetic, kinetic_axis_ops.z.overlap),
        (kinetic_axis_ops.x.overlap, kinetic_axis_ops.y.overlap, kinetic_axis_ops.z.kinetic),
    )
    function product_doside_support_local_separable_sum_reference(
        left_unit,
        right_unit,
        axis_factor_terms,
    )
        left_entries = CCPM._staged_unit_entries(left_unit)
        right_entries = CCPM._staged_unit_entries(right_unit)
        reference = zeros(Float64, length(left_unit.column_range), length(right_unit.column_range))
        for factors in axis_factor_terms
            reference .+= CCPM._contract_pair_block(
                left_entries,
                right_entries,
                factors[1],
                factors[2],
                factors[3],
            )
        end
        return reference
    end
    consistent_left_kinetic =
        CCPM._product_doside_retained_kinetic_block(
            consistent_left_unit,
            consistent_left_unit,
            kinetic_axis_ops,
        )
    consistent_left_kinetic_reference =
        product_doside_support_local_separable_sum_reference(
            consistent_left_unit,
            consistent_left_unit,
            kinetic_axis_factor_terms,
        )
    consistent_cross_kinetic =
        CCPM._product_doside_retained_kinetic_block(
            consistent_left_unit,
            consistent_right_unit,
            kinetic_axis_ops,
        )
    consistent_cross_separable_sum =
        CCPM._product_doside_retained_separable_sum_block(
            consistent_left_unit,
            consistent_right_unit,
            kinetic_axis_factor_terms,
        )
    consistent_cross_kinetic_reference =
        product_doside_support_local_separable_sum_reference(
            consistent_left_unit,
            consistent_right_unit,
            kinetic_axis_factor_terms,
        )
    @test size(consistent_left_kinetic) == (4, 4)
    @test size(consistent_cross_kinetic) == (4, 2)
    @test consistent_left_kinetic ≈ consistent_left_kinetic_reference atol = 1.0e-14 rtol = 1.0e-14
    @test consistent_left_kinetic ≈ transpose(consistent_left_kinetic) atol = 1.0e-14 rtol = 1.0e-14
    @test consistent_cross_kinetic ≈ consistent_cross_separable_sum atol = 1.0e-14 rtol = 1.0e-14
    @test consistent_cross_kinetic ≈ consistent_cross_kinetic_reference atol = 1.0e-14 rtol = 1.0e-14
    for term in (:overlap, :position_x, :position_y, :position_z)
        packet = CCP._cartesian_factorized_product_doside_raw_low_order_operator_packet(
            consistent_cross_pair;
            term,
            axis_metrics = distinct_axis_metrics,
        )
        retained = CCP._cartesian_retained_low_order_operator_block(
            packet,
            consistent_left_source_transform.retained_transform,
            consistent_right_source_transform.retained_transform,
        )
        direct_retained = CCPM._product_doside_retained_low_order_block(
            consistent_left_unit,
            consistent_right_unit,
            distinct_axis_metrics;
            term,
        )
        @test retained.retained_dimensions == (4, 2)
        @test size(direct_retained) == (4, 2)
        @test direct_retained ≈
              _product_staged_comparison_block_for_term(consistent_product_staged_blocks, term) atol = 1.0e-14 rtol = 1.0e-14
        @test retained.retained_operator_matrix ≈
              _product_staged_comparison_block_for_term(consistent_product_staged_blocks, term) atol = 1.0e-14 rtol = 1.0e-14
        @test packet.diagnostics.factorized_axis_path_used
        @test !packet.diagnostics.support_row_reference_used
        @test packet.diagnostics.fixture_only
        @test !packet.diagnostics.production_supported
        @test !packet.diagnostics.metric_execution_changed
        @test !retained.diagnostics.metric_execution_changed
        @test !retained.diagnostics.qwhamiltonian_consumes
        @test !retained.diagnostics.public_default_consumes
        @test !retained.diagnostics.backend_policy_changed
        @test !retained.diagnostics.quadrature_policy_changed
        @test !retained.diagnostics.cr2_science_status_changed
        @test !retained.diagnostics.ida_weight_division_allowed
    end
    @test_throws ArgumentError CCP._cartesian_physical_raw_low_order_operator_packet(
        distinct_position_omitting_pair;
        term = :position_x,
        axis_metrics = distinct_axis_metrics,
    )
    @test_throws ArgumentError CCP._cartesian_retained_low_order_operator_block(
        pqs_raw_overlap_packet,
        pqs_retained_transform,
    )
    mixed_cross_pair = only(
        pair for pair in mixed_resolved_pairs
        if pair.pair_key == (mixed_source_ids[1], mixed_source_ids[2])
    )
    @test Set((
        mixed_cross_pair.left_raw_source.source_dimension,
        mixed_cross_pair.right_raw_source.source_dimension,
    )) == Set((4, 5 * 5 * 5))
    @test Set((
        mixed_cross_pair.left_retained_transform.retained_dimension,
        mixed_cross_pair.right_retained_transform.retained_dimension,
    )) == Set((4, 98))
    @test mixed_cross_pair.diagnostics.pqs_factored_transform_present
    mixed_plan_audit = CCP._cartesian_raw_product_source_pair_plan_audit(mixed_pair_plan)
    @test length(mixed_plan_audit.resolved_pairs) == 3
    @test mixed_plan_audit.diagnostics.every_pair_resolves_raw_sources
    @test mixed_plan_audit.diagnostics.every_pair_resolves_retained_transforms
    @test mixed_plan_audit.diagnostics.every_pair_upper_triangular
    @test mixed_plan_audit.diagnostics.every_pair_placeholder_only
    @test mixed_plan_audit.diagnostics.raw_weight_roles_explicit
    @test mixed_plan_audit.diagnostics.retained_weight_roles_explicit
    @test !mixed_plan_audit.diagnostics.retained_ida_weight_division_allowed
    @test !mixed_plan_audit.diagnostics.raw_operator_matrices_built
    @test !mixed_plan_audit.diagnostics.retained_operator_blocks_built
    @test !mixed_plan_audit.diagnostics.metric_execution_changed
    @test !mixed_plan_audit.diagnostics.qwhamiltonian_consumes
    @test !mixed_plan_audit.diagnostics.public_default_consumes
    @test !mixed_plan_audit.diagnostics.backend_policy_changed
    @test !mixed_plan_audit.diagnostics.quadrature_policy_changed
    @test !mixed_plan_audit.diagnostics.cr2_science_status_changed
    pqs_product_policy = CCPM._pqs_product_mixed_block_policy()
    @test pqs_product_policy.pair_kind == :pqs_product_mixed
    @test pqs_product_policy.optimized_metric_path == :unsupported_pqs_product_optimized
    @test !pqs_product_policy.optimized_supported
    @test !pqs_product_policy.support_local_reference_allowed
    @test pqs_product_policy.support_local_reference_path ==
          :not_available_without_explicit_request
    @test pqs_product_policy.fixture_only
    @test !pqs_product_policy.production_supported
    @test pqs_product_policy.reason == :pqs_product_optimized_metric_not_implemented
    pqs_product_reference_policy = CCPM._pqs_product_mixed_block_policy(
        explicit_reference_requested = true,
    )
    @test pqs_product_reference_policy.support_local_reference_allowed
    @test pqs_product_reference_policy.support_local_reference_path ==
          :support_local_reference_explicit_only
    fake_product_resolved = CCP._CartesianResolvedContractionPayload3D(
        :product_staged_metric_contraction,
        true,
        :product_doside,
        1:1,
        Int[1],
        NTuple{3,Int}[(1, 1, 1)],
        cubic_pqs_payload,
        (),
        (
            source = :pqs_product_policy_test,
            rule_family = :product_owned_unit,
            rule_kind = :product_doside,
            metric_capability = :product_staged_metric_contraction,
            linear_vector_path = :product_staged_axis_projection,
            block_role = :product,
            unsupported = false,
            prototype = false,
        ),
        (source = :pqs_product_policy_test,),
    )
    pqs_product_dispatch = CCPM._metric_dispatch_plan_from_resolved_payloads(
        [cubic_pqs_resolved, fake_product_resolved],
    )
    @test !pqs_product_dispatch.plan_supported
    @test pqs_product_dispatch.pqs_product_unsupported_block_count == 1
    @test :unsupported_pqs_product_optimized in
          [path.path for path in pqs_product_dispatch.block_paths]
    pqs_product_error = try
        CCPM._resolved_payload_low_order_metric_block(
            cubic_pqs_resolved,
            fake_product_resolved,
            cubic_metrics,
        )
        nothing
    catch err
        err
    end
    @test pqs_product_error isa ArgumentError
    @test occursin(
        "PQS/product optimized metric blocks are explicitly unsupported",
        sprint(showerror, pqs_product_error),
    )

    rectangular_bundles = GaussletBases._CartesianNestedAxisBundles3D(
        bundle5,
        bundle5,
        bundle7,
    )
    rectangular_current = (1:5, 1:5, 1:7)
    rectangular_inner = (2:4, 2:4, 2:6)
    rectangular_expected = setdiff(
        GaussletBases._nested_box_support_indices(rectangular_current..., (5, 5, 7)),
        GaussletBases._nested_box_support_indices(rectangular_inner..., (5, 5, 7)),
    )
    sort!(rectangular_expected)
    rectangular = GaussletBases._nested_projected_q_shell_layer(
        rectangular_bundles,
        rectangular_current,
        rectangular_inner;
        bond_axis = :z,
        q = 5,
        L = 7,
        term_coefficients,
    )

    @test rectangular.support_indices == rectangular_expected
    @test rectangular.diagnostics.boundary_support_count == 130
    @test rectangular.diagnostics.full_block_column_count == 175
    @test rectangular.diagnostics.boundary_comx_product_mode_count == 130
    @test rectangular.diagnostics.boundary_comx_product_mode_selection_rule ==
          :any_axis_mode_index_first_or_last
    @test rectangular.diagnostics.retained_count == 130
    @test rectangular.diagnostics.rank_count == 130
    @test rectangular.diagnostics.rank_drop_count == 0
    @test rectangular.diagnostics.duplicate_count == 0
    @test rectangular.diagnostics.missing_count == 0
    @test rectangular.diagnostics.outside_count == 0
    @test rectangular.diagnostics.overlap_error < 1.0e-8
    @test rectangular.diagnostics.packet_overlap_error < 1.0e-8
    @test !rectangular.diagnostics.endcap_panel_stitching
    @test size(rectangular.coefficient_matrix) == (175, 130)
    @test !_pqs_one_hot_selector_columns(
        rectangular.coefficient_matrix[rectangular.support_indices, :],
    )
    @test all(isfinite, rectangular.packet.weights)
    @test minimum(rectangular.packet.weights) > -1.0e-10
    rectangular_descriptor =
        GaussletBases._nested_projected_q_shell_staged_unit_descriptor(rectangular)
    @test rectangular_descriptor.kind == :projected_q_shell
    @test rectangular_descriptor.current_box == rectangular_current
    @test rectangular_descriptor.inner_box == rectangular_inner
    @test rectangular_descriptor.bond_axis == :z
    @test rectangular_descriptor.q == 5
    @test rectangular_descriptor.L == 7
    @test rectangular_descriptor.support_count == 130
    @test rectangular_descriptor.mode_count == 130
    @test rectangular_descriptor.retained_count == 130
    rectangular_shared_raw_product_box_plan = GaussletBases._cartesian_raw_product_box_plan(
        rectangular_bundles,
        rectangular_descriptor.axis_intervals,
        (5, 5, 7);
        enforce_symmetric_odd = false,
    )
    @test all(
        mode -> _pqs_boundary_mode_rule_ok(mode, (5, 5, 7)),
        rectangular_descriptor.boundary_mode_indices,
    )
    @test rectangular_descriptor.boundary_column_indices ==
          rectangular.provenance.pqs_staged_unit_descriptor.boundary_column_indices
    @test rectangular_descriptor.cleanup_method == :projected_boundary_symmetric_lowdin
    @test rectangular_descriptor.cleanup_matrix_size == (130, 130)
    @test rectangular_descriptor.cleanup_rank_count == 130
    @test rectangular_descriptor.cleanup_rank_drop_count == 0
    @test rectangular_descriptor.support_local_coefficient_shape == (130, 130)
    @test :product_doside_unit in rectangular_descriptor.non_contracts
    @test :dense_full_parent_fallback in rectangular_descriptor.non_contracts
    @test rectangular_descriptor.diagnostics.metadata_only
    @test !rectangular_descriptor.active_consumption.fixed_block_sidecar_installed
    rectangular_rule = CCP.cartesian_contraction_rule(
        rectangular_descriptor;
        parent_dimension = 5 * 5 * 7,
    )
    @test CCP.contraction_rule_family(rectangular_rule) == :projected_q_shell_boundary_modes
    @test rectangular_rule.support_summary.entry_count == 130
    @test rectangular_rule.support_summary.unique_count == 130
    @test rectangular_rule.support_summary.missing_count == 45
    @test rectangular_rule.local_geometry.current_box == rectangular_current
    @test rectangular_rule.local_geometry.inner_box == rectangular_inner
    @test rectangular_rule.local_geometry.axis_local_coefficient_shapes == ((5, 5), (5, 5), (7, 7))
    @test CCP.contraction_rule_source_dimension(rectangular_rule) == 175
    @test CCP.contraction_rule_retained_dimension(rectangular_rule) == 130
    @test CCP.contraction_rule_transform_rule(rectangular_rule) ==
          :boundary_comx_product_mode_selection
    @test CCP.contraction_rule_cleanup_rule(rectangular_rule) == :shell_projection_lowdin
    @test CCP.contraction_rule_metric_capability(rectangular_rule) ==
          :pqs_low_order_product_metric_prototype
    rectangular_metric_prototype =
        GaussletBases._nested_projected_q_shell_descriptor_metric_prototype(
            rectangular,
            rectangular_bundles,
        )
    @test rectangular_metric_prototype.overlap ≈ rectangular.packet.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test rectangular_metric_prototype.weights ≈ rectangular.packet.weights atol = 1.0e-10 rtol = 1.0e-10
    @test rectangular_metric_prototype.diagnostics.descriptor_kind == :projected_q_shell
    @test !rectangular_metric_prototype.diagnostics.overlap_invariant_applied
    @test rectangular_metric_prototype.diagnostics.overlap_invariant_debug_check
    @test rectangular_metric_prototype.diagnostics.overlap_invariant_error < 1.0e-10
    @test !rectangular_metric_prototype.diagnostics.overlap_is_operator_target
    @test !rectangular_metric_prototype.diagnostics.dense_full_parent_matrix_used
    @test !rectangular_metric_prototype.diagnostics.product_doside_unit
    @test !rectangular_metric_prototype.diagnostics.fixed_block_sidecar_installed
    @test !rectangular_metric_prototype.diagnostics.optimized_sidecar_installed
    rectangular_metrics = _pqs_axis_metrics(rectangular_bundles)
    _check_pqs_raw_product_box_reference(
        CCPM,
        rectangular_descriptor,
        rectangular_metrics;
        expected_source_mode_dims = (5, 5, 7),
        expected_retained_count = 130,
        shared_raw_product_box_plan = rectangular_shared_raw_product_box_plan,
    )
    rectangular_pqs_pqs_source_box_plan =
        CCPM._pqs_raw_product_box_plan(
            rectangular_descriptor,
            rectangular_shared_raw_product_box_plan,
            rectangular_metrics,
        )
    _check_pqs_pqs_source_box_self_blocks(
        CCPM,
        rectangular_pqs_pqs_source_box_plan,
        rectangular_metrics;
        expected_source_mode_dims = (5, 5, 7),
        expected_retained_count = 130,
    )
    @test_throws ArgumentError CCPM._pqs_pqs_source_box_pair_plan(
        cubic_pqs_pqs_source_box_plan,
        rectangular_pqs_pqs_source_box_plan,
        rectangular_metrics,
    )
    rectangular_pqs_product_source_box_plan =
        CCPM._pqs_raw_product_box_plan(
            rectangular_descriptor,
            rectangular_metrics;
            shared_raw_product_box_plan = rectangular_shared_raw_product_box_plan,
        )
    _check_pqs_product_source_box_mixed_block(
        CCPM,
        rectangular_descriptor,
        rectangular_pqs_product_source_box_plan,
        product_unit,
        rectangular_metrics;
        expected_source_mode_dims = (5, 5, 7),
        expected_retained_count = 130,
        shared_raw_product_box_plan = rectangular_shared_raw_product_box_plan,
    )
    _check_pqs_product_source_box_shadow_blocks(
        CCPM,
        rectangular_descriptor,
        rectangular_pqs_product_source_box_plan,
        nonidentity_axis_product_unit,
        rectangular_metrics;
        expected_source_mode_dims = (5, 5, 7),
        expected_retained_count = 130,
        shared_raw_product_box_plan = rectangular_shared_raw_product_box_plan,
    )
    nuclear_axis_layers = (
        x = bundle5.basis,
        y = bundle5.basis,
        z = bundle7.basis,
    )
    nuclear_probe_expansion = CoulombGaussianExpansion(
        [0.65, 0.18],
        [0.35, 0.9];
        del = 0.0,
        s = 1.0,
        c = 1.0,
        maxu = 1.0,
    )
    nuclear_probe_centers = (
        (0.15, -0.2, 0.35),
        (-0.4, 0.1, -0.15),
    )
    nuclear_probe_charges = (2.0, 0.5)
    function _check_source_box_nuclear_attraction_by_center(
        result,
        positive_builder,
        expected_path::Symbol,
    )
        @test result.path == expected_path
        @test result.physical_operator == :electron_nuclear_attraction
        @test length(result.blocks_by_center) == length(nuclear_probe_centers)
        @test result.center_blocks === result.blocks_by_center
        expected_total = zeros(Float64, size(result.total_block))
        for center_index in eachindex(nuclear_probe_centers)
            contribution = result.blocks_by_center[center_index]
            positive = positive_builder(nuclear_probe_centers[center_index])
            @test contribution.center_index == center_index
            @test contribution.center == nuclear_probe_centers[center_index]
            @test contribution.nuclear_charge == nuclear_probe_charges[center_index]
            @test contribution.sign_charge_scale == -nuclear_probe_charges[center_index]
            @test contribution.gaussian_sum_block ≈ positive.block atol = 0.0 rtol = 0.0
            @test contribution.block ≈
                  (-nuclear_probe_charges[center_index]) .* positive.block atol = 1.0e-14 rtol = 1.0e-14
            @test !positive.diagnostics.nuclear_charge_applied
            @test !positive.diagnostics.nuclear_attraction_sign_applied
            @test contribution.positive_gaussian_sum.path == positive.path
            @test contribution.diagnostics.physical_operator ==
                  :electron_nuclear_attraction
            @test contribution.diagnostics.positive_gaussian_sum_component
            @test contribution.diagnostics.nuclear_charge_applied
            @test contribution.diagnostics.nuclear_attraction_sign_applied
            @test contribution.diagnostics.nuclear_charge_sign_applied
            @test contribution.diagnostics.center_contributions_preserved
            @test contribution.diagnostics.counterpoise_center_identity_preserved
            @test !contribution.diagnostics.ecp
            @test !contribution.diagnostics.packet_adoption
            @test !contribution.diagnostics.fixed_block_routing
            @test !contribution.diagnostics.qwhamiltonian_consumes
            @test !contribution.diagnostics.public_default_consumes
            @test !contribution.diagnostics.cr2_science_status_changed
            expected_total .+= contribution.block
        end
        @test result.total_block ≈ expected_total atol = 1.0e-14 rtol = 1.0e-14
        @test result.block ≈ result.total_block atol = 0.0 rtol = 0.0
        @test result.centers == nuclear_probe_centers
        @test result.nuclear_charges == nuclear_probe_charges
        @test result.diagnostics.source_box_first
        @test result.diagnostics.local_gaussian_source_box_terms
        @test result.diagnostics.physical_operator ==
              :electron_nuclear_attraction
        @test result.diagnostics.positive_gaussian_sum_component
        @test result.diagnostics.nuclear_charge_applied
        @test result.diagnostics.nuclear_attraction_sign_applied
        @test result.diagnostics.nuclear_charge_sign_applied
        @test result.diagnostics.center_contributions_preserved
        @test result.diagnostics.counterpoise_center_identity_preserved
        @test result.diagnostics.primary_result == :blocks_by_center
        @test result.diagnostics.total_block_is_derived_convenience
        @test result.diagnostics.block_field_is_derived_total
        @test result.diagnostics.center_count == length(nuclear_probe_centers)
        @test result.diagnostics.positive_gaussian_sum_convention
        @test !result.diagnostics.ecp
        @test !result.diagnostics.ecp_terms_implemented
        @test !result.diagnostics.electron_electron_terms_implemented
        @test !result.diagnostics.mwg_interaction_implemented
        @test !result.diagnostics.retained_pqs_weight_division_allowed
        @test !result.diagnostics.ida_weight_division_allowed
        @test !result.diagnostics.numerical_reference_fallback
        @test !result.diagnostics.shell_row_algorithm
        @test !result.diagnostics.packet_adoption
        @test !result.diagnostics.fixed_block_routing
        @test !result.diagnostics.qwhamiltonian_consumes
        @test !result.diagnostics.public_default_consumes
        @test !result.diagnostics.cr2_science_status_changed
        @test result.diagnostics.output_finite
    end
    product_nuclear_attraction =
        CCPM._product_doside_source_box_nuclear_attraction_by_center(
            product_unit,
            product_unit,
            nuclear_axis_layers,
            nuclear_probe_expansion;
            centers = nuclear_probe_centers,
            nuclear_charges = nuclear_probe_charges,
        )
    _check_source_box_nuclear_attraction_by_center(
        product_nuclear_attraction,
        center -> CCPM._product_doside_source_box_centered_local_gaussian_sum_block(
            product_unit,
            product_unit,
            nuclear_axis_layers,
            nuclear_probe_expansion;
            center,
        ),
        :product_doside_source_box_nuclear_attraction_by_center,
    )
    pqs_product_nuclear_attraction =
        CCPM._pqs_product_source_box_nuclear_attraction_by_center(
            rectangular_pqs_product_source_box_plan,
            product_unit,
            nuclear_axis_layers,
            nuclear_probe_expansion;
            centers = nuclear_probe_centers,
            nuclear_charges = nuclear_probe_charges,
        )
    _check_source_box_nuclear_attraction_by_center(
        pqs_product_nuclear_attraction,
        center -> CCPM._pqs_product_source_box_centered_local_gaussian_sum_block(
            rectangular_pqs_product_source_box_plan,
            product_unit,
            nuclear_axis_layers,
            nuclear_probe_expansion;
            center,
        ),
        :pqs_product_source_box_nuclear_attraction_by_center,
    )
    pqs_product_density_coefficients = [0.55, 0.21]
    function _pqs_product_density_term_tensor(source_count::Int, axis_scale::Float64)
        terms = Array{Float64,3}(undef, 2, source_count, source_count)
        for term in 1:2, col in 1:source_count, row in 1:source_count
            terms[term, row, col] =
                axis_scale / (1.0 + abs(row - col)) +
                0.03 * term +
                0.007 * row -
                0.005 * col
        end
        return terms
    end
    pqs_product_density_terms = (
        x = _pqs_product_density_term_tensor(5, 0.42),
        y = _pqs_product_density_term_tensor(5, 0.37),
        z = _pqs_product_density_term_tensor(7, 0.31),
    )
    pqs_product_density_weights = (
        x = [1.0 + 0.04 * index for index in 1:5],
        y = [0.9 + 0.03 * index for index in 1:5],
        z = [1.2 + 0.02 * index for index in 1:7],
    )
    pqs_product_density =
        CCPM._pqs_product_source_box_density_density_interaction_block(
            rectangular_pqs_product_source_box_plan,
            nonidentity_axis_product_unit;
            term_coefficients = pqs_product_density_coefficients,
            axis_pair_factor_terms = pqs_product_density_terms,
            axis_weights = pqs_product_density_weights,
        )
    pqs_product_density_reference =
        _pqs_product_density_density_explicit_reference(
            rectangular_descriptor,
            nonidentity_axis_product_unit,
            pqs_product_density_coefficients,
            pqs_product_density_terms,
        )
    @test pqs_product_density.path ==
          :pqs_product_source_box_density_density_interaction
    @test pqs_product_density.interaction_operator ==
          :electron_electron_density_density
    @test pqs_product_density.block ≈
          pqs_product_density_reference atol = 1.0e-13 rtol = 1.0e-13
    @test pqs_product_density.pair_plan.pair_kind ==
          :pqs_product_source_box_density_density_pair
    @test pqs_product_density.pair_plan.left_source_family ==
          :mode_selected_raw_product_box
    @test pqs_product_density.pair_plan.right_source_family == :product_doside
    @test pqs_product_density.pair_plan.left_source_dimensions == (5, 5, 7)
    @test pqs_product_density.pair_plan.right_source_dimensions == (2, 2, 1)
    @test pqs_product_density.pair_plan.left_retained_count == 130
    @test pqs_product_density.pair_plan.right_retained_count == 4
    @test pqs_product_density.pair_plan.supported_terms == (:pair_sum,)
    @test pqs_product_density.pair_plan.source_weights.pqs[3] ==
          pqs_product_density_weights.z
    @test pqs_product_density.pair_plan.source_weights.product[1] ==
          pqs_product_density_weights.x[1:2]
    @test pqs_product_density.pair_plan.product_retained_unit_plan.axis_coefficient_matrices[1] ==
          nonidentity_axis_product_unit.axes[1].coefficient_matrix
    @test pqs_product_density.diagnostics.interaction_operator ==
          :electron_electron_density_density
    @test pqs_product_density.diagnostics.output_representation ==
          :two_index_density_density
    @test !pqs_product_density.diagnostics.four_index_galerkin_tensor
    @test pqs_product_density.diagnostics.pqs_representation ==
          :mode_selected_raw_product_box
    @test pqs_product_density.diagnostics.product_representation ==
          :product_doside
    @test pqs_product_density.diagnostics.density_normalized_pair_factors
    @test !pqs_product_density.diagnostics.raw_weighted_pair_factors
    @test pqs_product_density.diagnostics.source_weight_division_owner ==
          :caller_supplied_density_normalized_pair_factors
    @test !pqs_product_density.diagnostics.source_weight_division_applied_by_helper
    @test !pqs_product_density.diagnostics.retained_pqs_weights_used
    @test !pqs_product_density.diagnostics.retained_weight_division_allowed
    @test !pqs_product_density.diagnostics.retained_pqs_weight_division_allowed
    @test !pqs_product_density.diagnostics.ida_weight_division_allowed
    @test !pqs_product_density.diagnostics.shell_projection_used
    @test !pqs_product_density.diagnostics.lowdin_cleanup_used
    @test !pqs_product_density.diagnostics.support_coefficient_matrix_used
    @test !pqs_product_density.diagnostics.support_local_pqs_oracle_used
    @test !pqs_product_density.diagnostics.packet_adoption
    @test !pqs_product_density.diagnostics.qwhamiltonian_consumes
    @test !pqs_product_density.diagnostics.mwg_ida_semantics_changed
    @test !pqs_product_density.diagnostics.numerical_reference_fallback
    @test pqs_product_density.diagnostics.electron_electron_terms_implemented
    @test pqs_product_density.diagnostics.projected_axis_term_dimensions.x == (2, 5, 2)
    @test pqs_product_density.diagnostics.projected_axis_term_dimensions.z == (2, 7, 1)
    pqs_product_raw_weighted_terms = (
        x = _raw_pair_terms_from_density_normalized(
            pqs_product_density_terms.x,
            pqs_product_density_weights.x,
        ),
        y = _raw_pair_terms_from_density_normalized(
            pqs_product_density_terms.y,
            pqs_product_density_weights.y,
        ),
        z = _raw_pair_terms_from_density_normalized(
            pqs_product_density_terms.z,
            pqs_product_density_weights.z,
        ),
    )
    pqs_product_raw_weighted_density =
        CCPM._pqs_product_source_box_raw_weighted_density_density_interaction_block(
            rectangular_pqs_product_source_box_plan,
            nonidentity_axis_product_unit;
            term_coefficients = pqs_product_density_coefficients,
            raw_axis_pair_factor_terms = pqs_product_raw_weighted_terms,
            axis_weights = pqs_product_density_weights,
        )
    @test pqs_product_raw_weighted_density.path ==
          :pqs_product_source_box_raw_weighted_density_density_interaction
    @test pqs_product_raw_weighted_density.block ≈
          pqs_product_density.block atol = 1.0e-13 rtol = 1.0e-13
    @test pqs_product_raw_weighted_density.normalized_axis_pair_factor_terms.x ≈
          pqs_product_density_terms.x atol = 1.0e-14 rtol = 1.0e-14
    @test pqs_product_raw_weighted_density.normalized_axis_pair_factor_terms.z ≈
          pqs_product_density_terms.z atol = 1.0e-14 rtol = 1.0e-14
    @test pqs_product_raw_weighted_density.diagnostics.pair_factor_normalization ==
          :raw_weighted
    @test pqs_product_raw_weighted_density.diagnostics.raw_weighted_pair_factors
    @test !pqs_product_raw_weighted_density.diagnostics.density_normalized_pair_factors
    @test pqs_product_raw_weighted_density.diagnostics.density_normalized_pair_factors_generated
    @test pqs_product_raw_weighted_density.diagnostics.source_weight_division_owner ==
          :source_box_raw_weights
    @test pqs_product_raw_weighted_density.diagnostics.source_weight_division_applied_by_helper
    @test pqs_product_raw_weighted_density.diagnostics.source_weight_division_shape ==
          :axis_pair_weight_outer
    @test !pqs_product_raw_weighted_density.diagnostics.retained_pqs_weights_used
    @test !pqs_product_raw_weighted_density.diagnostics.retained_weight_division_allowed
    @test !pqs_product_raw_weighted_density.diagnostics.retained_pqs_weight_division_allowed
    @test !pqs_product_raw_weighted_density.diagnostics.ida_weight_division_allowed
    @test !pqs_product_raw_weighted_density.diagnostics.packet_adoption
    @test !pqs_product_raw_weighted_density.diagnostics.qwhamiltonian_consumes
    @test !pqs_product_raw_weighted_density.diagnostics.mwg_ida_semantics_changed
    @test !pqs_product_raw_weighted_density.diagnostics.numerical_reference_fallback
    @test_throws ArgumentError CCPM._pqs_product_source_box_density_density_interaction_block(
        rectangular_pqs_product_source_box_plan,
        nonidentity_axis_product_unit;
        term_coefficients = pqs_product_density_coefficients,
        axis_pair_factor_terms = pqs_product_density_terms,
        axis_weights = pqs_product_density_weights,
        pair_factor_normalization = :raw_weighted,
    )
    @test_throws ArgumentError CCPM._pqs_product_source_box_density_density_interaction_block(
        rectangular_pqs_product_source_box_plan,
        nonidentity_axis_product_unit;
        term_coefficients = pqs_product_density_coefficients,
        axis_pair_factor_terms = pqs_product_density_terms,
        axis_weights = (
            x = [1.0, 1.1, 0.0, 1.3, 1.4],
            y = pqs_product_density_weights.y,
            z = pqs_product_density_weights.z,
        ),
    )
    @test_throws ArgumentError CCPM._pqs_product_source_box_density_density_interaction_block(
        rectangular_pqs_product_source_box_plan,
        nonidentity_axis_product_unit;
        term_coefficients = pqs_product_density_coefficients[1:1],
        axis_pair_factor_terms = pqs_product_density_terms,
        axis_weights = pqs_product_density_weights,
    )
    @test_throws ArgumentError CCPM._pqs_product_source_box_raw_weighted_density_density_interaction_block(
        rectangular_pqs_product_source_box_plan,
        nonidentity_axis_product_unit;
        term_coefficients = pqs_product_density_coefficients,
        raw_axis_pair_factor_terms = pqs_product_raw_weighted_terms,
        axis_weights = (
            x = pqs_product_density_weights.x,
            y = [0.93, 0.96, 0.0, 1.02, 1.05],
            z = pqs_product_density_weights.z,
        ),
    )
    pqs_pqs_density =
        CCPM._pqs_pqs_source_box_density_density_interaction_block(
            rectangular_pqs_pqs_source_box_plan,
            rectangular_pqs_pqs_source_box_plan;
            term_coefficients = pqs_product_density_coefficients,
            axis_pair_factor_terms = pqs_product_density_terms,
            axis_weights = pqs_product_density_weights,
        )
    pqs_pqs_density_reference =
        _pqs_pqs_density_density_explicit_reference(
            rectangular_descriptor,
            rectangular_descriptor,
            pqs_product_density_coefficients,
            pqs_product_density_terms,
        )
    @test pqs_pqs_density.path ==
          :pqs_pqs_source_box_density_density_interaction
    @test pqs_pqs_density.interaction_operator ==
          :electron_electron_density_density
    @test pqs_pqs_density.block ≈
          pqs_pqs_density_reference atol = 1.0e-13 rtol = 1.0e-13
    @test pqs_pqs_density.pair_plan.pair_kind ==
          :pqs_pqs_source_box_density_density_pair
    @test pqs_pqs_density.pair_plan.pair_family == :pqs_pqs
    @test pqs_pqs_density.pair_plan.left_source_family ==
          :mode_selected_raw_product_box
    @test pqs_pqs_density.pair_plan.right_source_family ==
          :mode_selected_raw_product_box
    @test pqs_pqs_density.pair_plan.left_source_dimensions == (5, 5, 7)
    @test pqs_pqs_density.pair_plan.right_source_dimensions == (5, 5, 7)
    @test pqs_pqs_density.pair_plan.left_retained_count == 130
    @test pqs_pqs_density.pair_plan.right_retained_count == 130
    @test pqs_pqs_density.pair_plan.supported_terms == (:pair_sum,)
    @test pqs_pqs_density.pair_plan.source_weights.left[3] ==
          pqs_product_density_weights.z
    @test pqs_pqs_density.pair_plan.source_weights.right[1] ==
          pqs_product_density_weights.x
    @test pqs_pqs_density.diagnostics.pair_family == :pqs_pqs
    @test pqs_pqs_density.diagnostics.interaction_operator ==
          :electron_electron_density_density
    @test pqs_pqs_density.diagnostics.output_representation ==
          :two_index_density_density
    @test !pqs_pqs_density.diagnostics.four_index_galerkin_tensor
    @test pqs_pqs_density.diagnostics.pqs_representation ==
          :mode_selected_raw_product_box
    @test pqs_pqs_density.diagnostics.density_normalized_pair_factors
    @test !pqs_pqs_density.diagnostics.raw_weighted_pair_factors
    @test pqs_pqs_density.diagnostics.source_weight_division_owner ==
          :caller_supplied_density_normalized_pair_factors
    @test !pqs_pqs_density.diagnostics.source_weight_division_applied_by_helper
    @test !pqs_pqs_density.diagnostics.retained_pqs_weights_used
    @test !pqs_pqs_density.diagnostics.retained_pqs_weights_positive_checked
    @test !pqs_pqs_density.diagnostics.retained_weight_division_allowed
    @test !pqs_pqs_density.diagnostics.retained_pqs_weight_division_allowed
    @test !pqs_pqs_density.diagnostics.ida_weight_division_allowed
    @test !pqs_pqs_density.diagnostics.shell_projection_used
    @test !pqs_pqs_density.diagnostics.lowdin_cleanup_used
    @test !pqs_pqs_density.diagnostics.support_coefficient_matrix_used
    @test !pqs_pqs_density.diagnostics.support_local_pqs_oracle_used
    @test !pqs_pqs_density.diagnostics.packet_adoption
    @test !pqs_pqs_density.diagnostics.qwhamiltonian_consumes
    @test !pqs_pqs_density.diagnostics.mwg_ida_semantics_changed
    @test !pqs_pqs_density.diagnostics.numerical_reference_fallback
    @test pqs_pqs_density.diagnostics.electron_electron_terms_implemented
    @test pqs_pqs_density.diagnostics.projected_axis_term_dimensions.x == (2, 5, 5)
    @test pqs_pqs_density.diagnostics.projected_axis_term_dimensions.z == (2, 7, 7)
    pqs_pqs_raw_weighted_density =
        CCPM._pqs_pqs_source_box_raw_weighted_density_density_interaction_block(
            rectangular_pqs_pqs_source_box_plan,
            rectangular_pqs_pqs_source_box_plan;
            term_coefficients = pqs_product_density_coefficients,
            raw_axis_pair_factor_terms = pqs_product_raw_weighted_terms,
            axis_weights = pqs_product_density_weights,
        )
    @test pqs_pqs_raw_weighted_density.path ==
          :pqs_pqs_source_box_raw_weighted_density_density_interaction
    @test pqs_pqs_raw_weighted_density.block ≈
          pqs_pqs_density.block atol = 1.0e-13 rtol = 1.0e-13
    @test pqs_pqs_raw_weighted_density.normalized_axis_pair_factor_terms.x ≈
          pqs_product_density_terms.x atol = 1.0e-14 rtol = 1.0e-14
    @test pqs_pqs_raw_weighted_density.normalized_axis_pair_factor_terms.z ≈
          pqs_product_density_terms.z atol = 1.0e-14 rtol = 1.0e-14
    @test pqs_pqs_raw_weighted_density.density_normalized_core.path ==
          :pqs_pqs_source_box_density_density_interaction
    @test pqs_pqs_raw_weighted_density.diagnostics.pair_factor_normalization ==
          :raw_weighted
    @test pqs_pqs_raw_weighted_density.diagnostics.raw_weighted_pair_factors
    @test !pqs_pqs_raw_weighted_density.diagnostics.density_normalized_pair_factors
    @test pqs_pqs_raw_weighted_density.diagnostics.density_normalized_pair_factors_generated
    @test pqs_pqs_raw_weighted_density.diagnostics.source_weight_division_owner ==
          :source_box_raw_weights
    @test pqs_pqs_raw_weighted_density.diagnostics.source_weight_division_applied_by_helper
    @test pqs_pqs_raw_weighted_density.diagnostics.source_weight_division_shape ==
          :axis_pair_weight_outer
    @test pqs_pqs_raw_weighted_density.diagnostics.density_normalized_core_helper ==
          :_pqs_pqs_source_box_density_density_interaction_block
    @test !pqs_pqs_raw_weighted_density.diagnostics.retained_pqs_weights_used
    @test !pqs_pqs_raw_weighted_density.diagnostics.retained_weight_division_allowed
    @test !pqs_pqs_raw_weighted_density.diagnostics.retained_pqs_weight_division_allowed
    @test !pqs_pqs_raw_weighted_density.diagnostics.ida_weight_division_allowed
    @test !pqs_pqs_raw_weighted_density.diagnostics.shell_projection_used
    @test !pqs_pqs_raw_weighted_density.diagnostics.lowdin_cleanup_used
    @test !pqs_pqs_raw_weighted_density.diagnostics.support_coefficient_matrix_used
    @test !pqs_pqs_raw_weighted_density.diagnostics.packet_adoption
    @test !pqs_pqs_raw_weighted_density.diagnostics.qwhamiltonian_consumes
    @test !pqs_pqs_raw_weighted_density.diagnostics.mwg_ida_semantics_changed
    @test !pqs_pqs_raw_weighted_density.diagnostics.numerical_reference_fallback
    @test_throws ArgumentError CCPM._pqs_pqs_source_box_density_density_interaction_block(
        rectangular_pqs_pqs_source_box_plan,
        rectangular_pqs_pqs_source_box_plan;
        term_coefficients = pqs_product_density_coefficients,
        axis_pair_factor_terms = pqs_product_density_terms,
        axis_weights = pqs_product_density_weights,
        pair_factor_normalization = :raw_weighted,
    )
    @test_throws ArgumentError CCPM._pqs_pqs_source_box_density_density_interaction_block(
        rectangular_pqs_pqs_source_box_plan,
        rectangular_pqs_pqs_source_box_plan;
        term_coefficients = pqs_product_density_coefficients,
        axis_pair_factor_terms = pqs_product_density_terms,
        axis_weights = (
            x = pqs_product_density_weights.x,
            y = pqs_product_density_weights.y,
            z = [1.22, 1.24, 1.26, 0.0, 1.3, 1.32, 1.34],
        ),
    )
    @test_throws ArgumentError CCPM._pqs_pqs_source_box_raw_weighted_density_density_interaction_block(
        rectangular_pqs_pqs_source_box_plan,
        rectangular_pqs_pqs_source_box_plan;
        term_coefficients = pqs_product_density_coefficients,
        raw_axis_pair_factor_terms = pqs_product_raw_weighted_terms,
        axis_weights = (
            x = [1.0, 1.1, 0.0, 1.3, 1.4],
            y = pqs_product_density_weights.y,
            z = pqs_product_density_weights.z,
        ),
    )
    @test_throws ArgumentError CCPM._pqs_pqs_source_box_density_density_interaction_block(
        rectangular_pqs_pqs_source_box_plan,
        rectangular_pqs_pqs_source_box_plan;
        term_coefficients = pqs_product_density_coefficients[1:1],
        axis_pair_factor_terms = pqs_product_density_terms,
        axis_weights = pqs_product_density_weights,
    )
    pqs_pqs_nuclear_attraction =
        CCPM._pqs_pqs_source_box_nuclear_attraction_by_center(
            rectangular_pqs_pqs_source_box_plan,
            rectangular_pqs_pqs_source_box_plan,
            nuclear_axis_layers,
            nuclear_probe_expansion;
            centers = nuclear_probe_centers,
            nuclear_charges = nuclear_probe_charges,
        )
    _check_source_box_nuclear_attraction_by_center(
        pqs_pqs_nuclear_attraction,
        center -> CCPM._pqs_pqs_source_box_centered_local_gaussian_sum_block(
            rectangular_pqs_pqs_source_box_plan,
            rectangular_pqs_pqs_source_box_plan,
            nuclear_axis_layers,
            nuclear_probe_expansion;
            center,
        ),
        :pqs_pqs_source_box_nuclear_attraction_by_center,
    )
    @test_throws ArgumentError CCPM._pqs_pqs_source_box_nuclear_attraction_by_center(
        rectangular_pqs_pqs_source_box_plan,
        rectangular_pqs_pqs_source_box_plan,
        nuclear_axis_layers,
        nuclear_probe_expansion;
        centers = nuclear_probe_centers,
        nuclear_charges = (1.0,),
    )
    @test_throws ArgumentError CCPM._pqs_pqs_source_box_nuclear_attraction_by_center(
        rectangular_pqs_pqs_source_box_plan,
        rectangular_pqs_pqs_source_box_plan,
        nuclear_axis_layers,
        nuclear_probe_expansion;
        centers = nuclear_probe_centers,
        nuclear_charges = (1.0, -0.5),
    )
    rectangular_parent_representation = GaussletBases._cartesian_direct_product_representation(
        (
            x = GaussletBases._cartesian_axis_representation(bundle5.basis),
            y = GaussletBases._cartesian_axis_representation(bundle5.basis),
            z = GaussletBases._cartesian_axis_representation(bundle7.basis),
        );
        route_metadata = (source = :pqs_source_box_gto_shadow_test,),
    )
    gto_probe = CartesianGaussianShellSupplementRepresentation3D(
        :pqs_source_box_gto_shadow_test,
        CartesianGaussianShellOrbitalRepresentation3D[
            CartesianGaussianShellOrbitalRepresentation3D(
                "shifted_s",
                (0, 0, 0),
                (0.15, -0.25, 0.35),
                [0.8, 0.35],
                [0.7, -0.2],
                :axiswise_normalized_cartesian_gaussian,
            ),
            CartesianGaussianShellOrbitalRepresentation3D(
                "shifted_px",
                (1, 0, 0),
                (-0.3, 0.1, -0.2),
                [0.6, 0.25],
                [0.9, 0.15],
                :axiswise_normalized_cartesian_gaussian,
            ),
            CartesianGaussianShellOrbitalRepresentation3D(
                "shifted_dz2",
                (0, 0, 2),
                (0.05, 0.2, -0.45),
                [0.7],
                [1.1],
                :axiswise_normalized_cartesian_gaussian,
            ),
        ],
        (source = :pqs_source_box_gto_shadow_test,),
    )
    _check_pqs_source_box_gto_cross_overlap_shadow(
        rectangular_descriptor,
        rectangular_parent_representation,
        gto_probe;
        expected_source_mode_dims = (5, 5, 7),
        expected_retained_count = 130,
        shared_raw_product_box_plan = rectangular_shared_raw_product_box_plan,
    )
    rectangular_product_metric =
        GaussletBases._nested_projected_q_shell_descriptor_metric_product_contraction(
            rectangular,
            rectangular_bundles,
        )
    @test rectangular_product_metric.overlap == Matrix{Float64}(I, 130, 130)
    @test rectangular_product_metric.weights ≈ rectangular_metric_prototype.weights atol = 1.0e-10 rtol = 1.0e-10
    @test rectangular_product_metric.weights ≈ rectangular.packet.weights atol = 1.0e-10 rtol = 1.0e-10
    @test rectangular_product_metric.position_x ≈ rectangular.packet.position_x atol = 1.0e-10 rtol = 1.0e-10
    @test rectangular_product_metric.position_y ≈ rectangular.packet.position_y atol = 1.0e-10 rtol = 1.0e-10
    @test rectangular_product_metric.position_z ≈ rectangular.packet.position_z atol = 1.0e-10 rtol = 1.0e-10
    @test rectangular_product_metric.coverage.piece_count == 6
    @test rectangular_product_metric.coverage.piece_roles == (:xlo, :xhi, :ylo, :yhi, :zlo, :zhi)
    @test rectangular_product_metric.coverage.support_count == 130
    @test rectangular_product_metric.coverage.unique_support_count == 130
    @test rectangular_product_metric.coverage.duplicate_count == 0
    @test rectangular_product_metric.coverage.missing_count == 0
    @test rectangular_product_metric.coverage.outside_count == 0
    @test rectangular_product_metric.coverage.coverage_ok
    @test rectangular_product_metric.diagnostics.descriptor_kind == :projected_q_shell
    @test rectangular_product_metric.diagnostics.slab_decomposed_product_contraction
    @test !rectangular_product_metric.diagnostics.support_local_boundary_matrix_used
    @test !rectangular_product_metric.diagnostics.dense_full_parent_matrix_used
    @test rectangular_product_metric.diagnostics.boundary_comx_product_modes_used
    @test rectangular_product_metric.diagnostics.raw_boundary_projection_used
    @test rectangular_product_metric.diagnostics.lowdin_cleanup_applied
    @test rectangular_product_metric.diagnostics.overlap_invariant_applied
    @test !rectangular_product_metric.diagnostics.overlap_invariant_debug_check
    @test rectangular_product_metric.diagnostics.overlap_invariant_error == 0.0
    @test !rectangular_product_metric.diagnostics.overlap_is_operator_target
    @test rectangular_product_metric.diagnostics.nontrivial_product_contracted_terms ==
          (:weights, :position_x, :position_y, :position_z)
    @test rectangular_product_metric.diagnostics.mode_axis_indices_cached
    @test rectangular_product_metric.diagnostics.symmetric_mode_matrix_assembly
    @test rectangular_product_metric.diagnostics.axis_piece_blocks_use_views
    @test !rectangular_product_metric.diagnostics.product_doside_unit
    @test !rectangular_product_metric.diagnostics.fixed_block_sidecar_installed
    @test !rectangular_product_metric.diagnostics.optimized_sidecar_installed
    @test rectangular_product_metric.diagnostics.prototype_only

    x_axis_bundles = GaussletBases._CartesianNestedAxisBundles3D(bundle7, bundle5, bundle5)
    x_axis = GaussletBases._nested_projected_q_shell_layer(
        x_axis_bundles,
        (1:7, 1:5, 1:5),
        (2:6, 2:4, 2:4);
        bond_axis = :x,
        q = 5,
        L = 7,
        term_coefficients,
    )
    @test x_axis.diagnostics.bond_axis == :x
    @test x_axis.diagnostics.boundary_support_count == 130
    @test x_axis.diagnostics.retained_count == 130
    @test x_axis.diagnostics.overlap_error < 1.0e-8
end

@testset "Cartesian nested endcap-panel owned shell producer" begin
    function _owned_unit_test_bundle(count::Int)
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
            backend = :numerical_reference,
            refinement_levels = 0,
        )
    end

    function _one_hot_selector_columns(matrix::AbstractMatrix{<:Real})
        dense = Matrix{Float64}(matrix)
        for column in axes(dense, 2)
            rows = findall(!iszero, dense[:, column])
            length(rows) == 1 || return false
            dense[only(rows), column] == 1.0 || return false
        end
        return true
    end

    function _owned_unit_overlap_gram(unit, bundles, dims)
        support_states = [
            GaussletBases._cartesian_unflat_index(index, dims) for index in unit.support_indices
        ]
        pgdg_x = bundles.bundle_x.pgdg_intermediate
        pgdg_y = bundles.bundle_y.pgdg_intermediate
        pgdg_z = bundles.bundle_z.pgdg_intermediate
        support_overlap = GaussletBases._nested_support_product_matrix(
            support_states,
            pgdg_x.overlap,
            pgdg_y.overlap,
            pgdg_z.overlap,
        )
        coefficients = Matrix{Float64}(unit.coefficient_matrix)
        return Matrix{Float64}(transpose(coefficients) * support_overlap * coefficients)
    end

    dims = (5, 5, 7)
    current_box = (1:5, 1:5, 1:7)
    inner_box = (2:4, 2:4, 2:6)
    bundle5 = _owned_unit_test_bundle(5)
    bundle7 = _owned_unit_test_bundle(7)
    bundles = GaussletBases._CartesianNestedAxisBundles3D(bundle5, bundle5, bundle7)
    expected_support = setdiff(
        GaussletBases._nested_box_support_indices(current_box..., dims),
        GaussletBases._nested_box_support_indices(inner_box..., dims),
    )
    sort!(expected_support)
    expected_roles = (
        :endcap_low,
        :endcap_high,
        :panel_y_low,
        :panel_x_high,
        :panel_y_high,
        :panel_x_low,
    )
    expected_support_counts = (25, 25, 20, 20, 20, 20)
    expected_unit_states = (
        sort([(ix, iy, 1) for ix in 1:5 for iy in 1:5]),
        sort([(ix, iy, 7) for ix in 1:5 for iy in 1:5]),
        sort([(ix, 1, iz) for ix in 1:4 for iz in 2:6]),
        sort([(5, iy, iz) for iy in 1:4 for iz in 2:6]),
        sort([(ix, 5, iz) for ix in 2:5 for iz in 2:6]),
        sort([(1, iy, iz) for iy in 2:5 for iz in 2:6]),
    )

    q, L = 4, 4
    retained_count = 2 * q^2 + 4 * q * L
    producer = GaussletBases._nested_endcap_panel_owned_units(
        bundles,
        current_box,
        inner_box;
        bond_axis = :z,
        q,
        L,
    )
    units = producer.units
    owned_support = reduce(vcat, (unit.support_indices for unit in units))
    unit_states = Tuple(
        sort([GaussletBases._cartesian_unflat_index(index, dims) for index in unit.support_indices])
        for unit in units
    )

    @test producer.support_contract == :thin_endcap_box_perimeter
    @test producer.coefficient_contract == :product_doside
    @test producer.expected_support_indices == expected_support
    @test length(expected_support) == 130
    @test producer.audit.expected_support_count == 130
    @test producer.audit.owned_support_count == 130
    @test producer.audit.duplicate_count == 0
    @test producer.audit.missing_count == 0
    @test producer.audit.outside_count == 0
    @test producer.audit.retained_count == 96
    @test producer.audit.retained_count == retained_count
    @test producer.audit.coverage_ok
    @test length(owned_support) == length(unique(owned_support))
    @test sort(owned_support) == expected_support
    @test getfield.(units, :role) == expected_roles
    @test unit_states == expected_unit_states
    @test length.(getfield.(units, :support_indices)) == expected_support_counts
    @test sum(size(unit.coefficient_matrix, 2) for unit in units) == retained_count
    @test all(unit.coefficient_matrix isa SparseMatrixCSC{Float64,Int} for unit in units)
    @test all(size(unit.coefficient_matrix, 1) == length(unit.support_indices) for unit in units)
    @test all(size(unit.coefficient_matrix, 2) == 16 for unit in units)
    @test all(all(isfinite, nonzeros(unit.coefficient_matrix)) for unit in units)
    @test all(!_one_hot_selector_columns(unit.coefficient_matrix) for unit in units)
    @test all(minimum(svdvals(Matrix{Float64}(unit.coefficient_matrix))) > 1.0e-10 for unit in units)
    for unit in units
        gram = _owned_unit_overlap_gram(unit, bundles, dims)
        @test norm(gram - I, Inf) < 1.0e-8
    end
    @test all(unit.metadata.q == q for unit in units)
    @test all(unit.metadata.L == L for unit in units)
    @test all(unit.metadata.bond_axis == :z for unit in units)
    @test all(unit.metadata.support_contract == :thin_endcap_box_perimeter for unit in units)
    @test all(unit.metadata.coefficient_contract == :product_doside for unit in units)
    @test first(units).metadata.retained_count == q * q
    @test last(units).metadata.retained_count == q * L

    direct_selector = GaussletBases._nested_endcap_panel_owned_units(
        dims,
        current_box,
        inner_box;
        bond_axis = :z,
        q,
        L,
        coefficient_contract = :direct_selector,
    )
    @test direct_selector.coefficient_contract == :direct_selector
    @test all(_one_hot_selector_columns(unit.coefficient_matrix) for unit in direct_selector.units)
    @test all(unit.metadata.coefficient_contract == :direct_selector for unit in direct_selector.units)
    @test_throws ArgumentError GaussletBases._nested_endcap_panel_owned_units(
        dims,
        current_box,
        inner_box;
        bond_axis = :z,
        q,
        L,
        coefficient_contract = :product_doside,
    )

    term_coefficients = Float64.(coulomb_gaussian_expansion(doacc = false).coefficients)
    layer = GaussletBases._nested_endcap_panel_shell_layer(
        bundles,
        current_box,
        inner_box;
        bond_axis = :z,
        q,
        L,
        packet_kernel = :support_reference,
        term_coefficients,
    )
    expected_ranges = [1:16, 17:32, 33:48, 49:64, 65:80, 81:96]
    expected_states = [
        GaussletBases._cartesian_unflat_index(index, dims) for index in expected_support
    ]
    outside_support = setdiff(collect(1:prod(dims)), layer.support_indices)

    @test layer isa GaussletBases._CartesianNestedEndcapPanelShellLayer3D
    @test layer.owned_units.coefficient_contract == :product_doside
    @test layer.unit_column_ranges == expected_ranges
    @test size(layer.coefficient_matrix) == (prod(dims), 96)
    @test layer.support_indices == expected_support
    @test layer.support_states == expected_states
    @test nnz(layer.coefficient_matrix[outside_support, :]) == 0
    @test all(all(isfinite, nonzeros(unit.coefficient_matrix)) for unit in layer.owned_units.units)
    @test layer.provenance.support_contract == :thin_endcap_box_perimeter
    @test layer.provenance.coefficient_contract == :product_doside
    @test layer.provenance.packet_kernel == :support_reference
    @test layer.provenance.q == q
    @test layer.provenance.L == L
    @test all(isfinite, layer.packet.overlap)
    @test all(isfinite, layer.packet.kinetic)
    @test all(isfinite, layer.packet.position_x)
    @test all(isfinite, layer.packet.position_y)
    @test all(isfinite, layer.packet.position_z)
    @test all(isfinite, layer.packet.x2_x)
    @test all(isfinite, layer.packet.x2_y)
    @test all(isfinite, layer.packet.x2_z)
    @test all(isfinite, layer.packet.weights)
    @test !isnothing(layer.packet.gaussian_sum)
    @test !isnothing(layer.packet.pair_sum)
    @test all(isfinite, layer.packet.gaussian_sum)
    @test all(isfinite, layer.packet.pair_sum)
    @test norm(layer.packet.overlap - I, Inf) < 1.0e-8

    factorized_layer = GaussletBases._nested_endcap_panel_shell_layer(
        bundles,
        current_box,
        inner_box;
        bond_axis = :z,
        q,
        L,
        packet_kernel = :factorized_direct,
        term_coefficients,
    )
    @test factorized_layer.provenance.packet_kernel == :factorized_direct
    for field in (
        :overlap,
        :kinetic,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :gaussian_sum,
        :pair_sum,
    )
        @test getfield(factorized_layer.packet, field) ≈ getfield(layer.packet, field) atol = 1.0e-9 rtol = 1.0e-9
    end
    @test factorized_layer.packet.weights ≈ layer.packet.weights atol = 1.0e-10 rtol = 1.0e-10
    @test_throws ArgumentError GaussletBases._nested_endcap_panel_shell_layer(
        direct_selector,
        bundles;
        term_coefficients,
    )

    x_axis = GaussletBases._nested_endcap_panel_owned_units(
        bundles,
        current_box,
        inner_box;
        bond_axis = :x,
        q = 3,
        L = 3,
    )
    @test getfield.(x_axis.units, :role) == (
        :endcap_low,
        :endcap_high,
        :panel_z_low,
        :panel_y_high,
        :panel_z_high,
        :panel_y_low,
    )
    @test x_axis.audit.coverage_ok
    @test x_axis.coefficient_contract == :product_doside

    @test_throws ArgumentError GaussletBases._nested_endcap_panel_owned_units(
        bundles,
        current_box,
        inner_box;
        bond_axis = :z,
        q = 0,
        L = 3,
    )
    @test_throws ArgumentError GaussletBases._nested_endcap_panel_owned_units(
        bundles,
        current_box,
        inner_box;
        bond_axis = :z,
        q = 3,
        L = 0,
    )
    @test_throws ArgumentError GaussletBases._nested_endcap_panel_owned_units(
        bundles,
        current_box,
        inner_box;
        bond_axis = :bad,
        q = 3,
        L = 3,
    )
    @test_throws ArgumentError GaussletBases._nested_endcap_panel_owned_units(
        bundles,
        current_box,
        current_box;
        bond_axis = :z,
        q = 3,
        L = 3,
    )
    @test_throws ArgumentError GaussletBases._nested_endcap_panel_owned_units(
        bundles,
        current_box,
        (3:4, 2:4, 2:6);
        bond_axis = :z,
        q = 3,
        L = 3,
    )
    @test_throws ArgumentError GaussletBases._nested_endcap_panel_owned_units(
        bundles,
        current_box,
        inner_box;
        bond_axis = :z,
        q = 5,
        L = 4,
    )
    @test_throws ArgumentError GaussletBases._nested_endcap_panel_owned_units(
        bundles,
        current_box,
        inner_box;
        bond_axis = :z,
        q = 4,
        L = 6,
    )
end

@testset "Bond-aligned diatomic atom-growth anatomy policy" begin
    parent29 = (1:29, 1:29, 1:29)
    recipe29 = GaussletBases._nested_bond_aligned_diatomic_atom_growth_recipe(
        parent29;
        bond_axis = :z,
        atom_axis_indices = (10, 20),
        protected_atom_side_count = 5,
    )
    anatomy29 = GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(recipe29)

    @test recipe29 isa GaussletBases._BondAlignedDiatomicAtomGrowthRecipe3D
    @test recipe29.contact_cap_policy == :single_shared_contact_cap
    @test recipe29.mismatch_absorption_policy == :outermost_shared_molecular_shell
    @test anatomy29.atom_side_count_ladder == [5, 7, 9]
    @test anatomy29.final_atom_side_count == 9
    @test anatomy29.contact_gap_count == 1
    @test anatomy29.contact_policy == :single_shared_contact_cap
    @test anatomy29.left_atom_box == (11:19, 11:19, 6:14)
    @test anatomy29.right_atom_box == (11:19, 11:19, 16:24)
    @test anatomy29.contact_box == (11:19, 11:19, 15:15)
    @test anatomy29.inner_atom_contact_box == (11:19, 11:19, 6:24)
    @test anatomy29.outer_regular_start_box == (6:24, 6:24, 1:29)
    @test anatomy29.regular_shared_shell_count == 5
    @test anatomy29.outer_mismatch_low_counts == (5, 5, 0)
    @test anatomy29.outer_mismatch_high_counts == (5, 5, 0)
    @test anatomy29.support_coverage.status == :full_parent_covered
    @test anatomy29.support_coverage.expected_support_count == 29^3
    @test anatomy29.support_coverage.atom_contact_support_count == 9 * 9 * 19
    @test anatomy29.support_coverage.shared_molecular_support_count == 29^3 - 9 * 9 * 19
    @test anatomy29.support_coverage.covered_support_count == 29^3
    @test anatomy29.support_coverage.duplicate_count == 0
    @test anatomy29.support_coverage.missing_count == 0
    @test anatomy29.support_coverage.outside_count == 0
    @test anatomy29.support_coverage.coverage_ok

    even_parent = (1:20, 1:20, 1:20)
    even_anatomy = GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(
        even_parent;
        bond_axis = :z,
        atom_axis_indices = (7, 14),
        protected_atom_side_count = 4,
    )

    @test even_anatomy.atom_side_count_ladder == [4, 6]
    @test even_anatomy.final_atom_side_count == 6
    @test even_anatomy.contact_gap_count == 0
    @test even_anatomy.contact_policy == :touching_atom_boxes
    @test isnothing(even_anatomy.contact_box)
    @test even_anatomy.left_atom_box == (8:13, 8:13, 5:10)
    @test even_anatomy.right_atom_box == (8:13, 8:13, 11:16)
    @test even_anatomy.left_atom_box[1:2] == even_anatomy.right_atom_box[1:2]
    @test 7 - first(even_anatomy.left_atom_box[3]) == 2
    @test last(even_anatomy.left_atom_box[3]) - 7 == 3
    @test 14 - first(even_anatomy.right_atom_box[3]) == 3
    @test last(even_anatomy.right_atom_box[3]) - 14 == 2
    @test even_anatomy.support_coverage.status == :full_parent_covered
    @test even_anatomy.support_coverage.coverage_ok

    asymmetric_parent = (1:11, 1:13, 1:25)
    asymmetric = GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(
        asymmetric_parent;
        bond_axis = :z,
        atom_axis_indices = (8, 17),
        protected_atom_side_count = 5,
    )

    @test asymmetric.contact_policy == :touching_atom_boxes
    @test asymmetric.inner_atom_contact_box == (2:10, 3:11, 4:21)
    @test asymmetric.outer_regular_start_box == (1:11, 2:12, 3:22)
    @test asymmetric.regular_shared_shell_count == 1
    @test asymmetric.outer_mismatch_low_counts == (0, 1, 2)
    @test asymmetric.outer_mismatch_high_counts == (0, 1, 3)
    @test asymmetric.recipe.mismatch_absorption_policy == :outermost_shared_molecular_shell
    @test asymmetric.support_coverage.expected_support_count == 11 * 13 * 25
    @test asymmetric.support_coverage.covered_support_count == 11 * 13 * 25
    @test asymmetric.support_coverage.coverage_ok

    expansion = coulomb_gaussian_expansion(doacc = false)
    real_basis = bond_aligned_homonuclear_qw_basis(
        family = :G10,
        bond_length = 3.0,
        core_spacing = 0.7,
        xmax_parallel = 8.0,
        xmax_transverse = 4.0,
        bond_axis = :z,
    )
    real_bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(real_basis, expansion)
    real_recipe = GaussletBases._nested_bond_aligned_diatomic_atom_growth_recipe(
        real_basis,
        real_bundles;
        protected_atom_side_count = 3,
    )
    real_anatomy = GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(real_recipe)
    @test real_recipe.parent_box == (1:7, 1:7, 1:13)
    @test real_recipe.atom_axis_indices == (5, 9)
    @test real_anatomy.contact_policy == :single_shared_contact_cap
    @test real_anatomy.support_coverage.coverage_ok

    @test_throws ArgumentError GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(
        (1:7, 1:7, 1:13);
        bond_axis = :z,
        atom_axis_indices = (3, 11),
        protected_atom_side_count = 3,
    )
end

@testset "Bond-aligned diatomic atom-growth construction plan" begin
    parent29 = (1:29, 1:29, 1:29)
    anatomy29 = GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(
        parent29;
        bond_axis = :z,
        atom_axis_indices = (10, 20),
        protected_atom_side_count = 5,
    )
    plan29 = GaussletBases._nested_bond_aligned_diatomic_atom_growth_construction_plan(
        anatomy29,
    )

    @test plan29 isa GaussletBases._BondAlignedDiatomicAtomGrowthConstructionPlan3D
    @test plan29.region_order == [
        :outer_mismatch_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :left_atom_box,
        :right_atom_box,
        :contact_cap,
    ]
    @test [region.order_index for region in plan29.regions] == collect(1:9)
    @test plan29.outer_mismatch_is_outermost
    @test plan29.middle_contact_clean
    @test plan29.support_coverage.expected_support_count == 29^3
    @test plan29.support_coverage.region_support_count == 29^3
    @test plan29.support_coverage.covered_support_count == 29^3
    @test plan29.support_coverage.duplicate_count == 0
    @test plan29.support_coverage.missing_count == 0
    @test plan29.support_coverage.outside_count == 0
    @test plan29.support_coverage.coverage_ok
    @test length(plan29.regions[1].support_indices) == 29^3 - 19 * 19 * 29
    @test [
        length(region.support_indices)
        for region in plan29.regions
        if region.role == :regular_shared_molecular_shell
    ] == [2666, 2178, 1738, 1346, 1002]
    @test length(plan29.regions[end - 2].support_indices) == 9^3
    @test length(plan29.regions[end - 1].support_indices) == 9^3
    @test length(plan29.regions[end].support_indices) == 9 * 9
    @test plan29.regions[1].metadata.mismatch_absorption_policy ==
        :outermost_shared_molecular_shell
    @test plan29.regions[2].metadata.shell_offset == 5
    @test plan29.regions[6].metadata.shell_offset == 1
    @test plan29.regions[end].metadata.contact_policy == :single_shared_contact_cap
    CCP = GaussletBases.CartesianContractedParents
    atom_region = CCP.cartesian_shell_region(
        plan29.regions[end - 2];
        parent_dimension = 29^3,
    )
    shared_region = CCP.cartesian_shell_region(
        plan29.regions[2];
        parent_dimension = 29^3,
    )
    contact_region = CCP.cartesian_shell_region(
        plan29.regions[end];
        parent_dimension = 29^3,
    )
    mismatch_region = CCP.cartesian_shell_region(
        plan29.regions[1];
        parent_dimension = 29^3,
    )
    @test atom_region isa CCP.CartesianShellRegion3D
    @test atom_region.region_family == :atom_core_cube
    @test atom_region.status == :clean
    @test atom_region.role == :left_atom_box
    @test atom_region.support_summary.entry_count == 9^3
    @test atom_region.support_summary.outside_count == 0
    @test atom_region.ownership_coverage_contract == :disjoint_partition_piece
    @test atom_region.retention.retention_rule == :protected_atom_cubic_shell
    @test atom_region.retention.preferred_contraction_rule == :complete_shell_sequence
    @test atom_region.retention.metric_capability == :support_local_product
    @test :coefficient_matrix in atom_region.retention.missing_payload_fields
    @test !atom_region.current_route_consumes
    @test !atom_region.descriptor_drives_builder
    @test atom_region.descriptor_only
    @test shared_region.region_family == :rectangular_molecular_shell
    @test shared_region.retention.retention_rule == :policy_selected_shared_exterior
    @test shared_region.retention.metric_capability == :metadata_only_policy_dependent
    @test contact_region.region_family == :shared_midpoint_slab_cap
    @test contact_region.status == :clean
    @test contact_region.retention.preferred_contraction_rule == :contact_cap_owned_slab
    @test !contact_region.current_route_consumes
    @test mismatch_region.region_family == :outer_boundary_shell
    @test mismatch_region.retention.retention_rule ==
          :outermost_mismatch_shared_molecular_shell
    @test mismatch_region.retention.preferred_contraction_rule ==
          :outer_mismatch_boundary_slab_set
    @test mismatch_region.geometry.mismatch_low_counts == (5, 5, 0)
    @test !mismatch_region.current_route_consumes

    even_plan = GaussletBases._nested_bond_aligned_diatomic_atom_growth_construction_plan(
        (1:20, 1:20, 1:20);
        bond_axis = :z,
        atom_axis_indices = (7, 14),
        protected_atom_side_count = 4,
    )
    @test even_plan.support_coverage.coverage_ok
    @test even_plan.middle_contact_clean
    @test !(:contact_cap in even_plan.region_order)
    @test even_plan.region_order[end - 1:end] == [:left_atom_box, :right_atom_box]

    expansion = coulomb_gaussian_expansion(doacc = false)
    real_basis = bond_aligned_homonuclear_qw_basis(
        family = :G10,
        bond_length = 5.0,
        core_spacing = 0.7,
        xmax_parallel = 8.0,
        xmax_transverse = 4.0,
        bond_axis = :z,
    )
    real_bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(real_basis, expansion)
    real_anatomy = GaussletBases._nested_bond_aligned_diatomic_atom_growth_anatomy(
        real_basis,
        real_bundles;
        protected_atom_side_count = 5,
    )
    real_plan = GaussletBases._nested_bond_aligned_diatomic_atom_growth_construction_plan(
        real_anatomy,
    )
    diagnostic_anatomy = bond_aligned_diatomic_nested_geometry_diagnostics(
        real_basis;
        nside = 5,
    ).atom_growth_anatomy.anatomy
    diagnostic_plan = GaussletBases._nested_bond_aligned_diatomic_atom_growth_construction_plan(
        diagnostic_anatomy,
    )
    @test real_plan.region_order == [
        :outer_mismatch_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :left_atom_box,
        :right_atom_box,
        :contact_cap,
    ]
    @test [length(region.support_indices) for region in real_plan.regions] ==
        [98, 362, 125, 125, 25]
    @test real_plan.support_coverage.expected_support_count == 7 * 7 * 15
    @test real_plan.support_coverage.region_support_count == 7 * 7 * 15
    @test real_plan.support_coverage.covered_support_count == 7 * 7 * 15
    @test real_plan.support_coverage.duplicate_count == 0
    @test real_plan.support_coverage.missing_count == 0
    @test real_plan.support_coverage.outside_count == 0
    @test real_plan.support_coverage.coverage_ok
    @test real_plan.outer_mismatch_is_outermost
    @test real_plan.middle_contact_clean
    @test real_plan.regions[1].metadata.low_counts == (0, 0, 1)
    @test real_plan.regions[1].metadata.high_counts == (0, 0, 1)
    @test diagnostic_plan.region_order == real_plan.region_order
    @test diagnostic_plan.support_coverage.coverage_ok
end

@testset "Bond-aligned diatomic high-order recipe policy metadata" begin
    plan29 = GaussletBases._nested_bond_aligned_diatomic_atom_growth_construction_plan(
        (1:29, 1:29, 1:29);
        bond_axis = :z,
        atom_axis_indices = (10, 20),
        protected_atom_side_count = 5,
    )
    policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        plan29;
        q_min = 4,
        atom_q = 4,
        shared_q = 4,
    )
    diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy_diagnostics(
            policy,
        )

    @test policy isa GaussletBases._BondAlignedDiatomicHighOrderRecipePolicy3D
    @test policy.recipe_label == :mixed_atom_cubic_shared_endcap_panel
    @test policy.q_min == 4
    @test policy.q_region_counts == Dict(4 => 9)
    @test policy.buildable_region_count == 7
    @test policy.planned_region_count == 2
    @test policy.experimental_region_count == 0
    @test policy.metadata.active_builder_uses_policy == false
    @test policy.metadata.source_builder_changed == false
    @test diagnostics.region_count == length(plan29.regions)
    @test diagnostics.active_builder_uses_policy == false
    @test diagnostics.support_coverage.coverage_ok
    @test diagnostics.outer_mismatch_is_outermost
    @test diagnostics.middle_contact_clean

    roles = [choice.region_role for choice in policy.region_choices]
    families = [choice.recipe_family for choice in policy.region_choices]
    categories = [choice.region_category for choice in policy.region_choices]
    statuses = [choice.buildability_status for choice in policy.region_choices]
    support_counts = [choice.support_count for choice in policy.region_choices]
    @test roles == plan29.region_order
    @test families[1] == :outermost_mismatch_shared_molecular_shell
    @test all(==(:shared_endcap_panel_exterior), families[2:6])
    @test all(==(:protected_atom_cubic_shell), families[7:8])
    @test families[9] == :shared_contact_cap
    @test categories[1] == :outer_mismatch
    @test all(==(:shared_exterior), categories[2:6])
    @test all(==(:atom_local), categories[7:8])
    @test categories[9] == :contact_cap
    @test statuses[1] == :planned_only
    @test all(==(:buildable_now), statuses[2:8])
    @test statuses[9] == :planned_only
    @test support_counts == [length(region.support_indices) for region in plan29.regions]
    @test all(choice.q == 4 && choice.order == 4 for choice in policy.region_choices)

    annulus_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        plan29;
        q_min = 4,
        atom_q = 4,
        shared_q = 5,
        contact_q = 4,
        outer_mismatch_q = 4,
        shared_exterior_family = :transverse_annulus_exterior,
    )
    annulus_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy_diagnostics(
            annulus_policy,
        )
    @test annulus_policy.recipe_label == :mixed_atom_cubic_shared_transverse_annulus
    @test annulus_policy.q_region_counts == Dict(4 => 4, 5 => 5)
    @test annulus_policy.buildable_region_count == 2
    @test annulus_policy.planned_region_count == 2
    @test annulus_policy.experimental_region_count == 5
    @test all(
        choice.recipe_family == :transverse_annulus_exterior &&
        choice.q == 5 &&
        choice.implementation_status == :promising_experimental &&
        choice.buildability_status == :planned_experimental
        for choice in annulus_policy.region_choices
        if choice.region_category == :shared_exterior
    )
    @test annulus_diagnostics.region_choices[2].recipe_family ==
        :transverse_annulus_exterior
    @test annulus_diagnostics.region_choices[2].q == 5

    @test_throws ArgumentError GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        plan29;
        q_min = 5,
        atom_q = 4,
    )
end

@testset "Bond-aligned diatomic high-order recipe realization audit" begin
    plan29 = GaussletBases._nested_bond_aligned_diatomic_atom_growth_construction_plan(
        (1:29, 1:29, 1:29);
        bond_axis = :z,
        atom_axis_indices = (10, 20),
        protected_atom_side_count = 5,
    )
    policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        plan29;
        q_min = 4,
    )
    audit =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_realization_audit(
            policy,
        )
    diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_realization_diagnostics(
            audit,
        )

    @test audit isa GaussletBases._BondAlignedDiatomicHighOrderRecipeRealizationAudit3D
    @test diagnostics.recipe_label == :mixed_atom_cubic_shared_endcap_panel
    @test diagnostics.region_count == 9
    @test diagnostics.mapped_region_count == 9
    @test diagnostics.missing_region_count == 0
    @test diagnostics.buildable_without_mapped_primitive_count == 0
    @test diagnostics.active_builder_consumed_region_count == 0
    @test diagnostics.descriptor_region_count == 2
    @test diagnostics.exact_descriptor_region_count == 2
    @test diagnostics.ready_for_opt_in_builder == true
    @test diagnostics.active_builder_uses_policy == false
    @test diagnostics.support_coverage.coverage_ok

    realizations = diagnostics.region_realizations
    @test realizations[1].role == :outer_mismatch_shared_molecular_shell
    @test realizations[1].mapped_primitive ==
        :_nested_diatomic_high_order_outer_mismatch_descriptor
    @test realizations[1].mapped_primitive_status == :construction_piece_descriptor
    @test isnothing(realizations[1].missing_implementation)
    @test realizations[1].realization_descriptor.parent_support_count ==
        length(plan29.regions[1].support_indices)
    @test realizations[1].realization_descriptor.owned_unit_count == 4
    @test realizations[1].realization_descriptor.primitive_family ==
        :outer_mismatch_boundary_slab_set
    @test realizations[1].realization_descriptor.support_coverage.coverage_ok
    @test realizations[1].realization_descriptor.exact_full_coverage
    @test all(
        piece.primitive_family == :outer_mismatch_boundary_slab
        for piece in realizations[1].realization_descriptor.pieces
    )
    @test all(
        realization.mapped_primitive == :_nested_endcap_panel_shell_layer &&
        realization.mapped_primitive_status == :existing_internal_primitive &&
        isnothing(realization.missing_implementation) &&
        realization.existing_opt_in_route == :shared_shell_layer_policy_endcap_panel_owned
        for realization in realizations
        if realization.region_category == :shared_exterior
    )
    shared_realizations = [
        realization for realization in realizations
        if realization.region_category == :shared_exterior
    ]
    @test all(
        !isnothing(realization.metadata.realization_notes.projected_q_shell_candidate)
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.mapped_primitive ==
        :_nested_projected_q_shell_layer
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.primitive_family ==
        :projected_q_shell
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.realization_primitive ==
        :projected_q_shell_boundary_comx_product_modes
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.support_contract ==
        :projected_q_shell_raw_boundary
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.coefficient_contract ==
        :full_block_boundary_comx_product_mode_projection
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.cleanup_contract ==
        :full_rank_symmetric_lowdin
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.mode_selection_rule ==
        :any_axis_mode_index_first_or_last
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.sidecar_status ==
        :not_yet_optimized_product_staged_for_pqs
        for realization in shared_realizations
    )
    @test all(
        !realization.metadata.realization_notes.projected_q_shell_candidate.active_builder_consumes &&
        !realization.metadata.realization_notes.projected_q_shell_candidate.source_builder_consumes &&
        !realization.metadata.realization_notes.projected_q_shell_candidate.fixed_block_consumes &&
        !realization.metadata.realization_notes.projected_q_shell_candidate.qw_consumes &&
        !realization.metadata.realization_notes.projected_q_shell_candidate.hamiltonian_consumes
        for realization in shared_realizations
    )
    @test all(
        realization.metadata.realization_notes.projected_q_shell_candidate.current_transitional_implementation ==
        :_nested_endcap_panel_shell_layer
        for realization in shared_realizations
    )
    @test all(
        realization.mapped_primitive == :_nested_bond_aligned_diatomic_sequence_for_box &&
        realization.mapped_primitive_status == :existing_internal_primitive &&
        isnothing(realization.missing_implementation)
        for realization in realizations
        if realization.region_category == :atom_local
    )
    @test realizations[end].role == :contact_cap
    @test realizations[end].mapped_primitive ==
        :_nested_diatomic_high_order_contact_cap_descriptor
    @test realizations[end].mapped_primitive_status == :construction_piece_descriptor
    @test isnothing(realizations[end].missing_implementation)
    @test realizations[end].realization_descriptor.parent_support_count ==
        length(plan29.regions[end].support_indices)
    @test realizations[end].realization_descriptor.owned_unit_count == 1
    @test realizations[end].realization_descriptor.primitive_family == :contact_cap_owned_slab
    @test realizations[end].realization_descriptor.support_coverage.coverage_ok
    @test realizations[end].realization_descriptor.exact_full_coverage
    @test all(
        !(
            realization.buildability_status == :buildable_now &&
            isnothing(realization.mapped_primitive)
        )
        for realization in realizations
    )
    @test all(!realization.active_builder_consumes for realization in realizations)

    annulus_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        plan29;
        q_min = 4,
        atom_q = 4,
        shared_q = 5,
        contact_q = 4,
        outer_mismatch_q = 4,
        shared_exterior_family = :transverse_annulus_exterior,
    )
    annulus_audit =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_realization_audit(
            annulus_policy,
        )
    annulus_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_realization_diagnostics(
            annulus_audit,
        )
    @test annulus_diagnostics.recipe_label == :mixed_atom_cubic_shared_transverse_annulus
    @test annulus_diagnostics.mapped_region_count == 4
    @test annulus_diagnostics.missing_region_count == 5
    @test annulus_diagnostics.buildable_without_mapped_primitive_count == 0
    @test annulus_diagnostics.descriptor_region_count == 2
    @test annulus_diagnostics.exact_descriptor_region_count == 2
    @test annulus_diagnostics.ready_for_opt_in_builder == false
    @test all(
        realization.recipe_family == :transverse_annulus_exterior &&
        realization.q == 5 &&
        realization.mapped_primitive_status == :missing_experimental_primitive &&
        realization.missing_implementation == :transverse_annulus_owned_unit_producer &&
        !(:projected_q_shell_candidate in propertynames(realization.metadata.realization_notes))
        for realization in annulus_diagnostics.region_realizations
        if realization.region_category == :shared_exterior
    )
end

@testset "Bond-aligned diatomic high-order recipe opt-in source construction" begin
    basis = bond_aligned_homonuclear_qw_basis(
        family = :G10,
        bond_length = 5.0,
        core_spacing = 0.7,
        xmax_parallel = 8.0,
        xmax_transverse = 4.0,
        bond_axis = :z,
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(basis, expansion)
    policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        basis,
        bundles;
        protected_atom_side_count = 5,
        q_min = 4,
    )
    realization_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_realization_diagnostics(
            GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_realization_audit(
                policy,
            ),
        )
    @test realization_diagnostics.ready_for_opt_in_builder
    @test !realization_diagnostics.active_builder_uses_policy
    @test realization_diagnostics.active_builder_consumed_region_count == 0

    construction =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :factorized_direct,
            build_sequence_packet = false,
        )
    diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction_diagnostics(
            construction,
        )

    @test construction isa
          GaussletBases._BondAlignedDiatomicHighOrderRecipeSourceConstruction3D
    @test diagnostics.recipe_label == :mixed_atom_cubic_shared_endcap_panel
    @test diagnostics.active_builder_consumes
    @test diagnostics.active_builder_uses_policy
    @test diagnostics.metadata.default_source_builder_changed == false
    @test diagnostics.consumed_region_count == diagnostics.region_count == 5
    @test diagnostics.unsupported_region_count == 0
    @test diagnostics.parent_dimension == 7 * 7 * 15
    @test diagnostics.fixed_dimension == 469
    @test diagnostics.support_coverage.expected_support_count == 7 * 7 * 15
    @test diagnostics.support_coverage.coverage_ok
    @test isnothing(construction.sequence.packet)

    region_builds = diagnostics.region_builds
    @test [build.role for build in region_builds] == [
        :outer_mismatch_shared_molecular_shell,
        :regular_shared_molecular_shell,
        :left_atom_box,
        :right_atom_box,
        :contact_cap,
    ]
    @test [build.primitive_family for build in region_builds] == [
        :outer_mismatch_boundary_slab_set,
        :shared_endcap_panel_shell_layer,
        :atom_local_complete_shell_sequence,
        :atom_local_complete_shell_sequence,
        :contact_cap_owned_slab,
    ]
    @test [build.built_support_count for build in region_builds] == [98, 362, 125, 125, 25]
    @test [build.retained_count for build in region_builds] == [98, 96, 125, 125, 25]
    @test [build.column_range for build in region_builds] ==
          [1:98, 374:469, 99:223, 224:348, 349:373]
    @test all(build.built && build.active_builder_consumes for build in region_builds)
    @test all(build.support_coverage.coverage_ok for build in region_builds)
    @test region_builds[2].metadata.support_contract == :thin_endcap_box_perimeter
    @test region_builds[2].metadata.coefficient_contract == :product_doside
    @test region_builds[5].metadata.descriptor_scope == :middle_contact_cap
    CCP = GaussletBases.CartesianContractedParents
    inventory = CCP.cartesian_shell_region_inventory(
        construction;
        parent_dimension = diagnostics.parent_dimension,
    )
    @test inventory isa CCP.CartesianShellRegionInventory3D
    @test inventory.region_count == diagnostics.region_count == length(region_builds)
    @test inventory.region_order == [build.role for build in region_builds]
    @test [region.role for region in inventory.regions] == inventory.region_order
    @test [region.status for region in inventory.regions] ==
          [:clean, :transitional, :clean, :clean, :clean]
    @test [region.retention.preferred_contraction_rule for region in inventory.regions] == [
        :outer_mismatch_boundary_slab_set,
        :old_endcap_panel_product_split,
        :complete_shell_sequence,
        :complete_shell_sequence,
        :contact_cap_owned_slab,
    ]
    @test inventory.status_counts.clean == 4
    @test inventory.status_counts.transitional == 1
    @test inventory.current_route_consumes_count == diagnostics.region_count
    @test inventory.descriptor_only_count == 0
    @test inventory.support_summary.region_support_entry_count ==
          diagnostics.support_coverage.covered_support_count
    @test inventory.support_summary.support_complete_by_region_counts
    @test inventory.support_summary.count_only_summaries_for_all_regions
    @test all(isnothing(region.support_indices) for region in inventory.regions)
    @test all(region.current_route_consumes for region in inventory.regions)
    @test all(!region.descriptor_drives_builder for region in inventory.regions)
    @test all(!region.descriptor_only for region in inventory.regions)
    @test inventory.regions[2].ownership_coverage_contract == :boundary_only
    @test inventory.regions[2].retention.metric_capability ==
          :product_staged_metric_contraction
    @test isempty(inventory.regions[2].retention.missing_payload_fields)
    @test inventory.regions[5].ownership_coverage_contract == :disjoint_partition_piece
    @test isempty(inventory.regions[5].retention.missing_payload_fields)
    @test diagnostics.metadata.q_policy == :atom_growth_endcap_panel_shared_q_variable
    @test diagnostics.metadata.q4_acceptance_fixture
    @test diagnostics.metadata.non_shared_q_policy == :fixed_q4_order4
    @test diagnostics.metadata.shared_q_values == (4,)
    @test diagnostics.metadata.shared_order_values == (4,)
    @test diagnostics.metadata.shared_shell_realization == :endcap_panel_owned
    @test !diagnostics.metadata.projected_q_shell_opt_in

    be2_basis = bond_aligned_homonuclear_qw_basis(
        family = :G10,
        bond_length = 5.0,
        core_spacing = 0.15,
        xmax_parallel = 10.5,
        xmax_transverse = 8.0,
        bond_axis = :z,
        nuclear_charge = 4.0,
    )
    be2_bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(be2_basis, expansion)
    @test GaussletBases._nested_axis_lengths(be2_bundles) == (15, 15, 27)
    be2_policy0 = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        be2_basis,
        be2_bundles;
        protected_atom_side_count = 5,
        q_min = 4,
    )
    be2_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        be2_policy0.construction_plan;
        q_min = 4,
        atom_q = 4,
        atom_order = 4,
        shared_q = 5,
        shared_order = 5,
        contact_q = 4,
        contact_order = 4,
        outer_mismatch_q = 4,
        outer_mismatch_order = 4,
    )
    be2_retention = GaussletBases._nested_resolve_complete_shell_retention(5)
    be2_shared_dimensions = [
        GaussletBases._nested_diatomic_projected_q_shell_adaptive_source_dimensions(
            be2_basis,
            be2_bundles,
            region,
            be2_retention;
            bond_axis = :z,
            nside = 5,
            selected_q = choice.q,
            shared_shell_angular_resolution_scale = 1.4,
        )
        for (choice, region) in zip(
            be2_policy.region_choices,
            be2_policy.construction_plan.regions,
        )
        if choice.recipe_family == :shared_endcap_panel_exterior
    ]
    be2_source_box_dimension_plans = [
        GaussletBases._nested_diatomic_source_box_dimension_plan(
            be2_basis,
            be2_bundles,
            region.box,
            something(region.inner_exclusion_box),
            be2_retention;
            bond_axis = :z,
            nside = 5,
            selected_q = choice.q,
            shared_shell_angular_resolution_scale = 1.4,
            support_count = length(region.support_indices),
        )
        for (choice, region) in zip(
            be2_policy.region_choices,
            be2_policy.construction_plan.regions,
        )
        if choice.recipe_family == :shared_endcap_panel_exterior
    ]
    be2_shared_regions = [
        region
        for (choice, region) in zip(
            be2_policy.region_choices,
            be2_policy.construction_plan.regions,
        )
        if choice.recipe_family == :shared_endcap_panel_exterior
    ]
    pqs_plan = GaussletBases._nested_projected_q_shell_source_mode_plan(
        (5, 5, 6);
        bond_axis = :z,
        selected_q = 5,
        physical_box_lengths = (11, 11, 21),
        support_count = 1002,
    )
    @test pqs_plan.raw_source_dims == (5, 5, 6)
    @test pqs_plan.source_mode_dims == (5, 5, 6)
    @test pqs_plan.axis_selector_retained_counts == (3, 3, 4)
    @test pqs_plan.raw_q == 5
    @test pqs_plan.raw_L == 6
    @test pqs_plan.raw_q_matches_selected_q
    @test pqs_plan.physical_box_lengths == (11, 11, 21)
    @test pqs_plan.support_count == 1002
    @test pqs_plan.pqs_retained_count == 114
    @test pqs_plan.decomposition_status == :adaptive_broad_support_q_local_modes
    @test !pqs_plan.broad_parent_boundary_reference
    @test !pqs_plan.excluded_from_mvp_gate
    pqs_mismatch_plan = GaussletBases._nested_projected_q_shell_source_mode_plan(
        (4, 4, 5);
        bond_axis = :z,
        selected_q = 5,
        physical_box_lengths = (9, 9, 9),
        support_count = 488,
    )
    @test !pqs_mismatch_plan.raw_q_matches_selected_q
    @test pqs_mismatch_plan.decomposition_status == :adaptive_raw_q_mismatch
    @test pqs_mismatch_plan.excluded_from_mvp_gate
    @test !(:adaptive_retain in propertynames(pqs_plan))
    @test [length.(region.box) for region in be2_shared_regions] ==
          [(15, 15, 25), (13, 13, 23), (11, 11, 21)]
    @test [dims.raw_source_dims for dims in be2_shared_dimensions] ==
          [(5, 5, 5), (5, 5, 5), (5, 5, 6)]
    @test [plan.dimension_policy for plan in be2_source_box_dimension_plans] ==
          fill(:diatomic_adaptive_angular_source_box, 3)
    @test [plan.source_mode_dims for plan in be2_source_box_dimension_plans] ==
          [(5, 5, 5), (5, 5, 5), (5, 5, 6)]
    @test [plan.side_dimensions for plan in be2_source_box_dimension_plans] ==
          [(5, 5, 5), (5, 5, 5), (5, 5, 6)]
    @test [plan.raw_L for plan in be2_source_box_dimension_plans] == [5, 5, 6]
    @test all(plan -> plan.total_source_dimensions_primary, be2_source_box_dimension_plans)
    @test all(
        plan -> plan.diagnostics.source_mode_dims_are_total_lengths,
        be2_source_box_dimension_plans,
    )
    @test all(
        plan -> plan.diagnostics.axis_selector_retained_counts_are_diagnostic,
        be2_source_box_dimension_plans,
    )
    @test [dims.pqs_retained_count for dims in be2_shared_dimensions] ==
          [98, 98, 114]
    @test [
        dims.source_box_dimension_plan.source_mode_dims
        for dims in be2_shared_dimensions
    ] == [(5, 5, 5), (5, 5, 5), (5, 5, 6)]
    @test all(
        dims -> !(:adaptive_retain in propertynames(dims)),
        be2_shared_dimensions,
    )
    @test all(dims -> dims.raw_q == 5, be2_shared_dimensions)
    @test all(dims -> dims.raw_q_matches_selected_q, be2_shared_dimensions)
    @test all(
        dims -> dims.decomposition_status == :adaptive_broad_support_q_local_modes,
        be2_shared_dimensions,
    )
    @test all(dims -> !dims.broad_parent_boundary_reference, be2_shared_dimensions)
    @test all(dims -> !dims.excluded_from_mvp_gate, be2_shared_dimensions)
    @test [
        GaussletBases._nested_diatomic_projected_q_shell_retained_count(length.(region.box))
        for region in be2_shared_regions
    ] == [1738, 1346, 1002]
    @test [
        dims.pqs_retained_count !=
        GaussletBases._nested_diatomic_projected_q_shell_retained_count(length.(region.box))
        for (dims, region) in zip(be2_shared_dimensions, be2_shared_regions)
    ] == [true, true, true]

    pqs_construction =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :support_reference,
            shared_shell_realization = :projected_q_shell,
        )
    pqs_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction_diagnostics(
            pqs_construction,
        )
    @test pqs_diagnostics.metadata.shared_shell_realization == :projected_q_shell
    @test pqs_diagnostics.metadata.projected_q_shell_opt_in
    @test pqs_diagnostics.metadata.default_source_builder_changed == false
    @test pqs_diagnostics.metadata.packet_kernel == :support_reference
    @test pqs_diagnostics.metadata.build_sequence_packet
    @test pqs_diagnostics.support_coverage.coverage_ok
    @test pqs_diagnostics.support_coverage.expected_support_count == 7 * 7 * 15
    @test pqs_diagnostics.support_coverage.covered_support_count == 7 * 7 * 15
    @test pqs_diagnostics.support_coverage.duplicate_count == 0
    @test pqs_diagnostics.support_coverage.missing_count == 0
    @test pqs_diagnostics.support_coverage.outside_count == 0
    @test [build.role for build in pqs_diagnostics.region_builds] ==
          [build.role for build in region_builds]
    @test [build.primitive_family for build in pqs_diagnostics.region_builds] == [
        :outer_mismatch_boundary_slab_set,
        :projected_q_shell,
        :atom_local_complete_shell_sequence,
        :atom_local_complete_shell_sequence,
        :contact_cap_owned_slab,
    ]
    @test [
        build.built_support_count for build in pqs_diagnostics.region_builds
    ] == [build.built_support_count for build in region_builds]
    @test pqs_diagnostics.region_builds[2].mapped_primitive ==
          :_nested_projected_q_shell_layer
    @test pqs_diagnostics.region_builds[2].metadata.support_contract ==
          :projected_q_shell_raw_boundary
    @test pqs_diagnostics.region_builds[2].metadata.coefficient_contract ==
          :full_block_boundary_comx_product_mode_projection
    @test pqs_diagnostics.region_builds[2].metadata.seed_contract ==
          :raw_boundary_projection_of_boundary_comx_product_modes_from_full_local_block_transform
    @test pqs_diagnostics.region_builds[2].metadata.cleanup_contract ==
          :full_rank_symmetric_lowdin
    @test pqs_diagnostics.region_builds[2].metadata.cleanup_method ==
          :projected_boundary_symmetric_lowdin
    @test !pqs_diagnostics.region_builds[2].metadata.pqs_product_staged_sidecar_available
    @test !pqs_diagnostics.region_builds[2].metadata.factorized_direct_allowed
    @test !pqs_diagnostics.region_builds[2].metadata.active_default_builder_changed
    @test pqs_diagnostics.region_builds[2].metadata.selected_q == 4
    @test pqs_diagnostics.region_builds[2].metadata.raw_source_dims == (5, 5, 6)
    @test pqs_diagnostics.region_builds[2].metadata.source_mode_dims == (5, 5, 6)
    @test !(
        :adaptive_retain in
        propertynames(pqs_diagnostics.region_builds[2].metadata)
    )
    @test pqs_diagnostics.region_builds[2].metadata.raw_q == 5
    @test pqs_diagnostics.region_builds[2].metadata.raw_L == 6
    @test !pqs_diagnostics.region_builds[2].metadata.raw_q_matches_selected_q
    @test pqs_diagnostics.region_builds[2].metadata.physical_box_lengths == (7, 7, 13)
    @test pqs_diagnostics.region_builds[2].metadata.support_count == 362
    @test pqs_diagnostics.region_builds[2].metadata.pqs_retained_count == 114
    @test pqs_diagnostics.region_builds[2].metadata.decomposition_status ==
          :adaptive_raw_q_mismatch
    @test !pqs_diagnostics.region_builds[2].metadata.broad_parent_boundary_reference
    @test pqs_diagnostics.region_builds[2].metadata.excluded_from_mvp_gate
    @test pqs_diagnostics.region_builds[2].metadata.policy_q == 4
    @test pqs_diagnostics.region_builds[2].metadata.policy_order == 4
    pqs_source_descriptor =
        pqs_diagnostics.region_builds[2].metadata.pqs_staged_unit_descriptor
    @test pqs_source_descriptor.kind == :projected_q_shell
    @test length.(pqs_source_descriptor.current_box) == (7, 7, 13)
    @test all(
        axis ->
            first(pqs_source_descriptor.inner_box[axis]) ==
            first(pqs_source_descriptor.current_box[axis]) + 1 &&
            last(pqs_source_descriptor.inner_box[axis]) ==
            last(pqs_source_descriptor.current_box[axis]) - 1,
        1:3,
    )
    @test pqs_source_descriptor.bond_axis == :z
    @test pqs_source_descriptor.q == 5
    @test pqs_source_descriptor.L == 6
    @test pqs_source_descriptor.support_count == 362
    @test pqs_source_descriptor.mode_count == 114
    @test pqs_source_descriptor.retained_count == 114
    @test pqs_source_descriptor.cleanup_method == :projected_boundary_symmetric_lowdin
    @test pqs_source_descriptor.cleanup_matrix_size == (114, 114)
    @test pqs_source_descriptor.cleanup_rank_count == 114
    @test pqs_source_descriptor.cleanup_rank_drop_count == 0
    @test pqs_source_descriptor.selection_rule == :any_axis_mode_index_first_or_last
    @test all(
        mode -> any(axis -> mode[axis] == 1 || mode[axis] == (5, 5, 6)[axis], 1:3),
        pqs_source_descriptor.boundary_mode_indices,
    )
    pqs_route_fact_diagnostic =
        CCPM._pqs_pqs_product_route_descriptor_diagnostic(
            pqs_construction,
            _pqs_axis_metrics(bundles),
        )
    @test pqs_route_fact_diagnostic.status == :descriptor_unavailable
    @test pqs_route_fact_diagnostic.descriptor === nothing
    @test :second_pqs_raw_plan in pqs_route_fact_diagnostic.missing
    @test :middle_product_doside_unit in pqs_route_fact_diagnostic.missing
    @test pqs_route_fact_diagnostic.diagnostics.pqs_descriptor_count == 1
    @test pqs_route_fact_diagnostic.diagnostics.pqs_raw_plan_convertible_count == 1
    @test pqs_route_fact_diagnostic.diagnostics.product_doside_unit_count == 0
    @test pqs_route_fact_diagnostic.diagnostics.direct_or_support_body_piece_count == 4
    @test !pqs_route_fact_diagnostic.diagnostics.descriptor_emitted
    @test pqs_route_fact_diagnostic.diagnostics.raw_plan_convertibility_checked
    @test pqs_route_fact_diagnostic.diagnostics.raw_plan_conversion_failure_count == 0
    @test pqs_route_fact_diagnostic.diagnostics.current_route_contains_pqs_descriptors
    @test !pqs_route_fact_diagnostic.diagnostics.current_route_contains_explicit_product_doside_body_unit
    @test !pqs_route_fact_diagnostic.diagnostics.packet_adoption
    @test !pqs_route_fact_diagnostic.diagnostics.fixed_block_construction_changed
    @test !pqs_route_fact_diagnostic.diagnostics.qwhamiltonian_changed
    @test !pqs_route_fact_diagnostic.diagnostics.shell_projection_used
    @test !pqs_route_fact_diagnostic.diagnostics.lowdin_cleanup_used
    @test !pqs_route_fact_diagnostic.diagnostics.support_local_pqs_oracle_used
    @test pqs_route_fact_diagnostic.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !pqs_route_fact_diagnostic.diagnostics.ida_weight_division_allowed
    @test !pqs_route_fact_diagnostic.diagnostics.direct_support_reinterpreted_as_product_doside
    @test any(
        mismatch ->
            hasproperty(mismatch, :reason) &&
            mismatch.reason == :direct_or_support_piece_not_product_doside &&
            hasproperty(mismatch, :primitive_family) &&
            mismatch.primitive_family == :contact_cap_owned_slab,
        pqs_route_fact_diagnostic.mismatches,
    )
    @test any(
        mismatch ->
            hasproperty(mismatch, :reason) &&
            mismatch.reason == :shared_pqs_descriptors_are_not_route_left_right_group,
        pqs_route_fact_diagnostic.mismatches,
    )
    pqs_retained_unit_audit =
        CCPM._pqs_route_retained_unit_fact_audit(pqs_construction)
    @test pqs_retained_unit_audit.object_kind ==
          :pqs_route_retained_unit_fact_audit
    @test pqs_retained_unit_audit.status == :audit_only
    @test length(pqs_retained_unit_audit.unit_facts) == 5
    @test pqs_retained_unit_audit.summary.classification_counts.product_box_constructible == 2
    @test pqs_retained_unit_audit.summary.classification_counts.needs_direct_support_retained_unit_kind == 2
    @test pqs_retained_unit_audit.summary.classification_counts.out_of_scope == 1
    @test pqs_retained_unit_audit.summary.product_box_constructible_slab_rule_count == 3
    pqs_retained_facts_by_role =
        Dict(fact.role => fact for fact in pqs_retained_unit_audit.unit_facts)
    contact_fact = pqs_retained_facts_by_role[:contact_cap]
    @test contact_fact.classification == :product_box_constructible
    @test contact_fact.primitive_family == :contact_cap_owned_slab
    @test !contact_fact.raw_product_box_operator_contract
    @test contact_fact.product_box_construction_rule_available
    @test !contact_fact.product_doside_unit
    @test contact_fact.safe_term_capability ==
          :product_box_rule_available_not_instantiated
    @test contact_fact.coefficient_scope == :support_local_direct_rows
    @test contact_fact.construction_rule.rule_kind ==
          :identity_selector_product_doside_slab
    @test contact_fact.construction_rule.fixed_axis == :z
    @test contact_fact.construction_rule.fixed_index == 8
    @test contact_fact.construction_rule.active_axes == (:x, :y)
    @test contact_fact.construction_rule.active_intervals == (2:6, 2:6)
    @test contact_fact.construction_rule.retained_count == 25
    outer_fact =
        pqs_retained_facts_by_role[:outer_mismatch_shared_molecular_shell]
    @test outer_fact.classification == :product_box_constructible
    @test outer_fact.primitive_family == :outer_mismatch_boundary_slab_set
    @test !outer_fact.raw_product_box_operator_contract
    @test outer_fact.product_box_construction_rule_available
    @test !outer_fact.product_doside_unit
    @test outer_fact.construction_rule.rule_kind ==
          :identity_selector_product_doside_boundary_slab_set
    @test outer_fact.construction_rule.boundary_slab_set
    @test outer_fact.construction_rule.slab_piece_count == 2
    @test all(
        piece_rule -> piece_rule.fixed_axis == :z,
        outer_fact.construction_rule.slab_piece_rules,
    )
    @test sum(
        piece_rule -> piece_rule.support_count,
        outer_fact.construction_rule.slab_piece_rules,
    ) == outer_fact.support_count
    outer_mismatch_product_units =
        CCPM._pqs_outer_mismatch_product_doside_units(pqs_construction)
    @test outer_mismatch_product_units.object_kind ==
          :pqs_outer_mismatch_product_doside_units_fixture
    @test outer_mismatch_product_units.status == :private_diagnostic_only
    @test outer_mismatch_product_units.fact.role ==
          :outer_mismatch_shared_molecular_shell
    @test outer_mismatch_product_units.fact.primitive_family ==
          :outer_mismatch_boundary_slab_set
    @test !outer_mismatch_product_units.fact.raw_product_box_operator_contract
    @test outer_mismatch_product_units.fact.product_box_construction_rule_available
    @test length(outer_mismatch_product_units.units) == 2
    @test map(unit -> unit.kind, outer_mismatch_product_units.units) ==
          (:product_doside, :product_doside)
    @test map(unit -> unit.role, outer_mismatch_product_units.units) ==
          (:outer_mismatch_z_low_slab, :outer_mismatch_z_high_slab)
    @test map(unit -> unit.column_range, outer_mismatch_product_units.units) ==
          (1:49, 50:98)
    @test map(unit -> length(unit.support_indices), outer_mismatch_product_units.units) ==
          (49, 49)
    outer_mismatch_piece_support_indices = vcat(
        outer_mismatch_product_units.units[1].support_indices,
        outer_mismatch_product_units.units[2].support_indices,
    )
    @test sort(outer_mismatch_piece_support_indices) ==
          sort(outer_mismatch_product_units.fact.support_indices)
    @test map(
        unit -> map(axis -> axis.kind, unit.axes),
        outer_mismatch_product_units.units,
    ) == ((:active, :active, :fixed), (:active, :active, :fixed))
    @test outer_mismatch_product_units.units[1].axes[1].interval == 1:7
    @test outer_mismatch_product_units.units[1].axes[2].interval == 1:7
    @test outer_mismatch_product_units.units[1].axes[3].fixed_index == 1
    @test outer_mismatch_product_units.units[2].axes[1].interval == 1:7
    @test outer_mismatch_product_units.units[2].axes[2].interval == 1:7
    @test outer_mismatch_product_units.units[2].axes[3].fixed_index == 15
    @test all(
        unit -> Matrix(unit.coefficient_matrix) == Matrix{Float64}(I, 49, 49),
        outer_mismatch_product_units.units,
    )
    @test all(
        unit -> unit.axis_function_indices ==
                GaussletBases._nested_product_axis_function_indices(3, 1, 7, 2, 7),
        outer_mismatch_product_units.units,
    )
    @test all(
        equivalence -> equivalence.support_indices_match,
        outer_mismatch_product_units.piece_equivalences,
    )
    @test all(
        equivalence -> equivalence.retained_count_match,
        outer_mismatch_product_units.piece_equivalences,
    )
    @test map(
        equivalence -> equivalence.column_range,
        outer_mismatch_product_units.piece_equivalences,
    ) == (1:49, 50:98)
    @test all(
        equivalence -> equivalence.coefficient_matrix_matches_direct_selector,
        outer_mismatch_product_units.piece_equivalences,
    )
    @test all(
        equivalence -> equivalence.max_parent_coefficient_error == 0.0,
        outer_mismatch_product_units.piece_equivalences,
    )
    @test outer_mismatch_product_units.aggregate_equivalence.support_indices_match
    @test outer_mismatch_product_units.aggregate_equivalence.audited_support_set_match
    @test outer_mismatch_product_units.aggregate_equivalence.retained_count_match
    @test outer_mismatch_product_units.aggregate_equivalence.column_range_partition
    @test outer_mismatch_product_units.aggregate_equivalence.coefficient_matrix_matches_direct_selector
    @test outer_mismatch_product_units.aggregate_equivalence.max_parent_coefficient_error == 0.0
    @test outer_mismatch_product_units.diagnostics.outer_mismatch_only
    @test outer_mismatch_product_units.diagnostics.boundary_slab_set
    @test outer_mismatch_product_units.diagnostics.product_doside_units_created
    @test outer_mismatch_product_units.diagnostics.unit_count == 2
    @test outer_mismatch_product_units.diagnostics.slab_piece_count == 2
    @test !outer_mismatch_product_units.diagnostics.route_descriptor_emitted
    @test !outer_mismatch_product_units.diagnostics.construction_mutated
    @test !outer_mismatch_product_units.diagnostics.sidecar_installation
    @test !outer_mismatch_product_units.diagnostics.packet_adoption
    @test !outer_mismatch_product_units.diagnostics.fixed_block_construction_changed
    @test !outer_mismatch_product_units.diagnostics.qwhamiltonian_changed
    @test !outer_mismatch_product_units.diagnostics.ida_weight_division_allowed
    @test outer_mismatch_product_units.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !outer_mismatch_product_units.diagnostics.input_fact_raw_product_box_operator_contract
    @test outer_mismatch_product_units.diagnostics.created_units_raw_product_box_operator_contract
    @test outer_mismatch_product_units.diagnostics.descriptor_piece_order_defines_columns
    @test outer_mismatch_product_units.diagnostics.audited_support_checked_as_set
    @test outer_mismatch_product_units.diagnostics.product_box_construction_rule_available
    @test !outer_mismatch_product_units.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    outer_mismatch_safe_term_metrics = _pqs_axis_metrics(bundles)
    outer_mismatch_safe_terms = (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    outer_mismatch_safe_term_comparison =
        CCPM._pqs_outer_mismatch_safe_term_operator_comparison(
            pqs_construction,
            outer_mismatch_safe_term_metrics,
        )
    @test outer_mismatch_safe_term_comparison.object_kind ==
          :pqs_outer_mismatch_safe_term_operator_comparison
    @test outer_mismatch_safe_term_comparison.status == :private_diagnostic_only
    @test outer_mismatch_safe_term_comparison.terms == outer_mismatch_safe_terms
    @test length(outer_mismatch_safe_term_comparison.fixture.units) == 2
    @test outer_mismatch_safe_term_comparison.max_block_error <= 1.0e-12
    @test outer_mismatch_safe_term_comparison.diagnostics.source ==
          :pqs_outer_mismatch_safe_term_operator_comparison
    @test outer_mismatch_safe_term_comparison.diagnostics.outer_mismatch_only
    @test outer_mismatch_safe_term_comparison.diagnostics.boundary_slab_set
    @test outer_mismatch_safe_term_comparison.diagnostics.private_diagnostic_only
    @test outer_mismatch_safe_term_comparison.diagnostics.terms_checked ==
          outer_mismatch_safe_terms
    @test outer_mismatch_safe_term_comparison.diagnostics.supported_terms ==
          outer_mismatch_safe_terms
    @test :weights in outer_mismatch_safe_term_comparison.diagnostics.unsupported_terms
    @test outer_mismatch_safe_term_comparison.diagnostics.product_path ==
          :_product_doside_source_box_reference_block
    @test outer_mismatch_safe_term_comparison.diagnostics.direct_oracle_path ==
          :support_local_direct_selector_contract_pair_block
    @test outer_mismatch_safe_term_comparison.diagnostics.product_doside_units_created
    @test outer_mismatch_safe_term_comparison.diagnostics.complete_slab_set_block_assembled
    @test outer_mismatch_safe_term_comparison.diagnostics.cross_slab_blocks_included
    @test outer_mismatch_safe_term_comparison.diagnostics.direct_support_oracle_compared
    @test !outer_mismatch_safe_term_comparison.diagnostics.route_descriptor_emitted
    @test !outer_mismatch_safe_term_comparison.diagnostics.construction_mutated
    @test !outer_mismatch_safe_term_comparison.diagnostics.sidecar_installation
    @test !outer_mismatch_safe_term_comparison.diagnostics.packet_adoption
    @test !outer_mismatch_safe_term_comparison.diagnostics.fixed_block_construction_changed
    @test !outer_mismatch_safe_term_comparison.diagnostics.qwhamiltonian_changed
    @test !outer_mismatch_safe_term_comparison.diagnostics.ida_weight_division_allowed
    @test outer_mismatch_safe_term_comparison.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !outer_mismatch_safe_term_comparison.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test outer_mismatch_safe_term_comparison.diagnostics.operator_factor_source ==
          :explicit_metric_operator_data
    @test outer_mismatch_safe_term_comparison.diagnostics.operator_metric_sources ==
          (:nested_pgdg_axis, :nested_pgdg_axis, :nested_pgdg_axis)
    @test !outer_mismatch_safe_term_comparison.diagnostics.input_metric_operator_data_pgdg_checked
    @test !outer_mismatch_safe_term_comparison.diagnostics.pgdg_analytic_operator_provenance_claimed
    @test !outer_mismatch_safe_term_comparison.diagnostics.numerical_reference_fallback
    @test outer_mismatch_safe_term_comparison.diagnostics.product_source_box_reference_compared
    @test outer_mismatch_safe_term_comparison.diagnostics.direct_support_oracle_entries_built
    @test outer_mismatch_safe_term_comparison.diagnostics.retained_count == 98
    @test outer_mismatch_safe_term_comparison.diagnostics.support_count == 98
    @test outer_mismatch_safe_term_comparison.diagnostics.unit_count == 2
    @test outer_mismatch_safe_term_comparison.diagnostics.max_block_error <= 1.0e-12
    @test outer_mismatch_safe_term_comparison.diagnostics.output_finite
    for term in outer_mismatch_safe_terms
        product_block = outer_mismatch_safe_term_comparison.product_blocks[term]
        oracle_block = outer_mismatch_safe_term_comparison.direct_oracle_blocks[term]
        @test size(product_block) == (98, 98)
        @test size(oracle_block) == (98, 98)
        @test all(isfinite, product_block)
        @test all(isfinite, oracle_block)
        @test product_block ≈ oracle_block atol = 1.0e-12 rtol = 0.0
        @test outer_mismatch_safe_term_comparison.term_errors[term] <= 1.0e-12
        @test outer_mismatch_safe_term_comparison.diagnostics.pair_block_counts[term] == 4
        @test outer_mismatch_safe_term_comparison.diagnostics.cross_slab_pair_counts[term] == 2
        @test length(outer_mismatch_safe_term_comparison.product_references[term].pair_references) == 4
    end
    @test_throws ArgumentError CCPM._pqs_outer_mismatch_safe_term_operator_comparison(
        pqs_construction,
        outer_mismatch_safe_term_metrics;
        terms = (:weights,),
    )
    left_atom_fact = pqs_retained_facts_by_role[:left_atom_box]
    right_atom_fact = pqs_retained_facts_by_role[:right_atom_box]
    @test left_atom_fact.classification ==
          :needs_direct_support_retained_unit_kind
    @test right_atom_fact.classification ==
          :needs_direct_support_retained_unit_kind
    @test !left_atom_fact.product_doside_unit
    @test !right_atom_fact.product_doside_unit
    @test !left_atom_fact.raw_product_box_operator_contract
    @test !right_atom_fact.raw_product_box_operator_contract
    @test !left_atom_fact.product_box_construction_rule_available
    @test !right_atom_fact.product_box_construction_rule_available
    @test left_atom_fact.safe_term_capability == :support_local_reference_only
    @test right_atom_fact.safe_term_capability == :support_local_reference_only
    atom_box_support_dense_units =
        CCPM._pqs_atom_box_support_dense_units(pqs_construction)
    @test atom_box_support_dense_units.object_kind ==
          :pqs_atom_box_support_dense_units_fixture
    @test atom_box_support_dense_units.status == :private_diagnostic_only
    @test length(atom_box_support_dense_units.units) == 2
    @test map(unit -> unit.role, atom_box_support_dense_units.units) ==
          (:left_atom_box, :right_atom_box)
    @test map(unit -> unit.kind, atom_box_support_dense_units.units) ==
          (:support_dense, :support_dense)
    @test map(unit -> unit.column_range, atom_box_support_dense_units.units) ==
          (99:223, 224:348)
    @test map(unit -> length(unit.support_indices), atom_box_support_dense_units.units) ==
          (125, 125)
    @test map(unit -> size(unit.coefficient_matrix), atom_box_support_dense_units.units) ==
          ((125, 125), (125, 125))
    @test all(
        unit -> map(axis -> axis.kind, unit.axes) == (:fixed, :fixed, :fixed),
        atom_box_support_dense_units.units,
    )
    @test all(
        unit -> all(index -> index == (1, 1, 1), unit.axis_function_indices),
        atom_box_support_dense_units.units,
    )
    @test all(
        equivalence -> equivalence.support_indices_match,
        atom_box_support_dense_units.equivalences,
    )
    @test all(
        equivalence -> equivalence.audited_support_set_match,
        atom_box_support_dense_units.equivalences,
    )
    @test all(
        equivalence -> equivalence.retained_count_match,
        atom_box_support_dense_units.equivalences,
    )
    @test all(
        equivalence -> equivalence.coefficient_matrix_matches_direct_support,
        atom_box_support_dense_units.equivalences,
    )
    @test all(
        equivalence -> equivalence.max_parent_coefficient_error == 0.0,
        atom_box_support_dense_units.equivalences,
    )
    @test all(
        equivalence -> equivalence.local_identity_error <= 1.0e-12,
        atom_box_support_dense_units.equivalences,
    )
    @test atom_box_support_dense_units.diagnostics.atom_box_only
    @test atom_box_support_dense_units.diagnostics.support_dense_direct_support_units_created
    @test !atom_box_support_dense_units.diagnostics.product_doside_units_created
    @test !atom_box_support_dense_units.diagnostics.raw_product_box_operator_contract
    @test atom_box_support_dense_units.diagnostics.support_local_reference_only
    @test !atom_box_support_dense_units.diagnostics.product_box_construction_rule_available
    @test !atom_box_support_dense_units.diagnostics.route_descriptor_emitted
    @test !atom_box_support_dense_units.diagnostics.construction_mutated
    @test !atom_box_support_dense_units.diagnostics.sidecar_installation
    @test !atom_box_support_dense_units.diagnostics.packet_adoption
    @test !atom_box_support_dense_units.diagnostics.fixed_block_construction_changed
    @test !atom_box_support_dense_units.diagnostics.qwhamiltonian_changed
    @test !atom_box_support_dense_units.diagnostics.ida_weight_division_allowed
    @test atom_box_support_dense_units.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !atom_box_support_dense_units.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test atom_box_support_dense_units.diagnostics.max_parent_coefficient_error == 0.0
    @test !atom_box_support_dense_units.diagnostics.local_identity_is_product_box_claim
    @test !atom_box_support_dense_units.diagnostics.safe_term_operator_comparison_added
    atom_box_safe_term_metrics = _pqs_axis_metrics(bundles)
    atom_box_safe_terms = (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    atom_box_safe_term_comparison =
        CCPM._pqs_atom_box_safe_term_operator_comparison(
            pqs_construction,
            atom_box_safe_term_metrics,
        )
    @test atom_box_safe_term_comparison.object_kind ==
          :pqs_atom_box_safe_term_operator_comparison
    @test atom_box_safe_term_comparison.status == :private_diagnostic_only
    @test atom_box_safe_term_comparison.terms == atom_box_safe_terms
    @test length(atom_box_safe_term_comparison.fixture.units) == 2
    @test atom_box_safe_term_comparison.max_block_error <= 1.0e-12
    @test atom_box_safe_term_comparison.diagnostics.source ==
          :pqs_atom_box_safe_term_operator_comparison
    @test atom_box_safe_term_comparison.diagnostics.atom_box_only
    @test atom_box_safe_term_comparison.diagnostics.support_dense_direct_support_units_created
    @test atom_box_safe_term_comparison.diagnostics.support_local_fallback_operator_comparison
    @test !atom_box_safe_term_comparison.diagnostics.product_doside_units_created
    @test !atom_box_safe_term_comparison.diagnostics.raw_product_box_operator_contract
    @test !atom_box_safe_term_comparison.diagnostics.product_box_construction_rule_available
    @test atom_box_safe_term_comparison.diagnostics.complete_atom_box_block_assembled
    @test atom_box_safe_term_comparison.diagnostics.cross_atom_blocks_included
    @test atom_box_safe_term_comparison.diagnostics.direct_support_oracle_compared
    @test !atom_box_safe_term_comparison.diagnostics.route_descriptor_emitted
    @test !atom_box_safe_term_comparison.diagnostics.construction_mutated
    @test !atom_box_safe_term_comparison.diagnostics.sidecar_installation
    @test !atom_box_safe_term_comparison.diagnostics.packet_adoption
    @test !atom_box_safe_term_comparison.diagnostics.fixed_block_construction_changed
    @test !atom_box_safe_term_comparison.diagnostics.qwhamiltonian_changed
    @test !atom_box_safe_term_comparison.diagnostics.ida_weight_division_allowed
    @test atom_box_safe_term_comparison.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !atom_box_safe_term_comparison.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test atom_box_safe_term_comparison.diagnostics.operator_factor_source ==
          :explicit_metric_operator_data
    @test atom_box_safe_term_comparison.diagnostics.operator_metric_sources ==
          (:nested_pgdg_axis, :nested_pgdg_axis, :nested_pgdg_axis)
    @test !atom_box_safe_term_comparison.diagnostics.input_metric_operator_data_pgdg_checked
    @test !atom_box_safe_term_comparison.diagnostics.pgdg_analytic_operator_provenance_claimed
    @test !atom_box_safe_term_comparison.diagnostics.numerical_reference_fallback
    @test atom_box_safe_term_comparison.diagnostics.retained_count == 250
    @test atom_box_safe_term_comparison.diagnostics.support_count == 250
    @test atom_box_safe_term_comparison.diagnostics.unit_count == 2
    @test atom_box_safe_term_comparison.diagnostics.max_block_error <= 1.0e-12
    @test atom_box_safe_term_comparison.diagnostics.output_finite
    for term in atom_box_safe_terms
        support_dense_block = atom_box_safe_term_comparison.support_dense_blocks[term]
        oracle_block = atom_box_safe_term_comparison.direct_oracle_blocks[term]
        @test size(support_dense_block) == (250, 250)
        @test size(oracle_block) == (250, 250)
        @test all(isfinite, support_dense_block)
        @test all(isfinite, oracle_block)
        @test support_dense_block ≈ oracle_block atol = 1.0e-12 rtol = 0.0
        @test atom_box_safe_term_comparison.term_errors[term] <= 1.0e-12
        @test atom_box_safe_term_comparison.diagnostics.pair_block_counts[term] == 4
        @test atom_box_safe_term_comparison.diagnostics.cross_atom_pair_counts[term] == 2
    end
    @test_throws ArgumentError CCPM._pqs_atom_box_safe_term_operator_comparison(
        pqs_construction,
        atom_box_safe_term_metrics;
        terms = (:weights,),
    )
    shared_pqs_fact =
        pqs_retained_facts_by_role[:regular_shared_molecular_shell]
    @test shared_pqs_fact.classification == :out_of_scope
    @test shared_pqs_fact.primitive_family == :projected_q_shell
    @test !shared_pqs_fact.product_doside_unit
    @test !shared_pqs_fact.raw_product_box_operator_contract
    @test !shared_pqs_fact.product_box_construction_rule_available
    @test shared_pqs_fact.safe_term_capability == :not_body_retained_unit
    @test shared_pqs_fact.notes.current_single_pqs_descriptor
    @test pqs_retained_unit_audit.diagnostics.private_diagnostic_only
    @test !pqs_retained_unit_audit.diagnostics.descriptor_emitted
    @test !pqs_retained_unit_audit.diagnostics.packet_adoption
    @test !pqs_retained_unit_audit.diagnostics.fixed_block_construction_changed
    @test !pqs_retained_unit_audit.diagnostics.qwhamiltonian_changed
    @test !pqs_retained_unit_audit.diagnostics.sidecar_mutation
    @test !pqs_retained_unit_audit.diagnostics.sidecar_installation
    @test !pqs_retained_unit_audit.diagnostics.direct_support_reinterpreted_as_product_doside
    @test pqs_retained_unit_audit.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !pqs_retained_unit_audit.diagnostics.ida_weight_division_allowed
    @test !pqs_retained_unit_audit.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    contact_cap_product_unit =
        CCPM._pqs_contact_cap_product_doside_unit(pqs_construction)
    @test contact_cap_product_unit.object_kind ==
          :pqs_contact_cap_product_doside_unit_fixture
    @test contact_cap_product_unit.status == :private_diagnostic_only
    @test contact_cap_product_unit.fact.role == :contact_cap
    @test contact_cap_product_unit.fact.classification ==
          :product_box_constructible
    @test !contact_cap_product_unit.fact.raw_product_box_operator_contract
    @test contact_cap_product_unit.fact.product_box_construction_rule_available
    @test contact_cap_product_unit.unit.role == :contact_cap_slab
    @test contact_cap_product_unit.unit.kind == :product_doside
    @test contact_cap_product_unit.unit.column_range == 349:373
    @test contact_cap_product_unit.unit.support_indices ==
          contact_cap_product_unit.fact.support_indices
    @test contact_cap_product_unit.unit.support_states == [
        GaussletBases._cartesian_unflat_index(index, (7, 7, 15))
        for index in contact_cap_product_unit.unit.support_indices
    ]
    @test Matrix(contact_cap_product_unit.unit.coefficient_matrix) ==
          Matrix{Float64}(I, 25, 25)
    @test map(axis -> axis.kind, contact_cap_product_unit.unit.axes) ==
          (:active, :active, :fixed)
    @test contact_cap_product_unit.unit.axes[1].interval == 2:6
    @test contact_cap_product_unit.unit.axes[2].interval == 2:6
    @test contact_cap_product_unit.unit.axes[3].fixed_index == 8
    @test Matrix(contact_cap_product_unit.unit.axes[1].coefficient_matrix) ==
          Matrix{Float64}(I, 5, 5)
    @test Matrix(contact_cap_product_unit.unit.axes[2].coefficient_matrix) ==
          Matrix{Float64}(I, 5, 5)
    @test contact_cap_product_unit.unit.axis_function_indices ==
          GaussletBases._nested_product_axis_function_indices(3, 1, 5, 2, 5)
    @test contact_cap_product_unit.equivalence.support_indices_match
    @test contact_cap_product_unit.equivalence.support_states_match
    @test contact_cap_product_unit.equivalence.retained_count_match
    @test contact_cap_product_unit.equivalence.column_range_match
    @test contact_cap_product_unit.equivalence.coefficient_matrix_matches_direct_selector
    @test contact_cap_product_unit.equivalence.max_parent_coefficient_error == 0.0
    @test contact_cap_product_unit.diagnostics.contact_cap_only
    @test contact_cap_product_unit.diagnostics.product_doside_unit_created
    @test !contact_cap_product_unit.diagnostics.route_descriptor_emitted
    @test !contact_cap_product_unit.diagnostics.construction_mutated
    @test !contact_cap_product_unit.diagnostics.sidecar_installation
    @test !contact_cap_product_unit.diagnostics.packet_adoption
    @test !contact_cap_product_unit.diagnostics.fixed_block_construction_changed
    @test !contact_cap_product_unit.diagnostics.qwhamiltonian_changed
    @test !contact_cap_product_unit.diagnostics.ida_weight_division_allowed
    @test contact_cap_product_unit.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test contact_cap_product_unit.diagnostics.product_box_construction_rule_available
    @test !contact_cap_product_unit.diagnostics.input_fact_raw_product_box_operator_contract
    @test contact_cap_product_unit.diagnostics.created_unit_raw_product_box_operator_contract
    contact_safe_term_metrics = _pqs_axis_metrics(bundles)
    contact_safe_term_comparison =
        CCPM._pqs_contact_cap_safe_term_operator_comparison(
            pqs_construction,
            contact_safe_term_metrics,
        )
    contact_safe_terms = (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    @test contact_safe_term_comparison.object_kind ==
          :pqs_contact_cap_safe_term_operator_comparison
    @test contact_safe_term_comparison.status == :private_diagnostic_only
    @test contact_safe_term_comparison.terms == contact_safe_terms
    @test contact_safe_term_comparison.fixture.unit.kind == :product_doside
    @test contact_safe_term_comparison.fixture.unit.column_range == 349:373
    @test contact_safe_term_comparison.max_block_error <= 1.0e-12
    @test contact_safe_term_comparison.diagnostics.source ==
          :pqs_contact_cap_safe_term_operator_comparison
    @test contact_safe_term_comparison.diagnostics.contact_cap_only
    @test contact_safe_term_comparison.diagnostics.private_diagnostic_only
    @test contact_safe_term_comparison.diagnostics.terms_checked ==
          contact_safe_terms
    @test contact_safe_term_comparison.diagnostics.supported_terms ==
          contact_safe_terms
    @test :weights in contact_safe_term_comparison.diagnostics.unsupported_terms
    @test contact_safe_term_comparison.diagnostics.product_path ==
          :_product_doside_source_box_reference_block
    @test contact_safe_term_comparison.diagnostics.direct_oracle_path ==
          :support_local_direct_selector_contract_pair_block
    @test contact_safe_term_comparison.diagnostics.current_direct_support_selector_compared
    @test contact_safe_term_comparison.diagnostics.product_doside_unit_created
    @test !contact_safe_term_comparison.diagnostics.input_fact_raw_product_box_operator_contract
    @test contact_safe_term_comparison.diagnostics.created_unit_raw_product_box_operator_contract
    @test contact_safe_term_comparison.diagnostics.product_box_construction_rule_available
    @test !contact_safe_term_comparison.diagnostics.route_descriptor_emitted
    @test !contact_safe_term_comparison.diagnostics.construction_mutated
    @test !contact_safe_term_comparison.diagnostics.sidecar_installation
    @test !contact_safe_term_comparison.diagnostics.packet_adoption
    @test !contact_safe_term_comparison.diagnostics.fixed_block_construction_changed
    @test !contact_safe_term_comparison.diagnostics.qwhamiltonian_changed
    @test !contact_safe_term_comparison.diagnostics.ida_weight_division_allowed
    @test contact_safe_term_comparison.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !contact_safe_term_comparison.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test contact_safe_term_comparison.diagnostics.operator_factor_source ==
          :explicit_metric_operator_data
    @test contact_safe_term_comparison.diagnostics.operator_metric_sources ==
          (:nested_pgdg_axis, :nested_pgdg_axis, :nested_pgdg_axis)
    @test !contact_safe_term_comparison.diagnostics.input_metric_operator_data_pgdg_checked
    @test !contact_safe_term_comparison.diagnostics.pgdg_analytic_operator_provenance_claimed
    @test !contact_safe_term_comparison.diagnostics.numerical_reference_fallback
    @test contact_safe_term_comparison.diagnostics.product_source_box_reference_compared
    @test contact_safe_term_comparison.diagnostics.direct_support_oracle_entries_built
    @test contact_safe_term_comparison.diagnostics.retained_count == 25
    @test contact_safe_term_comparison.diagnostics.support_count == 25
    @test contact_safe_term_comparison.diagnostics.column_range == 349:373
    @test contact_safe_term_comparison.diagnostics.max_block_error <= 1.0e-12
    @test contact_safe_term_comparison.diagnostics.output_finite
    for term in contact_safe_terms
        product_block = contact_safe_term_comparison.product_blocks[term]
        oracle_block = contact_safe_term_comparison.direct_oracle_blocks[term]
        @test size(product_block) == (25, 25)
        @test size(oracle_block) == (25, 25)
        @test all(isfinite, product_block)
        @test all(isfinite, oracle_block)
        @test product_block ≈ oracle_block atol = 1.0e-12 rtol = 0.0
        @test contact_safe_term_comparison.term_errors[term] <= 1.0e-12
        @test contact_safe_term_comparison.product_references[term].diagnostics.authoritative_block_compared
    end
    @test_throws ArgumentError CCPM._pqs_contact_cap_safe_term_operator_comparison(
        pqs_construction,
        contact_safe_term_metrics;
        terms = (:weights,),
    )
    current_route_inventory =
        CCPM._pqs_current_route_retained_unit_inventory(pqs_construction)
    @test current_route_inventory.object_kind ==
          :pqs_current_route_retained_unit_inventory_fixture
    @test current_route_inventory.status == :private_diagnostic_only
    @test length(current_route_inventory.units) == 6
    @test current_route_inventory.coverage.ordered_roles == (
        :outer_mismatch_z_low_slab,
        :outer_mismatch_z_high_slab,
        :left_atom_box,
        :right_atom_box,
        :contact_cap_slab,
        :regular_shared_molecular_shell,
    )
    @test current_route_inventory.coverage.ordered_column_ranges == (
        1:49,
        50:98,
        99:223,
        224:348,
        349:373,
        374:487,
    )
    @test current_route_inventory.coverage.first_column == 1
    @test current_route_inventory.coverage.last_column == 487
    @test current_route_inventory.coverage.represented_count == 487
    @test current_route_inventory.coverage.contiguous
    @test current_route_inventory.coverage.non_overlapping
    @test current_route_inventory.coverage.covers_every_column_once
    @test map(unit -> unit.category, current_route_inventory.units) == (
        :product_doside,
        :product_doside,
        :support_dense,
        :support_dense,
        :product_doside,
        :shell_realized_pqs_fixture,
    )
    @test map(unit -> unit.retained_count, current_route_inventory.units) ==
          (49, 49, 125, 125, 25, 114)
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].kind ==
          :projected_q_shell
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].active_representation_stage ==
          :shell_realized_pqs_fixture
    @test !current_route_inventory.by_role[:regular_shared_molecular_shell].raw_product_box_operator_contract
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].safe_term_capability ==
          :support_local_oracle_for_shell_realization
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].support_count == 362
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].retained_count == 114
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].raw_box_auxiliary_metadata.available
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].raw_box_auxiliary_metadata.reference_only
    @test !current_route_inventory.by_role[:regular_shared_molecular_shell].raw_box_auxiliary_metadata.active_current_route_contract
    shared_pqs_unit = current_route_inventory.by_role[:regular_shared_molecular_shell]
    shared_transform_fact = shared_pqs_unit.shell_realization_transform_fact
    @test shared_transform_fact.object_kind ==
          :pqs_current_route_shell_realization_transform_fact
    @test shared_transform_fact.status == :metadata_precursor
    @test shared_transform_fact.representation_stage == :shell_realized_pqs_fixture
    @test shared_transform_fact.source_box.source_mode_dims ==
          shared_pqs_unit.raw_box_auxiliary_metadata.source_mode_dims
    @test shared_transform_fact.source_box.source_mode_count ==
          prod(shared_transform_fact.source_box.source_mode_dims)
    @test shared_transform_fact.boundary_selection.mode_count == shared_pqs_unit.retained_count
    @test shared_transform_fact.shell_projection.support_count == shared_pqs_unit.support_count
    @test shared_transform_fact.shell_projection.matrix_shape ==
          (shared_pqs_unit.support_count, shared_pqs_unit.retained_count)
    @test shared_transform_fact.lowdin_cleanup.transform_shape ==
          (shared_pqs_unit.retained_count, shared_pqs_unit.retained_count)
    @test shared_transform_fact.retained_columns.support_local_coefficient_shape ==
          size(shared_pqs_unit.support_local_coefficient_matrix)
    @test shared_transform_fact.retained_columns.coefficient_matches_descriptor_realization
    @test shared_transform_fact.retained_columns.max_support_local_coefficient_error <=
          1.0e-12
    @test shared_transform_fact.shell_realization.shell_projection_used
    @test shared_transform_fact.shell_realization.lowdin_cleanup_used
    @test !shared_transform_fact.compact_source_space_transform.available
    @test !shared_transform_fact.source_box_operator_application_ready
    @test shared_transform_fact.diagnostics.support_local_oracle_used
    @test shared_transform_fact.diagnostics.shell_row_oracle_only
    @test shared_transform_fact.diagnostics.metadata_precursor
    shared_transform_fact_checked =
        CCPM._pqs_current_route_shell_realization_transform_fact(
            shared_pqs_unit;
            metrics = contact_safe_term_metrics,
        )
    @test shared_transform_fact_checked.shell_realization.isometry_checked
    @test shared_transform_fact_checked.shell_realization.isometry_error <= 1.0e-8
    @test shared_transform_fact_checked.shell_realization.isometric
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].equivalence.coefficient_matrix_matches_active_shell
    @test current_route_inventory.by_role[:regular_shared_molecular_shell].equivalence.max_parent_coefficient_error == 0.0
    @test current_route_inventory.by_role[:left_atom_box].category == :support_dense
    @test current_route_inventory.by_role[:left_atom_box].safe_term_capability ==
          :support_local_fallback_safe_terms
    @test current_route_inventory.by_role[:contact_cap_slab].category == :product_doside
    @test current_route_inventory.by_role[:outer_mismatch_z_low_slab].category ==
          :product_doside
    @test map(policy -> policy.pair_type, current_route_inventory.pair_policies) == (
        :product_product,
        :support_support,
        :support_product,
        :shell_realized_pqs_product,
        :shell_realized_pqs_support,
        :shell_realized_pqs_pqs,
        :raw_box_pqs_helpers,
    )
    @test current_route_inventory.pair_policies[1].policy ==
          :product_doside_source_box_path
    @test current_route_inventory.pair_policies[4].policy ==
          :support_local_oracle_for_shell_realization
    @test !current_route_inventory.pair_policies[4].active_current_route
    @test !current_route_inventory.pair_policies[4].active_algorithmic_policy
    @test !current_route_inventory.pair_policies[4].source_box_algorithm_available
    @test current_route_inventory.pair_policies[4].support_local_oracle_used
    @test current_route_inventory.pair_policies[4].shell_row_oracle_only
    @test !current_route_inventory.pair_policies[end].active_current_route
    @test current_route_inventory.diagnostics.private_diagnostic_only
    @test current_route_inventory.diagnostics.current_route_inventory
    @test !current_route_inventory.diagnostics.route_descriptor_emitted
    @test !current_route_inventory.diagnostics.construction_mutated
    @test !current_route_inventory.diagnostics.sidecar_installation
    @test !current_route_inventory.diagnostics.packet_adoption
    @test !current_route_inventory.diagnostics.fixed_block_construction_changed
    @test !current_route_inventory.diagnostics.qwhamiltonian_changed
    @test !current_route_inventory.diagnostics.ida_weight_division_allowed
    @test current_route_inventory.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !current_route_inventory.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test current_route_inventory.diagnostics.shared_pqs_active_representation ==
          :shell_realized_pqs_fixture
    @test !current_route_inventory.diagnostics.shared_pqs_raw_box_operator_contract
    @test current_route_inventory.diagnostics.shared_pqs_shell_realization_transform_fact_count == 1
    @test current_route_inventory.diagnostics.shared_pqs_source_box_operator_application_ready_count == 0
    @test current_route_inventory.diagnostics.raw_box_pqs_auxiliary_reference_available
    @test !current_route_inventory.diagnostics.whole_route_safe_term_matrix_consumer
    @test current_route_inventory.diagnostics.fixed_dimension == 487
    @test current_route_inventory.diagnostics.coverage_complete
    current_route_pair_inventory =
        CCPM._pqs_current_route_retained_pair_inventory(current_route_inventory)
    @test current_route_pair_inventory.object_kind ==
          :pqs_current_route_retained_pair_inventory_fixture
    @test current_route_pair_inventory.status == :private_diagnostic_only
    @test current_route_pair_inventory.unit_inventory === current_route_inventory
    @test length(current_route_pair_inventory.pairs) == 21
    expected_pair_roles = Tuple(
        (current_route_inventory.units[left].role, current_route_inventory.units[right].role)
        for left in 1:length(current_route_inventory.units)
        for right in left:length(current_route_inventory.units)
    )
    expected_pair_shapes = Tuple(
        (
            current_route_inventory.units[left].retained_count,
            current_route_inventory.units[right].retained_count,
        ) for left in 1:length(current_route_inventory.units)
        for right in left:length(current_route_inventory.units)
    )
    @test map(
        pair -> (pair.left_role, pair.right_role),
        current_route_pair_inventory.pairs,
    ) == expected_pair_roles
    @test map(pair -> pair.pair_shape, current_route_pair_inventory.pairs) ==
          expected_pair_shapes
    @test current_route_pair_inventory.counts.pair_count == 21
    @test current_route_pair_inventory.counts.product_product == 6
    @test current_route_pair_inventory.counts.support_support == 3
    @test current_route_pair_inventory.counts.support_product == 6
    @test current_route_pair_inventory.counts.shell_realized_pqs_product == 3
    @test current_route_pair_inventory.counts.shell_realized_pqs_support == 2
    @test current_route_pair_inventory.counts.shell_realized_pqs_pqs == 1
    @test current_route_pair_inventory.counts.raw_box_pqs_active == 0
    @test current_route_pair_inventory.counts.active_algorithmic_policy == 15
    @test current_route_pair_inventory.counts.source_box_algorithm_available == 6
    @test current_route_pair_inventory.counts.product_doside_source_box_path == 6
    @test current_route_pair_inventory.counts.support_local_fallback == 9
    @test current_route_pair_inventory.counts.support_local_oracle_for_shell_realization == 6
    @test current_route_pair_inventory.counts.support_local_oracle_pair_count == 6
    @test current_route_pair_inventory.counts.shell_row_oracle_pair_count == 6
    @test current_route_pair_inventory.pairs[1].pair_group == :product_product
    @test current_route_pair_inventory.pairs[1].policy ==
          :product_doside_source_box_path
    @test current_route_pair_inventory.pairs[3].pair_group == :support_product
    @test current_route_pair_inventory.pairs[3].policy == :support_local_fallback
    @test current_route_pair_inventory.pairs[6].pair_group ==
          :shell_realized_pqs_product
    @test current_route_pair_inventory.pairs[6].policy ==
          :support_local_oracle_for_shell_realization
    @test current_route_pair_inventory.pairs[6].shell_row_oracle_only
    @test current_route_pair_inventory.pairs[6].support_local_oracle_used
    @test !current_route_pair_inventory.pairs[6].active_algorithmic_policy
    @test !current_route_pair_inventory.pairs[6].source_box_algorithm_available
    @test current_route_pair_inventory.pairs[15].pair_group ==
          :shell_realized_pqs_support
    @test current_route_pair_inventory.pairs[15].policy ==
          :support_local_oracle_for_shell_realization
    @test current_route_pair_inventory.pairs[end].pair_group ==
          :shell_realized_pqs_pqs
    @test current_route_pair_inventory.pairs[end].policy ==
          :support_local_oracle_for_shell_realization
    @test current_route_pair_inventory.pairs[end].shell_row_oracle_only
    @test current_route_pair_inventory.pairs[end].support_local_oracle_used
    @test count(pair -> pair.active_current_route, current_route_pair_inventory.pairs) == 15
    @test count(pair -> pair.active_algorithmic_policy, current_route_pair_inventory.pairs) == 15
    @test count(pair -> pair.shell_row_oracle_only, current_route_pair_inventory.pairs) == 6
    @test all(
        pair -> !pair.raw_box_pqs_active_pair_policy,
        current_route_pair_inventory.pairs,
    )
    @test current_route_pair_inventory.diagnostics.private_diagnostic_only
    @test current_route_pair_inventory.diagnostics.current_route_pair_inventory
    @test current_route_pair_inventory.diagnostics.unit_inventory_complete
    @test current_route_pair_inventory.diagnostics.upper_triangular_pairs
    @test current_route_pair_inventory.diagnostics.pair_count == 21
    @test current_route_pair_inventory.diagnostics.raw_box_pqs_active_pair_policy_count == 0
    @test current_route_pair_inventory.diagnostics.active_algorithmic_policy_pair_count == 15
    @test current_route_pair_inventory.diagnostics.source_box_algorithm_available_pair_count == 6
    @test current_route_pair_inventory.diagnostics.support_local_oracle_for_shell_realization_pair_count == 6
    @test current_route_pair_inventory.diagnostics.shell_row_oracle_pair_count == 6
    @test current_route_pair_inventory.diagnostics.shell_realized_pqs_pairs_are_oracle_only
    @test !current_route_pair_inventory.diagnostics.route_descriptor_emitted
    @test !current_route_pair_inventory.diagnostics.construction_mutated
    @test !current_route_pair_inventory.diagnostics.sidecar_installation
    @test !current_route_pair_inventory.diagnostics.packet_adoption
    @test !current_route_pair_inventory.diagnostics.fixed_block_construction_changed
    @test !current_route_pair_inventory.diagnostics.qwhamiltonian_changed
    @test !current_route_pair_inventory.diagnostics.ida_weight_division_allowed
    @test current_route_pair_inventory.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !current_route_pair_inventory.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test !current_route_pair_inventory.diagnostics.whole_route_safe_term_matrix_consumer
    be2_inventory_timed = @timed begin
        be2_pqs_construction =
            GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
                be2_basis,
                be2_bundles,
                be2_policy;
                nside = 5,
                term_coefficients = Float64.(expansion.coefficients),
                packet_kernel = :support_reference,
                shared_shell_realization = :projected_q_shell,
            )
        be2_current_route_inventory =
            CCPM._pqs_current_route_retained_unit_inventory(be2_pqs_construction)
        be2_current_route_pair_inventory =
            CCPM._pqs_current_route_retained_pair_inventory(
                be2_current_route_inventory,
            )
        (
            construction = be2_pqs_construction,
            inventory = be2_current_route_inventory,
            pair_inventory = be2_current_route_pair_inventory,
        )
    end
    be2_inventory_payload = be2_inventory_timed.value
    be2_current_route_inventory = be2_inventory_payload.inventory
    be2_current_route_pair_inventory = be2_inventory_payload.pair_inventory
    be2_shared_pqs_units = Tuple(
        unit for unit in be2_current_route_inventory.units
        if unit.category == :shell_realized_pqs_fixture
    )
    @test be2_inventory_timed.time >= 0.0
    @test be2_inventory_timed.bytes >= 0
    @test be2_current_route_inventory.object_kind ==
          :pqs_current_route_retained_unit_inventory_fixture
    @test be2_current_route_inventory.status == :private_diagnostic_only
    @test length(be2_current_route_inventory.units) == 8
    @test be2_current_route_inventory.diagnostics.unit_count == 8
    @test be2_current_route_inventory.diagnostics.shared_pqs_unit_count == 3
    @test be2_current_route_inventory.diagnostics.shared_pqs_roles == (
        :regular_shared_molecular_shell_1,
        :regular_shared_molecular_shell_2,
        :regular_shared_molecular_shell_3,
    )
    @test be2_current_route_inventory.diagnostics.shared_pqs_original_roles ==
          (:regular_shared_molecular_shell, :regular_shared_molecular_shell, :regular_shared_molecular_shell)
    @test map(unit -> unit.role, be2_shared_pqs_units) ==
          be2_current_route_inventory.diagnostics.shared_pqs_roles
    @test map(unit -> unit.original_role, be2_shared_pqs_units) ==
          be2_current_route_inventory.diagnostics.shared_pqs_original_roles
    @test map(unit -> unit.retained_count, be2_shared_pqs_units) == (98, 98, 114)
    @test map(unit -> unit.support_count, be2_shared_pqs_units) ==
          (1738, 1346, 1002)
    @test map(unit -> unit.column_range, be2_shared_pqs_units) ==
          (1174:1271, 1272:1369, 1370:1483)
    @test map(unit -> unit.original_column_range, be2_shared_pqs_units) ==
          (1174:1271, 1272:1369, 1370:1483)
    @test be2_current_route_inventory.diagnostics.shared_pqs_column_ranges ==
          (1174:1271, 1272:1369, 1370:1483)
    @test be2_current_route_inventory.diagnostics.shared_pqs_original_column_ranges ==
          (1174:1271, 1272:1369, 1370:1483)
    @test be2_current_route_inventory.diagnostics.shared_pqs_shell_realization_transform_fact_count == 3
    @test be2_current_route_inventory.diagnostics.shared_pqs_source_box_operator_application_ready_count == 0
    @test be2_current_route_inventory.coverage.first_column == 1
    @test be2_current_route_inventory.coverage.last_column == 1483
    @test be2_current_route_inventory.coverage.represented_count == 1483
    @test be2_current_route_inventory.coverage.covers_every_column_once
    @test be2_current_route_inventory.diagnostics.fixed_dimension == 1483
    @test be2_current_route_inventory.diagnostics.coverage_complete
    @test !be2_current_route_inventory.diagnostics.q4_single_shared_role_order_preserved
    @test all(
        unit -> unit.active_representation_stage == :shell_realized_pqs_fixture,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> unit.safe_term_capability == :support_local_oracle_for_shell_realization,
        be2_shared_pqs_units,
    )
    @test all(unit -> !unit.raw_product_box_operator_contract, be2_shared_pqs_units)
    @test all(
        unit -> unit.raw_box_auxiliary_metadata.available,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> unit.raw_box_auxiliary_metadata.reference_only,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> !unit.raw_box_auxiliary_metadata.active_current_route_contract,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> unit.diagnostics.representation_stage == :shell_realized_pqs_fixture,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> unit.diagnostics.shell_projection_lowdin_realization,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> !unit.diagnostics.raw_product_box_operator_contract,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> !unit.diagnostics.ida_weight_division_allowed,
        be2_shared_pqs_units,
    )
    @test all(
        unit -> unit.diagnostics.retained_weight_semantics ==
                :not_positive_quadrature_weights,
        be2_shared_pqs_units,
    )
    @test length(be2_current_route_inventory.source_fixtures.shared_pqs) == 3
    @test be2_current_route_pair_inventory.object_kind ==
          :pqs_current_route_retained_pair_inventory_fixture
    @test be2_current_route_pair_inventory.unit_inventory ===
          be2_current_route_inventory
    @test length(be2_current_route_pair_inventory.pairs) == 36
    @test be2_current_route_pair_inventory.counts.pair_count == 36
    @test be2_current_route_pair_inventory.diagnostics.pair_count == 36
    @test be2_current_route_pair_inventory.diagnostics.expected_pair_count == 36
    @test be2_current_route_pair_inventory.diagnostics.unit_count == 8
    @test be2_current_route_pair_inventory.counts.raw_box_pqs_active == 0
    @test be2_current_route_pair_inventory.diagnostics.raw_box_pqs_active_pair_policy_count == 0
    @test be2_current_route_pair_inventory.counts.support_local_oracle_for_shell_realization == 21
    @test be2_current_route_pair_inventory.diagnostics.support_local_oracle_for_shell_realization_pair_count == 21
    @test be2_current_route_pair_inventory.diagnostics.shell_row_oracle_pair_count == 21
    @test be2_current_route_pair_inventory.diagnostics.shell_realized_pqs_pairs_are_oracle_only
    @test all(
        pair -> !pair.raw_box_pqs_active_pair_policy,
        be2_current_route_pair_inventory.pairs,
    )
    @test count(pair -> pair.shell_row_oracle_only, be2_current_route_pair_inventory.pairs) == 21
    @test count(pair -> pair.active_algorithmic_policy, be2_current_route_pair_inventory.pairs) == 15
    @test !be2_current_route_pair_inventory.diagnostics.whole_route_safe_term_matrix_consumer
    current_route_safe_terms = CCPM._pqs_current_route_safe_term_matrices(
        pqs_construction,
        contact_safe_term_metrics;
        inventory = current_route_inventory,
        pair_inventory = current_route_pair_inventory,
    )
    @test current_route_safe_terms.object_kind ==
          :pqs_current_route_safe_term_matrices_fixture
    @test current_route_safe_terms.status == :private_diagnostic_only
    @test current_route_safe_terms.terms == (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :kinetic,
    )
    @test current_route_safe_terms.global_max_error <= 1.0e-12
    for term in current_route_safe_terms.terms
        @test size(current_route_safe_terms.matrices[term]) == (487, 487)
        @test size(current_route_safe_terms.oracle_matrices[term]) == (487, 487)
        @test all(isfinite, current_route_safe_terms.matrices[term])
        @test all(isfinite, current_route_safe_terms.oracle_matrices[term])
        @test current_route_safe_terms.matrices[term] ≈
              current_route_safe_terms.oracle_matrices[term] atol = 1.0e-12 rtol = 0.0
        @test current_route_safe_terms.term_errors[term] <= 1.0e-12
    end
    @test current_route_safe_terms.diagnostics.private_diagnostic_only
    @test current_route_safe_terms.diagnostics.whole_route_safe_term_matrix_consumer
    @test current_route_safe_terms.diagnostics.retained_dimension == 487
    @test current_route_safe_terms.diagnostics.pair_count == 21
    @test current_route_safe_terms.diagnostics.product_source_box_pair_count == 6
    @test current_route_safe_terms.diagnostics.support_local_fallback_pair_count == 15
    @test current_route_safe_terms.diagnostics.support_local_oracle_for_shell_realization_pair_count == 6
    @test current_route_safe_terms.diagnostics.support_local_oracle_is_debug_validation
    @test current_route_safe_terms.diagnostics.shell_realized_pqs_pairs_use_oracle_not_algorithm
    @test !current_route_safe_terms.diagnostics.source_box_algorithm_available_for_shell_realized_pqs
    @test current_route_safe_terms.diagnostics.raw_box_pqs_active_pair_policy_count == 0
    @test current_route_safe_terms.diagnostics.global_max_error <= 1.0e-12
    @test current_route_safe_terms.diagnostics.finite_output
    @test current_route_safe_terms.diagnostics.support_local_oracle_compared
    @test !current_route_safe_terms.diagnostics.raw_box_pqs_active_policy_used
    @test !current_route_safe_terms.diagnostics.route_descriptor_emitted
    @test !current_route_safe_terms.diagnostics.construction_mutated
    @test !current_route_safe_terms.diagnostics.sidecar_installation
    @test !current_route_safe_terms.diagnostics.packet_adoption
    @test !current_route_safe_terms.diagnostics.fixed_block_construction_changed
    @test !current_route_safe_terms.diagnostics.qwhamiltonian_changed
    @test !current_route_safe_terms.diagnostics.ida_weight_division_allowed
    @test current_route_safe_terms.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !current_route_safe_terms.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test current_route_safe_terms.diagnostics.elapsed_seconds >= 0.0
    @test current_route_safe_terms.diagnostics.allocated_bytes >= 0
    @test_throws ArgumentError CCPM._pqs_current_route_safe_term_matrices(
        pqs_construction,
        contact_safe_term_metrics;
        inventory = current_route_inventory,
        pair_inventory = current_route_pair_inventory,
        terms = (:weights,),
    )
    @test :product_doside_unit in pqs_source_descriptor.non_contracts
    @test :dense_full_parent_fallback in pqs_source_descriptor.non_contracts
    @test pqs_source_descriptor.diagnostics.metadata_only
    @test !pqs_source_descriptor.active_consumption.fixed_block_sidecar_installed
    @test !pqs_source_descriptor.active_consumption.metric_packet_consumes
    @test !pqs_source_descriptor.active_consumption.by_center_consumes
    @test pqs_diagnostics.region_builds[2].retained_count == 114
    @test pqs_diagnostics.fixed_dimension == 487
    @test !isnothing(pqs_construction.sequence.packet)
    @test all(isfinite, pqs_construction.sequence.packet.overlap)
    @test all(isfinite, pqs_construction.sequence.packet.kinetic)
    @test all(isfinite, pqs_construction.sequence.packet.weights)
    @test all(isfinite, pqs_construction.sequence.packet.gaussian_sum)
    @test all(isfinite, pqs_construction.sequence.packet.pair_sum)
    @test norm(pqs_construction.sequence.packet.overlap - I, Inf) < 1.0e-8
    pqs_fixed_block =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_fixed_block(
            pqs_construction,
        )
    @test size(pqs_fixed_block.coefficient_matrix) == (7 * 7 * 15, 487)
    @test all(isfinite, pqs_fixed_block.overlap)
    @test all(isfinite, pqs_fixed_block.kinetic)
    @test all(isfinite, pqs_fixed_block.weights)
    @test all(isfinite, pqs_fixed_block.gaussian_sum)
    @test all(isfinite, pqs_fixed_block.pair_sum)
    @test norm(pqs_fixed_block.overlap - I, Inf) < 1.0e-8
    @test pqs_fixed_block.staged_by_center_sidecar[] === nothing
    current_route_authority_comparison =
        CCPM._pqs_current_route_safe_term_authority_comparison(
            pqs_construction,
            contact_safe_term_metrics;
            inventory = current_route_inventory,
            pair_inventory = current_route_pair_inventory,
            safe_terms = current_route_safe_terms,
            fixed_block = pqs_fixed_block,
        )
    @test current_route_authority_comparison.object_kind ==
          :pqs_current_route_safe_term_authority_comparison_fixture
    @test current_route_authority_comparison.status == :private_diagnostic_only
    @test current_route_authority_comparison.terms == current_route_safe_terms.terms
    @test current_route_authority_comparison.compared_terms ==
          current_route_safe_terms.terms
    @test isempty(current_route_authority_comparison.unavailable_terms)
    @test current_route_authority_comparison.max_authority_error <= 1.0e-8
    for term in current_route_authority_comparison.compared_terms
        @test current_route_authority_comparison.term_errors[term] <= 1.0e-8
        @test current_route_authority_comparison.authority_sources[term] ==
              :fixed_block
        @test current_route_authority_comparison.authority_fields[term] == term
        @test current_route_authority_comparison.authority_shapes[term] == (487, 487)
    end
    @test current_route_authority_comparison.diagnostics.private_diagnostic_only
    @test current_route_authority_comparison.diagnostics.current_route_safe_term_authority_comparison
    @test current_route_authority_comparison.diagnostics.authority_fixed_block_or_sequence_packet_only
    @test current_route_authority_comparison.diagnostics.support_local_oracle_secondary
    @test current_route_authority_comparison.diagnostics.support_local_oracle_global_max_error <=
          1.0e-12
    @test current_route_authority_comparison.diagnostics.compared_term_count ==
          length(current_route_safe_terms.terms)
    @test current_route_authority_comparison.diagnostics.unavailable_term_count == 0
    @test current_route_authority_comparison.diagnostics.max_authority_error <= 1.0e-8
    @test current_route_authority_comparison.diagnostics.finite_output
    @test !current_route_authority_comparison.diagnostics.construction_mutated
    @test !current_route_authority_comparison.diagnostics.sidecar_installation
    @test !current_route_authority_comparison.diagnostics.packet_adoption
    @test !current_route_authority_comparison.diagnostics.fixed_block_construction_changed
    @test !current_route_authority_comparison.diagnostics.qwhamiltonian_changed
    @test !current_route_authority_comparison.diagnostics.ida_weight_division_allowed
    @test current_route_authority_comparison.diagnostics.retained_weight_semantics ==
          :not_positive_quadrature_weights
    @test !current_route_authority_comparison.diagnostics.local_ecp_gaussian_mwg_interaction_changed
    @test current_route_authority_comparison.diagnostics.elapsed_seconds >= 0.0
    @test current_route_authority_comparison.diagnostics.allocated_bytes >= 0
    QWCS = GaussletBases.CartesianQWOperatorCarriedSpaces
    pqs_receipt = @test_logs min_level = Logging.Warn QWCS.cartesian_qw_operator_construction_receipt(
        pqs_fixed_block;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
        expansion = expansion,
    )
    pqs_operators = QWCS.qw_operator_construction_receipt_operators(pqs_receipt)
    pqs_receipt_diagnostics = QWCS.qw_operator_construction_receipt_diagnostics(
        pqs_receipt,
    )
    pqs_record_diagnostics = QWCS.qw_operator_construction_record_diagnostics(
        QWCS.qw_operator_construction_receipt_record(pqs_receipt),
    )
    @test pqs_receipt_diagnostics.delegated_to_existing_builder
    @test pqs_receipt_diagnostics.builder == :ordinary_cartesian_qiu_white_operators
    @test pqs_receipt_diagnostics.source_sidecar_agree
    @test isempty(pqs_receipt_diagnostics.mismatch_fields)
    @test pqs_receipt_diagnostics.operator_built
    @test pqs_receipt_diagnostics.gausslet_backend == :pgdg_localized_experimental
    @test pqs_receipt_diagnostics.interaction_treatment == :ggt_nearest
    @test pqs_receipt_diagnostics.nuclear_term_storage == :total_only
    @test pqs_receipt_diagnostics.dense_parent_matrix_used == false
    @test pqs_receipt_diagnostics.heavy_metric_packet_built == false
    @test pqs_receipt_diagnostics.new_hamiltonian_kernel_used == false
    @test pqs_receipt_diagnostics.numerical_outputs_changed == false
    @test pqs_receipt_diagnostics.operator_residual_count == 0
    @test pqs_record_diagnostics.source_sidecar_agree
    @test :carried_has_staged_sidecar in pqs_record_diagnostics.compared_fields
    @test !pqs_record_diagnostics.comparisons.carried_has_staged_sidecar.source_value
    @test !pqs_record_diagnostics.comparisons.carried_has_staged_sidecar.operator_value
    @test isnothing(
        pqs_record_diagnostics.comparisons.carried_staged_by_center_path.source_value,
    )
    @test isnothing(
        pqs_record_diagnostics.comparisons.carried_staged_by_center_path.operator_value,
    )
    @test pqs_record_diagnostics.source_parent_dimension == 7 * 7 * 15
    @test pqs_record_diagnostics.sidecar_parent_dimension == 7 * 7 * 15
    @test pqs_record_diagnostics.source_carried_dimension == 487
    @test pqs_record_diagnostics.sidecar_carried_dimension == 487
    @test pqs_record_diagnostics.source_carried_space_kind == :nested_fixed_block
    @test pqs_record_diagnostics.sidecar_input_kind == :nested_fixed_block_operator
    @test pqs_operators.gausslet_backend == :pgdg_localized_experimental
    @test pqs_operators.interaction_treatment == :ggt_nearest
    @test pqs_operators.nuclear_term_storage == :total_only
    @test pqs_operators.gausslet_count == 487
    @test pqs_operators.residual_count == 0
    @test size(pqs_operators.overlap) == (487, 487)
    @test size(pqs_operators.one_body_hamiltonian) == (487, 487)
    @test size(pqs_operators.interaction_matrix) == (487, 487)
    @test all(isfinite, pqs_operators.overlap)
    @test all(isfinite, pqs_operators.one_body_hamiltonian)
    @test all(isfinite, pqs_operators.interaction_matrix)
    @test norm(pqs_operators.overlap - I, Inf) < 1.0e-8
    @test norm(
        pqs_operators.one_body_hamiltonian - transpose(pqs_operators.one_body_hamiltonian),
        Inf,
    ) < 1.0e-12
    @test norm(
        pqs_operators.interaction_matrix - transpose(pqs_operators.interaction_matrix),
        Inf,
    ) < 1.0e-12
    @test_throws ArgumentError GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
        basis,
        bundles,
        policy;
        nside = 5,
        term_coefficients = Float64.(expansion.coefficients),
        packet_kernel = :factorized_direct,
        shared_shell_realization = :projected_q_shell,
    )

    shared_q5_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        policy.construction_plan;
        q_min = 4,
        atom_q = 4,
        atom_order = 4,
        shared_q = 5,
        shared_order = 5,
        contact_q = 4,
        contact_order = 4,
        outer_mismatch_q = 4,
        outer_mismatch_order = 4,
    )
    shared_q5_construction =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            shared_q5_policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :factorized_direct,
            build_sequence_packet = false,
        )
    shared_q5_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction_diagnostics(
            shared_q5_construction,
        )
    shared_q5_region_builds = shared_q5_diagnostics.region_builds
    @test shared_q5_diagnostics.recipe_label == :mixed_atom_cubic_shared_endcap_panel
    @test shared_q5_diagnostics.active_builder_consumes
    @test shared_q5_diagnostics.parent_dimension == 7 * 7 * 15
    @test shared_q5_diagnostics.fixed_dimension == 523
    @test shared_q5_diagnostics.support_coverage.coverage_ok
    @test [build.q for build in shared_q5_region_builds] == [4, 5, 4, 4, 4]
    @test [build.order for build in shared_q5_region_builds] == [4, 5, 4, 4, 4]
    @test [build.retained_count for build in shared_q5_region_builds] == [98, 150, 125, 125, 25]
    @test [build.column_range for build in shared_q5_region_builds] ==
          [1:98, 374:523, 99:223, 224:348, 349:373]
    @test shared_q5_region_builds[2].metadata.coefficient_contract == :product_doside
    @test shared_q5_region_builds[5].metadata.descriptor_scope == :middle_contact_cap
    @test shared_q5_diagnostics.metadata.q_policy == :atom_growth_endcap_panel_shared_q_variable
    @test !shared_q5_diagnostics.metadata.q4_acceptance_fixture
    @test shared_q5_diagnostics.metadata.non_shared_q_policy == :fixed_q4_order4
    @test shared_q5_diagnostics.metadata.shared_q_values == (5,)
    @test shared_q5_diagnostics.metadata.shared_order_values == (5,)
    @test isnothing(shared_q5_construction.sequence.packet)

    no_packet_readiness =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness(
            construction,
        )
    no_packet_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness_diagnostics(
            no_packet_readiness,
        )
    @test no_packet_diagnostics.parent_dimension == 7 * 7 * 15
    @test no_packet_diagnostics.fixed_dimension == 469
    @test no_packet_diagnostics.support_coverage.coverage_ok
    @test !no_packet_diagnostics.sequence_packet_available
    @test !no_packet_diagnostics.can_produce_fixed_block
    @test no_packet_diagnostics.fixed_block_missing_fields == [:sequence_packet]
    @test !no_packet_diagnostics.can_produce_nested_source
    @test :split_geometry in no_packet_diagnostics.nested_source_missing_fields
    @test no_packet_diagnostics.default_builders_unchanged

    packeted_construction =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :factorized_direct,
            build_sequence_packet = true,
        )
    readiness =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness(
            packeted_construction;
            build_fixed_block = true,
        )
    readiness_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness_diagnostics(
            readiness,
        )
    @test readiness isa
          GaussletBases._BondAlignedDiatomicHighOrderRecipeSourceReadiness3D
    @test readiness_diagnostics.sequence_packet_available
    @test readiness_diagnostics.overlap_available
    @test readiness_diagnostics.weights_available
    @test readiness_diagnostics.overlap_error < 1.0e-8
    @test readiness_diagnostics.can_produce_fixed_block
    @test isempty(readiness_diagnostics.fixed_block_missing_fields)
    @test readiness_diagnostics.fixed_block_built
    @test readiness_diagnostics.fixed_block_backend == :unknown
    @test readiness_diagnostics.fixed_block_dimension == 469
    @test readiness_diagnostics.fixed_block_support_count == 7 * 7 * 15
    @test readiness_diagnostics.staged_by_center_sidecar_available
    @test !readiness_diagnostics.can_produce_nested_source
    @test :child_sequences in readiness_diagnostics.nested_source_missing_fields
    fixed_block =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_fixed_block(
            packeted_construction,
        )
    @test fixed_block isa GaussletBases._NestedFixedBlock3D
    @test size(fixed_block.coefficient_matrix) == (7 * 7 * 15, 469)
    @test norm(fixed_block.overlap - I, Inf) < 1.0e-8
    @test all(isfinite, fixed_block.weights)
    @test fixed_block.staged_by_center_sidecar[] isa
          GaussletBases._CartesianNestedProductStagedByCenterSidecar3D

    receipt = @test_logs min_level = Logging.Warn QWCS.cartesian_qw_operator_construction_receipt(
        fixed_block;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
        expansion = expansion,
    )
    operators = QWCS.qw_operator_construction_receipt_operators(receipt)
    receipt_diagnostics = QWCS.qw_operator_construction_receipt_diagnostics(receipt)
    record_diagnostics = QWCS.qw_operator_construction_record_diagnostics(
        QWCS.qw_operator_construction_receipt_record(receipt),
    )
    @test receipt_diagnostics.delegated_to_existing_builder
    @test receipt_diagnostics.builder == :ordinary_cartesian_qiu_white_operators
    @test receipt_diagnostics.source_sidecar_agree
    @test isempty(receipt_diagnostics.mismatch_fields)
    @test receipt_diagnostics.operator_built
    @test receipt_diagnostics.gausslet_backend == :pgdg_localized_experimental
    @test receipt_diagnostics.interaction_treatment == :ggt_nearest
    @test receipt_diagnostics.nuclear_term_storage == :total_only
    @test receipt_diagnostics.dense_parent_matrix_used == false
    @test receipt_diagnostics.heavy_metric_packet_built == false
    @test receipt_diagnostics.new_hamiltonian_kernel_used == false
    @test receipt_diagnostics.numerical_outputs_changed == false
    @test record_diagnostics.source_sidecar_agree
    @test record_diagnostics.source_parent_dimension == 7 * 7 * 15
    @test record_diagnostics.sidecar_parent_dimension == 7 * 7 * 15
    @test record_diagnostics.source_carried_dimension == 469
    @test record_diagnostics.sidecar_carried_dimension == 469
    @test record_diagnostics.source_carried_space_kind == :nested_fixed_block
    @test record_diagnostics.sidecar_input_kind == :nested_fixed_block_operator
    @test operators.gausslet_backend == :pgdg_localized_experimental
    @test operators.interaction_treatment == :ggt_nearest
    @test operators.gausslet_count == 469
    @test operators.residual_count == 0
    @test size(operators.overlap) == (469, 469)
    @test size(operators.one_body_hamiltonian) == (469, 469)
    @test size(operators.interaction_matrix) == (469, 469)
    @test all(isfinite, operators.overlap)
    @test all(isfinite, operators.one_body_hamiltonian)
    @test all(isfinite, operators.interaction_matrix)
    @test norm(operators.overlap - I, Inf) < 1.0e-8
    @test norm(
        operators.one_body_hamiltonian - transpose(operators.one_body_hamiltonian),
        Inf,
    ) < 1.0e-12
    @test norm(
        operators.interaction_matrix - transpose(operators.interaction_matrix),
        Inf,
    ) < 1.0e-12

    shared_q5_packeted_construction =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            shared_q5_policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :factorized_direct,
            build_sequence_packet = true,
        )
    shared_q5_readiness =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness(
            shared_q5_packeted_construction;
            build_fixed_block = true,
        )
    shared_q5_readiness_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness_diagnostics(
            shared_q5_readiness,
        )
    @test shared_q5_readiness_diagnostics.sequence_packet_available
    @test shared_q5_readiness_diagnostics.overlap_available
    @test shared_q5_readiness_diagnostics.weights_available
    @test shared_q5_readiness_diagnostics.overlap_error < 1.0e-8
    @test shared_q5_readiness_diagnostics.can_produce_fixed_block
    @test isempty(shared_q5_readiness_diagnostics.fixed_block_missing_fields)
    @test shared_q5_readiness_diagnostics.fixed_block_built
    @test shared_q5_readiness_diagnostics.fixed_block_dimension == 523
    @test shared_q5_readiness_diagnostics.fixed_block_support_count == 7 * 7 * 15
    @test shared_q5_readiness_diagnostics.staged_by_center_sidecar_available
    shared_q5_fixed_block =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_fixed_block(
            shared_q5_packeted_construction,
        )
    @test shared_q5_fixed_block isa GaussletBases._NestedFixedBlock3D
    @test size(shared_q5_fixed_block.coefficient_matrix) == (7 * 7 * 15, 523)
    @test norm(shared_q5_fixed_block.overlap - I, Inf) < 1.0e-8
    @test all(isfinite, shared_q5_fixed_block.weights)
    @test shared_q5_fixed_block.staged_by_center_sidecar[] isa
          GaussletBases._CartesianNestedProductStagedByCenterSidecar3D

    shared_q5_receipt = @test_logs min_level = Logging.Warn QWCS.cartesian_qw_operator_construction_receipt(
        shared_q5_fixed_block;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
        expansion = expansion,
    )
    shared_q5_operators = QWCS.qw_operator_construction_receipt_operators(shared_q5_receipt)
    shared_q5_receipt_diagnostics =
        QWCS.qw_operator_construction_receipt_diagnostics(shared_q5_receipt)
    shared_q5_record_diagnostics = QWCS.qw_operator_construction_record_diagnostics(
        QWCS.qw_operator_construction_receipt_record(shared_q5_receipt),
    )
    @test shared_q5_receipt_diagnostics.delegated_to_existing_builder
    @test shared_q5_receipt_diagnostics.source_sidecar_agree
    @test isempty(shared_q5_receipt_diagnostics.mismatch_fields)
    @test shared_q5_receipt_diagnostics.operator_built
    @test shared_q5_receipt_diagnostics.gausslet_backend == :pgdg_localized_experimental
    @test shared_q5_receipt_diagnostics.interaction_treatment == :ggt_nearest
    @test shared_q5_receipt_diagnostics.nuclear_term_storage == :total_only
    @test shared_q5_receipt_diagnostics.dense_parent_matrix_used == false
    @test shared_q5_receipt_diagnostics.heavy_metric_packet_built == false
    @test shared_q5_receipt_diagnostics.new_hamiltonian_kernel_used == false
    @test shared_q5_receipt_diagnostics.numerical_outputs_changed == false
    @test shared_q5_record_diagnostics.source_sidecar_agree
    @test isempty(shared_q5_record_diagnostics.mismatch_fields)
    @test shared_q5_record_diagnostics.source_parent_dimension == 7 * 7 * 15
    @test shared_q5_record_diagnostics.sidecar_parent_dimension == 7 * 7 * 15
    @test shared_q5_record_diagnostics.source_carried_dimension == 523
    @test shared_q5_record_diagnostics.sidecar_carried_dimension == 523
    @test shared_q5_operators.gausslet_backend == :pgdg_localized_experimental
    @test shared_q5_operators.interaction_treatment == :ggt_nearest
    @test shared_q5_operators.gausslet_count == 523
    @test shared_q5_operators.residual_count == 0
    @test size(shared_q5_operators.overlap) == (523, 523)
    @test size(shared_q5_operators.one_body_hamiltonian) == (523, 523)
    @test size(shared_q5_operators.interaction_matrix) == (523, 523)
    @test all(isfinite, shared_q5_operators.overlap)
    @test all(isfinite, shared_q5_operators.one_body_hamiltonian)
    @test all(isfinite, shared_q5_operators.interaction_matrix)
    @test norm(shared_q5_operators.overlap - I, Inf) < 1.0e-8
    @test norm(
        shared_q5_operators.one_body_hamiltonian -
        transpose(shared_q5_operators.one_body_hamiltonian),
        Inf,
    ) < 1.0e-12
    @test norm(
        shared_q5_operators.interaction_matrix -
        transpose(shared_q5_operators.interaction_matrix),
        Inf,
    ) < 1.0e-12

    shared_q6_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        policy.construction_plan;
        q_min = 4,
        atom_q = 4,
        atom_order = 4,
        shared_q = 6,
        shared_order = 6,
        contact_q = 4,
        contact_order = 4,
        outer_mismatch_q = 4,
        outer_mismatch_order = 4,
    )
    shared_q6_packeted_construction =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            shared_q6_policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :factorized_direct,
            build_sequence_packet = true,
        )
    shared_q6_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction_diagnostics(
            shared_q6_packeted_construction,
        )
    shared_q6_region_builds = shared_q6_diagnostics.region_builds
    @test shared_q6_diagnostics.recipe_label == :mixed_atom_cubic_shared_endcap_panel
    @test shared_q6_diagnostics.active_builder_consumes
    @test shared_q6_diagnostics.parent_dimension == 7 * 7 * 15
    @test shared_q6_diagnostics.fixed_dimension == 589
    @test shared_q6_diagnostics.support_coverage.coverage_ok
    @test [build.q for build in shared_q6_region_builds] == [4, 6, 4, 4, 4]
    @test [build.order for build in shared_q6_region_builds] == [4, 6, 4, 4, 4]
    @test [build.retained_count for build in shared_q6_region_builds] ==
          [98, 216, 125, 125, 25]
    @test [build.column_range for build in shared_q6_region_builds] ==
          [1:98, 374:589, 99:223, 224:348, 349:373]
    @test shared_q6_region_builds[2].metadata.coefficient_contract == :product_doside
    @test all(
        build.retained_count == build.built_support_count for
        build in shared_q6_region_builds[[1, 3, 4, 5]]
    )
    @test shared_q6_diagnostics.metadata.non_shared_q_policy == :fixed_q4_order4
    @test shared_q6_diagnostics.metadata.shared_q_values == (6,)
    @test shared_q6_diagnostics.metadata.shared_order_values == (6,)

    shared_q6_readiness =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness(
            shared_q6_packeted_construction;
            build_fixed_block = true,
        )
    shared_q6_readiness_diagnostics =
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_readiness_diagnostics(
            shared_q6_readiness,
        )
    @test shared_q6_readiness_diagnostics.sequence_packet_available
    @test shared_q6_readiness_diagnostics.overlap_available
    @test shared_q6_readiness_diagnostics.weights_available
    @test shared_q6_readiness_diagnostics.overlap_error < 1.0e-8
    @test shared_q6_readiness_diagnostics.can_produce_fixed_block
    @test isempty(shared_q6_readiness_diagnostics.fixed_block_missing_fields)
    @test shared_q6_readiness_diagnostics.fixed_block_built
    @test shared_q6_readiness_diagnostics.fixed_block_dimension == 589
    @test shared_q6_readiness_diagnostics.fixed_block_support_count == 7 * 7 * 15
    @test shared_q6_readiness_diagnostics.staged_by_center_sidecar_available
    shared_q6_fixed_block = shared_q6_readiness.fixed_block
    @test shared_q6_fixed_block isa GaussletBases._NestedFixedBlock3D
    @test size(shared_q6_fixed_block.coefficient_matrix) == (7 * 7 * 15, 589)
    @test norm(shared_q6_fixed_block.overlap - I, Inf) < 1.0e-8
    @test all(isfinite, shared_q6_fixed_block.weights)
    @test shared_q6_fixed_block.staged_by_center_sidecar[] isa
          GaussletBases._CartesianNestedProductStagedByCenterSidecar3D

    shared_q6_receipt = @test_logs min_level = Logging.Warn QWCS.cartesian_qw_operator_construction_receipt(
        shared_q6_fixed_block;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
        expansion = expansion,
    )
    shared_q6_operators = QWCS.qw_operator_construction_receipt_operators(shared_q6_receipt)
    shared_q6_receipt_diagnostics =
        QWCS.qw_operator_construction_receipt_diagnostics(shared_q6_receipt)
    shared_q6_record_diagnostics = QWCS.qw_operator_construction_record_diagnostics(
        QWCS.qw_operator_construction_receipt_record(shared_q6_receipt),
    )
    @test shared_q6_receipt_diagnostics.delegated_to_existing_builder
    @test shared_q6_receipt_diagnostics.source_sidecar_agree
    @test isempty(shared_q6_receipt_diagnostics.mismatch_fields)
    @test shared_q6_receipt_diagnostics.operator_built
    @test shared_q6_receipt_diagnostics.gausslet_backend == :pgdg_localized_experimental
    @test shared_q6_receipt_diagnostics.interaction_treatment == :ggt_nearest
    @test shared_q6_receipt_diagnostics.nuclear_term_storage == :total_only
    @test shared_q6_receipt_diagnostics.dense_parent_matrix_used == false
    @test shared_q6_receipt_diagnostics.heavy_metric_packet_built == false
    @test shared_q6_receipt_diagnostics.new_hamiltonian_kernel_used == false
    @test shared_q6_receipt_diagnostics.numerical_outputs_changed == false
    @test shared_q6_record_diagnostics.source_sidecar_agree
    @test isempty(shared_q6_record_diagnostics.mismatch_fields)
    @test shared_q6_record_diagnostics.source_parent_dimension == 7 * 7 * 15
    @test shared_q6_record_diagnostics.sidecar_parent_dimension == 7 * 7 * 15
    @test shared_q6_record_diagnostics.source_carried_dimension == 589
    @test shared_q6_record_diagnostics.sidecar_carried_dimension == 589
    @test shared_q6_operators.gausslet_backend == :pgdg_localized_experimental
    @test shared_q6_operators.interaction_treatment == :ggt_nearest
    @test shared_q6_operators.gausslet_count == 589
    @test shared_q6_operators.residual_count == 0
    @test size(shared_q6_operators.overlap) == (589, 589)
    @test size(shared_q6_operators.one_body_hamiltonian) == (589, 589)
    @test size(shared_q6_operators.interaction_matrix) == (589, 589)
    @test all(isfinite, shared_q6_operators.overlap)
    @test all(isfinite, shared_q6_operators.one_body_hamiltonian)
    @test all(isfinite, shared_q6_operators.interaction_matrix)
    @test norm(shared_q6_operators.overlap - I, Inf) < 1.0e-8
    @test norm(
        shared_q6_operators.one_body_hamiltonian -
        transpose(shared_q6_operators.one_body_hamiltonian),
        Inf,
    ) < 1.0e-12
    @test norm(
        shared_q6_operators.interaction_matrix -
        transpose(shared_q6_operators.interaction_matrix),
        Inf,
    ) < 1.0e-12

    shared_q7_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        policy.construction_plan;
        q_min = 4,
        atom_q = 4,
        atom_order = 4,
        shared_q = 7,
        shared_order = 7,
        contact_q = 4,
        contact_order = 4,
        outer_mismatch_q = 4,
        outer_mismatch_order = 4,
    )
    shared_q7_error = try
        GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
            basis,
            bundles,
            shared_q7_policy;
            nside = 5,
            term_coefficients = Float64.(expansion.coefficients),
            packet_kernel = :factorized_direct,
            build_sequence_packet = false,
        )
        nothing
    catch err
        err
    end
    @test shared_q7_error isa ArgumentError
    @test occursin(
        "nested doside retained_count must not exceed the interval size",
        sprint(showerror, shared_q7_error),
    )

    q_row_expectations = (
        (q = 4, fixed_dimension = 469, retained = (98, 96, 125, 125, 25), shared = 96),
        (q = 5, fixed_dimension = 523, retained = (98, 150, 125, 125, 25), shared = 150),
        (q = 6, fixed_dimension = 589, retained = (98, 216, 125, 125, 25), shared = 216),
    )
    multi_shared_region_builds = (
        (role = :outer_mismatch_shared_molecular_shell, retained_count = 12),
        (role = :regular_shared_molecular_shell, retained_count = 24),
        (role = :regular_shared_molecular_shell, retained_count = 36),
        (role = :left_atom_box, retained_count = 64),
    )
    @test QWCS._nested_q_row_shared_retained_counts(multi_shared_region_builds) ==
          (24, 36)
    @test QWCS._nested_q_row_shared_retained_count(multi_shared_region_builds) ===
          nothing
    @test QWCS._nested_q_row_shared_retained_counts(()) == ()
    @test QWCS._nested_q_row_shared_retained_count(()) === nothing
    for expected in q_row_expectations
        q_row_receipt = @test_logs min_level = Logging.Warn QWCS._nested_bond_aligned_diatomic_high_order_q_row_route_receipt(
            basis;
            shared_q = expected.q,
            shared_order = expected.q,
            protected_atom_side_count = 5,
            q_min = 4,
            nside = 5,
            expansion = expansion,
            nuclear_charges = [1.0, 1.0],
            nuclear_term_storage = :total_only,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :pgdg_localized_experimental,
        )
        q_row_diagnostics =
            QWCS._nested_bond_aligned_diatomic_high_order_q_row_route_diagnostics(
                q_row_receipt,
            )
        @test q_row_diagnostics.route_label ==
              :bond_aligned_diatomic_high_order_q_row_route
        @test q_row_diagnostics.shared_q == expected.q
        @test q_row_diagnostics.shared_order == expected.q
        @test q_row_diagnostics.non_shared_q_policy == :fixed_q4_order4
        @test q_row_diagnostics.parent_dimension == 7 * 7 * 15
        @test q_row_diagnostics.fixed_dimension == expected.fixed_dimension
        @test q_row_diagnostics.retained_counts_by_region == expected.retained
        @test q_row_diagnostics.shared_retained_counts == (expected.shared,)
        @test q_row_diagnostics.shared_retained_count == expected.shared
        @test q_row_diagnostics.overlap_error < 1.0e-8
        @test q_row_diagnostics.staged_sidecar_available
        @test q_row_diagnostics.backend == :pgdg_localized_experimental
        @test q_row_diagnostics.residual_count == 0
        @test q_row_diagnostics.gausslet_count == expected.fixed_dimension
        @test q_row_diagnostics.source_sidecar_agree
        @test isempty(q_row_diagnostics.mismatch_fields)
        @test q_row_diagnostics.dense_parent_matrix_used == false
        @test q_row_diagnostics.heavy_metric_packet_built == false
        @test q_row_diagnostics.new_hamiltonian_kernel_used == false
        @test q_row_diagnostics.default_source_builder_changed == false
    end

    for expected in q_row_expectations
        fixture_receipt = @test_logs min_level = Logging.Warn QWCS._nested_bond_aligned_homonuclear_high_order_q_row_fixture_receipt(
            bond_length = 5.0,
            core_spacing = 0.7,
            xmax_parallel = 8.0,
            xmax_transverse = 4.0,
            shared_q = expected.q,
            shared_order = expected.q,
            family = :G10,
            bond_axis = :z,
            nuclear_charge = 1.0,
            reference_spacing = 1.0,
            tail_spacing = 10.0,
            protected_atom_side_count = 5,
            q_min = 4,
            nside = 5,
            expansion = expansion,
            nuclear_term_storage = :total_only,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :pgdg_localized_experimental,
        )
        fixture_diagnostics =
            QWCS._nested_bond_aligned_homonuclear_high_order_q_row_fixture_diagnostics(
                fixture_receipt,
            )
        fixture_provenance =
            QWCS._nested_bond_aligned_homonuclear_high_order_q_row_fixture_provenance(
                fixture_receipt,
            )
        route_diagnostics = fixture_diagnostics.q_row_route_diagnostics
        @test fixture_diagnostics.route_label ==
              :bond_aligned_homonuclear_high_order_q_row_fixture
        @test fixture_diagnostics.receipt_contract ==
              :construct_homonuclear_basis_then_delegate_q_row_route
        @test fixture_diagnostics.basis_constructor ==
              :bond_aligned_homonuclear_qw_basis
        @test fixture_diagnostics.family == :G10
        @test fixture_diagnostics.bond_length == 5.0
        @test fixture_diagnostics.core_spacing == 0.7
        @test fixture_diagnostics.xmax_parallel == 8.0
        @test fixture_diagnostics.xmax_transverse == 4.0
        @test fixture_diagnostics.bond_axis == :z
        @test fixture_diagnostics.reference_spacing == 1.0
        @test fixture_diagnostics.tail_spacing == 10.0
        @test fixture_diagnostics.nuclear_charge == 1.0
        @test fixture_diagnostics.basis_nuclear_charges == (1.0, 1.0)
        @test fixture_receipt.basis.nuclear_charges == [1.0, 1.0]
        @test fixture_diagnostics.basis_nuclei == ((0.0, 0.0, -2.5), (0.0, 0.0, 2.5))
        @test fixture_diagnostics.parent_axis_counts == (7, 7, 15)
        @test fixture_diagnostics.parent_dimension == 7 * 7 * 15
        @test fixture_diagnostics.flat_index_convention.one_based
        @test fixture_diagnostics.flat_index_convention.fastest_axis == :z
        @test fixture_diagnostics.flat_index_convention.slowest_axis == :x
        @test fixture_diagnostics.shared_q == expected.q
        @test fixture_diagnostics.shared_order == expected.q
        @test fixture_diagnostics.q_min == 4
        @test fixture_diagnostics.protected_atom_side_count == 5
        @test fixture_diagnostics.nside == 5
        @test fixture_diagnostics.fixed_dimension == expected.fixed_dimension
        @test fixture_diagnostics.retained_counts_by_region == expected.retained
        @test fixture_diagnostics.shared_retained_counts == (expected.shared,)
        @test fixture_diagnostics.shared_retained_count == expected.shared
        @test route_diagnostics.route_label ==
              :bond_aligned_diatomic_high_order_q_row_route
        @test route_diagnostics.parent_dimension == fixture_diagnostics.parent_dimension
        @test route_diagnostics.fixed_dimension == expected.fixed_dimension
        @test route_diagnostics.backend == :pgdg_localized_experimental
        @test route_diagnostics.source_sidecar_agree
        @test isempty(route_diagnostics.mismatch_fields)
        @test fixture_provenance.charge_policy == :basis_nuclear_charges_only
        @test fixture_provenance.homonuclear_only
        @test !fixture_provenance.heteronuclear_support
        @test !fixture_provenance.public_api
        @test !fixture_provenance.science_validation
    end

    @test_throws ArgumentError QWCS._nested_bond_aligned_diatomic_high_order_q_row_route_receipt(
        basis;
        shared_q = 4,
        expansion = expansion,
        gausslet_backend = :numerical_reference,
    )
    q_row_q7_error = try
        QWCS._nested_bond_aligned_diatomic_high_order_q_row_route_receipt(
            basis;
            shared_q = 7,
            shared_order = 7,
            expansion = expansion,
            gausslet_backend = :pgdg_localized_experimental,
        )
        nothing
    catch err
        err
    end
    @test q_row_q7_error isa ArgumentError
    @test occursin(
        "nested doside retained_count must not exceed the interval size",
        sprint(showerror, q_row_q7_error),
    )

    @test_throws ArgumentError QWCS._nested_bond_aligned_homonuclear_high_order_q_row_fixture_receipt(
        bond_length = 5.0,
        core_spacing = 0.7,
        xmax_parallel = 8.0,
        xmax_transverse = 4.0,
        shared_q = 4,
        expansion = expansion,
        gausslet_backend = :numerical_reference,
    )
    @test_throws UndefKeywordError QWCS._nested_bond_aligned_homonuclear_high_order_q_row_fixture_receipt(
        bond_length = 5.0,
        core_spacing = 0.7,
        xmax_transverse = 4.0,
        shared_q = 4,
        expansion = expansion,
    )
    fixture_q7_error = try
        QWCS._nested_bond_aligned_homonuclear_high_order_q_row_fixture_receipt(
            bond_length = 5.0,
            core_spacing = 0.7,
            xmax_parallel = 8.0,
            xmax_transverse = 4.0,
            shared_q = 7,
            shared_order = 7,
            expansion = expansion,
            gausslet_backend = :pgdg_localized_experimental,
        )
        nothing
    catch err
        err
    end
    @test fixture_q7_error isa ArgumentError
    @test occursin(
        "nested doside retained_count must not exceed the interval size",
        sprint(showerror, fixture_q7_error),
    )

    annulus_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        policy.construction_plan;
        q_min = 4,
        atom_q = 4,
        shared_q = 5,
        contact_q = 4,
        outer_mismatch_q = 4,
        shared_exterior_family = :transverse_annulus_exterior,
    )
    @test_throws ArgumentError GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
        basis,
        bundles,
        annulus_policy;
        nside = 5,
        term_coefficients = Float64.(expansion.coefficients),
        packet_kernel = :factorized_direct,
        build_sequence_packet = false,
    )

    atom_q5_policy = GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_policy(
        policy.construction_plan;
        q_min = 4,
        atom_q = 5,
        atom_order = 5,
        shared_q = 5,
        shared_order = 5,
        contact_q = 4,
        outer_mismatch_q = 4,
    )
    @test_throws ArgumentError GaussletBases._nested_bond_aligned_diatomic_high_order_recipe_source_construction(
        basis,
        bundles,
        atom_q5_policy;
        nside = 5,
        term_coefficients = Float64.(expansion.coefficients),
        packet_kernel = :factorized_direct,
        build_sequence_packet = false,
    )
end

@testset "Bond-aligned diatomic endcap-panel shared shell source policy" begin
    basis = bond_aligned_homonuclear_qw_basis(
        family = :G10,
        bond_length = 1.4,
        core_spacing = 0.7,
        xmax_parallel = 6.0,
        xmax_transverse = 4.0,
        bond_axis = :z,
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(basis, expansion)
    term_coefficients = Float64.(expansion.coefficients)

    default_source = GaussletBases._nested_bond_aligned_diatomic_source(
        basis,
        bundles;
        bond_axis = :z,
        nside = 5,
        term_coefficients,
        packet_kernel = :factorized_direct,
    )
    endcap_source = GaussletBases._nested_bond_aligned_diatomic_source(
        basis,
        bundles;
        bond_axis = :z,
        nside = 5,
        term_coefficients,
        packet_kernel = :factorized_direct,
        shared_shell_layer_policy = :endcap_panel_owned,
        shared_shell_endcap_panel_q = 4,
        shared_shell_endcap_panel_L = 4,
    )

    @test length(default_source.shared_shell_layers) == 1
    @test all(layer isa GaussletBases._CartesianNestedCompleteShell3D for layer in default_source.shared_shell_layers)
    @test length(endcap_source.shared_shell_layers) == length(default_source.shared_shell_layers)
    @test only(endcap_source.shared_shell_layers) isa GaussletBases._CartesianNestedEndcapPanelShellLayer3D
    @test length(endcap_source.child_sequences) == length(default_source.child_sequences) == 1
    @test size(only(endcap_source.child_sequences).coefficient_matrix, 2) ==
        size(only(default_source.child_sequences).coefficient_matrix, 2)

    default_shared_columns = [size(layer.coefficient_matrix, 2) for layer in default_source.shared_shell_layers]
    endcap_shared_columns = [size(layer.coefficient_matrix, 2) for layer in endcap_source.shared_shell_layers]
    @test default_shared_columns == [130]
    @test endcap_shared_columns == [96]
    @test size(endcap_source.sequence.coefficient_matrix, 2) ==
        size(default_source.sequence.coefficient_matrix, 2) - sum(default_shared_columns) +
        sum(endcap_shared_columns)

    layer = only(endcap_source.shared_shell_layers)
    @test layer.provenance.support_contract == :thin_endcap_box_perimeter
    @test layer.provenance.coefficient_contract == :product_doside
    @test layer.provenance.q == 4
    @test layer.provenance.L == 4
    @test layer.provenance.packet_kernel == :factorized_direct
    @test layer.owned_units.audit.coverage_ok
    @test layer.owned_units.audit.expected_support_count == 314
    @test layer.owned_units.audit.owned_support_count == 314
    @test layer.owned_units.audit.duplicate_count == 0
    @test layer.owned_units.audit.missing_count == 0
    @test layer.owned_units.audit.outside_count == 0
    @test layer.owned_units.audit.retained_count == 96
    @test length(layer.support_indices) == 314
    @test size(layer.coefficient_matrix) == (prod(GaussletBases._nested_axis_lengths(bundles)), 96)
    @test all(isfinite, layer.packet.overlap)
    @test norm(layer.packet.overlap - I, Inf) < 1.0e-8
    CCP = GaussletBases.CartesianContractedParents
    CP = GaussletBases.CartesianParentGaussletBases
    CCPM = GaussletBases.CartesianContractedParentMetrics
    pgdg_x = GaussletBases._nested_axis_pgdg(bundles, :x)
    pgdg_y = GaussletBases._nested_axis_pgdg(bundles, :y)
    pgdg_z = GaussletBases._nested_axis_pgdg(bundles, :z)
    axis_metrics = (
        x = (
            overlap = pgdg_x.overlap,
            position = pgdg_x.position,
            x2 = pgdg_x.x2,
            weights = pgdg_x.weights,
            centers = pgdg_x.centers,
            source = :nested_pgdg_axis,
        ),
        y = (
            overlap = pgdg_y.overlap,
            position = pgdg_y.position,
            x2 = pgdg_y.x2,
            weights = pgdg_y.weights,
            centers = pgdg_y.centers,
            source = :nested_pgdg_axis,
        ),
        z = (
            overlap = pgdg_z.overlap,
            position = pgdg_z.position,
            x2 = pgdg_z.x2,
            weights = pgdg_z.weights,
            centers = pgdg_z.centers,
            source = :nested_pgdg_axis,
        ),
    )
    safe_axis_data = (
        x = merge(axis_metrics.x, (kinetic = pgdg_x.kinetic,)),
        y = merge(axis_metrics.y, (kinetic = pgdg_y.kinetic,)),
        z = merge(axis_metrics.z, (kinetic = pgdg_z.kinetic,)),
    )
    dims = GaussletBases._nested_axis_lengths(bundles)
    pre_packet_source = CCP._cartesian_endcap_panel_pre_packet_build_source(
        layer.owned_units,
        layer.coefficient_matrix,
        layer.unit_column_ranges,
        layer.support_indices,
        dims,
    )
    post_layer_units = [
        GaussletBases._nested_product_staged_unit_from_owned_unit(
            owned_unit;
            column_range,
            dims,
        ) for (owned_unit, column_range) in
            zip(layer.owned_units.units, layer.unit_column_ranges)
    ]
    post_layer_sidecar = GaussletBases._CartesianNestedProductStagedByCenterSidecar3D(
        dims,
        post_layer_units,
        (; source = :post_layer_test_sidecar, support_contract = :product_owned_units),
        (
            parent_dimension = prod(dims),
            final_dimension = size(layer.coefficient_matrix, 2),
            unit_count = length(post_layer_units),
            product_unit_count = length(post_layer_units),
            generic_unit_count = 0,
            support_counts = Int[unit.diagnostics.support_count for unit in post_layer_units],
            max_support_count = maximum(unit.diagnostics.support_count for unit in post_layer_units),
        ),
    )
    post_layer_source = CCP._cartesian_packet_build_source(post_layer_sidecar)
    @test pre_packet_source.parent_dimension == post_layer_source.parent_dimension == prod(dims)
    @test pre_packet_source.contracted_dimension == post_layer_source.contracted_dimension == 96
    @test pre_packet_source.payload_kind_counts == post_layer_source.payload_kind_counts
    @test Dict(pre_packet_source.payload_kind_counts)[:product_doside] == 6
    @test [payload.payload_kind for payload in pre_packet_source.resolved_payloads] ==
          [payload.payload_kind for payload in post_layer_source.resolved_payloads]
    @test [payload.column_range for payload in pre_packet_source.resolved_payloads] ==
          [payload.column_range for payload in post_layer_source.resolved_payloads]
    @test [payload.support_indices for payload in pre_packet_source.resolved_payloads] ==
          [payload.support_indices for payload in post_layer_source.resolved_payloads]
    @test pre_packet_source.candidate_packet_fields == post_layer_source.candidate_packet_fields
    @test pre_packet_source.missing_packet_fields == post_layer_source.missing_packet_fields
    @test pre_packet_source.diagnostics.packet_construction_consumes_source == false
    @test pre_packet_source.diagnostics.source_object_builds_packet_matrices == false
    @test pre_packet_source.diagnostics.nested_shell_packet_remains_authoritative
    @test !pre_packet_source.diagnostics.numerical_packet_matrices_built
    @test !pre_packet_source.diagnostics.operator_data_available
    @test !pre_packet_source.diagnostics.packet_operator_data_checked
    pre_packet_shadow = CCPM._cartesian_packet_build_source_safe_field_shadow(
        pre_packet_source,
        safe_axis_data,
    )
    layer_parent = CP.cartesian_parent_gausslet_basis(basis)
    layer_contracted_parent = CCP.cartesian_contracted_parent(
        layer_parent,
        layer.coefficient_matrix;
        units = CCP._contracted_parent_units_from_staged_sidecar(post_layer_sidecar),
        metadata = (; source = :endcap_panel_layer_pre_packet_shadow_test),
    )
    layer_metric_packet = CCPM.cartesian_contracted_parent_metric_packet(
        layer_contracted_parent;
        axis_metrics,
        construction_path = :product_staged_metric_contraction,
    )
    @test pre_packet_shadow.diagnostics.source_driven_shadow_only
    @test !pre_packet_shadow.diagnostics.construction_adoption
    @test pre_packet_shadow.diagnostics.current_builder_authority ==
          :nested_shell_packet
    @test pre_packet_shadow.diagnostics.product_unit_count == 6
    @test pre_packet_shadow.diagnostics.support_dense_unit_count == 0
    @test pre_packet_shadow.overlap ≈ layer.packet.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test pre_packet_shadow.position_x ≈ layer.packet.position_x atol = 1.0e-10 rtol = 1.0e-10
    @test pre_packet_shadow.position_y ≈ layer.packet.position_y atol = 1.0e-10 rtol = 1.0e-10
    @test pre_packet_shadow.position_z ≈ layer.packet.position_z atol = 1.0e-10 rtol = 1.0e-10
    @test pre_packet_shadow.x2_x ≈ layer.packet.x2_x atol = 1.0e-10 rtol = 1.0e-10
    @test pre_packet_shadow.x2_y ≈ layer.packet.x2_y atol = 1.0e-10 rtol = 1.0e-10
    @test pre_packet_shadow.x2_z ≈ layer.packet.x2_z atol = 1.0e-10 rtol = 1.0e-10
    @test pre_packet_shadow.weights ≈ layer.packet.weights atol = 1.0e-10 rtol = 1.0e-10
    @test pre_packet_shadow.kinetic ≈ layer.packet.kinetic atol = 1.0e-10 rtol = 1.0e-10
    @test pre_packet_shadow.first_moments ≈ layer_metric_packet.first_moments atol = 1.0e-10 rtol = 1.0e-10
    sequence_core_coefficients =
        endcap_source.sequence.coefficient_matrix[:, endcap_source.sequence.core_column_range]
    pre_sequence_source = CCP._cartesian_nested_sequence_pre_packet_build_source(
        dims,
        endcap_source.sequence.core_indices,
        sequence_core_coefficients,
        endcap_source.sequence.core_column_range,
        endcap_source.sequence.shell_layers,
        endcap_source.sequence.layer_column_ranges,
        endcap_source.sequence.coefficient_matrix,
        endcap_source.sequence.support_indices,
    )
    @test pre_sequence_source.parent_dimension == prod(dims)
    @test pre_sequence_source.contracted_dimension == size(endcap_source.sequence.coefficient_matrix, 2)
    @test Dict(pre_sequence_source.payload_kind_counts)[:product_doside] == 6
    @test Dict(pre_sequence_source.payload_kind_counts)[:support_dense] == 1
    @test pre_sequence_source.diagnostics.packet_construction_consumes_source == false
    @test pre_sequence_source.diagnostics.source_object_builds_packet_matrices == false
    @test pre_sequence_source.diagnostics.nested_shell_packet_remains_authoritative
    @test !pre_sequence_source.diagnostics.numerical_packet_matrices_built
    @test !pre_sequence_source.diagnostics.operator_data_available
    @test !pre_sequence_source.diagnostics.packet_operator_data_checked
    layer_region = CCP.cartesian_shell_region(
        layer;
        parent_dimension = prod(GaussletBases._nested_axis_lengths(bundles)),
    )
    @test layer_region isa CCP.CartesianShellRegion3D
    @test layer_region.region_family == :endcap_panel_shared_exterior
    @test layer_region.role == :shared_endcap_panel_shell_layer
    @test layer_region.status == :transitional
    @test layer_region.box == layer.provenance.current_box
    @test layer_region.inner_exclusion_box == layer.provenance.inner_box
    @test layer_region.support_summary.entry_count == length(layer.support_indices)
    @test layer_region.support_summary.outside_count == 0
    @test layer_region.ownership_coverage_contract == :boundary_only
    @test layer_region.retention.retention_rule == :old_endcap_panel_product_split
    @test layer_region.retention.cleanup_rule == :locally_orthonormal_product_doside
    @test layer_region.retention.preferred_contraction_rule ==
          :old_endcap_panel_product_split
    @test layer_region.retention.expected_unit_family == :product_owned_unit
    @test layer_region.retention.metric_capability == :product_staged_metric_contraction
    @test isempty(layer_region.retention.missing_payload_fields)
    @test layer_region.current_route_consumes
    @test !layer_region.descriptor_drives_builder
    @test !layer_region.descriptor_only
    @test layer_region.geometry.q == 4
    @test layer_region.geometry.L == 4
    @test layer_region.geometry.unit_count == 6
    @test layer_region.diagnostics.transitional_current_active_implementation

    @test size(default_source.sequence.coefficient_matrix) == (539, 347)
    @test size(endcap_source.sequence.coefficient_matrix) == (539, 313)
    @test all(isfinite, endcap_source.sequence.packet.overlap)
    @test norm(endcap_source.sequence.packet.overlap - I, Inf) < 1.0e-8

    fixed_block = GaussletBases._nested_fixed_block(endcap_source)
    @test fixed_block isa GaussletBases._NestedFixedBlock3D
    @test size(fixed_block.coefficient_matrix) == (539, 313)
    @test fixed_block.staged_by_center_sidecar[] isa
          GaussletBases._CartesianNestedProductStagedByCenterSidecar3D
    @test fixed_block.staged_by_center_sidecar[].diagnostics.product_unit_count == 6
    @test fixed_block.staged_by_center_sidecar[].diagnostics.final_dimension == 313
    @test fixed_block.staged_by_center_sidecar[].diagnostics.max_support_count <= 225
    @test all(isfinite, fixed_block.overlap)
    @test norm(fixed_block.overlap - I, Inf) < 1.0e-8
    @test all(isfinite, fixed_block.weights)

    CCS = GaussletBases.CartesianCarriedSpaces
    contracted_parent = CCP.cartesian_contracted_parent(fixed_block)
    sidecar_units = fixed_block.staged_by_center_sidecar[].units
    parent_dim = CP.parent_dimension(CCP.contracted_parent_basis(contracted_parent))
    rule_built_units = [
        CCP.cartesian_contraction_unit_from_rule(
            CCP.cartesian_contraction_rule(sidecar_unit; parent_dimension = parent_dim),
            sidecar_unit,
        ) for sidecar_unit in sidecar_units
    ]
    @test CCP.contracted_parent_metadata(contracted_parent).staged_by_center_path ==
        :product_staged_factorized
    @test length(CCP.contracted_parent_units(contracted_parent)) == length(sidecar_units)
    @test length(rule_built_units) == length(sidecar_units)
    @test first(CCP.contracted_parent_units(contracted_parent)).metadata.staged_by_center_unit ===
        first(sidecar_units)
    product_sidecar_unit = first(unit for unit in sidecar_units if unit.kind == :product_doside)
    support_sidecar_unit = first(unit for unit in sidecar_units if unit.kind == :support_dense)
    product_resolved_payload = CCP._cartesian_resolved_contraction_payload(
        product_sidecar_unit;
        parent_dimension = parent_dim,
    )
    support_resolved_payload = CCP._cartesian_resolved_contraction_payload(
        support_sidecar_unit;
        parent_dimension = parent_dim,
    )
    @test product_resolved_payload isa CCP._CartesianResolvedContractionPayload3D
    @test support_resolved_payload isa CCP._CartesianResolvedContractionPayload3D
    @test product_resolved_payload.metric_path == :product_staged_metric_contraction
    @test support_resolved_payload.metric_path == :support_local_product
    @test product_resolved_payload.ready_for_metric_execution
    @test support_resolved_payload.ready_for_metric_execution
    @test product_resolved_payload.payload_kind == :product_doside
    @test support_resolved_payload.payload_kind == :support_dense
    @test product_resolved_payload.column_range == product_sidecar_unit.column_range
    @test support_resolved_payload.column_range == support_sidecar_unit.column_range
    @test product_resolved_payload.support_indices == product_sidecar_unit.support_indices
    @test support_resolved_payload.support_indices == support_sidecar_unit.support_indices
    @test product_resolved_payload.support_states == product_sidecar_unit.support_states
    @test support_resolved_payload.support_states == support_sidecar_unit.support_states
    @test product_resolved_payload.payload === product_sidecar_unit
    @test support_resolved_payload.payload === support_sidecar_unit
    @test isempty(product_resolved_payload.missing_fields)
    @test isempty(support_resolved_payload.missing_fields)
    product_sidecar_entries = CCPM._staged_unit_entries(product_sidecar_unit)
    support_sidecar_entries = CCPM._staged_unit_entries(support_sidecar_unit)
    @test length(product_sidecar_entries) == length(product_sidecar_unit.column_range)
    @test length(support_sidecar_entries) == length(support_sidecar_unit.column_range)
    @test sum(length, product_sidecar_entries) ==
          count(!iszero, Matrix{Float64}(product_sidecar_unit.coefficient_matrix))
    @test sum(length, support_sidecar_entries) ==
          count(!iszero, Matrix{Float64}(support_sidecar_unit.coefficient_matrix))
    @test product_resolved_payload.diagnostics.linear_vector_path ==
          :product_staged_axis_projection
    @test support_resolved_payload.diagnostics.linear_vector_path ==
          :support_local_fallback
    @test product_resolved_payload.diagnostics.block_role == :product
    @test support_resolved_payload.diagnostics.block_role == :fallback
    for (sidecar_unit, rule_unit, adapted_unit) in zip(
        sidecar_units,
        rule_built_units,
        CCP.contracted_parent_units(contracted_parent),
    )
        @test rule_unit.role == sidecar_unit.role
        @test rule_unit.support_indices == sidecar_unit.support_indices
        @test rule_unit.column_range == sidecar_unit.column_range
        @test adapted_unit.role == rule_unit.role
        @test adapted_unit.support_indices == rule_unit.support_indices
        @test adapted_unit.column_range == rule_unit.column_range
        @test adapted_unit.metadata.staged_by_center_unit === sidecar_unit
        @test adapted_unit.metadata.contraction_rule.kind == sidecar_unit.kind
        @test adapted_unit.metadata.rule_driven_unit_creation
        @test adapted_unit.metadata.rule_family ==
              adapted_unit.metadata.contraction_rule.rule_family
        @test adapted_unit.metadata.rule_kind ==
              adapted_unit.metadata.contraction_rule.kind
        @test adapted_unit.metadata.rule_metric_capability ==
              adapted_unit.metadata.contraction_rule.metric_capability
        adapted_resolved_payload = CCP._cartesian_resolved_contraction_payload(
            adapted_unit;
            parent_dimension = parent_dim,
        )
        @test adapted_resolved_payload.payload === sidecar_unit
        @test adapted_resolved_payload.column_range == sidecar_unit.column_range
        @test adapted_resolved_payload.ready_for_metric_execution
    end
    contraction_rules = [
        CCP.contraction_unit_rule(
            unit;
            parent_dimension = parent_dim,
        ) for unit in CCP.contracted_parent_units(contracted_parent)
    ]
    product_rules = filter(
        rule -> CCP.contraction_rule_family(rule) == :product_owned_unit,
        contraction_rules,
    )
    support_rules = filter(
        rule -> CCP.contraction_rule_family(rule) == :support_dense_fallback,
        contraction_rules,
    )
    @test length(product_rules) == 6
    @test length(support_rules) >= 1
    first_product_rule = first(product_rules)
    @test first_product_rule.kind == :product_doside
    @test first_product_rule.support_summary.parent_dimension == 539
    @test first_product_rule.support_summary.entry_count ==
          length(first_product_rule.support_indices)
    @test first_product_rule.support_summary.duplicate_count == 0
    @test first_product_rule.support_summary.outside_count == 0
    @test CCP.contraction_rule_retained_dimension(first_product_rule) ==
          length(first_product_rule.column_range)
    @test CCP.contraction_rule_transform_rule(first_product_rule) ==
          :two_active_axis_product_doside
    @test CCP.contraction_rule_cleanup_rule(first_product_rule) ==
          :locally_orthonormal_product_doside
    @test CCP.contraction_rule_metric_capability(first_product_rule) ==
          :product_staged_metric_contraction
    @test first_product_rule.local_geometry.axis_function_index_count ==
          CCP.contraction_rule_retained_dimension(first_product_rule)
    @test first_product_rule.diagnostics.coefficient_contract == :product_doside
    first_support_rule = first(support_rules)
    @test first_support_rule.kind == :support_dense
    @test first_support_rule.support_summary.outside_count == 0
    @test CCP.contraction_rule_transform_rule(first_support_rule) ==
          :explicit_support_dense_coefficients
    @test CCP.contraction_rule_cleanup_rule(first_support_rule) ==
          :external_or_already_cleaned
    @test CCP.contraction_rule_metric_capability(first_support_rule) ==
          :support_local_product
    rule_inventory = CCP.contracted_parent_rule_inventory(contracted_parent)
    rule_family_counts = Dict(rule_inventory.rule_family_counts)
    @test rule_inventory.rule_count == length(CCP.contracted_parent_units(contracted_parent))
    @test rule_inventory.unit_count == length(CCP.contracted_parent_units(contracted_parent))
    @test rule_family_counts[:product_owned_unit] == 6
    @test rule_family_counts[:support_dense_fallback] == length(support_rules)
    @test rule_inventory.parent_dimension == 539
    @test rule_inventory.contracted_dimension == 313
    @test rule_inventory.total_retained_dimension == 313
    @test rule_inventory.support_summary.parent_dimension == 539
    @test rule_inventory.support_summary.outside_count == 0
    @test rule_inventory.support_summary.missing_count == 0
    @test rule_inventory.support_summary.support_complete
    @test Set(rule_inventory.metric_capabilities) ==
          Set([:product_staged_metric_contraction, :support_local_product])
    @test rule_inventory.every_unit_has_rule_metadata
    @test rule_inventory.every_unit_rule_derivable
    @test !rule_inventory.any_metadata_only_rule
    @test !rule_inventory.any_prototype_rule
    @test rule_inventory.diagnostics.parent_level_unit_inventory
    @test rule_inventory.diagnostics.all_rules_have_column_ranges
    dispatch_shadow = CCPM._contracted_parent_metric_dispatch_shadow_plan(contracted_parent)
    resolved_payloads = [
        CCP._cartesian_resolved_contraction_payload(unit; parent_dimension = parent_dim)
        for unit in CCP.contracted_parent_units(contracted_parent)
    ]
    @test dispatch_shadow.comparison.agree
    @test isempty(dispatch_shadow.comparison.mismatch_fields)
    @test dispatch_shadow.payload_plan.plan_supported
    @test dispatch_shadow.rule_plan.plan_supported
    @test dispatch_shadow.payload_plan.unit_count == rule_inventory.rule_count
    @test dispatch_shadow.rule_plan.unit_count == rule_inventory.rule_count
    @test dispatch_shadow.resolved_plan.unit_count == rule_inventory.rule_count
    @test dispatch_shadow.payload_plan.product_unit_count == 6
    @test dispatch_shadow.rule_plan.product_unit_count == 6
    @test dispatch_shadow.resolved_plan.product_unit_count == 6
    @test dispatch_shadow.payload_plan.support_fallback_unit_count == length(support_rules)
    @test dispatch_shadow.rule_plan.support_fallback_unit_count == length(support_rules)
    @test dispatch_shadow.resolved_plan.support_fallback_unit_count == length(support_rules)
    expected_product_blocks = 6 * (6 + 1) ÷ 2
    expected_total_blocks = rule_inventory.rule_count * (rule_inventory.rule_count + 1) ÷ 2
    @test dispatch_shadow.payload_plan.product_product_block_count == expected_product_blocks
    @test dispatch_shadow.rule_plan.product_product_block_count == expected_product_blocks
    @test dispatch_shadow.payload_plan.fallback_block_count ==
          expected_total_blocks - expected_product_blocks
    @test dispatch_shadow.rule_plan.fallback_block_count ==
          expected_total_blocks - expected_product_blocks
    @test dispatch_shadow.resolved_plan.fallback_block_count ==
          expected_total_blocks - expected_product_blocks
    @test dispatch_shadow.payload_plan.unsupported_unit_count == 0
    @test dispatch_shadow.rule_plan.unsupported_unit_count == 0
    @test dispatch_shadow.resolved_plan.unsupported_unit_count == 0
    @test dispatch_shadow.rule_plan.prototype_rule_count == 0
    @test dispatch_shadow.resolved_plan.prototype_rule_count == 0
    @test dispatch_shadow.resolved_comparison.agree
    @test isempty(dispatch_shadow.resolved_comparison.mismatch_fields)
    @test all(payload -> payload.ready_for_metric_execution, resolved_payloads)
    @test [payload.diagnostics.block_role for payload in resolved_payloads] ==
          [path.block_role for path in dispatch_shadow.payload_plan.unit_paths]
    @test [payload.diagnostics.linear_vector_path for payload in resolved_payloads] ==
          [path.linear_vector_path for path in dispatch_shadow.payload_plan.unit_paths]
    @test [payload.diagnostics.metric_capability for payload in resolved_payloads] ==
          [path.metric_capability for path in dispatch_shadow.payload_plan.unit_paths]
    @test [path.path for path in dispatch_shadow.payload_plan.block_paths] ==
          [path.path for path in dispatch_shadow.rule_plan.block_paths]
    packet_build_plan = CCP._cartesian_packet_build_plan(
        fixed_block.staged_by_center_sidecar[],
    )
    packet_build_source = packet_build_plan.source
    packet_payload_counts = Dict(packet_build_source.payload_kind_counts)
    pre_sequence_payload_counts = Dict(pre_sequence_source.payload_kind_counts)
    @test packet_build_plan isa CCP._CartesianPacketBuildPlan3D
    @test packet_build_source isa CCP._CartesianPacketBuildSource3D
    @test packet_build_source.parent_dimension == parent_dim
    @test packet_build_source.contracted_dimension == 313
    @test pre_sequence_source isa CCP._CartesianPacketBuildSource3D
    @test pre_sequence_source.parent_dimension == packet_build_source.parent_dimension
    @test pre_sequence_source.contracted_dimension == packet_build_source.contracted_dimension
    @test pre_sequence_source.payload_kind_counts == packet_build_source.payload_kind_counts
    @test pre_sequence_payload_counts[:product_doside] == 6
    @test pre_sequence_payload_counts[:support_dense] == length(support_rules)
    @test length(packet_build_source.resolved_payloads) == length(resolved_payloads)
    @test length(pre_sequence_source.resolved_payloads) ==
          length(packet_build_source.resolved_payloads)
    @test [payload.payload_kind for payload in pre_sequence_source.resolved_payloads] ==
          [payload.payload_kind for payload in packet_build_source.resolved_payloads]
    @test [payload.column_range for payload in pre_sequence_source.resolved_payloads] ==
          [payload.column_range for payload in packet_build_source.resolved_payloads]
    @test [payload.support_indices for payload in pre_sequence_source.resolved_payloads] ==
          [payload.support_indices for payload in packet_build_source.resolved_payloads]
    @test [payload.payload for payload in packet_build_source.resolved_payloads] ==
          [payload.payload for payload in resolved_payloads]
    @test [payload.payload_kind for payload in packet_build_source.resolved_payloads] ==
          [payload.payload_kind for payload in resolved_payloads]
    @test [payload.column_range for payload in packet_build_source.resolved_payloads] ==
          [payload.column_range for payload in resolved_payloads]
    @test all(payload -> payload.ready_for_metric_execution, packet_build_source.resolved_payloads)
    @test packet_build_source.column_ranges ==
          [unit.column_range for unit in sidecar_units]
    @test packet_build_source.column_coverage.entry_count == 313
    @test packet_build_source.column_coverage.unique_count == 313
    @test packet_build_source.column_coverage.duplicate_count == 0
    @test packet_build_source.column_coverage.missing_count == 0
    @test packet_build_source.column_coverage.outside_count == 0
    @test packet_build_source.support_union_summary.parent_dimension == parent_dim
    @test packet_build_source.support_union_summary.outside_count == 0
    @test pre_sequence_source.column_coverage.entry_count ==
          packet_build_source.column_coverage.entry_count
    @test pre_sequence_source.column_coverage.unique_count ==
          packet_build_source.column_coverage.unique_count
    @test pre_sequence_source.column_coverage.duplicate_count ==
          packet_build_source.column_coverage.duplicate_count
    @test pre_sequence_source.column_coverage.missing_count ==
          packet_build_source.column_coverage.missing_count
    @test pre_sequence_source.column_coverage.outside_count ==
          packet_build_source.column_coverage.outside_count
    @test pre_sequence_source.support_union_summary.entry_count ==
          packet_build_source.support_union_summary.entry_count
    @test pre_sequence_source.support_union_summary.unique_count ==
          packet_build_source.support_union_summary.unique_count
    @test pre_sequence_source.support_union_summary.outside_count ==
          packet_build_source.support_union_summary.outside_count
    @test packet_payload_counts[:product_doside] == 6
    @test packet_payload_counts[:support_dense] == length(support_rules)
    @test packet_build_source.candidate_packet_fields == (
        :overlap,
        :position_x,
        :position_y,
        :position_z,
        :x2_x,
        :x2_y,
        :x2_z,
        :weights,
        :first_moments,
        :kinetic,
    )
    @test !(:x2_x in packet_build_source.missing_packet_fields)
    @test !(:x2_y in packet_build_source.missing_packet_fields)
    @test !(:x2_z in packet_build_source.missing_packet_fields)
    @test :nuclear_one_body in packet_build_source.missing_packet_fields
    @test :local_coulomb_one_body in packet_build_source.missing_packet_fields
    @test :local_ecp_one_body in packet_build_source.missing_packet_fields
    @test :gaussian_local_terms in packet_build_source.missing_packet_fields
    @test :mwg_interaction in packet_build_source.missing_packet_fields
    @test packet_build_source.axis_operator_requirements.kinetic ==
          (:overlap, :kinetic)
    @test packet_build_source.diagnostics.metadata_only
    @test packet_build_source.diagnostics.descriptor_does_not_drive_builder
    @test packet_build_source.diagnostics.current_builder_authority ==
          :nested_shell_packet
    @test !packet_build_source.diagnostics.numerical_packet_matrices_built
    @test !packet_build_source.diagnostics.operator_data_available
    @test !packet_build_source.diagnostics.packet_operator_data_checked
    @test !packet_build_source.diagnostics.overlap_matrix_built
    @test !packet_build_source.diagnostics.kinetic_matrix_built
    @test packet_build_source.diagnostics.column_ranges_cover_contract
    @test packet_build_source.diagnostics.column_layout_ready_for_candidate_fields
    @test packet_build_source.diagnostics.support_union_summary_informational
    @test !packet_build_source.diagnostics.parent_support_complete_required
    @test packet_build_source.diagnostics.overlapping_payload_support_allowed
    @test pre_sequence_source.candidate_packet_fields ==
          packet_build_source.candidate_packet_fields
    @test pre_sequence_source.missing_packet_fields ==
          packet_build_source.missing_packet_fields
    @test pre_sequence_source.diagnostics.packet_construction_consumes_source == false
    @test pre_sequence_source.diagnostics.source_object_builds_packet_matrices == false
    @test pre_sequence_source.diagnostics.nested_shell_packet_remains_authoritative
    @test !pre_sequence_source.diagnostics.numerical_packet_matrices_built
    @test !pre_sequence_source.diagnostics.operator_data_available
    @test !pre_sequence_source.diagnostics.packet_operator_data_checked
    @test !packet_build_source.diagnostics.full_packet_builder_ready
    @test packet_build_plan.current_builder_authority == :nested_shell_packet
    @test !packet_build_plan.descriptor_drives_builder
    @test !packet_build_plan.numerical_packet_matrices_built
    @test packet_build_plan.diagnostics.metadata_only
    @test packet_build_plan.diagnostics.current_nested_shell_packet_authoritative
    @test !packet_build_plan.diagnostics.fixed_block_construction_changed
    @test !packet_build_plan.diagnostics.metric_packet_execution_changed
    @test !packet_build_plan.diagnostics.qwhamiltonian_changed
    @test !packet_build_plan.diagnostics.backend_policy_changed
    @test !packet_build_plan.diagnostics.quadrature_policy_changed
    @test !packet_build_plan.diagnostics.ida_positive_weight_semantics_changed
    @test !packet_build_plan.diagnostics.cr2_science_status_changed
    @test !packet_build_plan.diagnostics.operator_data_available
    @test !packet_build_plan.diagnostics.packet_operator_data_checked
    support_metric_packet = CCPM.cartesian_contracted_parent_metric_packet(
        contracted_parent;
        axis_metrics,
    )
    product_metric_packet = CCPM.cartesian_contracted_parent_metric_packet(
        contracted_parent;
        axis_metrics,
        construction_path = :product_staged_metric_contraction,
    )
    resolved_metric_packet = CCPM._resolved_payload_product_staged_metric_packet(
        contracted_parent;
        axis_metrics,
    )
    @test CP.parent_dimension(CCP.contracted_parent_basis(contracted_parent)) == 539
    @test support_metric_packet.diagnostics.construction_path == :support_local_product
    @test product_metric_packet.diagnostics.construction_path == :product_staged_metric_contraction
    @test !(:resolved_payload_count in propertynames(product_metric_packet.diagnostics))
    @test !(:default_metric_execution_changed in propertynames(product_metric_packet.diagnostics))
    @test resolved_metric_packet.diagnostics.construction_path ==
          :resolved_payload_product_staged_metric_contraction
    @test resolved_metric_packet.diagnostics.source == :resolved_payload_metric_shadow
    @test !resolved_metric_packet.diagnostics.default_metric_execution_changed
    @test resolved_metric_packet.diagnostics.resolved_payload_count ==
          rule_inventory.rule_count
    @test product_metric_packet.diagnostics.dense_parent_matrix_used == false
    @test resolved_metric_packet.diagnostics.dense_parent_matrix_used == false
    @test product_metric_packet.diagnostics.product_unit_count == 6
    @test resolved_metric_packet.diagnostics.product_unit_count ==
          product_metric_packet.diagnostics.product_unit_count
    @test product_metric_packet.diagnostics.generic_unit_count >= 1
    @test resolved_metric_packet.diagnostics.generic_unit_count ==
          product_metric_packet.diagnostics.generic_unit_count
    @test product_metric_packet.diagnostics.product_block_count > 0
    @test resolved_metric_packet.diagnostics.product_block_count ==
          product_metric_packet.diagnostics.product_block_count
    @test product_metric_packet.diagnostics.fallback_block_count > 0
    @test resolved_metric_packet.diagnostics.fallback_block_count ==
          product_metric_packet.diagnostics.fallback_block_count
    @test resolved_metric_packet.overlap ≈ product_metric_packet.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test resolved_metric_packet.weights ≈ product_metric_packet.weights atol = 1.0e-12 rtol = 1.0e-12
    @test resolved_metric_packet.centers ≈ product_metric_packet.centers atol = 1.0e-12 rtol = 1.0e-12
    @test resolved_metric_packet.first_moments ≈ product_metric_packet.first_moments atol = 1.0e-12 rtol = 1.0e-12
    @test product_metric_packet.overlap ≈ support_metric_packet.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test product_metric_packet.weights ≈ support_metric_packet.weights atol = 1.0e-10 rtol = 1.0e-10
    @test product_metric_packet.centers ≈ support_metric_packet.centers atol = 1.0e-10 rtol = 1.0e-10
    @test product_metric_packet.first_moments ≈ support_metric_packet.first_moments atol = 1.0e-10 rtol = 1.0e-10
    @test product_metric_packet.overlap ≈ fixed_block.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test product_metric_packet.weights ≈ fixed_block.weights atol = 1.0e-10 rtol = 1.0e-10
    @test product_metric_packet.centers ≈ fixed_block.fixed_centers atol = 1.0e-10 rtol = 1.0e-10
    pre_sequence_shadow = CCPM._cartesian_packet_build_source_safe_field_shadow(
        pre_sequence_source,
        safe_axis_data,
    )
    @test pre_sequence_shadow.diagnostics.source_driven_shadow_only
    @test !pre_sequence_shadow.diagnostics.construction_adoption
    @test pre_sequence_shadow.diagnostics.current_builder_authority ==
          :nested_shell_packet
    @test pre_sequence_shadow.diagnostics.product_unit_count == 6
    @test pre_sequence_shadow.diagnostics.support_dense_unit_count == length(support_rules)
    @test pre_sequence_shadow.overlap ≈ fixed_block.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test pre_sequence_shadow.position_x ≈ fixed_block.position_x atol = 1.0e-10 rtol = 1.0e-10
    @test pre_sequence_shadow.position_y ≈ fixed_block.position_y atol = 1.0e-10 rtol = 1.0e-10
    @test pre_sequence_shadow.position_z ≈ fixed_block.position_z atol = 1.0e-10 rtol = 1.0e-10
    @test pre_sequence_shadow.x2_x ≈ fixed_block.x2_x atol = 1.0e-10 rtol = 1.0e-10
    @test pre_sequence_shadow.x2_y ≈ fixed_block.x2_y atol = 1.0e-10 rtol = 1.0e-10
    @test pre_sequence_shadow.x2_z ≈ fixed_block.x2_z atol = 1.0e-10 rtol = 1.0e-10
    @test pre_sequence_shadow.weights ≈ fixed_block.weights atol = 1.0e-10 rtol = 1.0e-10
    @test pre_sequence_shadow.kinetic ≈ fixed_block.kinetic atol = 1.0e-10 rtol = 1.0e-10
    @test pre_sequence_shadow.first_moments ≈ product_metric_packet.first_moments atol = 1.0e-10 rtol = 1.0e-10
    safe_field_shadow = CCPM._cartesian_packet_build_source_safe_field_shadow(
        packet_build_source,
        safe_axis_data,
    )
    @test safe_field_shadow.diagnostics.source ==
          :cartesian_packet_build_source_safe_field_shadow
    @test safe_field_shadow.diagnostics.source_driven_shadow_only
    @test !safe_field_shadow.diagnostics.construction_adoption
    @test safe_field_shadow.diagnostics.current_builder_authority ==
          :nested_shell_packet
    @test safe_field_shadow.diagnostics.nested_shell_packet_remains_authoritative
    @test !safe_field_shadow.diagnostics.descriptor_drives_builder
    @test !safe_field_shadow.diagnostics.numerical_packet_authority_changed
    @test !safe_field_shadow.diagnostics.fixed_block_construction_changed
    @test !safe_field_shadow.diagnostics.metric_packet_execution_changed
    @test !safe_field_shadow.diagnostics.qwhamiltonian_changed
    @test !safe_field_shadow.diagnostics.backend_policy_changed
    @test !safe_field_shadow.diagnostics.quadrature_policy_changed
    @test !safe_field_shadow.diagnostics.ida_positive_weight_semantics_changed
    @test !safe_field_shadow.diagnostics.cr2_science_status_changed
    @test safe_field_shadow.diagnostics.requested_fields ==
          packet_build_source.candidate_packet_fields
    @test safe_field_shadow.diagnostics.missing_packet_fields ==
          packet_build_source.missing_packet_fields
    @test safe_field_shadow.diagnostics.x2_built
    @test !safe_field_shadow.diagnostics.gaussian_terms_built
    @test !safe_field_shadow.diagnostics.nuclear_one_body_built
    @test !safe_field_shadow.diagnostics.local_coulomb_one_body_built
    @test !safe_field_shadow.diagnostics.local_ecp_one_body_built
    @test !safe_field_shadow.diagnostics.pair_mwg_interaction_built
    @test safe_field_shadow.diagnostics.product_unit_count == 6
    @test safe_field_shadow.diagnostics.support_dense_unit_count == length(support_rules)
    @test safe_field_shadow.diagnostics.low_order_product_block_count == expected_product_blocks
    @test safe_field_shadow.diagnostics.kinetic_product_block_count == expected_product_blocks
    @test safe_field_shadow.diagnostics.low_order_fallback_block_count ==
          expected_total_blocks - expected_product_blocks
    @test safe_field_shadow.diagnostics.x2_product_block_count == expected_product_blocks
    @test safe_field_shadow.diagnostics.x2_fallback_block_count ==
          expected_total_blocks - expected_product_blocks
    @test safe_field_shadow.diagnostics.kinetic_fallback_block_count ==
          expected_total_blocks - expected_product_blocks
    @test safe_field_shadow.overlap ≈ fixed_block.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test safe_field_shadow.position_x ≈ fixed_block.position_x atol = 1.0e-10 rtol = 1.0e-10
    @test safe_field_shadow.position_y ≈ fixed_block.position_y atol = 1.0e-10 rtol = 1.0e-10
    @test safe_field_shadow.position_z ≈ fixed_block.position_z atol = 1.0e-10 rtol = 1.0e-10
    @test safe_field_shadow.x2_x ≈ fixed_block.x2_x atol = 1.0e-10 rtol = 1.0e-10
    @test safe_field_shadow.x2_y ≈ fixed_block.x2_y atol = 1.0e-10 rtol = 1.0e-10
    @test safe_field_shadow.x2_z ≈ fixed_block.x2_z atol = 1.0e-10 rtol = 1.0e-10
    @test safe_field_shadow.weights ≈ fixed_block.weights atol = 1.0e-10 rtol = 1.0e-10
    @test safe_field_shadow.kinetic ≈ fixed_block.kinetic atol = 1.0e-10 rtol = 1.0e-10
    @test safe_field_shadow.first_moments ≈ product_metric_packet.first_moments atol = 1.0e-10 rtol = 1.0e-10
    @test_throws ArgumentError CCPM._cartesian_packet_build_source_safe_field_shadow(
        packet_build_source,
        safe_axis_data;
        fields = (:gaussian_sum,),
    )
    kinetic_axis_ops = (
        x = (overlap = pgdg_x.overlap, kinetic = pgdg_x.kinetic),
        y = (overlap = pgdg_y.overlap, kinetic = pgdg_y.kinetic),
        z = (overlap = pgdg_z.overlap, kinetic = pgdg_z.kinetic),
    )
    kinetic_product_units = [
        unit for unit in sidecar_units if unit.kind == :product_doside
    ]
    @test length(kinetic_product_units) ==
          fixed_block.staged_by_center_sidecar[].diagnostics.product_unit_count
    kinetic_shadow_packet = CCPM._product_doside_retained_kinetic_shadow_matrix(
        sidecar_units,
        kinetic_axis_ops;
        final_dimension = size(fixed_block.kinetic, 1),
    )
    @test kinetic_shadow_packet.product_units == kinetic_product_units
    @test size(kinetic_shadow_packet.kinetic) == size(fixed_block.kinetic)
    @test kinetic_shadow_packet.diagnostics.product_only_shadow
    @test !kinetic_shadow_packet.diagnostics.full_kinetic_packet
    @test kinetic_shadow_packet.diagnostics.support_dense_blocks_absent
    @test kinetic_shadow_packet.diagnostics.mixed_blocks_absent
    @test !kinetic_shadow_packet.diagnostics.default_execution_changed
    @test !kinetic_shadow_packet.diagnostics.qwhamiltonian_consumes
    @test !kinetic_shadow_packet.diagnostics.backend_policy_changed
    @test !kinetic_shadow_packet.diagnostics.quadrature_policy_changed
    @test !kinetic_shadow_packet.diagnostics.cr2_science_status_changed
    @test kinetic_shadow_packet.diagnostics.product_unit_count ==
          length(kinetic_product_units)
    @test kinetic_shadow_packet.diagnostics.final_dimension ==
          size(fixed_block.kinetic, 1)
    kinetic_shadow_block_count = 0
    kinetic_shadow_errors = Float64[]
    for right_index in eachindex(kinetic_product_units)
        right = kinetic_product_units[right_index]
        for left_index in 1:right_index
            left = kinetic_product_units[left_index]
            kinetic_shadow_block =
                kinetic_shadow_packet.kinetic[left.column_range, right.column_range]
            kinetic_reference_block =
                fixed_block.kinetic[left.column_range, right.column_range]
            push!(
                kinetic_shadow_errors,
                maximum(abs.(kinetic_shadow_block .- kinetic_reference_block)),
            )
            @test kinetic_shadow_block ≈ kinetic_reference_block atol = 1.0e-10 rtol = 1.0e-10
            if left_index != right_index
                kinetic_mirror_reference =
                    fixed_block.kinetic[right.column_range, left.column_range]
                @test transpose(kinetic_shadow_block) ≈ kinetic_mirror_reference atol = 1.0e-10 rtol = 1.0e-10
            end
            kinetic_shadow_block_count += 1
        end
    end
    @test kinetic_shadow_block_count ==
          length(kinetic_product_units) * (length(kinetic_product_units) + 1) ÷ 2
    @test kinetic_shadow_packet.diagnostics.product_block_count ==
          kinetic_shadow_block_count
    @test maximum(kinetic_shadow_errors) <= 1.0e-10
    product_columns = Int[]
    for unit in kinetic_product_units
        append!(product_columns, unit.column_range)
    end
    product_columns = sort(unique(product_columns))
    nonproduct_columns = setdiff(1:size(fixed_block.kinetic, 1), product_columns)
    @test !isempty(nonproduct_columns)
    @test maximum(abs.(kinetic_shadow_packet.kinetic[nonproduct_columns, :])) == 0.0
    @test maximum(abs.(kinetic_shadow_packet.kinetic[:, nonproduct_columns])) == 0.0
    full_kinetic_shadow_packet = CCPM._staged_retained_kinetic_shadow_matrix(
        sidecar_units,
        kinetic_axis_ops;
        final_dimension = size(fixed_block.kinetic, 1),
    )
    expected_kinetic_total_blocks = length(sidecar_units) * (length(sidecar_units) + 1) ÷ 2
    @test size(full_kinetic_shadow_packet.kinetic) == size(fixed_block.kinetic)
    @test full_kinetic_shadow_packet.diagnostics.full_private_shadow_matrix
    @test !full_kinetic_shadow_packet.diagnostics.product_only_shadow
    @test !full_kinetic_shadow_packet.diagnostics.production_adoption
    @test !full_kinetic_shadow_packet.diagnostics.default_execution_changed
    @test !full_kinetic_shadow_packet.diagnostics.metric_packet_execution_changed
    @test !full_kinetic_shadow_packet.diagnostics.fixed_block_construction_changed
    @test !full_kinetic_shadow_packet.diagnostics.qwhamiltonian_consumes
    @test !full_kinetic_shadow_packet.diagnostics.backend_policy_changed
    @test !full_kinetic_shadow_packet.diagnostics.quadrature_policy_changed
    @test !full_kinetic_shadow_packet.diagnostics.ida_positive_weight_semantics_changed
    @test !full_kinetic_shadow_packet.diagnostics.cr2_science_status_changed
    @test full_kinetic_shadow_packet.diagnostics.product_unit_count ==
          length(kinetic_product_units)
    @test full_kinetic_shadow_packet.diagnostics.generic_unit_count ==
          length(sidecar_units) - length(kinetic_product_units)
    @test full_kinetic_shadow_packet.diagnostics.product_block_count ==
          kinetic_shadow_block_count
    @test full_kinetic_shadow_packet.diagnostics.fallback_block_count ==
          expected_kinetic_total_blocks - kinetic_shadow_block_count
    @test full_kinetic_shadow_packet.diagnostics.total_block_count ==
          expected_kinetic_total_blocks
    @test full_kinetic_shadow_packet.diagnostics.final_dimension ==
          size(fixed_block.kinetic, 1)
    @test maximum(abs.(full_kinetic_shadow_packet.kinetic .- fixed_block.kinetic)) <= 1.0e-10

    carried = CCS.cartesian_carried_space(fixed_block)
    carried_parent = CCS.carried_space_parent(carried)
    carried_contracted_parent = CCS.carried_space_contracted_parent(carried)
    carried_representation = CCS.carried_space_representation(carried)
    carried_diagnostics = CCS.carried_space_diagnostics(carried)
    @test carried_parent isa CP.CartesianParentGaussletBasis3D
    @test carried_contracted_parent isa CCP.CartesianContractedParent3D
    @test carried_representation isa CartesianBasisRepresentation3D
    @test CP.parent_dimension(carried_parent) == 539
    @test carried_diagnostics.parent_axis_counts == CP.parent_axis_counts(carried_parent)
    @test carried_diagnostics.has_contracted_parent
    @test carried_diagnostics.contracted_dimension == 313
    @test carried_diagnostics.contracted_dimension_matches_representation
    @test carried_diagnostics.contracted_parent_dimension_matches_parent
    @test carried_diagnostics.has_staged_sidecar
    @test carried_diagnostics.staged_by_center_path == :product_staged_factorized
    @test carried_diagnostics.dense_parent_matrix_used == false
    @test carried_diagnostics.heavy_metric_packet_built == false
    @test CCS.carried_space_provenance(carried).input_kind == :nested_fixed_block
    @test CCP.contracted_parent_metadata(carried_contracted_parent).staged_by_center_sidecar ===
        fixed_block.staged_by_center_sidecar[]
    @test first(CCP.contracted_parent_units(carried_contracted_parent)).metadata.staged_by_center_unit ===
        first(sidecar_units)
    carried_metric_packet = CCPM.cartesian_contracted_parent_metric_packet(
        carried_contracted_parent;
        axis_metrics,
        construction_path = :product_staged_metric_contraction,
    )
    @test carried_metric_packet.diagnostics.construction_path ==
        :product_staged_metric_contraction
    @test carried_metric_packet.diagnostics.dense_parent_matrix_used == false
    @test carried_metric_packet.overlap ≈ product_metric_packet.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test carried_metric_packet.weights ≈ product_metric_packet.weights atol = 1.0e-10 rtol = 1.0e-10
    @test carried_metric_packet.centers ≈ product_metric_packet.centers atol = 1.0e-10 rtol = 1.0e-10

    @test_throws ArgumentError GaussletBases._nested_bond_aligned_diatomic_source(
        basis,
        bundles;
        bond_axis = :z,
        nside = 5,
        term_coefficients,
        packet_kernel = :factorized_direct,
        shared_shell_layer_policy = :bad_policy,
    )
end

@testset "Cartesian nested shell first packet" begin
    function _fixed_a_nested_shell_basis(count::Int; a::Float64 = 0.25, xmax::Float64 = 10.0, tail_spacing::Float64 = 10.0)
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
    basis, s = _fixed_a_nested_shell_basis(13)
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    pgdg = bundle.pgdg_intermediate
    term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
    interval = 2:(length(basis) - 1)
    shell = GaussletBases._nested_xy_shell_pair(
        bundle,
        interval,
        interval;
        retain_x = 4,
        retain_y = 3,
        term_coefficients,
    )
    packet = shell.packet
    face_low, face_high = shell.faces
    nface = size(face_low.coefficient_matrix, 2)
    direct_overlap = transpose(shell.coefficient_matrix) * shell.coefficient_matrix
    low_z_mean = sum(diag(packet.position_z)[1:nface]) / nface
    high_z_mean = sum(diag(packet.position_z)[(nface + 1):end]) / nface

    @test s > 0.0
    @test shell isa GaussletBases._CartesianNestedXYShell3D
    @test shell.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test size(shell.coefficient_matrix) == (length(basis)^3, 2 * nface)
    @test nface == 9
    @test length(shell.support_indices) == 2 * length(interval)^2
    @test length(shell.support_states) == length(shell.support_indices)
    @test isempty(intersect(face_low.support_indices, face_high.support_indices))
    @test norm(packet.overlap - I, Inf) < 1.0e-10
    @test packet.overlap ≈ direct_overlap atol = 1.0e-12 rtol = 1.0e-12
    @test packet.kinetic ≈ transpose(packet.kinetic) atol = 1.0e-10 rtol = 1.0e-10
    @test packet.position_x ≈ transpose(packet.position_x) atol = 1.0e-10 rtol = 1.0e-10
    @test packet.position_y ≈ transpose(packet.position_y) atol = 1.0e-10 rtol = 1.0e-10
    @test packet.position_z ≈ transpose(packet.position_z) atol = 1.0e-10 rtol = 1.0e-10
    @test packet.x2_x ≈ transpose(packet.x2_x) atol = 1.0e-10 rtol = 1.0e-10
    @test packet.x2_y ≈ transpose(packet.x2_y) atol = 1.0e-10 rtol = 1.0e-10
    @test packet.x2_z ≈ transpose(packet.x2_z) atol = 1.0e-10 rtol = 1.0e-10
    @test !hasproperty(packet, :gaussian_terms)
    @test !hasproperty(packet, :pair_terms)
    @test !hasproperty(packet, :term_storage)
    @test !isnothing(packet.gaussian_sum)
    @test !isnothing(packet.pair_sum)
    support_coefficients = GaussletBases._nested_support_coefficient_slice(
        shell.coefficient_matrix,
        shell.support_indices,
    )
    @test support_coefficients isa SparseMatrixCSC{Float64,Int}
    support_workspace, contraction_scratch = GaussletBases._nested_support_reference_workspaces(
        support_coefficients,
        length(shell.support_indices),
        size(shell.coefficient_matrix, 2),
    )
    @test size(support_workspace) == (0, 0)
    @test size(contraction_scratch) == (0, 0)
    support_weights = GaussletBases._nested_support_weights(shell.support_states, pgdg.weights)
    fixed_weights = vec(transpose(support_coefficients) * support_weights)
    weighted_support_coefficients = support_coefficients .* reshape(1.0 ./ fixed_weights, 1, :)
    @test weighted_support_coefficients isa SparseMatrixCSC{Float64,Int}
    support_axes = GaussletBases._nested_support_axes(shell.support_states)
    gaussian_reference = GaussletBases._nested_support_reference_gaussian_sum(
        support_axes,
        support_coefficients,
        support_workspace,
        contraction_scratch,
        term_coefficients,
        pgdg.gaussian_factor_terms,
        pgdg.gaussian_factor_terms,
        pgdg.gaussian_factor_terms,
    )
    pair_reference = GaussletBases._nested_support_reference_pair_sum(
        support_axes,
        weighted_support_coefficients,
        support_workspace,
        contraction_scratch,
        term_coefficients,
        pgdg.pair_factor_terms_raw,
        pgdg.pair_factor_terms_raw,
        pgdg.pair_factor_terms_raw,
    )
    @test packet.gaussian_sum ≈ gaussian_reference atol = 1.0e-10 rtol = 1.0e-10
    @test packet.pair_sum ≈ pair_reference atol = 1.0e-10 rtol = 1.0e-10
    @test low_z_mean < 0.0
    @test high_z_mean > 0.0
end

@testset "Cartesian nested support immediate contraction helpers" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    basis = build_basis(MappedUniformBasisSpec(:G10;
        count = 13,
        mapping = AsinhMapping(a = 0.25, s = asinh(10.0 / 0.25) / (6.0 - 1.0), tail_spacing = 10.0),
        reference_spacing = 1.0,
    ))
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    pgdg = bundle.pgdg_intermediate
    term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
    interval = 2:(length(basis) - 1)
    shell = GaussletBases._nested_rectangular_shell(
        bundle,
        interval,
        interval,
        interval;
        retain_xy = (4, 3),
        retain_xz = (4, 3),
        retain_yz = (4, 3),
        term_coefficients,
    )
    support_states = shell.support_states
    support_coefficients = Matrix{Float64}(shell.coefficient_matrix[shell.support_indices, :])
    nsupport = length(support_states)
    nfixed = size(support_coefficients, 2)
    workspace = Matrix{Float64}(undef, nsupport, nsupport)
    scratch = Matrix{Float64}(undef, nfixed, nsupport)

    overlap_support = GaussletBases._nested_support_product_matrix(
        support_states,
        pgdg.overlap,
        pgdg.overlap,
        pgdg.overlap,
    )
    overlap_workspace = Matrix{Float64}(undef, nsupport, nsupport)
    GaussletBases._nested_fill_support_product_matrix!(
        overlap_workspace,
        support_states,
        pgdg.overlap,
        pgdg.overlap,
        pgdg.overlap,
    )
    overlap_reference = transpose(support_coefficients) * overlap_support * support_coefficients
    overlap_contracted = Matrix{Float64}(undef, nfixed, nfixed)
    GaussletBases._nested_contract_support_product!(
        overlap_contracted,
        workspace,
        scratch,
        support_states,
        support_coefficients,
        pgdg.overlap,
        pgdg.overlap,
        pgdg.overlap;
        beta = 0.0,
    )

    kinetic_reference_support = GaussletBases._nested_sum_of_support_products(
        support_states,
        (
            (pgdg.kinetic, pgdg.overlap, pgdg.overlap),
            (pgdg.overlap, pgdg.kinetic, pgdg.overlap),
            (pgdg.overlap, pgdg.overlap, pgdg.kinetic),
        ),
    )
    kinetic_reference = transpose(support_coefficients) * kinetic_reference_support * support_coefficients
    kinetic_contracted = Matrix{Float64}(undef, nfixed, nfixed)
    GaussletBases._nested_contract_sum_of_support_products!(
        kinetic_contracted,
        workspace,
        scratch,
        support_states,
        support_coefficients,
        (
            (pgdg.kinetic, pgdg.overlap, pgdg.overlap),
            (pgdg.overlap, pgdg.kinetic, pgdg.overlap),
            (pgdg.overlap, pgdg.overlap, pgdg.kinetic),
        );
        beta = 0.0,
    )

    @test overlap_workspace ≈ overlap_support atol = 0.0 rtol = 0.0
    @test overlap_contracted ≈ overlap_reference atol = 1.0e-10 rtol = 1.0e-10
    @test kinetic_contracted ≈ kinetic_reference atol = 1.0e-10 rtol = 1.0e-10
end

@testset "Cartesian nested shell interface" begin
    function _fixed_a_multi_face_basis(count::Int; a::Float64 = 0.25, xmax::Float64 = 10.0, tail_spacing::Float64 = 10.0)
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
    basis, s = _fixed_a_multi_face_basis(13)
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    interval = 2:(length(basis) - 1)
    term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
    shell = GaussletBases._nested_rectangular_shell(
        bundle,
        interval,
        interval,
        interval;
        retain_xy = (4, 3),
        retain_xz = (4, 3),
        retain_yz = (4, 3),
        term_coefficients,
    )
    packet = shell.packet
    direct_overlap = transpose(shell.coefficient_matrix) * shell.coefficient_matrix

    @test s > 0.0
    @test shell isa GaussletBases._CartesianNestedShell3D
    @test shell.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test length(shell.faces) == 6
    @test length(shell.face_column_ranges) == 6
    @test size(shell.coefficient_matrix) == (length(basis)^3, 54)
    @test length(shell.support_indices) == 6 * length(interval)^2
    @test length(shell.support_states) == length(shell.support_indices)
    @test all(length(face.support_indices) == length(interval)^2 for face in shell.faces)
    for left in 1:length(shell.faces), right in (left + 1):length(shell.faces)
        @test isempty(intersect(shell.faces[left].support_indices, shell.faces[right].support_indices))
    end
    @test norm(packet.overlap - I, Inf) < 1.0e-10
    @test packet.overlap ≈ direct_overlap atol = 1.0e-12 rtol = 1.0e-12
    @test packet.kinetic ≈ transpose(packet.kinetic) atol = 1.0e-10 rtol = 1.0e-10
    @test packet.position_x ≈ transpose(packet.position_x) atol = 1.0e-10 rtol = 1.0e-10
    @test packet.position_y ≈ transpose(packet.position_y) atol = 1.0e-10 rtol = 1.0e-10
    @test packet.position_z ≈ transpose(packet.position_z) atol = 1.0e-10 rtol = 1.0e-10
    @test !hasproperty(packet, :gaussian_terms)
    @test !hasproperty(packet, :pair_terms)
    @test !hasproperty(packet, :term_storage)
    @test !isnothing(packet.gaussian_sum)
    @test !isnothing(packet.pair_sum)
    support_coefficients = GaussletBases._nested_support_coefficient_slice(
        shell.coefficient_matrix,
        shell.support_indices,
    )
    support_workspace, contraction_scratch = GaussletBases._nested_support_reference_workspaces(
        support_coefficients,
        length(shell.support_indices),
        size(shell.coefficient_matrix, 2),
    )
    support_axes = GaussletBases._nested_support_axes(shell.support_states)
    support_weights = GaussletBases._nested_support_weights(shell.support_states, bundle.pgdg_intermediate.weights)
    fixed_weights = vec(transpose(support_coefficients) * support_weights)
    weighted_support_coefficients = support_coefficients .* reshape(1.0 ./ fixed_weights, 1, :)
    gaussian_reference = GaussletBases._nested_support_reference_gaussian_sum(
        support_axes,
        support_coefficients,
        support_workspace,
        contraction_scratch,
        term_coefficients,
        bundle.pgdg_intermediate.gaussian_factor_terms,
        bundle.pgdg_intermediate.gaussian_factor_terms,
        bundle.pgdg_intermediate.gaussian_factor_terms,
    )
    pair_reference = GaussletBases._nested_support_reference_pair_sum(
        support_axes,
        weighted_support_coefficients,
        support_workspace,
        contraction_scratch,
        term_coefficients,
        bundle.pgdg_intermediate.pair_factor_terms_raw,
        bundle.pgdg_intermediate.pair_factor_terms_raw,
        bundle.pgdg_intermediate.pair_factor_terms_raw,
    )
    @test packet.gaussian_sum ≈ gaussian_reference atol = 1.0e-10 rtol = 1.0e-10
    @test packet.pair_sum ≈ pair_reference atol = 1.0e-10 rtol = 1.0e-10

    for (face, columns) in zip(shell.faces, shell.face_column_ranges)
        mean_value =
            face.fixed_axis == :x ? sum(diag(packet.position_x)[columns]) / length(columns) :
            face.fixed_axis == :y ? sum(diag(packet.position_y)[columns]) / length(columns) :
            sum(diag(packet.position_z)[columns]) / length(columns)
        if face.fixed_side == :low
            @test mean_value < 0.0
        else
            @test mean_value > 0.0
        end
    end
end

@testset "Cartesian nested fixed-block QW-PGDG adapter" begin
    CCS = GaussletBases.CartesianCarriedSpaces
    QWCS = GaussletBases.CartesianQWOperatorCarriedSpaces
    (
        basis,
        bundle,
        shell,
        fixed_block,
        shell_plus_core,
        fixed_block_shell_plus_core,
        legacy,
        baseline,
        nested,
        nested_shell_plus_core,
        baseline_check,
        nested_check,
        nested_shell_plus_core_check,
    ) = _nested_qiu_white_nearest_fixture()

    @test fixed_block isa GaussletBases._NestedFixedBlock3D
    @test fixed_block.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test fixed_block.parent_basis === basis
    @test fixed_block.shell === shell
    @test size(fixed_block.coefficient_matrix) == size(shell.coefficient_matrix)
    @test size(fixed_block.overlap) == (54, 54)
    @test size(fixed_block.kinetic) == (54, 54)
    @test size(fixed_block.fixed_centers) == (54, 3)
    @test length(fixed_block.support_indices) == length(shell.support_indices)
    @test norm(fixed_block.overlap - I, Inf) < 1.0e-10
    @test fixed_block.overlap ≈ shell.packet.overlap atol = 0.0 rtol = 0.0
    @test !hasproperty(fixed_block, :gaussian_terms)
    @test !hasproperty(fixed_block, :pair_terms)
    @test !hasproperty(fixed_block, :term_storage)
    @test fixed_block.gaussian_sum ≈ shell.packet.gaussian_sum atol = 0.0 rtol = 0.0
    @test fixed_block.pair_sum ≈ shell.packet.pair_sum atol = 0.0 rtol = 0.0

    @test baseline.interaction_treatment == :ggt_nearest
    @test nested.interaction_treatment == :ggt_nearest
    @test nested.basis === fixed_block
    @test nested.gausslet_count == size(fixed_block.overlap, 1)
    @test nested.residual_count >= 1
    @test baseline.gausslet_count == length(bundle.pgdg_intermediate.centers)^3
    @test norm(nested.overlap - I, Inf) < 1.0e-10
    @test nested_check.overlap_error < 1.0e-10
    @test isfinite(nested_check.orbital_energy)
    @test isfinite(nested_check.vee_expectation)
    @test nested_check.orbital_energy < 0.0
    @test nested_check.vee_expectation > 0.0
    @test any(orbital.kind == :nested_fixed for orbital in orbitals(nested))
    @test all(startswith(orbital.label, "nf") for orbital in orbitals(nested)[1:nested.gausslet_count])

    overlap_before = copy(nested.overlap)
    one_body_before = copy(nested.one_body_hamiltonian)
    interaction_before = copy(nested.interaction_matrix)
    nested_sidecar = QWCS.cartesian_qw_operator_carried_space_sidecar(nested)
    nested_carried = QWCS.qw_operator_carried_space(nested_sidecar)
    nested_representation = QWCS.qw_operator_basis_representation(nested_sidecar)
    nested_diagnostics = QWCS.qw_operator_carried_space_diagnostics(nested_sidecar)
    @test QWCS.qw_operator_carried_space_provenance(nested_sidecar).input_kind ==
        :nested_fixed_block_operator
    @test nested_carried isa CCS.CartesianCarriedSpace3D
    @test nested_representation isa CartesianBasisRepresentation3D
    @test nested_diagnostics.operator_dimension == size(nested.overlap, 1)
    @test nested_diagnostics.operator_gausslet_count == nested.gausslet_count
    @test nested_diagnostics.operator_residual_count == nested.residual_count
    @test nested_diagnostics.raw_parent_dimension == size(nested.raw_to_final, 1)
    @test nested_diagnostics.carried_dimension == size(fixed_block.overlap, 1)
    @test nested_diagnostics.carried_dimension_matches_operator_gausslet_count
    @test nested_diagnostics.operator_representation_matches_operator_dimension
    @test nested_diagnostics.carried_has_contracted_parent
    @test nested_diagnostics.carried_has_staged_sidecar == false
    @test nested_diagnostics.dense_parent_matrix_used == false
    @test nested_diagnostics.heavy_metric_packet_built == false
    @test nested.overlap == overlap_before
    @test nested.one_body_hamiltonian == one_body_before
    @test nested.interaction_matrix == interaction_before

    nested_build_source = QWCS.cartesian_qw_operator_build_source(
        fixed_block,
        legacy;
        Z = 2.0,
        interaction_treatment = nested.interaction_treatment,
        gausslet_backend = nested.gausslet_backend,
    )
    nested_build_diagnostics =
        QWCS.operator_build_source_diagnostics(nested_build_source)
    @test QWCS.operator_build_source_provenance(nested_build_source).input_kind ==
        :atomic_nested_fixed_block_input
    @test nested_build_source.basis_family == :one_center_atomic
    @test nested_build_source.carried_space_kind == :nested_fixed_block
    @test nested_build_source.gausslet_backend == nested.gausslet_backend
    @test nested_build_source.interaction_treatment == nested.interaction_treatment
    @test nested_build_source.nuclear_term_storage == nested.nuclear_term_storage
    @test nested_build_diagnostics.carried_dimension ==
        nested_diagnostics.carried_dimension
    @test nested_build_diagnostics.carried_has_contracted_parent ==
        nested_diagnostics.carried_has_contracted_parent
    @test nested_build_diagnostics.carried_has_staged_sidecar ==
        nested_diagnostics.carried_has_staged_sidecar
    @test nested_build_diagnostics.dense_parent_matrix_used == false
    @test nested_build_diagnostics.heavy_metric_packet_built == false
    @test nested_build_diagnostics.operator_built == false

    nested_record =
        QWCS.cartesian_qw_operator_construction_record(nested_build_source, nested)
    nested_record_diagnostics =
        QWCS.qw_operator_construction_record_diagnostics(nested_record)
    @test nested_record_diagnostics.source_sidecar_agree
    @test isempty(nested_record_diagnostics.mismatch_fields)
    @test isempty(nested_record_diagnostics.ambiguous_mismatch_fields)
    @test :gausslet_backend in nested_record_diagnostics.compared_fields
    @test :nuclear_term_storage in nested_record_diagnostics.compared_fields
    @test :nuclear_charges in nested_record_diagnostics.compared_fields
    @test :carried_dimension in nested_record_diagnostics.compared_fields
    @test :carried_parent_axis_counts in nested_record_diagnostics.compared_fields
    @test :carried_parent_dimension in nested_record_diagnostics.compared_fields
    @test :carried_representation_basis_kind in nested_record_diagnostics.compared_fields
    @test :carried_representation_parent_kind in nested_record_diagnostics.compared_fields
    @test :carried_representation_final_dimension in nested_record_diagnostics.compared_fields
    @test :carried_axis_sharing in nested_record_diagnostics.compared_fields
    @test :carried_provenance_input_kind in nested_record_diagnostics.compared_fields
    @test :carried_provenance_route_metadata in nested_record_diagnostics.compared_fields
    @test nested_record_diagnostics.source_basis_family == :one_center_atomic
    @test nested_record_diagnostics.source_carried_space_kind == :nested_fixed_block
    @test nested_record_diagnostics.sidecar_input_kind == :nested_fixed_block_operator
    @test nested_record_diagnostics.source_parent_axis_counts ==
        nested_record_diagnostics.sidecar_parent_axis_counts
    @test nested_record_diagnostics.source_parent_dimension ==
        nested_record_diagnostics.sidecar_parent_dimension
    @test :coefficient_matrix_values in
        nested_record_diagnostics.intentionally_not_compared
    @test :overlap_matrix_values in
        nested_record_diagnostics.intentionally_not_compared
    @test nested_record_diagnostics.numerical_outputs_changed == false
    @test nested_record_diagnostics.dense_parent_matrix_used == false
    @test nested_record_diagnostics.heavy_metric_packet_built == false
    @test nested_record_diagnostics.operator_built == false
    @test QWCS.qw_operator_construction_record_sidecar(nested_record) isa
        QWCS.CartesianQWOperatorCarriedSpaceSidecar
    @test QWCS.qw_operator_construction_record_provenance(nested_record).source ==
        :cartesian_qw_operator_construction_record
    @test nested.overlap == overlap_before
    @test nested.one_body_hamiltonian == one_body_before
    @test nested.interaction_matrix == interaction_before

    @test shell_plus_core isa GaussletBases._CartesianNestedShellPlusCore3D
    @test fixed_block_shell_plus_core isa GaussletBases._NestedFixedBlock3D
    @test fixed_block_shell_plus_core.parent_basis === basis
    @test fixed_block_shell_plus_core.shell === shell_plus_core
    inner_len = length(basis) - 2
    @test first(shell_plus_core.core_column_range) == 1
    @test last(shell_plus_core.core_column_range) == length(shell_plus_core.core_indices)
    @test length(shell_plus_core.core_indices) == inner_len^3
    @test isempty(intersect(shell_plus_core.core_indices, shell.support_indices))
    @test size(fixed_block_shell_plus_core.overlap, 1) == length(shell_plus_core.core_indices) + size(shell.coefficient_matrix, 2)
    @test norm(fixed_block_shell_plus_core.overlap - I, Inf) < 1.0e-10
    @test nested_shell_plus_core.interaction_treatment == :ggt_nearest
    @test nested_shell_plus_core.basis === fixed_block_shell_plus_core
    @test nested_shell_plus_core.gausslet_count == size(fixed_block_shell_plus_core.overlap, 1)
    @test nested_shell_plus_core.residual_count >= 1
    @test nested_shell_plus_core_check.overlap_error < 1.0e-10
    @test isfinite(nested_shell_plus_core_check.orbital_energy)
    @test isfinite(nested_shell_plus_core_check.vee_expectation)
    @test nested_shell_plus_core_check.orbital_energy < 0.0
    @test nested_shell_plus_core_check.vee_expectation > 0.0
    @test abs(nested_shell_plus_core_check.vee_expectation - baseline_check.vee_expectation) < abs(nested_check.vee_expectation - baseline_check.vee_expectation)
    @test abs(nested_shell_plus_core_check.orbital_energy - baseline_check.orbital_energy) < abs(nested_check.orbital_energy - baseline_check.orbital_energy)
    @test abs(nested_shell_plus_core_check.vee_expectation - baseline_check.vee_expectation) < 1.0e-4
    @test abs(nested_shell_plus_core_check.orbital_energy - baseline_check.orbital_energy) < 1.0e-4
end

function _one_center_atomic_full_parent_contract_fixture(;
    Z::Float64 = 2.0,
    d::Float64 = 0.15,
    count::Int = 19,
    nside::Int = 7,
    tail_spacing::Float64 = 10.0,
)
    key = Symbol(
        :one_center_atomic_full_parent_contract_fixture,
        round(Int, 1000 * Z),
        round(Int, 1000 * d),
        count,
        nside,
    )
    return _cached_fixture(key, () -> begin
        basis = build_basis(MappedUniformBasisSpec(:G10;
            count,
            mapping = white_lindsey_atomic_mapping(; Z, d, tail_spacing),
            reference_spacing = 1.0,
        ))
        expansion = coulomb_gaussian_expansion(doacc = false)
        sequence = build_one_center_atomic_full_parent_shell_sequence(
            basis;
            exponents = expansion.exponents,
            gausslet_backend = :numerical_reference,
            refinement_levels = 0,
            nside,
        )
        audit = GaussletBases._nested_shell_sequence_contract_audit(sequence, (count, count, count))
        (basis, sequence, audit)
    end)
end

function _one_center_atomic_legacy_profile_contract_fixture(;
    Z::Float64 = 2.0,
    d::Float64 = 0.2,
    count::Int = 15,
    nside::Int = 5,
    working_box::UnitRange{Int} = 2:14,
    tail_spacing::Float64 = 10.0,
)
    key = Symbol(
        :one_center_atomic_legacy_profile_contract_fixture,
        round(Int, 1000 * Z),
        round(Int, 1000 * d),
        count,
        nside,
        first(working_box),
        last(working_box),
    )
    return _cached_fixture(key, () -> begin
        basis = build_basis(MappedUniformBasisSpec(:G10;
            count,
            mapping = white_lindsey_atomic_mapping(; Z, d, tail_spacing),
            reference_spacing = 1.0,
        ))
        expansion = coulomb_gaussian_expansion(doacc = false)
        sequence = build_one_center_atomic_legacy_profile_shell_sequence(
            basis;
            exponents = expansion.exponents,
            gausslet_backend = :numerical_reference,
            refinement_levels = 0,
            working_box,
            nside,
        )
        diagnostics = one_center_atomic_nested_structure_diagnostics(
            sequence;
            parent_side_count = count,
            nside,
        )
        ownership = GaussletBases._nested_shell_sequence_piece_ownership_audit(sequence)
        (basis, sequence, diagnostics, ownership)
    end)
end

function _ne_repo_v6z_sp_basis_text()
    return "#BASIS SET: Ne repo-v6z-sp\n" *
           "Ne    S\n" *
           "      9.024000e+05           5.510000e-06\n" *
           "      1.351000e+05           4.282000e-05\n" *
           "      3.075000e+04           2.251400e-04\n" *
           "      8.710000e+03           9.501600e-04\n" *
           "      2.842000e+03           3.447190e-03\n" *
           "      1.026000e+03           1.112545e-02\n" *
           "      4.001000e+02           3.220568e-02\n" *
           "      1.659000e+02           8.259891e-02\n" *
           "      7.221000e+01           1.799056e-01\n" *
           "      3.266000e+01           3.060521e-01\n" *
           "      1.522000e+01           3.401256e-01\n" *
           "      7.149000e+00           1.761682e-01\n" *
           "      2.957000e+00           2.101527e-02\n" *
           "      1.335000e+00          -5.074500e-04\n" *
           "      5.816000e-01           1.057850e-03\n" *
           "      2.463000e-01          -5.988000e-05\n" *
           "Ne    S\n" *
           "      7.149000e+00           1.000000e+00\n" *
           "Ne    S\n" *
           "      2.957000e+00           1.000000e+00\n" *
           "Ne    S\n" *
           "      1.335000e+00           1.000000e+00\n" *
           "Ne    S\n" *
           "      9.024000e+05          -1.290000e-06\n" *
           "      1.351000e+05          -1.005000e-05\n" *
           "      3.075000e+04          -5.293000e-05\n" *
           "      8.710000e+03          -2.231200e-04\n" *
           "      2.842000e+03          -8.133800e-04\n" *
           "      1.026000e+03          -2.632300e-03\n" *
           "      4.001000e+02          -7.759100e-03\n" *
           "      1.659000e+02          -2.045277e-02\n" *
           "      7.221000e+01          -4.797505e-02\n" *
           "      3.266000e+01          -9.340086e-02\n" *
           "      1.522000e+01          -1.427721e-01\n" *
           "      7.149000e+00          -1.022908e-01\n" *
           "      2.957000e+00           1.587858e-01\n" *
           "      1.335000e+00           4.494079e-01\n" *
           "      5.816000e-01           4.334854e-01\n" *
           "      2.463000e-01           1.215757e-01\n" *
           "Ne    S\n" *
           "      5.816000e-01           1.000000e+00\n" *
           "Ne    S\n" *
           "      2.463000e-01           1.000000e+00\n" *
           "Ne    P\n" *
           "      4.281000e+00           1.000000e+00\n" *
           "Ne    P\n" *
           "      1.915000e+00           1.000000e+00\n" *
           "Ne    P\n" *
           "      8.156000e+02           1.837600e-04\n" *
           "      1.933000e+02           1.585090e-03\n" *
           "      6.260000e+01           8.414640e-03\n" *
           "      2.361000e+01           3.220033e-02\n" *
           "      9.762000e+00           9.396390e-02\n" *
           "      4.281000e+00           2.004808e-01\n" *
           "      1.915000e+00           3.031137e-01\n" *
           "      8.476000e-01           3.297578e-01\n" *
           "      3.660000e-01           2.366743e-01\n" *
           "      1.510000e-01           6.911689e-02\n" *
           "Ne    P\n" *
           "      8.476000e-01           1.000000e+00\n" *
           "Ne    P\n" *
           "      3.660000e-01           1.000000e+00\n" *
           "Ne    P\n" *
           "      1.510000e-01           1.000000e+00\n" *
           "END\n"
end

function _one_center_atomic_legacy_profile_ne_residual_completion_fixture()
    return _cached_fixture(:one_center_atomic_legacy_profile_ne_residual_completion_fixture, () -> begin
        overlap_only_expansion = CoulombGaussianExpansion(
            [0.0],
            [1.0];
            del = 1.0,
            s = 1.0,
            c = 1.0,
            maxu = 1.0,
        )
        mktemp() do path, io
            write(io, _ne_repo_v6z_sp_basis_text())
            close(io)

            basis = build_basis(MappedUniformBasisSpec(:G10;
                count = 29,
                mapping = white_lindsey_atomic_mapping(Z = 10.0, d = 0.03, tail_spacing = 10.0),
                reference_spacing = 1.0,
            ))
            bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
                basis;
                exponents = overlap_only_expansion.exponents,
                center = 0.0,
                backend = :numerical_reference,
            )
            fixed_block = one_center_atomic_legacy_profile_fixed_block(
                bundle;
                working_box = 2:28,
                nside = 7,
            )
            supplement = legacy_atomic_gaussian_supplement(
                "Ne",
                "repo-v6z-sp";
                lmax = 1,
                basisfile = path,
            )
            supplement3d = GaussletBases._atomic_cartesian_shell_supplement_3d(supplement)
            blocks = GaussletBases._qwrg_atomic_cartesian_blocks_3d(
                bundle,
                supplement3d,
                overlap_only_expansion,
            )
            overlap_fg = GaussletBases._qwrg_contract_parent_ga_matrix(
                fixed_block.coefficient_matrix,
                blocks.overlap_ga,
            )
            near_null = diagnose_qwrg_residual_space(
                fixed_block.overlap,
                overlap_fg,
                blocks.overlap_aa;
                keep_policy = :near_null_only,
                keep_abs_tol = GaussletBases._qwrg_atomic_residual_keep_tol(),
                accept_tol = GaussletBases._qwrg_atomic_residual_accept_tol(),
            )
            near_null_data = GaussletBases._qwrg_residual_space(
                fixed_block.overlap,
                overlap_fg,
                blocks.overlap_aa;
                keep_policy = :near_null_only,
                keep_abs_tol = GaussletBases._qwrg_atomic_residual_keep_tol(),
                accept_tol = GaussletBases._qwrg_atomic_residual_accept_tol(),
            )
            near_null_total_basis = size(near_null_data.raw_to_final, 2)
            legacy_alias = diagnose_qwrg_residual_space(
                fixed_block.overlap,
                overlap_fg,
                blocks.overlap_aa;
                keep_policy = :legacy_profile,
                keep_abs_tol = GaussletBases._qwrg_atomic_residual_keep_tol(),
                accept_tol = GaussletBases._qwrg_atomic_residual_accept_tol(),
            )
            return (
                fixed_gausslet_count = size(fixed_block.overlap, 1),
                supplement_count = length(supplement3d.orbitals),
                near_null = near_null,
                near_null_data = near_null_data,
                near_null_total_basis = near_null_total_basis,
                legacy_alias = legacy_alias,
            )
        end
    end)
end

function _one_center_atomic_ns9_legacy_profile_qw_fixture()
    return _cached_fixture(:one_center_atomic_ns9_legacy_profile_qw_fixture, () -> begin
        mktemp() do path, io
            write(io, _ne_repo_v6z_sp_basis_text())
            close(io)

            basis = build_basis(MappedUniformBasisSpec(
                :G10;
                count = 29,
                mapping = white_lindsey_atomic_mapping(Z = 10.0, d = 0.03, tail_spacing = 10.0),
                reference_spacing = 1.0,
            ))
            expansion = coulomb_gaussian_expansion(doacc = false)
            fixed_block = one_center_atomic_legacy_profile_fixed_block(
                basis;
                expansion = expansion,
                working_box = 2:28,
                nside = 9,
            )
            supplement = legacy_atomic_gaussian_supplement(
                "Ne",
                "repo-v6z-sp";
                lmax = 1,
                basisfile = path,
            )
            supplement3d = GaussletBases._atomic_cartesian_shell_supplement_3d(supplement)
            bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
                basis;
                exponents = expansion.exponents,
                center = 0.0,
                backend = :numerical_reference,
            )
            blocks = GaussletBases._qwrg_atomic_cartesian_blocks_3d(
                bundle,
                supplement3d,
                expansion,
            )
            overlap_fg = GaussletBases._qwrg_contract_parent_ga_matrix(
                fixed_block.coefficient_matrix,
                blocks.overlap_ga,
            )
            residual_data = GaussletBases._qwrg_residual_space(
                fixed_block.overlap,
                overlap_fg,
                blocks.overlap_aa;
                keep_policy = :near_null_only,
                keep_abs_tol = GaussletBases._qwrg_atomic_residual_keep_tol(),
                accept_tol = GaussletBases._qwrg_atomic_residual_accept_tol(),
            )
            operators = ordinary_cartesian_qiu_white_operators(
                fixed_block,
                supplement;
                expansion,
                Z = 10.0,
                interaction_treatment = :ggt_nearest,
                residual_keep_policy = :near_null_only,
            )
            return (
                fixed_block = fixed_block,
                residual_data = residual_data,
                operators = operators,
            )
        end
    end)
end

@testset "One-center atomic full-parent nested contract" begin
    basis, sequence, audit = _one_center_atomic_full_parent_contract_fixture()
    count = length(basis)
    diagnostics = one_center_atomic_nested_structure_diagnostics(
        sequence;
        parent_side_count = count,
        nside = 7,
    )
    common_contract = GaussletBases._nested_glass_box_contract(diagnostics)
    count_only_27 = one_center_atomic_nested_structure_diagnostics(27; nside = 7)
    count_only_29 = one_center_atomic_nested_structure_diagnostics(29; nside = 7)
    count_only_5 = one_center_atomic_nested_structure_diagnostics(15; nside = 5)

    @test sequence isa GaussletBases._CartesianNestedShellSequence3D
    @test sequence.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test sequence.working_box == (1:count, 1:count, 1:count)
    @test audit.full_parent_working_box
    @test audit.support_count == count^3
    @test audit.expected_support_count == count^3
    @test audit.missing_row_count == 0
    @test audit.ownership_group_count_min == 1
    @test audit.ownership_group_count_max == 1
    @test audit.ownership_unowned_row_count == 0
    @test audit.ownership_multi_owned_row_count == 0

    @test diagnostics.parent_side_count == count
    @test diagnostics.working_box_side_count == count
    @test diagnostics.nside == 7
    @test diagnostics.core_side_count == 7
    @test diagnostics.shell_layer_count == 6
    @test diagnostics.expected_shell_increment == 7^3 - 5^3
    @test diagnostics.expected_shell_increment == 218
    @test diagnostics.expected_face_retained_count == 6 * 5^2
    @test diagnostics.expected_edge_retained_count == 12 * 5
    @test diagnostics.expected_corner_retained_count == 8
    @test diagnostics.total_face_retained_count == 6 * (6 * 5^2)
    @test diagnostics.total_edge_retained_count == 6 * (12 * 5)
    @test diagnostics.total_corner_retained_count == 6 * 8
    @test diagnostics.total_expected_gausslet_count == 7^3 + 6 * 218
    @test diagnostics.total_actual_gausslet_count == 7^3 + 6 * 218
    @test diagnostics.total_actual_gausslet_count == size(sequence.coefficient_matrix, 2)
    @test diagnostics.layers_match_expected
    @test length(diagnostics.layer_structures) == 6
    @test all(layer.face_retained_count == 150 for layer in diagnostics.layer_structures)
    @test all(layer.edge_retained_count == 60 for layer in diagnostics.layer_structures)
    @test all(layer.corner_retained_count == 8 for layer in diagnostics.layer_structures)
    @test all(layer.retained_dimension == 218 for layer in diagnostics.layer_structures)
    @test common_contract.fixed_dimension == diagnostics.total_actual_gausslet_count
    @test common_contract.contract_audit === diagnostics
    @test common_contract.layer_dimensions ==
          Int[layer.retained_dimension for layer in diagnostics.layer_structures]
    @test common_contract.layer_provenance ==
          GaussletBases._CartesianNestedShellLayerProvenance3D[
              shell.provenance for shell in sequence.shell_layers
          ]
    @test common_contract.leaf_count === nothing
    @test common_contract.layer_provenance[1].source_box == sequence.working_box
    @test common_contract.layer_provenance[1].source_point_count == 19^3 - 17^3
    @test common_contract.layer_provenance[end].next_inner_box == ((7:13), (7:13), (7:13))
    report_text = one_center_atomic_nested_structure_report(diagnostics)
    @test occursin("layer_1_source_box = ", report_text)
    @test occursin("layer_1_next_inner_box = ", report_text)
    @test occursin("layer_1_source_point_count = ", report_text)
    @test !occursin("leaf_count", report_text)

    @test count_only_5.expected_shell_increment == 5^3 - 3^3
    @test count_only_5.expected_shell_increment == 98
    @test count_only_5.layer_structures[1].provenance.source_point_count == 15^3 - 13^3

    @test count_only_27.parent_side_count == 27
    @test count_only_27.working_box_side_count == 27
    @test count_only_27.nside == 7
    @test count_only_27.core_side_count == 7
    @test count_only_27.shell_layer_count == 10
    @test count_only_27.expected_shell_increment == 218
    @test count_only_27.total_expected_gausslet_count == 343 + 10 * 218
    @test count_only_27.total_actual_gausslet_count == 343 + 10 * 218
    @test count_only_27.total_actual_gausslet_count == 2523
    @test count_only_27.layers_match_expected

    @test count_only_29.parent_side_count == 29
    @test count_only_29.working_box_side_count == 29
    @test count_only_29.shell_layer_count == 11
    @test count_only_29.expected_shell_increment == 218
    @test count_only_29.total_actual_gausslet_count == 343 + 11 * 218
end

@testset "One-center atomic legacy-profile nested contract" begin
    basis, sequence, diagnostics, ownership = _one_center_atomic_legacy_profile_contract_fixture()
    common_contract = GaussletBases._nested_glass_box_contract(diagnostics)
    range_groups = UnitRange{Int}[sequence.core_column_range]
    append!(range_groups, sequence.layer_column_ranges)
    support_group_counts = Int[]
    for row in sequence.support_indices
        nzcols = findall(!iszero, @view sequence.coefficient_matrix[row, :])
        touched_groups = 0
        for range in range_groups
            any(col -> col in range, nzcols) && (touched_groups += 1)
        end
        push!(support_group_counts, touched_groups)
    end

    @test sequence isa GaussletBases._CartesianNestedShellSequence3D
    @test sequence.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test sequence.working_box == (2:14, 2:14, 2:14)
    @test length(sequence.support_indices) == 13^3
    @test ownership.min_group_count == 0
    @test ownership.max_group_count == 1
    @test ownership.unowned_row_count == length(basis)^3 - 13^3
    @test ownership.multi_owned_row_count == 0
    @test minimum(support_group_counts) == 1
    @test maximum(support_group_counts) == 1
    @test diagnostics.parent_side_count == length(basis)
    @test diagnostics.working_box_side_count == 13
    @test diagnostics.nside == 5
    @test diagnostics.core_side_count == 5
    @test diagnostics.shell_layer_count == 4
    @test diagnostics.expected_shell_increment == 98
    @test diagnostics.total_actual_gausslet_count == 5^3 + 4 * 98
    @test diagnostics.layers_match_expected
    @test common_contract.fixed_dimension == diagnostics.total_actual_gausslet_count
    @test common_contract.contract_audit === diagnostics
    @test common_contract.layer_dimensions == [98, 98, 98, 98]
    @test common_contract.layer_provenance ==
          GaussletBases._CartesianNestedShellLayerProvenance3D[
              shell.provenance for shell in sequence.shell_layers
          ]
    @test common_contract.leaf_count === nothing
    @test common_contract.layer_provenance[1].source_box == ((2:14), (2:14), (2:14))
    @test common_contract.layer_provenance[end].next_inner_box == ((6:10), (6:10), (6:10))

    count_only_legacy_ne = one_center_atomic_nested_structure_diagnostics(
        29;
        working_box_side_count = 27,
        nside = 7,
    )
    count_only_modern_ne = one_center_atomic_nested_structure_diagnostics(29; nside = 7)
    @test count_only_legacy_ne.parent_side_count == 29
    @test count_only_legacy_ne.working_box_side_count == 27
    @test count_only_legacy_ne.shell_layer_count == 10
    @test count_only_legacy_ne.expected_shell_increment == 218
    @test count_only_legacy_ne.total_actual_gausslet_count == 2523
    @test count_only_modern_ne.working_box_side_count == 29
    @test count_only_modern_ne.shell_layer_count == 11
    @test count_only_modern_ne.total_actual_gausslet_count == 2741
end

@testset "One-center atomic fixed-block timing surface" begin
    function _timing_labels(report::GaussletBases.TimeG.TimingReport)
        labels = String[]
        function _visit(node)
            push!(labels, node.label)
            foreach(_visit, node.children)
        end
        foreach(_visit, report.roots)
        return labels
    end

    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 13,
            mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)

    timed_full = one_center_atomic_full_parent_fixed_block(
        basis;
        exponents = expansion.exponents,
        nside = 5,
        timing = true,
    )
    timed_legacy = one_center_atomic_legacy_profile_fixed_block(
        basis;
        exponents = expansion.exponents,
        working_box = 2:12,
        nside = 5,
        timing = true,
    )

    @test timed_full isa GaussletBases.TimedNestedFixedBlockBuild
    @test timed_legacy isa GaussletBases.TimedNestedFixedBlockBuild
    @test timed_full.timings isa GaussletBases.TimeG.TimingReport
    @test timed_legacy.timings isa GaussletBases.TimeG.TimingReport
    @test timed_full.fixed_block.shell.working_box == (1:13, 1:13, 1:13)
    @test timed_legacy.fixed_block.shell.working_box == (2:12, 2:12, 2:12)
    @test !hasproperty(timed_full.fixed_block, :gaussian_terms)
    @test !hasproperty(timed_full.fixed_block, :pair_terms)
    @test !hasproperty(timed_full.fixed_block, :term_storage)
    @test !isnothing(timed_full.fixed_block.gaussian_sum)
    @test !isnothing(timed_full.fixed_block.pair_sum)
    @test !hasproperty(timed_legacy.fixed_block, :gaussian_terms)
    @test !hasproperty(timed_legacy.fixed_block, :pair_terms)
    @test !hasproperty(timed_legacy.fixed_block, :term_storage)
    @test !isnothing(timed_legacy.fixed_block.gaussian_sum)
    @test !isnothing(timed_legacy.fixed_block.pair_sum)
    @test norm(timed_full.fixed_block.overlap - I, Inf) < 1.0e-10
    @test norm(timed_legacy.fixed_block.overlap - I, Inf) < 1.0e-10
    full_labels = _timing_labels(timed_full.timings)
    legacy_labels = _timing_labels(timed_legacy.timings)
    @test "fixed_block.total" in full_labels
    @test "fixed_block.parent_bundle" in full_labels
    @test "fixed_block.sequence_build" in full_labels
    @test "fixed_block.adapter" in full_labels
    @test "diatomic.packet.total" in full_labels
    @test "diatomic.packet.gaussian_terms" in full_labels
    @test "diatomic.packet.pair_terms" in full_labels
    @test "diatomic.packet.total" in legacy_labels

    full_report = nested_fixed_block_timing_report(timed_full)
    legacy_report = nested_fixed_block_timing_report(timed_legacy.timings)
    @test occursin("fixed_block.total", full_report)
    @test occursin("shell_layer.nonpacket", full_report)
    @test occursin("sequence_merge.nonpacket", full_report)
    @test occursin("diatomic.packet.gaussian_terms", full_report)
    @test occursin("diatomic.packet.pair_terms", full_report)
    @test occursin("diatomic.packet.total", legacy_report)
end

@testset "Global timing macro surface" begin
    old_config = GaussletBases.TimeG._TIMING_CONFIG[]
    try
        @test timing_enabled() == GaussletBases.TimeG.timing_enabled()
        @test timing_live_enabled() == GaussletBases.TimeG.timing_live_enabled()

        reset_timing_report!()
        set_timing!(true)
        set_timing_live!(false)
        set_timing_thresholds!(expand = 0.0, drop = 0.0)

        @timeg "outer" begin
            sleep(0.002)
            @timeg "inner" begin
                sleep(0.001)
            end
        end

        report = current_timing_report()
        @test report isa GaussletBases.TimeG.TimingReport
        @test length(report.roots) == 1
        root = only(report.roots)
        @test root.label == "outer"
        @test root.elapsed_seconds > 0.0
        @test root.self_seconds >= 0.0
        @test root.call_count == 1
        @test length(root.children) == 1
        child = only(root.children)
        @test child.label == "inner"
        @test child.elapsed_seconds > 0.0

        rendered = timing_report(report)
        @test occursin("GaussletBases timing report", rendered)
        @test occursin("outer", rendered)
        @test occursin("inner", rendered)

        live_path, live_io = mktemp()
        close(live_io)
        try
            open(live_path, "w") do io
                redirect_stdout(io) do
                    reset_timing_report!()
                    set_timing!(true)
                    set_timing_live!(true)
                    set_timing_thresholds!(expand = 0.0, drop = 0.0)
                    @timeg "live outer" begin
                        @timeg "live inner" begin
                            sleep(0.001)
                        end
                    end
                end
            end
            live_output = read(live_path, String)
            @test occursin("live outer: ", live_output)
            @test occursin("live inner: ", live_output)
            @test occursin("seconds", live_output)
        finally
            rm(live_path; force = true)
        end

        reset_timing_report!()
        set_timing!(false)
        set_timing_live!(false)
        @timeg "disabled" begin
            sleep(0.001)
        end
        disabled_report = current_timing_report()
        @test isempty(disabled_report.roots)
    finally
        GaussletBases.TimeG._TIMING_CONFIG[] = old_config
        reset_timing_report!()
    end
end

@testset "One-center atomic compact fixed-block term storage" begin
    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 13,
            mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)

    compact_full = one_center_atomic_full_parent_fixed_block(
        basis;
        expansion,
        nside = 5,
    )
    compact_legacy = one_center_atomic_legacy_profile_fixed_block(
        basis;
        expansion,
        working_box = 2:12,
        nside = 5,
    )

    @test !hasproperty(compact_full, :term_storage)
    @test !hasproperty(compact_full, :gaussian_terms)
    @test !hasproperty(compact_full, :pair_terms)
    @test !isnothing(compact_full.gaussian_sum)
    @test !isnothing(compact_full.pair_sum)

    @test !hasproperty(compact_legacy, :term_storage)
    @test !hasproperty(compact_legacy, :gaussian_terms)
    @test !hasproperty(compact_legacy, :pair_terms)
    @test !isnothing(compact_legacy.gaussian_sum)
    @test !isnothing(compact_legacy.pair_sum)

    full_contract = GaussletBases._nested_glass_box_contract(
        one_center_atomic_nested_structure_diagnostics(compact_full; nside = 5),
    )
    legacy_contract = GaussletBases._nested_glass_box_contract(
        one_center_atomic_nested_structure_diagnostics(compact_legacy; nside = 5),
    )
    @test full_contract.fixed_dimension == size(compact_full.overlap, 1)
    @test legacy_contract.fixed_dimension == size(compact_legacy.overlap, 1)
    @test full_contract.leaf_count === nothing
    @test legacy_contract.leaf_count === nothing

end

@testset "Cartesian basis representation for direct-product QW bases" begin
    CP = GaussletBases.CartesianParentGaussletBases
    CCS = GaussletBases.CartesianCarriedSpaces
    QWCS = GaussletBases.CartesianQWOperatorCarriedSpaces
    basis, operators, _check = _bond_aligned_diatomic_qw_fixture()
    representation = basis_representation(basis)
    metadata = basis_metadata(representation)
    chain_basis, _chain_ops, _chain_diagnostics = _bond_aligned_homonuclear_chain_qw_fixture()
    square_basis, _square_ops, _square_diagnostics, _square_check =
        _axis_aligned_homonuclear_square_lattice_qw_fixture()
    chain_representation = basis_representation(chain_basis)
    square_representation = basis_representation(square_basis)

    @test representation isa CartesianBasisRepresentation3D
    @test metadata.basis_kind == :direct_product
    @test metadata.parent_kind == :cartesian_product_basis
    @test metadata.parent_axis_counts == (length(basis.basis_x), length(basis.basis_y), length(basis.basis_z))
    @test metadata.parent_dimension == prod(metadata.parent_axis_counts)
    @test metadata.final_dimension == prod(metadata.parent_axis_counts)
    @test metadata.axis_sharing == :shared_xy
    @test metadata.route_metadata.basis_family == :bond_aligned_diatomic
    @test metadata.route_metadata.bond_axis == basis.bond_axis
    @test metadata.route_metadata.nuclei == basis.nuclei
    @test representation.contraction_kind == :identity
    @test isnothing(representation.coefficient_matrix)
    @test isnothing(representation.support_indices)
    @test isnothing(representation.support_states)
    @test metadata.basis_labels == representation.parent_labels
    @test metadata.basis_centers == representation.parent_centers
    @test size(metadata.basis_centers, 1) == metadata.final_dimension
    @test size(metadata.basis_centers, 2) == 3
    @test chain_representation.metadata.basis_kind == :direct_product
    @test chain_representation.metadata.route_metadata.basis_family ==
        :bond_aligned_homonuclear_chain
    @test square_representation.metadata.basis_kind == :direct_product
    @test square_representation.metadata.route_metadata.basis_family ==
        :axis_aligned_homonuclear_square_lattice

    carried = CCS.cartesian_carried_space(basis)
    chain_carried = CCS.cartesian_carried_space(chain_basis)
    square_carried = CCS.cartesian_carried_space(square_basis)
    @test CCS.carried_space_parent(carried) isa CP.CartesianParentGaussletBasis3D
    @test isnothing(CCS.carried_space_contracted_parent(carried))
    @test CCS.carried_space_representation(carried) isa CartesianBasisRepresentation3D
    @test CCS.carried_space_diagnostics(carried).parent_axis_counts ==
        metadata.parent_axis_counts
    @test CCS.carried_space_diagnostics(carried).parent_dimension ==
        metadata.parent_dimension
    @test CCS.carried_space_diagnostics(carried).representation_final_dimension ==
        metadata.final_dimension
    @test CCS.carried_space_diagnostics(carried).has_contracted_parent == false
    @test CCS.carried_space_diagnostics(carried).dense_parent_matrix_used == false
    @test CCS.carried_space_diagnostics(carried).heavy_metric_packet_built == false
    @test CCS.carried_space_provenance(carried).input_kind ==
        :bond_aligned_diatomic_qw_basis
    @test CCS.carried_space_diagnostics(chain_carried).parent_axis_counts ==
        chain_representation.metadata.parent_axis_counts
    @test CCS.carried_space_provenance(chain_carried).input_kind ==
        :bond_aligned_homonuclear_chain_qw_basis
    @test CCS.carried_space_diagnostics(square_carried).parent_axis_counts ==
        square_representation.metadata.parent_axis_counts
    @test CCS.carried_space_provenance(square_carried).input_kind ==
        :axis_aligned_homonuclear_square_lattice_qw_basis

    overlap_before = copy(operators.overlap)
    one_body_before = copy(operators.one_body_hamiltonian)
    interaction_before = copy(operators.interaction_matrix)
    operator_sidecar = QWCS.cartesian_qw_operator_carried_space_sidecar(operators)
    operator_carried = QWCS.qw_operator_carried_space(operator_sidecar)
    operator_representation = QWCS.qw_operator_basis_representation(operator_sidecar)
    operator_diagnostics = QWCS.qw_operator_carried_space_diagnostics(operator_sidecar)
    @test QWCS.qw_operator_carried_space_provenance(operator_sidecar).input_kind ==
        :bond_aligned_direct_product_operator
    @test operator_carried isa CCS.CartesianCarriedSpace3D
    @test operator_representation isa CartesianBasisRepresentation3D
    @test operator_diagnostics.operator_dimension == size(operators.overlap, 1)
    @test operator_diagnostics.operator_gausslet_count == operators.gausslet_count
    @test operator_diagnostics.operator_residual_count == operators.residual_count
    @test operator_diagnostics.carried_dimension == metadata.final_dimension
    @test operator_diagnostics.carried_dimension_matches_operator_gausslet_count
    @test operator_diagnostics.operator_representation_matches_operator_dimension
    @test operator_diagnostics.carried_has_contracted_parent == false
    @test operator_diagnostics.carried_has_staged_sidecar == false
    @test operator_diagnostics.dense_parent_matrix_used == false
    @test operator_diagnostics.heavy_metric_packet_built == false
    @test operators.overlap == overlap_before
    @test operators.one_body_hamiltonian == one_body_before
    @test operators.interaction_matrix == interaction_before

    build_source = QWCS.cartesian_qw_operator_build_source(
        basis;
        nuclear_charges = operators.nuclear_charges,
        nuclear_term_storage = :auto,
        interaction_treatment = operators.interaction_treatment,
        gausslet_backend = operators.gausslet_backend,
    )
    build_diagnostics = QWCS.operator_build_source_diagnostics(build_source)
    @test QWCS.operator_build_source_carried_space(build_source) isa
        CCS.CartesianCarriedSpace3D
    @test QWCS.operator_build_source_provenance(build_source).input_kind ==
        :bond_aligned_direct_product_input
    @test build_source.basis_family == :bond_aligned_diatomic
    @test build_source.carried_space_kind == :direct_product
    @test build_source.nuclear_charges == operators.nuclear_charges
    @test build_source.gausslet_backend == operators.gausslet_backend
    @test build_source.interaction_treatment == operators.interaction_treatment
    @test build_source.nuclear_term_storage == operators.nuclear_term_storage
    @test build_diagnostics.carried_dimension == operator_diagnostics.carried_dimension
    @test build_diagnostics.carried_has_contracted_parent ==
        operator_diagnostics.carried_has_contracted_parent
    @test build_diagnostics.carried_has_staged_sidecar ==
        operator_diagnostics.carried_has_staged_sidecar
    @test build_diagnostics.dense_parent_matrix_used == false
    @test build_diagnostics.heavy_metric_packet_built == false
    @test build_diagnostics.operator_built == false

    construction_record =
        QWCS.cartesian_qw_operator_construction_record(build_source, operators)
    record_diagnostics =
        QWCS.qw_operator_construction_record_diagnostics(construction_record)
    @test record_diagnostics.source_sidecar_agree
    @test isempty(record_diagnostics.mismatch_fields)
    @test isempty(record_diagnostics.ambiguous_mismatch_fields)
    @test :operator_input_kind in record_diagnostics.compared_fields
    @test :gausslet_backend in record_diagnostics.compared_fields
    @test :interaction_treatment in record_diagnostics.compared_fields
    @test :nuclear_charges in record_diagnostics.compared_fields
    @test :carried_parent_axis_counts in record_diagnostics.compared_fields
    @test :carried_parent_dimension in record_diagnostics.compared_fields
    @test :carried_representation_basis_kind in record_diagnostics.compared_fields
    @test :carried_representation_parent_kind in record_diagnostics.compared_fields
    @test :carried_representation_final_dimension in record_diagnostics.compared_fields
    @test :carried_axis_sharing in record_diagnostics.compared_fields
    @test :carried_provenance_input_kind in record_diagnostics.compared_fields
    @test :carried_provenance_route_metadata in record_diagnostics.compared_fields
    @test :carried_has_staged_sidecar in record_diagnostics.compared_fields
    @test record_diagnostics.source_basis_family == :bond_aligned_diatomic
    @test record_diagnostics.source_carried_space_kind == :direct_product
    @test record_diagnostics.sidecar_input_kind == :bond_aligned_direct_product_operator
    @test record_diagnostics.source_parent_axis_counts ==
        record_diagnostics.sidecar_parent_axis_counts
    @test record_diagnostics.source_parent_dimension ==
        record_diagnostics.sidecar_parent_dimension
    @test :coefficient_matrix_values in
        record_diagnostics.intentionally_not_compared
    @test :interaction_matrix_values in
        record_diagnostics.intentionally_not_compared
    @test record_diagnostics.numerical_outputs_changed == false
    @test record_diagnostics.dense_parent_matrix_used == false
    @test record_diagnostics.heavy_metric_packet_built == false
    @test record_diagnostics.operator_built == false
    @test QWCS.qw_operator_construction_record_sidecar(construction_record) isa
        QWCS.CartesianQWOperatorCarriedSpaceSidecar
    @test QWCS.qw_operator_construction_record_provenance(construction_record).source ==
        :cartesian_qw_operator_construction_record

    construction_receipt = QWCS.cartesian_qw_operator_construction_receipt(
        basis;
        nuclear_charges = operators.nuclear_charges,
        nuclear_term_storage = operators.nuclear_term_storage,
        interaction_treatment = operators.interaction_treatment,
        gausslet_backend = operators.gausslet_backend,
    )
    receipt_operators =
        QWCS.qw_operator_construction_receipt_operators(construction_receipt)
    receipt_record =
        QWCS.qw_operator_construction_receipt_record(construction_receipt)
    receipt_diagnostics =
        QWCS.qw_operator_construction_receipt_diagnostics(construction_receipt)
    @test QWCS.qw_operator_construction_receipt_source(construction_receipt) isa
        QWCS.CartesianOperatorBuildSource3D
    @test receipt_record isa QWCS.CartesianQWOperatorConstructionRecord3D
    @test receipt_diagnostics.delegated_to_existing_builder
    @test receipt_diagnostics.builder == :ordinary_cartesian_qiu_white_operators
    @test receipt_diagnostics.source_sidecar_agree
    @test isempty(receipt_diagnostics.mismatch_fields)
    @test receipt_diagnostics.operator_built
    @test receipt_diagnostics.new_hamiltonian_kernel_used == false
    @test receipt_diagnostics.dense_parent_matrix_used == false
    @test receipt_diagnostics.heavy_metric_packet_built == false
    @test receipt_diagnostics.numerical_outputs_changed == false
    @test QWCS.qw_operator_construction_receipt_provenance(construction_receipt).source ==
        :cartesian_qw_operator_construction_receipt
    @test receipt_operators.overlap == operators.overlap
    @test receipt_operators.one_body_hamiltonian == operators.one_body_hamiltonian
    @test receipt_operators.interaction_matrix == operators.interaction_matrix
    @test receipt_operators.gausslet_count == operators.gausslet_count
    @test receipt_operators.residual_count == operators.residual_count
    @test receipt_operators.gausslet_backend == operators.gausslet_backend
    @test receipt_operators.interaction_treatment == operators.interaction_treatment
    @test receipt_operators.nuclear_term_storage == operators.nuclear_term_storage

    mismatched_source = QWCS.cartesian_qw_operator_build_source(
        basis;
        nuclear_charges = [
            operators.nuclear_charges[1] + 0.5,
            operators.nuclear_charges[2],
        ],
        nuclear_term_storage = :auto,
        interaction_treatment = operators.interaction_treatment,
        gausslet_backend = operators.gausslet_backend,
    )
    mismatch_record = QWCS.cartesian_qw_operator_construction_record(
        mismatched_source,
        operators;
        throw_on_mismatch = false,
    )
    mismatch_diagnostics =
        QWCS.qw_operator_construction_record_diagnostics(mismatch_record)
    @test !mismatch_diagnostics.source_sidecar_agree
    @test :nuclear_charges in mismatch_diagnostics.mismatch_fields
    @test_throws ArgumentError QWCS.cartesian_qw_operator_construction_record(
        mismatched_source,
        operators,
    )
    mismatched_carried_source = QWCS.CartesianOperatorBuildSource3D(
        chain_carried,
        build_source.basis_family,
        build_source.carried_space_kind,
        build_source.nuclei,
        build_source.nuclear_charges,
        build_source.gausslet_backend,
        build_source.requested_gausslet_backend,
        build_source.interaction_treatment,
        build_source.nuclear_term_storage,
        build_source.requested_nuclear_term_storage,
        build_source.capabilities,
        build_source.diagnostics,
        build_source.provenance,
    )
    carried_mismatch_record = QWCS.cartesian_qw_operator_construction_record(
        mismatched_carried_source,
        operators;
        throw_on_mismatch = false,
    )
    carried_mismatch_diagnostics =
        QWCS.qw_operator_construction_record_diagnostics(carried_mismatch_record)
    @test !carried_mismatch_diagnostics.source_sidecar_agree
    @test !isempty(
        intersect(
            carried_mismatch_diagnostics.mismatch_fields,
            [
                :carried_parent_axis_counts,
                :carried_parent_dimension,
                :carried_representation_final_dimension,
                :carried_provenance_input_kind,
            ],
        ),
    )
    @test_throws ArgumentError QWCS.cartesian_qw_operator_construction_record(
        mismatched_carried_source,
        operators,
    )
end

@testset "Cartesian parent gausslet basis identity" begin
    CP = GaussletBases.CartesianParentGaussletBases
    atomic_axis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 7,
            mapping = IdentityMapping(),
            reference_spacing = 1.0,
        ),
    )
    atomic_parent = CP.CartesianParentGaussletBasis3D(atomic_axis)
    construction_allocated = @allocated CP.CartesianParentGaussletBasis3D(atomic_axis)

    diatomic_basis = bond_aligned_homonuclear_qw_basis(
        bond_length = 1.4;
        core_spacing = 0.5,
        xmax_parallel = 3.0,
        xmax_transverse = 2.0,
    )
    chain_basis = bond_aligned_homonuclear_chain_qw_basis(
        natoms = 3;
        spacing = 1.2,
        core_spacing = 0.5,
        xmax_parallel = 2.5,
        xmax_transverse = 2.0,
    )
    square_basis = axis_aligned_homonuclear_square_lattice_qw_basis(
        n = 2;
        spacing = 1.2,
        core_spacing = 0.5,
        xmax_in_plane = 2.5,
        xmax_transverse = 2.0,
    )
    diatomic_parent = CP.CartesianParentGaussletBasis3D(diatomic_basis)
    chain_parent = CP.CartesianParentGaussletBasis3D(chain_basis)
    square_parent = CP.CartesianParentGaussletBasis3D(square_basis)

    function _check_parent_index_contract(parent)
        dims = CP.parent_axis_counts(parent)
        states = (
            (1, 1, 1),
            (min(2, dims[1]), min(3, dims[2]), min(4, dims[3])),
            dims,
        )
        axes = CP.parent_axes(parent)
        for state in states
            flat = CP.parent_flat_index(parent, state...)
            @test flat == GaussletBases._cartesian_flat_index(state..., dims)
            @test CP.parent_unflat_index(parent, flat) ==
                GaussletBases._cartesian_unflat_index(flat, dims)
            @test CP.parent_center(parent, state) == (
                centers(axes.x)[state[1]],
                centers(axes.y)[state[2]],
                centers(axes.z)[state[3]],
            )
        end
    end

    @test CP.cartesian_parent_gausslet_basis(atomic_parent) === atomic_parent
    @test CP.parent_axis_counts(CP.cartesian_parent_gausslet_basis(atomic_axis)) ==
        CP.parent_axis_counts(atomic_parent)
    @test CP.parent_axis_counts(CP.cartesian_parent_gausslet_basis(diatomic_basis)) ==
        CP.parent_axis_counts(diatomic_parent)
    @test CP.parent_axis_counts(CP.cartesian_parent_gausslet_basis(chain_basis)) ==
        CP.parent_axis_counts(chain_parent)
    @test CP.parent_axis_counts(CP.cartesian_parent_gausslet_basis(square_basis)) ==
        CP.parent_axis_counts(square_parent)

    @test CP.parent_axes(atomic_parent).x === atomic_axis
    @test CP.axis_basis(atomic_parent, :y) === atomic_axis
    @test CP.parent_box(atomic_parent) == (1:7, 1:7, 1:7)
    @test CP.parent_axis_counts(atomic_parent) == (7, 7, 7)
    @test CP.parent_dimension(atomic_parent) == 7^3
    @test atomic_parent.axis_sharing == :shared_xyz
    @test atomic_parent.metadata.basis_family == :mapped_uniform_same_axis
    @test construction_allocated < 50_000
    @test fieldnames(typeof(atomic_parent)) == (:axes, :parent_box, :axis_sharing, :metadata)
    @test !hasproperty(atomic_parent, :gausslet_backend)
    @test !hasproperty(atomic_parent, :backend)
    @test !hasproperty(atomic_parent.metadata, :basis_centers)
    @test !hasproperty(atomic_parent.metadata, :parent_centers)

    @test CP.parent_axes(diatomic_parent).x === diatomic_basis.basis_x
    @test CP.parent_axes(diatomic_parent).z === diatomic_basis.basis_z
    @test diatomic_parent.axis_sharing == :shared_xy
    @test diatomic_parent.metadata.basis_family == :bond_aligned_diatomic
    @test diatomic_parent.metadata.bond_axis == diatomic_basis.bond_axis
    @test diatomic_parent.metadata.nuclei == diatomic_basis.nuclei
    @test diatomic_parent.metadata.nuclear_charges == diatomic_basis.nuclear_charges
    @test CP.parent_box(diatomic_parent) == (
        1:length(diatomic_basis.basis_x),
        1:length(diatomic_basis.basis_y),
        1:length(diatomic_basis.basis_z),
    )

    @test chain_parent.axis_sharing == :shared_xy
    @test chain_parent.metadata.basis_family == :bond_aligned_homonuclear_chain
    @test chain_parent.metadata.chain_axis == chain_basis.chain_axis
    @test chain_parent.metadata.chain_coordinates == chain_basis.chain_coordinates
    @test chain_parent.metadata.nuclei == chain_basis.nuclei

    @test square_parent.axis_sharing == :shared_xy
    @test square_parent.metadata.basis_family == :axis_aligned_homonuclear_square_lattice
    @test square_parent.metadata.lattice_size == square_basis.lattice_size
    @test square_parent.metadata.x_coordinates == square_basis.x_coordinates
    @test square_parent.metadata.y_coordinates == square_basis.y_coordinates

    for parent in (atomic_parent, diatomic_parent, chain_parent, square_parent)
        _check_parent_index_contract(parent)
        @test !hasproperty(parent, :gausslet_backend)
        @test !hasproperty(parent, :backend)
    end

    @test_throws ArgumentError CP.axis_basis(atomic_parent, :q)
    @test_throws ArgumentError CP.parent_flat_index(atomic_parent, 0, 1, 1)
    @test_throws ArgumentError CP.parent_unflat_index(atomic_parent, CP.parent_dimension(atomic_parent) + 1)
    @test_throws ArgumentError CP.cartesian_parent_gausslet_basis((not_a_parent = true,))
end

@testset "Cartesian contracted parent scaffold" begin
    CP = GaussletBases.CartesianParentGaussletBases
    CCP = GaussletBases.CartesianContractedParents
    axis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 3,
            mapping = IdentityMapping(),
            reference_spacing = 1.0,
        ),
    )
    parent = CP.cartesian_parent_gausslet_basis(axis)
    parent_dim = CP.parent_dimension(parent)

    coefficients = zeros(Float64, parent_dim, 4)
    coefficients[1, 1] = 1.0
    coefficients[2, 2] = 1.0
    coefficients[2, 3] = 2.0
    coefficients[5, 4] = 1.0
    coefficients[6, 4] = -1.0
    unit_a = CCP.CartesianContractionUnit3D(
        :cube_a,
        [1, 2, 3],
        1:2;
        metadata = (shape = :cube,),
    )
    unit_b = CCP.CartesianContractionUnit3D(
        :cube_b,
        [2, 5, 6],
        3:4;
        metadata = (shape = :cube,),
    )
    contracted = CCP.CartesianContractedParent3D(
        parent,
        coefficients;
        units = [unit_a, unit_b],
        metadata = (source = :synthetic,),
    )
    audit = CCP.contracted_parent_structural_audit(contracted)

    @test CCP.contracted_parent_basis(contracted) === parent
    @test CCP.contracted_parent_coefficients(contracted) == coefficients
    @test CCP.contracted_parent_parent_dimension(contracted) == parent_dim
    @test CCP.contracted_parent_dimension(contracted) == 4
    @test CCP.contracted_parent_metadata(contracted).source == :synthetic
    @test CCP.contracted_parent_units(contracted) == [unit_a, unit_b]
    @test CCP.contracted_parent_unit_column_ranges(contracted) == [1:2, 3:4]
    @test CCP.contracted_parent_unit_support_indices(contracted) == [[1, 2, 3], [2, 5, 6]]
    @test CCP.contracted_parent_support_indices(contracted) == [1, 2, 3, 2, 5, 6]
    @test CCP.contraction_unit_role(unit_a) == :cube_a
    @test CCP.contraction_unit_support_indices(unit_a) == [1, 2, 3]
    @test CCP.contraction_unit_column_range(unit_a) == 1:2
    @test CCP.contraction_unit_metadata(unit_a).shape == :cube

    @test coefficients[:, 3] == 2.0 .* coefficients[:, 2]
    @test audit.parent_dimension == parent_dim
    @test audit.contracted_dimension == 4
    @test audit.unit_count == 2
    @test audit.support_entry_count == 6
    @test audit.unique_support_count == 5
    @test audit.duplicate_support_count == 1
    @test audit.missing_support_count == parent_dim - 5
    @test audit.outside_support_count == 0
    @test !audit.support_complete
    @test audit.column_entry_count == 4
    @test audit.unique_column_count == 4
    @test audit.duplicate_column_count == 0
    @test audit.missing_column_count == 0
    @test audit.outside_column_count == 0
    @test audit.column_ranges_cover_contract
    @test audit.structural_ok
    @test !hasproperty(contracted, :gausslet_backend)
    @test !hasproperty(contracted, :backend)
    @test !hasproperty(contracted, :overlap)
    @test !hasproperty(contracted, :interaction_matrix)

    sparse_coefficients = sparse(coefficients)
    sparse_contracted = CCP.CartesianContractedParent3D(
        parent,
        sparse_coefficients;
        units = [unit_a, unit_b],
    )
    @test CCP.contracted_parent_coefficients(sparse_contracted) isa SparseMatrixCSC{Float64,Int}
    @test CCP.contracted_parent_coefficients(sparse_contracted) == sparse_coefficients

    outside_unit = CCP.CartesianContractionUnit3D(:outside, [1, parent_dim + 1], 1:1)
    outside = CCP.CartesianContractedParent3D(
        parent,
        coefficients[:, 1:1];
        units = [outside_unit],
    )
    outside_audit = CCP.contracted_parent_structural_audit(outside)
    @test outside_audit.outside_support_count == 1
    @test !outside_audit.support_complete
    @test !outside_audit.structural_ok

    overlapping_columns = CCP.CartesianContractedParent3D(
        parent,
        coefficients;
        units = [
            CCP.CartesianContractionUnit3D(:left, [1], 1:2),
            CCP.CartesianContractionUnit3D(:right, [2], 2:4),
        ],
    )
    overlapping_column_audit = CCP.contracted_parent_structural_audit(overlapping_columns)
    @test overlapping_column_audit.duplicate_column_count == 1
    @test !overlapping_column_audit.column_ranges_cover_contract
    @test !overlapping_column_audit.structural_ok

    missing_columns = CCP.CartesianContractedParent3D(
        parent,
        coefficients;
        units = [
            CCP.CartesianContractionUnit3D(:first, [1], 1:1),
            CCP.CartesianContractionUnit3D(:last, [2], 3:4),
        ],
    )
    missing_column_audit = CCP.contracted_parent_structural_audit(missing_columns)
    @test missing_column_audit.missing_column_count == 1
    @test !missing_column_audit.column_ranges_cover_contract
    @test !missing_column_audit.structural_ok

    @test_throws ArgumentError CCP.CartesianContractionUnit3D(:empty, [1], 1:0)
    @test_throws ArgumentError CCP.CartesianContractedParent3D(
        parent,
        coefficients;
        units = [CCP.CartesianContractionUnit3D(:bad_columns, [1], 4:5)],
    )
    @test_throws DimensionMismatch CCP.CartesianContractedParent3D(
        parent,
        coefficients[1:(end - 1), :],
    )
end

@testset "Cartesian contracted parent metric packet" begin
    CP = GaussletBases.CartesianParentGaussletBases
    CCP = GaussletBases.CartesianContractedParents
    CCPM = GaussletBases.CartesianContractedParentMetrics
    axis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 3,
            mapping = IdentityMapping(),
            reference_spacing = 1.0,
        ),
    )
    parent = CP.cartesian_parent_gausslet_basis(axis)
    parent_dim = CP.parent_dimension(parent)
    coefficients = zeros(Float64, parent_dim, 3)
    coefficients[1, 1] = 1.0
    coefficients[2, 1] = 0.25
    coefficients[5, 2] = 1.0
    coefficients[14, 2] = -0.5
    coefficients[parent_dim, 3] = 1.0
    units = [
        CCP.CartesianContractionUnit3D(:left, [1, 2], 1:1),
        CCP.CartesianContractionUnit3D(:middle, [5, 14], 2:2),
        CCP.CartesianContractionUnit3D(:right, [parent_dim], 3:3),
    ]
    contracted = CCP.CartesianContractedParent3D(parent, coefficients; units)
    packet = CCPM.cartesian_contracted_parent_metric_packet(contracted)
    reference = CCPM.cartesian_contracted_parent_metric_packet_dense_reference(contracted)

    @test packet isa CCPM.CartesianContractedParentMetricPacket3D
    @test CCPM.contracted_parent_metric_packet_parent(packet) === contracted
    @test packet.diagnostics.construction_path == :support_local_product
    @test packet.diagnostics.dense_parent_matrix_used == false
    @test reference.diagnostics.construction_path == :dense_reference_oracle
    @test reference.diagnostics.dense_parent_matrix_used == true
    @test size(packet.overlap) == (3, 3)
    @test size(packet.centers) == (3, 3)
    @test length(packet.weights) == 3
    @test packet.overlap ≈ reference.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test packet.position_x ≈ reference.position_x atol = 1.0e-12 rtol = 1.0e-12
    @test packet.position_y ≈ reference.position_y atol = 1.0e-12 rtol = 1.0e-12
    @test packet.position_z ≈ reference.position_z atol = 1.0e-12 rtol = 1.0e-12
    @test packet.weights ≈ reference.weights atol = 1.0e-12 rtol = 1.0e-12
    @test packet.first_moments ≈ reference.first_moments atol = 1.0e-12 rtol = 1.0e-12
    @test packet.centers ≈ reference.centers atol = 1.0e-12 rtol = 1.0e-12
    @test isfinite(packet.diagnostics.overlap_symmetry_error)
    @test isfinite(packet.diagnostics.overlap_identity_error)

    sparse_coefficients = sparse([1, 7, 13, 27], [1, 1, 2, 2], [0.5, -0.25, 1.0, 0.125], parent_dim, 2)
    sparse_contracted = CCP.CartesianContractedParent3D(parent, sparse_coefficients)
    sparse_packet = CCPM.cartesian_contracted_parent_metric_packet(sparse_contracted)
    sparse_reference = CCPM.cartesian_contracted_parent_metric_packet_dense_reference(sparse_contracted)
    @test CCP.contracted_parent_coefficients(sparse_contracted) isa SparseMatrixCSC{Float64,Int}
    @test sparse_packet.diagnostics.dense_parent_matrix_used == false
    @test sparse_packet.diagnostics.coefficient_storage == :SparseMatrixCSC
    @test sparse_packet.overlap ≈ sparse_reference.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test sparse_packet.weights ≈ sparse_reference.weights atol = 1.0e-12 rtol = 1.0e-12

    larger_axis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 5,
            mapping = IdentityMapping(),
            reference_spacing = 1.0,
        ),
    )
    larger_parent = CP.cartesian_parent_gausslet_basis(larger_axis)
    larger_coefficients = sparse([1, 125], [1, 2], [1.0, 1.0], 125, 2)
    larger_contracted = CCP.CartesianContractedParent3D(larger_parent, larger_coefficients)
    larger_packet = CCPM.cartesian_contracted_parent_metric_packet(larger_contracted)
    @test larger_packet.diagnostics.dense_parent_matrix_used == false
    @test_throws ArgumentError CCPM.cartesian_contracted_parent_metric_packet_dense_reference(
        larger_contracted;
        max_parent_dimension = 64,
    )

    distorted_axis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 3,
            mapping = AsinhMapping(a = 0.25, s = 0.5, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    distorted_parent = CP.cartesian_parent_gausslet_basis(distorted_axis)
    distorted_coefficients = sparse([1, 27], [1, 2], [1.0, 1.0], 27, 2)
    distorted_contracted = CCP.CartesianContractedParent3D(
        distorted_parent,
        distorted_coefficients,
    )
    explicit_axis_metric = (
        overlap = Matrix{Float64}(I, 3, 3),
        position = Matrix(Diagonal(centers(distorted_axis))),
        weights = ones(Float64, 3),
        centers = Float64.(centers(distorted_axis)),
        source = :test_explicit_no_quadrature,
    )
    explicit_packet = CCPM.cartesian_contracted_parent_metric_packet(
        distorted_contracted;
        axis_metrics = (
            x = explicit_axis_metric,
            y = explicit_axis_metric,
            z = explicit_axis_metric,
        ),
    )
    @test_throws ArgumentError CCPM.cartesian_contracted_parent_metric_packet(distorted_contracted)
    @test explicit_packet.diagnostics.dense_parent_matrix_used == false
    @test explicit_packet.diagnostics.axis_metric_sources.x == :test_explicit_no_quadrature
    @test_throws ArgumentError CCPM.cartesian_contracted_parent_metric_packet(
        sparse_contracted;
        construction_path = :product_staged_metric_contraction,
    )
end

@testset "Cartesian basis representation for nested fixed blocks" begin
    CP = GaussletBases.CartesianParentGaussletBases
    CCP = GaussletBases.CartesianContractedParents
    CCS = GaussletBases.CartesianCarriedSpaces
    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 13,
            mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    fixed_block = one_center_atomic_full_parent_fixed_block(
        basis;
        expansion,
        nside = 5,
    )
    fixed_parent = CP.cartesian_parent_gausslet_basis(fixed_block)
    fixed_contracted_parent = CCP.cartesian_contracted_parent(fixed_block)
    fixed_contracted_audit = CCP.contracted_parent_structural_audit(fixed_contracted_parent)
    representation = basis_representation(fixed_block)
    metadata = basis_metadata(representation)

    @test representation isa CartesianBasisRepresentation3D
    @test metadata.basis_kind == :nested_fixed_block
    @test metadata.parent_kind == :cartesian_product_basis
    @test metadata.parent_axis_counts == (13, 13, 13)
    @test metadata.parent_axis_counts == CP.parent_axis_counts(fixed_parent)
    @test metadata.parent_dimension == 13^3
    @test metadata.parent_dimension == CP.parent_dimension(fixed_parent)
    @test metadata.final_dimension == size(fixed_block.coefficient_matrix, 2)
    @test metadata.working_box == (1:13, 1:13, 1:13)
    @test metadata.route_metadata.shell_kind == :shell_sequence
    @test metadata.route_metadata.working_box_profile == :full_parent
    @test metadata.route_metadata.nside == 5
    @test metadata.route_metadata.support_count == length(fixed_block.support_indices)
    @test size(representation.coefficient_matrix) == size(fixed_block.coefficient_matrix)
    @test representation.support_indices == fixed_block.support_indices
    @test length(representation.support_states) == length(fixed_block.support_indices)
    @test size(metadata.basis_centers) == size(fixed_block.fixed_centers)
    @test CP.parent_center(fixed_parent, (1, 1, 1)) == (
        centers(basis)[1],
        centers(basis)[1],
        centers(basis)[1],
    )
    @test !hasproperty(fixed_parent, :gausslet_backend)
    @test !hasproperty(fixed_parent, :backend)
    @test CCP.contracted_parent_basis(fixed_contracted_parent).parent_box ==
        fixed_parent.parent_box
    @test CCP.contracted_parent_parent_dimension(fixed_contracted_parent) ==
        metadata.parent_dimension
    @test CCP.contracted_parent_dimension(fixed_contracted_parent) ==
        size(fixed_block.coefficient_matrix, 2)
    @test CCP.contracted_parent_coefficients(fixed_contracted_parent) isa
        SparseMatrixCSC{Float64,Int}
    @test CCP.contracted_parent_coefficients(fixed_contracted_parent) ==
        Matrix{Float64}(fixed_block.coefficient_matrix)
    @test only(CCP.contracted_parent_units(fixed_contracted_parent)).role ==
        :nested_fixed_block
    @test only(CCP.contracted_parent_unit_column_ranges(fixed_contracted_parent)) ==
        1:size(fixed_block.coefficient_matrix, 2)
    @test fixed_contracted_audit.outside_support_count == 0
    @test fixed_contracted_audit.column_ranges_cover_contract
    @test fixed_contracted_audit.structural_ok
    @test !hasproperty(fixed_contracted_parent, :gausslet_backend)
    @test !hasproperty(fixed_contracted_parent, :backend)
    @test !hasproperty(fixed_contracted_parent, :interaction_matrix)
    @test !isnothing(fixed_block.factorized_cartesian_parent_basis[])
    @test hasproperty(representation.parent_data, :factorized_cartesian_parent_basis)
    @test representation.parent_data.factorized_cartesian_parent_basis ===
          fixed_block.factorized_cartesian_parent_basis[]
    fixed_carried = CCS.cartesian_carried_space(fixed_block)
    @test CCS.carried_space_parent(fixed_carried) isa CP.CartesianParentGaussletBasis3D
    @test CCS.carried_space_contracted_parent(fixed_carried) isa CCP.CartesianContractedParent3D
    @test CCS.carried_space_representation(fixed_carried) isa CartesianBasisRepresentation3D
    @test CCS.carried_space_diagnostics(fixed_carried).parent_dimension ==
        metadata.parent_dimension
    @test CCS.carried_space_diagnostics(fixed_carried).contracted_dimension ==
        metadata.final_dimension
    @test CCS.carried_space_diagnostics(fixed_carried).contracted_dimension_matches_representation
    @test CCS.carried_space_diagnostics(fixed_carried).contracted_parent_dimension_matches_parent
    @test CCS.carried_space_provenance(fixed_carried).input_kind == :nested_fixed_block

    square_basis, _source, square_fixed_block, _diagnostics =
        _axis_aligned_homonuclear_square_lattice_nested_fixture()
    square_representation = basis_representation(square_fixed_block)
    square_metadata = basis_metadata(square_representation)
    @test square_metadata.basis_kind == :nested_fixed_block
    @test square_metadata.parent_axis_counts == (
        length(square_basis.basis_x),
        length(square_basis.basis_y),
        length(square_basis.basis_z),
    )
    @test square_metadata.parent_dimension == prod(square_metadata.parent_axis_counts)
    @test square_metadata.final_dimension == size(square_fixed_block.coefficient_matrix, 2)
    @test size(square_representation.coefficient_matrix) == size(square_fixed_block.coefficient_matrix)
    @test square_metadata.working_box == square_fixed_block.shell.working_box
    @test square_metadata.route_metadata.support_count == length(square_fixed_block.support_indices)
    @test !isnothing(square_fixed_block.factorized_cartesian_parent_basis[])
    @test hasproperty(square_representation.parent_data, :factorized_cartesian_parent_basis)
    @test square_representation.parent_data.factorized_cartesian_parent_basis ===
          square_fixed_block.factorized_cartesian_parent_basis[]
    square_carried = CCS.cartesian_carried_space(square_fixed_block)
    @test CCS.carried_space_parent(square_carried) isa CP.CartesianParentGaussletBasis3D
    @test CCS.carried_space_contracted_parent(square_carried) isa CCP.CartesianContractedParent3D
    @test CCS.carried_space_diagnostics(square_carried).parent_axis_counts ==
        square_metadata.parent_axis_counts
    @test CCS.carried_space_diagnostics(square_carried).contracted_dimension ==
        square_metadata.final_dimension
    @test CCS.carried_space_diagnostics(square_carried).contracted_dimension_matches_representation
    @test CCS.carried_space_provenance(square_carried).input_kind == :nested_fixed_block
end

function _with_sparse_nested_coefficients(fixed_block::GaussletBases._NestedFixedBlock3D)
    return GaussletBases._NestedFixedBlock3D(
        fixed_block.parent_basis,
        fixed_block.shell,
        fixed_block.gausslet_backend,
        sparse(fixed_block.coefficient_matrix),
        fixed_block.support_indices,
        fixed_block.overlap,
        fixed_block.kinetic,
        fixed_block.position_x,
        fixed_block.position_y,
        fixed_block.position_z,
        fixed_block.x2_x,
        fixed_block.x2_y,
        fixed_block.x2_z,
        fixed_block.weights,
        fixed_block.gaussian_sum,
        fixed_block.pair_sum,
        fixed_block.fixed_centers,
        GaussletBases._nested_factorized_basis_cache(
            fixed_block.factorized_cartesian_parent_basis[],
        ),
        GaussletBases._nested_staged_by_center_sidecar_cache(
            fixed_block.staged_by_center_sidecar[],
        ),
    )
end

@testset "Nested coefficient maps support sparse storage" begin
    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 13,
            mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    direct_fixed_block = one_center_atomic_full_parent_fixed_block(
        basis;
        expansion,
        nside = 5,
    )
    sparse_fixed_block = _with_sparse_nested_coefficients(direct_fixed_block)

    direct_representation = basis_representation(direct_fixed_block)
    sparse_representation = basis_representation(sparse_fixed_block)
    support_coefficients = GaussletBases._nested_support_coefficient_slice(
        direct_fixed_block.shell.coefficient_matrix,
        direct_fixed_block.shell.support_indices,
    )
    support_workspace, contraction_scratch = GaussletBases._nested_support_reference_workspaces(
        support_coefficients,
        length(direct_fixed_block.shell.support_indices),
        size(direct_fixed_block.shell.coefficient_matrix, 2),
    )

    @test direct_fixed_block.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test direct_representation.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test sparse_fixed_block.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test sparse_representation.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
    @test support_coefficients isa SparseMatrixCSC{Float64,Int}
    @test size(support_workspace) == (0, 0)
    @test size(contraction_scratch) == (0, 0)
    @test size(support_coefficients) == (
        length(direct_fixed_block.shell.support_indices),
        size(direct_fixed_block.shell.coefficient_matrix, 2),
    )
    @test Matrix(sparse_representation.coefficient_matrix) ≈
        Matrix(direct_representation.coefficient_matrix) atol = 1.0e-12 rtol = 1.0e-12
    @test cross_overlap(sparse_representation, sparse_representation) ≈
        cross_overlap(direct_representation, direct_representation) atol = 1.0e-10 rtol = 1.0e-10

    mktemp() do sparse_path, sparse_io
        close(sparse_io)
        sparse_matrix = sparse(direct_fixed_block.coefficient_matrix)
        jldopen(sparse_path, "w") do file
            file["matrix"] = sparse_matrix
        end
        restored = jldopen(sparse_path, "r") do file
            file["matrix"]
        end
        @test restored isa SparseMatrixCSC{Float64,Int}
        @test restored == sparse_matrix
    end
end

function _atomic_hybrid_cartesian_representation_fixture()
    return _cached_fixture(:atomic_hybrid_cartesian_representation_fixture, () -> begin
        basis = build_basis(
            MappedUniformBasisSpec(
                :G10;
                count = 13,
                mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
                reference_spacing = 1.0,
            ),
        )
        expansion = coulomb_gaussian_expansion(doacc = false)
        fixed_full = one_center_atomic_full_parent_fixed_block(
            basis;
            expansion,
            nside = 5,
        )
        fixed_legacy = one_center_atomic_legacy_profile_fixed_block(
            basis;
            expansion,
            working_box = 2:12,
            nside = 5,
        )
        supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 1)
        full_ops = ordinary_cartesian_qiu_white_operators(
            fixed_full,
            supplement;
            expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
            residual_keep_policy = :near_null_only,
        )
        legacy_ops = ordinary_cartesian_qiu_white_operators(
            fixed_legacy,
            supplement;
            expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
            residual_keep_policy = :near_null_only,
        )
        (
            basis = basis,
            expansion = expansion,
            fixed_full = fixed_full,
            fixed_legacy = fixed_legacy,
            fixed_full_rep = basis_representation(fixed_full),
            fixed_legacy_rep = basis_representation(fixed_legacy),
            supplement = supplement,
            full_ops = full_ops,
            legacy_ops = legacy_ops,
            full_rep = basis_representation(full_ops),
            legacy_rep = basis_representation(legacy_ops),
        )
    end)
end

function _metric_normalize_orbital(
    coefficients::AbstractVector,
    overlap::AbstractMatrix{<:Real},
)
    orbital = Float64[Float64(real(value)) for value in coefficients]
    norm2 = Float64(real(dot(orbital, overlap * orbital)))
    norm2 > 0.0 || throw(ArgumentError("orbital must have nonzero target-metric norm"))
    return orbital ./ sqrt(norm2)
end

function _metric_orbital_overlap(
    left::AbstractVector,
    right::AbstractVector,
    overlap::AbstractMatrix{<:Real},
)
    normalized_left = _metric_normalize_orbital(left, overlap)
    normalized_right = _metric_normalize_orbital(right, overlap)
    return Float64(real(dot(normalized_left, overlap * normalized_right)))
end

function _ordinary_cartesian_hybrid_orbital_observables(
    operators::OrdinaryCartesianOperators3D,
    orbital::AbstractVector;
    overlap_tol::Real = 1.0e-7,
)
    overlap = Matrix{Float64}(operators.overlap)
    normalized = _metric_normalize_orbital(orbital, overlap)
    one_body = Float64(real(dot(normalized, operators.one_body_hamiltonian * normalized)))
    vee = GaussletBases.ordinary_cartesian_vee_expectation(
        operators,
        normalized;
        overlap_tol = overlap_tol,
    )
    return (
        orbital = normalized,
        metric_norm_error = abs(Float64(real(dot(normalized, overlap * normalized))) - 1.0),
        one_body = one_body,
        vee = vee,
        total = 2.0 * one_body + vee,
    )
end

function _atomic_direct_product_he_extent_change_contract_fixture(;
    source_count::Int = 3,
    target_count::Int = 5,
)
    key = Symbol(:atomic_direct_product_he_extent_change_contract, source_count, target_count)
    return _cached_fixture(key, () -> begin
        mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0)
        source_basis = build_basis(
            MappedUniformBasisSpec(
                :G10;
                count = source_count,
                mapping,
                reference_spacing = 1.0,
            ),
        )
        target_basis = build_basis(
            MappedUniformBasisSpec(
                :G10;
                count = target_count,
                mapping,
                reference_spacing = 1.0,
            ),
        )

        source_rep = basis_representation(source_basis)
        target_rep = basis_representation(target_basis)
        offset = (target_count - source_count) ÷ 2
        shared_slice = (offset + 1):(offset + source_count)

        return (
            source_count = source_count,
            target_count = target_count,
            shared_slice = shared_slice,
            source_rep = source_rep,
            target_rep = target_rep,
            centers_subset =
                source_rep.metadata.center_data == target_rep.metadata.center_data[shared_slice],
            weights_subset =
                source_rep.metadata.integral_weight_data ==
                target_rep.metadata.integral_weight_data[shared_slice],
            coefficient_core_match =
                source_rep.coefficient_matrix ==
                target_rep.coefficient_matrix[shared_slice, shared_slice],
        )
    end)
end

function _atomic_hybrid_he_same_parent_stress_fixture(;
    parent_count::Int = 7,
    source_working_box::UnitRange{Int} = 2:6,
    supplement_lmax::Int = 1,
)
    key = Symbol(
        :atomic_hybrid_he_same_parent_stress_fixture,
        parent_count,
        first(source_working_box),
        last(source_working_box),
        supplement_lmax,
    )
    return _cached_fixture(key, () -> begin
        expansion = coulomb_gaussian_expansion(doacc = false)
        mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0)
        supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = supplement_lmax)

        parent_basis = build_basis(
            MappedUniformBasisSpec(
                :G10;
                count = parent_count,
                mapping,
                reference_spacing = 1.0,
            ),
        )

        source_fixed = one_center_atomic_legacy_profile_fixed_block(
            parent_basis;
            expansion,
            working_box = source_working_box,
            nside = 5,
        )
        target_fixed = one_center_atomic_full_parent_fixed_block(
            parent_basis;
            expansion,
            nside = 5,
        )

        source_ops = ordinary_cartesian_qiu_white_operators(
            source_fixed,
            supplement;
            expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
            residual_keep_policy = :near_null_only,
        )
        target_ops = ordinary_cartesian_qiu_white_operators(
            target_fixed,
            supplement;
            expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
            residual_keep_policy = :near_null_only,
        )

        source_rep = basis_representation(source_ops)
        target_rep = basis_representation(target_ops)
        source_check = GaussletBases.ordinary_cartesian_1s2_check(
            source_ops;
            overlap_tol = 1.0e-7,
        )
        target_check = GaussletBases.ordinary_cartesian_1s2_check(
            target_ops;
            overlap_tol = 1.0e-7,
        )
        source_observables = _ordinary_cartesian_hybrid_orbital_observables(
            source_ops,
            source_check.orbital;
            overlap_tol = 1.0e-7,
        )
        target_observables = _ordinary_cartesian_hybrid_orbital_observables(
            target_ops,
            target_check.orbital;
            overlap_tol = 1.0e-7,
        )

        transfer = transfer_orbitals(source_observables.orbital, source_rep, target_rep)
        transferred_observables = _ordinary_cartesian_hybrid_orbital_observables(
            target_ops,
            transfer.coefficients;
            overlap_tol = 1.0e-7,
        )
        target_overlap = Matrix{Float64}(target_ops.overlap)
        overlap_with_target = _metric_orbital_overlap(
            transferred_observables.orbital,
            target_observables.orbital,
            target_overlap,
        )
        sign = overlap_with_target < 0.0 ? -1.0 : 1.0
        aligned_transferred_observables = _ordinary_cartesian_hybrid_orbital_observables(
            target_ops,
            sign .* transferred_observables.orbital;
            overlap_tol = 1.0e-7,
        )
        aligned_overlap_to_target = abs(
            _metric_orbital_overlap(
                aligned_transferred_observables.orbital,
                target_observables.orbital,
                target_overlap,
            ),
        )

        return (
            parent_basis = parent_basis,
            source_fixed = source_fixed,
            target_fixed = target_fixed,
            supplement = supplement,
            source_working_box = source_working_box,
            target_working_box = target_fixed.shell.working_box,
            source_ops = source_ops,
            target_ops = target_ops,
            source_rep = source_rep,
            target_rep = target_rep,
            source_check = source_check,
            target_check = target_check,
            source_observables = source_observables,
            target_observables = target_observables,
            transfer = transfer,
            transferred_observables = transferred_observables,
            aligned_transferred_observables = aligned_transferred_observables,
            aligned_overlap_to_target = aligned_overlap_to_target,
        )
    end)
end

@testset "Cartesian basis representation for atomic QW residual bases" begin
    fixture = _atomic_hybrid_cartesian_representation_fixture()
    operators = fixture.full_ops
    representation = fixture.full_rep
    metadata = basis_metadata(representation)
    supplement_representation = representation.parent_data.supplement_representation

    @test representation isa CartesianBasisRepresentation3D
    @test metadata.basis_kind == :hybrid_residual
    @test metadata.parent_kind == :cartesian_plus_supplement_raw
    @test metadata.final_dimension == length(operators.orbital_data)
    @test metadata.final_dimension == size(operators.raw_to_final, 2)
    @test metadata.parent_dimension == size(operators.raw_to_final, 1)
    @test metadata.route_metadata.gausslet_count == operators.gausslet_count
    @test metadata.route_metadata.residual_count == operators.residual_count
    @test metadata.route_metadata.supplement_kind == :atomic_cartesian_shell
    @test metadata.route_metadata.supplement_lmax == fixture.supplement.lmax
    @test size(representation.coefficient_matrix) == size(operators.raw_to_final)
    @test length(representation.parent_labels) == size(operators.raw_to_final, 1)
    @test size(representation.parent_centers, 1) == size(operators.raw_to_final, 1)
    @test hasproperty(representation.parent_data, :cartesian_parent_representation)
    @test representation.parent_data.cartesian_parent_representation.metadata.basis_kind ==
        :nested_fixed_block
    @test representation.parent_data.cartesian_parent_representation.metadata.final_dimension ==
        operators.gausslet_count
    @test hasproperty(representation.parent_data, :supplement_representation)
    @test hasproperty(representation.parent_data, :factorized_cartesian_parent_basis)
    @test hasproperty(representation.parent_data, :cartesian_supplement_axis_tables)
    @test supplement_representation isa CartesianGaussianShellSupplementRepresentation3D
    @test supplement_representation.supplement_kind == :atomic_cartesian_shell
    @test length(supplement_representation.orbitals) ==
        size(operators.raw_to_final, 1) - operators.gausslet_count
    @test size(representation.parent_data.cartesian_supplement_axis_tables.x, 2) ==
        length(supplement_representation.orbitals)
    @test size(representation.parent_data.cartesian_supplement_axis_tables.y, 2) ==
        length(supplement_representation.orbitals)
    @test size(representation.parent_data.cartesian_supplement_axis_tables.z, 2) ==
        length(supplement_representation.orbitals)
    @test any(
        orbital -> sum(orbital.angular_powers) > 0,
        supplement_representation.orbitals,
    )
end

@testset "Cartesian basis representation cross overlap" begin
    diatomic_basis14, diatomic_ops14, _check14 = _bond_aligned_diatomic_qw_fixture(; bond_length = 1.4)
    diatomic_basis20, _diatomic_ops20, _check20 = _bond_aligned_diatomic_qw_fixture(; bond_length = 2.0)
    diatomic_rep14 = basis_representation(diatomic_basis14)
    diatomic_rep20 = basis_representation(diatomic_basis20)
    direct_self = cross_overlap(diatomic_rep14, diatomic_rep14)
    direct_cross = cross_overlap(diatomic_rep14, diatomic_rep20)
    direct_cross_reverse = cross_overlap(diatomic_rep20, diatomic_rep14)

    @test size(direct_self) == size(diatomic_ops14.overlap)
    @test direct_self ≈ diatomic_ops14.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test size(direct_cross) == (
        diatomic_rep14.metadata.final_dimension,
        diatomic_rep20.metadata.final_dimension,
    )
    @test direct_cross ≈ transpose(direct_cross_reverse) atol = 1.0e-10 rtol = 1.0e-10

    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 13,
            mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    fixed_full = one_center_atomic_full_parent_fixed_block(
        basis;
        expansion = expansion,
        nside = 5,
    )
    fixed_legacy = one_center_atomic_legacy_profile_fixed_block(
        basis;
        expansion = expansion,
        working_box = 2:12,
        nside = 5,
    )
    fixed_full_rep = basis_representation(fixed_full)
    fixed_legacy_rep = basis_representation(fixed_legacy)
    S1d = basis_representation(basis).basis_matrices.overlap
    Sparent = kron(S1d, kron(S1d, S1d))

    fixed_self = cross_overlap(fixed_full_rep, fixed_full_rep)
    fixed_cross = cross_overlap(fixed_full_rep, fixed_legacy_rep)
    fixed_cross_reverse = cross_overlap(fixed_legacy_rep, fixed_full_rep)
    fixed_cross_expected =
        transpose(fixed_full.coefficient_matrix) * Sparent * fixed_legacy.coefficient_matrix

    @test size(fixed_self) == size(fixed_full.overlap)
    @test fixed_self ≈ fixed_full.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test size(fixed_cross) == (
        size(fixed_full.coefficient_matrix, 2),
        size(fixed_legacy.coefficient_matrix, 2),
    )
    @test fixed_cross ≈ fixed_cross_expected atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_cross ≈ transpose(fixed_cross_reverse) atol = 1.0e-10 rtol = 1.0e-10

    square_basis, _square_source, square_fixed_block, _square_diagnostics =
        _axis_aligned_homonuclear_square_lattice_nested_fixture()
    square_basis_rep = basis_representation(square_basis)
    square_fixed_rep = basis_representation(square_fixed_block)
    square_parent_x = basis_representation(square_basis.basis_x).basis_matrices.overlap
    square_parent_y = basis_representation(square_basis.basis_y).basis_matrices.overlap
    square_parent_z = basis_representation(square_basis.basis_z).basis_matrices.overlap
    square_parent_overlap = kron(square_parent_x, kron(square_parent_y, square_parent_z))
    square_cross = cross_overlap(square_basis_rep, square_fixed_rep)
    square_cross_expected = square_parent_overlap * square_fixed_block.coefficient_matrix

    @test size(square_cross) == (
        square_basis_rep.metadata.final_dimension,
        square_fixed_rep.metadata.final_dimension,
    )
    @test square_cross ≈ square_cross_expected atol = 1.0e-10 rtol = 1.0e-10

    hybrid_fixture = _atomic_hybrid_cartesian_representation_fixture()
    hybrid_self = cross_overlap(hybrid_fixture.full_rep, hybrid_fixture.full_rep)
    hybrid_cross = cross_overlap(hybrid_fixture.full_rep, hybrid_fixture.legacy_rep)
    hybrid_cross_reverse = cross_overlap(hybrid_fixture.legacy_rep, hybrid_fixture.full_rep)
    hybrid_parent_cross = cross_overlap(hybrid_fixture.fixed_full_rep, hybrid_fixture.full_rep)
    hybrid_parent_cross_reverse = cross_overlap(hybrid_fixture.full_rep, hybrid_fixture.fixed_full_rep)
    hybrid_self_dense = GaussletBases._cartesian_mixed_raw_cross_overlap_dense_reference(
        hybrid_fixture.full_rep,
        hybrid_fixture.full_rep,
    )

    @test size(hybrid_self) == size(hybrid_fixture.full_ops.overlap)
    @test hybrid_self ≈ hybrid_self_dense atol = 1.0e-10 rtol = 1.0e-10
    @test size(hybrid_cross) == (
        hybrid_fixture.full_rep.metadata.final_dimension,
        hybrid_fixture.legacy_rep.metadata.final_dimension,
    )
    @test hybrid_cross ≈ transpose(hybrid_cross_reverse) atol = 1.0e-10 rtol = 1.0e-10
    @test size(hybrid_parent_cross) == (
        hybrid_fixture.fixed_full_rep.metadata.final_dimension,
        hybrid_fixture.full_rep.metadata.final_dimension,
    )
    @test hybrid_parent_cross ≈ transpose(hybrid_parent_cross_reverse) atol = 1.0e-10 rtol = 1.0e-10
end

@testset "Cartesian basis projector and orbital transfer" begin
    square_basis, _square_source, square_fixed_block, _square_diagnostics =
        _axis_aligned_homonuclear_square_lattice_nested_fixture()
    square_basis_rep = basis_representation(square_basis)
    square_fixed_rep = basis_representation(square_fixed_block)

    direct_self_projector = basis_projector(square_basis_rep, square_basis_rep)
    direct_self_coefficients =
        reshape(sin.(Float64.(1:(2 * square_basis_rep.metadata.final_dimension))), :, 2)
    direct_self_transfer =
        transfer_orbitals(direct_self_coefficients, square_basis_rep, square_basis_rep)

    @test direct_self_projector.matrix ≈ Matrix{Float64}(
        LinearAlgebra.I,
        square_basis_rep.metadata.final_dimension,
        square_basis_rep.metadata.final_dimension,
    ) atol = 1.0e-10 rtol = 1.0e-10
    @test direct_self_transfer.coefficients ≈ direct_self_coefficients atol = 1.0e-10 rtol = 1.0e-10
    @test direct_self_transfer.diagnostics.transfer_path == :same_parent_cross_overlap_transfer
    @test direct_self_transfer.diagnostics.projector_residual_inf < 1.0e-10
    @test direct_self_transfer.diagnostics.transferred_residual_inf < 1.0e-10
    @test direct_self_transfer.diagnostics.source_metric_trace ≈ direct_self_transfer.diagnostics.target_metric_trace atol =
          1.0e-10 rtol = 1.0e-10

    nested_coefficients = cos.(Float64.(1:square_fixed_rep.metadata.final_dimension))
    nested_to_direct = transfer_orbitals(nested_coefficients, square_fixed_rep, square_basis_rep)
    nested_embedded = square_fixed_block.coefficient_matrix * nested_coefficients
    direct_back_to_nested =
        transfer_orbitals(nested_to_direct.coefficients, square_basis_rep, square_fixed_rep)

    @test nested_to_direct.coefficients ≈ nested_embedded atol = 1.0e-10 rtol = 1.0e-10
    @test direct_back_to_nested.coefficients ≈ nested_coefficients atol = 1.0e-10 rtol = 1.0e-10
    @test direct_back_to_nested.diagnostics.transferred_residual_inf < 1.0e-10

    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 13,
            mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    fixed_full = one_center_atomic_full_parent_fixed_block(
        basis;
        expansion = expansion,
        nside = 5,
    )
    fixed_legacy = one_center_atomic_legacy_profile_fixed_block(
        basis;
        expansion = expansion,
        working_box = 2:12,
        nside = 5,
    )
    fixed_full_rep = basis_representation(fixed_full)
    fixed_legacy_rep = basis_representation(fixed_legacy)

    full_to_legacy = basis_projector(fixed_full_rep, fixed_legacy_rep)
    legacy_to_full = basis_projector(fixed_legacy_rep, fixed_full_rep)
    full_coefficients =
        reshape(cos.(Float64.(1:(2 * fixed_full_rep.metadata.final_dimension))), :, 2)
    legacy_coefficients =
        reshape(sin.(Float64.(1:(2 * fixed_legacy_rep.metadata.final_dimension))), :, 2)
    transferred_full_to_legacy =
        transfer_orbitals(full_coefficients, fixed_full_rep, fixed_legacy_rep)
    transferred_legacy_to_full =
        transfer_orbitals(legacy_coefficients, fixed_legacy_rep, fixed_full_rep)

    @test full_to_legacy.matrix ≈ cross_overlap(fixed_legacy_rep, fixed_full_rep) atol =
          1.0e-10 rtol = 1.0e-10
    @test legacy_to_full.matrix ≈ cross_overlap(fixed_full_rep, fixed_legacy_rep) atol =
          1.0e-10 rtol = 1.0e-10
    @test transferred_full_to_legacy.diagnostics.transferred_residual_inf < 1.0e-10
    @test transferred_legacy_to_full.diagnostics.transferred_residual_inf < 1.0e-10

    hybrid_fixture = _atomic_hybrid_cartesian_representation_fixture()
    hybrid_self_projector = basis_projector(hybrid_fixture.full_rep, hybrid_fixture.full_rep)
    hybrid_self_dense = GaussletBases._cartesian_mixed_raw_cross_overlap_dense_reference(
        hybrid_fixture.full_rep,
        hybrid_fixture.full_rep,
    )
    hybrid_cross_dense = GaussletBases._cartesian_mixed_raw_cross_overlap_dense_reference(
        hybrid_fixture.legacy_rep,
        hybrid_fixture.full_rep,
    )
    hybrid_self_coefficients = reshape(
        sin.(Float64.(1:(2 * hybrid_fixture.full_rep.metadata.final_dimension))),
        :,
        2,
    )
    hybrid_self_transfer = transfer_orbitals(
        hybrid_self_coefficients,
        hybrid_fixture.full_rep,
        hybrid_fixture.full_rep,
    )
    hybrid_full_to_legacy = basis_projector(hybrid_fixture.full_rep, hybrid_fixture.legacy_rep)
    hybrid_parent_to_full =
        basis_projector(hybrid_fixture.fixed_full_rep, hybrid_fixture.full_rep)
    hybrid_full_to_parent =
        basis_projector(hybrid_fixture.full_rep, hybrid_fixture.fixed_full_rep)
    hybrid_transfer_from_projector =
        transfer_orbitals(hybrid_self_coefficients, hybrid_self_projector)

    @test hybrid_self_projector.matrix ≈ Matrix{Float64}(
        LinearAlgebra.I,
        hybrid_fixture.full_rep.metadata.final_dimension,
        hybrid_fixture.full_rep.metadata.final_dimension,
    ) atol = 1.0e-10 rtol = 1.0e-10
    @test cross_overlap(hybrid_fixture.full_rep, hybrid_fixture.full_rep) ≈
          hybrid_self_dense atol = 1.0e-10 rtol = 1.0e-10
    @test cross_overlap(hybrid_fixture.legacy_rep, hybrid_fixture.full_rep) ≈
          hybrid_cross_dense atol = 1.0e-10 rtol = 1.0e-10
    @test hybrid_full_to_legacy.matrix ≈ hybrid_cross_dense atol = 1.0e-10 rtol = 1.0e-10
    @test hybrid_self_transfer.coefficients ≈ hybrid_self_coefficients atol = 1.0e-10 rtol = 1.0e-10
    @test hybrid_self_transfer.projector !== nothing
    @test hybrid_self_transfer.projector.matrix ≈ hybrid_self_projector.matrix atol = 1.0e-10 rtol = 1.0e-10
    @test hybrid_self_transfer.diagnostics.transfer_path == :hybrid_mixed_raw_cross_overlap_transfer
    @test hybrid_self_transfer.diagnostics.projector_residual_inf < 1.0e-10
    @test hybrid_self_transfer.diagnostics.transferred_residual_inf < 1.0e-10
    @test hybrid_transfer_from_projector.coefficients ≈ hybrid_self_coefficients atol = 1.0e-10 rtol = 1.0e-10
    @test hybrid_transfer_from_projector.projector === hybrid_self_projector
    @test hybrid_full_to_legacy.matrix ≈
          cross_overlap(hybrid_fixture.legacy_rep, hybrid_fixture.full_rep) atol = 1.0e-10 rtol = 1.0e-10
    @test hybrid_parent_to_full.matrix ≈
          cross_overlap(hybrid_fixture.full_rep, hybrid_fixture.fixed_full_rep) atol = 1.0e-10 rtol = 1.0e-10
    @test hybrid_full_to_parent.matrix ≈
          cross_overlap(hybrid_fixture.fixed_full_rep, hybrid_fixture.full_rep) atol = 1.0e-10 rtol = 1.0e-10
end

@testset "Cartesian basis bundle export" begin
    square_basis, _square_source, square_fixed_block, _square_diagnostics =
        _axis_aligned_homonuclear_square_lattice_nested_fixture()
    square_basis_rep = basis_representation(square_basis)

    square_bundle = cartesian_basis_bundle_payload(
        square_basis;
        meta = (example = "test_cartesian_basis_bundle_basis_only",),
    )

    @test square_bundle.basis["format"] == "cartesian_basis_bundle_v1"
    @test square_bundle.basis["version"] == 1
    @test square_bundle.basis["basis_kind"] == "direct_product"
    @test square_bundle.basis["parent_kind"] == "cartesian_product_basis"
    @test square_bundle.basis["contraction_kind"] == "identity"
    @test size(square_bundle.basis["basis_centers"]) == size(square_basis_rep.metadata.basis_centers)
    @test length(square_bundle.basis["final_integral_weights"]) == square_basis_rep.metadata.final_dimension
    @test square_bundle.ham === nothing
    @test !square_bundle.meta["has_ham"]
    @test square_bundle.meta["example"] == "test_cartesian_basis_bundle_basis_only"

    fixed_bundle = cartesian_basis_bundle_payload(square_fixed_block)
    @test fixed_bundle.basis["basis_kind"] == "nested_fixed_block"
    @test fixed_bundle.basis["support_indices_present"]
    @test size(fixed_bundle.basis["support_states"], 2) == 3
    @test fixed_bundle.basis["final_integral_weights"] ≈ square_fixed_block.weights atol = 1.0e-12 rtol = 1.0e-12
    @test fixed_bundle.ham === nothing

    sparse_square_fixed_rep = basis_representation(_with_sparse_nested_coefficients(square_fixed_block))
    sparse_fixed_bundle = cartesian_basis_bundle_payload(sparse_square_fixed_rep)
    @test sparse_fixed_bundle.basis["coefficient_matrix"] isa SparseMatrixCSC{Float64,Int}
    @test Matrix(sparse_fixed_bundle.basis["coefficient_matrix"]) ≈
        Matrix(square_fixed_block.coefficient_matrix) atol = 1.0e-12 rtol = 1.0e-12

    diatomic_basis, diatomic_ops, _diatomic_check = _bond_aligned_diatomic_qw_fixture()
    operator_bundle = cartesian_basis_bundle_payload(
        diatomic_ops;
        meta = (example = "test_cartesian_basis_bundle_with_ham",),
    )
    operator_basis_only_bundle = cartesian_basis_bundle_payload(
        diatomic_ops;
        include_ham = false,
        meta = (example = "test_cartesian_basis_bundle_basis_only_from_ops",),
    )

    @test operator_bundle.basis["basis_kind"] == "direct_product"
    @test operator_bundle.ham !== nothing
    @test operator_bundle.ham["format"] == "cartesian_hamiltonian_bundle_v1"
    @test operator_bundle.ham["model_kind"] == "ordinary_cartesian_operators"
    @test size(operator_bundle.ham["overlap"]) == size(diatomic_ops.overlap)
    @test size(operator_bundle.ham["one_body_hamiltonian"]) == size(diatomic_ops.one_body_hamiltonian)
    @test size(operator_bundle.ham["interaction_matrix"]) == size(diatomic_ops.interaction_matrix)
    @test operator_bundle.ham["nuclear_term_storage"] == "by_center"
    @test operator_bundle.ham["default_nuclear_charges"] == [1.0, 1.0]
    @test operator_bundle.ham["nuclear_one_body_by_center/count"] == 2
    @test size(operator_bundle.ham["kinetic_one_body"]) == size(diatomic_ops.one_body_hamiltonian)
    @test operator_bundle.ham["basis_integral_weights"] == operator_bundle.basis["final_integral_weights"]
    @test operator_bundle.meta["has_ham"]

    mktempdir() do dir
        basis_only_path = joinpath(dir, "square_basis_only.jld2")
        sparse_fixed_path = joinpath(dir, "square_sparse_fixed.jld2")
        ops_path = joinpath(dir, "diatomic_ops_bundle.jld2")
        ops_basis_only_path = joinpath(dir, "diatomic_ops_basis_only_bundle.jld2")

        @test write_cartesian_basis_bundle_jld2(
            basis_only_path,
            square_basis;
            meta = (example = "test_cartesian_basis_bundle_basis_only",),
        ) == basis_only_path
        @test write_cartesian_basis_bundle_jld2(sparse_fixed_path, sparse_square_fixed_rep) ==
            sparse_fixed_path
        @test write_cartesian_basis_bundle_jld2(
            ops_path,
            diatomic_ops;
            meta = (example = "test_cartesian_basis_bundle_with_ham",),
        ) == ops_path
        @test write_cartesian_basis_bundle_jld2(
            ops_basis_only_path,
            diatomic_ops;
            include_ham = false,
            meta = (example = "test_cartesian_basis_bundle_basis_only_from_ops",),
        ) == ops_basis_only_path

        jldopen(basis_only_path, "r") do file
            top_keys = Set(
                key isa AbstractVector ? join(string.(key), "/") : string(key) for key in keys(file)
            )
            @test "basis" in top_keys
            @test "meta" in top_keys
            @test !("ham" in top_keys)
            @test String(file["basis/format"]) == "cartesian_basis_bundle_v1"
            @test Int(file["basis/version"]) == 1
            @test String(file["basis/basis_kind"]) == "direct_product"
            @test size(file["basis/basis_centers"]) == size(square_basis_rep.metadata.basis_centers)
            @test size(file["basis/final_integral_weights"]) == (square_basis_rep.metadata.final_dimension,)
            @test String(file["basis/axes/x/format"]) == "basis_representation_1d_v1"
            @test String(file["meta/producer"]) ==
                "GaussletBases.write_cartesian_basis_bundle_jld2"
        end

        jldopen(sparse_fixed_path, "r") do file
            basis_values = GaussletBases._cartesian_jld_group_values(file["basis"])
            meta_values = GaussletBases._cartesian_jld_group_values(file["meta"])
            @test file["basis/coefficient_matrix"] isa SparseMatrixCSC{Float64,Int}
            @test Set(keys(basis_values)) == Set(keys(sparse_fixed_bundle.basis))
            @test Set(keys(meta_values)) == Set(keys(sparse_fixed_bundle.meta))
            @test basis_values["final_integral_weights"] ≈
                sparse_fixed_bundle.basis["final_integral_weights"] atol = 1.0e-12 rtol = 1.0e-12
        end

        sparse_fixed_bundle_roundtrip = read_cartesian_basis_bundle(sparse_fixed_path)
        @test sparse_fixed_bundle_roundtrip.basis.coefficient_matrix isa SparseMatrixCSC{Float64,Int}
        @test sparse_fixed_bundle_roundtrip.basis.coefficient_matrix ==
            sparse_square_fixed_rep.coefficient_matrix
        @test cross_overlap(sparse_fixed_bundle_roundtrip, sparse_fixed_bundle_roundtrip) ≈
            cross_overlap(sparse_square_fixed_rep, sparse_square_fixed_rep) atol = 1.0e-10 rtol = 1.0e-10

        jldopen(ops_path, "r") do file
            top_keys = Set(
                key isa AbstractVector ? join(string.(key), "/") : string(key) for key in keys(file)
            )
            ham_values = GaussletBases._cartesian_jld_group_values(file["ham"])
            @test "basis" in top_keys
            @test "ham" in top_keys
            @test "meta" in top_keys
            @test Set(keys(ham_values)) == Set(keys(operator_bundle.ham))
            @test String(file["ham/format"]) == "cartesian_hamiltonian_bundle_v1"
            @test String(file["ham/model_kind"]) == "ordinary_cartesian_operators"
            @test size(file["ham/overlap"]) == size(diatomic_ops.overlap)
            @test size(file["ham/one_body_hamiltonian"]) == size(diatomic_ops.one_body_hamiltonian)
            @test size(file["ham/interaction_matrix"]) == size(diatomic_ops.interaction_matrix)
            @test String(file["ham/nuclear_term_storage"]) == "by_center"
            @test Int(file["ham/nuclear_one_body_by_center/count"]) == 2
            @test size(file["ham/kinetic_one_body"]) == size(diatomic_ops.one_body_hamiltonian)
            @test String(file["meta/manifest/contract/format"]) == "cartesian_basis_bundle_v1"
            @test Bool(file["meta/has_ham"])
        end

        jldopen(ops_basis_only_path, "r") do file
            top_keys = Set(
                key isa AbstractVector ? join(string.(key), "/") : string(key) for key in keys(file)
            )
            basis_values = GaussletBases._cartesian_jld_group_values(file["basis"])
            meta_values = GaussletBases._cartesian_jld_group_values(file["meta"])
            @test "basis" in top_keys
            @test "meta" in top_keys
            @test !("ham" in top_keys)
            @test Set(keys(basis_values)) == Set(keys(operator_basis_only_bundle.basis))
            @test Set(keys(meta_values)) == Set(keys(operator_basis_only_bundle.meta))
            @test !Bool(file["meta/has_ham"])
            @test String(file["meta/example"]) == "test_cartesian_basis_bundle_basis_only_from_ops"
        end

        ops_basis_only_bundle_roundtrip = read_cartesian_basis_bundle(ops_basis_only_path)
        @test ops_basis_only_bundle_roundtrip.ham === nothing
        @test cross_overlap(ops_basis_only_bundle_roundtrip, ops_basis_only_bundle_roundtrip) ≈
            diatomic_ops.overlap atol = 1.0e-10 rtol = 1.0e-10
    end

    hybrid_fixture = _atomic_hybrid_cartesian_representation_fixture()
    hybrid_bundle = cartesian_basis_bundle_payload(
        hybrid_fixture.full_ops;
        meta = (example = "test_cartesian_hybrid_bundle",),
    )

    @test hybrid_bundle.basis["basis_kind"] == "hybrid_residual"
    @test hybrid_bundle.basis["parent_kind"] == "cartesian_plus_supplement_raw"
    @test hybrid_bundle.basis["parent/format"] == "cartesian_plus_supplement_raw_v1"
    @test hybrid_bundle.basis["parent/cartesian/format"] == "cartesian_basis_bundle_v1"
    @test hybrid_bundle.basis["parent/supplement/format"] ==
        "cartesian_gaussian_shell_supplement_v1"
    @test haskey(hybrid_bundle.basis, "parent/cartesian_supplement_axis_tables/x")
    @test haskey(hybrid_bundle.basis, "parent/cartesian_supplement_axis_tables/y")
    @test haskey(hybrid_bundle.basis, "parent/cartesian_supplement_axis_tables/z")
    @test haskey(hybrid_bundle.basis, "parent/exact_cartesian_supplement_overlap")
    @test haskey(hybrid_bundle.basis, "parent/exact_supplement_overlap")
    @test hybrid_bundle.basis["parent/supplement/orbital_count"] ==
        size(hybrid_fixture.full_ops.raw_to_final, 1) - hybrid_fixture.full_ops.gausslet_count
    @test hybrid_bundle.ham !== nothing
    @test hybrid_bundle.ham["model_kind"] == "ordinary_cartesian_operators"
    @test size(hybrid_bundle.ham["overlap"]) == size(hybrid_fixture.full_ops.overlap)
    @test hybrid_bundle.meta["example"] == "test_cartesian_hybrid_bundle"

    mktempdir() do dir
        hybrid_path = joinpath(dir, "atomic_hybrid_ops_bundle.jld2")

        @test write_cartesian_basis_bundle_jld2(
            hybrid_path,
            hybrid_fixture.full_ops;
            meta = (example = "test_cartesian_hybrid_bundle",),
        ) == hybrid_path

        jldopen(hybrid_path, "r") do file
            @test String(file["basis/parent_kind"]) == "cartesian_plus_supplement_raw"
            @test String(file["basis/parent/format"]) == "cartesian_plus_supplement_raw_v1"
            @test String(file["basis/parent/cartesian/format"]) == "cartesian_basis_bundle_v1"
            @test String(file["basis/parent/supplement/format"]) ==
                "cartesian_gaussian_shell_supplement_v1"
            @test size(file["basis/parent/cartesian_supplement_axis_tables/x"], 2) ==
                size(hybrid_fixture.full_ops.raw_to_final, 1) - hybrid_fixture.full_ops.gausslet_count
            @test size(file["basis/parent/cartesian_supplement_axis_tables/y"], 2) ==
                size(hybrid_fixture.full_ops.raw_to_final, 1) - hybrid_fixture.full_ops.gausslet_count
            @test size(file["basis/parent/cartesian_supplement_axis_tables/z"], 2) ==
                size(hybrid_fixture.full_ops.raw_to_final, 1) - hybrid_fixture.full_ops.gausslet_count
            @test size(file["basis/parent/exact_cartesian_supplement_overlap"]) ==
                (hybrid_fixture.full_ops.gausslet_count,
                 size(hybrid_fixture.full_ops.raw_to_final, 1) - hybrid_fixture.full_ops.gausslet_count)
            @test size(file["basis/parent/exact_supplement_overlap"]) ==
                (size(hybrid_fixture.full_ops.raw_to_final, 1) - hybrid_fixture.full_ops.gausslet_count,
                 size(hybrid_fixture.full_ops.raw_to_final, 1) - hybrid_fixture.full_ops.gausslet_count)
            @test Int(file["basis/parent/supplement/orbital_count"]) ==
                size(hybrid_fixture.full_ops.raw_to_final, 1) - hybrid_fixture.full_ops.gausslet_count
            @test size(file["ham/overlap"]) == size(hybrid_fixture.full_ops.overlap)
            @test size(file["ham/one_body_hamiltonian"]) ==
                size(hybrid_fixture.full_ops.one_body_hamiltonian)
        end
    end

    bond_aligned_hybrid_fixture = _bond_aligned_diatomic_nested_hybrid_bundle_fixture()
    bond_aligned_hybrid_trimmed_fixture =
        _bond_aligned_diatomic_nested_hybrid_bundle_fixture(; max_width = 1.0)
    bond_aligned_hybrid_supplement3d =
        GaussletBases._bond_aligned_diatomic_cartesian_shell_supplement_3d(
            bond_aligned_hybrid_fixture.supplement,
        )
    bond_aligned_hybrid_bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(
        bond_aligned_hybrid_fixture.basis,
        bond_aligned_hybrid_fixture.hybrid_ops.expansion;
        gausslet_backend = bond_aligned_hybrid_fixture.hybrid_ops.gausslet_backend,
    )
    bond_aligned_hybrid_overlap_blocks =
        GaussletBases._qwrg_diatomic_cartesian_shell_overlap_blocks_3d(
            bond_aligned_hybrid_bundles,
            bond_aligned_hybrid_supplement3d,
            bond_aligned_hybrid_fixture.basis,
            bond_aligned_hybrid_fixture.hybrid_ops.expansion,
        )
    bond_aligned_hybrid_full_blocks = GaussletBases._qwrg_diatomic_cartesian_shell_blocks_3d(
        bond_aligned_hybrid_bundles,
        bond_aligned_hybrid_supplement3d,
        bond_aligned_hybrid_fixture.basis,
        bond_aligned_hybrid_fixture.hybrid_ops.expansion,
        bond_aligned_hybrid_fixture.hybrid_ops.nuclear_charges,
    )
    bond_aligned_hybrid_bundle = cartesian_basis_bundle_payload(
        bond_aligned_hybrid_fixture.hybrid_ops;
        meta = (example = "test_cartesian_bond_aligned_diatomic_hybrid_bundle",),
    )
    bond_aligned_hybrid_trimmed_bundle = cartesian_basis_bundle_payload(
        bond_aligned_hybrid_trimmed_fixture.hybrid_ops;
        meta = (example = "test_cartesian_bond_aligned_diatomic_hybrid_bundle_trimmed",),
    )

    @test bond_aligned_hybrid_bundle.basis["basis_kind"] == "hybrid_residual"
    @test bond_aligned_hybrid_bundle.basis["parent_kind"] == "cartesian_plus_supplement_raw"
    @test bond_aligned_hybrid_bundle.basis["parent/format"] == "cartesian_plus_supplement_raw_v1"
    @test bond_aligned_hybrid_bundle.basis["parent/supplement/format"] ==
        "cartesian_gaussian_shell_supplement_v1"
    @test haskey(bond_aligned_hybrid_bundle.basis, "parent/cartesian_supplement_axis_tables/x")
    @test haskey(bond_aligned_hybrid_bundle.basis, "parent/exact_cartesian_supplement_overlap")
    @test haskey(bond_aligned_hybrid_bundle.basis, "parent/exact_supplement_overlap")
    @test bond_aligned_hybrid_overlap_blocks.overlap_ga ≈
        bond_aligned_hybrid_full_blocks.overlap_ga atol = 1.0e-12 rtol = 1.0e-12
    @test bond_aligned_hybrid_overlap_blocks.overlap_aa ≈
        bond_aligned_hybrid_full_blocks.overlap_aa atol = 1.0e-12 rtol = 1.0e-12
    @test bond_aligned_hybrid_trimmed_bundle.basis["parent/supplement/metadata/max_width"] == 1.0
    @test Int(bond_aligned_hybrid_trimmed_bundle.basis["parent/supplement/orbital_count"]) <
        Int(bond_aligned_hybrid_bundle.basis["parent/supplement/orbital_count"])

    mktempdir() do dir
        hybrid_path = joinpath(dir, "bond_aligned_diatomic_hybrid_ops_bundle.jld2")
        hybrid_trimmed_path =
            joinpath(dir, "bond_aligned_diatomic_hybrid_ops_bundle_trimmed.jld2")

        @test write_cartesian_basis_bundle_jld2(
            hybrid_path,
            bond_aligned_hybrid_fixture.hybrid_ops;
            meta = (example = "test_cartesian_bond_aligned_diatomic_hybrid_bundle",),
        ) == hybrid_path
        @test write_cartesian_basis_bundle_jld2(
            hybrid_trimmed_path,
            bond_aligned_hybrid_trimmed_fixture.hybrid_ops;
            meta = (example = "test_cartesian_bond_aligned_diatomic_hybrid_bundle_trimmed",),
        ) == hybrid_trimmed_path

        jldopen(hybrid_path, "r") do file
            @test String(file["basis/parent_kind"]) == "cartesian_plus_supplement_raw"
            @test String(file["basis/parent/format"]) == "cartesian_plus_supplement_raw_v1"
            @test String(file["basis/parent/supplement/format"]) ==
                "cartesian_gaussian_shell_supplement_v1"
            @test Int(file["basis/parent/supplement/orbital_count"]) ==
                size(bond_aligned_hybrid_fixture.hybrid_ops.raw_to_final, 1) -
                bond_aligned_hybrid_fixture.hybrid_ops.gausslet_count
            @test size(file["basis/parent/exact_cartesian_supplement_overlap"]) ==
                (bond_aligned_hybrid_fixture.hybrid_ops.gausslet_count,
                 size(bond_aligned_hybrid_fixture.hybrid_ops.raw_to_final, 1) -
                 bond_aligned_hybrid_fixture.hybrid_ops.gausslet_count)
            @test size(file["basis/parent/exact_supplement_overlap"]) ==
                (size(bond_aligned_hybrid_fixture.hybrid_ops.raw_to_final, 1) -
                 bond_aligned_hybrid_fixture.hybrid_ops.gausslet_count,
                 size(bond_aligned_hybrid_fixture.hybrid_ops.raw_to_final, 1) -
                 bond_aligned_hybrid_fixture.hybrid_ops.gausslet_count)
            @test String(file["meta/example"]) ==
                "test_cartesian_bond_aligned_diatomic_hybrid_bundle"
        end

        jldopen(hybrid_trimmed_path, "r") do file
            @test Float64(file["basis/parent/supplement/metadata/max_width"]) == 1.0
            @test Int(file["basis/parent/supplement/orbital_count"]) <
                Int(bond_aligned_hybrid_bundle.basis["parent/supplement/orbital_count"])
            @test String(file["meta/example"]) ==
                "test_cartesian_bond_aligned_diatomic_hybrid_bundle_trimmed"
        end
    end
end

@testset "Cartesian basis bundle overlap and projector" begin
    square_basis, _square_source, square_fixed_block, _square_diagnostics =
        _axis_aligned_homonuclear_square_lattice_nested_fixture()
    square_basis_rep = basis_representation(square_basis)
    square_fixed_rep = basis_representation(square_fixed_block)

    diatomic_basis14, diatomic_ops14, _diatomic_check14 =
        _bond_aligned_diatomic_qw_fixture(; bond_length = 1.4)
    diatomic_basis20, _diatomic_ops20, _diatomic_check20 =
        _bond_aligned_diatomic_qw_fixture(; bond_length = 2.0)
    diatomic_rep14 = basis_representation(diatomic_basis14)
    diatomic_rep20 = basis_representation(diatomic_basis20)
    bond_aligned_hybrid_fixture = _bond_aligned_diatomic_nested_hybrid_bundle_fixture()

    dir = mktempdir()
    try
        square_path = joinpath(dir, "square_basis.jld2")
        square_fixed_path = joinpath(dir, "square_fixed.jld2")
        diatomic14_path = joinpath(dir, "diatomic14.jld2")
        diatomic20_path = joinpath(dir, "diatomic20.jld2")
        diatomic_ops_path = joinpath(dir, "diatomic_ops.jld2")
        atomic_fixed_full_path = joinpath(dir, "atomic_fixed_full.jld2")
        atomic_hybrid_full_path = joinpath(dir, "atomic_hybrid_full.jld2")
        atomic_hybrid_legacy_path = joinpath(dir, "atomic_hybrid_legacy.jld2")
        bond_aligned_hybrid_fixed_path = joinpath(dir, "bond_aligned_hybrid_fixed.jld2")
        bond_aligned_hybrid_path = joinpath(dir, "bond_aligned_hybrid_ops.jld2")

        write_cartesian_basis_bundle_jld2(square_path, square_basis)
        write_cartesian_basis_bundle_jld2(square_fixed_path, square_fixed_block)
        write_cartesian_basis_bundle_jld2(diatomic14_path, diatomic_basis14)
        write_cartesian_basis_bundle_jld2(diatomic20_path, diatomic_basis20)
        write_cartesian_basis_bundle_jld2(diatomic_ops_path, diatomic_ops14)
        hybrid_fixture = _atomic_hybrid_cartesian_representation_fixture()
        write_cartesian_basis_bundle_jld2(atomic_fixed_full_path, hybrid_fixture.fixed_full)
        write_cartesian_basis_bundle_jld2(atomic_hybrid_full_path, hybrid_fixture.full_ops)
        write_cartesian_basis_bundle_jld2(atomic_hybrid_legacy_path, hybrid_fixture.legacy_ops)
        write_cartesian_basis_bundle_jld2(
            bond_aligned_hybrid_fixed_path,
            bond_aligned_hybrid_fixture.fixed_block,
        )
        write_cartesian_basis_bundle_jld2(
            bond_aligned_hybrid_path,
            bond_aligned_hybrid_fixture.hybrid_ops,
        )

        square_bundle = read_cartesian_basis_bundle(square_path)
        square_fixed_bundle = read_cartesian_basis_bundle(square_fixed_path)
        diatomic14_bundle = read_cartesian_basis_bundle(diatomic14_path)
        diatomic20_bundle = read_cartesian_basis_bundle(diatomic20_path)
        diatomic_ops_bundle = read_cartesian_basis_bundle(diatomic_ops_path)
        atomic_hybrid_full_bundle = read_cartesian_basis_bundle(atomic_hybrid_full_path)
        atomic_hybrid_legacy_bundle = read_cartesian_basis_bundle(atomic_hybrid_legacy_path)
        bond_aligned_hybrid_bundle = read_cartesian_basis_bundle(bond_aligned_hybrid_path)

        @test square_bundle.path == abspath(square_path)
        @test square_bundle.diagnostics.basis_kind == :direct_product
        @test square_bundle.diagnostics.final_dimension == square_basis_rep.metadata.final_dimension
        @test square_bundle.ham === nothing
        @test diatomic_ops_bundle.ham !== nothing
        @test diatomic_ops_bundle.ham["model_kind"] == "ordinary_cartesian_operators"
        @test diatomic_ops_bundle.diagnostics.has_ham

        diatomic_atom_a_ops = ordinary_cartesian_qiu_white_operators(
            diatomic_basis14;
            nuclear_charges = [1.0, 0.0],
            nuclear_term_storage = :total_only,
            interaction_treatment = :ggt_nearest,
        )
        diatomic_atom_b_ops = ordinary_cartesian_qiu_white_operators(
            diatomic_basis14;
            nuclear_charges = [0.0, 1.0],
            nuclear_term_storage = :total_only,
            interaction_treatment = :ggt_nearest,
        )
        @test assembled_one_body_hamiltonian(diatomic_ops14) ≈
              diatomic_ops14.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
        @test assembled_one_body_hamiltonian(diatomic_ops_bundle) ≈
              diatomic_ops14.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
        @test assembled_one_body_hamiltonian(diatomic_ops14; nuclear_charges = [1.0, 0.0]) ≈
              diatomic_atom_a_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(diatomic_ops_bundle; nuclear_charges = [1.0, 0.0]) ≈
              diatomic_atom_a_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(diatomic_ops14; nuclear_charges = [0.0, 1.0]) ≈
              diatomic_atom_b_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(diatomic_ops_bundle; nuclear_charges = [0.0, 1.0]) ≈
              diatomic_atom_b_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10

        hybrid_atom_a_ops = ordinary_cartesian_qiu_white_operators(
            bond_aligned_hybrid_fixture.fixed_block,
            bond_aligned_hybrid_fixture.supplement;
            nuclear_charges = [1.0, 0.0],
            nuclear_term_storage = :total_only,
            interaction_treatment = :ggt_nearest,
        )
        hybrid_atom_b_ops = ordinary_cartesian_qiu_white_operators(
            bond_aligned_hybrid_fixture.fixed_block,
            bond_aligned_hybrid_fixture.supplement;
            nuclear_charges = [0.0, 1.0],
            nuclear_term_storage = :total_only,
            interaction_treatment = :ggt_nearest,
        )
        @test assembled_one_body_hamiltonian(bond_aligned_hybrid_fixture.hybrid_ops) ≈
              bond_aligned_hybrid_fixture.hybrid_ops.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
        @test assembled_one_body_hamiltonian(bond_aligned_hybrid_bundle) ≈
              bond_aligned_hybrid_fixture.hybrid_ops.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
        @test assembled_one_body_hamiltonian(
                  bond_aligned_hybrid_fixture.hybrid_ops;
                  nuclear_charges = [1.0, 0.0],
              ) ≈ hybrid_atom_a_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(
                  bond_aligned_hybrid_bundle;
                  nuclear_charges = [1.0, 0.0],
              ) ≈ hybrid_atom_a_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(
                  bond_aligned_hybrid_fixture.hybrid_ops;
                  nuclear_charges = [0.0, 1.0],
              ) ≈ hybrid_atom_b_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
        @test assembled_one_body_hamiltonian(
                  bond_aligned_hybrid_bundle;
                  nuclear_charges = [0.0, 1.0],
              ) ≈ hybrid_atom_b_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10

        loaded_square_rep = load_cartesian_basis_representation(square_path)
        @test loaded_square_rep.metadata.basis_kind == square_basis_rep.metadata.basis_kind
        @test loaded_square_rep.metadata.final_dimension == square_basis_rep.metadata.final_dimension
        @test loaded_square_rep.metadata.parent_kind == square_basis_rep.metadata.parent_kind

        direct_self_disk = cross_overlap(square_bundle, square_bundle)
        direct_cross_disk = cross_overlap(diatomic14_bundle, diatomic20_bundle)
        nested_cross_disk = cross_overlap(square_path, square_fixed_path)

        @test direct_self_disk ≈ cross_overlap(square_basis_rep, square_basis_rep) atol = 1.0e-10 rtol = 1.0e-10
        @test direct_cross_disk ≈ cross_overlap(diatomic_rep14, diatomic_rep20) atol = 1.0e-10 rtol = 1.0e-10
        @test nested_cross_disk ≈ cross_overlap(square_basis_rep, square_fixed_rep) atol = 1.0e-10 rtol = 1.0e-10

        disk_projector = basis_projector(square_fixed_path, square_path)
        memory_projector = basis_projector(square_fixed_rep, square_basis_rep)
        @test disk_projector.matrix ≈ memory_projector.matrix atol = 1.0e-10 rtol = 1.0e-10
        @test disk_projector.diagnostics.transfer_path == memory_projector.diagnostics.transfer_path

        fixed_coefficients = cos.(Float64.(1:square_fixed_rep.metadata.final_dimension))
        disk_transfer = transfer_orbitals(fixed_coefficients, square_fixed_path, square_path)
        memory_transfer = transfer_orbitals(fixed_coefficients, square_fixed_rep, square_basis_rep)
        @test disk_transfer.coefficients ≈ memory_transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test disk_transfer.diagnostics.transferred_residual_inf < 1.0e-10

        @test atomic_hybrid_full_bundle.diagnostics.parent_kind == :cartesian_plus_supplement_raw
        @test atomic_hybrid_full_bundle.ham !== nothing
        @test hasproperty(
            atomic_hybrid_full_bundle.basis.parent_data,
            :cartesian_supplement_axis_tables,
        )
        @test bond_aligned_hybrid_bundle.diagnostics.parent_kind == :cartesian_plus_supplement_raw
        @test bond_aligned_hybrid_bundle.ham !== nothing
        @test hasproperty(
            bond_aligned_hybrid_bundle.basis.parent_data,
            :cartesian_supplement_axis_tables,
        )

        disk_hybrid_self = cross_overlap(atomic_hybrid_full_path, atomic_hybrid_full_path)
        memory_hybrid_self = cross_overlap(hybrid_fixture.full_rep, hybrid_fixture.full_rep)
        @test disk_hybrid_self ≈ memory_hybrid_self atol = 1.0e-10 rtol = 1.0e-10

        disk_hybrid_cross = cross_overlap(atomic_hybrid_full_path, atomic_hybrid_legacy_path)
        memory_hybrid_cross = cross_overlap(hybrid_fixture.full_rep, hybrid_fixture.legacy_rep)
        @test disk_hybrid_cross ≈ memory_hybrid_cross atol = 1.0e-10 rtol = 1.0e-10

        disk_hybrid_parent = cross_overlap(atomic_fixed_full_path, atomic_hybrid_full_path)
        memory_hybrid_parent = cross_overlap(hybrid_fixture.fixed_full_rep, hybrid_fixture.full_rep)
        @test disk_hybrid_parent ≈ memory_hybrid_parent atol = 1.0e-10 rtol = 1.0e-10

        disk_hybrid_projector =
            basis_projector(atomic_hybrid_full_path, atomic_hybrid_legacy_path)
        memory_hybrid_projector =
            basis_projector(hybrid_fixture.full_rep, hybrid_fixture.legacy_rep)
        @test disk_hybrid_projector.matrix ≈ memory_hybrid_projector.matrix atol = 1.0e-10 rtol = 1.0e-10
        @test disk_hybrid_projector.diagnostics.transfer_path ==
            memory_hybrid_projector.diagnostics.transfer_path

        hybrid_coefficients = cos.(Float64.(1:hybrid_fixture.full_rep.metadata.final_dimension))
        disk_hybrid_transfer = transfer_orbitals(
            hybrid_coefficients,
            atomic_hybrid_full_path,
            atomic_hybrid_legacy_path,
        )
        memory_hybrid_transfer = transfer_orbitals(
            hybrid_coefficients,
            hybrid_fixture.full_rep,
            hybrid_fixture.legacy_rep,
        )
        disk_hybrid_fast_transfer = transfer_orbitals(
            hybrid_coefficients,
            atomic_hybrid_full_path,
            atomic_hybrid_legacy_path;
            materialize_projector = false,
        )
        bundle_hybrid_fast_transfer = transfer_orbitals(
            hybrid_coefficients,
            atomic_hybrid_full_bundle,
            atomic_hybrid_legacy_bundle;
            materialize_projector = false,
        )
        memory_hybrid_fast_transfer = transfer_orbitals(
            hybrid_coefficients,
            hybrid_fixture.full_rep,
            hybrid_fixture.legacy_rep;
            materialize_projector = false,
        )
        @test disk_hybrid_transfer.coefficients ≈
            memory_hybrid_transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test disk_hybrid_transfer.projector !== nothing
        @test memory_hybrid_transfer.projector !== nothing
        @test disk_hybrid_transfer.projector.matrix ≈
            memory_hybrid_transfer.projector.matrix atol = 1.0e-10 rtol = 1.0e-10
        @test disk_hybrid_transfer.diagnostics.transferred_residual_inf < 1.0e-10
        @test memory_hybrid_fast_transfer.coefficients ≈
            memory_hybrid_transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test disk_hybrid_fast_transfer.coefficients ≈
            memory_hybrid_fast_transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test bundle_hybrid_fast_transfer.coefficients ≈
            memory_hybrid_fast_transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test disk_hybrid_fast_transfer.projector === nothing
        @test bundle_hybrid_fast_transfer.projector === nothing
        @test memory_hybrid_fast_transfer.projector === nothing
        @test disk_hybrid_fast_transfer.diagnostics.transfer_path == :hybrid_mixed_raw_cross_overlap_transfer
        @test bundle_hybrid_fast_transfer.diagnostics.transfer_path == :hybrid_mixed_raw_cross_overlap_transfer
        @test isnan(disk_hybrid_fast_transfer.diagnostics.projector_residual_inf)
        @test isnan(bundle_hybrid_fast_transfer.diagnostics.projector_residual_inf)
        @test isnan(disk_hybrid_fast_transfer.diagnostics.transferred_residual_inf)
        @test isnan(bundle_hybrid_fast_transfer.diagnostics.transferred_residual_inf)

        disk_bond_aligned_hybrid_self =
            cross_overlap(bond_aligned_hybrid_path, bond_aligned_hybrid_path)
        memory_bond_aligned_hybrid_self =
            cross_overlap(
                bond_aligned_hybrid_fixture.hybrid_rep,
                bond_aligned_hybrid_fixture.hybrid_rep,
            )
        @test disk_bond_aligned_hybrid_self ≈
            memory_bond_aligned_hybrid_self atol = 1.0e-10 rtol = 1.0e-10

        disk_bond_aligned_hybrid_cross =
            cross_overlap(bond_aligned_hybrid_fixed_path, bond_aligned_hybrid_path)
        memory_bond_aligned_hybrid_cross =
            cross_overlap(
                bond_aligned_hybrid_fixture.fixed_rep,
                bond_aligned_hybrid_fixture.hybrid_rep,
            )
        @test disk_bond_aligned_hybrid_cross ≈
            memory_bond_aligned_hybrid_cross atol = 1.0e-10 rtol = 1.0e-10

        disk_bond_aligned_hybrid_projector =
            basis_projector(bond_aligned_hybrid_fixed_path, bond_aligned_hybrid_path)
        memory_bond_aligned_hybrid_projector =
            basis_projector(
                bond_aligned_hybrid_fixture.fixed_rep,
                bond_aligned_hybrid_fixture.hybrid_rep,
            )
        @test disk_bond_aligned_hybrid_projector.matrix ≈
            memory_bond_aligned_hybrid_projector.matrix atol = 1.0e-10 rtol = 1.0e-10
        @test disk_bond_aligned_hybrid_projector.diagnostics.transfer_path ==
            memory_bond_aligned_hybrid_projector.diagnostics.transfer_path

        bond_aligned_coefficients =
            cos.(Float64.(1:bond_aligned_hybrid_fixture.fixed_rep.metadata.final_dimension))
        disk_bond_aligned_hybrid_transfer = transfer_orbitals(
            bond_aligned_coefficients,
            bond_aligned_hybrid_fixed_path,
            bond_aligned_hybrid_path,
        )
        memory_bond_aligned_hybrid_transfer = transfer_orbitals(
            bond_aligned_coefficients,
            bond_aligned_hybrid_fixture.fixed_rep,
            bond_aligned_hybrid_fixture.hybrid_rep,
        )
        @test disk_bond_aligned_hybrid_transfer.coefficients ≈
            memory_bond_aligned_hybrid_transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test disk_bond_aligned_hybrid_transfer.projector.matrix ≈
            memory_bond_aligned_hybrid_transfer.projector.matrix atol = 1.0e-10 rtol = 1.0e-10
        @test disk_bond_aligned_hybrid_transfer.diagnostics.transferred_residual_inf < 1.0e-10
    finally
        rm(dir; recursive = true, force = true)
    end
end

@testset "Atomic direct-product He extent change is not an outer-only identity" begin
    fixture = _atomic_direct_product_he_extent_change_contract_fixture()

    @test fixture.source_count == 3
    @test fixture.target_count == 5
    @test fixture.shared_slice == 2:4
    @test fixture.centers_subset
    @test !fixture.weights_subset
    @test !fixture.coefficient_core_match
end

@testset "Atomic hybrid He orbital transfer remains stable across same-parent different-final-contraction change" begin
    fixture = _atomic_hybrid_he_same_parent_stress_fixture()

    @test fixture.source_ops isa OrdinaryCartesianOperators3D
    @test fixture.target_ops isa OrdinaryCartesianOperators3D
    @test fixture.source_rep.metadata.parent_kind == :cartesian_plus_supplement_raw
    @test fixture.target_rep.metadata.parent_kind == :cartesian_plus_supplement_raw
    @test fixture.source_working_box == 2:6
    @test fixture.target_working_box == (1:7, 1:7, 1:7)
    @test length(fixture.source_ops.orbital_data) == 134
    @test length(fixture.target_ops.orbital_data) == 232
    @test fixture.transfer.diagnostics.transfer_path == :hybrid_mixed_raw_cross_overlap_transfer
    @test fixture.transfer.diagnostics.transferred_residual_inf < 1.0e-10

    @test fixture.source_observables.metric_norm_error < 1.0e-12
    @test fixture.target_observables.metric_norm_error < 1.0e-12
    @test fixture.aligned_transferred_observables.metric_norm_error < 1.0e-12

    source_self_overlap = cross_overlap(fixture.source_rep, fixture.source_rep)
    target_self_overlap = cross_overlap(fixture.target_rep, fixture.target_rep)
    cross_overlap_source_target = cross_overlap(fixture.source_rep, fixture.target_rep)
    @test source_self_overlap ≈ fixture.source_ops.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test target_self_overlap ≈ fixture.target_ops.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test maximum(abs, cross_overlap_source_target) <= 1.0 + 1.0e-10
    @test maximum(svdvals(cross_overlap_source_target)) <= 1.0 + 1.0e-10

    @test fixture.target_observables.total < fixture.source_observables.total
    @test fixture.aligned_overlap_to_target > 0.999995
    @test abs(
        fixture.aligned_transferred_observables.one_body - fixture.target_observables.one_body,
    ) < 1.0e-4
    @test abs(
        fixture.aligned_transferred_observables.vee - fixture.target_observables.vee,
    ) < 2.0e-4
    @test abs(
        fixture.aligned_transferred_observables.total - fixture.target_observables.total,
    ) < 5.0e-4

    mktempdir() do dir
        source_path = joinpath(dir, "he_source_hybrid.jld2")
        target_path = joinpath(dir, "he_target_hybrid.jld2")

        write_cartesian_basis_bundle_jld2(source_path, fixture.source_ops)
        write_cartesian_basis_bundle_jld2(target_path, fixture.target_ops)

        source_bundle = read_cartesian_basis_bundle(source_path)
        target_bundle = read_cartesian_basis_bundle(target_path)
        @test hasproperty(source_bundle.basis.parent_data, :exact_cartesian_supplement_overlap)
        @test hasproperty(source_bundle.basis.parent_data, :exact_supplement_overlap)
        @test hasproperty(target_bundle.basis.parent_data, :exact_cartesian_supplement_overlap)
        @test hasproperty(target_bundle.basis.parent_data, :exact_supplement_overlap)

        disk_source_self = cross_overlap(source_path, source_path)
        disk_target_self = cross_overlap(target_path, target_path)
        disk_cross = cross_overlap(source_path, target_path)
        @test disk_source_self ≈ fixture.source_ops.overlap atol = 1.0e-10 rtol = 1.0e-10
        @test disk_target_self ≈ fixture.target_ops.overlap atol = 1.0e-10 rtol = 1.0e-10
        @test disk_source_self ≈ source_self_overlap atol = 1.0e-10 rtol = 1.0e-10
        @test disk_target_self ≈ target_self_overlap atol = 1.0e-10 rtol = 1.0e-10
        @test disk_cross ≈ cross_overlap_source_target atol = 1.0e-10 rtol = 1.0e-10
        @test maximum(abs, disk_cross) <= 1.0 + 1.0e-10
        @test maximum(svdvals(disk_cross)) <= 1.0 + 1.0e-10

        disk_transfer = transfer_orbitals(
            fixture.source_observables.orbital,
            source_path,
            target_path,
        )

        @test disk_transfer.coefficients ≈ fixture.transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test disk_transfer.diagnostics.transfer_path ==
            fixture.transfer.diagnostics.transfer_path
    end
end

@testset "Mapped ordinary Cartesian 1D working representation uses localized Gaussian contract" begin
    mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0)
    basis_a = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 5,
            mapping = mapping,
            reference_spacing = 1.0,
        ),
    )
    basis_b = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 7,
            mapping = mapping,
            reference_spacing = 1.0,
        ),
    )

    rep_a = GaussletBases._mapped_ordinary_working_basis_representation(basis_a)
    rep_b = GaussletBases._mapped_ordinary_working_basis_representation(basis_b)
    S_AA = cross_overlap(rep_a, rep_a)
    S_BB = cross_overlap(rep_b, rep_b)
    S_AB = cross_overlap(rep_a, rep_b)
    I_A = Matrix{Float64}(I, size(S_AA, 1), size(S_AA, 2))
    I_B = Matrix{Float64}(I, size(S_BB, 1), size(S_BB, 2))

    @test all(primitive -> primitive isa Gaussian, primitives(primitive_set(rep_a)))
    @test all(primitive -> primitive isa Gaussian, primitives(primitive_set(rep_b)))
    @test norm(S_AA - I_A, Inf) < 1.0e-12
    @test norm(S_BB - I_B, Inf) < 1.0e-12
    @test maximum(svdvals(S_AB)) <= 1.0 + 1.0e-10

    fixture = _atomic_hybrid_he_same_parent_stress_fixture()
    @test all(
        primitive -> primitive isa Gaussian,
        primitives(primitive_set(fixture.source_rep.axis_representations.x)),
    )
    @test all(
        primitive -> primitive isa Gaussian,
        primitives(primitive_set(fixture.target_rep.axis_representations.x)),
    )
end

@testset "One-center atomic factorized direct packet kernel" begin
    basis = build_basis(
        MappedUniformBasisSpec(
            :G10;
            count = 13,
            mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.2, tail_spacing = 10.0),
            reference_spacing = 1.0,
        ),
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )

    shell = GaussletBases._nested_rectangular_shell(
        bundle,
        2:12,
        2:12,
        2:12;
        retain_xy = (4, 3),
        retain_xz = (4, 3),
        retain_yz = (4, 3),
        packet_kernel = :factorized_direct,
        term_coefficients = term_coefficients,
    )
    factorized = GaussletBases._nested_extract_factorized_basis(
        shell.coefficient_matrix,
        (13, 13, 13),
    )
    reconstructed = GaussletBases._nested_reconstruct_factorized_coefficients(factorized)
    @test factorized.reconstruction_max_error < 1.0e-10
    @test reconstructed ≈ shell.coefficient_matrix atol = 1.0e-10 rtol = 1.0e-10

    dims = (3, 2, 2)
    x1 = [1.0, 0.5, 0.0]
    x2 = [1.0, -0.25, 0.2]
    y1 = [1.0, 0.0]
    y2 = [1.0, 0.3]
    z1 = [1.0, -0.4]
    z2 = [1.0, 0.25]
    amplitudes = [2.0, -1.5, 0.75]
    factorable = zeros(Float64, prod(dims), 3)
    factors = ((x1, y1, z1), (x1, y2, z1), (x2, y1, z2))
    for column in 1:3
        xvec, yvec, zvec = factors[column]
        amplitude = amplitudes[column]
        for ix in 1:dims[1], iy in 1:dims[2], iz in 1:dims[3]
            factorable[GaussletBases._cartesian_flat_index(ix, iy, iz, dims), column] =
                amplitude * xvec[ix] * yvec[iy] * zvec[iz]
        end
    end
    hand_factorized = GaussletBases._nested_extract_factorized_basis(factorable, dims)
    @test hand_factorized.reconstruction_max_error < 1.0e-12
    @test hand_factorized.basis_triplets == [(1, 1, 1), (1, 2, 1), (2, 1, 2)]
    @test hand_factorized.basis_amplitudes ≈ amplitudes atol = 1.0e-12 rtol = 1.0e-12
    @test hand_factorized.x_functions[:, 1] ≈ x1 atol = 1.0e-12 rtol = 1.0e-12
    @test hand_factorized.x_functions[:, 2] ≈ x2 atol = 1.0e-12 rtol = 1.0e-12
    @test hand_factorized.y_functions[:, 1] ≈ y1 atol = 1.0e-12 rtol = 1.0e-12
    @test hand_factorized.y_functions[:, 2] ≈ y2 atol = 1.0e-12 rtol = 1.0e-12
    @test hand_factorized.z_functions[:, 1] ≈ z1 atol = 1.0e-12 rtol = 1.0e-12
    @test hand_factorized.z_functions[:, 2] ≈ z2 atol = 1.0e-12 rtol = 1.0e-12
    @test GaussletBases._nested_reconstruct_factorized_coefficients(hand_factorized) ≈
          factorable atol = 1.0e-12 rtol = 1.0e-12

    broken_factorable = copy(factorable)
    broken_factorable[GaussletBases._cartesian_flat_index(2, 2, 2, dims), 2] += 1.0e-4
    @test_throws ArgumentError GaussletBases._nested_extract_factorized_basis(
        broken_factorable,
        dims,
    )

    full_reference = GaussletBases._build_one_center_atomic_shell_sequence(
        bundle,
        (1:13, 1:13, 1:13);
        nside = 5,
        packet_kernel = :support_reference,
        term_coefficients = term_coefficients,
    )
    full_direct = GaussletBases._build_one_center_atomic_shell_sequence(
        bundle,
        (1:13, 1:13, 1:13);
        nside = 5,
        packet_kernel = :factorized_direct,
        term_coefficients = term_coefficients,
    )
    legacy_reference = GaussletBases._build_one_center_atomic_shell_sequence(
        bundle,
        (2:12, 2:12, 2:12);
        nside = 5,
        packet_kernel = :support_reference,
        term_coefficients = term_coefficients,
    )
    legacy_direct = GaussletBases._build_one_center_atomic_shell_sequence(
        bundle,
        (2:12, 2:12, 2:12);
        nside = 5,
        packet_kernel = :factorized_direct,
        term_coefficients = term_coefficients,
    )

    fixed_full_reference = GaussletBases._nested_fixed_block(full_reference, bundle)
    fixed_full_direct = GaussletBases._nested_fixed_block(full_direct, bundle)
    fixed_legacy_reference = GaussletBases._nested_fixed_block(legacy_reference, bundle)
    fixed_legacy_direct = GaussletBases._nested_fixed_block(legacy_direct, bundle)

    carried_legacy = fixed_legacy_reference.factorized_cartesian_parent_basis[]
    @test !isnothing(carried_legacy)
    extracted_legacy = GaussletBases._nested_extract_factorized_basis(
        fixed_legacy_reference.coefficient_matrix,
        (length(basis), length(basis), length(basis)),
    )
    @test GaussletBases._nested_factorized_parent_basis(fixed_legacy_reference) === carried_legacy
    @test carried_legacy.basis_triplets == extracted_legacy.basis_triplets
    @test carried_legacy.basis_amplitudes ≈
          extracted_legacy.basis_amplitudes atol = 1.0e-12 rtol = 1.0e-12
    @test GaussletBases._nested_reconstruct_factorized_coefficients(carried_legacy) ≈
          fixed_legacy_reference.coefficient_matrix atol = 1.0e-10 rtol = 1.0e-10
    fixed_legacy_reference.factorized_cartesian_parent_basis[] = nothing
    lazy_legacy = GaussletBases._nested_factorized_parent_basis(fixed_legacy_reference)
    @test fixed_legacy_reference.factorized_cartesian_parent_basis[] === lazy_legacy
    @test lazy_legacy.basis_triplets == carried_legacy.basis_triplets
    @test lazy_legacy.basis_amplitudes ≈
          carried_legacy.basis_amplitudes atol = 1.0e-12 rtol = 1.0e-12

    @test fixed_full_direct.overlap ≈ fixed_full_reference.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.kinetic ≈ fixed_full_reference.kinetic atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.position_x ≈ fixed_full_reference.position_x atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.position_y ≈ fixed_full_reference.position_y atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.position_z ≈ fixed_full_reference.position_z atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.x2_x ≈ fixed_full_reference.x2_x atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.x2_y ≈ fixed_full_reference.x2_y atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.x2_z ≈ fixed_full_reference.x2_z atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.weights ≈ fixed_full_reference.weights atol = 1.0e-10 rtol = 1.0e-10
    @test !hasproperty(fixed_full_direct, :gaussian_terms)
    @test !hasproperty(fixed_full_direct, :pair_terms)
    @test !hasproperty(fixed_full_direct, :term_storage)
    @test !hasproperty(fixed_full_reference, :gaussian_terms)
    @test !hasproperty(fixed_full_reference, :pair_terms)
    @test !hasproperty(fixed_full_reference, :term_storage)
    @test fixed_full_direct.gaussian_sum ≈ fixed_full_reference.gaussian_sum atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.pair_sum ≈ fixed_full_reference.pair_sum atol = 1.0e-10 rtol = 1.0e-10

    @test fixed_legacy_direct.overlap ≈ fixed_legacy_reference.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_legacy_direct.kinetic ≈ fixed_legacy_reference.kinetic atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_legacy_direct.position_x ≈ fixed_legacy_reference.position_x atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_legacy_direct.position_y ≈ fixed_legacy_reference.position_y atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_legacy_direct.position_z ≈ fixed_legacy_reference.position_z atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_legacy_direct.x2_x ≈ fixed_legacy_reference.x2_x atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_legacy_direct.x2_y ≈ fixed_legacy_reference.x2_y atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_legacy_direct.x2_z ≈ fixed_legacy_reference.x2_z atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_legacy_direct.weights ≈ fixed_legacy_reference.weights atol = 1.0e-10 rtol = 1.0e-10
    @test !hasproperty(fixed_legacy_direct, :gaussian_terms)
    @test !hasproperty(fixed_legacy_direct, :pair_terms)
    @test !hasproperty(fixed_legacy_direct, :term_storage)
    @test !hasproperty(fixed_legacy_reference, :gaussian_terms)
    @test !hasproperty(fixed_legacy_reference, :pair_terms)
    @test !hasproperty(fixed_legacy_reference, :term_storage)
    @test fixed_legacy_direct.gaussian_sum ≈ fixed_legacy_reference.gaussian_sum atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_legacy_direct.pair_sum ≈ fixed_legacy_reference.pair_sum atol = 1.0e-10 rtol = 1.0e-10
end

@testset "QW residual-space keep policy is near-null-only and stabilized" begin
    # Literal residual-overlap spectrum observed on the anchored one-center
    # Ne legacy-profile case:
    # parent side = 29, working box = 2:28, nside = 7, supplement lmax = 1.
    residual_overlap_eigenvalues = Float64[
        6.486197469e-08,
        3.165964397e-06,
        3.165964398e-06,
        3.165964398e-06,
        5.681904400e-06,
        1.681965647e-05,
        3.337404514e-05,
        5.805312472e-05,
        5.805312472e-05,
        5.805312472e-05,
        7.256172691e-05,
        1.406818079e-04,
        1.406818079e-04,
        1.406818079e-04,
        1.927015773e-04,
        1.927015773e-04,
        1.927015773e-04,
        1.995510583e-04,
        4.261857498e-04,
        4.261857498e-04,
        4.261857498e-04,
        5.359433116e-04,
        1.945893481e-03,
        1.945893481e-03,
        1.945893481e-03,
    ]
    gausslet_overlap = Matrix{Float64}(I, 1, 1)
    overlap_ga = zeros(Float64, 1, length(residual_overlap_eigenvalues))
    overlap_aa = Matrix(Diagonal(residual_overlap_eigenvalues))
    near_null_diagnostics = diagnose_qwrg_residual_space(
        gausslet_overlap,
        overlap_ga,
        overlap_aa;
        keep_policy = :near_null_only,
        keep_abs_tol = GaussletBases._qwrg_atomic_residual_keep_tol(),
        accept_tol = GaussletBases._qwrg_atomic_residual_accept_tol(),
    )
    near_null_data = GaussletBases._qwrg_residual_space(
        gausslet_overlap,
        overlap_ga,
        overlap_aa;
        keep_policy = :near_null_only,
        keep_abs_tol = GaussletBases._qwrg_atomic_residual_keep_tol(),
        accept_tol = GaussletBases._qwrg_atomic_residual_accept_tol(),
    )
    legacy_alias_diagnostics = diagnose_qwrg_residual_space(
        gausslet_overlap,
        overlap_ga,
        overlap_aa;
        keep_policy = :legacy_profile,
        keep_abs_tol = GaussletBases._qwrg_atomic_residual_keep_tol(),
        accept_tol = GaussletBases._qwrg_atomic_residual_accept_tol(),
    )

    @test near_null_diagnostics.keep_policy == :near_null_only
    @test near_null_diagnostics.gaussian_count == 25
    @test near_null_diagnostics.supplement_numerical_rank == 25
    @test near_null_diagnostics.residual_numerical_rank == 25
    @test near_null_diagnostics.kept_count == 24
    @test near_null_diagnostics.discarded_count == 1
    @test near_null_diagnostics.keep_tol ≈ 1.0e-7 atol = 0.0 rtol = 0.0
    @test near_null_diagnostics.accept_tol ≈ 1.0e-7 atol = 0.0 rtol = 0.0
    @test near_null_diagnostics.kept_block_stabilization_null_tol ≈ 1.0e-12 atol = 1.0e-15 rtol = 0.0
    @test near_null_diagnostics.kept_block_stabilization_correction_passes >= 1
    @test near_null_diagnostics.kept_block_stabilization_clipped_count == 0
    @test near_null_diagnostics.kept_block_stabilization_dropped_count == 0
    @test near_null_diagnostics.kept_block_pre_stabilization_overlap_error < 1.0e-12
    @test near_null_diagnostics.kept_block_post_stabilization_overlap_error < 1.0e-12
    @test near_null_diagnostics.kept_block_pre_stabilization_symmetry_defect < 1.0e-12
    @test near_null_diagnostics.kept_block_post_stabilization_symmetry_defect < 1.0e-12
    @test near_null_diagnostics.kept_block_pre_stabilization_negative_count == 0
    @test near_null_diagnostics.kept_block_post_stabilization_negative_count == 0
    @test near_null_diagnostics.kept_block_pre_stabilization_near_null_count == 0
    @test near_null_diagnostics.kept_block_post_stabilization_near_null_count == 0
    @test norm(near_null_data.final_overlap - I, Inf) < 1.0e-10
    @test legacy_alias_diagnostics.keep_policy == :near_null_only
    @test legacy_alias_diagnostics.kept_count == near_null_diagnostics.kept_count
    @test legacy_alias_diagnostics.keep_tol == near_null_diagnostics.keep_tol
    @test legacy_alias_diagnostics.kept_block_post_stabilization_overlap_error ==
        near_null_diagnostics.kept_block_post_stabilization_overlap_error

    nsynthetic = 69
    synthetic_raw_overlap = Matrix{Float64}(I, nsynthetic, nsynthetic)
    synthetic_coefficients = Matrix{Float64}(I, nsynthetic, nsynthetic)
    @inbounds for i in 1:nsynthetic, j in 1:nsynthetic
        synthetic_coefficients[i, j] += 8.0e-9 * sin(Float64(i + 2 * j))
    end
    synthetic_stabilization = GaussletBases._qwrg_stabilize_residual_coefficients(
        synthetic_raw_overlap,
        synthetic_coefficients,
    )
    @test synthetic_stabilization.pre_error > 1.0e-8
    @test synthetic_stabilization.post_error < 1.0e-10
    @test synthetic_stabilization.post_symmetry_defect < 1.0e-12
    @test synthetic_stabilization.pre_negative_count == 0
    @test synthetic_stabilization.post_negative_count == 0
    @test synthetic_stabilization.dropped_count == 0
    @test synthetic_stabilization.correction_passes >= 1
end

@testset "One-center atomic legacy-profile residual completion contract" begin
    if !_RUN_SLOW_TESTS
        @test true
    else
        data = _one_center_atomic_legacy_profile_ne_residual_completion_fixture()

        @test data.fixed_gausslet_count == 2523
        @test data.supplement_count == 25

        @test data.near_null.keep_policy == :near_null_only
        @test data.near_null.residual_numerical_rank == 25
        @test data.near_null.kept_count == 24
        @test data.near_null.discarded_count == 1
        @test data.near_null.keep_tol ≈ 1.0e-7 atol = 0.0 rtol = 0.0
        @test data.near_null.accept_tol ≈ 1.0e-7 atol = 0.0 rtol = 0.0
        @test data.near_null.kept_block_pre_stabilization_overlap_error > 0.0
        @test data.near_null.kept_block_post_stabilization_overlap_error <
            data.near_null.kept_block_pre_stabilization_overlap_error
        @test data.near_null.kept_block_post_stabilization_overlap_error < 1.0e-9
        @test data.near_null.kept_block_post_stabilization_symmetry_defect < 1.0e-9
        @test data.near_null.kept_block_stabilization_dropped_count == 0
        @test norm(data.near_null_data.final_overlap - I, Inf) < 1.0e-7
        @test data.near_null_total_basis == 2547
        @test data.legacy_alias.keep_policy == :near_null_only
        @test data.legacy_alias.kept_count == data.near_null.kept_count
    end
end

@testset "Atomic residual keep policy rejects relative_case_scale on public QW routes" begin
    if !_legacy_basisfile_available()
        @test true
    else
        source_basis_qw, _legacy_qw, _ordinary_l0, _ordinary_l0_check = _qiu_white_full_nearest_fixture()
        supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 1)
        err = try
            ordinary_cartesian_qiu_white_operators(
                source_basis_qw,
                supplement;
                expansion = coulomb_gaussian_expansion(doacc = false),
                Z = 2.0,
                interaction_treatment = :ggt_nearest,
                residual_keep_policy = :relative_case_scale,
            )
            nothing
        catch caught
            caught
        end
        @test err isa ArgumentError
        @test occursin(":near_null_only", sprint(showerror, err))
        @test occursin(":legacy_profile", sprint(showerror, err))
    end
end

@testset "One-center atomic ns=9 legacy-profile residual stabilization closes center-extraction failure" begin
    if !_RUN_SLOW_TESTS
        @test true
    else
        data = _one_center_atomic_ns9_legacy_profile_qw_fixture()
        @test data.residual_data.diagnostics.kept_count == 24
        @test data.residual_data.diagnostics.keep_policy == :near_null_only
        @test data.residual_data.diagnostics.keep_tol ≈ 1.0e-7 atol = 0.0 rtol = 0.0
        @test data.residual_data.diagnostics.accept_tol ≈ 1.0e-7 atol = 0.0 rtol = 0.0
        @test data.residual_data.diagnostics.kept_block_pre_stabilization_overlap_error > 0.0
        @test data.residual_data.diagnostics.kept_block_post_stabilization_overlap_error <=
            data.residual_data.diagnostics.kept_block_pre_stabilization_overlap_error
        @test data.residual_data.diagnostics.kept_block_post_stabilization_overlap_error < 1.0e-9
        @test data.residual_data.diagnostics.kept_block_post_stabilization_symmetry_defect < 1.0e-9
        @test data.residual_data.diagnostics.kept_block_stabilization_dropped_count == 0
        @test norm(data.residual_data.final_overlap - I, Inf) < 1.0e-7
        @test data.operators.residual_count == 24
        @test norm(data.operators.overlap - I, Inf) < 1.0e-7
        check = GaussletBases.ordinary_cartesian_1s2_check(data.operators)
        @test isfinite(check.orbital_energy)
        @test check.overlap_error < 1.0e-7
    end
end

@testset "Cartesian nested shell sequence fixed-block" begin
    (
        basis,
        bundle,
        shell1,
        shell2,
        shell_plus_core,
        shell_sequence,
        fixed_shell_plus_core,
        fixed_sequence,
        legacy,
        baseline,
        shell_plus_core_ops,
        shell_sequence_ops,
        baseline_check,
        shell_plus_core_check,
        shell_sequence_check,
    ) = _nested_qiu_white_shell_sequence_fixture()

    @test shell_sequence isa GaussletBases._CartesianNestedShellSequence3D
    @test length(shell_sequence.shell_layers) == 2
    @test shell_sequence.shell_layers[1] === shell1
    @test shell_sequence.shell_layers[2] === shell2
    @test first(shell_sequence.core_column_range) == 1
    @test last(shell_sequence.core_column_range) == length(shell_sequence.core_indices)
    @test length(shell_sequence.layer_column_ranges) == 2
    @test length(shell_sequence.core_indices) == 11^3
    @test isempty(intersect(shell_sequence.core_indices, shell1.support_indices))
    @test isempty(intersect(shell_sequence.core_indices, shell2.support_indices))
    @test isempty(intersect(shell1.support_indices, shell2.support_indices))

    @test fixed_sequence isa GaussletBases._NestedFixedBlock3D
    @test fixed_sequence.parent_basis === basis
    @test fixed_sequence.shell === shell_sequence
    @test shell_plus_core_ops.gausslet_count == 1385
    @test shell_sequence_ops.gausslet_count == 1439
    @test baseline.gausslet_count == 17^3
    @test norm(fixed_shell_plus_core.overlap - I, Inf) < 1.0e-10
    @test norm(fixed_sequence.overlap - I, Inf) < 1.0e-10
    @test shell_plus_core_check.overlap_error < 1.0e-10
    @test shell_sequence_check.overlap_error < 1.0e-10
    @test shell_plus_core_check.orbital_energy < 0.0
    @test shell_sequence_check.orbital_energy < 0.0
    @test shell_plus_core_check.vee_expectation > 0.0
    @test shell_sequence_check.vee_expectation > 0.0
    @test abs(shell_plus_core_check.orbital_energy - baseline_check.orbital_energy) < 1.0e-4
    @test abs(shell_plus_core_check.vee_expectation - baseline_check.vee_expectation) < 1.0e-4
    @test abs(shell_sequence_check.orbital_energy - baseline_check.orbital_energy) < 1.0e-4
    @test abs(shell_sequence_check.vee_expectation - baseline_check.vee_expectation) < 1.0e-4
    @test abs(shell_sequence_check.orbital_energy - shell_plus_core_check.orbital_energy) < 1.0e-4
    @test abs(shell_sequence_check.vee_expectation - shell_plus_core_check.vee_expectation) < 1.0e-4
end

@testset "Cartesian nested fixed-nside compression policy" begin
    (
        basis,
        bundle,
        shell1,
        shell2,
        shell3,
        grow_sequence,
        shrinking_sequence,
        fixed_grow,
        fixed_shrink,
        legacy,
        baseline,
        grow_ops,
        shrink_ops,
        baseline_check,
        grow_check,
        shrink_check,
    ) = _nested_qiu_white_nside_sequence_fixture()

    @test GaussletBases._nested_shrunk_interval(4:14, 0; nside = 5) == 4:14
    @test GaussletBases._nested_shrunk_interval(4:14, 1; nside = 5) == 5:13
    @test GaussletBases._nested_shrunk_interval(4:14, 2; nside = 5) == 6:12
    @test GaussletBases._nested_shrunk_interval(4:14, 3; nside = 5) == 7:11
    @test GaussletBases._nested_shrunk_interval(4:14, 4; nside = 5) == 7:11

    @test shrinking_sequence isa GaussletBases._CartesianNestedShellSequence3D
    @test length(shrinking_sequence.shell_layers) == 3
    @test shrinking_sequence.shell_layers[1] === shell1
    @test shrinking_sequence.shell_layers[2] === shell2
    @test shrinking_sequence.shell_layers[3] === shell3
    @test length(grow_sequence.core_indices) == 5^3
    @test length(shrinking_sequence.core_indices) == 5^3
    @test first(shrinking_sequence.core_column_range) == 1
    @test last(shrinking_sequence.core_column_range) == 3^3
    @test isempty(intersect(shrinking_sequence.core_indices, shell1.support_indices))
    @test isempty(intersect(shrinking_sequence.core_indices, shell2.support_indices))
    @test isempty(intersect(shrinking_sequence.core_indices, shell3.support_indices))
    @test isempty(intersect(shell1.support_indices, shell2.support_indices))
    @test isempty(intersect(shell1.support_indices, shell3.support_indices))
    @test isempty(intersect(shell2.support_indices, shell3.support_indices))

    @test fixed_grow isa GaussletBases._NestedFixedBlock3D
    @test fixed_shrink isa GaussletBases._NestedFixedBlock3D
    @test grow_ops.gausslet_count == 287
    @test shrink_ops.gausslet_count == 189
    @test baseline.gausslet_count == 17^3
    @test norm(fixed_grow.overlap - I, Inf) < 1.0e-10
    @test norm(fixed_shrink.overlap - I, Inf) < 1.0e-10
    @test grow_check.overlap_error < 1.0e-10
    @test shrink_check.overlap_error < 1.0e-10
    @test isfinite(grow_check.orbital_energy)
    @test isfinite(grow_check.vee_expectation)
    @test isfinite(shrink_check.orbital_energy)
    @test isfinite(shrink_check.vee_expectation)
    @test grow_check.orbital_energy < 0.0
    @test shrink_check.orbital_energy < 0.0
    @test grow_check.vee_expectation > 0.0
    @test shrink_check.vee_expectation > 0.0
    @test shrink_ops.gausslet_count < grow_ops.gausslet_count
end

@testset "Cartesian nested complete shell layer" begin
    (
        basis,
        bundle,
        shell1_complete,
        shell2_complete,
        shell3_complete,
        shell4_complete,
        interval1,
        interval2,
        interval3,
        interval4,
        core5,
        complete_sequence,
        fixed_complete_sequence,
        legacy,
        baseline,
        complete_sequence_ops,
        baseline_check,
        complete_sequence_check,
    ) = _nested_qiu_white_complete_shell_sequence_fixture()

    @test shell1_complete isa GaussletBases._CartesianNestedCompleteShell3D
    @test shell2_complete isa GaussletBases._CartesianNestedCompleteShell3D
    @test shell3_complete isa GaussletBases._CartesianNestedCompleteShell3D
    @test shell4_complete isa GaussletBases._CartesianNestedCompleteShell3D
    @test length(shell1_complete.faces) == 6
    @test length(shell1_complete.edges) == 12
    @test length(shell1_complete.corners) == 8
    @test length(shell2_complete.faces) == 6
    @test length(shell2_complete.edges) == 12
    @test length(shell2_complete.corners) == 8
    @test length(shell3_complete.faces) == 6
    @test length(shell3_complete.edges) == 12
    @test length(shell3_complete.corners) == 8
    @test length(shell4_complete.faces) == 6
    @test length(shell4_complete.edges) == 12
    @test length(shell4_complete.corners) == 8

    @test length(shell1_complete.support_indices) == 13^3 - 11^3
    @test length(shell2_complete.support_indices) == 11^3 - 9^3
    @test length(shell3_complete.support_indices) == 9^3 - 7^3
    @test length(shell4_complete.support_indices) == 7^3 - 5^3
    @test shell1_complete.provenance.source_box == (
        (first(interval1) - 1):(last(interval1) + 1),
        (first(interval1) - 1):(last(interval1) + 1),
        (first(interval1) - 1):(last(interval1) + 1),
    )
    @test shell1_complete.provenance.next_inner_box == (interval1, interval1, interval1)
    @test shell1_complete.provenance.source_point_count == 13^3 - 11^3
    @test shell1_complete.provenance.retained_fixed_count == size(shell1_complete.coefficient_matrix, 2)
    @test shell2_complete.provenance.source_box == (
        (first(interval2) - 1):(last(interval2) + 1),
        (first(interval2) - 1):(last(interval2) + 1),
        (first(interval2) - 1):(last(interval2) + 1),
    )
    @test shell2_complete.provenance.next_inner_box == (interval2, interval2, interval2)
    @test shell2_complete.provenance.source_point_count == 11^3 - 9^3
    @test shell2_complete.provenance.retained_fixed_count == size(shell2_complete.coefficient_matrix, 2)
    @test shell3_complete.provenance.source_box == (
        (first(interval3) - 1):(last(interval3) + 1),
        (first(interval3) - 1):(last(interval3) + 1),
        (first(interval3) - 1):(last(interval3) + 1),
    )
    @test shell3_complete.provenance.next_inner_box == (interval3, interval3, interval3)
    @test shell3_complete.provenance.source_point_count == 9^3 - 7^3
    @test shell3_complete.provenance.retained_fixed_count == size(shell3_complete.coefficient_matrix, 2)
    @test shell4_complete.provenance.source_box == (
        (first(interval4) - 1):(last(interval4) + 1),
        (first(interval4) - 1):(last(interval4) + 1),
        (first(interval4) - 1):(last(interval4) + 1),
    )
    @test shell4_complete.provenance.next_inner_box == (interval4, interval4, interval4)
    @test shell4_complete.provenance.source_point_count == 7^3 - 5^3
    @test shell4_complete.provenance.retained_fixed_count == size(shell4_complete.coefficient_matrix, 2)
    @test sum(length(face.support_indices) for face in shell1_complete.faces) == 6 * 11^2
    @test sum(length(edge.support_indices) for edge in shell1_complete.edges) == 12 * 11
    @test sum(length(corner.support_indices) for corner in shell1_complete.corners) == 8

    @test complete_sequence isa GaussletBases._CartesianNestedShellSequence3D
    @test length(complete_sequence.core_indices) == 5^3
    @test complete_sequence.working_box == (3:15, 3:15, 3:15)
    @test baseline.gausslet_count == 17^3
    @test norm(fixed_complete_sequence.overlap - I, Inf) < 1.0e-10
    @test all(isfinite, fixed_complete_sequence.weights)
    @test minimum(fixed_complete_sequence.weights) > 0.0
    @test complete_sequence_check.overlap_error < 1.0e-10
    @test isfinite(complete_sequence_check.orbital_energy)
    @test isfinite(complete_sequence_check.vee_expectation)
    # Post-hardening residual-space route check for the complete-shell candidate only.
    @test abs(complete_sequence_check.vee_expectation - baseline_check.vee_expectation) < 3.0e-4

    expansion = coulomb_gaussian_expansion(doacc = false)
    overlap_parent, one_body_parent, interaction_parent = _nested_parent_fixed_problem(bundle, expansion; Z = 2.0)
    parent_modes = eigen(Hermitian(one_body_parent), Hermitian(overlap_parent))
    parent_ground = parent_modes.vectors[:, 1]
    parent_ground_vee = _nested_vee_from_orbital(interaction_parent, parent_ground)
    projected_complete = _nested_fixed_projected_orbital(overlap_parent, fixed_complete_sequence, parent_ground)
    projected_complete_vee = _nested_vee_from_orbital(
        GaussletBases._qwrg_fixed_block_interaction_matrix(fixed_complete_sequence, expansion),
        projected_complete,
    )
    term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
    @test abs(projected_complete_vee - parent_ground_vee) < 5.0e-4
    @test_throws ArgumentError GaussletBases._nested_shell_sequence(
        bundle,
        core5,
        core5,
        core5,
        [shell1_complete, shell2_complete, shell3_complete],
        term_coefficients = term_coefficients,
    )
end
