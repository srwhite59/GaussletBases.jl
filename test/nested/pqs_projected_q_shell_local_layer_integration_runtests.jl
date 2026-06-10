# Integration/slow test. Do not include in default nested runner.

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
    # Legacy smoke only: CPB provider tests own detailed one-body product math.
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
    @test product_source_box_shadow.diagnostics.output_finite
    for term in product_source_box_shadow.terms
        block = product_source_box_shadow.blocks[term]
        @test all(isfinite, block)
        @test size(block) == (8, 8)
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
