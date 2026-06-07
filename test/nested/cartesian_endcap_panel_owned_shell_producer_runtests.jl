# Integration/slow test. Do not include in default nested runner.

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
