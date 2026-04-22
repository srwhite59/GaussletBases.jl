@testset "Uniform basis" begin
    ub = build_basis(UniformBasisSpec(:G10; xmin = -2.0, xmax = 2.0, spacing = 1.0))
    primitive_data = primitives(ub)
    coefficient_matrix = stencil_matrix(ub)
    x = 0.2

    @test ub isa UniformBasis
    @test length(ub) == 5
    @test ub[2] isa Gausslet
    @test centers(ub) == [center(ub[i]) for i in 1:length(ub)]
    @test reference_centers(ub) == centers(ub)
    @test length(primitive_data) == size(coefficient_matrix, 1)
    @test size(coefficient_matrix, 2) == length(ub)
    @test sum(coefficient_matrix[mu, 3] * primitive_data[mu](x) for mu in eachindex(primitive_data)) ≈
          ub[3](x) atol = 1.0e-12 rtol = 1.0e-12
end

@testset "Gausslet construction and evaluation" begin
    g = Gausslet(:G10; center = 0.0, spacing = 1.0)
    x = 0.2

    @test g isa Gausslet
    @test g(x) == value(g, x)
    @test center(g) == 0.0

    st = stencil(g)
    explicit_sum = sum(coefficients(st)[i] * primitives(st)[i](x) for i in eachindex(coefficients(st)))

    @test direct_value(g, x) == explicit_sum
    @test st(x) == explicit_sum
    @test g(x) == explicit_sum
    @test length(coefficients(st)) == length(primitives(st))
    @test integral_weight(g) ≈ 1.0 atol = 1.0e-12
end

@testset "Stencil consistency" begin
    g = Gausslet(:G10; center = 0.0, spacing = 1.0)
    st = stencil(g)

    @test coefficients(st)[1] == coefficients(st)[end]
    @test center(primitives(st)[div(length(primitives(st)) + 1, 2)]) == 0.0
end

@testset "Coordinate mappings" begin
    map = AsinhMapping(c = 0.15, s = 0.15)
    x = 3.0

    @test map(x) == uofx(map, x)
    @test xofu(map, uofx(map, x)) ≈ x atol = 1.0e-12 rtol = 1.0e-12
    @test dudx(map, x) > 0.0
    @test uofx(map, x) > asinh(x / 1.0) / 0.15
end

@testset "AsinhMapping constructor semantics" begin
    c0 = 0.15
    s0 = 0.15
    t0 = 10.0
    map_from_c = AsinhMapping(c = c0, s = s0, tail_spacing = t0)
    map_from_a = AsinhMapping(a = c0 / s0, s = s0, tail_spacing = t0)

    for x in (-2.0, -0.5, 0.0, 0.75, 3.0)
        @test uofx(map_from_c, x) ≈ uofx(map_from_a, x) atol = 1.0e-12 rtol = 1.0e-12
        @test dudx(map_from_c, x) ≈ dudx(map_from_a, x) atol = 1.0e-12 rtol = 1.0e-12
        @test du2dx2(map_from_c, x) ≈ du2dx2(map_from_a, x) atol = 1.0e-12 rtol = 1.0e-12
    end

    @test xofu(map_from_c, 1.25) ≈ xofu(map_from_a, 1.25) atol = 1.0e-12 rtol = 1.0e-12
    @test_throws ArgumentError AsinhMapping(c = c0, a = c0 / s0, s = s0)
    @test_throws ArgumentError AsinhMapping(s = s0)
end

@testset "AsinhMapping keeps the linear tail term" begin
    a0 = 1.0
    s0 = 0.15
    t0 = 10.0
    x = 3.0
    map = AsinhMapping(a = a0, s = s0, tail_spacing = t0)

    @test uofx(map, x) ≈ x / t0 + asinh(x / a0) / s0 atol = 1.0e-12 rtol = 1.0e-12
    @test dudx(map, x) ≈ 1.0 / t0 + 1.0 / (s0 * sqrt(x * x + a0 * a0)) atol = 1.0e-12 rtol = 1.0e-12
end

@testset "White-Lindsey atomic mapping matches the legacy one-center formula" begin
    function legacy_atomic_wl_u(x::Float64, Z::Float64, d::Float64, wi::Float64)
        a = sqrt(d / Z)
        s = sqrt(d * Z)
        return x / wi + asinh(x / a) / s
    end

    function legacy_atomic_wl_dudx(x::Float64, Z::Float64, d::Float64, wi::Float64)
        a = sqrt(d / Z)
        s = sqrt(d * Z)
        return 1.0 / wi + 1.0 / (s * sqrt(x * x + a * a))
    end

    for (Z, d, wi, expected_s) in ((10.0, 0.02, 6.0, sqrt(0.2)), (10.0, 0.03, 6.0, sqrt(0.3)), (2.0, 0.2, 10.0, sqrt(0.4)))
        map = white_lindsey_atomic_mapping(Z = Z, d = d, tail_spacing = wi)

        @test map isa AsinhMapping
        @test map.a ≈ sqrt(d / Z) atol = 1.0e-14 rtol = 0.0
        @test map.s ≈ expected_s atol = 1.0e-14 rtol = 0.0
        @test map.a * map.s ≈ d atol = 1.0e-14 rtol = 0.0
        @test map.tail_spacing ≈ wi atol = 1.0e-14 rtol = 0.0

        for x in (-12.0, -1.25, 0.0, 2.5, 11.0)
            xval = Float64(x)
            @test uofx(map, xval) ≈ legacy_atomic_wl_u(xval, Z, d, wi) atol = 1.0e-12 rtol = 1.0e-12
            @test dudx(map, xval) ≈ legacy_atomic_wl_dudx(xval, Z, d, wi) atol = 1.0e-12 rtol = 1.0e-12
        end
    end
end

@testset "CombinedInvsqrtMapping supports bond-aligned diatomic symmetry" begin
    basis = bond_aligned_homonuclear_qw_basis(
        bond_length = 1.4,
        core_spacing = 0.5,
        xmax_parallel = 6.0,
        xmax_transverse = 4.0,
        bond_axis = :z,
    )

    @test basis isa BondAlignedDiatomicQWBasis3D
    @test mapping(basis.basis_x) === mapping(basis.basis_y)
    @test mapping(basis.basis_z) isa CombinedInvsqrtMapping
    @test centers(basis.basis_z) ≈ -reverse(centers(basis.basis_z)) atol = 1.0e-12 rtol = 1.0e-12
    @test xofu(mapping(basis.basis_z), uofx(mapping(basis.basis_z), 0.7)) ≈ 0.7 atol = 1.0e-10 rtol = 0.0
    @test xofu(mapping(basis.basis_z), uofx(mapping(basis.basis_z), -0.7)) ≈ -0.7 atol = 1.0e-10 rtol = 0.0
    @test 1.0 / dudx(mapping(basis.basis_z), -0.7) ≈ 0.5 atol = 1.0e-10 rtol = 0.0
    @test 1.0 / dudx(mapping(basis.basis_z), 0.7) ≈ 0.5 atol = 1.0e-10 rtol = 0.0
    @test mapping(basis.basis_x).centers == [0.0]
    @test 1.0 / dudx(mapping(basis.basis_x), 0.0) ≈ 0.5 atol = 1.0e-10 rtol = 0.0
end

@testset "CombinedInvsqrtMapping supports first bond-aligned heteronuclear rule" begin
    bond_length = 1.45
    basis = bond_aligned_heteronuclear_qw_basis(
        atoms = ("He", "H"),
        bond_length = bond_length,
        core_spacings = (0.25, 0.5),
        nuclear_charges = (2.0, 1.0),
        xmax_parallel = 6.0,
        xmax_transverse = 4.0,
        bond_axis = :z,
    )
    zmap = mapping(basis.basis_z)
    xmap = mapping(basis.basis_x)

    @test basis isa BondAlignedDiatomicQWBasis3D
    @test mapping(basis.basis_x) === mapping(basis.basis_y)
    @test zmap isa CombinedInvsqrtMapping
    @test xmap isa CombinedInvsqrtMapping
    @test 1.0 / dudx(zmap, -0.5 * bond_length) ≈ 0.25 atol = 1.0e-10 rtol = 0.0
    @test 1.0 / dudx(zmap, 0.5 * bond_length) ≈ 0.5 atol = 1.0e-10 rtol = 0.0
    @test xmap.centers == [0.0]
    @test 1.0 / dudx(xmap, 0.0) ≈ 0.25 atol = 1.0e-10 rtol = 0.0
end

@testset "CombinedInvsqrtMapping supports experimental homonuclear chain geometry" begin
    basis4 = bond_aligned_homonuclear_chain_qw_basis(
        natoms = 4,
        spacing = 1.4,
        core_spacing = 0.5,
        xmax_parallel = 4.0,
        xmax_transverse = 3.5,
        chain_axis = :z,
    )
    explicit_basis = bond_aligned_homonuclear_chain_qw_basis(
        chain_coordinates = [-1.4, 0.0, 1.4],
        core_spacing = 0.5,
        xmax_parallel = 4.0,
        xmax_transverse = 3.5,
        chain_axis = :x,
    )
    diagnostics4 = bond_aligned_homonuclear_chain_geometry_diagnostics(basis4)
    explicit_diagnostics = bond_aligned_homonuclear_chain_geometry_diagnostics(explicit_basis)

    @test basis4 isa BondAlignedHomonuclearChainQWBasis3D
    @test explicit_basis isa BondAlignedHomonuclearChainQWBasis3D
    @test mapping(basis4.basis_x) === mapping(basis4.basis_y)
    @test mapping(basis4.basis_z) isa CombinedInvsqrtMapping
    @test mapping(explicit_basis.basis_x) isa CombinedInvsqrtMapping
    @test diagnostics4.axis_monotone
    @test diagnostics4.chain_coordinates ≈ [-2.1, -0.7, 0.7, 2.1] atol = 1.0e-12 rtol = 0.0
    @test centers(basis4.basis_z) ≈ -reverse(centers(basis4.basis_z)) atol = 1.0e-12 rtol = 1.0e-12
    @test diagnostics4.axis_center_symmetry_error < 1.0e-12
    @test diagnostics4.local_spacings_at_nuclei ≈ fill(0.5, 4) atol = 1.0e-10 rtol = 0.0
    @test all(diagnostics4.local_spacings_at_midpoints .> 0.45)
    @test explicit_diagnostics.chain_coordinates ≈ [-1.4, 0.0, 1.4] atol = 1.0e-12 rtol = 0.0
    @test explicit_diagnostics.local_spacings_at_nuclei ≈ fill(0.5, 3) atol = 1.0e-10 rtol = 0.0
    @test xofu(mapping(explicit_basis.basis_x), uofx(mapping(explicit_basis.basis_x), 0.0)) ≈ 0.0 atol = 1.0e-10 rtol = 0.0
end

@testset "CombinedInvsqrtMapping supports experimental homonuclear square-lattice geometry" begin
    basis2 = axis_aligned_homonuclear_square_lattice_qw_basis(
        n = 2,
        spacing = 1.4,
        core_spacing = 0.5,
        xmax_in_plane = 3.5,
        xmax_transverse = 3.0,
    )
    basis3 = axis_aligned_homonuclear_square_lattice_qw_basis(
        n = 3,
        spacing = 1.2,
        core_spacing = 0.5,
        xmax_in_plane = 3.5,
        xmax_transverse = 3.0,
    )
    explicit_basis = axis_aligned_homonuclear_square_lattice_qw_basis(
        x_coordinates = [-1.2, 0.0, 1.2],
        y_coordinates = [-1.2, 0.0, 1.2],
        core_spacing = 0.5,
        xmax_in_plane = 3.5,
        xmax_transverse = 3.0,
    )

    diagnostics2 = axis_aligned_homonuclear_square_lattice_geometry_diagnostics(basis2)
    diagnostics3 = axis_aligned_homonuclear_square_lattice_geometry_diagnostics(basis3)
    explicit_diagnostics = axis_aligned_homonuclear_square_lattice_geometry_diagnostics(explicit_basis)

    @test basis2 isa AxisAlignedHomonuclearSquareLatticeQWBasis3D
    @test basis3 isa AxisAlignedHomonuclearSquareLatticeQWBasis3D
    @test explicit_basis isa AxisAlignedHomonuclearSquareLatticeQWBasis3D
    @test length(basis2.nuclei) == 4
    @test length(basis3.nuclei) == 9
    @test mapping(basis2.basis_x) isa CombinedInvsqrtMapping
    @test mapping(basis2.basis_y) isa CombinedInvsqrtMapping
    @test mapping(basis2.basis_z) isa CombinedInvsqrtMapping
    @test diagnostics2.x_axis_monotone
    @test diagnostics2.y_axis_monotone
    @test diagnostics2.z_axis_monotone
    @test diagnostics2.x_coordinates ≈ [-0.7, 0.7] atol = 1.0e-12 rtol = 0.0
    @test diagnostics2.y_coordinates ≈ [-0.7, 0.7] atol = 1.0e-12 rtol = 0.0
    @test diagnostics3.x_coordinates ≈ [-1.2, 0.0, 1.2] atol = 1.0e-12 rtol = 0.0
    @test diagnostics3.y_coordinates ≈ [-1.2, 0.0, 1.2] atol = 1.0e-12 rtol = 0.0
    @test diagnostics2.x_axis_center_symmetry_error < 1.0e-12
    @test diagnostics2.y_axis_center_symmetry_error < 1.0e-12
    @test diagnostics3.x_axis_center_symmetry_error < 1.0e-12
    @test diagnostics3.y_axis_center_symmetry_error < 1.0e-12
    @test diagnostics2.xy_axis_center_match_error < 1.0e-12
    @test diagnostics3.xy_axis_center_match_error < 1.0e-12
    @test diagnostics2.local_spacings_at_x_coordinates ≈ fill(0.5, 2) atol = 1.0e-10 rtol = 0.0
    @test diagnostics2.local_spacings_at_y_coordinates ≈ fill(0.5, 2) atol = 1.0e-10 rtol = 0.0
    @test diagnostics3.local_spacings_at_x_coordinates ≈ fill(0.5, 3) atol = 1.0e-10 rtol = 0.0
    @test diagnostics3.local_spacings_at_y_coordinates ≈ fill(0.5, 3) atol = 1.0e-10 rtol = 0.0
    @test diagnostics2.local_spacing_at_plane_center_x > 0.45
    @test diagnostics2.local_spacing_at_plane_center_y > 0.45
    @test diagnostics3.local_spacing_at_plane_center_x ≥ 0.5 - 1.0e-10
    @test diagnostics3.local_spacing_at_plane_center_y ≥ 0.5 - 1.0e-10
    @test all(diagnostics3.representative_midpoint_spacings_x .> 0.45)
    @test all(diagnostics3.representative_midpoint_spacings_y .> 0.45)
    @test diagnostics2.local_spacing_at_plane_center_z ≈ 0.5 atol = 1.0e-10 rtol = 0.0
    @test explicit_diagnostics.x_coordinates ≈ [-1.2, 0.0, 1.2] atol = 1.0e-12 rtol = 0.0
    @test explicit_diagnostics.y_coordinates ≈ [-1.2, 0.0, 1.2] atol = 1.0e-12 rtol = 0.0
    @test explicit_diagnostics.xy_axis_center_match_error < 1.0e-12
end

@testset "XGaussian center" begin
    g = XGaussian(alpha = 0.23)
    @test center(g) == 0.23
end

if _RUN_SLOW_TESTS
    @testset "Half-line basis" begin
        spec = HalfLineBasisSpec(:G10;
            xmax = 2.0,
            reference_spacing = 1.0,
            tails = 2,
            mapping = AsinhMapping(a = 1.0, s = 0.25),
        )
        hb = build_basis(spec)
        primitive_data = primitives(hb)
        coefficient_matrix = stencil_matrix(hb)
        st = stencil(hb[2])

        @test hb isa HalfLineBasis
        @test length(hb) >= 3
        @test hb[1] isa BoundaryGausslet
        @test hb[1](0.2) == value(hb[1], 0.2)
        @test hb[1](-0.5) ≈ 0.0 atol = 1.0e-12
        @test sum(coefficient_matrix[mu, 2] * primitive_data[mu](0.3) for mu in eachindex(primitive_data)) ≈
              hb[2](0.3) atol = 1.0e-12 rtol = 1.0e-12
        @test st(0.3) ≈ direct_value(hb[2], 0.3) atol = 1.0e-12 rtol = 1.0e-12
        @test length(coefficients(st)) == length(primitive_data)
        @test collect(coefficients(st)) ≈ coefficient_matrix[:, 2] atol = 1.0e-12 rtol = 1.0e-12
        @test all(primitives(st)[i] === primitive_data[i] for i in eachindex(primitive_data))
        @test all(primitive -> primitive isa Distorted && primitive.primitive isa HalfLineGaussian, primitive_data)
        @test center(hb[2]) ≈ xofu(mapping(hb), reference_center(hb[2])) atol = 1.0e-12 rtol = 1.0e-12
        @test issorted(reference_centers(hb))
        @test issorted(centers(hb))
    end
end

@testset "Radial basis with count" begin
    spec = RadialBasisSpec(:G10;
        count = 6,
        mapping = AsinhMapping(a = 1.0, s = 0.2),
        reference_spacing = 1.0,
        tails = 3,
        odd_even_kmax = 2,
        xgaussians = [XGaussian(alpha = 0.2)],
    )
    rb = build_basis(spec)
    primitive_data = primitives(rb)
    coefficient_matrix = stencil_matrix(rb)
    st = stencil(rb[2])

    @test rb isa RadialBasis
    @test length(rb) == 6
    @test rb[1] isa RadialGausslet
    @test abs(rb[1](0.0)) ≤ 1.0e-8
    @test issorted(reference_centers(rb))
    @test issorted(centers(rb))
    @test sum(coefficient_matrix[mu, 2] * primitive_data[mu](0.3) for mu in eachindex(primitive_data)) ≈
          rb[2](0.3) atol = 1.0e-12 rtol = 1.0e-12
    @test st(0.3) ≈ direct_value(rb[2], 0.3) atol = 1.0e-12 rtol = 1.0e-12
    @test length(coefficients(st)) == length(primitive_data)
    @test collect(coefficients(st)) ≈ coefficient_matrix[:, 2] atol = 1.0e-12 rtol = 1.0e-12
    @test all(primitives(st)[i] === primitive_data[i] for i in eachindex(primitive_data))
    @test any(primitive -> primitive isa Distorted && primitive.primitive isa XGaussian, primitive_data)
end

if _RUN_SLOW_TESTS
    @testset "Radial basis with rmax" begin
        spec = RadialBasisSpec(:G10;
            rmax = 5.0,
            mapping = AsinhMapping(c = 0.15, s = 0.15),
            reference_spacing = 1.0,
            tails = 3,
            odd_even_kmax = 2,
            xgaussians = XGaussian[],
        )
        rb = build_basis(spec)

        @test rb isa RadialBasis
        @test length(rb) >= 3
    end
end

if _RUN_SLOW_TESTS
    @testset "Construction grid controls" begin
        rspec = RadialBasisSpec(:G10;
            count = 4,
            mapping = AsinhMapping(c = 0.15, s = 0.15),
            reference_spacing = 1.0,
            tails = 2,
            odd_even_kmax = 1,
            xgaussians = XGaussian[],
        )
        rb_fixed = build_basis(rspec; grid_h = 0.04, refine_grid_h = false)
        rb_refined = build_basis(rspec; grid_h = 0.04, refine_grid_h = true)
        rdata_fixed = GaussletBases._build_radial_coefficients(rspec; grid_h = 0.04)
        rdata_refined = GaussletBases._select_construction_data(
            h -> GaussletBases._build_radial_coefficients(rspec; grid_h = h),
            GaussletBases._radial_overlap_deviation,
            0.04;
            refine_grid_h = true,
        )

        @test rb_fixed isa RadialBasis
        @test rb_refined isa RadialBasis
        @test GaussletBases._radial_overlap_deviation(rdata_refined) <= GaussletBases._radial_overlap_deviation(rdata_fixed) + 1.0e-12

        hspec = HalfLineBasisSpec(:G10;
            xmax = 2.0,
            reference_spacing = 1.0,
            tails = 2,
            mapping = AsinhMapping(a = 1.0, s = 0.2),
        )
        hb_fixed = build_basis(hspec; grid_h = 0.04, refine_grid_h = false)
        hb_refined = build_basis(hspec; grid_h = 0.04, refine_grid_h = true)
        hdata_fixed = GaussletBases._build_halfline_coefficients(hspec; grid_h = 0.04)
        hdata_refined = GaussletBases._select_construction_data(
            h -> GaussletBases._build_halfline_coefficients(hspec; grid_h = h),
            GaussletBases._halfline_overlap_deviation,
            0.04;
            refine_grid_h = true,
        )

        @test hb_fixed isa HalfLineBasis
        @test hb_refined isa HalfLineBasis
        @test GaussletBases._halfline_overlap_deviation(hdata_refined) <= GaussletBases._halfline_overlap_deviation(hdata_fixed) + 1.0e-12
    end
end

@testset "Primitive contractions" begin
    rb = build_basis(RadialBasisSpec(:G10;
        count = 6,
        mapping = AsinhMapping(c = 0.15, s = 0.15),
        reference_spacing = 1.0,
        tails = 3,
        odd_even_kmax = 2,
        xgaussians = [XGaussian(alpha = 0.2)],
    ))
    coefficient_matrix = stencil_matrix(rb)
    nprimitive = size(coefficient_matrix, 1)
    vmu = collect(1.0:nprimitive)
    dmu = collect(range(1.0, step = 0.1, length = nprimitive))
    Amunu = [1.0 / (i + j) for i in 1:nprimitive, j in 1:nprimitive]

    @test contract_primitive_vector(rb, vmu) ≈ coefficient_matrix' * vmu atol = 1.0e-12 rtol = 1.0e-12
    @test contract_primitive_diagonal(rb, dmu) ≈ coefficient_matrix' * Diagonal(dmu) * coefficient_matrix atol = 1.0e-12 rtol = 1.0e-12
    @test contract_primitive_matrix(rb, Amunu) ≈ coefficient_matrix' * Amunu * coefficient_matrix atol = 1.0e-12 rtol = 1.0e-12
end

@testset "Primitive sets and metadata" begin
    g = Gausslet(:G10; center = 0.0, spacing = 1.0)
    plain_set = PrimitiveSet1D(primitives(stencil(g)); name = :plain_gausslet_stencil)
    overlap_analytic = overlap_matrix(plain_set)
    overlap_numerical = GaussletBases._primitive_overlap_matrix(
        plain_set,
        GaussletBases._NumericalPrimitiveMatrixBackend(),
    )
    position_analytic = position_matrix(plain_set)
    position_numerical = GaussletBases._primitive_position_matrix(
        plain_set,
        GaussletBases._NumericalPrimitiveMatrixBackend(),
    )
    kinetic_analytic = kinetic_matrix(plain_set)
    kinetic_numerical = GaussletBases._primitive_kinetic_matrix(
        plain_set,
        GaussletBases._NumericalPrimitiveMatrixBackend(),
    )

    map = AsinhMapping(c = 0.15, s = 0.15)
    distorted_set = PrimitiveSet1D(
        [Distorted(primitive, map) for primitive in primitives(plain_set)];
        name = :distorted_gausslet_stencil,
    )
    distorted_overlap = overlap_matrix(distorted_set)

    ub = build_basis(UniformBasisSpec(:G10; xmin = -2.0, xmax = 2.0, spacing = 1.0))
    metadata = basis_metadata(ub)

    @test length(plain_set) == length(primitives(stencil(g)))
    @test overlap_analytic ≈ overlap_numerical atol = 1.0e-9 rtol = 1.0e-9
    @test position_analytic ≈ position_numerical atol = 1.0e-9 rtol = 1.0e-9
    @test kinetic_analytic ≈ kinetic_numerical atol = 1.0e-9 rtol = 1.0e-9
    @test all(isfinite, distorted_overlap)
    @test distorted_overlap ≈ transpose(distorted_overlap) atol = 1.0e-10 rtol = 1.0e-10
    @test metadata.basis_kind == :uniform
    @test metadata.family_name == :G10
    @test size(metadata.coefficient_matrix, 1) == length(metadata.primitive_set)
    @test size(metadata.coefficient_matrix, 2) == length(ub)
    @test metadata.center_data == centers(ub)
end

@testset "Basis contraction from primitive layer" begin
    ub = build_basis(UniformBasisSpec(:G10; xmin = -1.0, xmax = 1.0, spacing = 1.0))
    P = primitive_set(ub)
    Smu = overlap_matrix(P)
    Xmu = position_matrix(P)
    Tmu = kinetic_matrix(P)
    Sb = contract_primitive_matrix(ub, Smu)
    Xb = contract_primitive_matrix(ub, Xmu)
    Tb = contract_primitive_matrix(ub, Tmu)
    Sref, Tref = _midpoint_reference_matrices(ub)
    Xref = _midpoint_reference_position_matrix(ub)

    @test Sb ≈ Sref atol = 1.0e-12 rtol = 1.0e-12
    @test Xb ≈ Xref atol = 1.0e-12 rtol = 1.0e-12
    @test Tb ≈ Tref atol = 1.0e-12 rtol = 1.0e-12
end

@testset "Ordinary Coulomb expansion and Gaussian factors" begin
    expansion = coulomb_gaussian_expansion()
    sample_points = [1.0e-3, 1.0e-2, 0.1, 1.0, 5.0, 20.0]

    @test expansion isa CoulombGaussianExpansion
    @test length(expansion) == length(expansion.coefficients) == length(expansion.exponents)
    @test all(expansion.exponents .> 0.0)
    @test sprint(show, expansion) |> x -> occursin("CoulombGaussianExpansion", x)
    @test maximum(abs.(expansion.(sample_points) .* sample_points .- 1.0)) ≤ 1.0e-6

    ub = build_basis(UniformBasisSpec(:G10; xmin = -1.0, xmax = 1.0, spacing = 1.0))
    gaussian_basis_matrix = gaussian_factor_matrix(ub; exponent = 0.7, center = 0.25)
    gaussian_reference = _midpoint_reference_gaussian_factor_matrix(ub; exponent = 0.7, center = 0.25)
    @test gaussian_basis_matrix ≈ gaussian_reference atol = 1.0e-10 rtol = 1.0e-10

    map = AsinhMapping(c = 0.15, s = 0.15)
    distorted_set = PrimitiveSet1D(
        [
            Distorted(Gaussian(center = -0.5, width = 0.35), map),
            Distorted(Gaussian(center = 0.0, width = 0.35), map),
            Distorted(Gaussian(center = 0.5, width = 0.35), map),
        ];
        name = :distorted_gaussian_triplet,
    )
    distorted_matrix = gaussian_factor_matrix(distorted_set; exponent = 0.7, center = 0.25)
    @test distorted_matrix ≈ transpose(distorted_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test all(isfinite, distorted_matrix)
end

@testset "Mapped uniform basis" begin
    ub = build_basis(UniformBasisSpec(:G10; xmin = -2.0, xmax = 2.0, spacing = 1.0))
    identity_mapped = build_basis(MappedUniformBasisSpec(:G10;
        count = 5,
        mapping = IdentityMapping(),
        reference_spacing = 1.0,
    ))
    fitted_map = fit_asinh_mapping_for_extent(npoints = 9, xmax = 6.0)
    strength_map = fit_asinh_mapping_for_strength(s = 0.5, npoints = 5, xmax = 6.0)
    mapped = build_basis(MappedUniformBasisSpec(:G10;
        count = 5,
        mapping = fit_asinh_mapping_for_extent(npoints = 5, xmax = 6.0),
        reference_spacing = 1.0,
    ))
    mapped_representation = basis_representation(mapped; operators = (:overlap, :kinetic))
    @test length(identity_mapped) == 5
    @test reference_centers(identity_mapped) == centers(ub)
    @test centers(identity_mapped) == centers(ub)
    @test stencil_matrix(identity_mapped) ≈ stencil_matrix(ub) atol = 1.0e-12 rtol = 1.0e-12
    @test primitive_set(identity_mapped).primitive_data == primitive_set(ub).primitive_data

    @test fitted_map isa AsinhMapping
    @test xofu(fitted_map, 4.0) ≈ 6.0 atol = 1.0e-10 rtol = 0.0
    @test xofu(fitted_map, -4.0) ≈ -6.0 atol = 1.0e-10 rtol = 0.0
    @test strength_map isa AsinhMapping
    @test xofu(strength_map, 2.0) ≈ 6.0 atol = 1.0e-10 rtol = 0.0
    @test xofu(strength_map, -2.0) ≈ -6.0 atol = 1.0e-10 rtol = 0.0

    @test mapped isa MappedUniformBasis
    @test length(mapped) == 5
    @test length(primitives(mapped)) == size(stencil_matrix(mapped), 1)
    @test reference_centers(mapped) == [-2.0, -1.0, 0.0, 1.0, 2.0]
    @test first(centers(mapped)) ≈ -6.0 atol = 1.0e-10 rtol = 0.0
    @test last(centers(mapped)) ≈ 6.0 atol = 1.0e-10 rtol = 0.0
    @test mapped_representation.basis_matrices.overlap ≈
          transpose(mapped_representation.basis_matrices.overlap) atol = 1.0e-10 rtol = 1.0e-10
    @test mapped_representation.basis_matrices.kinetic ≈
          transpose(mapped_representation.basis_matrices.kinetic) atol = 1.0e-10 rtol = 1.0e-10
end

@testset "Mapped PGDG prototype" begin
    uniform = build_basis(UniformBasisSpec(:G10; xmin = -2.0, xmax = 2.0, spacing = 1.0))
    uniform_representation = basis_representation(uniform; operators = (:overlap, :kinetic))
    identity_mapped = build_basis(MappedUniformBasisSpec(:G10;
        count = 5,
        mapping = IdentityMapping(),
        reference_spacing = 1.0,
    ))
    identity_representation = basis_representation(identity_mapped; operators = (:overlap, :kinetic))
    identity_proto = mapped_pgdg_prototype(identity_mapped)
    identity_refined = GaussletBases.mapped_pgdg_logfit_prototype(identity_mapped)

    basis, mapping_value, representation, prototype, localized, refined, refined_localized, exponent, factor_numeric, factor_pgdg, factor_localized, factor_refined, factor_refined_localized =
        _mapped_pgdg_1d_fixture()

    plain_gaussian = x -> exp(-0.5 * ((x - 0.4) / 0.7)^2)
    shifted_gaussian = x -> exp(-0.5 * ((x + 1.3) / 0.5)^2)
    xgaussian_like = x -> x * exp(-0.5 * (x / 0.8)^2)
    span_pre = _subspace_overlap_metric(basis, prototype)
    span_refined = _subspace_overlap_metric(basis, refined)
    projector_pre = _projector_difference_metric(basis, prototype)
    projector_refined = _projector_difference_metric(basis, refined)
    plain_error_pre = _projection_error(prototype, plain_gaussian)
    shifted_error_pre = _projection_error(prototype, shifted_gaussian)
    xgaussian_error_pre = _projection_error(prototype, xgaussian_like)
    plain_error_refined = _projection_error(refined, plain_gaussian)
    shifted_error_refined = _projection_error(refined, shifted_gaussian)
    xgaussian_error_refined = _projection_error(refined, xgaussian_like)

    @test identity_proto isa MappedPGDGPrototype1D
    @test identity_refined isa GaussletBases.MappedPGDGLogFitPrototype1D
    @test length(identity_proto) == length(uniform)
    @test overlap_matrix(identity_proto) ≈ uniform_representation.basis_matrices.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test kinetic_matrix(identity_proto) ≈ uniform_representation.basis_matrices.kinetic atol = 1.0e-12 rtol = 1.0e-12
    @test overlap_matrix(identity_proto) ≈ identity_representation.basis_matrices.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test kinetic_matrix(identity_proto) ≈ identity_representation.basis_matrices.kinetic atol = 1.0e-12 rtol = 1.0e-12
    @test overlap_matrix(identity_refined) ≈ uniform_representation.basis_matrices.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test kinetic_matrix(identity_refined) ≈ uniform_representation.basis_matrices.kinetic atol = 1.0e-12 rtol = 1.0e-12

    @test mapping_value isa AsinhMapping
    @test prototype isa MappedPGDGPrototype1D
    @test localized isa MappedPGDGLocalized1D
    @test refined isa GaussletBases.MappedPGDGLogFitPrototype1D
    @test refined_localized isa MappedPGDGLocalized1D
    @test occursin("experimental", sprint(show, prototype))
    @test occursin("experimental", sprint(show, localized))
    @test occursin("experimental", sprint(show, refined))
    @test length(primitives(prototype)) == size(stencil_matrix(prototype), 1)
    @test overlap_matrix(prototype) ≈ transpose(overlap_matrix(prototype)) atol = 1.0e-10 rtol = 1.0e-10
    @test kinetic_matrix(prototype) ≈ transpose(kinetic_matrix(prototype)) atol = 1.0e-10 rtol = 1.0e-10
    @test factor_pgdg ≈ transpose(factor_pgdg) atol = 1.0e-10 rtol = 1.0e-10
    @test norm(representation.basis_matrices.overlap - overlap_matrix(prototype), Inf) < 0.06
    @test norm(representation.basis_matrices.kinetic - kinetic_matrix(prototype), Inf) < 0.06
    @test norm(factor_numeric - factor_pgdg, Inf) < 0.07
    @test overlap_matrix(localized) ≈ I atol = 1.0e-10 rtol = 1.0e-10
    @test position_matrix(localized) ≈ Diagonal(centers(localized)) atol = 1.0e-10 rtol = 1.0e-10
    @test factor_localized ≈ transpose(factor_localized) atol = 1.0e-10 rtol = 1.0e-10
    @test overlap_matrix(refined) ≈ transpose(overlap_matrix(refined)) atol = 1.0e-10 rtol = 1.0e-10
    @test kinetic_matrix(refined) ≈ transpose(kinetic_matrix(refined)) atol = 1.0e-10 rtol = 1.0e-10
    @test factor_refined ≈ transpose(factor_refined) atol = 1.0e-10 rtol = 1.0e-10
    @test overlap_matrix(refined_localized) ≈ I atol = 1.0e-10 rtol = 1.0e-10
    @test position_matrix(refined_localized) ≈ Diagonal(centers(refined_localized)) atol = 1.0e-10 rtol = 1.0e-10
    @test factor_refined_localized ≈ transpose(factor_refined_localized) atol = 1.0e-10 rtol = 1.0e-10
    @test span_refined.min_sv > span_pre.min_sv
    @test projector_refined.frob < projector_pre.frob
    @test projector_refined.op < projector_pre.op
    @test plain_error_refined ≤ plain_error_pre + 0.01
    @test shifted_error_refined ≤ shifted_error_pre + 0.01
    @test xgaussian_error_refined ≤ xgaussian_error_pre + 0.02
end

@testset "Ternary PGDG refinement mask" begin
    mask = GaussletBases._default_ternary_gaussian_refinement_mask()
    offsets = GaussletBases._refinement_mask_offsets(mask)
    residue_sums = GaussletBases._refinement_mask_residue_sums(mask)
    residue_ripple = GaussletBases._refinement_mask_residue_ripple(mask)

    @test mask.factor == 3
    @test mask.rho ≈ 1.2 atol = 0.0 rtol = 0.0
    @test mask.support_radius == 24
    @test mask.half_window ≈ 8.0 atol = 0.0 rtol = 0.0
    @test length(mask) == 49
    @test offsets == collect(-24:24)
    @test mask.coefficients ≈ reverse(mask.coefficients) atol = 0.0 rtol = 0.0
    @test all(>(0.0), mask.coefficients)
    @test maximum(
        abs(
            coefficient - GaussletBases._analytic_ternary_refinement_coefficient(offset, mask.rho),
        ) for (offset, coefficient) in zip(offsets, mask.coefficients)
    ) == 0.0
    @test abs(sum(mask.coefficients) - 3.0) < 1.0e-10
    @test maximum(abs.(residue_sums .- 1.0)) < 3.0e-11
    @test residue_ripple < 3.0e-11

    one_step = GaussletBases._apply_gaussian_refinement_mask(mask, [1.0])
    two_step_direct = GaussletBases._apply_gaussian_refinement_mask(mask, one_step)
    two_step_repeat = GaussletBases._apply_gaussian_refinement_mask_repeated(mask, [1.0]; levels = 2)
    three_step_repeat = GaussletBases._apply_gaussian_refinement_mask_repeated(mask, [1.0]; levels = 3)

    @test one_step.offset == -mask.support_radius
    @test one_step.coefficients ≈ mask.coefficients atol = 0.0 rtol = 0.0
    @test two_step_direct.offset == two_step_repeat.offset
    @test two_step_direct.coefficients ≈ two_step_repeat.coefficients atol = 1.0e-14 rtol = 1.0e-14
    @test abs(sum(two_step_repeat.coefficients) - sum(mask.coefficients)^2) < 1.0e-9
    @test abs(sum(three_step_repeat.coefficients) - sum(mask.coefficients)^3) < 1.0e-7
    @test three_step_repeat.offset == -(mask.factor^3 - 1) ÷ (mask.factor - 1) * mask.support_radius
end

@testset "Mapped ordinary PGDG intermediate layer" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    basis = build_basis(MappedUniformBasisSpec(:G10;
        count = 5,
        mapping = fit_asinh_mapping_for_strength(s = 0.5, npoints = 5, xmax = 6.0),
        reference_spacing = 1.0,
    ))

    intermediate = GaussletBases._mapped_ordinary_pgdg_intermediate_1d(
        basis;
        exponents = expansion.exponents[1:3],
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents[1:3],
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    experimental_bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents[1:3],
        backend = :pgdg_experimental,
        refinement_levels = 0,
    )
    localized_bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents[1:3],
        backend = :pgdg_localized_experimental,
        refinement_levels = 0,
    )

    @test intermediate isa GaussletBases._MappedOrdinaryPGDGIntermediate1D
    @test intermediate.refinement_levels == 0
    @test intermediate.refinement_mask.factor == 3
    @test intermediate.refinement_mask.rho ≈ 1.2 atol = 0.0 rtol = 0.0
    @test intermediate.base_layer !== intermediate.auxiliary_layer
    @test size(intermediate.overlap) == (length(basis), length(basis))
    @test size(intermediate.kinetic) == (length(basis), length(basis))
    @test size(intermediate.position) == (length(basis), length(basis))
    @test size(intermediate.x2) == (length(basis), length(basis))
    @test length(intermediate.gaussian_factors) == 3
    @test size(intermediate.gaussian_factor_terms) == (3, length(basis), length(basis))
    @test length(intermediate.pair_factors) == 3
    @test size(intermediate.pair_factor_terms) == (3, length(basis), length(basis))
    @test length(intermediate.weights) == length(basis)
    @test length(intermediate.centers) == length(basis)
    @test intermediate.overlap ≈ transpose(intermediate.overlap) atol = 1.0e-10 rtol = 1.0e-10
    @test norm(intermediate.overlap - I, Inf) < 1.0e-10
    @test intermediate.kinetic ≈ transpose(intermediate.kinetic) atol = 1.0e-10 rtol = 1.0e-10
    @test maximum(abs.(intermediate.pair_factor_terms[1, :, :] .- transpose(intermediate.pair_factor_terms[1, :, :]))) < 1.0e-10
    @test bundle.pgdg_intermediate.refinement_levels == 0
    @test bundle.pgdg_intermediate.gaussian_factor_terms ≈ intermediate.gaussian_factor_terms atol = 0.0 rtol = 0.0
    @test GaussletBases._supports_analytic_gaussian_backend(primitive_set(bundle.pgdg_intermediate.base_layer))
    @test GaussletBases._supports_analytic_gaussian_backend(primitive_set(bundle.pgdg_intermediate.auxiliary_layer))
    for candidate in (experimental_bundle, localized_bundle)
        @test GaussletBases._supports_analytic_gaussian_backend(primitive_set(candidate.layer))
        @test GaussletBases._supports_analytic_gaussian_backend(primitive_set(candidate.pgdg_intermediate.base_layer))
        @test GaussletBases._supports_analytic_gaussian_backend(primitive_set(candidate.pgdg_intermediate.auxiliary_layer))
    end
    @test_throws ArgumentError GaussletBases._mapped_ordinary_pgdg_intermediate_1d(
        basis;
        exponents = expansion.exponents[1:3],
        backend = :numerical_reference,
        refinement_levels = 1,
    )
    @test_throws ArgumentError GaussletBases._mapped_ordinary_pgdg_intermediate_1d(
        basis;
        exponents = expansion.exponents[1:3],
        backend = :pgdg_experimental,
        working_layer = basis,
        refinement_levels = 0,
    )
    @test_throws ArgumentError GaussletBases._mapped_ordinary_pgdg_intermediate_1d(
        basis;
        exponents = expansion.exponents[1:3],
        backend = :pgdg_localized_experimental,
        working_layer = basis,
        refinement_levels = 0,
    )
end

@testset "Basis representation" begin
    ub = build_basis(UniformBasisSpec(:G10; xmin = -1.0, xmax = 1.0, spacing = 1.0))
    rep = basis_representation(ub)
    primitive_layer = primitive_set(ub)
    metadata = basis_metadata(ub)

    overlap_public = overlap_matrix(primitive_layer)
    position_public = position_matrix(primitive_layer)
    kinetic_public = kinetic_matrix(primitive_layer)

    @test rep.metadata.basis_kind == metadata.basis_kind
    @test rep.metadata.family_name == metadata.family_name
    @test rep.coefficient_matrix ≈ stencil_matrix(ub) atol = 1.0e-12 rtol = 1.0e-12
    @test length(rep.primitive_set) == length(primitive_layer)
    @test rep.primitive_matrices.overlap ≈ overlap_public atol = 1.0e-12 rtol = 1.0e-12
    @test rep.primitive_matrices.position ≈ position_public atol = 1.0e-12 rtol = 1.0e-12
    @test rep.primitive_matrices.kinetic ≈ kinetic_public atol = 1.0e-12 rtol = 1.0e-12
    @test rep.basis_matrices.overlap ≈ contract_primitive_matrix(ub, overlap_public) atol = 1.0e-12 rtol = 1.0e-12
    @test rep.basis_matrices.position ≈ contract_primitive_matrix(ub, position_public) atol = 1.0e-12 rtol = 1.0e-12
    @test rep.basis_matrices.kinetic ≈ contract_primitive_matrix(ub, kinetic_public) atol = 1.0e-12 rtol = 1.0e-12
    @test rep.basis_matrices.position ≈ _midpoint_reference_position_matrix(ub) atol = 1.0e-12 rtol = 1.0e-12
end

@testset "Basis partitions" begin
    ub = build_basis(UniformBasisSpec(:G10; xmin = -2.0, xmax = 2.0, spacing = 1.0))
    rep = basis_representation(ub)
    partition = basis_partition(rep, [-2.5, -0.5, 0.5, 2.5])
    overlap = rep.basis_matrices.overlap
    kinetic = rep.basis_matrices.kinetic

    assigned = sort!(vcat([box_indices(partition, i) for i in 1:length(boxes(partition))]...))
    @test assigned == collect(1:length(ub))
    @test length(unique(assigned)) == length(ub)

    assembled = zeros(Float64, length(ub), length(ub))
    for i in 1:length(boxes(partition))
        for j in 1:length(boxes(partition))
            assembled[box_indices(partition, i), box_indices(partition, j)] .=
                box_coupling(overlap, partition, i, j)
        end
    end

    @test assembled ≈ overlap atol = 1.0e-12 rtol = 1.0e-12
    @test box_block(rep, partition, :overlap, 1) ≈ overlap[1:2, 1:2] atol = 1.0e-12 rtol = 1.0e-12
    @test box_coupling(rep, partition, :kinetic, 1, 2) ≈ kinetic[1:2, 3:3] atol = 1.0e-12 rtol = 1.0e-12
end

@testset "Hierarchical basis partitions" begin
    ub = build_basis(UniformBasisSpec(:G10; xmin = -2.0, xmax = 2.0, spacing = 1.0))
    rep = basis_representation(ub)
    base_partition = basis_partition(rep, [-2.5, -0.5, 0.5, 2.5])
    hierarchy = hierarchical_partition(base_partition)
    refined = refine_partition(hierarchy, 1)
    overlap = rep.basis_matrices.overlap

    leaves = leaf_boxes(refined)
    leaf_indices = sort!(vcat([box.basis_indices for box in leaves]...))
    @test leaf_indices == collect(1:length(ub))
    @test length(unique(leaf_indices)) == length(ub)

    assembled = zeros(Float64, length(ub), length(ub))
    for box_i in leaves
        for box_j in leaves
            assembled[box_i.basis_indices, box_j.basis_indices] .=
                box_coupling(overlap, refined, box_i.index, box_j.index)
        end
    end

    @test assembled ≈ overlap atol = 1.0e-12 rtol = 1.0e-12
    @test box_parent(refined, 4) == 1
    @test box_parent(refined, 5) == 1
    @test box_children(refined, 1) == [4, 5]
    @test box_level(refined, 1) == 0
    @test box_level(refined, 4) == 1
    @test isempty(box_children(refined, 2))
    @test box_indices(refined, 2) == box_indices(base_partition, 2)
    @test box_indices(refined, 3) == box_indices(base_partition, 3)
    @test box_block(rep, refined, :overlap, 4) ≈ overlap[1:1, 1:1] atol = 1.0e-12 rtol = 1.0e-12
    @test box_coupling(rep, refined, :kinetic, 4, 2) ≈ rep.basis_matrices.kinetic[1:1, 3:3] atol = 1.0e-12 rtol = 1.0e-12
end

@testset "Leaf-local PGDG generation" begin
    ub = build_basis(UniformBasisSpec(:G10; xmin = -2.0, xmax = 2.0, spacing = 1.0))
    rep = basis_representation(ub)
    coarse_hierarchy = hierarchical_partition(rep, [-2.5, -0.5, 0.5, 2.5])
    refined_hierarchy = refine_partition(coarse_hierarchy, 1)

    coarse_pgdg = build_leaf_pgdg(coarse_hierarchy; primitives_per_leaf = 2, width_scale = 0.75)
    refined_pgdg = build_leaf_pgdg(refined_hierarchy; primitives_per_leaf = 2, width_scale = 0.75)
    refined_rep = basis_representation(refined_pgdg)

    @test coarse_pgdg isa LeafLocalPGDG1D
    @test refined_pgdg isa LeafLocalPGDG1D
    @test size(stencil_matrix(refined_pgdg), 1) == length(primitive_set(refined_pgdg))
    @test size(stencil_matrix(refined_pgdg), 2) == length(primitive_set(refined_pgdg))
    @test primitive_set(refined_pgdg).name_value == :leaf_pgdg_1d

    @test length(leaf_primitive_indices(coarse_pgdg, 1)) == 2
    @test length(leaf_primitive_indices(refined_pgdg, 4)) == 2
    @test length(leaf_primitive_indices(refined_pgdg, 5)) == 2
    @test length(leaf_primitive_indices(refined_pgdg, 4)) + length(leaf_primitive_indices(refined_pgdg, 5)) >
          length(leaf_primitive_indices(coarse_pgdg, 1))

    centers_coarse = coarse_pgdg.metadata.center_data
    centers_refined = refined_pgdg.metadata.center_data
    @test centers_refined[leaf_primitive_indices(refined_pgdg, 2)] ≈ centers_coarse[leaf_primitive_indices(coarse_pgdg, 2)] atol = 1.0e-12 rtol = 1.0e-12
    @test centers_refined[leaf_primitive_indices(refined_pgdg, 3)] ≈ centers_coarse[leaf_primitive_indices(coarse_pgdg, 3)] atol = 1.0e-12 rtol = 1.0e-12

    @test all(isfinite, refined_rep.basis_matrices.overlap)
    @test all(isfinite, refined_rep.basis_matrices.kinetic)
    @test refined_rep.basis_matrices.overlap ≈ transpose(refined_rep.basis_matrices.overlap) atol = 1.0e-12 rtol = 1.0e-12
    @test refined_rep.basis_matrices.kinetic ≈ transpose(refined_rep.basis_matrices.kinetic) atol = 1.0e-12 rtol = 1.0e-12
    @test length(refined_pgdg.leaf_box_ids) == length(leaf_boxes(refined_hierarchy))
end

@testset "Leaf-local PGDG augmentation" begin
    ub = build_basis(UniformBasisSpec(:G10; xmin = -2.0, xmax = 2.0, spacing = 1.0))
    rep = basis_representation(ub)
    hierarchy = refine_partition(hierarchical_partition(rep, [-2.5, -0.5, 0.5, 2.5]), 1)
    base_pgdg = build_leaf_pgdg(hierarchy; primitives_per_leaf = 2, width_scale = 0.75)
    augmented_pgdg = augment_leaf_pgdg(
        base_pgdg;
        by_leaf = Dict(
            4 => [LeafGaussianSpec1D(relative_position = 0.5, width_scale = 0.2)],
        ),
    )
    augmented_rep = basis_representation(augmented_pgdg)

    @test length(leaf_primitive_indices(augmented_pgdg, 4)) == length(leaf_primitive_indices(base_pgdg, 4)) + 1
    @test length(leaf_primitive_indices(augmented_pgdg, 5)) == length(leaf_primitive_indices(base_pgdg, 5))
    @test length(leaf_primitive_indices(augmented_pgdg, 2)) == length(leaf_primitive_indices(base_pgdg, 2))
    @test length(primitive_set(augmented_pgdg)) == length(primitive_set(base_pgdg)) + 1

    base_centers = base_pgdg.metadata.center_data
    augmented_centers = augmented_pgdg.metadata.center_data
    @test augmented_centers[leaf_primitive_indices(augmented_pgdg, 5)] ≈ base_centers[leaf_primitive_indices(base_pgdg, 5)] atol = 1.0e-12 rtol = 1.0e-12
    @test augmented_centers[leaf_primitive_indices(augmented_pgdg, 2)] ≈ base_centers[leaf_primitive_indices(base_pgdg, 2)] atol = 1.0e-12 rtol = 1.0e-12
    @test augmented_centers[leaf_primitive_indices(augmented_pgdg, 3)] ≈ base_centers[leaf_primitive_indices(base_pgdg, 3)] atol = 1.0e-12 rtol = 1.0e-12

    @test count(==(:augmented), primitive_origins(augmented_pgdg)) == 1
    @test count(==(:generated), primitive_origins(augmented_pgdg)) == length(primitive_set(base_pgdg))
    @test primitive_leaf_boxes(augmented_pgdg)[last(leaf_primitive_indices(augmented_pgdg, 4))] == 4

    @test all(isfinite, augmented_rep.basis_matrices.overlap)
    @test all(isfinite, augmented_rep.basis_matrices.position)
    @test all(isfinite, augmented_rep.basis_matrices.kinetic)
    @test augmented_rep.basis_matrices.overlap ≈ transpose(augmented_rep.basis_matrices.overlap) atol = 1.0e-12 rtol = 1.0e-12
    @test augmented_rep.basis_matrices.position ≈ transpose(augmented_rep.basis_matrices.position) atol = 1.0e-12 rtol = 1.0e-12
    @test augmented_rep.basis_matrices.kinetic ≈ transpose(augmented_rep.basis_matrices.kinetic) atol = 1.0e-12 rtol = 1.0e-12
end

@testset "Global mapped layer and leaf contraction" begin
    mapping = AsinhMapping(c = 0.15, s = 0.15)
    global_layer = build_global_mapped_primitive_layer(
        xmin = -2.0,
        xmax = 2.0,
        mapping = mapping,
        reference_spacing = 0.5,
        width_scale = 1.0,
    )
    global_rep = basis_representation(global_layer)
    coarse_hierarchy = hierarchical_partition(global_layer, [-2.5, -0.5, 0.5, 2.5])
    refined_hierarchy = refine_partition(coarse_hierarchy, 1)
    coarse_contracted = contract_leaf_boxes(global_layer, coarse_hierarchy; retained_per_leaf = 1)
    refined_contracted = contract_leaf_boxes(global_layer, refined_hierarchy; retained_per_leaf = 1)
    refined_rep = basis_representation(refined_contracted)

    @test global_layer isa GlobalMappedPrimitiveLayer1D
    @test size(stencil_matrix(global_layer), 1) == length(primitive_set(global_layer))
    @test size(stencil_matrix(global_layer), 2) == length(primitive_set(global_layer))
    @test all(isfinite, global_rep.basis_matrices.overlap)
    @test all(isfinite, global_rep.basis_matrices.kinetic)
    @test global_rep.basis_matrices.overlap ≈ transpose(global_rep.basis_matrices.overlap) atol = 1.0e-12 rtol = 1.0e-12
    @test global_rep.basis_matrices.kinetic ≈ transpose(global_rep.basis_matrices.kinetic) atol = 1.0e-12 rtol = 1.0e-12

    @test refined_contracted isa LeafBoxContractionLayer1D
    @test length(leaf_contractions(refined_contracted)) == length(leaf_boxes(refined_hierarchy))
    @test size(stencil_matrix(refined_contracted), 1) == length(primitive_set(global_layer))
    @test size(stencil_matrix(refined_contracted), 2) == length(leaf_boxes(refined_hierarchy))
    @test size(stencil_matrix(refined_contracted), 2) == size(stencil_matrix(coarse_contracted), 2) + 1
    @test leaf_contractions(refined_contracted)[1].leaf_box_index == 4
    @test leaf_contractions(refined_contracted)[1].primitive_indices == box_indices(refined_hierarchy, 4)

    untouched_coarse = Dict(contraction.leaf_box_index => contraction.retained_centers for contraction in leaf_contractions(coarse_contracted))
    untouched_refined = Dict(contraction.leaf_box_index => contraction.retained_centers for contraction in leaf_contractions(refined_contracted))
    @test untouched_refined[2] ≈ untouched_coarse[2] atol = 1.0e-12 rtol = 1.0e-12
    @test untouched_refined[3] ≈ untouched_coarse[3] atol = 1.0e-12 rtol = 1.0e-12

    @test all(isfinite, refined_rep.basis_matrices.overlap)
    @test all(isfinite, refined_rep.basis_matrices.kinetic)
    @test refined_rep.basis_matrices.overlap ≈ transpose(refined_rep.basis_matrices.overlap) atol = 1.0e-12 rtol = 1.0e-12
    @test refined_rep.basis_matrices.kinetic ≈ transpose(refined_rep.basis_matrices.kinetic) atol = 1.0e-12 rtol = 1.0e-12
end

