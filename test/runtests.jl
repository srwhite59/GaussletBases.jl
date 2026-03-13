using Test

using Gausslets

@testset "Uniform basis" begin
    ub = build_basis(UniformBasisSpec(:G10; xmin = -2.0, xmax = 2.0, spacing = 1.0))

    @test ub isa UniformBasis
    @test length(ub) == 5
    @test ub[2] isa Gausslet
    @test centers(ub) == [center(ub[i]) for i in 1:length(ub)]
    @test reference_centers(ub) == centers(ub)
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

@testset "XGaussian center" begin
    g = XGaussian(alpha = 0.23)
    @test center(g) == 0.23
end

@testset "Half-line basis" begin
    spec = HalfLineBasisSpec(:G10;
        xmax = 4.0,
        reference_spacing = 1.0,
        tails = 3,
        mapping = AsinhMapping(a = 1.0, s = 0.2),
    )
    hb = build_basis(spec)

    @test hb isa HalfLineBasis
    @test length(hb) >= 4
    @test hb[1] isa BoundaryGausslet
    @test hb[1](0.2) == value(hb[1], 0.2)
    @test hb[1](-0.5) ≈ 0.0 atol = 1.0e-12
    @test stencil(hb[2])(0.3) ≈ direct_value(hb[2], 0.3) atol = 1.0e-12 rtol = 1.0e-12
    @test center(hb[2]) ≈ xofu(mapping(hb), reference_center(hb[2])) atol = 1.0e-12 rtol = 1.0e-12
    @test issorted(reference_centers(hb))
    @test issorted(centers(hb))
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

    @test rb isa RadialBasis
    @test length(rb) == 6
    @test rb[1] isa RadialGausslet
    @test abs(rb[1](0.0)) ≤ 1.0e-8
    @test issorted(reference_centers(rb))
    @test issorted(centers(rb))
    @test stencil(rb[2])(0.3) ≈ direct_value(rb[2], 0.3) atol = 1.0e-12 rtol = 1.0e-12
end

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

@testset "README example slice" begin
    ub = build_basis(UniformBasisSpec(:G10; xmin = -2.0, xmax = 2.0, spacing = 1.0))
    @test ub[3] isa Gausslet

    rb = build_basis(RadialBasisSpec(:G10;
        count = 6,
        mapping = AsinhMapping(c = 0.15, s = 0.15),
        reference_spacing = 1.0,
        tails = 3,
        odd_even_kmax = 2,
        xgaussians = [XGaussian(alpha = 0.2)],
    ))

    g = Gausslet(:G10; center = 0.0, spacing = 1.0)
    st = stencil(g)
    map = AsinhMapping(c = 0.15, s = 0.15)
    x = 0.2

    @test g(x) == value(g, x)
    @test direct_value(g, x) == st(x)
    @test xofu(map, map(3.0)) ≈ 3.0 atol = 1.0e-12 rtol = 1.0e-12
    @test rb[2](0.2) == value(rb[2], 0.2)
    @test stencil(rb[2])(0.2) ≈ direct_value(rb[2], 0.2) atol = 1.0e-12 rtol = 1.0e-12
end
