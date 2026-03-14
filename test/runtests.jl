using Test
using LinearAlgebra

using GaussletBases

const _PROJECT_ROOT = dirname(@__DIR__)

function _radial_operator_fixture(; refine = 8, rmax = 12.0)
    rb = build_basis(RadialBasisSpec(:G10;
        count = 6,
        mapping = AsinhMapping(c = 0.15, s = 0.15),
        reference_spacing = 1.0,
        tails = 3,
        odd_even_kmax = 2,
        xgaussians = [XGaussian(alpha = 0.2)],
    ))
    grid = radial_quadrature(rb; refine = refine, rmax = rmax)
    return rb, grid
end

function _run_example_script(name::AbstractString)
    example_path = joinpath(_PROJECT_ROOT, "examples", name)
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(_PROJECT_ROOT) $example_path`
    return success(cmd)
end

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
    primitive_data = primitives(hb)
    coefficient_matrix = stencil_matrix(hb)
    st = stencil(hb[2])

    @test hb isa HalfLineBasis
    @test length(hb) >= 4
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

@testset "Radial quadrature and diagnostics" begin
    rb, grid = _radial_operator_fixture()
    @test_throws ArgumentError radial_quadrature(rb; refine = 8)

    points = quadrature_points(grid)
    weights = quadrature_weights(grid)
    diag_rb = basis_diagnostics(rb, grid)
    diag_ub = basis_diagnostics(build_basis(UniformBasisSpec(:G10; xmin = -2.0, xmax = 2.0, spacing = 1.0)))
    diag_hb = basis_diagnostics(build_basis(HalfLineBasisSpec(:G10;
        xmax = 4.0,
        reference_spacing = 1.0,
        tails = 3,
        mapping = AsinhMapping(a = 1.0, s = 0.2),
    )))
    mc = moment_center(rb[2], grid)

    @test length(points) == length(weights)
    @test length(points) > length(rb)
    @test all(weight -> weight > 0.0, weights)
    @test issorted(points)
    @test uofx(mapping(rb), points[end]) >= maximum(reference_centers(rb))
    @test isfinite(mc)
    @test mc >= 0.0
    @test :overlap_error in propertynames(diag_rb)
    @test :moment_centers in propertynames(diag_rb)
    @test :center_mismatches in propertynames(diag_rb)
    @test :D in propertynames(diag_rb)
    @test isfinite(diag_rb.overlap_error)
    @test isfinite(diag_rb.D)
    @test isfinite(diag_ub.overlap_error)
    @test isfinite(diag_ub.D)
    @test isfinite(diag_hb.overlap_error)
    @test isfinite(diag_hb.D)
end

@testset "Radial operator matrices" begin
    rb, grid = _radial_operator_fixture(; refine = 24)
    points = quadrature_points(grid)
    weights = quadrature_weights(grid)

    overlap = overlap_matrix(rb, grid)
    kinetic = kinetic_matrix(rb, grid)
    nuclear = nuclear_matrix(rb, grid; Z = 2.0)
    centr0 = centrifugal_matrix(rb, grid; l = 0)
    centr2 = centrifugal_matrix(rb, grid; l = 2)
    multipole0 = multipole_matrix(rb, grid; L = 0)
    multipole1 = multipole_matrix(rb, grid; L = 1)

    values = [rb[j](points[i]) for i in eachindex(points), j in 1:length(rb)]
    wchi = vec(transpose(values) * weights)
    kernel1 = [
        weights[i] * weights[j] * min(points[i], points[j]) / max(points[i], points[j])^2
        for i in eachindex(points), j in eachindex(points)
    ]
    multipole1_explicit =
        Diagonal(1.0 ./ wchi) * transpose(values) * kernel1 * values * Diagonal(1.0 ./ wchi)

    @test all(isfinite, overlap)
    @test overlap ≈ transpose(overlap) atol = 1.0e-12 rtol = 1.0e-12
    @test norm(overlap - I, Inf) ≤ 2.0e-3

    @test all(isfinite, kinetic)
    @test kinetic ≈ transpose(kinetic) atol = 1.0e-12 rtol = 1.0e-12

    @test all(isfinite, nuclear)
    @test nuclear ≈ transpose(nuclear) atol = 1.0e-12 rtol = 1.0e-12

    @test all(isfinite, centr0)
    @test centr0 ≈ zeros(Float64, length(rb), length(rb)) atol = 1.0e-12 rtol = 1.0e-12
    @test all(isfinite, centr2)
    @test centr2 ≈ transpose(centr2) atol = 1.0e-12 rtol = 1.0e-12
    @test norm(centr2, Inf) > 1.0e-8

    @test multipole0 isa Matrix{Float64}
    @test all(isfinite, multipole0)
    @test multipole0 ≈ transpose(multipole0) atol = 1.0e-12 rtol = 1.0e-12
    @test all(isfinite, multipole1)
    @test multipole1 ≈ transpose(multipole1) atol = 1.0e-12 rtol = 1.0e-12
    @test multipole1 ≈ multipole1_explicit atol = 1.0e-10 rtol = 1.0e-10
end

@testset "Radial atomic operators" begin
    rb, grid = _radial_operator_fixture(; refine = 24)
    ops = atomic_operators(rb, grid; Z = 2.0, lmax = 2)

    @test ops isa RadialAtomicOperators
    @test ops.overlap ≈ overlap_matrix(rb, grid) atol = 1.0e-12 rtol = 1.0e-12
    @test ops.kinetic ≈ kinetic_matrix(rb, grid) atol = 1.0e-12 rtol = 1.0e-12
    @test ops.nuclear ≈ nuclear_matrix(rb, grid; Z = 2.0) atol = 1.0e-12 rtol = 1.0e-12
    @test centrifugal(ops, 2) ≈ centrifugal_matrix(rb, grid; l = 2) atol = 1.0e-12 rtol = 1.0e-12
    @test multipole(ops, 1) ≈ multipole_matrix(rb, grid; L = 1) atol = 1.0e-12 rtol = 1.0e-12
    @test size(multipole(ops, 4)) == (length(rb), length(rb))
    @test_throws BoundsError multipole(ops, 5)
end

@testset "REPL displays" begin
    family = GaussletFamily(:G10)
    map = AsinhMapping(c = 0.15, s = 0.15)
    ub_spec = UniformBasisSpec(:G10; xmin = -2.0, xmax = 2.0, spacing = 1.0)
    hb_spec = HalfLineBasisSpec(:G10;
        xmax = 4.0,
        reference_spacing = 1.0,
        tails = 3,
        mapping = map,
    )
    rb_spec = RadialBasisSpec(:G10;
        count = 6,
        mapping = map,
        reference_spacing = 1.0,
        tails = 3,
        odd_even_kmax = 2,
        xgaussians = [XGaussian(alpha = 0.2)],
    )
    ub = build_basis(ub_spec)
    hb = build_basis(hb_spec)
    rb, grid = _radial_operator_fixture(; refine = 24)
    ops = atomic_operators(rb, grid; Z = 2.0, lmax = 2)

    @test sprint(show, family) == "GaussletFamily(:G10)"
    @test occursin("AsinhMapping(", sprint(show, map))
    @test occursin("UniformBasisSpec(", sprint(show, ub_spec))
    @test occursin("HalfLineBasisSpec(", sprint(show, hb_spec))
    @test occursin("RadialBasisSpec(", sprint(show, rb_spec))
    @test occursin("UniformBasis(length=5", sprint(show, ub))
    @test occursin("HalfLineBasis(length=", sprint(show, hb))
    @test occursin("RadialBasis(length=6", sprint(show, rb))
    @test occursin("RadialQuadratureGrid(length=", sprint(show, grid))
    @test occursin("RadialAtomicOperators(size=(6, 6)", sprint(show, ops))
end

@testset "Documentation consistency" begin
    design = read(joinpath(_PROJECT_ROOT, "DESIGN.md"), String)
    readme = read(joinpath(_PROJECT_ROOT, "README.md"), String)
    status = read(joinpath(_PROJECT_ROOT, "STATUS.md"), String)

    @test !occursin("primitive_kinetic_matrix", design)
    @test !occursin("CombinedMapping", design)
    @test !occursin("ScaledMapping", design)
    @test !occursin("NoDiagonalApproximation", design)
    @test occursin("atomic_operators", readme)
    @test occursin("examples/", readme)
    @test occursin("Exact non-diagonal electron-electron API", status)
end

@testset "Example scripts" begin
    @test _run_example_script("01_first_gausslet.jl")
    @test _run_example_script("02_radial_basis.jl")
    @test _run_example_script("03_radial_operators.jl")
end

@testset "README example slice" begin
    rb, grid = _radial_operator_fixture(; refine = 24)
    rf = rb[2]
    primitive_data = primitives(rb)
    coefficient_matrix = stencil_matrix(rb)
    Amunu = Matrix{Float64}(I, length(primitive_data), length(primitive_data))
    A = contract_primitive_matrix(rb, Amunu)
    diag = basis_diagnostics(rb, grid)
    ops = atomic_operators(rb, grid; Z = 2.0, lmax = 2)

    @test sum(coefficient_matrix[mu, 2] * primitive_data[mu](0.2) for mu in eachindex(primitive_data)) ≈
          rf(0.2) atol = 1.0e-12 rtol = 1.0e-12
    @test size(A) == (length(rb), length(rb))
    @test isfinite(moment_center(rf, grid))
    @test isfinite(diag.D)
    @test quadrature_points(grid)[end] >= center(rf)
    @test quadrature_weights(grid)[1] > 0.0
    @test size(ops.overlap) == (length(rb), length(rb))
    @test size(centrifugal(ops, 2)) == (length(rb), length(rb))
    @test size(multipole(ops, 1)) == (length(rb), length(rb))
end
