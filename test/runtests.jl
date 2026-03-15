using Test
using LinearAlgebra

using GaussletBases

const _PROJECT_ROOT = dirname(@__DIR__)

function _radial_operator_fixture(; accuracy = :medium, refine = nothing, quadrature_rmax = 12.0)
    rb = build_basis(RadialBasisSpec(:G10;
        count = 6,
        mapping = AsinhMapping(c = 0.15, s = 0.15),
        reference_spacing = 1.0,
        tails = 3,
        odd_even_kmax = 2,
        xgaussians = [XGaussian(alpha = 0.2)],
    ))
    grid = radial_quadrature(
        rb;
        accuracy = accuracy,
        refine = refine,
        quadrature_rmax = quadrature_rmax,
    )
    return rb, grid
end

function _run_example_script(name::AbstractString)
    example_path = joinpath(_PROJECT_ROOT, "examples", name)
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(_PROJECT_ROOT) $example_path`
    return success(cmd)
end

function _midpoint_reference_matrices(basis; xmin = -20.0, xmax = 20.0, h = 0.02)
    points = collect((xmin + 0.5 * h):h:(xmax - 0.5 * h))
    weights = fill(h, length(points))
    values = [basis[j](x) for x in points, j in 1:length(basis)]
    derivatives = [derivative(basis[j], x) for x in points, j in 1:length(basis)]
    overlap = transpose(values) * (weights .* values)
    kinetic = 0.5 .* (transpose(derivatives) * (weights .* derivatives))
    return overlap, kinetic
end

function _midpoint_reference_position_matrix(basis; xmin = -20.0, xmax = 20.0, h = 0.02)
    points = collect((xmin + 0.5 * h):h:(xmax - 0.5 * h))
    weights = fill(h, length(points))
    values = [basis[j](x) for x in points, j in 1:length(basis)]
    return transpose(values) * (((weights .* points)) .* values)
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

@testset "Radial quadrature and diagnostics" begin
    rb, grid = _radial_operator_fixture()
    @test radial_quadrature(rb) isa RadialQuadratureGrid
    @test radial_quadrature(rb; accuracy = :medium) isa RadialQuadratureGrid
    @test_throws ArgumentError radial_quadrature(rb; accuracy = :low)

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

@testset "Quadrature accuracy profiles" begin
    Z = 10.0
    s = 0.2
    rb = build_basis(RadialBasisSpec(:G10;
        rmax = 30.0,
        mapping = AsinhMapping(c = s / (2Z), s = s),
        reference_spacing = 1.0,
        tails = 6,
        odd_even_kmax = 6,
        xgaussians = XGaussian[],
    ))

    grid_medium = radial_quadrature(rb; accuracy = :medium)
    grid_high = radial_quadrature(rb; accuracy = :high)
    grid_veryhigh = radial_quadrature(rb; accuracy = :veryhigh)

    function hydrogenic_error(basis, grid, charge)
        overlap = overlap_matrix(basis, grid)
        hamiltonian = kinetic_matrix(basis, grid) +
                      nuclear_matrix(basis, grid; Z = charge) +
                      centrifugal_matrix(basis, grid; l = 0)
        eigenvalues = eigen(Hermitian(hamiltonian), Hermitian(overlap)).values
        return abs(minimum(real(eigenvalues)) + 0.5 * charge * charge)
    end

    err_medium = hydrogenic_error(rb, grid_medium, Z)
    err_high = hydrogenic_error(rb, grid_high, Z)
    err_veryhigh = hydrogenic_error(rb, grid_veryhigh, Z)

    @test length(quadrature_points(grid_medium)) <= length(quadrature_points(grid_high))
    @test length(quadrature_points(grid_high)) <= length(quadrature_points(grid_veryhigh))
    @test err_high <= err_medium + 1.0e-12
    @test err_veryhigh <= err_high + 1.0e-12
    @test basis_diagnostics(rb; accuracy = :veryhigh).overlap_error <=
          basis_diagnostics(rb; accuracy = :medium).overlap_error + 1.0e-9
end

@testset "Recommended radial diagnostics cutoff" begin
    for Z in (2.0, 10.0)
        s = 0.2
        rb = build_basis(RadialBasisSpec(:G10;
            rmax = 30.0,
            mapping = AsinhMapping(c = s / (2Z), s = s),
            reference_spacing = 1.0,
            tails = 6,
            odd_even_kmax = 6,
            xgaussians = XGaussian[],
        ))

        Z == 2.0 && @test_logs (:warn, r"truncating basis tails") radial_quadrature(rb; refine = 24, quadrature_rmax = 30.0)

        diag = basis_diagnostics(rb)
        @test diag.overlap_error < 1.0e-5
        @test diag.D < 1.0e-3
    end
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
    rep = basis_representation(ub)

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
    @test occursin("BasisRepresentation1D(kind=:uniform", sprint(show, rep))
end

@testset "Documentation consistency" begin
    design = read(joinpath(_PROJECT_ROOT, "DESIGN.md"), String)
    readme = read(joinpath(_PROJECT_ROOT, "README.md"), String)
    atomic_setup = read(joinpath(_PROJECT_ROOT, "docs", "recommended_atomic_setup.md"), String)
    radial_workflow = read(joinpath(_PROJECT_ROOT, "docs", "first_radial_workflow.md"), String)
    primitive_layer_note = read(joinpath(_PROJECT_ROOT, "docs", "intermediate_primitive_layer.md"), String)
    terminology = read(joinpath(_PROJECT_ROOT, "docs", "terminology.md"), String)
    roadmap = read(joinpath(_PROJECT_ROOT, "ROADMAP.md"), String)
    status = read(joinpath(_PROJECT_ROOT, "STATUS.md"), String)

    @test !occursin("primitive_kinetic_matrix", design)
    @test !occursin("CombinedMapping", design)
    @test !occursin("ScaledMapping", design)
    @test !occursin("NoDiagonalApproximation", design)
    @test occursin("atomic_operators", readme)
    @test occursin("examples/04_hydrogen_ground_state.jl", readme)
    @test occursin("recommended_atomic_setup.md", readme)
    @test occursin("first_radial_workflow.md", readme)
    @test occursin("terminology.md", readme)
    @test occursin("ROADMAP.md", readme)
    @test occursin("radial_quadrature(rb)", readme)
    @test occursin("accuracy = :high", readme)
    @test occursin("quadrature_rmax", readme)
    @test occursin("Pkg.add(url = \"https://github.com/srwhite59/GaussletBases.jl\")", readme)
    @test occursin("s = 0.2", atomic_setup)
    @test occursin("odd_even_kmax = 6", atomic_setup)
    @test occursin("accuracy = :veryhigh", atomic_setup)
    @test occursin("examples/04_hydrogen_ground_state.jl", radial_workflow)
    @test occursin("accuracy = :medium", radial_workflow)
    @test occursin("PrimitiveSet1D", primitive_layer_note)
    @test occursin("BasisRepresentation1D", primitive_layer_note)
    @test occursin("primitive_set(basis)", primitive_layer_note)
    @test occursin("position_matrix(set::PrimitiveSet1D)", primitive_layer_note)
    @test occursin("basis_representation(basis; operators = (:overlap, :position, :kinetic))", primitive_layer_note)
    @test occursin("analytic path", primitive_layer_note)
    @test occursin("numerical path", primitive_layer_note)
    @test occursin("the basis is not the quadrature grid", terminology)
    @test occursin("more conventional gausslet functionality", roadmap)
    @test occursin("accuracy = :medium", roadmap)
    @test !occursin("Gausslets.jl", readme)
    @test occursin("Exact non-diagonal electron-electron API", status)
end

@testset "Example scripts" begin
    @test _run_example_script("01_first_gausslet.jl")
    @test _run_example_script("02_radial_basis.jl")
    @test _run_example_script("03_radial_operators.jl")
    @test _run_example_script("04_hydrogen_ground_state.jl")
    @test _run_example_script("05_primitive_sets.jl")
    @test _run_example_script("06_basis_contraction.jl")
    @test _run_example_script("07_position_contraction.jl")
    @test _run_example_script("08_basis_representation.jl")
end

@testset "README example slice" begin
    Z = 1.0
    s = 0.2
    map = AsinhMapping(c = s / (2Z), s = s)
    rb = build_basis(RadialBasisSpec(:G10;
        rmax = 30.0,
        mapping = map,
        reference_spacing = 1.0,
        tails = 6,
        odd_even_kmax = 6,
        xgaussians = XGaussian[],
    ))
    f = rb[4]
    grid = radial_quadrature(rb)
    primitive_data = primitives(rb)
    coefficient_matrix = stencil_matrix(rb)
    diag = basis_diagnostics(rb)
    S = overlap_matrix(rb, grid)
    H = kinetic_matrix(rb, grid) +
        nuclear_matrix(rb, grid; Z = Z) +
        centrifugal_matrix(rb, grid; l = 0)
    eig = eigen(Hermitian(H), Hermitian(S))
    E0 = minimum(real(eig.values))
    ops = atomic_operators(rb, grid; Z = 2.0, lmax = 2)

    @test sum(coefficient_matrix[mu, 4] * primitive_data[mu](0.3) for mu in eachindex(primitive_data)) ≈
          f(0.3) atol = 1.0e-12 rtol = 1.0e-12
    @test isfinite(moment_center(f, grid))
    @test diag.overlap_error < 1.0e-5
    @test diag.D < 1.0e-3
    @test quadrature_points(grid)[end] >= center(f)
    @test quadrature_weights(grid)[1] > 0.0
    @test isfinite(E0)
    @test E0 < 0.0
    @test norm(S - I, Inf) < 1.0e-5
    @test size(ops.overlap) == (length(rb), length(rb))
    @test size(centrifugal(ops, 2)) == (length(rb), length(rb))
    @test size(multipole(ops, 1)) == (length(rb), length(rb))
end
