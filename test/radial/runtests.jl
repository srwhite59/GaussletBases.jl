@testset "Recommended xgaussian presets" begin
    @test isempty(recommended_xgaussians(0))
    @test [g.alpha for g in recommended_xgaussians(1)] == [0.0936]
    @test [g.alpha for g in recommended_xgaussians(2)] == [0.0936, 0.0236]
    @test_throws ArgumentError recommended_xgaussians(-1)
    @test_throws ArgumentError recommended_xgaussians(3)

    spec_default = RadialBasisSpec(:G10; rmax = 8.0, mapping = AsinhMapping(c = 0.1, s = 0.2))
    spec_none = RadialBasisSpec(:G10; rmax = 8.0, mapping = AsinhMapping(c = 0.1, s = 0.2), xgaussian_count = 0)
    spec_explicit = RadialBasisSpec(:G10; rmax = 8.0, mapping = AsinhMapping(c = 0.1, s = 0.2), xgaussians = XGaussian[])

    @test [g.alpha for g in spec_default.xgaussians] == [0.0936, 0.0236]
    @test isempty(spec_none.xgaussians)
    @test isempty(spec_explicit.xgaussians)
    @test spec_default.rmax_count_policy == :ceil_reference
    @test_throws ArgumentError RadialBasisSpec(
        :G10;
        rmax = 8.0,
        mapping = AsinhMapping(c = 0.1, s = 0.2),
        xgaussian_count = 0,
        xgaussians = recommended_xgaussians(),
    )
    @test_throws ArgumentError RadialBasisSpec(
        :G10;
        rmax = 8.0,
        mapping = AsinhMapping(c = 0.1, s = 0.2),
        rmax_count_policy = :bogus,
    )
end

@testset "Radial rmax count policy can imitate legacy strict trim" begin
    cases = (
        (2.0, 10.0, 31, 30),
        (10.0, 30.0, 47, 46),
    )
    for (Z, rmax, nceil, nlegacy) in cases
        mapping = AsinhMapping(c = 0.2 / (2Z), s = 0.2)
        spec_current = RadialBasisSpec(
            :G10;
            rmax = rmax,
            mapping = mapping,
            tails = 6,
            odd_even_kmax = 6,
            xgaussian_count = 2,
        )
        spec_legacy = RadialBasisSpec(
            :G10;
            rmax = rmax,
            mapping = mapping,
            tails = 6,
            odd_even_kmax = 6,
            xgaussian_count = 2,
            rmax_count_policy = :legacy_strict_trim,
        )
        basis_current = build_basis(spec_current)
        basis_legacy = build_basis(spec_legacy)
        @test length(basis_current) == nceil
        @test length(basis_legacy) == nlegacy
        @test maximum(reference_centers(basis_current)) ≈ maximum(reference_centers(basis_legacy)) + 1.0 atol = 1.0e-10 rtol = 1.0e-10
    end
end

@testset "Named paper-parity radial prototype" begin
    @test radial_boundary_prototype_names() == [:paper_parity_g10_k6_x2]
    @test isfile(joinpath(_PROJECT_ROOT, "data", "radial", "paper_parity_g10_k6_x2.jld2"))
    @test_throws ArgumentError radial_boundary_prototype(:bogus_named_prototype)

    prototype = radial_boundary_prototype()
    @test prototype.name == :paper_parity_g10_k6_x2
    @test prototype.family_value.name == :G10
    @test prototype.reference_spacing == 1.0
    @test prototype.odd_seed_half_width == 24
    @test prototype.even_tail_kmax == 6
    @test [gaussian.alpha for gaussian in prototype.xgaussians] ==
          [0.09358986806, 0.02357750369]
    @test prototype.stage_dimensions ==
          (
              seed_count = 49,
              raw_odd_count = 24,
              raw_even_count = 7,
              cleaned_even_count = 6,
              xgaussian_count = 2,
              canonical_base_count = 51,
              final_dimension = 32,
              runtime_primitive_count = 411,
          )
    @test size(prototype.canonical_coefficients_big) == (51, 32)
    @test size(prototype.canonical_coefficients) == (51, 32)
    @test size(prototype.runtime_coefficients) == (411, 32)
    @test length(prototype.runtime_primitives) == 411
    @test prototype.diagnostics.expected_final_dimension == 32
    @test prototype.diagnostics.retained_dimension == 32
    @test prototype.diagnostics.mode_drop_count == 0
    @test prototype.diagnostics.smallest_overlap_eigenvalue > 1.0e-12
    @test prototype.diagnostics.overlap_identity_error < 1.0e-25
    @test prototype.diagnostics.evaluation_overlap_identity_error < 1.0e-8
    @test prototype.diagnostics.D < 5.0e-5
    @test prototype.diagnostics.centers_monotone
    @test prototype.checksums.reference_centers ==
          "1a62d84b4572b588dbc3bb0d872dbc62804aa481"
    @test prototype.checksums.canonical_coefficients ==
          "5b29d6bdcaf6b6001b80d1b328b47d250ec57490"
    @test prototype.checksums.runtime_coefficients ==
          "8e5a2d1d14cc704db3f65aabb19b4f38825fe436"
    @test prototype.checksums.sampled_basis_u0_10_du0_01 ==
          "d5394fabb5f9376897175b806a1262e2cb7a536d"

    analytic_basis = build_basis(prototype)
    diag_grid = GaussletBases._paper_parity_prototype_grid(IdentityMapping())
    analytic_diag = basis_diagnostics(analytic_basis, diag_grid)
    @test reference_centers(analytic_basis) ≈ prototype.reference_centers atol = 1.0e-12 rtol = 1.0e-12
    @test analytic_diag.overlap_error ≈
          prototype.diagnostics.evaluation_overlap_identity_error atol = 1.0e-12 rtol = 1.0e-8
    @test analytic_diag.D ≈ prototype.diagnostics.D atol = 1.0e-12 rtol = 1.0e-8

    numerical_basis = GaussletBases._paper_parity_numerical_reference_basis()
    @test maximum(abs.(reference_centers(analytic_basis) .- reference_centers(numerical_basis))) <
          5.0e-10
    points = Float64[0.01 * i for i in 0:1000]
    analytic_sample = _basis_sample_matrix(analytic_basis, points)
    numerical_sample = _basis_sample_matrix(numerical_basis, points)
    @test maximum(abs.(analytic_sample .- numerical_sample)) < 5.0e-9
    @test norm(analytic_sample .- numerical_sample) / sqrt(length(analytic_sample)) < 1.0e-9

    paper_map = AsinhMapping(c = 0.2 / (2 * 10.0), s = 0.2)
    analytic_mapped = build_basis(prototype; mapping = paper_map)
    numerical_mapped = GaussletBases._paper_parity_numerical_reference_basis(mapping = paper_map)
    grid_points, grid_weights = GaussletBases._make_erf_grid(
        ;
        h = 0.001,
        rmax = 80.0,
        sigma = 3.0,
        s0 = 6.5,
    )
    paper_grid = RadialQuadratureGrid(grid_points, grid_weights; mapping = paper_map)
    analytic_ops = atomic_operators(analytic_mapped, paper_grid; Z = 10.0, lmax = 2)
    numerical_ops = atomic_operators(numerical_mapped, paper_grid; Z = 10.0, lmax = 2)
    @test opnorm(analytic_ops.overlap - numerical_ops.overlap, Inf) /
          opnorm(numerical_ops.overlap, Inf) < 1.0e-9
    @test opnorm(analytic_ops.kinetic - numerical_ops.kinetic, Inf) /
          opnorm(numerical_ops.kinetic, Inf) < 1.0e-9
    @test opnorm(analytic_ops.nuclear - numerical_ops.nuclear, Inf) /
          opnorm(numerical_ops.nuclear, Inf) < 1.0e-9
    @test opnorm(centrifugal(analytic_ops, 1) - centrifugal(numerical_ops, 1), Inf) /
          opnorm(centrifugal(numerical_ops, 1), Inf) < 2.0e-7
    @test opnorm(multipole(analytic_ops, 0) - multipole(numerical_ops, 0), Inf) /
          opnorm(multipole(numerical_ops, 0), Inf) < 1.0e-9
end

@testset "Paper-parity radial prototype plus tail extension" begin
    prototype = radial_boundary_prototype()

    ne_map = AsinhMapping(c = 0.2 / (2 * 10.0), s = 0.2)
    ne_extended = build_paper_parity_radial_basis(prototype; rmax = 30.0, mapping = ne_map)
    ne_overload = build_basis(prototype; rmax = 30.0, mapping = ne_map)
    ne_reference = build_basis(
        RadialBasisSpec(
            :G10;
            rmax = 30.0,
            mapping = ne_map,
            reference_spacing = 1.0,
            tails = 6,
            odd_even_kmax = 6,
            xgaussians = prototype.xgaussians,
            rmax_count_policy = :legacy_strict_trim,
        ),
    )
    @test length(ne_extended) == 46
    @test last(centers(ne_extended)) ≈ 28.34082360929913 atol = 1.0e-9 rtol = 1.0e-12
    @test reference_centers(ne_overload) ≈ reference_centers(ne_extended) atol = 1.0e-12 rtol = 1.0e-12
    @test maximum(abs.(reference_centers(ne_extended) .- reference_centers(ne_reference))) < 5.0e-8
    ne_points = collect(0.0:0.05:30.0)
    ne_sample_extended = _basis_sample_matrix(ne_extended, ne_points)
    ne_sample_reference = _basis_sample_matrix(ne_reference, ne_points)
    @test maximum(abs.(ne_sample_extended .- ne_sample_reference)) < 1.0e-6
    @test norm(ne_sample_extended .- ne_sample_reference) / sqrt(length(ne_sample_extended)) < 1.0e-8

    he_map = AsinhMapping(c = 0.2 / (2 * 2.0), s = 0.2)
    he_extended = build_paper_parity_radial_basis(prototype; rmax = 10.0, mapping = he_map)
    he_reference = build_basis(
        RadialBasisSpec(
            :G10;
            rmax = 10.0,
            mapping = he_map,
            reference_spacing = 1.0,
            tails = 6,
            odd_even_kmax = 6,
            xgaussians = prototype.xgaussians,
            rmax_count_policy = :legacy_strict_trim,
        ),
    )
    @test length(he_extended) == 30
    @test length(he_extended) < length(build_basis(prototype; mapping = he_map))
    @test maximum(abs.(reference_centers(he_extended) .- reference_centers(he_reference))) < 5.0e-8
    he_points = collect(0.0:0.05:10.0)
    he_sample_extended = _basis_sample_matrix(he_extended, he_points)
    he_sample_reference = _basis_sample_matrix(he_reference, he_points)
    @test maximum(abs.(he_sample_extended .- he_sample_reference)) < 1.0e-6
    @test norm(he_sample_extended .- he_sample_reference) / sqrt(length(he_sample_extended)) < 1.0e-8
end

@testset "Runtime family tables use trimmed machine-significant tails" begin
    high_prec_path = joinpath(_PROJECT_ROOT, "src", "internal", "families_high_prec.jl")
    @test isfile(high_prec_path)

    expected_runtime_radii = Dict(:G4 => 54, :G6 => 48, :G8 => 67, :G10 => 75)
    expected_high_prec_radii = Dict(:G4 => 91, :G6 => 104, :G8 => 118, :G10 => 132)
    for family_name in (:G4, :G6, :G8, :G10)
        runtime_radius = length(GaussletFamily(family_name).positive_coefficients) - 1
        high_prec_radius = length(_HIGH_PREC_FAMILY_REFERENCE[family_name]) - 1
        @test runtime_radius == expected_runtime_radii[family_name]
        @test high_prec_radius == expected_high_prec_radii[family_name]
        @test high_prec_radius > runtime_radius
    end

    gausslet = Gausslet(:G10; center = 0.0, spacing = 1.0)
    sample_points = collect(-8.0:0.25:8.0)
    runtime_values = Float64[direct_value(gausslet, point) for point in sample_points]
    high_prec_values = Float64[
        _high_prec_gausslet_value(:G10, point; center = 0.0, spacing = 1.0)
        for point in sample_points
    ]
    @test maximum(abs.(runtime_values .- high_prec_values)) < 1.0e-12

    overlap_points = -10.0:0.02:10.0
    h = step(overlap_points)
    runtime_basis = [
        Float64[direct_value(Gausslet(:G10; center = center_value, spacing = 1.0), point) for point in overlap_points]
        for center_value in (-1.0, 0.0, 1.0)
    ]
    high_prec_basis = [
        Float64[_high_prec_gausslet_value(:G10, point; center = center_value, spacing = 1.0) for point in overlap_points]
        for center_value in (-1.0, 0.0, 1.0)
    ]
    runtime_overlap = h .* reduce(hcat, runtime_basis)' * reduce(hcat, runtime_basis)
    high_prec_overlap = h .* reduce(hcat, high_prec_basis)' * reduce(hcat, high_prec_basis)
    @test runtime_overlap ≈ high_prec_overlap atol = 1.0e-11 rtol = 1.0e-11

    Z = 10.0
    s = 0.2
    radial_basis = build_basis(RadialBasisSpec(:G10;
        rmax = 30.0,
        mapping = AsinhMapping(c = s / (2Z), s = s),
        reference_spacing = 1.0,
        tails = 6,
        odd_even_kmax = 6,
        xgaussians = XGaussian[],
    ))
    runtime_tail_bound = GaussletBases._radial_quadrature_tail_bound(radial_basis)
    high_prec_tail_bound = _high_prec_radial_tail_bound(radial_basis)

    @test runtime_tail_bound < high_prec_tail_bound
    @test runtime_tail_bound / high_prec_tail_bound < 0.8
end

@testset "Radial quadrature extent policy" begin
    Z = 2.0
    s = 0.2
    map = AsinhMapping(c = s / (2Z), s = s)

    basis_tails4 = build_basis(RadialBasisSpec(:G10;
        count = 11,
        mapping = map,
        reference_spacing = 1.0,
        tails = 4,
        odd_even_kmax = 2,
        xgaussian_count = 0,
    ))
    basis_tails6 = build_basis(RadialBasisSpec(:G10;
        count = 11,
        mapping = map,
        reference_spacing = 1.0,
        tails = 6,
        odd_even_kmax = 2,
        xgaussian_count = 0,
    ))

    primitive_spacing = basis_tails4.spec.reference_spacing / 3.0
    family_radius = length(family(basis_tails4).positive_coefficients) - 1
    expected_quadrature_umax =
        maximum(reference_centers(basis_tails4)) +
        primitive_spacing * (family_radius + 12.0)

    @test GaussletBases._radial_quadrature_umax(basis_tails4) ≈ expected_quadrature_umax atol = 1.0e-12 rtol = 1.0e-12
    @test GaussletBases._radial_quadrature_umax(basis_tails6) ≈ expected_quadrature_umax atol = 1.0e-12 rtol = 1.0e-12
    @test GaussletBases._radial_build_umax(basis_tails6) > GaussletBases._radial_build_umax(basis_tails4)
    @test abs(GaussletBases._radial_quadrature_umax(basis_tails6) - GaussletBases._radial_quadrature_umax(basis_tails4)) < 1.0e-9
    @test abs(GaussletBases._radial_quadrature_tail_bound(basis_tails6) - GaussletBases._radial_quadrature_tail_bound(basis_tails4)) < 1.0e-7
end

@testset "Interval-sampled radial setup layer" begin
    family = GaussletFamily(:G10)
    spacing = 1.0
    seed = Gausslet(family; center = 0.0, spacing = spacing)
    js = collect(-2:2)
    xgrid, weights = GaussletBases._make_erf_grid(; h = 0.05, rmax = 6.0)
    phi = (j::Int, x::Float64) -> seed(x - j * spacing)

    _, _, sampled_dense = GaussletBases._seed_scalar_integrals(phi, js, xgrid, weights)
    overlap_dense = sampled_dense' * (weights .* sampled_dense)
    position_dense = sampled_dense' * ((xgrid .* weights) .* sampled_dense)

    _, sampled_intervals = GaussletBases._sample_shifted_gausslets(family, js, xgrid, spacing)
    sampled_interval = GaussletBases._interval_sample_matrix(sampled_intervals, length(xgrid))
    overlap_interval = GaussletBases._interval_gram_matrix(sampled_intervals, weights)
    position_interval = GaussletBases._interval_gram_matrix(sampled_intervals, xgrid .* weights)

    @test sampled_interval ≈ sampled_dense atol = 1.0e-12 rtol = 1.0e-12
    @test overlap_interval ≈ overlap_dense atol = 1.0e-12 rtol = 1.0e-12
    @test position_interval ≈ position_dense atol = 1.0e-12 rtol = 1.0e-12

    xgaussians = recommended_xgaussians(2)
    sampled_x_dense = GaussletBases._xgaussian_sample_matrix(xgaussians, xgrid)
    sampled_x_interval = GaussletBases._interval_sample_matrix(
        GaussletBases._sample_xgaussian_intervals(xgaussians, xgrid),
        length(xgrid),
    )

    @test sampled_x_interval ≈ sampled_x_dense atol = 1.0e-12 rtol = 1.0e-12
end

@testset "Recommended radial front-door hydrogen" begin
    Z = 1.0
    s = 0.2
    rb = build_basis(RadialBasisSpec(:G10;
        rmax = 30.0,
        mapping = AsinhMapping(c = s / (2Z), s = s),
    ))
    diag = @test_logs min_level = Logging.Warn basis_diagnostics(rb)
    grid = @test_logs min_level = Logging.Warn radial_quadrature(rb)
    hamiltonian = kinetic_matrix(rb, grid) +
                  nuclear_matrix(rb, grid; Z = Z) +
                  centrifugal_matrix(rb, grid; l = 0)
    ground_energy = minimum(real(eigen(Hermitian(hamiltonian)).values))

    @test length(rb) == 35
    @test diag.overlap_error < 1.5e-5
    @test diag.D < 1.0e-5
    @test abs(ground_energy + 0.5) < 1.0e-6
end

@testset "Radial quadrature and diagnostics" begin
    rb, grid = _quick_radial_operator_fixture()
    @test radial_quadrature(rb) isa RadialQuadratureGrid
    @test radial_quadrature(rb; accuracy = :medium) isa RadialQuadratureGrid
    @test_throws ArgumentError radial_quadrature(rb; accuracy = :low)

    points = quadrature_points(grid)
    weights = quadrature_weights(grid)
    diag_rb = basis_diagnostics(rb, grid)
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
end

if _RUN_SLOW_TESTS
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
            @test norm(overlap - I, Inf) ≤ 1.0e-5
            eigenvalues = eigen(Hermitian(hamiltonian)).values
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
end

@testset "Radial operator matrices" begin
    rb, grid = _quick_radial_operator_fixture()
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

@testset "Stabilized radial multipole builder" begin
    rb, grid = _quick_radial_operator_fixture()
    points = quadrature_points(grid)
    weights = quadrature_weights(grid)
    values = GaussletBases._basis_values_matrix(rb, points)

    raw_benign = GaussletBases._integral_diagonal_kernel_matrix_raw(values, points, weights, 3)
    stable_benign = GaussletBases._integral_diagonal_kernel_matrix(values, points, weights, 3)
    @test opnorm(stable_benign - raw_benign, Inf) / opnorm(raw_benign, Inf) < 1.0e-12

    risky_points = exp.(range(log(1.0e-6), log(1.0e6), length = 401))
    risky_weights = fill((risky_points[end] - risky_points[1]) / (length(risky_points) - 1), length(risky_points))
    risky_values = hcat(
        exp.(-risky_points ./ 10.0),
        sqrt.(risky_points) .* exp.(-risky_points ./ 20.0),
    )

    raw_risky =
        GaussletBases._integral_diagonal_kernel_matrix_raw(risky_values, risky_points, risky_weights, 120)
    stable_risky =
        GaussletBases._integral_diagonal_kernel_matrix(risky_values, risky_points, risky_weights, 120)

    @test !all(isfinite, raw_risky)
    @test all(isfinite, stable_risky)
    @test stable_risky ≈ transpose(stable_risky) atol = 1.0e-12 rtol = 1.0e-12
end

@testset "Radial primitive operator contraction" begin
    rb, grid = _quick_radial_operator_fixture()
    P = primitive_set(rb)

    overlap_mu = overlap_matrix(P, grid)
    kinetic_mu = kinetic_matrix(P, grid)
    nuclear_mu = nuclear_matrix(P, grid; Z = 2.0)
    centr2_mu = centrifugal_matrix(P, grid; l = 2)

    overlap_b = contract_primitive_matrix(rb, overlap_mu)
    kinetic_b = contract_primitive_matrix(rb, kinetic_mu)
    nuclear_b = contract_primitive_matrix(rb, nuclear_mu)
    centr2_b = contract_primitive_matrix(rb, centr2_mu)

    @test all(isfinite, overlap_mu)
    @test all(isfinite, kinetic_mu)
    @test all(isfinite, nuclear_mu)
    @test all(isfinite, centr2_mu)

    @test overlap_mu ≈ transpose(overlap_mu) atol = 1.0e-12 rtol = 1.0e-12
    @test kinetic_mu ≈ transpose(kinetic_mu) atol = 1.0e-12 rtol = 1.0e-12
    @test nuclear_mu ≈ transpose(nuclear_mu) atol = 1.0e-12 rtol = 1.0e-12
    @test centr2_mu ≈ transpose(centr2_mu) atol = 1.0e-12 rtol = 1.0e-12

    @test norm(overlap_b - overlap_matrix(rb, grid), Inf) ≤ 1.0e-10
    @test norm(kinetic_b - kinetic_matrix(rb, grid), Inf) ≤ 1.0e-10
    @test norm(nuclear_b - nuclear_matrix(rb, grid; Z = 2.0), Inf) ≤ 1.0e-10
    @test norm(centr2_b - centrifugal_matrix(rb, grid; l = 2), Inf) /
          norm(centrifugal_matrix(rb, grid; l = 2), Inf) ≤ 1.0e-9
end

@testset "Radial atomic operators" begin
    rb, grid, ops, _, _, _ = _quick_radial_atomic_fixture()

    @test ops isa RadialAtomicOperators
    @test ops.overlap ≈ overlap_matrix(rb, grid) atol = 1.0e-12 rtol = 1.0e-12
    @test ops.kinetic ≈ kinetic_matrix(rb, grid) atol = 1.0e-12 rtol = 1.0e-12
    @test ops.nuclear ≈ nuclear_matrix(rb, grid; Z = 2.0) atol = 1.0e-12 rtol = 1.0e-12
    @test centrifugal(ops, 2) ≈ centrifugal_matrix(rb, grid; l = 2) atol = 1.0e-12 rtol = 1.0e-12
    @test multipole(ops, 1) ≈ multipole_matrix(rb, grid; L = 1) atol = 1.0e-12 rtol = 1.0e-12
    @test size(multipole(ops, 4)) == (length(rb), length(rb))
    @test_throws BoundsError multipole(ops, 5)
end

