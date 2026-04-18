using Test
using Logging
using LinearAlgebra
using JLD2
using SparseArrays

using GaussletBases

# The old 1D COMX-cleaned hybrid route is intentionally unexported from the
# public API, but we keep a few internal regression checks for it here.
const hybrid_mapped_ordinary_basis = GaussletBases.hybrid_mapped_ordinary_basis
const HybridMappedOrdinaryBasis1D = GaussletBases.HybridMappedOrdinaryBasis1D

const _PROJECT_ROOT = dirname(@__DIR__)
const _RUN_SLOW_TESTS = get(ENV, "GAUSSLETBASES_SLOW_TESTS", "0") == "1"
const _FIXTURE_CACHE = Dict{Symbol,Any}()
const _TEST_GROUP_ENV = strip(get(ENV, "GAUSSLETBASES_TEST_GROUPS", "all"))
const _AVAILABLE_TEST_GROUPS = (
    :radial,
    :core,
    :nested,
    :ordinary,
    :diatomic,
    :angular,
    :ida,
    :docs,
    :examples,
    :misc,
)
const _FAST_TEST_GROUPS = Set((
    :radial,
    :core,
    :angular,
    :ida,
    :docs,
))

function _parse_selected_test_groups(spec::AbstractString)
    spec_clean = strip(spec)
    isempty(spec_clean) && return Set(_AVAILABLE_TEST_GROUPS)
    lowercase(spec_clean) == "all" && return Set(_AVAILABLE_TEST_GROUPS)

    selected = Set{Symbol}()
    tokens = filter(!isempty, strip.(split(spec_clean, ',')))
    for token in tokens
        lowered = lowercase(token)
        if lowered == "all"
            union!(selected, _AVAILABLE_TEST_GROUPS)
        elseif lowered == "fast"
            union!(selected, _FAST_TEST_GROUPS)
        else
            group = Symbol(lowered)
            group in _AVAILABLE_TEST_GROUPS ||
                throw(
                    ArgumentError(
                        "unknown GAUSSLETBASES_TEST_GROUPS token $(repr(token)); available groups are $(join(string.(_AVAILABLE_TEST_GROUPS), ", ")) plus aliases all, fast",
                    ),
                )
            push!(selected, group)
        end
    end
    return selected
end

const _SELECTED_TEST_GROUPS = _parse_selected_test_groups(_TEST_GROUP_ENV)
_test_group_enabled(group::Symbol) = group in _SELECTED_TEST_GROUPS

if _TEST_GROUP_ENV != "all"
    @info(
        "GaussletBases test selection active",
        groups = sort!(collect(_SELECTED_TEST_GROUPS)),
        slow_tests = _RUN_SLOW_TESTS,
    )
end

module _HighPrecFamilyReference
include(joinpath(dirname(@__DIR__), "src", "internal", "families_high_prec.jl"))
end

const _HIGH_PREC_FAMILY_REFERENCE = _HighPrecFamilyReference._HIGH_PREC_POSITIVE_COEFFICIENTS

_cached_fixture(key::Symbol, builder::Function) = get!(_FIXTURE_CACHE, key) do
    builder()
end

function _high_prec_gausslet_value(
    family_name::Symbol,
    x::Real;
    center::Real = 0.0,
    spacing::Real = 1.0,
)
    positive_coefficients = _HIGH_PREC_FAMILY_REFERENCE[family_name]
    coefficients_full = GaussletBases._full_family_coefficients(positive_coefficients)
    radius = length(positive_coefficients) - 1
    primitive_spacing = Float64(spacing) / 3.0
    scale = 1.0 / sqrt(2.0 * pi * Float64(spacing))

    total = 0.0
    for (offset, coefficient) in zip((-radius):radius, coefficients_full)
        primitive_center = Float64(center) + primitive_spacing * offset
        total += coefficient * exp(-0.5 * ((Float64(x) - primitive_center) / primitive_spacing)^2)
    end
    return scale * total
end

function _high_prec_radial_tail_bound(basis::RadialBasis)
    positive_coefficients = _HIGH_PREC_FAMILY_REFERENCE[basis.spec.family_value.name]
    raw_radius = length(positive_coefficients) - 1
    primitive_spacing = basis.spec.reference_spacing / 3.0
    reference_upper =
        maximum(reference_centers(basis)) +
        primitive_spacing * (raw_radius + 12.0)
    mapped_upper = xofu(mapping(basis), reference_upper)

    xgaussian_upper = -Inf
    for primitive in primitives(basis)
        if primitive isa XGaussian || (primitive isa Distorted && primitive.primitive isa XGaussian)
            _, xhi = GaussletBases._primitive_physical_bounds(primitive)
            xgaussian_upper = max(xgaussian_upper, xhi)
        end
    end

    return max(mapped_upper, xgaussian_upper)
end

function _orthonormalize_1d(overlap::AbstractMatrix, operators::AbstractVector{<:AbstractMatrix})
    decomposition = eigen(Symmetric(Matrix{Float64}(overlap)))
    invhalf =
        decomposition.vectors *
        Diagonal(1.0 ./ sqrt.(decomposition.values)) *
        transpose(decomposition.vectors)
    return Tuple(invhalf * operator * invhalf for operator in operators)
end

function _cleanup_comx_transform(
    overlap::AbstractMatrix,
    position::AbstractMatrix,
    sign_vector::AbstractVector,
)
    vectors, invhalf = GaussletBases._s_invsqrt_reduced(overlap)
    localizer, center_values = GaussletBases._comx_reduced(vectors, invhalf, position)
    transform = Matrix{Float64}(vectors * (invhalf * localizer))
    ordering = sortperm(center_values)
    transform = transform[:, ordering]
    for column in 1:size(transform, 2)
        sign_probe = dot(sign_vector, view(transform, :, column))
        if abs(sign_probe) <= 1.0e-12
            pivot = argmax(abs.(view(transform, :, column)))
            sign_probe = transform[pivot, column]
        end
        if sign_probe < 0.0
            transform[:, column] .*= -1.0
        end
    end
    return transform, Float64[Float64(center_values[index]) for index in ordering]
end

function _basis_sample_matrix(basis_like, points::AbstractVector{<:Real})
    primitive_values = GaussletBases._primitive_sample_matrix(primitive_set(basis_like), Float64[Float64(point) for point in points])
    return primitive_values * Matrix{Float64}(stencil_matrix(basis_like))
end

function _projection_error(basis_like, target::Function; xmin = -10.0, xmax = 10.0, h = 0.05)
    points = collect((xmin + 0.5 * h):h:(xmax - 0.5 * h))
    weights = fill(h, length(points))
    values = _basis_sample_matrix(basis_like, points)
    target_values = Float64[target(point) for point in points]
    gram = transpose(values) * (weights .* values)
    rhs = transpose(values) * (weights .* target_values)
    coefficients = gram \ rhs
    residual = target_values - values * coefficients
    residual_norm = sqrt(sum(weights .* (residual .^ 2)))
    target_norm = sqrt(sum(weights .* (target_values .^ 2)))
    return residual_norm / target_norm
end

function _subspace_overlap_metric(basis_a, basis_b; xmin = -10.0, xmax = 10.0, h = 0.05)
    points = collect((xmin + 0.5 * h):h:(xmax - 0.5 * h))
    weights_sqrt = sqrt.(fill(h, length(points)))
    values_a = _basis_sample_matrix(basis_a, points) .* weights_sqrt
    values_b = _basis_sample_matrix(basis_b, points) .* weights_sqrt
    qa = Matrix(qr(values_a).Q[:, 1:size(values_a, 2)])
    qb = Matrix(qr(values_b).Q[:, 1:size(values_b, 2)])
    singular_values = svdvals(transpose(qa) * qb)
    min_sv = minimum(singular_values)
    return (min_sv = min_sv, max_angle = acos(clamp(min_sv, -1.0, 1.0)))
end

function _projector_difference_metric(basis_a, basis_b; xmin = -10.0, xmax = 10.0, h = 0.05)
    points = collect((xmin + 0.5 * h):h:(xmax - 0.5 * h))
    weights_sqrt = sqrt.(fill(h, length(points)))
    values_a = _basis_sample_matrix(basis_a, points) .* weights_sqrt
    values_b = _basis_sample_matrix(basis_b, points) .* weights_sqrt
    qa = Matrix(qr(values_a).Q[:, 1:size(values_a, 2)])
    qb = Matrix(qr(values_b).Q[:, 1:size(values_b, 2)])
    pa = qa * transpose(qa)
    pb = qb * transpose(qb)
    return (frob = norm(pa - pb), op = opnorm(pa - pb))
end

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

if _test_group_enabled(:radial)
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
end

function _quick_radial_operator_fixture()
    return _cached_fixture(:quick_radial_operator_fixture, () -> begin
        _radial_operator_fixture()
    end)
end

function _quick_cartesian_hydrogen_fixture()
    return _cached_fixture(:quick_cartesian_hydrogen_fixture, () -> begin
        Z = 1.0
        basis = build_basis(UniformBasisSpec(:G10; xmin = -2.0, xmax = 2.0, spacing = 1.0))
        representation = basis_representation(basis; operators = (:overlap, :kinetic))
        overlap_1d = representation.basis_matrices.overlap
        kinetic_1d = representation.basis_matrices.kinetic
        expansion = coulomb_gaussian_expansion(doacc = false)
        gaussian_factors = gaussian_factor_matrices(
            basis;
            exponents = expansion.exponents,
            center = 0.0,
        )
        overlap_3d = kron(overlap_1d, kron(overlap_1d, overlap_1d))
        kinetic_3d =
            kron(kinetic_1d, kron(overlap_1d, overlap_1d)) +
            kron(overlap_1d, kron(kinetic_1d, overlap_1d)) +
            kron(overlap_1d, kron(overlap_1d, kinetic_1d))
        nuclear_3d = zeros(Float64, size(overlap_3d))
        for term in eachindex(expansion.coefficients)
            factor = gaussian_factors[term]
            nuclear_3d .-= Z * expansion.coefficients[term] .* kron(factor, kron(factor, factor))
        end
        hamiltonian = kinetic_3d + nuclear_3d
        energy = minimum(eigen(Hermitian(hamiltonian)).values)
        (
            basis,
            representation,
            expansion,
            overlap_1d,
            kinetic_1d,
            gaussian_factors,
            overlap_3d,
            nuclear_3d,
            hamiltonian,
            energy,
        )
    end)
end

function _quick_mapped_cartesian_hydrogen_fixture()
    return _cached_fixture(:quick_mapped_cartesian_hydrogen_fixture, () -> begin
        Z = 1.0
        mapping = fit_asinh_mapping_for_extent(npoints = 5, xmax = 6.0)
        basis = build_basis(MappedUniformBasisSpec(:G10;
            count = 5,
            mapping = mapping,
            reference_spacing = 1.0,
        ))
        representation = basis_representation(basis; operators = (:overlap, :kinetic))
        overlap_1d = representation.basis_matrices.overlap
        kinetic_1d = representation.basis_matrices.kinetic
        expansion = coulomb_gaussian_expansion(doacc = false)
        gaussian_factors = gaussian_factor_matrices(
            basis;
            exponents = expansion.exponents,
            center = 0.0,
        )
        transformed = _orthonormalize_1d(overlap_1d, [kinetic_1d, gaussian_factors...])
        kinetic_orth = first(transformed)
        gaussian_orth = collect(transformed[2:end])
        identity_1d = Matrix{Float64}(I, length(basis), length(basis))
        hamiltonian =
            kron(kinetic_orth, kron(identity_1d, identity_1d)) +
            kron(identity_1d, kron(kinetic_orth, identity_1d)) +
            kron(identity_1d, kron(identity_1d, kinetic_orth))
        for term in eachindex(expansion.coefficients)
            factor = gaussian_orth[term]
            hamiltonian .-= Z * expansion.coefficients[term] .* kron(factor, kron(factor, factor))
        end
        energy = minimum(eigen(Hermitian(hamiltonian)).values)
        (
            basis,
            mapping,
            representation,
            overlap_1d,
            kinetic_1d,
            expansion,
            hamiltonian,
            energy,
        )
    end)
end

function _truncate_coulomb_expansion(expansion::CoulombGaussianExpansion, nterms::Int)
    return CoulombGaussianExpansion(
        expansion.coefficients[1:nterms],
        expansion.exponents[1:nterms];
        del = expansion.del,
        s = expansion.s,
        c = expansion.c,
        maxu = expansion.maxu,
    )
end

function _quick_ordinary_cartesian_ida_fixture(; backend = :pgdg_experimental, mapped = true, s = 0.5)
    key = Symbol(:ordinary_cartesian_ida, backend, mapped ? :mapped : :identity, s)
    return _cached_fixture(key, () -> begin
        full_expansion = coulomb_gaussian_expansion(doacc = false)
        expansion = _truncate_coulomb_expansion(full_expansion, 3)
        mapping_value = mapped ?
            fit_asinh_mapping_for_strength(s = s, npoints = 5, xmax = 6.0) :
            IdentityMapping()
        basis = build_basis(MappedUniformBasisSpec(:G10;
            count = 5,
            mapping = mapping_value,
            reference_spacing = 1.0,
        ))
        operators = ordinary_cartesian_ida_operators(
            basis;
            expansion = expansion,
            Z = 2.0,
            backend = backend,
        )
        (basis, expansion, operators)
    end)
end

function _quick_ordinary_cartesian_1s2_vee_fixture()
    return _cached_fixture(:ordinary_cartesian_1s2_vee, () -> begin
        basis = build_basis(MappedUniformBasisSpec(:G10;
            count = 9,
            mapping = IdentityMapping(),
            reference_spacing = 0.5,
        ))
        operators = ordinary_cartesian_ida_operators(
            basis;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
            backend = :numerical_reference,
        )
        decomposition = eigen(Hermitian(operators.one_body_hamiltonian))
        orbital_energy = decomposition.values[1]
        orbital = decomposition.vectors[:, 1]
        vee = ordinary_cartesian_vee_expectation(operators, orbital)
        (basis, operators, orbital_energy, orbital, vee)
    end)
end

function _quick_hybrid_cartesian_1s2_vee_fixture()
    return _cached_fixture(:hybrid_cartesian_1s2_vee, () -> begin
        source_basis = build_basis(MappedUniformBasisSpec(:G10;
            count = 7,
            mapping = fit_asinh_mapping_for_strength(s = 0.2, npoints = 7, xmax = 6.0),
            reference_spacing = 1.0,
        ))
        core_gaussians = [
            Gaussian(center = 0.0, width = 0.2),
            Gaussian(center = 0.0, width = 0.6),
        ]
        pure_operators = ordinary_cartesian_ida_operators(
            source_basis;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
            backend = :pgdg_localized_experimental,
        )
        hybrid_basis = hybrid_mapped_ordinary_basis(
            source_basis;
            core_gaussians = core_gaussians,
            backend = :pgdg_localized_experimental,
        )
        hybrid_operators = ordinary_cartesian_ida_operators(
            hybrid_basis;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
        )
        pure_check = GaussletBases.ordinary_cartesian_1s2_check(pure_operators)
        hybrid_check = GaussletBases.ordinary_cartesian_1s2_check(hybrid_operators)
        (
            source_basis,
            core_gaussians,
            pure_operators,
            hybrid_basis,
            hybrid_operators,
            pure_check,
            hybrid_check,
        )
    end)
end

function _friendly_hybrid_residual_vee_fixture(count::Int, s::Float64)
    key = Symbol(:friendly_hybrid_residual_vee, count, round(Int, 1000 * s))
    return _cached_fixture(key, () -> begin
        source_basis = build_basis(MappedUniformBasisSpec(:G10;
            count = count,
            mapping = fit_asinh_mapping_for_strength(s = s, npoints = count, xmax = 6.0),
            reference_spacing = 1.0,
        ))
        core_gaussians = [
            Gaussian(center = 0.0, width = 0.2),
            Gaussian(center = 0.0, width = 0.6),
        ]
        pure_operators = ordinary_cartesian_ida_operators(
            source_basis;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
            backend = :pgdg_localized_experimental,
        )
        hybrid_basis = hybrid_mapped_ordinary_basis(
            source_basis;
            core_gaussians = core_gaussians,
            backend = :pgdg_localized_experimental,
        )
        combined_operators = ordinary_cartesian_ida_operators(
            hybrid_basis;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
        )
        residual_operators = ordinary_cartesian_ida_operators(
            hybrid_basis;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
            interaction_treatment = :residual_gaussian_nearest,
        )
        (
            source_basis,
            hybrid_basis,
            pure_operators,
            combined_operators,
            residual_operators,
            GaussletBases.ordinary_cartesian_1s2_check(pure_operators),
            GaussletBases.ordinary_cartesian_1s2_check(combined_operators),
            GaussletBases.ordinary_cartesian_1s2_check(residual_operators),
        )
    end)
end

function _quick_hybrid_mapped_ordinary_fixture()
    return _cached_fixture(:quick_hybrid_mapped_ordinary_fixture, () -> begin
        full_expansion = coulomb_gaussian_expansion(doacc = false)
        expansion = _truncate_coulomb_expansion(full_expansion, 3)
        basis = build_basis(MappedUniformBasisSpec(:G10;
            count = 5,
            mapping = fit_asinh_mapping_for_strength(s = 0.2, npoints = 5, xmax = 6.0),
            reference_spacing = 1.0,
        ))
        core_gaussians = [
            Gaussian(center = 0.0, width = 0.2),
            Gaussian(center = 0.0, width = 0.6),
        ]
        hybrid_reference = hybrid_mapped_ordinary_basis(
            basis;
            core_gaussians = core_gaussians,
            backend = :numerical_reference,
        )
        hybrid_analytic = hybrid_mapped_ordinary_basis(
            basis;
            core_gaussians = core_gaussians,
            backend = :pgdg_localized_experimental,
        )
        one_body_reference = mapped_ordinary_one_body_operators(
            hybrid_reference;
            exponents = expansion.exponents,
        )
        one_body_analytic = mapped_ordinary_one_body_operators(
            hybrid_analytic;
            exponents = expansion.exponents,
        )
        hard_basis = build_basis(MappedUniformBasisSpec(:G10;
            count = 5,
            mapping = fit_asinh_mapping_for_strength(s = 2.0, npoints = 5, xmax = 6.0),
            reference_spacing = 1.0,
        ))
        hard_reference = mapped_ordinary_one_body_operators(
            hard_basis;
            exponents = expansion.exponents,
            backend = :numerical_reference,
        )
        hard_analytic = mapped_ordinary_one_body_operators(
            hard_basis;
            exponents = expansion.exponents,
            backend = :pgdg_localized_experimental,
        )
        energy_reference = mapped_cartesian_hydrogen_energy(
            hybrid_reference;
            expansion = expansion,
            Z = 1.0,
        )
        energy_analytic = mapped_cartesian_hydrogen_energy(
            hybrid_analytic;
            expansion = expansion,
            Z = 1.0,
        )
        hard_energy_reference = mapped_cartesian_hydrogen_energy(
            hard_basis;
            expansion = expansion,
            Z = 1.0,
            backend = :numerical_reference,
        )
        hard_energy_analytic = mapped_cartesian_hydrogen_energy(
            hard_basis;
            expansion = expansion,
            Z = 1.0,
            backend = :pgdg_localized_experimental,
        )
        (
            basis,
            core_gaussians,
            expansion,
            hybrid_reference,
            hybrid_analytic,
            one_body_reference,
            one_body_analytic,
            hard_reference,
            hard_analytic,
            energy_reference,
            energy_analytic,
            hard_energy_reference,
            hard_energy_analytic,
        )
    end)
end

_legacy_basisfile_path_for_tests() = GaussletBases._legacy_basisfile_path()
_legacy_basisfile_available() = isfile(_legacy_basisfile_path_for_tests())

function _legacy_he_s_hybrid_fixture(basis_name::String)
    key = Symbol(:legacy_he_s_hybrid_fixture, Symbol(lowercase(basis_name)))
    return _cached_fixture(key, () -> begin
        source_basis = build_basis(MappedUniformBasisSpec(:G10;
            count = 11,
            mapping = fit_asinh_mapping_for_strength(s = 0.6, npoints = 11, xmax = 6.0),
            reference_spacing = 1.0,
        ))
        expansion = coulomb_gaussian_expansion(doacc = false)
        pure_operators = ordinary_cartesian_ida_operators(
            source_basis;
            expansion = expansion,
            Z = 2.0,
            backend = :pgdg_localized_experimental,
        )
        toy_basis = hybrid_mapped_ordinary_basis(
            source_basis;
            core_gaussians = [Gaussian(center = 0.0, width = 0.2), Gaussian(center = 0.0, width = 0.6)],
            backend = :pgdg_localized_experimental,
        )
        toy_operators = ordinary_cartesian_ida_operators(
            toy_basis;
            expansion = expansion,
            Z = 2.0,
        )
        legacy = legacy_atomic_gaussian_supplement("He", basis_name; lmax = 0)
        legacy_basis = hybrid_mapped_ordinary_basis(
            source_basis;
            core_gaussians = legacy,
            backend = :pgdg_localized_experimental,
        )
        legacy_operators = ordinary_cartesian_ida_operators(
            legacy_basis;
            expansion = expansion,
            Z = 2.0,
        )
        (
            source_basis,
            legacy,
            GaussletBases.ordinary_cartesian_1s2_check(pure_operators),
            GaussletBases.ordinary_cartesian_1s2_check(toy_operators),
            GaussletBases.ordinary_cartesian_1s2_check(legacy_operators),
        )
    end)
end

function _legacy_he_s_mwg_fixture(basis_name::String; s::Float64 = 0.6)
    key = Symbol(:legacy_he_s_mwg_fixture, Symbol(lowercase(basis_name)), round(Int, 1000 * s))
    return _cached_fixture(key, () -> begin
        source_basis = build_basis(MappedUniformBasisSpec(:G10;
            count = 11,
            mapping = fit_asinh_mapping_for_strength(s = s, npoints = 11, xmax = 6.0),
            reference_spacing = 1.0,
        ))
        expansion = coulomb_gaussian_expansion(doacc = false)
        pure_operators = ordinary_cartesian_ida_operators(
            source_basis;
            expansion = expansion,
            Z = 2.0,
            backend = :pgdg_localized_experimental,
        )
        legacy = legacy_atomic_gaussian_supplement("He", basis_name; lmax = 0)
        hybrid_basis = hybrid_mapped_ordinary_basis(
            source_basis;
            core_gaussians = legacy,
            backend = :pgdg_localized_experimental,
        )
        combined_operators = ordinary_cartesian_ida_operators(
            hybrid_basis;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :combined_basis,
        )
        nearest_operators = ordinary_cartesian_ida_operators(
            hybrid_basis;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :residual_gaussian_nearest,
        )
        mwg_operators = ordinary_cartesian_ida_operators(
            hybrid_basis;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :residual_gaussian_mwg,
        )
        mwg_data = GaussletBases._hybrid_residual_gaussian_mwg_data(hybrid_basis)
        (
            source_basis,
            legacy,
            hybrid_basis,
            pure_operators,
            combined_operators,
            nearest_operators,
            mwg_operators,
            GaussletBases.ordinary_cartesian_1s2_check(pure_operators),
            GaussletBases.ordinary_cartesian_1s2_check(combined_operators),
            GaussletBases.ordinary_cartesian_1s2_check(nearest_operators),
            GaussletBases.ordinary_cartesian_1s2_check(mwg_operators),
            mwg_data,
        )
    end)
end

function _qiu_white_reference_fixture(; basis_name::String = "cc-pVTZ", count::Int = 9, s::Float64 = 0.8, nterms::Int = 3)
    key = Symbol(:qiu_white_reference_fixture, Symbol(lowercase(basis_name)), count, round(Int, 1000 * s), nterms)
    return _cached_fixture(key, () -> begin
        source_basis = build_basis(MappedUniformBasisSpec(:G10;
            count = count,
            mapping = fit_asinh_mapping_for_strength(s = s, npoints = count, xmax = 6.0),
            reference_spacing = 1.0,
        ))
        expansion = _truncate_coulomb_expansion(coulomb_gaussian_expansion(doacc = false), nterms)
        legacy = legacy_atomic_gaussian_supplement("He", basis_name; lmax = 0)
        hybrid_basis = hybrid_mapped_ordinary_basis(
            source_basis;
            core_gaussians = legacy,
            backend = :pgdg_localized_experimental,
        )
        surrogate_mwg = ordinary_cartesian_ida_operators(
            hybrid_basis;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :residual_gaussian_mwg,
        )
        qiu_nearest = ordinary_cartesian_qiu_white_operators(
            source_basis,
            legacy;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        qiu_mwg = ordinary_cartesian_qiu_white_operators(
            source_basis,
            legacy;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :mwg,
        )
        (
            source_basis,
            legacy,
            surrogate_mwg,
            qiu_nearest,
            qiu_mwg,
            GaussletBases.ordinary_cartesian_1s2_check(surrogate_mwg),
            GaussletBases.ordinary_cartesian_1s2_check(qiu_nearest),
            GaussletBases.ordinary_cartesian_1s2_check(qiu_mwg),
        )
    end)
end

function _qiu_white_full_nearest_fixture(; basis_name::String = "cc-pVTZ", count::Int = 9, s::Float64 = 0.8)
    key = Symbol(:qiu_white_full_nearest_fixture, Symbol(lowercase(basis_name)), count, round(Int, 1000 * s))
    return _cached_fixture(key, () -> begin
        source_basis = build_basis(MappedUniformBasisSpec(:G10;
            count = count,
            mapping = fit_asinh_mapping_for_strength(s = s, npoints = count, xmax = 6.0),
            reference_spacing = 1.0,
        ))
        legacy = legacy_atomic_gaussian_supplement("He", basis_name; lmax = 0)
        operators = ordinary_cartesian_qiu_white_operators(
            source_basis,
            legacy;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        (
            source_basis,
            legacy,
            operators,
            GaussletBases.ordinary_cartesian_1s2_check(operators),
        )
    end)
end

function _nested_qiu_white_nearest_fixture(; basis_name::String = "cc-pVTZ", count::Int = 13, a::Float64 = 0.25, xmax::Float64 = 10.0, tail_spacing::Float64 = 10.0)
    key = Symbol(
        :nested_qiu_white_nearest_fixture,
        Symbol(lowercase(basis_name)),
        count,
        round(Int, 1000 * a),
        round(Int, 1000 * xmax),
    )
    return _cached_fixture(key, () -> begin
        endpoint = (count - 1) / 2
        s = asinh(xmax / a) / (endpoint - xmax / tail_spacing)
        source_basis = build_basis(MappedUniformBasisSpec(:G10;
            count = count,
            mapping = AsinhMapping(a = a, s = s, tail_spacing = tail_spacing),
            reference_spacing = 1.0,
        ))
        expansion = coulomb_gaussian_expansion(doacc = false)
        term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
        bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
            source_basis;
            exponents = expansion.exponents,
            backend = :numerical_reference,
            refinement_levels = 0,
        )
        interval = 2:(length(source_basis) - 1)
        shell = GaussletBases._nested_rectangular_shell(
            bundle,
            interval,
            interval,
            interval;
            retain_xy = (4, 3),
            retain_xz = (4, 3),
            retain_yz = (4, 3),
            term_coefficients = term_coefficients,
        )
        fixed_block = GaussletBases._nested_fixed_block(shell, bundle)
        shell_plus_core = GaussletBases._nested_shell_plus_core(
            bundle,
            shell,
            interval,
            interval,
            interval,
            term_coefficients = term_coefficients,
        )
        fixed_block_shell_plus_core = GaussletBases._nested_fixed_block(shell_plus_core, bundle)
        legacy = legacy_atomic_gaussian_supplement("He", basis_name; lmax = 0)
        baseline = ordinary_cartesian_qiu_white_operators(
            source_basis,
            legacy;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        nested = ordinary_cartesian_qiu_white_operators(
            fixed_block,
            legacy;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        nested_shell_plus_core = ordinary_cartesian_qiu_white_operators(
            fixed_block_shell_plus_core,
            legacy;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        (
            source_basis,
            bundle,
            shell,
            fixed_block,
            shell_plus_core,
            fixed_block_shell_plus_core,
            legacy,
            baseline,
            nested,
            nested_shell_plus_core,
            GaussletBases.ordinary_cartesian_1s2_check(baseline),
            GaussletBases.ordinary_cartesian_1s2_check(nested),
            GaussletBases.ordinary_cartesian_1s2_check(nested_shell_plus_core),
        )
    end)
end

function _nested_qiu_white_shell_sequence_fixture(; basis_name::String = "cc-pVTZ", count::Int = 17, a::Float64 = 0.25, xmax::Float64 = 10.0, tail_spacing::Float64 = 10.0)
    key = Symbol(
        :nested_qiu_white_shell_sequence_fixture,
        Symbol(lowercase(basis_name)),
        count,
        round(Int, 1000 * a),
        round(Int, 1000 * xmax),
    )
    return _cached_fixture(key, () -> begin
        endpoint = (count - 1) / 2
        s = asinh(xmax / a) / (endpoint - xmax / tail_spacing)
        source_basis = build_basis(MappedUniformBasisSpec(:G10;
            count = count,
            mapping = AsinhMapping(a = a, s = s, tail_spacing = tail_spacing),
            reference_spacing = 1.0,
        ))
        expansion = coulomb_gaussian_expansion(doacc = false)
        term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
        bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
            source_basis;
            exponents = expansion.exponents,
            backend = :numerical_reference,
            refinement_levels = 0,
        )
        core_interval = 4:14
        shell1 = GaussletBases._nested_rectangular_shell(
            bundle,
            core_interval,
            core_interval,
            core_interval;
            retain_xy = (4, 3),
            retain_xz = (4, 3),
            retain_yz = (4, 3),
            x_fixed = (3, 15),
            y_fixed = (3, 15),
            z_fixed = (3, 15),
            term_coefficients = term_coefficients,
        )
        shell_plus_core = GaussletBases._nested_shell_plus_core(
            bundle,
            shell1,
            core_interval,
            core_interval,
            core_interval,
            term_coefficients = term_coefficients,
        )
        fixed_shell_plus_core = GaussletBases._nested_fixed_block(shell_plus_core, bundle)
        outer_interval = 3:15
        shell2 = GaussletBases._nested_rectangular_shell(
            bundle,
            outer_interval,
            outer_interval,
            outer_interval;
            retain_xy = (4, 3),
            retain_xz = (4, 3),
            retain_yz = (4, 3),
            x_fixed = (2, 16),
            y_fixed = (2, 16),
            z_fixed = (2, 16),
            term_coefficients = term_coefficients,
        )
        shell_sequence = GaussletBases._nested_shell_sequence(
            bundle,
            core_interval,
            core_interval,
            core_interval,
            [shell1, shell2];
            enforce_coverage = false,
            term_coefficients = term_coefficients,
        )
        fixed_sequence = GaussletBases._nested_fixed_block(shell_sequence, bundle)
        legacy = legacy_atomic_gaussian_supplement("He", basis_name; lmax = 0)
        baseline = ordinary_cartesian_qiu_white_operators(
            source_basis,
            legacy;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        shell_plus_core_ops = ordinary_cartesian_qiu_white_operators(
            fixed_shell_plus_core,
            legacy;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        shell_sequence_ops = ordinary_cartesian_qiu_white_operators(
            fixed_sequence,
            legacy;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        (
            source_basis,
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
            GaussletBases.ordinary_cartesian_1s2_check(baseline),
            GaussletBases.ordinary_cartesian_1s2_check(shell_plus_core_ops),
            GaussletBases.ordinary_cartesian_1s2_check(shell_sequence_ops),
        )
    end)
end

function _nested_qiu_white_nside_sequence_fixture(; basis_name::String = "cc-pVTZ", count::Int = 17, a::Float64 = 0.25, xmax::Float64 = 10.0, tail_spacing::Float64 = 10.0, nside::Int = 5)
    key = Symbol(
        :nested_qiu_white_nside_sequence_fixture,
        Symbol(lowercase(basis_name)),
        count,
        round(Int, 1000 * a),
        round(Int, 1000 * xmax),
        nside,
    )
    return _cached_fixture(key, () -> begin
        endpoint = (count - 1) / 2
        s = asinh(xmax / a) / (endpoint - xmax / tail_spacing)
        source_basis = build_basis(MappedUniformBasisSpec(:G10;
            count = count,
            mapping = AsinhMapping(a = a, s = s, tail_spacing = tail_spacing),
            reference_spacing = 1.0,
        ))
        expansion = coulomb_gaussian_expansion(doacc = false)
        term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
        bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
            source_basis;
            exponents = expansion.exponents,
            backend = :numerical_reference,
            refinement_levels = 0,
        )
        core_interval = 4:14
        shell1 = GaussletBases._nested_rectangular_shell(
            bundle,
            core_interval,
            core_interval,
            core_interval;
            retain_xy = (4, 3),
            retain_xz = (4, 3),
            retain_yz = (4, 3),
            x_fixed = (3, 15),
            y_fixed = (3, 15),
            z_fixed = (3, 15),
            term_coefficients = term_coefficients,
        )
        middle_interval = 5:13
        shell2 = GaussletBases._nested_rectangular_shell(
            bundle,
            middle_interval,
            middle_interval,
            middle_interval;
            retain_xy = (4, 3),
            retain_xz = (4, 3),
            retain_yz = (4, 3),
            x_fixed = (4, 14),
            y_fixed = (4, 14),
            z_fixed = (4, 14),
            term_coefficients = term_coefficients,
        )
        inner_shell_interval = 6:12
        shell3 = GaussletBases._nested_rectangular_shell(
            bundle,
            inner_shell_interval,
            inner_shell_interval,
            inner_shell_interval;
            retain_xy = (4, 3),
            retain_xz = (4, 3),
            retain_yz = (4, 3),
            x_fixed = (5, 13),
            y_fixed = (5, 13),
            z_fixed = (5, 13),
            term_coefficients = term_coefficients,
        )
        inner_direct_interval = 7:11
        grow_sequence = GaussletBases._nested_shell_sequence(
            bundle,
            inner_direct_interval,
            inner_direct_interval,
            inner_direct_interval,
            [shell1, shell2, shell3];
            enforce_coverage = false,
            term_coefficients = term_coefficients,
        )
        shrinking_sequence = GaussletBases._nested_nside_shell_sequence(
            bundle,
            core_interval,
            core_interval,
            core_interval,
            [shell1, shell2, shell3];
            nside = nside,
            term_coefficients = term_coefficients,
        )
        fixed_grow = GaussletBases._nested_fixed_block(grow_sequence, bundle)
        fixed_shrink = GaussletBases._nested_fixed_block(shrinking_sequence, bundle)
        legacy = legacy_atomic_gaussian_supplement("He", basis_name; lmax = 0)
        baseline = ordinary_cartesian_qiu_white_operators(
            source_basis,
            legacy;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        grow_ops = ordinary_cartesian_qiu_white_operators(
            fixed_grow,
            legacy;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        shrink_ops = ordinary_cartesian_qiu_white_operators(
            fixed_shrink,
            legacy;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        (
            source_basis,
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
            GaussletBases.ordinary_cartesian_1s2_check(baseline),
            GaussletBases.ordinary_cartesian_1s2_check(grow_ops),
            GaussletBases.ordinary_cartesian_1s2_check(shrink_ops),
        )
    end)
end

function _bond_aligned_diatomic_qw_fixture(; bond_length::Float64 = 1.4)
    key = Symbol(:bond_aligned_diatomic_qw_fixture, round(Int, 1000 * bond_length))
    return _cached_fixture(key, () -> begin
        basis = bond_aligned_homonuclear_qw_basis(
            bond_length = bond_length,
            core_spacing = 0.5,
            xmax_parallel = 6.0,
            xmax_transverse = 4.0,
            bond_axis = :z,
        )
        operators = ordinary_cartesian_qiu_white_operators(
            basis;
            nuclear_charges = [1.0, 1.0],
            interaction_treatment = :ggt_nearest,
        )
        (
            basis,
            operators,
            GaussletBases.ordinary_cartesian_1s2_check(operators),
        )
    end)
end

function _bond_aligned_homonuclear_chain_qw_fixture(;
    natoms::Int = 3,
    spacing::Float64 = 1.2,
    chain_axis::Symbol = :z,
)
    key = Symbol(
        :bond_aligned_homonuclear_chain_qw_fixture,
        natoms,
        round(Int, 1000 * spacing),
        chain_axis,
    )
    return _cached_fixture(key, () -> begin
        basis = bond_aligned_homonuclear_chain_qw_basis(
            natoms = natoms,
            spacing = spacing,
            core_spacing = 0.5,
            xmax_parallel = 2.0,
            xmax_transverse = 2.0,
            chain_axis = chain_axis,
        )
        operators = ordinary_cartesian_qiu_white_operators(
            basis;
            nuclear_charges = fill(1.0, natoms),
            interaction_treatment = :ggt_nearest,
        )
        (
            basis,
            operators,
            bond_aligned_homonuclear_chain_geometry_diagnostics(basis),
        )
    end)
end

function _axis_aligned_homonuclear_square_lattice_qw_fixture(;
    n::Int = 2,
    spacing::Float64 = 1.2,
)
    key = Symbol(
        :axis_aligned_homonuclear_square_lattice_qw_fixture,
        n,
        round(Int, 1000 * spacing),
    )
    return _cached_fixture(key, () -> begin
        basis = axis_aligned_homonuclear_square_lattice_qw_basis(
            n = n,
            spacing = spacing,
            core_spacing = 0.5,
            xmax_in_plane = 2.0,
            xmax_transverse = 2.0,
        )
        operators = ordinary_cartesian_qiu_white_operators(
            basis;
            nuclear_charges = fill(1.0, length(basis.nuclei)),
            interaction_treatment = :ggt_nearest,
        )
        (
            basis,
            operators,
            axis_aligned_homonuclear_square_lattice_geometry_diagnostics(basis),
            GaussletBases.ordinary_cartesian_1s2_check(operators),
        )
    end)
end

function _axis_aligned_homonuclear_square_lattice_nested_fixture(;
    n::Int = 2,
    spacing::Float64 = 1.2,
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
)
    key = Symbol(
        :axis_aligned_homonuclear_square_lattice_nested_fixture,
        n,
        round(Int, 1000 * spacing),
        nside,
        round(Int, 1000 * min_in_plane_aspect_ratio),
    )
    return _cached_fixture(key, () -> begin
        basis = axis_aligned_homonuclear_square_lattice_qw_basis(
            n = n,
            spacing = spacing,
            core_spacing = 0.5,
            xmax_in_plane = 2.0,
            xmax_transverse = 2.0,
        )
        source = GaussletBases._axis_aligned_homonuclear_square_lattice_nested_fixed_source(
            basis;
            nside = nside,
            min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
        )
        (
            basis,
            source,
            GaussletBases._nested_fixed_block(source),
            axis_aligned_homonuclear_square_lattice_nested_geometry_diagnostics(
                basis;
                nside = nside,
                min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
            ),
        )
    end)
end

function _axis_aligned_homonuclear_square_lattice_nested_qw_fixture(;
    n::Int = 2,
    spacing::Float64 = 1.2,
    nside::Int = 5,
    min_in_plane_aspect_ratio::Float64 = 0.15,
)
    key = Symbol(
        :axis_aligned_homonuclear_square_lattice_nested_qw_fixture,
        n,
        round(Int, 1000 * spacing),
        nside,
        round(Int, 1000 * min_in_plane_aspect_ratio),
    )
    return _cached_fixture(key, () -> begin
        basis = axis_aligned_homonuclear_square_lattice_qw_basis(
            n = n,
            spacing = spacing,
            core_spacing = 0.5,
            xmax_in_plane = 2.0,
            xmax_transverse = 2.0,
        )
        path = experimental_axis_aligned_homonuclear_square_lattice_nested_qw_operators(
            basis;
            nuclear_charges = fill(1.0, length(basis.nuclei)),
            nside = nside,
            min_in_plane_aspect_ratio = min_in_plane_aspect_ratio,
        )
        return (
            basis,
            path,
            GaussletBases.ordinary_cartesian_1s2_check(path.operators),
        )
    end)
end

function _bond_aligned_homonuclear_chain_nested_fixture(;
    natoms::Int = 4,
    spacing::Float64 = 1.2,
    chain_axis::Symbol = :z,
    nside::Int = 5,
    odd_chain_policy::Symbol = :strict_current,
)
    key = Symbol(
        :bond_aligned_homonuclear_chain_nested_fixture,
        natoms,
        round(Int, 1000 * spacing),
        chain_axis,
        nside,
        odd_chain_policy,
    )
    return _cached_fixture(key, () -> begin
        basis = bond_aligned_homonuclear_chain_qw_basis(
            natoms = natoms,
            spacing = spacing,
            core_spacing = 0.5,
            xmax_parallel = 2.0,
            xmax_transverse = 2.0,
            chain_axis = chain_axis,
        )
        source = GaussletBases._bond_aligned_homonuclear_chain_nested_fixed_source(
            basis;
            nside = nside,
            odd_chain_policy = odd_chain_policy,
        )
        (
            basis,
            source,
            GaussletBases._nested_fixed_block(source),
            bond_aligned_homonuclear_chain_nested_geometry_diagnostics(
                basis;
                nside = nside,
                odd_chain_policy = odd_chain_policy,
            ),
        )
    end)
end

function _bond_aligned_homonuclear_chain_nested_qw_fixture(;
    natoms::Int = 3,
    spacing::Float64 = 1.2,
    chain_axis::Symbol = :z,
    nside::Int = 5,
    odd_chain_policy::Symbol = :central_ternary_relaxed,
)
    key = Symbol(
        :bond_aligned_homonuclear_chain_nested_qw_fixture,
        natoms,
        round(Int, 1000 * spacing),
        chain_axis,
        nside,
        odd_chain_policy,
    )
    return _cached_fixture(key, () -> begin
        basis = bond_aligned_homonuclear_chain_qw_basis(
            natoms = natoms,
            spacing = spacing,
            core_spacing = 0.5,
            xmax_parallel = 2.0,
            xmax_transverse = 2.0,
            chain_axis = chain_axis,
        )
        path = experimental_bond_aligned_homonuclear_chain_nested_qw_operators(
            basis;
            nuclear_charges = fill(1.0, natoms),
            nside = nside,
            odd_chain_policy = odd_chain_policy,
        )
        return (
            basis,
            path,
            GaussletBases.ordinary_cartesian_1s2_check(path.operators),
        )
    end)
end

function _bond_aligned_diatomic_nested_fixed_block_fixture(; bond_length::Float64 = 1.4)
    key = Symbol(:bond_aligned_diatomic_nested_fixed_block_fixture, round(Int, 1000 * bond_length))
    return _cached_fixture(key, () -> begin
        basis, operators, check = _bond_aligned_diatomic_qw_fixture(; bond_length = bond_length)
        expansion = coulomb_gaussian_expansion(doacc = false)
        nested = bond_aligned_diatomic_nested_fixed_block(
            basis;
            expansion = expansion,
        )
        source = nested.source
        fixed_block = nested.fixed_block
        parent_modes = eigen(Hermitian(operators.one_body_hamiltonian), Hermitian(operators.overlap))
        parent_ground = parent_modes.vectors[:, 1]
        projected = _nested_fixed_projected_orbital(operators.overlap, fixed_block, parent_ground)
        projected_vee = _nested_vee_from_orbital(
            GaussletBases._qwrg_fixed_block_interaction_matrix(fixed_block, expansion),
            projected,
        )
        capture, projected_energy = _nested_projector_stats(
            operators.overlap,
            operators.one_body_hamiltonian,
            fixed_block,
            parent_ground,
        )
        (
            basis,
            operators,
            check,
            expansion,
            source,
            fixed_block,
            parent_modes,
            parent_ground,
            projected,
            projected_vee,
            capture,
            projected_energy,
        )
    end)
end

function _bond_aligned_diatomic_nested_qw_fixture(; bond_length::Float64 = 1.4)
    key = Symbol(:bond_aligned_diatomic_nested_qw_fixture, round(Int, 1000 * bond_length))
    return _cached_fixture(key, () -> begin
        (
            basis,
            parent_ops,
            parent_check,
            expansion,
            source,
            fixed_block,
            _parent_modes,
            _parent_ground,
            _projected,
            _projected_vee,
            _capture,
            _projected_energy,
        ) = _bond_aligned_diatomic_nested_fixed_block_fixture(; bond_length = bond_length)
        nested_ops = ordinary_cartesian_qiu_white_operators(
            fixed_block;
            nuclear_charges = [1.0, 1.0],
            expansion = expansion,
            interaction_treatment = :ggt_nearest,
        )
        (
            basis,
            parent_ops,
            parent_check,
            expansion,
            source,
            fixed_block,
            nested_ops,
            GaussletBases.ordinary_cartesian_1s2_check(nested_ops),
        )
    end)
end

function _bond_aligned_diatomic_hybrid_qw_fixture(; bond_length::Float64 = 1.4)
    key = Symbol(:bond_aligned_diatomic_hybrid_qw_fixture, round(Int, 1000 * bond_length))
    return _cached_fixture(key, () -> begin
        basis, parent_ops, parent_check = _bond_aligned_diatomic_qw_fixture(; bond_length = bond_length)
        supplement = legacy_bond_aligned_diatomic_gaussian_supplement(
            "H",
            "cc-pVTZ",
            basis.nuclei;
            lmax = 1,
        )
        ordinary_ops = ordinary_cartesian_qiu_white_operators(
            basis,
            supplement;
            nuclear_charges = [1.0, 1.0],
            interaction_treatment = :ggt_nearest,
        )
        (
            basis,
            parent_ops,
            parent_check,
            supplement,
            ordinary_ops,
            GaussletBases.ordinary_cartesian_1s2_check(ordinary_ops),
        )
    end)
end

function _bond_aligned_diatomic_nested_hybrid_bundle_fixture(;
    atom::String = "H",
    basis_name::String = "cc-pVTZ",
    bond_length::Float64 = 1.4,
    core_spacing::Float64 = 0.5,
    xmax_parallel::Float64 = 8.0,
    xmax_transverse::Float64 = 5.0,
    bond_axis::Symbol = :z,
    nuclear_charge::Float64 = 1.0,
    nside::Int = 5,
    supplement_lmax::Int = 1,
    max_width::Union{Nothing,Float64} = nothing,
)
    key = Symbol(
        :bond_aligned_diatomic_nested_hybrid_bundle_fixture,
        atom,
        basis_name,
        round(Int, 1000 * bond_length),
        round(Int, 1000 * core_spacing),
        round(Int, 1000 * xmax_parallel),
        round(Int, 1000 * xmax_transverse),
        bond_axis,
        round(Int, 1000 * nuclear_charge),
        nside,
        supplement_lmax,
        isnothing(max_width) ? :none : round(Int, 1000 * max_width),
    )
    return _cached_fixture(key, () -> begin
        basis = bond_aligned_homonuclear_qw_basis(
            bond_length = bond_length,
            core_spacing = core_spacing,
            xmax_parallel = xmax_parallel,
            xmax_transverse = xmax_transverse,
            bond_axis = bond_axis,
            nuclear_charge = nuclear_charge,
        )
        source = bond_aligned_diatomic_nested_fixed_source(basis; nside = nside)
        fixed_block = GaussletBases._nested_fixed_block(source)
        supplement = legacy_bond_aligned_diatomic_gaussian_supplement(
            atom,
            basis_name,
            basis.nuclei;
            lmax = supplement_lmax,
            max_width = max_width,
        )
        hybrid_ops = ordinary_cartesian_qiu_white_operators(
            fixed_block,
            supplement;
            nuclear_charges = fill(nuclear_charge, length(basis.nuclei)),
            interaction_treatment = :ggt_nearest,
        )
        return (
            basis = basis,
            source = source,
            fixed_block = fixed_block,
            fixed_rep = basis_representation(fixed_block),
            supplement = supplement,
            hybrid_ops = hybrid_ops,
            hybrid_rep = basis_representation(hybrid_ops),
        )
    end)
end

function _bond_aligned_heteronuclear_hybrid_qw_fixture(; bond_length::Float64 = 1.45)
    key = Symbol(:bond_aligned_heteronuclear_hybrid_qw_fixture, round(Int, 1000 * bond_length))
    return _cached_fixture(key, () -> begin
        basis = bond_aligned_heteronuclear_qw_basis(
            atoms = ("He", "H"),
            bond_length = bond_length,
            core_spacings = (0.25, 0.5),
            nuclear_charges = (2.0, 1.0),
            xmax_parallel = 6.0,
            xmax_transverse = 4.0,
            bond_axis = :z,
        )
        parent_ops = ordinary_cartesian_qiu_white_operators(
            basis;
            nuclear_charges = [2.0, 1.0],
            interaction_treatment = :ggt_nearest,
        )
        supplement = legacy_bond_aligned_heteronuclear_gaussian_supplement(
            "He",
            "cc-pVTZ",
            "H",
            "cc-pVTZ",
            basis.nuclei;
            lmax = 1,
        )
        ordinary_ops = ordinary_cartesian_qiu_white_operators(
            basis,
            supplement;
            nuclear_charges = [2.0, 1.0],
            interaction_treatment = :ggt_nearest,
        )
        (
            basis,
            parent_ops,
            GaussletBases.ordinary_cartesian_1s2_check(parent_ops),
            supplement,
            ordinary_ops,
            GaussletBases.ordinary_cartesian_1s2_check(ordinary_ops),
        )
    end)
end

function _bond_aligned_heteronuclear_nested_fixed_block_fixture(; bond_length::Float64 = 1.45)
    key = Symbol(:bond_aligned_heteronuclear_nested_fixed_block_fixture, round(Int, 1000 * bond_length))
    return _cached_fixture(key, () -> begin
        (
            basis,
            parent_ops,
            parent_check,
            supplement,
            ordinary_ops,
            ordinary_check,
        ) = _bond_aligned_heteronuclear_hybrid_qw_fixture(; bond_length = bond_length)
        expansion = coulomb_gaussian_expansion(doacc = false)
        nested = bond_aligned_diatomic_nested_fixed_block(
            basis;
            expansion = expansion,
        )
        source = nested.source
        fixed_block = nested.fixed_block
        parent_modes = eigen(Hermitian(parent_ops.one_body_hamiltonian), Hermitian(parent_ops.overlap))
        parent_ground = parent_modes.vectors[:, 1]
        projected = _nested_fixed_projected_orbital(parent_ops.overlap, fixed_block, parent_ground)
        projected_vee = _nested_vee_from_orbital(
            GaussletBases._qwrg_fixed_block_interaction_matrix(fixed_block, expansion),
            projected,
        )
        capture, projected_energy = _nested_projector_stats(
            parent_ops.overlap,
            parent_ops.one_body_hamiltonian,
            fixed_block,
            parent_ground,
        )
        (
            basis,
            parent_ops,
            parent_check,
            supplement,
            ordinary_ops,
            ordinary_check,
            expansion,
            source,
            fixed_block,
            parent_modes,
            parent_ground,
            projected,
            projected_vee,
            capture,
            projected_energy,
        )
    end)
end

function _bond_aligned_heteronuclear_nested_hybrid_qw_fixture(; bond_length::Float64 = 1.45)
    key = Symbol(:bond_aligned_heteronuclear_nested_hybrid_qw_fixture, round(Int, 1000 * bond_length))
    return _cached_fixture(key, () -> begin
        (
            basis,
            parent_ops,
            parent_check,
            supplement,
            _ordinary_ops,
            _ordinary_check,
            expansion,
            source,
            fixed_block,
            _parent_modes,
            _parent_ground,
            _projected,
            _projected_vee,
            _capture,
            _projected_energy,
        ) = _bond_aligned_heteronuclear_nested_fixed_block_fixture(; bond_length = bond_length)
        nested_ops = ordinary_cartesian_qiu_white_operators(
            fixed_block,
            supplement;
            nuclear_charges = [2.0, 1.0],
            expansion = expansion,
            interaction_treatment = :ggt_nearest,
        )
        (
            basis,
            parent_ops,
            parent_check,
            source,
            fixed_block,
            supplement,
            nested_ops,
            GaussletBases.ordinary_cartesian_1s2_check(nested_ops),
        )
    end)
end

function _bond_aligned_diatomic_nested_hybrid_qw_fixture(; bond_length::Float64 = 1.4)
    key = Symbol(:bond_aligned_diatomic_nested_hybrid_qw_fixture, round(Int, 1000 * bond_length))
    return _cached_fixture(key, () -> begin
        (
            basis,
            parent_ops,
            parent_check,
            expansion,
            source,
            fixed_block,
            _parent_modes,
            _parent_ground,
            _projected,
            _projected_vee,
            _capture,
            _projected_energy,
        ) = _bond_aligned_diatomic_nested_fixed_block_fixture(; bond_length = bond_length)
        supplement = legacy_bond_aligned_diatomic_gaussian_supplement(
            "H",
            "cc-pVTZ",
            basis.nuclei;
            lmax = 1,
        )
        nested_ops = ordinary_cartesian_qiu_white_operators(
            fixed_block,
            supplement;
            nuclear_charges = [1.0, 1.0],
            expansion = expansion,
            interaction_treatment = :ggt_nearest,
        )
        (
            basis,
            parent_ops,
            parent_check,
            source,
            fixed_block,
            supplement,
            nested_ops,
            GaussletBases.ordinary_cartesian_1s2_check(nested_ops),
        )
    end)
end

function _bond_aligned_diatomic_nested_hybrid_qw_shared_shell_experiment_fixture(
    ;
    bond_length::Float64 = 1.4,
    shared_shell_retain_xy::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_xz::Union{Nothing,Tuple{Int,Int}} = nothing,
    shared_shell_retain_yz::Union{Nothing,Tuple{Int,Int}} = nothing,
)
    retain_xy_key = isnothing(shared_shell_retain_xy) ? "default" : string(shared_shell_retain_xy[1], "_", shared_shell_retain_xy[2])
    retain_xz_key = isnothing(shared_shell_retain_xz) ? "default" : string(shared_shell_retain_xz[1], "_", shared_shell_retain_xz[2])
    retain_yz_key = isnothing(shared_shell_retain_yz) ? "default" : string(shared_shell_retain_yz[1], "_", shared_shell_retain_yz[2])
    key = Symbol(
        :bond_aligned_diatomic_nested_hybrid_qw_shared_shell_experiment_fixture,
        round(Int, 1000 * bond_length),
        :_,
        retain_xy_key,
        :_,
        retain_xz_key,
        :_,
        retain_yz_key,
    )
    return _cached_fixture(key, () -> begin
        basis, parent_ops, parent_check = _bond_aligned_diatomic_qw_fixture(; bond_length = bond_length)
        expansion = coulomb_gaussian_expansion(doacc = false)
        nested = bond_aligned_diatomic_nested_fixed_block(
            basis;
            expansion = expansion,
            shared_shell_retain_xy = shared_shell_retain_xy,
            shared_shell_retain_xz = shared_shell_retain_xz,
            shared_shell_retain_yz = shared_shell_retain_yz,
        )
        source = nested.source
        fixed_block = nested.fixed_block
        parent_modes = eigen(Hermitian(parent_ops.one_body_hamiltonian), Hermitian(parent_ops.overlap))
        parent_ground = parent_modes.vectors[:, 1]
        projected = _nested_fixed_projected_orbital(parent_ops.overlap, fixed_block, parent_ground)
        projected_vee = _nested_vee_from_orbital(
            GaussletBases._qwrg_fixed_block_interaction_matrix(fixed_block, expansion),
            projected,
        )
        capture, projected_energy = _nested_projector_stats(
            parent_ops.overlap,
            parent_ops.one_body_hamiltonian,
            fixed_block,
            parent_ground,
        )
        supplement = legacy_bond_aligned_diatomic_gaussian_supplement(
            "H",
            "cc-pVTZ",
            basis.nuclei;
            lmax = 1,
        )
        nested_ops = ordinary_cartesian_qiu_white_operators(
            fixed_block,
            supplement;
            nuclear_charges = [1.0, 1.0],
            expansion = expansion,
            interaction_treatment = :ggt_nearest,
        )
        (
            basis,
            parent_ops,
            parent_check,
            expansion,
            source,
            fixed_block,
            parent_modes,
            parent_ground,
            projected,
            projected_vee,
            capture,
            projected_energy,
            supplement,
            nested_ops,
            GaussletBases.ordinary_cartesian_1s2_check(nested_ops),
        )
    end)
end

function _nested_complete_shell_intervals(count::Int)
    count >= 15 || throw(ArgumentError("complete-shell fixture expects count >= 15"))
    outer_start = div(count - 13, 2) + 1
    outer_stop = outer_start + 12
    interval1 = (outer_start + 1):(outer_stop - 1)
    interval2 = (outer_start + 2):(outer_stop - 2)
    interval3 = (outer_start + 3):(outer_stop - 3)
    interval4 = (outer_start + 4):(outer_stop - 4)
    return interval1, interval2, interval3, interval4, interval4
end

function _nested_qiu_white_complete_shell_sequence_fixture(; basis_name::String = "cc-pVTZ", count::Int = 17, a::Float64 = 0.25, xmax::Float64 = 10.0, tail_spacing::Float64 = 10.0)
    key = Symbol(
        :nested_qiu_white_complete_shell_sequence_fixture,
        Symbol(lowercase(basis_name)),
        count,
        round(Int, 1000 * a),
        round(Int, 1000 * xmax),
    )
    return _cached_fixture(key, () -> begin
        endpoint = (count - 1) / 2
        s = asinh(xmax / a) / (endpoint - xmax / tail_spacing)
        source_basis = build_basis(MappedUniformBasisSpec(:G10;
            count = count,
            mapping = AsinhMapping(a = a, s = s, tail_spacing = tail_spacing),
            reference_spacing = 1.0,
        ))
        expansion = coulomb_gaussian_expansion(doacc = false)
        term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
        bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
            source_basis;
            exponents = expansion.exponents,
            backend = :numerical_reference,
            refinement_levels = 0,
        )

        interval1, interval2, interval3, interval4, core5 = _nested_complete_shell_intervals(count)

        shell1_face = GaussletBases._nested_rectangular_shell(
            bundle,
            interval1,
            interval1,
            interval1;
            retain_xy = (4, 3),
            retain_xz = (4, 3),
            retain_yz = (4, 3),
            x_fixed = (first(interval1) - 1, last(interval1) + 1),
            y_fixed = (first(interval1) - 1, last(interval1) + 1),
            z_fixed = (first(interval1) - 1, last(interval1) + 1),
            term_coefficients = term_coefficients,
        )
        shell2_face = GaussletBases._nested_rectangular_shell(
            bundle,
            interval2,
            interval2,
            interval2;
            retain_xy = (4, 3),
            retain_xz = (4, 3),
            retain_yz = (4, 3),
            x_fixed = (first(interval2) - 1, last(interval2) + 1),
            y_fixed = (first(interval2) - 1, last(interval2) + 1),
            z_fixed = (first(interval2) - 1, last(interval2) + 1),
            term_coefficients = term_coefficients,
        )
        shell3_face = GaussletBases._nested_rectangular_shell(
            bundle,
            interval3,
            interval3,
            interval3;
            retain_xy = (4, 3),
            retain_xz = (4, 3),
            retain_yz = (4, 3),
            x_fixed = (first(interval3) - 1, last(interval3) + 1),
            y_fixed = (first(interval3) - 1, last(interval3) + 1),
            z_fixed = (first(interval3) - 1, last(interval3) + 1),
            term_coefficients = term_coefficients,
        )
        shell4_face = GaussletBases._nested_rectangular_shell(
            bundle,
            interval4,
            interval4,
            interval4;
            retain_xy = (4, 3),
            retain_xz = (4, 3),
            retain_yz = (4, 3),
            x_fixed = (first(interval4) - 1, last(interval4) + 1),
            y_fixed = (first(interval4) - 1, last(interval4) + 1),
            z_fixed = (first(interval4) - 1, last(interval4) + 1),
            term_coefficients = term_coefficients,
        )

        shell1_complete = GaussletBases._nested_complete_rectangular_shell(
            bundle,
            interval1,
            interval1,
            interval1;
            retain_xy = (4, 3),
            retain_xz = (4, 3),
            retain_yz = (4, 3),
            retain_x_edge = 3,
            retain_y_edge = 3,
            retain_z_edge = 3,
            x_fixed = (first(interval1) - 1, last(interval1) + 1),
            y_fixed = (first(interval1) - 1, last(interval1) + 1),
            z_fixed = (first(interval1) - 1, last(interval1) + 1),
            term_coefficients = term_coefficients,
        )
        shell2_complete = GaussletBases._nested_complete_rectangular_shell(
            bundle,
            interval2,
            interval2,
            interval2;
            retain_xy = (4, 3),
            retain_xz = (4, 3),
            retain_yz = (4, 3),
            retain_x_edge = 3,
            retain_y_edge = 3,
            retain_z_edge = 3,
            x_fixed = (first(interval2) - 1, last(interval2) + 1),
            y_fixed = (first(interval2) - 1, last(interval2) + 1),
            z_fixed = (first(interval2) - 1, last(interval2) + 1),
            term_coefficients = term_coefficients,
        )
        shell3_complete = GaussletBases._nested_complete_rectangular_shell(
            bundle,
            interval3,
            interval3,
            interval3;
            retain_xy = (4, 3),
            retain_xz = (4, 3),
            retain_yz = (4, 3),
            retain_x_edge = 3,
            retain_y_edge = 3,
            retain_z_edge = 3,
            x_fixed = (first(interval3) - 1, last(interval3) + 1),
            y_fixed = (first(interval3) - 1, last(interval3) + 1),
            z_fixed = (first(interval3) - 1, last(interval3) + 1),
            term_coefficients = term_coefficients,
        )
        shell4_complete = GaussletBases._nested_complete_rectangular_shell(
            bundle,
            interval4,
            interval4,
            interval4;
            retain_xy = (4, 3),
            retain_xz = (4, 3),
            retain_yz = (4, 3),
            retain_x_edge = 3,
            retain_y_edge = 3,
            retain_z_edge = 3,
            x_fixed = (first(interval4) - 1, last(interval4) + 1),
            y_fixed = (first(interval4) - 1, last(interval4) + 1),
            z_fixed = (first(interval4) - 1, last(interval4) + 1),
            term_coefficients = term_coefficients,
        )

        shell_plus_core = GaussletBases._nested_shell_plus_core(
            bundle,
            shell1_face,
            interval1,
            interval1,
            interval1,
            term_coefficients = term_coefficients,
        )
        face_sequence = GaussletBases._nested_shell_sequence(
            bundle,
            core5,
            core5,
            core5,
            [shell1_face, shell2_face, shell3_face, shell4_face];
            enforce_coverage = false,
            term_coefficients = term_coefficients,
        )
        complete_sequence = GaussletBases._nested_shell_sequence(
            bundle,
            core5,
            core5,
            core5,
            [shell1_complete, shell2_complete, shell3_complete, shell4_complete],
            term_coefficients = term_coefficients,
        )

        fixed_shell_plus_core = GaussletBases._nested_fixed_block(shell_plus_core, bundle)
        fixed_face_sequence = GaussletBases._nested_fixed_block(face_sequence, bundle)
        fixed_complete_sequence = GaussletBases._nested_fixed_block(complete_sequence, bundle)

        legacy = legacy_atomic_gaussian_supplement("He", basis_name; lmax = 0)
        baseline = ordinary_cartesian_qiu_white_operators(
            source_basis,
            legacy;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        shell_plus_core_ops = ordinary_cartesian_qiu_white_operators(
            fixed_shell_plus_core,
            legacy;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        face_sequence_ops = ordinary_cartesian_qiu_white_operators(
            fixed_face_sequence,
            legacy;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        complete_sequence_ops = ordinary_cartesian_qiu_white_operators(
            fixed_complete_sequence,
            legacy;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )

        (
            source_basis,
            bundle,
            shell1_face,
            shell2_face,
            shell3_face,
            shell4_face,
            shell1_complete,
            shell2_complete,
            shell3_complete,
            shell4_complete,
            interval1,
            interval2,
            interval3,
            interval4,
            core5,
            shell_plus_core,
            face_sequence,
            complete_sequence,
            fixed_shell_plus_core,
            fixed_face_sequence,
            fixed_complete_sequence,
            legacy,
            baseline,
            shell_plus_core_ops,
            face_sequence_ops,
            complete_sequence_ops,
            GaussletBases.ordinary_cartesian_1s2_check(baseline),
            GaussletBases.ordinary_cartesian_1s2_check(shell_plus_core_ops),
            GaussletBases.ordinary_cartesian_1s2_check(face_sequence_ops),
            GaussletBases.ordinary_cartesian_1s2_check(complete_sequence_ops),
        )
    end)
end

function _nested_qiu_white_hierarchical_core_fixture(; basis_name::String = "cc-pVTZ", count::Int = 17, a::Float64 = 0.25, xmax::Float64 = 10.0, tail_spacing::Float64 = 10.0)
    key = Symbol(
        :nested_qiu_white_hierarchical_core_fixture,
        Symbol(lowercase(basis_name)),
        count,
        round(Int, 1000 * a),
        round(Int, 1000 * xmax),
    )
    return _cached_fixture(key, () -> begin
        (
            source_basis,
            bundle,
            _shell1_face,
            _shell2_face,
            _shell3_face,
            _shell4_face,
            shell1_complete,
            shell2_complete,
            shell3_complete,
            shell4_complete,
            _interval1,
            _interval2,
            _interval3,
            _interval4,
            core5,
            _shell_plus_core,
            _face_sequence,
            complete_sequence,
            _fixed_shell_plus_core,
            _fixed_face_sequence,
            fixed_complete_sequence,
            legacy,
            baseline,
            _shell_plus_core_ops,
            _face_sequence_ops,
            complete_sequence_ops,
            baseline_check,
            _shell_plus_core_check,
            _face_sequence_check,
            complete_sequence_check,
        ) = _nested_qiu_white_complete_shell_sequence_fixture(
            basis_name = basis_name,
            count = count,
            a = a,
            xmax = xmax,
            tail_spacing = tail_spacing,
        )

        expansion = coulomb_gaussian_expansion(doacc = false)
        term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
        refined_data = GaussletBases._nested_shell_sequence_with_hierarchical_core_refinement(
            bundle,
            core5,
            core5,
            core5,
            [shell1_complete, shell2_complete, shell3_complete, shell4_complete];
            retain_xy = (2, 2),
            retain_xz = (2, 2),
            retain_yz = (2, 2),
            retain_x_edge = 2,
            retain_y_edge = 2,
            retain_z_edge = 2,
            term_coefficients = term_coefficients,
        )
        fixed_refined_sequence = GaussletBases._nested_fixed_block(refined_data.sequence, bundle)
        refined_sequence_ops = ordinary_cartesian_qiu_white_operators(
            fixed_refined_sequence,
            legacy;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )

        (
            source_basis,
            bundle,
            core5,
            complete_sequence,
            fixed_complete_sequence,
            complete_sequence_ops,
            complete_sequence_check,
            refined_data.refined_core,
            refined_data.sequence,
            fixed_refined_sequence,
            refined_sequence_ops,
            GaussletBases.ordinary_cartesian_1s2_check(refined_sequence_ops),
            legacy,
            baseline,
            baseline_check,
        )
    end)
end

function _nested_parent_fixed_problem(bundle, expansion; Z::Float64 = 2.0)
    gg = GaussletBases._qwrg_gausslet_1d_blocks(bundle)
    n1d = size(gg.overlap_gg, 1)
    overlap = zeros(Float64, n1d^3, n1d^3)
    GaussletBases._qwrg_fill_product_matrix!(overlap, gg.overlap_gg, gg.overlap_gg, gg.overlap_gg)
    one_body = GaussletBases._qwrg_gausslet_one_body_matrix(gg, expansion; Z = Z)
    interaction = GaussletBases._qwrg_gausslet_interaction_matrix(gg, expansion)
    return overlap, one_body, interaction
end

function _nested_fixed_projected_orbital(
    overlap_parent::AbstractMatrix{<:Real},
    fixed_block,
    parent_vector::AbstractVector{<:Real},
)
    rhs = transpose(fixed_block.coefficient_matrix) * (overlap_parent * parent_vector)
    coeffs = fixed_block.overlap \ rhs
    norm2 = real(dot(coeffs, fixed_block.overlap * coeffs))
    return coeffs ./ sqrt(norm2)
end

function _nested_vee_from_orbital(
    interaction::AbstractMatrix{<:Real},
    orbital::AbstractVector{<:Real},
)
    weights = Float64[abs2(coefficient) for coefficient in orbital]
    weights ./= sum(weights)
    return Float64(real(dot(weights, interaction * weights)))
end

function _nested_projector_stats(
    overlap_parent::AbstractMatrix{<:Real},
    one_body_parent::AbstractMatrix{<:Real},
    fixed_block,
    parent_vector::AbstractVector{<:Real},
)
    rhs = transpose(fixed_block.coefficient_matrix) * (overlap_parent * parent_vector)
    projected_coords = fixed_block.overlap \ rhs
    retained =
        Float64(real(dot(rhs, projected_coords) / dot(parent_vector, overlap_parent * parent_vector)))
    projected = fixed_block.coefficient_matrix * projected_coords
    projected_norm = real(dot(projected, overlap_parent * projected))
    projected_energy = Float64(real(dot(projected, one_body_parent * projected) / projected_norm))
    return retained, projected_energy
end

function _quick_ordinary_sho_smoke_fixture()
    return _cached_fixture(:quick_ordinary_sho_smoke_fixture, () -> begin
        basis = build_basis(MappedUniformBasisSpec(:G10;
            count = 3,
            mapping = fit_asinh_mapping_for_strength(s = 0.2, npoints = 3, xmax = 4.0),
            reference_spacing = 1.0,
        ))
        core_gaussians = [
            Gaussian(center = 0.0, width = 0.2),
            Gaussian(center = 0.0, width = 0.6),
        ]
        hybrid_reference = hybrid_mapped_ordinary_basis(
            basis;
            core_gaussians = core_gaussians,
            backend = :numerical_reference,
        )
        hybrid_analytic = hybrid_mapped_ordinary_basis(
            basis;
            core_gaussians = core_gaussians,
            backend = :pgdg_localized_experimental,
        )
        hybrid_analytic_core = GaussletBases._ordinary_sho_core(
            mapped_ordinary_one_body_operators(hybrid_analytic),
        )
        centered_analytic = GaussletBases._ordinary_sho_spectrum_from_core(
            hybrid_analytic_core;
            omega = 0.25,
            center = 0.0,
            nev = 3,
        )
        centered_hamiltonian = ordinary_sho_hamiltonian(
            hybrid_analytic;
            omega = 0.25,
            center = 0.0,
        )
        (
            hybrid_analytic,
            centered_analytic,
            centered_hamiltonian,
        )
    end)
end

function _slow_ordinary_sho_fixture()
    return _cached_fixture(:slow_ordinary_sho_fixture, () -> begin
        basis = build_basis(MappedUniformBasisSpec(:G10;
            count = 5,
            mapping = fit_asinh_mapping_for_strength(s = 0.2, npoints = 5, xmax = 6.0),
            reference_spacing = 1.0,
        ))
        core_gaussians = [
            Gaussian(center = 0.0, width = 0.2),
            Gaussian(center = 0.0, width = 0.6),
        ]
        hybrid_reference = hybrid_mapped_ordinary_basis(
            basis;
            core_gaussians = core_gaussians,
            backend = :numerical_reference,
        )
        hybrid_analytic = hybrid_mapped_ordinary_basis(
            basis;
            core_gaussians = core_gaussians,
            backend = :pgdg_localized_experimental,
        )
        hybrid_reference_core = GaussletBases._ordinary_sho_core(
            mapped_ordinary_one_body_operators(hybrid_reference),
        )
        hybrid_analytic_core = GaussletBases._ordinary_sho_core(
            mapped_ordinary_one_body_operators(hybrid_analytic),
        )
        centered_reference = GaussletBases._ordinary_sho_spectrum_from_core(
            hybrid_reference_core;
            omega = 0.25,
            center = 0.0,
            nev = 3,
        )
        centered_analytic = GaussletBases._ordinary_sho_spectrum_from_core(
            hybrid_analytic_core;
            omega = 0.25,
            center = 0.0,
            nev = 3,
        )
        shifted_reference = GaussletBases._ordinary_sho_spectrum_from_core(
            hybrid_reference_core;
            omega = 0.25,
            center = 1.0,
            nev = 3,
        )
        shifted_analytic = GaussletBases._ordinary_sho_spectrum_from_core(
            hybrid_analytic_core;
            omega = 0.25,
            center = 1.0,
            nev = 3,
        )
        stress_basis = build_basis(MappedUniformBasisSpec(:G10;
            count = 5,
            mapping = fit_asinh_mapping_for_strength(s = 2.0, npoints = 5, xmax = 6.0),
            reference_spacing = 1.0,
        ))
        stress_reference_core = GaussletBases._ordinary_sho_core(
            mapped_ordinary_one_body_operators(stress_basis; backend = :numerical_reference),
        )
        stress_analytic_core = GaussletBases._ordinary_sho_core(
            mapped_ordinary_one_body_operators(stress_basis; backend = :pgdg_localized_experimental),
        )
        stress_reference = GaussletBases._ordinary_sho_spectrum_from_core(
            stress_reference_core;
            omega = 0.5,
            center = 0.0,
            nev = 3,
        )
        stress_analytic = GaussletBases._ordinary_sho_spectrum_from_core(
            stress_analytic_core;
            omega = 0.5,
            center = 0.0,
            nev = 3,
        )
        (
            centered_reference,
            centered_analytic,
            shifted_reference,
            shifted_analytic,
            stress_reference,
            stress_analytic,
        )
    end)
end

function _cartesian_hydrogen_energy(
    overlap_1d::AbstractMatrix,
    kinetic_1d::AbstractMatrix,
    gaussian_factors::AbstractVector{<:AbstractMatrix},
    expansion::CoulombGaussianExpansion;
    Z::Real = 1.0,
)
    transformed = _orthonormalize_1d(overlap_1d, [kinetic_1d, gaussian_factors...])
    kinetic_orth = first(transformed)
    gaussian_orth = collect(transformed[2:end])
    n1d = size(overlap_1d, 1)
    identity_1d = Matrix{Float64}(I, n1d, n1d)
    hamiltonian =
        kron(kinetic_orth, kron(identity_1d, identity_1d)) +
        kron(identity_1d, kron(kinetic_orth, identity_1d)) +
        kron(identity_1d, kron(identity_1d, kinetic_orth))

    for term in eachindex(expansion.coefficients)
        factor = gaussian_orth[term]
        hamiltonian .-= Float64(Z) * expansion.coefficients[term] .* kron(factor, kron(factor, factor))
    end

    energy = minimum(eigen(Hermitian(hamiltonian)).values)
    return hamiltonian, energy
end

function _localized_numerical_reference_1d(
    basis::MappedUniformBasis,
    exponents::AbstractVector{<:Real};
    center::Real = 0.0,
)
    representation = basis_representation(basis; operators = (:overlap, :position, :kinetic))
    transform, centers_value = _cleanup_comx_transform(
        representation.basis_matrices.overlap,
        representation.basis_matrices.position,
        integral_weights(basis),
    )
    overlap = Matrix{Float64}(transpose(transform) * representation.basis_matrices.overlap * transform)
    kinetic = Matrix{Float64}(transpose(transform) * representation.basis_matrices.kinetic * transform)
    gaussian_factors = Matrix{Float64}[
        Matrix{Float64}(transpose(transform) * gaussian_factor_matrix(basis; exponent = exponent, center = center) * transform)
        for exponent in exponents
    ]
    return (transform, centers_value, overlap, kinetic, gaussian_factors)
end

function _mapped_pgdg_1d_fixture()
    return _cached_fixture(:mapped_pgdg_1d_fixture, () -> begin
        mapping = fit_asinh_mapping_for_extent(npoints = 5, xmax = 6.0)
        basis = build_basis(MappedUniformBasisSpec(:G10;
            count = 5,
            mapping = mapping,
            reference_spacing = 1.0,
        ))
        representation = basis_representation(basis; operators = (:overlap, :kinetic))
        prototype = mapped_pgdg_prototype(basis)
        localized = mapped_pgdg_localized(prototype)
        refined = GaussletBases.mapped_pgdg_logfit_prototype(basis)
        refined_localized = mapped_pgdg_localized(refined)
        exponent = 0.35
        factor_numeric = gaussian_factor_matrix(basis; exponent = exponent, center = 0.0)
        factor_pgdg = gaussian_factor_matrix(prototype; exponent = exponent, center = 0.0)
        factor_localized = gaussian_factor_matrix(localized; exponent = exponent, center = 0.0)
        factor_refined = gaussian_factor_matrix(refined; exponent = exponent, center = 0.0)
        factor_refined_localized = gaussian_factor_matrix(refined_localized; exponent = exponent, center = 0.0)
        (
            basis,
            mapping,
            representation,
            prototype,
            localized,
            refined,
            refined_localized,
            exponent,
            factor_numeric,
            factor_pgdg,
            factor_localized,
            factor_refined,
            factor_refined_localized,
        )
    end)
end

function _mapped_cartesian_hydrogen_comparison(expansion::CoulombGaussianExpansion)
    basis = build_basis(MappedUniformBasisSpec(:G10;
        count = 5,
        mapping = fit_asinh_mapping_for_extent(npoints = 5, xmax = 6.0),
        reference_spacing = 1.0,
    ))
    representation = basis_representation(basis; operators = (:overlap, :kinetic))
    overlap_numeric = representation.basis_matrices.overlap
    kinetic_numeric = representation.basis_matrices.kinetic
    factors_numeric = gaussian_factor_matrices(
        basis;
        exponents = expansion.exponents,
        center = 0.0,
    )
    hamiltonian_numeric, energy_numeric = _cartesian_hydrogen_energy(
        overlap_numeric,
        kinetic_numeric,
        factors_numeric,
        expansion;
        Z = 1.0,
    )

    prototype = mapped_pgdg_prototype(basis)
    overlap_pgdg = overlap_matrix(prototype)
    position_pgdg = position_matrix(prototype)
    kinetic_pgdg = kinetic_matrix(prototype)
    factors_pgdg = gaussian_factor_matrices(
        prototype;
        exponents = expansion.exponents,
        center = 0.0,
    )
    hamiltonian_pgdg, energy_pgdg = _cartesian_hydrogen_energy(
        overlap_pgdg,
        kinetic_pgdg,
        factors_pgdg,
        expansion;
        Z = 1.0,
    )

    localized = mapped_pgdg_localized(prototype)
    overlap_localized = overlap_matrix(localized)
    position_localized = position_matrix(localized)
    kinetic_localized = kinetic_matrix(localized)
    factors_localized = gaussian_factor_matrices(
        localized;
        exponents = expansion.exponents,
        center = 0.0,
    )
    hamiltonian_localized, energy_localized = _cartesian_hydrogen_energy(
        overlap_localized,
        kinetic_localized,
        factors_localized,
        expansion;
        Z = 1.0,
    )

    refined = GaussletBases.mapped_pgdg_logfit_prototype(basis)
    overlap_refined = overlap_matrix(refined)
    position_refined = position_matrix(refined)
    kinetic_refined = kinetic_matrix(refined)
    factors_refined = gaussian_factor_matrices(
        refined;
        exponents = expansion.exponents,
        center = 0.0,
    )
    hamiltonian_refined, energy_refined = _cartesian_hydrogen_energy(
        overlap_refined,
        kinetic_refined,
        factors_refined,
        expansion;
        Z = 1.0,
    )

    refined_localized = mapped_pgdg_localized(refined)
    overlap_refined_localized = overlap_matrix(refined_localized)
    position_refined_localized = position_matrix(refined_localized)
    kinetic_refined_localized = kinetic_matrix(refined_localized)
    factors_refined_localized = gaussian_factor_matrices(
        refined_localized;
        exponents = expansion.exponents,
        center = 0.0,
    )
    hamiltonian_refined_localized, energy_refined_localized = _cartesian_hydrogen_energy(
        overlap_refined_localized,
        kinetic_refined_localized,
        factors_refined_localized,
        expansion;
        Z = 1.0,
    )

    return (
        basis,
        prototype,
        localized,
        refined,
        refined_localized,
        overlap_numeric,
        kinetic_numeric,
        hamiltonian_numeric,
        energy_numeric,
        overlap_pgdg,
        position_pgdg,
        kinetic_pgdg,
        hamiltonian_pgdg,
        energy_pgdg,
        overlap_localized,
        position_localized,
        kinetic_localized,
        hamiltonian_localized,
        energy_localized,
        overlap_refined,
        position_refined,
        kinetic_refined,
        hamiltonian_refined,
        energy_refined,
        overlap_refined_localized,
        position_refined_localized,
        kinetic_refined_localized,
        hamiltonian_refined_localized,
        energy_refined_localized,
    )
end

function _quick_radial_atomic_fixture()
    return _cached_fixture(:quick_radial_atomic_fixture, () -> begin
        rb, grid = _quick_radial_operator_fixture()
        radial_ops = atomic_operators(rb, grid; Z = 2.0, lmax = 2)
        channels = ylm_channels(2)
        atom = atomic_one_body_operators(radial_ops, channels)
        ida = atomic_ida_operators(radial_ops, channels)
        (rb, grid, radial_ops, channels, atom, ida)
    end)
end

function _quick_hydrogen_ylm_fixture()
    return _cached_fixture(:quick_hydrogen_ylm_fixture, () -> begin
        Z = 1.0
        s = 0.2
        lmax = 2
        rb = build_basis(RadialBasisSpec(:G10;
            rmax = 30.0,
            mapping = AsinhMapping(c = s / (2Z), s = s),
            reference_spacing = 1.0,
            tails = 6,
            odd_even_kmax = 6,
            xgaussians = XGaussian[],
        ))
        grid = radial_quadrature(rb)
        radial_ops = atomic_operators(rb, grid; Z = Z, lmax = lmax)
        atom = atomic_one_body_operators(radial_ops; lmax = lmax)
        (rb, grid, radial_ops, atom)
    end)
end

function _shell_local_injected_angular_fixture(order::Int)
    key = Symbol("shell_local_injected_angular_fixture_", order)
    return _cached_fixture(key, () -> begin
        build_shell_local_injected_angular_basis(order)
    end)
end

function _shell_local_angular_profile_fixture(order::Int)
    key = Symbol("shell_local_angular_profile_fixture_", order)
    return _cached_fixture(key, () -> begin
        shell_local_angular_profile(order)
    end)
end

function _shell_local_angular_profile_overlap_fixture(source_order::Int, target_order::Int)
    key = Symbol("shell_local_angular_profile_overlap_fixture_$(source_order)_$(target_order)")
    return _cached_fixture(key, () -> begin
        adjacent_shell_local_angular_profile_overlap(source_order, target_order)
    end)
end

function _atomic_shell_local_angular_fixture()
    return _cached_fixture(:atomic_shell_local_angular_fixture, () -> begin
        shell_radii = [0.3, 0.8, 1.5, 3.5]
        shell_orders = [15, 32, 51, 32]
        build_atomic_shell_local_angular_assembly(shell_radii; shell_orders = shell_orders)
    end)
end

function _atomic_injected_angular_one_body_benchmark_fixture()
    return _cached_fixture(:atomic_injected_angular_one_body_benchmark_fixture, () -> begin
        _, _, radial_ops, _, _, _ = _quick_radial_atomic_fixture()
        build_atomic_injected_angular_one_body_benchmark(radial_ops)
    end)
end

function _atomic_injected_angular_hf_style_benchmark_fixture()
    return _cached_fixture(:atomic_injected_angular_hf_style_benchmark_fixture, () -> begin
        _, _, radial_ops, _, _, _ = _quick_radial_atomic_fixture()
        build_atomic_injected_angular_hf_style_benchmark(radial_ops)
    end)
end

function _atomic_injected_angular_small_ed_benchmark_fixture()
    return _cached_fixture(:atomic_injected_angular_small_ed_benchmark_fixture, () -> begin
        _, _, radial_ops, _, _, _ = _quick_radial_atomic_fixture()
        build_atomic_injected_angular_small_ed_benchmark(radial_ops)
    end)
end

function _atomic_injected_angular_hfdmrg_hf_adapter_fixture()
    return _cached_fixture(:atomic_injected_angular_hfdmrg_hf_adapter_fixture, () -> begin
        _, _, radial_ops, _, _, _ = _quick_radial_atomic_fixture()
        build_atomic_injected_angular_hfdmrg_hf_adapter(radial_ops)
    end)
end

function _paper_style_angular_anchor_radial_fixture(; Z::Float64 = 2.0, lmax::Int = 2)
    ztag = replace(string(round(Z; digits = 2)), "." => "p")
    key = Symbol("paper_style_angular_anchor_radial_fixture_Z", ztag, "_lmax", lmax)
    return _cached_fixture(key, () -> begin
        s = 0.2
        rb = build_basis(RadialBasisSpec(:G10;
            rmax = 30.0,
            mapping = AsinhMapping(c = s / (2Z), s = s),
            reference_spacing = 1.0,
            tails = 6,
            odd_even_kmax = 6,
            xgaussian_count = 2,
            rmax_count_policy = :legacy_strict_trim,
        ))
        grid = radial_quadrature(rb)
        radial_ops = atomic_operators(rb, grid; Z = Z, lmax = lmax)
        (rb, grid, radial_ops)
    end)
end

function _paper_style_angular_hfdmrg_payload_fixture(order::Int; Z::Float64 = 2.0, lmax::Int = 2)
    ztag = replace(string(round(Z; digits = 2)), "." => "p")
    key = Symbol("paper_style_angular_hfdmrg_payload_fixture_Z", ztag, "_order", order, "_lmax", lmax)
    return _cached_fixture(key, () -> begin
        _, _, radial_ops = _paper_style_angular_anchor_radial_fixture(; Z = Z, lmax = lmax)
        build_atomic_injected_angular_hfdmrg_payload(
            radial_ops;
            shell_orders = fill(order, length(radial_ops.shell_centers_r)),
            nelec = Int(round(Z)),
        )
    end)
end

function _paper_style_fixed_radial_angular_sequence_fixture(;
    Z::Float64 = 2.0,
    lmax::Int = 2,
    N_sph_values::AbstractVector{<:Integer} = [10, 15, 32],
)
    ztag = replace(string(round(Z; digits = 2)), "." => "p")
    ntag = join(Int.(N_sph_values), "_")
    key = Symbol(
        "paper_style_fixed_radial_angular_sequence_fixture_Z",
        ztag,
        "_Nsph_",
        ntag,
        "_lmax_",
        lmax,
    )
    return _cached_fixture(key, () -> begin
        _, _, radial_ops = _paper_style_angular_anchor_radial_fixture(; Z = Z, lmax = lmax)
        build_atomic_fixed_radial_angular_sequence(
            radial_ops,
            Int[Int(value) for value in N_sph_values];
            nelec = Int(round(Z)),
        )
    end)
end

function _paper_style_angular_one_body_benchmark_fixture(
    order::Int;
    Z::Float64 = 4.0,
    lmax::Int = 2,
)
    ztag = replace(string(round(Z; digits = 2)), "." => "p")
    key = Symbol(
        "paper_style_angular_one_body_benchmark_fixture_Z",
        ztag,
        "_order",
        order,
        "_lmax",
        lmax,
    )
    return _cached_fixture(key, () -> begin
        _, _, radial_ops = _paper_style_angular_anchor_radial_fixture(; Z = Z, lmax = lmax)
        build_atomic_injected_angular_one_body_benchmark(
            radial_ops;
            shell_orders = fill(order, length(radial_ops.shell_centers_r)),
        )
    end)
end

function _paper_style_angular_cartesian_moments_fixture(
    order::Int;
    Z::Float64 = 2.0,
    lmax::Int = 2,
    l_inject::Union{Int,Symbol} = :auto,
)
    ztag = replace(string(round(Z; digits = 2)), "." => "p")
    ltag = l_inject isa Symbol ? String(l_inject) : string(l_inject)
    key = Symbol(
        "paper_style_angular_cartesian_moments_fixture_Z",
        ztag,
        "_order",
        order,
        "_lmax",
        lmax,
        "_linject_",
        ltag,
    )
    return _cached_fixture(key, () -> begin
        rb, grid, radial_ops = _paper_style_angular_anchor_radial_fixture(; Z = Z, lmax = lmax)
        benchmark = build_atomic_injected_angular_one_body_benchmark(
            radial_ops;
            shell_orders = fill(order, length(radial_ops.shell_centers_r)),
            l_inject = l_inject,
        )
        bundle_explicit = build_atomic_injected_angular_cartesian_moments(
            benchmark;
            radial_basis = rb,
            radial_grid = grid,
        )
        bundle_reconstructed = build_atomic_injected_angular_cartesian_moments(benchmark)
        return (rb, grid, radial_ops, benchmark, bundle_explicit, bundle_reconstructed)
    end)
end

function _paper_style_angular_hf_style_benchmark_fixture(
    order::Int;
    Z::Float64 = 10.0,
    lmax::Int = 6,
)
    ztag = replace(string(round(Z; digits = 2)), "." => "p")
    key = Symbol(
        "paper_style_angular_hf_style_benchmark_fixture_Z",
        ztag,
        "_order",
        order,
        "_lmax",
        lmax,
    )
    return _cached_fixture(key, () -> begin
        _, _, radial_ops = _paper_style_angular_anchor_radial_fixture(; Z = Z, lmax = lmax)
        build_atomic_injected_angular_hf_style_benchmark(
            radial_ops;
            shell_orders = fill(order, length(radial_ops.shell_centers_r)),
            nelec = Int(round(Z)),
        )
    end)
end

function _solve_hfdmrg_from_payload_direct(
    payload::AtomicInjectedAngularHFDMRGHFAdapter,
    hfdmrg;
    nblockcenter::Int = 2,
    blocksize::Int = 100,
    maxiter::Int = 100,
    cutoff::Real = 1.0e-8,
    scf_cutoff::Real = 1.0e-9,
    verbose::Bool = false,
)
    if payload.solver_mode == :restricted_closed_shell
        psiup, psidn, energy = hfdmrg.solve_hfdmrg(
            payload.hamiltonian,
            payload.interaction,
            payload.psiup0;
            nblockcenter = nblockcenter,
            blocksize = blocksize,
            maxiter = maxiter,
            cutoff = cutoff,
            scf_cutoff = scf_cutoff,
            verbose = verbose,
        )
    else
        psiup, psidn, energy = hfdmrg.solve_hfdmrg(
            payload.hamiltonian,
            payload.interaction,
            payload.psiup0,
            payload.psidn0;
            nblockcenter = nblockcenter,
            blocksize = blocksize,
            maxiter = maxiter,
            cutoff = cutoff,
            scf_cutoff = scf_cutoff,
            verbose = verbose,
        )
    end
    return (psiup = psiup, psidn = psidn, energy = energy)
end

function _local_hfdmrg_module()
    return _cached_fixture(:local_hfdmrg_module, () -> begin
        path = "/Users/srw/Dropbox/codexhome/work/hfdmrg/src"
        isdir(path) || return nothing
        path in LOAD_PATH || push!(LOAD_PATH, path)
        try
            @eval Main using HFDMRG
            return Base.invokelatest(getfield, Main, :HFDMRG)
        catch
            return nothing
        end
    end)
end

function _tiny_atomic_ida_two_electron_fixture()
    return _cached_fixture(:tiny_atomic_ida_two_electron_fixture, () -> begin
        Z = 2.0
        s = 2.0
        lmax = 1
        rb = build_basis(RadialBasisSpec(:G10;
            rmax = 5.0,
            mapping = AsinhMapping(c = s / (2Z), s = s),
            reference_spacing = 1.0,
            tails = 6,
            odd_even_kmax = 6,
            xgaussians = XGaussian[],
        ))
        grid = radial_quadrature(rb; accuracy = :medium)
        radial_ops = atomic_operators(rb, grid; Z = Z, lmax = lmax)
        ida = atomic_ida_operators(radial_ops; lmax = lmax)
        problem = atomic_ida_two_electron_problem(ida)
        (rb, grid, radial_ops, ida, problem)
    end)
end

function _tiny_atomic_ida_lanczos_fixture()
    return _cached_fixture(:tiny_atomic_ida_lanczos_fixture, () -> begin
        Z = 2.0
        s = 0.5
        lmax = 0
        rb = build_basis(RadialBasisSpec(:G10;
            rmax = 8.0,
            mapping = AsinhMapping(c = s / (2Z), s = s),
            reference_spacing = 1.0,
            tails = 6,
            odd_even_kmax = 6,
            xgaussians = XGaussian[],
        ))
        grid = radial_quadrature(rb; accuracy = :high)
        radial_ops = atomic_operators(rb, grid; Z = Z, lmax = lmax)
        ida = atomic_ida_operators(radial_ops; lmax = lmax)
        problem = atomic_ida_two_electron_problem(ida)
        (rb, grid, radial_ops, ida, problem)
    end)
end

function _tiny_atomic_ida_uhf_fixture()
    return _cached_fixture(:tiny_atomic_ida_uhf_fixture, () -> begin
        rb, grid, radial_ops, ida, problem = _tiny_atomic_ida_lanczos_fixture()
        scf = uhf_scf(ida; nalpha = 1, nbeta = 1, maxiter = 80, damping = 0.25, tol = 1.0e-10)
        (rb, grid, radial_ops, ida, problem, scf)
    end)
end

function _direct_dense_angular_kernel(table, channels, L)
    nchannels = length(channels)
    prefactor = 4 * pi / (2 * L + 1)
    kernel = zeros(Float64, nchannels, nchannels, nchannels, nchannels)

    for alpha in 1:nchannels, alphap in 1:nchannels, beta in 1:nchannels, betap in 1:nchannels
        total = 0.0
        for M in -L:L
            total += prefactor *
                     (isodd(M) ? -1.0 : 1.0) *
                     GaussletBases.gaunt_value(
                        table,
                        L,
                        channels[alpha].l,
                        channels[alpha].m,
                        channels[alphap].l,
                        channels[alphap].m,
                        M,
                    ) *
                     GaussletBases.gaunt_value(
                        table,
                        L,
                        channels[beta].l,
                        channels[beta].m,
                        channels[betap].l,
                        channels[betap].m,
                        -M,
                    )
        end
        kernel[alpha, alphap, beta, betap] = total
    end

    return kernel
end

function _dense_direct_reference(ops::AtomicIDAOperators, density::AbstractMatrix{<:Real})
    radial_dim = size(radial_multipole(ops, 0), 1)
    nchannels = length(ops.one_body.channels)
    norbitals = length(orbitals(ops))
    size(density) == (norbitals, norbitals) || throw(DimensionMismatch("density matrix has the wrong size"))

    orbital_index(channel_index, radial_index) = (channel_index - 1) * radial_dim + radial_index
    reference = zeros(Float64, norbitals, norbitals)
    radial_multipoles = [radial_multipole(ops, L) for L in 0:GaussletBases.gaunt_Lmax(ops.gaunt_table)]
    angular_kernels = [angular_kernel(ops, L) for L in 0:GaussletBases.gaunt_Lmax(ops.gaunt_table)]

    for left_radial in 1:radial_dim, left_channel in 1:nchannels, right_channel in 1:nchannels
        total = 0.0
        for right_radial in 1:radial_dim, (level_index, multipole) in enumerate(radial_multipoles), beta in 1:nchannels, betap in 1:nchannels
            total += multipole[left_radial, right_radial] *
                     angular_kernels[level_index][left_channel, right_channel, beta, betap] *
                     density[orbital_index(beta, right_radial), orbital_index(betap, right_radial)]
        end
        reference[orbital_index(left_channel, left_radial), orbital_index(right_channel, left_radial)] = total
    end

    return 0.5 .* (reference .+ transpose(reference))
end

function _dense_exchange_reference(ops::AtomicIDAOperators, density::AbstractMatrix{<:Real})
    radial_dim = size(radial_multipole(ops, 0), 1)
    nchannels = length(ops.one_body.channels)
    norbitals = length(orbitals(ops))
    size(density) == (norbitals, norbitals) || throw(DimensionMismatch("density matrix has the wrong size"))

    orbital_index(channel_index, radial_index) = (channel_index - 1) * radial_dim + radial_index
    reference = zeros(Float64, norbitals, norbitals)
    radial_multipoles = [radial_multipole(ops, L) for L in 0:GaussletBases.gaunt_Lmax(ops.gaunt_table)]
    angular_kernels = [angular_kernel(ops, L) for L in 0:GaussletBases.gaunt_Lmax(ops.gaunt_table)]

    for left_radial in 1:radial_dim, right_radial in 1:radial_dim, left_channel in 1:nchannels, right_channel in 1:nchannels
        total = 0.0
        for (level_index, multipole) in enumerate(radial_multipoles), left_source in 1:nchannels, right_source in 1:nchannels
            total += multipole[left_radial, right_radial] *
                     angular_kernels[level_index][left_channel, left_source, right_channel, right_source] *
                     density[orbital_index(left_source, left_radial), orbital_index(right_source, right_radial)]
        end
        reference[orbital_index(left_channel, left_radial), orbital_index(right_channel, right_radial)] = total
    end

    return 0.5 .* (reference .+ transpose(reference))
end

function _quick_display_fixture()
    return _cached_fixture(:quick_display_fixture, () -> begin
        family = GaussletFamily(:G10)
        map = AsinhMapping(c = 0.15, s = 0.15)
        ub_spec = UniformBasisSpec(:G10; xmin = -1.0, xmax = 1.0, spacing = 1.0)
        rb_spec = RadialBasisSpec(:G10;
            count = 4,
            mapping = map,
            reference_spacing = 1.0,
            tails = 2,
            odd_even_kmax = 1,
            xgaussians = XGaussian[],
        )
        ub = build_basis(ub_spec)
        rb, grid, ops, channels, atom, ida = _quick_radial_atomic_fixture()
        tiny_rb, tiny_grid, tiny_radial_ops, tiny_ida, tiny_problem = _tiny_atomic_ida_two_electron_fixture()
        rep = basis_representation(ub)
        (
            family,
            map,
            ub_spec,
            rb_spec,
            ub,
            rb,
            grid,
            ops,
            rep,
            channels,
            atom,
            ida,
            tiny_rb,
            tiny_grid,
            tiny_radial_ops,
            tiny_ida,
            tiny_problem,
        )
    end)
end

function _slow_display_fixture()
    return _cached_fixture(:slow_display_fixture, () -> begin
        family = GaussletFamily(:G10)
        map = AsinhMapping(c = 0.15, s = 0.15)
        ub_spec = UniformBasisSpec(:G10; xmin = -1.0, xmax = 1.0, spacing = 1.0)
        hb_spec = HalfLineBasisSpec(:G10;
            xmax = 2.0,
            reference_spacing = 1.0,
            tails = 2,
            mapping = map,
        )
        rb_spec = RadialBasisSpec(:G10;
            count = 4,
            mapping = map,
            reference_spacing = 1.0,
            tails = 2,
            odd_even_kmax = 1,
            xgaussians = XGaussian[],
        )
        ub = build_basis(ub_spec)
        hb = build_basis(hb_spec)
        rb = build_basis(rb_spec)
        grid = radial_quadrature(rb; accuracy = :medium, quadrature_rmax = 6.0)
        ops = atomic_operators(rb, grid; Z = 2.0, lmax = 1)
        rep = basis_representation(ub)
        partition = basis_partition(rep, [-1.5, -0.5, 0.5, 1.5])
        hierarchy = refine_partition(hierarchical_partition(partition), 1)
        pgdg = build_leaf_pgdg(hierarchy)
        augmented_pgdg = augment_leaf_pgdg(
            pgdg;
            by_leaf = Dict(
                4 => [LeafGaussianSpec1D(relative_position = 0.5, width_scale = 0.2)],
            ),
        )
        spec = LeafGaussianSpec1D(relative_position = 0.5, width_scale = 0.2)
        global_layer = build_global_mapped_primitive_layer(
            xmin = -2.0,
            xmax = 2.0,
            mapping = map,
            reference_spacing = 0.5,
            width_scale = 1.0,
        )
        contracted_layer = contract_leaf_boxes(
            global_layer,
            refine_partition(hierarchical_partition(global_layer, [-2.5, -0.5, 0.5, 2.5]), 1);
            retained_per_leaf = 1,
        )
        channels = ylm_channels(1)
        atom = atomic_one_body_operators(ops; lmax = 1)
        ida = atomic_ida_operators(ops; lmax = 1)
        (
            family,
            map,
            ub_spec,
            hb_spec,
            rb_spec,
            ub,
            hb,
            rb,
            grid,
            ops,
            rep,
            partition,
            hierarchy,
            pgdg,
            augmented_pgdg,
            spec,
            global_layer,
            contracted_layer,
            channels,
            atom,
            ida,
        )
    end)
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

function _midpoint_reference_gaussian_factor_matrix(
    basis;
    exponent::Real,
    center::Real = 0.0,
    xmin = -20.0,
    xmax = 20.0,
    h = 0.02,
)
    points = collect((xmin + 0.5 * h):h:(xmax - 0.5 * h))
    weights = fill(h, length(points))
    values = [basis[j](x) for x in points, j in 1:length(basis)]
    factor = exp.(-Float64(exponent) .* ((points .- Float64(center)) .^ 2))
    return transpose(values) * ((weights .* factor) .* values)
end

if _test_group_enabled(:core)
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

end

if _test_group_enabled(:nested)
@testset "Cartesian nested face first primitive" begin
    function _fixed_a_nested_test_basis(count::Int; a::Float64 = 0.25, xmax::Float64 = 10.0, tail_spacing::Float64 = 10.0)
        endpoint = (count - 1) / 2
        s = asinh(xmax / a) / (endpoint - xmax / tail_spacing)
        basis = build_basis(MappedUniformBasisSpec(:G10;
            count = count,
            mapping = AsinhMapping(a = a, s = s, tail_spacing = tail_spacing),
            reference_spacing = 1.0,
        ))
        return basis, s
    end

    expansion = coulomb_gaussian_expansion(doacc = false)
    basis, s = _fixed_a_nested_test_basis(13)
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents[1:3],
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    pgdg = bundle.pgdg_intermediate
    interval = 2:(length(basis) - 1)
    side = GaussletBases._nested_doside_1d(bundle, interval, 4)

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

@testset "Cartesian nested shell first packet" begin
    function _fixed_a_nested_shell_basis(count::Int; a::Float64 = 0.25, xmax::Float64 = 10.0, tail_spacing::Float64 = 10.0)
        endpoint = (count - 1) / 2
        s = asinh(xmax / a) / (endpoint - xmax / tail_spacing)
        basis = build_basis(MappedUniformBasisSpec(:G10;
            count = count,
            mapping = AsinhMapping(a = a, s = s, tail_spacing = tail_spacing),
            reference_spacing = 1.0,
        ))
        return basis, s
    end

    expansion = coulomb_gaussian_expansion(doacc = false)
    basis, s = _fixed_a_nested_shell_basis(13)
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents[1:3],
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    pgdg = bundle.pgdg_intermediate
    term_coefficients = Float64[Float64(value) for value in expansion.coefficients[1:3]]
    interval = 2:(length(basis) - 1)
    shell = GaussletBases._nested_xy_shell_pair(
        bundle,
        interval,
        interval;
        retain_x = 4,
        retain_y = 3,
        term_coefficients = term_coefficients,
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
    @test isnothing(packet.gaussian_terms)
    @test isnothing(packet.pair_terms)
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
        exponents = expansion.exponents[1:3],
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    pgdg = bundle.pgdg_intermediate
    term_coefficients = Float64[Float64(value) for value in expansion.coefficients[1:3]]
    interval = 2:(length(basis) - 1)
    shell = GaussletBases._nested_rectangular_shell(
        bundle,
        interval,
        interval,
        interval;
        retain_xy = (4, 3),
        retain_xz = (4, 3),
        retain_yz = (4, 3),
        term_coefficients = term_coefficients,
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
            count = count,
            mapping = AsinhMapping(a = a, s = s, tail_spacing = tail_spacing),
            reference_spacing = 1.0,
        ))
        return basis, s
    end

    expansion = coulomb_gaussian_expansion(doacc = false)
    basis, s = _fixed_a_multi_face_basis(13)
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents[1:3],
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    interval = 2:(length(basis) - 1)
    term_coefficients = Float64[Float64(value) for value in expansion.coefficients[1:3]]
    shell = GaussletBases._nested_rectangular_shell(
        bundle,
        interval,
        interval,
        interval;
        retain_xy = (4, 3),
        retain_xz = (4, 3),
        retain_yz = (4, 3),
        term_coefficients = term_coefficients,
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
    @test isnothing(packet.gaussian_terms)
    @test isnothing(packet.pair_terms)
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
    @test isnothing(fixed_block.gaussian_terms)
    @test isnothing(fixed_block.pair_terms)
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
            count = count,
            mapping = white_lindsey_atomic_mapping(Z = Z, d = d, tail_spacing = tail_spacing),
            reference_spacing = 1.0,
        ))
        expansion = coulomb_gaussian_expansion(doacc = false)
        sequence = build_one_center_atomic_full_parent_shell_sequence(
            basis;
            exponents = expansion.exponents,
            gausslet_backend = :numerical_reference,
            refinement_levels = 0,
            nside = nside,
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
            count = count,
            mapping = white_lindsey_atomic_mapping(Z = Z, d = d, tail_spacing = tail_spacing),
            reference_spacing = 1.0,
        ))
        expansion = coulomb_gaussian_expansion(doacc = false)
        sequence = build_one_center_atomic_legacy_profile_shell_sequence(
            basis;
            exponents = expansion.exponents,
            gausslet_backend = :numerical_reference,
            refinement_levels = 0,
            working_box = working_box,
            nside = nside,
        )
        diagnostics = one_center_atomic_nested_structure_diagnostics(
            sequence;
            parent_side_count = count,
            nside = nside,
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
                expansion = expansion,
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

    @test count_only_5.expected_shell_increment == 5^3 - 3^3
    @test count_only_5.expected_shell_increment == 98

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
    @test isnothing(timed_full.fixed_block.gaussian_terms)
    @test isnothing(timed_full.fixed_block.pair_terms)
    @test !isnothing(timed_full.fixed_block.gaussian_sum)
    @test !isnothing(timed_full.fixed_block.pair_sum)
    @test isnothing(timed_legacy.fixed_block.gaussian_terms)
    @test isnothing(timed_legacy.fixed_block.pair_terms)
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
        expansion = expansion,
        nside = 5,
    )
    compact_legacy = one_center_atomic_legacy_profile_fixed_block(
        basis;
        expansion = expansion,
        working_box = 2:12,
        nside = 5,
    )

    @test compact_full.term_storage == :compact_production
    @test isnothing(compact_full.gaussian_terms)
    @test isnothing(compact_full.pair_terms)
    @test !isnothing(compact_full.gaussian_sum)
    @test !isnothing(compact_full.pair_sum)

    @test compact_legacy.term_storage == :compact_production
    @test isnothing(compact_legacy.gaussian_terms)
    @test isnothing(compact_legacy.pair_terms)
    @test !isnothing(compact_legacy.gaussian_sum)
    @test !isnothing(compact_legacy.pair_sum)

    @test_throws MethodError one_center_atomic_full_parent_fixed_block(
        basis;
        expansion = expansion,
        nside = 5,
        retain_term_tensors = true,
    )
    @test_throws MethodError one_center_atomic_legacy_profile_fixed_block(
        basis;
        expansion = expansion,
        working_box = 2:12,
        nside = 5,
        retain_term_tensors = true,
    )
end

@testset "Cartesian basis representation for direct-product QW bases" begin
    basis, _operators, _check = _bond_aligned_diatomic_qw_fixture()
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
end

@testset "Cartesian basis representation for nested fixed blocks" begin
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
        expansion = expansion,
        nside = 5,
    )
    representation = basis_representation(fixed_block)
    metadata = basis_metadata(representation)

    @test representation isa CartesianBasisRepresentation3D
    @test metadata.basis_kind == :nested_fixed_block
    @test metadata.parent_kind == :cartesian_product_basis
    @test metadata.parent_axis_counts == (13, 13, 13)
    @test metadata.parent_dimension == 13^3
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
end

function _with_sparse_nested_coefficients(fixed_block::GaussletBases._NestedFixedBlock3D)
    return GaussletBases._NestedFixedBlock3D(
        fixed_block.parent_basis,
        fixed_block.shell,
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
        fixed_block.gaussian_terms,
        fixed_block.pair_terms,
        fixed_block.term_storage,
        fixed_block.fixed_centers,
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
        expansion = expansion,
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
            expansion = expansion,
            nside = 5,
        )
        fixed_legacy = one_center_atomic_legacy_profile_fixed_block(
            basis;
            expansion = expansion,
            working_box = 2:12,
            nside = 5,
        )
        supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 1)
        full_ops = ordinary_cartesian_qiu_white_operators(
            fixed_full,
            supplement;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
            residual_keep_policy = :near_null_only,
        )
        legacy_ops = ordinary_cartesian_qiu_white_operators(
            fixed_legacy,
            supplement;
            expansion = expansion,
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
                mapping = mapping,
                reference_spacing = 1.0,
            ),
        )
        target_basis = build_basis(
            MappedUniformBasisSpec(
                :G10;
                count = target_count,
                mapping = mapping,
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
                mapping = mapping,
                reference_spacing = 1.0,
            ),
        )

        source_fixed = one_center_atomic_legacy_profile_fixed_block(
            parent_basis;
            expansion = expansion,
            working_box = source_working_box,
            nside = 5,
        )
        target_fixed = one_center_atomic_full_parent_fixed_block(
            parent_basis;
            expansion = expansion,
            nside = 5,
        )

        source_ops = ordinary_cartesian_qiu_white_operators(
            source_fixed,
            supplement;
            expansion = expansion,
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
            residual_keep_policy = :near_null_only,
        )
        target_ops = ordinary_cartesian_qiu_white_operators(
            target_fixed,
            supplement;
            expansion = expansion,
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
            target_working_box = target_fixed.working_box,
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

    mktempdir() do dir
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
        @test disk_hybrid_transfer.coefficients ≈
            memory_hybrid_transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test disk_hybrid_transfer.projector !== nothing
        @test memory_hybrid_transfer.projector !== nothing
        @test disk_hybrid_transfer.projector.matrix ≈
            memory_hybrid_transfer.projector.matrix atol = 1.0e-10 rtol = 1.0e-10
        @test disk_hybrid_transfer.diagnostics.transferred_residual_inf < 1.0e-10

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

    @test fixed_full_direct.overlap ≈ fixed_full_reference.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.kinetic ≈ fixed_full_reference.kinetic atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.position_x ≈ fixed_full_reference.position_x atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.position_y ≈ fixed_full_reference.position_y atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.position_z ≈ fixed_full_reference.position_z atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.x2_x ≈ fixed_full_reference.x2_x atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.x2_y ≈ fixed_full_reference.x2_y atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.x2_z ≈ fixed_full_reference.x2_z atol = 1.0e-10 rtol = 1.0e-10
    @test fixed_full_direct.weights ≈ fixed_full_reference.weights atol = 1.0e-10 rtol = 1.0e-10
    @test isnothing(fixed_full_direct.gaussian_terms)
    @test isnothing(fixed_full_direct.pair_terms)
    @test isnothing(fixed_full_reference.gaussian_terms)
    @test isnothing(fixed_full_reference.pair_terms)
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
    @test isnothing(fixed_legacy_direct.gaussian_terms)
    @test isnothing(fixed_legacy_direct.pair_terms)
    @test isnothing(fixed_legacy_reference.gaussian_terms)
    @test isnothing(fixed_legacy_reference.pair_terms)
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
        shell1_face,
        shell2_face,
        shell3_face,
        shell4_face,
        shell1_complete,
        shell2_complete,
        shell3_complete,
        shell4_complete,
        interval1,
        interval2,
        interval3,
        interval4,
        core5,
        shell_plus_core,
        face_sequence,
        complete_sequence,
        fixed_shell_plus_core,
        fixed_face_sequence,
        fixed_complete_sequence,
        legacy,
        baseline,
        shell_plus_core_ops,
        face_sequence_ops,
        complete_sequence_ops,
        baseline_check,
        shell_plus_core_check,
        face_sequence_check,
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
    @test sum(length(face.support_indices) for face in shell1_complete.faces) == 6 * 11^2
    @test sum(length(edge.support_indices) for edge in shell1_complete.edges) == 12 * 11
    @test sum(length(corner.support_indices) for corner in shell1_complete.corners) == 8
    @test isempty(intersect(shell1_face.support_indices, shell2_face.support_indices))
    @test isempty(intersect(shell2_face.support_indices, shell3_face.support_indices))
    @test isempty(intersect(shell3_face.support_indices, shell4_face.support_indices))

    @test face_sequence isa GaussletBases._CartesianNestedShellSequence3D
    @test complete_sequence isa GaussletBases._CartesianNestedShellSequence3D
    @test length(face_sequence.core_indices) == 5^3
    @test length(complete_sequence.core_indices) == 5^3
    @test complete_sequence.working_box == (3:15, 3:15, 3:15)
    @test shell_plus_core_ops.gausslet_count == 1385
    @test face_sequence_ops.gausslet_count < complete_sequence_ops.gausslet_count < shell_plus_core_ops.gausslet_count
    @test baseline.gausslet_count == 17^3
    @test norm(fixed_face_sequence.overlap - I, Inf) < 1.0e-10
    @test norm(fixed_complete_sequence.overlap - I, Inf) < 1.0e-10
    @test all(isfinite, fixed_shell_plus_core.weights)
    @test minimum(fixed_shell_plus_core.weights) > 0.0
    @test all(isfinite, fixed_complete_sequence.weights)
    @test minimum(fixed_complete_sequence.weights) > 0.0
    @test complete_sequence_check.overlap_error < 1.0e-10
    @test isfinite(complete_sequence_check.orbital_energy)
    @test isfinite(complete_sequence_check.vee_expectation)
    @test abs(complete_sequence_check.orbital_energy - baseline_check.orbital_energy) <
        abs(face_sequence_check.orbital_energy - baseline_check.orbital_energy)
    @test abs(complete_sequence_check.vee_expectation - baseline_check.vee_expectation) <
        abs(face_sequence_check.vee_expectation - baseline_check.vee_expectation)
    @test abs(complete_sequence_check.vee_expectation - baseline_check.vee_expectation) < 2.0e-4
    @test abs(complete_sequence_check.vee_expectation - shell_plus_core_check.vee_expectation) < 2.0e-4

    expansion = coulomb_gaussian_expansion(doacc = false)
    overlap_parent, one_body_parent, interaction_parent = _nested_parent_fixed_problem(bundle, expansion; Z = 2.0)
    parent_modes = eigen(Hermitian(one_body_parent), Hermitian(overlap_parent))
    parent_ground = parent_modes.vectors[:, 1]
    parent_ground_vee = _nested_vee_from_orbital(interaction_parent, parent_ground)
    projected_shell_plus_core = _nested_fixed_projected_orbital(overlap_parent, fixed_shell_plus_core, parent_ground)
    projected_complete = _nested_fixed_projected_orbital(overlap_parent, fixed_complete_sequence, parent_ground)
    projected_shell_plus_core_vee = _nested_vee_from_orbital(
        GaussletBases._qwrg_fixed_block_interaction_matrix(fixed_shell_plus_core, expansion),
        projected_shell_plus_core,
    )
    projected_complete_vee = _nested_vee_from_orbital(
        GaussletBases._qwrg_fixed_block_interaction_matrix(fixed_complete_sequence, expansion),
        projected_complete,
    )
    term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
    @test abs(projected_shell_plus_core_vee - parent_ground_vee) < 1.0e-4
    @test abs(projected_complete_vee - parent_ground_vee) < 5.0e-4
    @test_throws ArgumentError GaussletBases._nested_shell_sequence(
        bundle,
        core5,
        core5,
        core5,
        [shell1_face, shell2_face, shell3_face, shell4_face],
        term_coefficients = term_coefficients,
    )
    @test_throws ArgumentError GaussletBases._nested_shell_sequence(
        bundle,
        core5,
        core5,
        core5,
        [shell1_complete, shell2_complete, shell3_complete],
        term_coefficients = term_coefficients,
    )
end

end

if _test_group_enabled(:ordinary)
@testset "Bond-aligned homonuclear chain ordinary QW reference path" begin
    basis2, operators2, diagnostics2 = _bond_aligned_homonuclear_chain_qw_fixture(; natoms = 2, spacing = 1.4)
    basis3, operators3, diagnostics3 = _bond_aligned_homonuclear_chain_qw_fixture(; natoms = 3, spacing = 1.2)

    @test basis2 isa BondAlignedHomonuclearChainQWBasis3D
    @test basis3 isa BondAlignedHomonuclearChainQWBasis3D
    @test operators2 isa QiuWhiteResidualGaussianOperators
    @test operators3 isa QiuWhiteResidualGaussianOperators
    @test operators2.gaussian_data === nothing
    @test operators3.gaussian_data === nothing
    @test operators2.residual_count == 0
    @test operators3.residual_count == 0
    @test operators2.gausslet_count == length(basis2.basis_x) * length(basis2.basis_y) * length(basis2.basis_z)
    @test operators3.gausslet_count == length(basis3.basis_x) * length(basis3.basis_y) * length(basis3.basis_z)
    @test size(operators2.overlap, 1) < 1000
    @test size(operators3.overlap, 1) < 1000
    @test norm(operators2.overlap - I, Inf) < 1.0e-8
    @test norm(operators3.overlap - I, Inf) < 1.0e-8
    @test operators2.one_body_hamiltonian ≈ transpose(operators2.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test operators3.one_body_hamiltonian ≈ transpose(operators3.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test operators2.interaction_matrix ≈ transpose(operators2.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test operators3.interaction_matrix ≈ transpose(operators3.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test all(isfinite, operators2.one_body_hamiltonian)
    @test all(isfinite, operators3.one_body_hamiltonian)
    @test all(isfinite, operators2.interaction_matrix)
    @test all(isfinite, operators3.interaction_matrix)
    @test minimum(diag(operators2.interaction_matrix)) > 0.0
    @test minimum(diag(operators3.interaction_matrix)) > 0.0
    @test diagnostics2.axis_monotone
    @test diagnostics3.axis_monotone
    @test all(diagnostics2.local_spacings_at_midpoints .> 0.45)
    @test all(diagnostics3.local_spacings_at_midpoints .> 0.45)
    @test length(basis2.basis_z) >= length(basis2.basis_x)
    @test length(basis3.basis_z) > length(basis3.basis_x)
end

@testset "Axis-aligned homonuclear square-lattice ordinary QW reference path" begin
    basis2, operators2, diagnostics2, check2 =
        _axis_aligned_homonuclear_square_lattice_qw_fixture(; n = 2, spacing = 1.4)
    basis3, operators3, diagnostics3, check3 =
        _axis_aligned_homonuclear_square_lattice_qw_fixture(; n = 3, spacing = 1.2)

    @test basis2 isa AxisAlignedHomonuclearSquareLatticeQWBasis3D
    @test basis3 isa AxisAlignedHomonuclearSquareLatticeQWBasis3D
    @test operators2 isa QiuWhiteResidualGaussianOperators
    @test operators3 isa QiuWhiteResidualGaussianOperators
    @test operators2.gaussian_data === nothing
    @test operators3.gaussian_data === nothing
    @test operators2.residual_count == 0
    @test operators3.residual_count == 0
    @test operators2.gausslet_count == length(basis2.basis_x) * length(basis2.basis_y) * length(basis2.basis_z)
    @test operators3.gausslet_count == length(basis3.basis_x) * length(basis3.basis_y) * length(basis3.basis_z)
    @test size(operators2.overlap, 1) < 1500
    @test size(operators3.overlap, 1) < 4000
    @test norm(operators2.overlap - I, Inf) < 1.0e-8
    @test norm(operators3.overlap - I, Inf) < 1.0e-8
    @test operators2.one_body_hamiltonian ≈ transpose(operators2.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test operators3.one_body_hamiltonian ≈ transpose(operators3.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test operators2.interaction_matrix ≈ transpose(operators2.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test operators3.interaction_matrix ≈ transpose(operators3.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test all(isfinite, operators2.one_body_hamiltonian)
    @test all(isfinite, operators3.one_body_hamiltonian)
    @test all(isfinite, operators2.interaction_matrix)
    @test all(isfinite, operators3.interaction_matrix)
    @test minimum(diag(operators2.interaction_matrix)) > 0.0
    @test minimum(diag(operators3.interaction_matrix)) > 0.0
    @test diagnostics2.x_axis_monotone
    @test diagnostics2.y_axis_monotone
    @test diagnostics3.x_axis_monotone
    @test diagnostics3.y_axis_monotone
    @test diagnostics2.xy_axis_center_match_error < 1.0e-12
    @test diagnostics3.xy_axis_center_match_error < 1.0e-12
    @test diagnostics2.local_spacing_at_plane_center_x > 0.45
    @test diagnostics2.local_spacing_at_plane_center_y > 0.45
    @test all(diagnostics3.representative_midpoint_spacings_x .> 0.45)
    @test all(diagnostics3.representative_midpoint_spacings_y .> 0.45)
    @test isfinite(check2.orbital_energy)
    @test isfinite(check3.orbital_energy)
    @test isfinite(check2.vee_expectation)
    @test isfinite(check3.vee_expectation)
end

@testset "Axis-aligned homonuclear square-lattice nested geometry diagnostics" begin
    basis2, source2, fixed2, diagnostics2 =
        _axis_aligned_homonuclear_square_lattice_nested_fixture(; n = 2, spacing = 1.4)
    basis3, source3, fixed3, diagnostics3 =
        _axis_aligned_homonuclear_square_lattice_nested_fixture(; n = 3, spacing = 1.2)

    @test basis2 isa AxisAlignedHomonuclearSquareLatticeQWBasis3D
    @test basis3 isa AxisAlignedHomonuclearSquareLatticeQWBasis3D
    @test diagnostics2.retention_contract.retain_xy == (3, 3)
    @test diagnostics3.retention_contract.retain_xy == (3, 3)
    @test diagnostics2.retention_contract.retain_x_edge == 3
    @test diagnostics3.retention_contract.retain_x_edge == 3
    @test diagnostics2.retention_contract.shell_increment == 98
    @test diagnostics3.retention_contract.shell_increment == 98
    @test diagnostics2.retention_contract.matches_nside_default
    @test diagnostics3.retention_contract.matches_nside_default
    @test diagnostics2.contract_audit.full_parent_working_box
    @test diagnostics3.contract_audit.full_parent_working_box
    @test diagnostics2.contract_audit.missing_row_count == 0
    @test diagnostics3.contract_audit.missing_row_count == 0
    @test diagnostics2.contract_audit.ownership_unowned_row_count == 0
    @test diagnostics3.contract_audit.ownership_unowned_row_count == 0
    @test diagnostics2.contract_audit.ownership_multi_owned_row_count == 0
    @test diagnostics3.contract_audit.ownership_multi_owned_row_count == 0
    @test diagnostics2.shared_shells_match_contract
    @test diagnostics3.shared_shells_match_contract
    @test diagnostics2.leaf_count == 4
    @test diagnostics3.leaf_count == 9
    @test diagnostics2.root_node.did_split
    @test diagnostics3.root_node.did_split
    @test diagnostics2.root_node.accepted_candidate_index == 1
    @test diagnostics3.root_node.accepted_candidate_index == 1
    @test diagnostics2.root_node.child_count == 2
    @test diagnostics3.root_node.child_count == 3
    @test !diagnostics2.root_node.local_resolution_warning
    @test !diagnostics3.root_node.local_resolution_warning
    @test diagnostics2.root_node.candidate_summaries[1].split_family == :split_x_binary
    @test diagnostics2.root_node.candidate_summaries[2].split_family == :split_y_binary
    @test diagnostics3.root_node.candidate_summaries[1].split_family == :split_x_ternary
    @test diagnostics3.root_node.candidate_summaries[2].split_family == :split_y_ternary
    @test diagnostics2.root_node.candidate_summaries[1].symmetry_preserving
    @test diagnostics2.root_node.candidate_summaries[2].symmetry_preserving
    @test diagnostics3.root_node.candidate_summaries[1].symmetry_preserving
    @test diagnostics3.root_node.candidate_summaries[2].symmetry_preserving
    @test all(
        ratio >= diagnostics2.root_node.min_in_plane_aspect_ratio for
        ratio in diagnostics2.root_node.candidate_summaries[1].child_in_plane_aspect_ratios
    )
    @test all(
        ratio >= diagnostics3.root_node.min_in_plane_aspect_ratio for
        ratio in diagnostics3.root_node.candidate_summaries[1].child_in_plane_aspect_ratios
    )
    @test length(source2.root_geometry.child_nodes) == 2
    @test length(source3.root_geometry.child_nodes) == 3
    @test source2.root_geometry.child_nodes[1].accepted_candidate_index == 1
    @test source2.root_geometry.child_nodes[2].accepted_candidate_index == 1
    @test all(child.accepted_candidate_index == 1 for child in source3.root_geometry.child_nodes)
    @test size(fixed2.overlap, 1) == diagnostics2.fixed_dimension
    @test size(fixed3.overlap, 1) == diagnostics3.fixed_dimension
    @test norm(fixed2.overlap - I, Inf) < 1.0e-8
    @test norm(fixed3.overlap - I, Inf) < 1.0e-8
    @test all(isfinite, fixed2.weights)
    @test all(isfinite, fixed3.weights)
    @test minimum(fixed2.weights) > 0.0
    @test minimum(fixed3.weights) > 0.0
    @test all(isfinite, fixed2.fixed_centers)
    @test all(isfinite, fixed3.fixed_centers)

    mktempdir() do dir
        report_path = joinpath(dir, "square_lattice_nested_geometry_report.txt")
        written = write_axis_aligned_homonuclear_square_lattice_nested_geometry_report(
            report_path,
            basis3;
            nside = 5,
            min_in_plane_aspect_ratio = 0.15,
        )
        @test written.leaf_count == diagnostics3.leaf_count
        report_text = read(report_path, String)
        @test occursin("# retain_xy = (3, 3)", report_text)
        @test occursin("# shell_increment = 98", report_text)
        @test occursin("# ownership_multi_owned_row_count = 0", report_text)
        @test occursin("candidate[1].split_family = split_x_ternary", report_text)
        @test occursin("candidate[2].split_family = split_y_ternary", report_text)
        @test occursin("[node root_1]", report_text)
    end
end

@testset "Experimental axis-aligned homonuclear square-lattice nested QW consumer path" begin
    basis2, path2, check2 =
        _axis_aligned_homonuclear_square_lattice_nested_qw_fixture(; n = 2, spacing = 1.4)
    basis3, path3, check3 =
        _axis_aligned_homonuclear_square_lattice_nested_qw_fixture(; n = 3, spacing = 1.2)

    for (basis, path, check, n, expected_dim) in (
        (basis2, path2, check2, 2, 773),
        (basis3, path3, check3, 3, 703),
    )
        operators = path.operators
        fixed_block = path.fixed_block
        diagnostics = path.diagnostics
        @test basis isa AxisAlignedHomonuclearSquareLatticeQWBasis3D
        @test path isa GaussletBases.ExperimentalAxisAlignedHomonuclearSquareLatticeNestedQWPath
        @test path.basis === basis
        @test path.source === diagnostics.source
        @test path.fixed_block === fixed_block
        @test path.min_in_plane_aspect_ratio == 0.15
        @test operators isa GaussletBases.QiuWhiteResidualGaussianOperators
        @test operators.basis === fixed_block
        @test operators.gaussian_data === nothing
        @test operators.interaction_treatment == :ggt_nearest
        @test operators.residual_count == 0
        @test operators.gausslet_count == size(fixed_block.overlap, 1)
        @test size(operators.residual_centers) == (0, 3)
        @test size(operators.residual_widths) == (0, 3)
        @test norm(operators.overlap - I, Inf) < 1.0e-8
        @test norm(operators.one_body_hamiltonian - transpose(operators.one_body_hamiltonian), Inf) < 1.0e-10
        @test norm(operators.interaction_matrix - transpose(operators.interaction_matrix), Inf) < 1.0e-10
        @test all(isfinite, operators.one_body_hamiltonian)
        @test all(isfinite, operators.interaction_matrix)
        @test all(isfinite, fixed_block.weights)
        @test minimum(fixed_block.weights) > 0.0
        @test all(isfinite, fixed_block.fixed_centers)
        @test size(fixed_block.overlap, 1) == diagnostics.fixed_dimension
        @test size(fixed_block.overlap, 1) == expected_dim
        @test check.overlap_error < 1.0e-8
        @test isfinite(check.orbital_energy)
        @test isfinite(check.vee_expectation)
        @test length(basis.nuclei) == n * n
    end

    @test path2.diagnostics.root_node.did_split
    @test path2.diagnostics.root_node.accepted_candidate_index == 1
    @test path2.diagnostics.leaf_count == 4
    @test path2.diagnostics.root_node.child_count == 2
    @test path2.diagnostics.root_node.min_in_plane_aspect_ratio == 0.15
    @test path2.diagnostics.root_node.candidate_summaries[1].split_family == :split_x_binary
    @test path2.diagnostics.root_node.candidate_summaries[1].did_split

    @test path3.diagnostics.root_node.did_split
    @test path3.diagnostics.root_node.accepted_candidate_index == 1
    @test path3.diagnostics.leaf_count == 9
    @test path3.diagnostics.root_node.child_count == 3
    @test path3.diagnostics.root_node.min_in_plane_aspect_ratio == 0.15
    @test path3.diagnostics.root_node.candidate_summaries[1].split_family == :split_x_ternary
    @test path3.diagnostics.root_node.candidate_summaries[2].split_family == :split_y_ternary
    @test path3.diagnostics.root_node.candidate_summaries[1].did_split
    @test path3.diagnostics.root_node.candidate_summaries[1].child_in_plane_aspect_ratios[2] < 0.2
    @test path3.diagnostics.root_node.candidate_summaries[1].child_in_plane_aspect_ratios[2] >=
        path3.min_in_plane_aspect_ratio
    @test !path3.diagnostics.root_node.local_resolution_warning
end

@testset "Experimental homonuclear square-lattice nested dense export" begin
    basis2, path2, _check2 =
        _axis_aligned_homonuclear_square_lattice_nested_qw_fixture(; n = 2, spacing = 1.4)
    basis3, path3, _check3 =
        _axis_aligned_homonuclear_square_lattice_nested_qw_fixture(; n = 3, spacing = 1.2)

    for (basis, path, n, natoms, expected_family) in (
        (basis2, path2, 2, 4, "split_x_binary"),
        (basis3, path3, 3, 9, "split_x_ternary"),
    )
        payload_data = experimental_homonuclear_square_lattice_nested_dense_payload(
            path;
            meta = (example = "test_experimental_square_lattice_nested_dense_export",),
        )
        @test payload_data.bridge_meta["format"] ==
            "experimental_homonuclear_square_lattice_nested_dense_v1"
        @test payload_data.bridge_meta["experimental"]
        @test payload_data.bridge_meta["lattice_size"] == n
        @test payload_data.bridge_meta["fixed_dimension"] == size(path.fixed_block.overlap, 1)
        @test payload_data.bridge_meta["leaf_count"] == path.diagnostics.leaf_count
        @test payload_data.bridge_meta["root_accepted_candidate_index"] ==
            something(path.diagnostics.root_node.accepted_candidate_index, 0)
        @test payload_data.bridge_meta["root_accepted_split_family"] == expected_family
        @test payload_data.bridge_meta["root_child_count"] == path.diagnostics.root_node.child_count
        @test payload_data.bridge_meta["min_in_plane_aspect_ratio"] == path.min_in_plane_aspect_ratio
        @test payload_data.bridge_meta["residual_sector_empty"]
        @test payload_data.bridge_meta["coordinate_provenance"] == "uniform_square_spacing"
        @test size(payload_data.payload["S"]) == size(path.operators.overlap)
        @test size(payload_data.payload["H1"]) == size(path.operators.one_body_hamiltonian)
        @test size(payload_data.payload["Vee"]) == size(path.operators.interaction_matrix)
        @test size(payload_data.payload["basis_centers"]) == size(path.fixed_block.fixed_centers)
        @test size(payload_data.payload["nuclear_coordinates_xyz"]) == (natoms, 3)
        @test payload_data.payload["nuclear_charges"] == fill(1.0, natoms)
        @test payload_data.payload["lattice_x_coordinates"] == basis.x_coordinates
        @test payload_data.payload["lattice_y_coordinates"] == basis.y_coordinates
        @test payload_data.payload["orbital_labels"] ==
            String[orbital.label for orbital in orbitals(path.operators)]
        @test payload_data.payload["geometry_report_text"] isa String
        @test occursin("min_in_plane_aspect_ratio = 0.15", payload_data.payload["geometry_report_text"])
        @test occursin("[node root]", payload_data.payload["geometry_report_text"])
        @test size(payload_data.payload["root_accepted_child_planar_counts"], 2) == 2
        @test size(payload_data.payload["root_accepted_child_physical_widths"], 2) == 3
        @test payload_data.meta_values["manifest/contract/status"] == "experimental"
        @test payload_data.meta_values["manifest/source/lattice_size"] == n
        @test payload_data.meta_values["manifest/source/min_in_plane_aspect_ratio"] ==
            path.min_in_plane_aspect_ratio
        @test payload_data.meta_values["manifest/source/fixed_dimension"] ==
            path.diagnostics.fixed_dimension
    end

    payload3 = experimental_homonuclear_square_lattice_nested_dense_payload(
        path3;
        meta = (example = "test_experimental_square_lattice_nested_dense_export",),
    )
    @test payload3.payload["root_accepted_child_in_plane_aspect_ratios"][2] < 0.2

    mktempdir() do dir
        export_path = joinpath(dir, "h3_square_lattice_nested_dense.jld2")
        @test write_experimental_homonuclear_square_lattice_nested_dense_jld2(
            export_path,
            path3;
            meta = (example = "test_experimental_square_lattice_nested_dense_export",),
        ) == export_path
        jldopen(export_path, "r") do file
            @test String(file["bridge/format"]) ==
                "experimental_homonuclear_square_lattice_nested_dense_v1"
            @test Bool(file["bridge/experimental"])
            @test Int(file["bridge/lattice_size"]) == 3
            @test String(file["bridge/root_accepted_split_family"]) == "split_x_ternary"
            @test String(file["bridge/root_accepted_split_axis"]) == "x"
            @test Int(file["bridge/root_accepted_candidate_index"]) == 1
            @test Int(file["bridge/root_child_count"]) == 3
            @test Float64(file["bridge/min_in_plane_aspect_ratio"]) == 0.15
            @test Bool(file["bridge/residual_sector_empty"])
            @test size(file["S"]) == size(path3.operators.overlap)
            @test size(file["H1"]) == size(path3.operators.one_body_hamiltonian)
            @test size(file["Vee"]) == size(path3.operators.interaction_matrix)
            @test size(file["nuclear_coordinates_xyz"]) == (9, 3)
            @test size(file["root_accepted_child_planar_counts"]) == (3, 2)
            @test Float64(file["root_accepted_child_in_plane_aspect_ratios"][2]) < 0.2
            @test String(file["meta/producer"]) ==
                "GaussletBases.write_experimental_homonuclear_square_lattice_nested_dense_jld2"
            @test String(file["meta/manifest/contract/status"]) == "experimental"
            @test Int(file["meta/manifest/source/lattice_size"]) == 3
            @test Float64(file["meta/manifest/source/min_in_plane_aspect_ratio"]) == 0.15
            @test String(file["meta/example"]) ==
                "test_experimental_square_lattice_nested_dense_export"
        end
    end
end

@testset "Bond-aligned homonuclear chain nested geometry diagnostics" begin
    basis3, source3, fixed3, diagnostics3 = _bond_aligned_homonuclear_chain_nested_fixture(; natoms = 3, spacing = 1.2)
    basis4, source4, fixed4, diagnostics4 = _bond_aligned_homonuclear_chain_nested_fixture(; natoms = 4, spacing = 1.2)
    basis5, source5, fixed5, diagnostics5 = _bond_aligned_homonuclear_chain_nested_fixture(; natoms = 5, spacing = 1.2)
    basis3r, source3r, fixed3r, diagnostics3r = _bond_aligned_homonuclear_chain_nested_fixture(; natoms = 3, spacing = 1.2, odd_chain_policy = :central_ternary_relaxed)
    basis5r, source5r, fixed5r, diagnostics5r = _bond_aligned_homonuclear_chain_nested_fixture(; natoms = 5, spacing = 1.2, odd_chain_policy = :central_ternary_relaxed)

    @test basis3 isa BondAlignedHomonuclearChainQWBasis3D
    @test basis4 isa BondAlignedHomonuclearChainQWBasis3D
    @test basis5 isa BondAlignedHomonuclearChainQWBasis3D
    @test diagnostics3.retention_contract.retain_xy == (3, 3)
    @test diagnostics4.retention_contract.retain_xy == (3, 3)
    @test diagnostics5.retention_contract.retain_xy == (3, 3)
    @test diagnostics3.retention_contract.retain_x_edge == 3
    @test diagnostics4.retention_contract.retain_x_edge == 3
    @test diagnostics5.retention_contract.retain_x_edge == 3
    @test diagnostics3.retention_contract.shell_increment == 98
    @test diagnostics4.retention_contract.shell_increment == 98
    @test diagnostics5.retention_contract.shell_increment == 98
    @test diagnostics3.retention_contract.matches_nside_default
    @test diagnostics4.retention_contract.matches_nside_default
    @test diagnostics5.retention_contract.matches_nside_default
    @test diagnostics3.contract_audit.full_parent_working_box
    @test diagnostics4.contract_audit.full_parent_working_box
    @test diagnostics5.contract_audit.full_parent_working_box
    @test diagnostics3.contract_audit.missing_row_count == 0
    @test diagnostics4.contract_audit.missing_row_count == 0
    @test diagnostics5.contract_audit.missing_row_count == 0
    @test diagnostics3.contract_audit.ownership_unowned_row_count == 0
    @test diagnostics4.contract_audit.ownership_unowned_row_count == 0
    @test diagnostics5.contract_audit.ownership_unowned_row_count == 0
    @test diagnostics3.contract_audit.ownership_multi_owned_row_count == 0
    @test diagnostics4.contract_audit.ownership_multi_owned_row_count == 0
    @test diagnostics5.contract_audit.ownership_multi_owned_row_count == 0
    @test diagnostics3.shared_shells_match_contract
    @test diagnostics4.shared_shells_match_contract
    @test diagnostics5.shared_shells_match_contract
    @test diagnostics3.leaf_count == 1
    @test diagnostics4.leaf_count == 2
    @test diagnostics5.leaf_count == 1
    @test !diagnostics3.root_node.did_split
    @test diagnostics4.root_node.did_split
    @test !diagnostics5.root_node.did_split
    @test diagnostics3.root_node.child_count == 0
    @test diagnostics4.root_node.child_count == 2
    @test diagnostics5.root_node.child_count == 0
    @test diagnostics3.root_node.odd_chain_policy == :strict_current
    @test diagnostics5.root_node.odd_chain_policy == :strict_current
    @test diagnostics3.root_node.accepted_candidate_index === nothing
    @test diagnostics4.root_node.accepted_candidate_index == 2
    @test diagnostics5.root_node.accepted_candidate_index === nothing
    @test diagnostics3.root_node.candidate_summaries[1].split_kind == :ternary
    @test diagnostics4.root_node.candidate_summaries[2].split_kind == :binary
    @test diagnostics5.root_node.candidate_summaries[2].split_kind == :ternary
    @test diagnostics3.root_node.candidate_summaries[1].nucleus_ranges == [1:1, 2:2, 3:3]
    @test diagnostics4.root_node.candidate_summaries[2].nucleus_ranges == [1:2, 3:4]
    @test diagnostics5.root_node.candidate_summaries[2].nucleus_ranges == [1:2, 3:3, 4:5]
    @test diagnostics3.root_node.local_resolution_warning
    @test !diagnostics4.root_node.local_resolution_warning
    @test diagnostics5.root_node.local_resolution_warning
    @test length(source3.root_geometry.child_nodes) == 0
    @test length(source4.root_geometry.child_nodes) == 2
    @test length(source5.root_geometry.child_nodes) == 0
    @test !diagnostics3.root_node.candidate_summaries[1].did_split
    @test diagnostics4.root_node.candidate_summaries[2].did_split
    @test diagnostics5.root_node.candidate_summaries[2].count_eligible
    @test !diagnostics5.root_node.candidate_summaries[2].shape_eligible
    @test diagnostics5.root_node.candidate_summaries[2].child_parallel_counts == [6, 3, 6]
    @test diagnostics5.root_node.candidate_summaries[2].child_parallel_to_transverse_ratios[2] <
        diagnostics5.root_node.odd_chain_policy_thresholds.center_parallel_to_transverse_ratio_min
    @test diagnostics4.root_node.candidate_summaries[2].midpoint_values == [0.0]
    @test abs(
        length(diagnostics4.root_node.candidate_summaries[2].child_boxes[1][3]) -
        length(diagnostics4.root_node.candidate_summaries[2].child_boxes[2][3]),
    ) <= 1
    @test all(
        widths[3] > 0.0 for
        widths in diagnostics4.root_node.candidate_summaries[2].child_physical_widths
    )
    @test size(fixed3.overlap, 1) == diagnostics3.fixed_dimension
    @test size(fixed4.overlap, 1) == diagnostics4.fixed_dimension
    @test size(fixed5.overlap, 1) == diagnostics5.fixed_dimension
    @test norm(fixed3.overlap - I, Inf) < 1.0e-8
    @test norm(fixed4.overlap - I, Inf) < 1.0e-8
    @test norm(fixed5.overlap - I, Inf) < 1.0e-8
    @test all(isfinite, fixed3.weights)
    @test all(isfinite, fixed4.weights)
    @test all(isfinite, fixed5.weights)
    @test all(isfinite, fixed3.fixed_centers)
    @test all(isfinite, fixed4.fixed_centers)
    @test all(isfinite, fixed5.fixed_centers)

    @test diagnostics3r.root_node.odd_chain_policy == :central_ternary_relaxed
    @test diagnostics5r.root_node.odd_chain_policy == :central_ternary_relaxed
    @test diagnostics3r.root_node.did_split
    @test diagnostics5r.root_node.did_split
    @test diagnostics3r.root_node.accepted_candidate_index == 1
    @test diagnostics5r.root_node.accepted_candidate_index == 2
    @test diagnostics3r.root_node.child_count == 3
    @test diagnostics5r.root_node.child_count == 3
    @test diagnostics3r.root_node.candidate_summaries[1].did_split
    @test diagnostics5r.root_node.candidate_summaries[2].did_split
    @test diagnostics3r.root_node.candidate_summaries[1].child_parallel_counts == [4, 3, 4]
    @test diagnostics5r.root_node.candidate_summaries[2].child_parallel_counts == [6, 3, 6]
    @test diagnostics3r.root_node.odd_chain_policy_thresholds.center_parallel_count_min == 3
    @test diagnostics5r.root_node.odd_chain_policy_thresholds.center_parallel_to_transverse_ratio_min == 0.35
    @test size(fixed3r.overlap, 1) == diagnostics3r.fixed_dimension
    @test size(fixed5r.overlap, 1) == diagnostics5r.fixed_dimension
    @test norm(fixed3r.overlap - I, Inf) < 1.0e-8
    @test norm(fixed5r.overlap - I, Inf) < 1.0e-8
    @test all(isfinite, fixed3r.weights)
    @test all(isfinite, fixed5r.weights)
    @test all(isfinite, fixed3r.fixed_centers)
    @test all(isfinite, fixed5r.fixed_centers)
    @test length(source3r.root_geometry.child_nodes) == 3
    @test length(source5r.root_geometry.child_nodes) == 3

    mktempdir() do dir
        report_path = joinpath(dir, "chain_nested_report.txt")
        written = write_bond_aligned_homonuclear_chain_nested_geometry_report(
            report_path,
            basis5r;
            nside = 5,
            odd_chain_policy = :central_ternary_relaxed,
        )
        report_text = read(report_path, String)
        @test written.leaf_count == diagnostics5r.leaf_count
        @test occursin("[node root]", report_text)
        @test occursin("# retain_xy = (3, 3)", report_text)
        @test occursin("# shell_increment = 98", report_text)
        @test occursin("# ownership_multi_owned_row_count = 0", report_text)
        @test occursin("odd_chain_policy = central_ternary_relaxed", report_text)
        @test occursin("candidate[2].child_parallel_counts = [6, 3, 6]", report_text)
        @test occursin("candidate[2].accepted = true", report_text)
    end
end

@testset "Non-atomic nested routes default to nside-driven complete-shell retention" begin
    diatomic_basis = bond_aligned_homonuclear_qw_basis(
        bond_length = 1.4,
        core_spacing = 0.5,
        xmax_parallel = 4.0,
        xmax_transverse = 3.0,
        bond_axis = :z,
    )
    diatomic_diagnostics = bond_aligned_diatomic_nested_geometry_diagnostics(
        diatomic_basis;
        nside = 7,
    )
    @test diatomic_diagnostics.nside == 7
    @test diatomic_diagnostics.child_shell_retention_contract.retain_xy == (5, 5)
    @test diatomic_diagnostics.child_shell_retention_contract.retain_x_edge == 5
    @test diatomic_diagnostics.child_shell_retention_contract.shell_increment == 218
    @test diatomic_diagnostics.child_shell_retention_contract.matches_nside_default
    @test diatomic_diagnostics.shared_shell_retention_contract.matches_nside_default
    @test diatomic_diagnostics.contract_audit.support_count ==
        diatomic_diagnostics.contract_audit.expected_support_count
    @test diatomic_diagnostics.contract_audit.missing_row_count == 0
    @test diatomic_diagnostics.contract_audit.ownership_unowned_row_count == 0
    @test diatomic_diagnostics.contract_audit.ownership_multi_owned_row_count == 0

    chain_basis = bond_aligned_homonuclear_chain_qw_basis(
        natoms = 4,
        spacing = 1.2,
        core_spacing = 0.5,
        xmax_parallel = 2.0,
        xmax_transverse = 3.0,
        chain_axis = :z,
    )
    chain_diagnostics = bond_aligned_homonuclear_chain_nested_geometry_diagnostics(
        chain_basis;
        nside = 7,
    )
    @test chain_diagnostics.nside == 7
    @test chain_diagnostics.retention_contract.retain_xy == (5, 5)
    @test chain_diagnostics.retention_contract.retain_x_edge == 5
    @test chain_diagnostics.retention_contract.shell_increment == 218
    @test chain_diagnostics.retention_contract.matches_nside_default
    @test chain_diagnostics.shared_shell_dimensions == [218]
    @test chain_diagnostics.shared_shells_match_contract
    @test chain_diagnostics.contract_audit.full_parent_working_box
    @test chain_diagnostics.contract_audit.support_count ==
        chain_diagnostics.contract_audit.expected_support_count
    @test chain_diagnostics.contract_audit.missing_row_count == 0
    @test chain_diagnostics.contract_audit.ownership_unowned_row_count == 0
    @test chain_diagnostics.contract_audit.ownership_multi_owned_row_count == 0

    square_basis = axis_aligned_homonuclear_square_lattice_qw_basis(
        n = 3,
        spacing = 1.2,
        core_spacing = 0.5,
        xmax_in_plane = 3.0,
        xmax_transverse = 3.0,
    )
    square_diagnostics = axis_aligned_homonuclear_square_lattice_nested_geometry_diagnostics(
        square_basis;
        nside = 7,
    )
    @test square_diagnostics.nside == 7
    @test square_diagnostics.retention_contract.retain_xy == (5, 5)
    @test square_diagnostics.retention_contract.retain_x_edge == 5
    @test square_diagnostics.retention_contract.shell_increment == 218
    @test square_diagnostics.retention_contract.matches_nside_default
    @test square_diagnostics.shared_shell_dimensions == [218]
    @test square_diagnostics.shared_shells_match_contract
    @test square_diagnostics.contract_audit.full_parent_working_box
    @test square_diagnostics.contract_audit.support_count ==
        square_diagnostics.contract_audit.expected_support_count
    @test square_diagnostics.contract_audit.missing_row_count == 0
    @test square_diagnostics.contract_audit.ownership_unowned_row_count == 0
    @test square_diagnostics.contract_audit.ownership_multi_owned_row_count == 0
end

@testset "Bond-aligned diatomic nested source reuse path" begin
    basis = bond_aligned_homonuclear_qw_basis(
        bond_length = 1.4,
        core_spacing = 0.5,
        xmax_parallel = 4.0,
        xmax_transverse = 3.0,
        bond_axis = :z,
    )

    diagnostics_via_basis = bond_aligned_diatomic_nested_geometry_diagnostics(
        basis;
        nside = 5,
    )
    source = bond_aligned_diatomic_nested_fixed_source(
        basis;
        nside = 5,
    )
    diagnostics_via_source = bond_aligned_diatomic_nested_geometry_diagnostics(source)
    fixed_via_basis = bond_aligned_diatomic_nested_fixed_block(
        basis;
        nside = 5,
    )
    fixed_via_source = bond_aligned_diatomic_nested_fixed_block(source)
    source_payload = bond_aligned_diatomic_source_geometry_payload(source)
    source_slice = bond_aligned_diatomic_plane_slice(
        source_payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-5,
    )

    @test diagnostics_via_source.nside == diagnostics_via_basis.nside
    @test diagnostics_via_source.geometry.did_split == diagnostics_via_basis.geometry.did_split
    @test diagnostics_via_source.geometry.count_eligible ==
        diagnostics_via_basis.geometry.count_eligible
    @test diagnostics_via_source.geometry.unsplit_aspect_eligible ==
        diagnostics_via_basis.geometry.unsplit_aspect_eligible
    @test diagnostics_via_source.geometry.shape_eligible ==
        diagnostics_via_basis.geometry.shape_eligible
    @test diagnostics_via_source.geometry.split_index ==
        diagnostics_via_basis.geometry.split_index
    @test diagnostics_via_source.geometry.working_box ==
        diagnostics_via_basis.geometry.working_box
    @test diagnostics_via_source.geometry.shared_midpoint_box ==
        diagnostics_via_basis.geometry.shared_midpoint_box
    @test diagnostics_via_source.geometry.child_boxes ==
        diagnostics_via_basis.geometry.child_boxes
    @test maximum(
        abs,
        reduce(vcat, (
            collect(widths_source .- widths_basis) for
            (widths_source, widths_basis) in zip(
                diagnostics_via_source.geometry.child_physical_widths,
                diagnostics_via_basis.geometry.child_physical_widths,
            )
        );
        init = Float64[]),
    ) < 1.0e-12
    @test diagnostics_via_source.child_shell_retention_contract ==
        diagnostics_via_basis.child_shell_retention_contract
    @test diagnostics_via_source.shared_shell_retention_contract ==
        diagnostics_via_basis.shared_shell_retention_contract
    @test diagnostics_via_source.shared_shell_dimensions ==
        diagnostics_via_basis.shared_shell_dimensions
    @test diagnostics_via_source.child_sequence_dimensions ==
        diagnostics_via_basis.child_sequence_dimensions
    @test diagnostics_via_source.fixed_dimension == diagnostics_via_basis.fixed_dimension
    @test diagnostics_via_source.contract_audit.support_count ==
        diagnostics_via_basis.contract_audit.support_count
    @test diagnostics_via_source.contract_audit.expected_support_count ==
        diagnostics_via_basis.contract_audit.expected_support_count
    @test diagnostics_via_source.contract_audit.missing_row_count ==
        diagnostics_via_basis.contract_audit.missing_row_count
    @test source_payload.bond_axis == basis.bond_axis
    @test length(source_payload.points) == prod(length.(source.geometry.parent_box))
    @test source_slice.selected_count > 0
    @test source.sequence.packet.term_storage == :compact_production
    @test isnothing(source.sequence.packet.gaussian_terms)
    @test isnothing(source.sequence.packet.pair_terms)
    @test !isnothing(source.sequence.packet.gaussian_sum)
    @test !isnothing(source.sequence.packet.pair_sum)
    @test all(isnothing(sequence.support_states) for sequence in source.child_sequences)
    @test all(isnothing(sequence.packet) for sequence in source.child_sequences)
    @test !isnothing(source.sequence.support_states)
    @test !isnothing(source.sequence.packet)

    @test fixed_via_source.source === source
    @test norm(fixed_via_source.fixed_block.overlap - fixed_via_basis.fixed_block.overlap, Inf) <
        1.0e-12
    @test norm(
        fixed_via_source.fixed_block.coefficient_matrix -
        fixed_via_basis.fixed_block.coefficient_matrix,
        Inf,
    ) < 1.0e-12
    @test fixed_via_source.fixed_block.support_indices == fixed_via_basis.fixed_block.support_indices
    @test fixed_via_source.fixed_block.working_box == fixed_via_basis.fixed_block.working_box
end

@testset "Experimental bond-aligned homonuclear chain nested QW consumer path" begin
    basis3, path3, check3 = _bond_aligned_homonuclear_chain_nested_qw_fixture(; natoms = 3, spacing = 1.2)
    basis4, path4, check4 = _bond_aligned_homonuclear_chain_nested_qw_fixture(; natoms = 4, spacing = 1.2)
    basis5, path5, check5 = _bond_aligned_homonuclear_chain_nested_qw_fixture(; natoms = 5, spacing = 1.2)

    for (basis, path, check, natoms) in ((basis3, path3, check3, 3), (basis4, path4, check4, 4), (basis5, path5, check5, 5))
        operators = path.operators
        fixed_block = path.fixed_block
        diagnostics = path.diagnostics
        @test basis isa BondAlignedHomonuclearChainQWBasis3D
        @test path isa GaussletBases.ExperimentalBondAlignedHomonuclearChainNestedQWPath
        @test path.basis === basis
        @test path.source === diagnostics.source
        @test path.fixed_block === fixed_block
        @test operators isa GaussletBases.QiuWhiteResidualGaussianOperators
        @test operators.basis === fixed_block
        @test operators.gaussian_data === nothing
        @test operators.interaction_treatment == :ggt_nearest
        @test operators.residual_count == 0
        @test operators.gausslet_count == size(fixed_block.overlap, 1)
        @test size(operators.residual_centers) == (0, 3)
        @test size(operators.residual_widths) == (0, 3)
        @test norm(operators.overlap - I, Inf) < 1.0e-8
        @test norm(operators.one_body_hamiltonian - transpose(operators.one_body_hamiltonian), Inf) < 1.0e-10
        @test norm(operators.interaction_matrix - transpose(operators.interaction_matrix), Inf) < 1.0e-10
        @test all(isfinite, operators.one_body_hamiltonian)
        @test all(isfinite, operators.interaction_matrix)
        @test all(isfinite, fixed_block.weights)
        @test minimum(fixed_block.weights) > 0.0
        @test all(isfinite, fixed_block.fixed_centers)
        @test diagnostics.root_node.odd_chain_policy == :central_ternary_relaxed
        @test path.odd_chain_policy == :central_ternary_relaxed
        @test size(fixed_block.overlap, 1) == diagnostics.fixed_dimension
        @test size(fixed_block.overlap, 1) < 1000
        @test check.overlap_error < 1.0e-8
        @test isfinite(check.orbital_energy)
        @test isfinite(check.vee_expectation)
        @test length(basis.nuclei) == natoms
    end

    @test path3.diagnostics.root_node.did_split
    @test path3.diagnostics.root_node.accepted_candidate_index == 1
    @test path3.diagnostics.leaf_count == 3
    @test path3.diagnostics.root_node.child_count == 3

    @test path4.diagnostics.root_node.did_split
    @test path4.diagnostics.root_node.accepted_candidate_index == 2
    @test path4.diagnostics.leaf_count == 2
    @test path4.diagnostics.root_node.child_count == 2

    @test path5.diagnostics.root_node.did_split
    @test path5.diagnostics.root_node.accepted_candidate_index == 2
    @test path5.diagnostics.leaf_count == 3
    @test path5.diagnostics.root_node.child_count == 3
    @test path5.diagnostics.root_node.candidate_summaries[2].child_parallel_counts == [6, 3, 6]
    @test path5.diagnostics.root_node.candidate_summaries[2].child_parallel_to_transverse_ratios[2] > 0.35
end

@testset "Experimental homonuclear chain nested dense export" begin
    basis3, path3, _check3 = _bond_aligned_homonuclear_chain_nested_qw_fixture(; natoms = 3, spacing = 1.2)
    basis4, path4, _check4 = _bond_aligned_homonuclear_chain_nested_qw_fixture(; natoms = 4, spacing = 1.2)
    basis5, path5, _check5 = _bond_aligned_homonuclear_chain_nested_qw_fixture(; natoms = 5, spacing = 1.2)

    for (basis, path, natoms) in ((basis3, path3, 3), (basis4, path4, 4), (basis5, path5, 5))
        payload_data = experimental_homonuclear_chain_nested_dense_payload(
            path;
            meta = (example = "test_experimental_chain_nested_dense_export",),
        )
        @test payload_data.bridge_meta["format"] == "experimental_homonuclear_chain_nested_dense_v1"
        @test payload_data.bridge_meta["experimental"]
        @test payload_data.bridge_meta["odd_chain_policy"] == "central_ternary_relaxed"
        @test payload_data.bridge_meta["fixed_dimension"] == size(path.fixed_block.overlap, 1)
        @test payload_data.bridge_meta["leaf_count"] == path.diagnostics.leaf_count
        @test payload_data.bridge_meta["root_accepted_candidate_index"] ==
            something(path.diagnostics.root_node.accepted_candidate_index, 0)
        @test payload_data.bridge_meta["residual_sector_empty"]
        @test size(payload_data.payload["S"]) == size(path.operators.overlap)
        @test size(payload_data.payload["H1"]) == size(path.operators.one_body_hamiltonian)
        @test size(payload_data.payload["Vee"]) == size(path.operators.interaction_matrix)
        @test size(payload_data.payload["basis_centers"]) == size(path.fixed_block.fixed_centers)
        @test size(payload_data.payload["nuclear_coordinates_xyz"]) == (natoms, 3)
        @test payload_data.payload["nuclear_charges"] == fill(1.0, natoms)
        @test payload_data.payload["orbital_labels"] ==
            String[orbital.label for orbital in orbitals(path.operators)]
        @test payload_data.payload["geometry_report_text"] isa String
        @test occursin("odd_chain_policy = central_ternary_relaxed", payload_data.payload["geometry_report_text"])
        @test occursin("[node root]", payload_data.payload["geometry_report_text"])
        @test payload_data.meta_values["manifest/contract/status"] == "experimental"
        @test payload_data.meta_values["manifest/source/odd_chain_policy"] == "central_ternary_relaxed"
        @test payload_data.meta_values["manifest/source/fixed_dimension"] == path.diagnostics.fixed_dimension
    end

    mktempdir() do dir
        export_path = joinpath(dir, "h5_chain_nested_dense.jld2")
        @test write_experimental_homonuclear_chain_nested_dense_jld2(
            export_path,
            path5;
            meta = (example = "test_experimental_chain_nested_dense_export",),
        ) == export_path
        jldopen(export_path, "r") do file
            @test String(file["bridge/format"]) == "experimental_homonuclear_chain_nested_dense_v1"
            @test Bool(file["bridge/experimental"])
            @test String(file["bridge/odd_chain_policy"]) == "central_ternary_relaxed"
            @test Bool(file["bridge/root_did_split"])
            @test Int(file["bridge/root_accepted_candidate_index"]) == 2
            @test Int(file["bridge/fixed_dimension"]) == size(path5.fixed_block.overlap, 1)
            @test Bool(file["bridge/residual_sector_empty"])
            @test size(file["S"]) == size(path5.operators.overlap)
            @test size(file["H1"]) == size(path5.operators.one_body_hamiltonian)
            @test size(file["Vee"]) == size(path5.operators.interaction_matrix)
            @test size(file["nuclear_coordinates_xyz"]) == (5, 3)
            @test String(file["meta/producer"]) ==
                "GaussletBases.write_experimental_homonuclear_chain_nested_dense_jld2"
            @test String(file["meta/manifest/contract/status"]) == "experimental"
            @test String(file["meta/manifest/source/odd_chain_policy"]) == "central_ternary_relaxed"
            @test String(file["meta/example"]) == "test_experimental_chain_nested_dense_export"
        end
    end
end

@testset "Atomic hybrid anchor comparison" begin
    if !_legacy_basisfile_available()
        @test true
    else
        (
            _source_basis,
            bundle,
            _shell1_face,
            _shell2_face,
            _shell3_face,
            _shell4_face,
            _shell1_complete,
            _shell2_complete,
            _shell3_complete,
            _shell4_complete,
            _interval1,
            _interval2,
            _interval3,
            _interval4,
            _core5,
            _shell_plus_core,
            _face_sequence,
            complete_sequence,
            fixed_shell_plus_core,
            _fixed_face_sequence,
            fixed_complete_sequence,
            legacy,
            baseline,
            shell_plus_core_ops,
            _face_sequence_ops,
            complete_sequence_ops,
            baseline_check,
            shell_plus_core_check,
            _face_sequence_check,
            complete_sequence_check,
        ) = _nested_qiu_white_complete_shell_sequence_fixture(count = 15)

        expansion = coulomb_gaussian_expansion(doacc = false)
        overlap_parent, one_body_parent, interaction_parent =
            _nested_parent_fixed_problem(bundle, expansion; Z = 2.0)
        parent_modes = eigen(Hermitian(one_body_parent), Hermitian(overlap_parent))
        parent_ground = parent_modes.vectors[:, 1]
        parent_ground_vee = _nested_vee_from_orbital(interaction_parent, parent_ground)
        projected_complete = _nested_fixed_projected_orbital(
            overlap_parent,
            fixed_complete_sequence,
            parent_ground,
        )
        projected_complete_vee = _nested_vee_from_orbital(
            GaussletBases._qwrg_fixed_block_interaction_matrix(fixed_complete_sequence, expansion),
            projected_complete,
        )
        complete_ground_capture, complete_ground_energy = _nested_projector_stats(
            overlap_parent,
            one_body_parent,
            fixed_complete_sequence,
            parent_ground,
        )
        complete_average4_capture =
            sum(
                _nested_projector_stats(
                    overlap_parent,
                    one_body_parent,
                    fixed_complete_sequence,
                    parent_modes.vectors[:, index],
                )[1] for index in 1:4
            ) / 4

        @test legacy isa LegacyAtomicGaussianSupplement
        @test legacy.lmax == 0
        @test complete_sequence.working_box == (2:14, 2:14, 2:14)
        @test baseline.gausslet_count == 15^3
        @test shell_plus_core_ops.gausslet_count == 1385
        @test complete_sequence_ops.gausslet_count < shell_plus_core_ops.gausslet_count
        @test norm(fixed_complete_sequence.overlap - I, Inf) < 1.0e-10
        @test complete_sequence_check.overlap_error < 1.0e-10
        @test all(isfinite, fixed_complete_sequence.weights)
        @test minimum(fixed_complete_sequence.weights) > 0.0
        @test maximum(fixed_complete_sequence.weights) < 10.0
        @test abs(complete_sequence_check.orbital_energy - baseline_check.orbital_energy) < 2.0e-4
        @test abs(complete_sequence_check.vee_expectation - baseline_check.vee_expectation) < 1.0e-4
        @test abs(complete_sequence_check.vee_expectation - shell_plus_core_check.vee_expectation) < 1.0e-4
        @test abs(projected_complete_vee - parent_ground_vee) < 1.5e-4
        @test complete_ground_capture > 0.99999
        @test complete_average4_capture > 0.999
        @test complete_ground_energy - parent_modes.values[1] < 1.0e-4
    end
end

@testset "Atomic hierarchical core-only refinement" begin
    if !_legacy_basisfile_available()
        @test true
    else
        for count in (17, 15)
            (
                _source_basis,
                bundle,
                core5,
                complete_sequence,
                fixed_complete_sequence,
                complete_sequence_ops,
                complete_sequence_check,
                refined_core,
                refined_sequence,
                fixed_refined_sequence,
                refined_sequence_ops,
                refined_sequence_check,
                legacy,
                _baseline,
                _baseline_check,
            ) = _nested_qiu_white_hierarchical_core_fixture(count = count)

            expansion = coulomb_gaussian_expansion(doacc = false)
            overlap_parent, one_body_parent, interaction_parent =
                _nested_parent_fixed_problem(bundle, expansion; Z = 2.0)
            parent_modes = eigen(Hermitian(one_body_parent), Hermitian(overlap_parent))
            parent_ground = parent_modes.vectors[:, 1]
            parent_ground_vee = _nested_vee_from_orbital(interaction_parent, parent_ground)
            projected_refined = _nested_fixed_projected_orbital(
                overlap_parent,
                fixed_refined_sequence,
                parent_ground,
            )
            projected_refined_vee = _nested_vee_from_orbital(
                GaussletBases._qwrg_fixed_block_interaction_matrix(fixed_refined_sequence, expansion),
                projected_refined,
            )
            refined_ground_capture, refined_ground_energy = _nested_projector_stats(
                overlap_parent,
                one_body_parent,
                fixed_refined_sequence,
                parent_ground,
            )

            inner_core = (first(core5) + 1):(last(core5) - 1)
            @test legacy isa LegacyAtomicGaussianSupplement
            @test legacy.lmax == 0
            @test refined_core isa GaussletBases._CartesianNestedShellSequence3D
            @test refined_sequence isa GaussletBases._CartesianNestedShellSequence3D
            @test refined_core.working_box == (core5, core5, core5)
            @test refined_sequence.working_box == complete_sequence.working_box
            @test length(refined_core.core_indices) == length(inner_core)^3
            @test size(refined_core.coefficient_matrix, 2) == 53
            @test refined_sequence_ops.gausslet_count == 445
            @test refined_sequence_ops.gausslet_count < complete_sequence_ops.gausslet_count
            @test norm(fixed_refined_sequence.overlap - I, Inf) < 1.0e-10
            @test refined_sequence_check.overlap_error < 1.0e-10
            @test all(isfinite, fixed_refined_sequence.weights)
            @test minimum(fixed_refined_sequence.weights) > 0.0
            @test refined_ground_capture > 0.999
            @test abs(refined_sequence_check.vee_expectation - complete_sequence_check.vee_expectation) > 5.0e-4
            @test abs(projected_refined_vee - parent_ground_vee) > 1.0e-3
            @test refined_ground_energy - parent_modes.values[1] > 1.0e-2
        end
    end
end

@testset "Mapped ordinary one-body backends" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    mild_basis = build_basis(MappedUniformBasisSpec(:G10;
        count = 5,
        mapping = fit_asinh_mapping_for_strength(s = 0.5, npoints = 5, xmax = 6.0),
        reference_spacing = 1.0,
    ))
    mild_reference = mapped_ordinary_one_body_operators(
        mild_basis;
        exponents = expansion.exponents[1:3],
        backend = :numerical_reference,
    )
    mild_analytic = mapped_ordinary_one_body_operators(
        mild_basis;
        exponents = expansion.exponents[1:3],
        backend = :pgdg_experimental,
    )
    mild_localized = mapped_ordinary_one_body_operators(
        mild_basis;
        exponents = expansion.exponents[1:3],
        backend = :pgdg_localized_experimental,
    )
    mild_oracle = GaussletBases._mapped_ordinary_localized_oracle_operators(
        mild_basis;
        exponents = expansion.exponents[1:3],
    )
    (_, _, overlap_reference_localized, kinetic_reference_localized, _) =
        _localized_numerical_reference_1d(mild_basis, expansion.exponents[1:3])
    raw_localized = mapped_pgdg_localized(GaussletBases.mapped_pgdg_logfit_prototype(mild_basis))
    raw_localized_kinetic = kinetic_matrix(raw_localized)

    @test mild_reference isa MappedOrdinaryOneBody1D
    @test mild_analytic isa MappedOrdinaryOneBody1D
    @test mild_localized isa MappedOrdinaryOneBody1D
    @test mild_oracle isa MappedOrdinaryOneBody1D
    @test mild_reference.backend == :numerical_reference
    @test mild_analytic.backend == :pgdg_experimental
    @test mild_localized.backend == :pgdg_localized_experimental
    @test mild_oracle.backend == :pgdg_localized_oracle
    @test occursin("experimental=true", sprint(show, mild_analytic))
    @test occursin("experimental=true", sprint(show, mild_localized))
    @test occursin("experimental=true", sprint(show, mild_oracle))
    @test !occursin("experimental=true", sprint(show, mild_reference))
    @test mild_reference.overlap ≈ transpose(mild_reference.overlap) atol = 1.0e-10 rtol = 1.0e-10
    @test mild_analytic.overlap ≈ transpose(mild_analytic.overlap) atol = 1.0e-10 rtol = 1.0e-10
    @test mild_localized.overlap ≈ transpose(mild_localized.overlap) atol = 1.0e-10 rtol = 1.0e-10
    @test mild_reference.kinetic ≈ transpose(mild_reference.kinetic) atol = 1.0e-10 rtol = 1.0e-10
    @test mild_analytic.kinetic ≈ transpose(mild_analytic.kinetic) atol = 1.0e-10 rtol = 1.0e-10
    @test mild_localized.kinetic ≈ transpose(mild_localized.kinetic) atol = 1.0e-10 rtol = 1.0e-10
    @test length(mild_reference.gaussian_factors) == 3
    @test length(mild_analytic.gaussian_factors) == 3
    @test length(mild_localized.gaussian_factors) == 3
    @test mild_basis isa MappedUniformBasis
    @test norm(mild_reference.overlap - mild_analytic.overlap, Inf) < 0.05
    @test norm(mild_reference.kinetic - mild_analytic.kinetic, Inf) < 0.05
    @test norm(mild_reference.gaussian_factors[1] - mild_analytic.gaussian_factors[1], Inf) < 0.05
    @test norm(mild_localized.overlap - I, Inf) < 1.0e-10
    @test norm(mild_localized.overlap - I, Inf) < norm(mild_analytic.overlap - I, Inf)
    @test norm(mild_localized.overlap - overlap_reference_localized, Inf) < 1.0e-10
    @test norm(mild_localized.kinetic - kinetic_reference_localized, Inf) <
          norm(raw_localized_kinetic - kinetic_reference_localized, Inf)
    @test norm(mild_oracle.kinetic - kinetic_reference_localized, Inf) <
          norm(mild_localized.kinetic - kinetic_reference_localized, Inf)
end

@testset "Public ordinary and nested backend contract" begin
    function _reference_only_backend_error(route_builder)
        err = try
            route_builder()
            nothing
        catch err
            err
        end
        @test err isa ArgumentError
        text = sprint(showerror, err)
        @test occursin("numerical-reference-only route", text)
        @test occursin("PGDG production-contract support is not yet implemented here", text)
        return text
    end

    function _check_direct_product_backend_pair(
        basis,
        nuclear_charges::AbstractVector{<:Real},
    )
        reference = ordinary_cartesian_qiu_white_operators(
            basis;
            nuclear_charges = nuclear_charges,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :numerical_reference,
        )
        localized = ordinary_cartesian_qiu_white_operators(
            basis;
            nuclear_charges = nuclear_charges,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :pgdg_localized_experimental,
        )
        reference_check = GaussletBases.ordinary_cartesian_1s2_check(reference)
        localized_check = GaussletBases.ordinary_cartesian_1s2_check(localized)

        @test reference.gausslet_backend == :numerical_reference
        @test localized.gausslet_backend == :pgdg_localized_experimental
        @test reference.gaussian_data === nothing
        @test localized.gaussian_data === nothing
        @test reference.residual_count == 0
        @test localized.residual_count == 0
        @test norm(reference.overlap - I, Inf) < 1.0e-8
        @test norm(localized.overlap - I, Inf) < 1.0e-8
        @test localized.one_body_hamiltonian ≈ transpose(localized.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
        @test localized.interaction_matrix ≈ transpose(localized.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
        @test norm(localized.overlap - reference.overlap, Inf) < 0.05
        @test norm(localized.one_body_hamiltonian - reference.one_body_hamiltonian, Inf) < 0.35
        @test norm(localized.interaction_matrix - reference.interaction_matrix, Inf) < 0.3
        @test abs(localized_check.orbital_energy - reference_check.orbital_energy) < 0.05
        @test abs(localized_check.vee_expectation - reference_check.vee_expectation) < 0.05
        return reference, localized, reference_check, localized_check
    end

    expansion = coulomb_gaussian_expansion(doacc = false)
    mapped_basis = build_basis(MappedUniformBasisSpec(:G10;
        count = 5,
        mapping = fit_asinh_mapping_for_strength(s = 0.5, npoints = 5, xmax = 6.0),
        reference_spacing = 1.0,
    ))
    diatomic_basis = bond_aligned_homonuclear_qw_basis(
        bond_length = 1.4,
        core_spacing = 0.5,
        xmax_parallel = 6.0,
        xmax_transverse = 4.0,
        bond_axis = :z,
    )
    chain_basis = bond_aligned_homonuclear_chain_qw_basis(
        natoms = 3,
        spacing = 1.2,
        core_spacing = 0.5,
        xmax_parallel = 2.0,
        xmax_transverse = 2.0,
        chain_axis = :z,
    )
    square_basis = axis_aligned_homonuclear_square_lattice_qw_basis(
        n = 2,
        spacing = 1.4,
        core_spacing = 0.5,
        xmax_in_plane = 2.0,
        xmax_transverse = 2.0,
    )

    @test mapped_ordinary_one_body_operators(
        mapped_basis;
        exponents = expansion.exponents[1:3],
        backend = :pgdg_experimental,
    ).backend == :pgdg_experimental
    @test mapped_ordinary_one_body_operators(
        mapped_basis;
        exponents = expansion.exponents[1:3],
        backend = :pgdg_localized_experimental,
    ).backend == :pgdg_localized_experimental

    diatomic_reference, diatomic_localized, diatomic_reference_check, diatomic_localized_check =
        _check_direct_product_backend_pair(diatomic_basis, [1.0, 1.0])
    chain_reference, chain_localized, chain_reference_check, chain_localized_check =
        _check_direct_product_backend_pair(
            chain_basis,
            fill(1.0, length(chain_basis.nuclei)),
        )
    square_reference, square_localized, square_reference_check, square_localized_check =
        _check_direct_product_backend_pair(
            square_basis,
            fill(1.0, length(square_basis.nuclei)),
        )

    @test diatomic_localized.gausslet_count == diatomic_reference.gausslet_count
    @test chain_localized.gausslet_count == chain_reference.gausslet_count
    @test square_localized.gausslet_count == square_reference.gausslet_count
    @test isfinite(diatomic_localized_check.orbital_energy)
    @test isfinite(chain_localized_check.orbital_energy)
    @test isfinite(square_localized_check.orbital_energy)
    @test isfinite(diatomic_localized_check.vee_expectation)
    @test isfinite(chain_localized_check.vee_expectation)
    @test isfinite(square_localized_check.vee_expectation)

    direct_product_error_text = _reference_only_backend_error(() ->
        ordinary_cartesian_qiu_white_operators(
            diatomic_basis;
            nuclear_charges = [1.0, 1.0],
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :pgdg_experimental,
        )
    )
    @test occursin("bond-aligned ordinary_cartesian_qiu_white_operators", direct_product_error_text)
    @test occursin(":pgdg_localized_experimental", direct_product_error_text)

    diatomic_nested_source_text = _reference_only_backend_error(() ->
        bond_aligned_diatomic_nested_fixed_source(
            diatomic_basis;
            expansion = expansion,
            gausslet_backend = :pgdg_localized_experimental,
        )
    )
    @test occursin("bond-aligned diatomic nested fixed source", diatomic_nested_source_text)

    diatomic_nested_fixed_text = _reference_only_backend_error(() ->
        bond_aligned_diatomic_nested_fixed_block(
            diatomic_basis;
            expansion = expansion,
            gausslet_backend = :pgdg_experimental,
        )
    )
    @test occursin("bond-aligned diatomic nested fixed source", diatomic_nested_fixed_text)

    chain_text = _reference_only_backend_error(() ->
        experimental_bond_aligned_homonuclear_chain_nested_qw_operators(
            chain_basis;
            nuclear_charges = fill(1.0, length(chain_basis.nuclei)),
            gausslet_backend = :pgdg_experimental,
        )
    )
    @test occursin("bond-aligned homonuclear chain nested fixed source", chain_text)

    square_text = _reference_only_backend_error(() ->
        experimental_axis_aligned_homonuclear_square_lattice_nested_qw_operators(
            square_basis;
            nuclear_charges = fill(1.0, length(square_basis.nuclei)),
            gausslet_backend = :pgdg_localized_experimental,
        )
    )
    @test occursin("axis-aligned homonuclear square-lattice nested fixed source", square_text)

    if !_legacy_basisfile_available()
        @test true
    else
        supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 0)
        supplement_text = _reference_only_backend_error(() ->
            ordinary_cartesian_qiu_white_operators(
                mapped_basis,
                supplement;
                expansion = expansion,
                Z = 2.0,
                interaction_treatment = :ggt_nearest,
                gausslet_backend = :pgdg_localized_experimental,
            )
        )
        @test occursin("ordinary_cartesian_qiu_white_operators", supplement_text)
    end
end

@testset "Ordinary Cartesian IDA operators" begin
    mild_basis, mild_expansion, mild_analytic = _quick_ordinary_cartesian_ida_fixture(
        backend = :pgdg_experimental,
        mapped = true,
        s = 0.5,
    )
    (_, _, identity_analytic) = _quick_ordinary_cartesian_ida_fixture(
        backend = :pgdg_experimental,
        mapped = false,
    )

    function reconstruct_interaction(expansion::CoulombGaussianExpansion, factors::AbstractVector{<:AbstractMatrix})
        interaction = zeros(Float64, size(first(factors), 1)^3, size(first(factors), 1)^3)
        for term in eachindex(expansion.coefficients)
            factor = factors[term]
            interaction .+= expansion.coefficients[term] .* kron(factor, kron(factor, factor))
        end
        return interaction
    end

    @test mild_basis isa MappedUniformBasis
    @test mild_analytic isa OrdinaryCartesianIDAOperators
    @test identity_analytic isa OrdinaryCartesianIDAOperators
    @test mild_analytic.backend == :pgdg_experimental
    @test identity_analytic.backend == :pgdg_experimental
    @test occursin("experimental=true", sprint(show, mild_analytic))
    @test occursin("experimental=true", sprint(show, identity_analytic))
    @test length(orbitals(mild_analytic)) == length(mild_basis)^3
    @test orbitals(mild_analytic)[1].ix == 1
    @test orbitals(mild_analytic)[1].iy == 1
    @test orbitals(mild_analytic)[1].iz == 1
    @test orbitals(mild_analytic)[end].ix == length(mild_basis)
    @test orbitals(mild_analytic)[end].iy == length(mild_basis)
    @test orbitals(mild_analytic)[end].iz == length(mild_basis)
    @test mild_analytic.overlap_3d ≈ transpose(mild_analytic.overlap_3d) atol = 1.0e-10 rtol = 1.0e-10
    @test mild_analytic.one_body_hamiltonian ≈ transpose(mild_analytic.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test mild_analytic.interaction_matrix ≈ transpose(mild_analytic.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test all(
        factor -> isapprox(factor, transpose(factor); atol = 1.0e-10, rtol = 1.0e-10),
        mild_analytic.pair_factors_1d,
    )
    @test mild_analytic.interaction_matrix ≈ reconstruct_interaction(mild_expansion, mild_analytic.pair_factors_1d) atol = 1.0e-10 rtol = 1.0e-10
    @test minimum(diag(mild_analytic.interaction_matrix)) > 0.0
    @test minimum(diag(identity_analytic.interaction_matrix)) > 0.0
    @test size(identity_analytic.one_body_hamiltonian) == (length(identity_analytic.orbital_data), length(identity_analytic.orbital_data))
    @test opnorm(mild_analytic.interaction_matrix, Inf) > 0.0
end

@testset "Ordinary Cartesian 1s^2 Vee check" begin
    basis, operators, orbital_energy, orbital, vee = _quick_ordinary_cartesian_1s2_vee_fixture()
    reference_value = 1.25

    @test basis isa MappedUniformBasis
    @test operators isa OrdinaryCartesianIDAOperators
    @test operators.backend == :numerical_reference
    @test norm(operators.overlap_3d - I, Inf) < 1.0e-10
    @test operators.interaction_matrix ≈ transpose(operators.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test isfinite(orbital_energy)
    @test orbital_energy < 0.0
    @test isfinite(vee)
    @test vee > 0.0
    @test abs(sum(abs2, orbital) - 1.0) < 1.0e-10
    @test abs(vee - reference_value) < 0.02
    @test ordinary_cartesian_vee_expectation(operators, 2.0 .* orbital) ≈ vee atol = 1.0e-12 rtol = 0.0
end

@testset "Legacy 1D hybrid ordinary Cartesian 1s^2 Vee check" begin
    (
        source_basis,
        core_gaussians,
        pure_operators,
        hybrid_basis,
        hybrid_operators,
        pure_check,
        hybrid_check,
    ) = _quick_hybrid_cartesian_1s2_vee_fixture()
    reference_value = 1.25

    @test source_basis isa MappedUniformBasis
    @test hybrid_basis isa HybridMappedOrdinaryBasis1D
    @test hybrid_operators isa OrdinaryCartesianIDAOperators
    @test hybrid_operators.backend == :pgdg_localized_experimental
    @test length(core_gaussians) == 2
    @test length(hybrid_basis) == length(source_basis) + length(core_gaussians)
    @test length(orbitals(hybrid_operators)) == length(hybrid_basis)^3
    @test norm(hybrid_operators.overlap_3d - I, Inf) < 1.0e-10
    @test hybrid_operators.one_body_hamiltonian ≈ transpose(hybrid_operators.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test hybrid_operators.interaction_matrix ≈ transpose(hybrid_operators.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test isfinite(pure_check.orbital_energy)
    @test isfinite(hybrid_check.orbital_energy)
    @test isfinite(pure_check.vee_expectation)
    @test isfinite(hybrid_check.vee_expectation)
    @test pure_check.vee_expectation > 0.0
    @test hybrid_check.vee_expectation > 0.0
    @test hybrid_check.vee_expectation > pure_check.vee_expectation
    @test abs(hybrid_check.vee_expectation - reference_value) <
          abs(pure_check.vee_expectation - reference_value)
    @test abs(hybrid_check.orbital_energy + 2.0) < abs(pure_check.orbital_energy + 2.0)
@test ordinary_cartesian_vee_expectation(hybrid_operators, 3.0 .* hybrid_check.orbital) ≈
          hybrid_check.vee_expectation atol = 1.0e-12 rtol = 0.0
end

@testset "Legacy 1D hybrid residual Gaussian interaction treatment" begin
    (
        source_basis,
        hybrid_basis,
        pure_operators,
        combined_operators,
        residual_operators,
        pure_check,
        combined_check,
        residual_check,
    ) = _friendly_hybrid_residual_vee_fixture(11, 0.6)
    reference_value = 1.25

    @test source_basis isa MappedUniformBasis
    @test hybrid_basis isa HybridMappedOrdinaryBasis1D
    @test combined_operators.interaction_treatment == :combined_basis
    @test residual_operators.interaction_treatment == :residual_gaussian_nearest
    @test norm(combined_operators.overlap_3d - residual_operators.overlap_3d, Inf) < 1.0e-12
    @test norm(combined_operators.one_body_hamiltonian - residual_operators.one_body_hamiltonian, Inf) < 1.0e-12
    @test combined_operators.interaction_matrix ≈ transpose(combined_operators.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test residual_operators.interaction_matrix ≈ transpose(residual_operators.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test pure_check.vee_expectation > 0.0
    @test combined_check.vee_expectation > 0.0
    @test residual_check.vee_expectation > 0.0
    @test abs(combined_check.vee_expectation - residual_check.vee_expectation) > 1.0e-2
    @test abs(combined_check.vee_expectation - reference_value) <
          abs(residual_check.vee_expectation - reference_value)
    @test abs(combined_check.orbital_energy - residual_check.orbital_energy) < 1.0e-12
    @test abs(pure_check.vee_expectation - reference_value) < 0.02
    @test abs(combined_check.orbital_energy + 2.0) < abs(pure_check.orbital_energy + 2.0)
end

@testset "Legacy atomic Gaussian supplement" begin
    if !_legacy_basisfile_available()
        @test true
    else
        vtz0 = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 0)
        vtz1 = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 1)
        vtz = legacy_s_gaussian_data("He", "cc-pVTZ")
        vqz = legacy_s_gaussian_data("He", "cc-pVQZ")
        vtz_uncontracted = legacy_s_gaussian_data("He", "cc-pVTZ"; uncontracted = true)

        @test vtz0 isa LegacyAtomicGaussianSupplement
        @test vtz isa LegacySGaussianData
        @test vqz isa LegacySGaussianData
        @test vtz0.lmax == 0
        @test vtz1.lmax == 1
        @test length(vtz0.shells) == 3
        @test any(shell -> shell.l == 1, vtz1.shells)
        @test length(vtz1.shells) > length(vtz0.shells)
        @test vtz.primitive_exponents == vtz0.primitive_exponents
        @test vtz.primitive_widths == vtz0.primitive_widths
        @test vtz.contraction_matrix ≈ vtz0.contraction_matrix atol = 0.0 rtol = 0.0
        @test vtz.widths ≈ vtz0.widths atol = 0.0 rtol = 0.0
        @test vtz1.primitive_exponents == vtz0.primitive_exponents
        @test vtz1.primitive_widths == vtz0.primitive_widths
        @test vtz1.contraction_matrix ≈ vtz0.contraction_matrix atol = 0.0 rtol = 0.0
        @test vtz1.widths ≈ vtz0.widths atol = 0.0 rtol = 0.0
        @test !vtz.uncontracted
        @test !vqz.uncontracted
        @test vtz.max_width === nothing
        @test vqz.max_width === nothing
        @test vtz.primitive_exponents == [234.0, 35.16, 7.989, 2.212, 0.6669, 0.2089]
        @test vqz.primitive_exponents == [528.5, 79.31, 18.05, 5.085, 1.609, 0.5363, 0.1833]
        @test length(vtz.primitive_gaussians) == 6
        @test length(vqz.primitive_gaussians) == 7
        @test size(vtz.contraction_matrix) == (6, 3)
        @test size(vqz.contraction_matrix) == (7, 4)
        @test length(vtz.gaussians) == 3
        @test length(vqz.gaussians) == 4
        @test all(isfinite, vtz.widths)
        @test all(isfinite, vqz.widths)
        @test all(>(0.0), vtz.widths)
        @test all(>(0.0), vqz.widths)
        @test all(abs(gaussian.center_value) < 1.0e-12 for gaussian in vtz.gaussians)
        @test all(abs(gaussian.center_value) < 1.0e-12 for gaussian in vqz.gaussians)
        @test vtz.primitive_widths[1] ≈ inv(sqrt(2 * 234.0)) atol = 1.0e-12 rtol = 0.0
        @test vqz.primitive_widths[end] ≈ inv(sqrt(2 * 0.1833)) atol = 1.0e-12 rtol = 0.0
        @test vtz_uncontracted.uncontracted
        @test size(vtz_uncontracted.contraction_matrix) == (6, 6)
        @test vtz_uncontracted.contraction_matrix ≈ Matrix{Float64}(I, 6, 6) atol = 0.0 rtol = 0.0
        @test length(vtz_uncontracted.gaussians) == 6
    end
end

@testset "Legacy bond-aligned diatomic Gaussian supplement width cutoff" begin
    if !_legacy_basisfile_available()
        @test true
    else
        nuclei = [(-0.7, 0.0, 0.0), (0.7, 0.0, 0.0)]
        full = legacy_bond_aligned_diatomic_gaussian_supplement(
            "H",
            "cc-pVTZ",
            nuclei;
            lmax = 1,
        )
        trimmed = legacy_bond_aligned_diatomic_gaussian_supplement(
            "H",
            "cc-pVTZ",
            nuclei;
            lmax = 1,
            max_width = 1.0,
        )
        hetero_full = legacy_bond_aligned_heteronuclear_gaussian_supplement(
            "He",
            "cc-pVTZ",
            "H",
            "cc-pVTZ",
            nuclei;
            lmax = 1,
        )
        hetero_trimmed = legacy_bond_aligned_heteronuclear_gaussian_supplement(
            "He",
            "cc-pVTZ",
            "H",
            "cc-pVTZ",
            nuclei;
            lmax = 1,
            max_width = 1.0,
        )
        full_cart = GaussletBases._bond_aligned_diatomic_cartesian_shell_supplement_3d(full)
        trimmed_cart = GaussletBases._bond_aligned_diatomic_cartesian_shell_supplement_3d(trimmed)
        hetero_full_cart =
            GaussletBases._bond_aligned_diatomic_cartesian_shell_supplement_3d(hetero_full)
        hetero_trimmed_cart =
            GaussletBases._bond_aligned_diatomic_cartesian_shell_supplement_3d(hetero_trimmed)

        @test full.max_width === nothing
        @test trimmed.max_width == 1.0
        @test hetero_full.max_width === nothing
        @test hetero_trimmed.max_width == 1.0
        @test full.atomic_source.max_width === nothing
        @test trimmed.atomic_source.max_width === nothing
        @test length(full_cart.orbitals) == 18
        @test length(trimmed_cart.orbitals) == 8
        @test length(hetero_trimmed_cart.orbitals) < length(hetero_full_cart.orbitals)
    end
end

@testset "Vendored legacy BasisSets lookup and overrides" begin
    vendored = GaussletBases._vendored_legacy_basisfile_path()
    @test isfile(vendored)
    @test occursin(joinpath("data", "legacy", "BasisSets"), vendored)

    withenv("GAUSSLETBASES_BASISSETS_PATH" => nothing) do
        @test GaussletBases._legacy_basisfile_path() == vendored
        he_supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 0)
        h_supplement = legacy_atomic_gaussian_supplement("H", "cc-pVTZ"; lmax = 1)
        @test he_supplement.basisfile == vendored
        @test h_supplement.basisfile == vendored
    end

    mktemp() do path, io
        write(
            io,
            "#BASIS SET: He repo-test\n" *
            "He    S\n" *
            "      1.0000000              1.0000000\n" *
            "END\n",
        )
        close(io)

        withenv("GAUSSLETBASES_BASISSETS_PATH" => path) do
            @test GaussletBases._legacy_basisfile_path() == path
            supplement = legacy_atomic_gaussian_supplement("He", "repo-test"; lmax = 0)
            @test supplement.basisfile == path
            @test supplement.primitive_exponents == [1.0]
        end

        supplement = legacy_atomic_gaussian_supplement("He", "repo-test"; lmax = 0, basisfile = path)
        @test supplement.basisfile == path
        @test supplement.primitive_exponents == [1.0]
    end
end

@testset "Atomic lmax=0 supplement uses the explicit 3D shell route in QW consumers" begin
    if !_legacy_basisfile_available()
        @test true
    else
        source_basis, _legacy_old, baseline_ops, baseline_check = _qiu_white_full_nearest_fixture()
        supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 0)
        @test !GaussletBases._legacy_atomic_has_nonseparable_shells(supplement)
        ordinary_ops = ordinary_cartesian_qiu_white_operators(
            source_basis,
            supplement;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        explicit_ordinary = GaussletBases._ordinary_cartesian_qiu_white_operators_atomic_shell_3d(
            source_basis,
            supplement;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
            gausslet_backend = :numerical_reference,
            timing = false,
        )
        ordinary_check = GaussletBases.ordinary_cartesian_1s2_check(ordinary_ops)
        explicit_ordinary_check = GaussletBases.ordinary_cartesian_1s2_check(explicit_ordinary)

        (
            _source_basis_nested,
            _bundle_nested,
            _shell_nested,
            _fixed_block_nested,
            _shell_plus_core_nested,
            fixed_block_shell_plus_core,
            _legacy_nested,
            _baseline_nested,
            _nested_shell_only,
            nested_shell_plus_core,
            _baseline_nested_check,
            _nested_shell_only_check,
            nested_shell_plus_core_check,
        ) = _nested_qiu_white_nearest_fixture()
        nested_ops = ordinary_cartesian_qiu_white_operators(
            fixed_block_shell_plus_core,
            supplement;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        explicit_nested = GaussletBases._ordinary_cartesian_qiu_white_operators_nested_atomic_shell_3d(
            fixed_block_shell_plus_core,
            supplement;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
            gausslet_backend = :numerical_reference,
            timing = false,
        )
        nested_check = GaussletBases.ordinary_cartesian_1s2_check(nested_ops)
        explicit_nested_check = GaussletBases.ordinary_cartesian_1s2_check(explicit_nested)

        @test ordinary_ops.gaussian_data isa LegacyAtomicGaussianSupplement
        @test nested_ops.gaussian_data isa LegacyAtomicGaussianSupplement
        @test ordinary_ops.residual_count > 0
        @test nested_ops.residual_count > 0
        @test ordinary_check.overlap_error < 1.0e-8
        @test nested_check.overlap_error < 1.0e-8
        @test ordinary_check.orbital_energy ≈ explicit_ordinary_check.orbital_energy atol = 1.0e-12 rtol = 1.0e-12
        @test ordinary_check.vee_expectation ≈ explicit_ordinary_check.vee_expectation atol = 1.0e-12 rtol = 1.0e-12
        @test nested_check.orbital_energy ≈ explicit_nested_check.orbital_energy atol = 1.0e-12 rtol = 1.0e-12
        @test nested_check.vee_expectation ≈ explicit_nested_check.vee_expectation atol = 1.0e-12 rtol = 1.0e-12
        @test norm(ordinary_ops.overlap - explicit_ordinary.overlap, Inf) < 1.0e-12
        @test norm(ordinary_ops.one_body_hamiltonian - explicit_ordinary.one_body_hamiltonian, Inf) < 1.0e-12
        @test norm(nested_ops.overlap - explicit_nested.overlap, Inf) < 1.0e-12
        @test norm(nested_ops.one_body_hamiltonian - explicit_nested.one_body_hamiltonian, Inf) < 1.0e-12
        @test ordinary_check.orbital_energy ≈ baseline_check.orbital_energy atol = 1.0e-12 rtol = 1.0e-12
        @test nested_check.orbital_energy ≈ nested_shell_plus_core_check.orbital_energy atol = 1.0e-12 rtol = 1.0e-12
    end
end

@testset "Active atomic lmax=1 supplement is explicit and physical in QW routes" begin
    if !_legacy_basisfile_available()
        @test true
    else
        supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 1)
        supplement3d = GaussletBases._atomic_cartesian_shell_supplement_3d(supplement)

        @test GaussletBases._legacy_atomic_has_nonseparable_shells(supplement)
        @test any(orbital -> orbital.label == "px1", supplement3d.orbitals)
        @test any(orbital -> orbital.label == "py1", supplement3d.orbitals)
        @test any(orbital -> orbital.label == "pz1", supplement3d.orbitals)

        source_basis_hybrid, _legacy_old, _pure_check, _toy_check, _legacy_check =
            _legacy_he_s_hybrid_fixture("cc-pVTZ")
        hybrid_err = try
            hybrid_mapped_ordinary_basis(
                source_basis_hybrid;
                core_gaussians = supplement,
                backend = :pgdg_localized_experimental,
            )
            nothing
        catch err
            err
        end
        @test hybrid_err isa ArgumentError
        @test occursin("l > 0", sprint(showerror, hybrid_err))
        @test occursin("explicit 3D", sprint(showerror, hybrid_err))

        source_basis_qw, _legacy_qw, ordinary_l0, ordinary_l0_check = _qiu_white_full_nearest_fixture()
        ordinary_l1 = ordinary_cartesian_qiu_white_operators(
            source_basis_qw,
            supplement;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        ordinary_l1_check = GaussletBases.ordinary_cartesian_1s2_check(ordinary_l1)
        @test ordinary_l1.gaussian_data isa LegacyAtomicGaussianSupplement
        @test ordinary_l1.residual_count > 0
        @test ordinary_l1_check.overlap_error < 1.0e-8
        @test isfinite(ordinary_l1_check.orbital_energy)
        @test isfinite(ordinary_l1_check.vee_expectation)
        @test ordinary_l1_check.vee_expectation > 0.0
        @test abs(ordinary_l1_check.orbital_energy - ordinary_l0_check.orbital_energy) > 1.0e-6
        @test abs(ordinary_l1_check.vee_expectation - ordinary_l0_check.vee_expectation) > 1.0e-4
        @test size(ordinary_l1.one_body_hamiltonian) != size(ordinary_l0.one_body_hamiltonian) ||
              norm(ordinary_l1.one_body_hamiltonian - ordinary_l0.one_body_hamiltonian, Inf) > 1.0e-6

        (
            _source_basis_nested,
            _bundle_nested,
            _shell_nested,
            _fixed_block_nested,
            _shell_plus_core_nested,
            fixed_block_shell_plus_core,
            _legacy_nested,
            _baseline_nested,
            _nested_shell_only,
            nested_l0,
            _baseline_nested_check,
            _nested_shell_only_check,
            nested_l0_check,
        ) = _nested_qiu_white_nearest_fixture()
        nested_l1 = ordinary_cartesian_qiu_white_operators(
            fixed_block_shell_plus_core,
            supplement;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        nested_l1_check = GaussletBases.ordinary_cartesian_1s2_check(nested_l1)
        @test nested_l1.gaussian_data isa LegacyAtomicGaussianSupplement
        @test nested_l1.residual_count > 0
        @test nested_l1_check.overlap_error < 1.0e-8
        @test isfinite(nested_l1_check.orbital_energy)
        @test isfinite(nested_l1_check.vee_expectation)
        @test nested_l1_check.vee_expectation > 0.0
        @test abs(nested_l1_check.orbital_energy - nested_l0_check.orbital_energy) > 1.0e-6
        @test abs(nested_l1_check.vee_expectation - nested_l0_check.vee_expectation) > 1.0e-4
        @test size(nested_l1.one_body_hamiltonian) != size(nested_l0.one_body_hamiltonian) ||
              norm(nested_l1.one_body_hamiltonian - nested_l0.one_body_hamiltonian, Inf) > 1.0e-6
    end
end

@testset "Ne cc-pV6Z atomic shell referee for lmax=0 and lmax=1" begin
    mktemp() do path, io
        write(io, _ne_repo_v6z_sp_basis_text())
        close(io)

        basis = build_basis(MappedUniformBasisSpec(:G10;
            count = 9,
            mapping = fit_asinh_mapping_for_strength(s = 0.8, npoints = 9, xmax = 6.0),
            reference_spacing = 1.0,
        ))
        expansion = coulomb_gaussian_expansion(doacc = false)
        bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
            basis,
            exponents = expansion.exponents,
            center = 0.0,
            backend = :numerical_reference,
        )

        for lmax in (0, 1)
            supplement = legacy_atomic_gaussian_supplement("Ne", "repo-v6z-sp"; lmax = lmax, basisfile = path)
            supplement3d = GaussletBases._atomic_cartesian_shell_supplement_3d(supplement)
            public_ops = ordinary_cartesian_qiu_white_operators(
                basis,
                supplement;
                expansion = expansion,
                Z = 10.0,
                interaction_treatment = :ggt_nearest,
            )
            explicit_ops = GaussletBases._ordinary_cartesian_qiu_white_operators_atomic_shell_3d(
                basis,
                supplement;
                expansion = expansion,
                Z = 10.0,
                interaction_treatment = :ggt_nearest,
                gausslet_backend = :numerical_reference,
                timing = false,
            )
            supplement3d = GaussletBases._atomic_cartesian_shell_supplement_3d(supplement)
            blocks = GaussletBases._qwrg_atomic_cartesian_blocks_3d(bundle, supplement3d, expansion)
            shell_h1 = GaussletBases._qwrg_atomic_cartesian_one_body_aa(blocks, expansion; Z = 10.0)
            shell_vals = eigen(Hermitian(shell_h1), Hermitian(blocks.overlap_aa)).values

            @test norm(public_ops.overlap - explicit_ops.overlap, Inf) < 1.0e-12
            @test norm(public_ops.one_body_hamiltonian - explicit_ops.one_body_hamiltonian, Inf) < 1.0e-12
            @test abs(shell_vals[1] + 50.0) < 2.0e-3
            @test abs(shell_vals[2] + 12.5) < 2.0e-3

            if lmax == 0
                @test abs(shell_vals[3] + (50.0 / 9.0)) < 2.0e-3
            else
                @test length(supplement3d.orbitals) == 25
                @test maximum(abs.(shell_vals[3:5] .+ 12.5)) < 2.0e-3
                @test maximum(abs.(shell_vals[3:5] .- shell_vals[3])) < 1.0e-10
                @test abs(shell_vals[6] + (50.0 / 9.0)) < 2.0e-3
            end
        end
    end
end

@testset "Atomic lmax=2 supplement is explicit in QW routes but not yet molecular" begin
    mktemp() do path, io
        write(
            io,
            "#BASIS SET: He repo-spd\n" *
            "He    S\n" *
            "      1.0000000              1.0000000\n" *
            "He    P\n" *
            "      0.8000000              1.0000000\n" *
            "He    D\n" *
            "      0.6000000              1.0000000\n" *
            "END\n",
        )
        close(io)

        supplement = legacy_atomic_gaussian_supplement("He", "repo-spd"; lmax = 2, basisfile = path)
        supplement3d = GaussletBases._atomic_cartesian_shell_supplement_3d(supplement)

        @test any(orbital -> orbital.label == "dxx1", supplement3d.orbitals)
        @test any(orbital -> orbital.label == "dyy1", supplement3d.orbitals)
        @test any(orbital -> orbital.label == "dzz1", supplement3d.orbitals)
        @test any(orbital -> orbital.label == "dxy1", supplement3d.orbitals)
        @test any(orbital -> orbital.label == "dxz1", supplement3d.orbitals)
        @test any(orbital -> orbital.label == "dyz1", supplement3d.orbitals)

        source_basis = build_basis(MappedUniformBasisSpec(:G10;
            count = 9,
            mapping = fit_asinh_mapping_for_strength(s = 0.8, npoints = 9, xmax = 6.0),
            reference_spacing = 1.0,
        ))
        ordinary_ops = ordinary_cartesian_qiu_white_operators(
            source_basis,
            supplement;
            expansion = coulomb_gaussian_expansion(doacc = false),
            Z = 2.0,
            interaction_treatment = :ggt_nearest,
        )
        ordinary_check = GaussletBases.ordinary_cartesian_1s2_check(ordinary_ops)
        @test ordinary_ops.gaussian_data isa LegacyAtomicGaussianSupplement
        @test ordinary_ops.residual_count > 0
        @test ordinary_check.overlap_error < 1.0e-8
        @test isfinite(ordinary_check.orbital_energy)
        @test isfinite(ordinary_check.vee_expectation)
        @test ordinary_check.vee_expectation > 0.0

        diatomic = legacy_bond_aligned_diatomic_gaussian_supplement(
            "He",
            "repo-spd",
            [(-0.7, 0.0, 0.0), (0.7, 0.0, 0.0)];
            lmax = 2,
            basisfile = path,
        )
        diatomic_err = try
            GaussletBases._bond_aligned_diatomic_cartesian_shell_supplement_3d(diatomic)
            nothing
        catch err
            err
        end
        @test diatomic_err isa ArgumentError
        @test occursin("lmax <= 1", sprint(showerror, diatomic_err))
    end
end

@testset "Legacy He s hybrid supplement check" begin
    if !_legacy_basisfile_available()
        @test true
    else
        (
            source_basis,
            legacy,
            pure_check,
            toy_check,
            legacy_check,
        ) = _legacy_he_s_hybrid_fixture("cc-pVTZ")

        @test source_basis isa MappedUniformBasis
        @test legacy isa LegacySGaussianData
        @test isfinite(pure_check.orbital_energy)
        @test isfinite(toy_check.orbital_energy)
        @test isfinite(legacy_check.orbital_energy)
        @test isfinite(pure_check.vee_expectation)
        @test isfinite(toy_check.vee_expectation)
        @test isfinite(legacy_check.vee_expectation)
        @test legacy_check.vee_expectation > 0.0
        @test abs(legacy_check.orbital_energy + 2.0) < abs(toy_check.orbital_energy + 2.0)
        @test abs(legacy_check.orbital_energy + 2.0) < abs(pure_check.orbital_energy + 2.0)
    end
end

@testset "Legacy He s MWG residual interaction" begin
    if !_legacy_basisfile_available()
        @test true
    else
        (
            source_basis,
            legacy,
            hybrid_basis,
            pure_operators,
            combined_operators,
            nearest_operators,
            mwg_operators,
            pure_check,
            combined_check,
            nearest_check,
            mwg_check,
            mwg_data,
        ) = _legacy_he_s_mwg_fixture("cc-pVTZ")

        @test source_basis isa MappedUniformBasis
        @test legacy isa LegacySGaussianData
        @test hybrid_basis isa HybridMappedOrdinaryBasis1D
        @test combined_operators.interaction_treatment == :combined_basis
        @test nearest_operators.interaction_treatment == :residual_gaussian_nearest
        @test mwg_operators.interaction_treatment == :residual_gaussian_mwg
        @test norm(combined_operators.overlap_3d - nearest_operators.overlap_3d, Inf) < 1.0e-12
        @test norm(combined_operators.overlap_3d - mwg_operators.overlap_3d, Inf) < 1.0e-12
        @test norm(combined_operators.one_body_hamiltonian - nearest_operators.one_body_hamiltonian, Inf) < 1.0e-12
        @test norm(combined_operators.one_body_hamiltonian - mwg_operators.one_body_hamiltonian, Inf) < 1.0e-12
        @test mwg_operators.interaction_matrix ≈ transpose(mwg_operators.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
        @test all(isfinite, mwg_data.residual_centers)
        @test all(isfinite, mwg_data.residual_widths)
        @test all(>(0.0), mwg_data.residual_widths)
        @test maximum(abs.(mwg_data.residual_centers)) < 1.0e-8
        @test pure_check.vee_expectation > 0.0
        @test combined_check.vee_expectation > 0.0
        @test nearest_check.vee_expectation > 0.0
        @test mwg_check.vee_expectation > 0.0
        @test abs(mwg_check.vee_expectation - nearest_check.vee_expectation) > 1.0e-3
        @test abs(mwg_check.orbital_energy - combined_check.orbital_energy) < 1.0e-12
    end
end

@testset "Qiu-White residual Gaussian reference path" begin
    if !_legacy_basisfile_available() || !_RUN_SLOW_TESTS
        @test true
    else
        (
            source_basis,
            legacy,
            surrogate_mwg,
            qiu_nearest,
            qiu_mwg,
            surrogate_check,
            nearest_check,
            mwg_check,
        ) = _qiu_white_reference_fixture()

        @test source_basis isa MappedUniformBasis
        @test legacy isa LegacySGaussianData
        @test surrogate_mwg isa OrdinaryCartesianIDAOperators
        @test qiu_nearest isa QiuWhiteResidualGaussianOperators
        @test qiu_mwg isa QiuWhiteResidualGaussianOperators
        @test qiu_nearest.interaction_treatment == :ggt_nearest
        @test qiu_mwg.interaction_treatment == :mwg
        @test qiu_mwg.gausslet_backend == :numerical_reference
        @test qiu_nearest.residual_count > 0
        @test qiu_mwg.gausslet_count == length(source_basis)^3
        @test qiu_mwg.residual_count > 0
        @test size(qiu_mwg.raw_to_final, 1) == qiu_mwg.gausslet_count + length(legacy.gaussians)
        @test size(qiu_mwg.raw_to_final, 2) == qiu_mwg.gausslet_count + qiu_mwg.residual_count
        @test norm(qiu_nearest.overlap - I, Inf) < 1.0e-8
        @test norm(qiu_mwg.overlap - I, Inf) < 1.0e-8
        @test qiu_nearest.one_body_hamiltonian ≈ transpose(qiu_nearest.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
        @test qiu_mwg.one_body_hamiltonian ≈ transpose(qiu_mwg.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
        @test qiu_nearest.interaction_matrix ≈ transpose(qiu_nearest.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
        @test qiu_mwg.interaction_matrix ≈ transpose(qiu_mwg.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
        @test size(qiu_mwg.interaction_matrix) == (qiu_mwg.gausslet_count + qiu_mwg.residual_count, qiu_mwg.gausslet_count + qiu_mwg.residual_count)
        @test all(isfinite, qiu_nearest.residual_centers)
        @test all(isnan, qiu_nearest.residual_widths)
        @test all(isfinite, qiu_mwg.residual_centers)
        @test all(isfinite, qiu_mwg.residual_widths)
        @test all(>(0.0), vec(qiu_mwg.residual_widths))
        @test maximum(abs.(qiu_mwg.residual_centers)) < 1.0e-6
        @test isfinite(surrogate_check.orbital_energy)
        @test isfinite(nearest_check.orbital_energy)
        @test isfinite(mwg_check.orbital_energy)
        @test isfinite(surrogate_check.vee_expectation)
        @test isfinite(nearest_check.vee_expectation)
        @test isfinite(mwg_check.vee_expectation)
        @test surrogate_check.vee_expectation > 0.0
        @test nearest_check.vee_expectation > 0.0
        @test mwg_check.vee_expectation > 0.0

        (
            _full_source_basis,
            _full_legacy,
            qiu_full_nearest,
            qiu_full_nearest_check,
        ) = _qiu_white_full_nearest_fixture()

        @test norm(qiu_full_nearest.overlap - I, Inf) < 1.0e-8
        @test qiu_full_nearest.residual_count > 0
        @test -3.5 < qiu_full_nearest_check.orbital_energy < -1.8
        @test 1.0 < qiu_full_nearest_check.vee_expectation < 1.4
        @test abs(qiu_full_nearest_check.vee_expectation - 1.25) < 0.05
    end
end

end

if _test_group_enabled(:diatomic)
@testset "Bond-aligned diatomic QW reference path" begin
    basis14, operators14, check14 = _bond_aligned_diatomic_qw_fixture(; bond_length = 1.4)
    basis20, operators20, check20 = _bond_aligned_diatomic_qw_fixture(; bond_length = 2.0)

    @test basis14 isa BondAlignedDiatomicQWBasis3D
    @test operators14 isa QiuWhiteResidualGaussianOperators
    @test operators14.gaussian_data === nothing
    @test operators14.residual_count == 0
    @test operators14.gausslet_count == length(basis14.basis_x) * length(basis14.basis_y) * length(basis14.basis_z)
    @test norm(operators14.overlap - I, Inf) < 1.0e-8
    @test norm(operators20.overlap - I, Inf) < 1.0e-8
    @test operators14.one_body_hamiltonian ≈ transpose(operators14.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test operators14.interaction_matrix ≈ transpose(operators14.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test isfinite(check14.orbital_energy)
    @test isfinite(check14.vee_expectation)
    @test isfinite(check20.orbital_energy)
    @test isfinite(check20.vee_expectation)
    @test check14.orbital_energy < -1.0
    @test check20.orbital_energy < -1.0
    @test 0.5 < check14.vee_expectation < 1.0
    @test 0.5 < check20.vee_expectation < 1.0
    @test check14.orbital_energy < check20.orbital_energy
    @test length(basis14.basis_x) == length(basis14.basis_y)
    @test length(basis14.basis_z) > length(basis14.basis_x)
    @test length(basis20.basis_z) >= length(basis14.basis_z)
end

@testset "Per-center nuclear one-body reassembly on diatomic routes" begin
    basis, _operators, _check = _bond_aligned_diatomic_qw_fixture(; bond_length = 1.4)
    full_ops = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :by_center,
        interaction_treatment = :ggt_nearest,
    )
    atom_a_ops = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = [1.0, 0.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
    )
    atom_b_ops = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = [0.0, 1.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
    )

    @test full_ops.nuclear_term_storage == :by_center
    @test full_ops.nuclear_charges == [1.0, 1.0]
    @test !isnothing(full_ops.kinetic_one_body)
    @test !isnothing(full_ops.nuclear_one_body_by_center)
    @test length(full_ops.nuclear_one_body_by_center) == 2
    @test assembled_one_body_hamiltonian(full_ops) ≈
          full_ops.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test assembled_one_body_hamiltonian(full_ops; nuclear_charges = [1.0, 0.0]) ≈
          atom_a_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
    @test assembled_one_body_hamiltonian(full_ops; nuclear_charges = [0.0, 1.0]) ≈
          atom_b_ops.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
end

@testset "Nuclear term storage auto stays lightweight on longer-center routes" begin
    basis, _operators, _diagnostics = _bond_aligned_homonuclear_chain_qw_fixture(; natoms = 3)
    auto_ops = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = fill(1.0, length(basis.nuclei)),
        nuclear_term_storage = :auto,
        interaction_treatment = :ggt_nearest,
    )

    @test auto_ops.nuclear_term_storage == :total_only
    @test isnothing(auto_ops.kinetic_one_body)
    @test isnothing(auto_ops.nuclear_one_body_by_center)
    @test_throws ArgumentError assembled_one_body_hamiltonian(
        auto_ops;
        nuclear_charges = [1.0, 0.0, 1.0],
    )
end

@testset "Ordinary Cartesian naming surface distinguishes geometry from supplement" begin
    @test OrdinaryCartesianOrbital3D === QiuWhiteHybridOrbital3D
    @test OrdinaryCartesianOperators3D === QiuWhiteResidualGaussianOperators

    (
        direct_basis,
        _direct_parent_ops,
        _direct_parent_check,
        direct_supplement,
        direct_hybrid_ops,
        _direct_hybrid_check,
    ) = _bond_aligned_diatomic_hybrid_qw_fixture(; bond_length = 1.4)

    direct_via_clear_name = ordinary_cartesian_product_operators(
        direct_basis,
        direct_supplement;
        nuclear_charges = [1.0, 1.0],
        interaction_treatment = :ggt_nearest,
    )

    @test direct_via_clear_name isa OrdinaryCartesianOperators3D
    @test direct_via_clear_name.overlap ≈ direct_hybrid_ops.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test direct_via_clear_name.one_body_hamiltonian ≈ direct_hybrid_ops.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test direct_via_clear_name.interaction_matrix ≈ direct_hybrid_ops.interaction_matrix atol = 1.0e-12 rtol = 1.0e-12
    @test direct_via_clear_name.residual_count == direct_hybrid_ops.residual_count

    (
        _nested_basis,
        _nested_parent_ops,
        _nested_parent_check,
        nested_expansion,
        _nested_source,
        nested_fixed_block,
        nested_pure_ops,
        _nested_pure_check,
    ) = _bond_aligned_diatomic_nested_qw_fixture(; bond_length = 1.4)

    nested_via_clear_name = nested_cartesian_operators(
        nested_fixed_block;
        nuclear_charges = [1.0, 1.0],
        expansion = nested_expansion,
        interaction_treatment = :ggt_nearest,
    )

    @test nested_via_clear_name isa OrdinaryCartesianOperators3D
    @test nested_via_clear_name.gaussian_data === nothing
    @test nested_via_clear_name.overlap ≈ nested_pure_ops.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test nested_via_clear_name.one_body_hamiltonian ≈ nested_pure_ops.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test nested_via_clear_name.interaction_matrix ≈ nested_pure_ops.interaction_matrix atol = 1.0e-12 rtol = 1.0e-12
    @test nested_via_clear_name.residual_count == 0
end

@testset "Bond-aligned diatomic split geometry" begin
    basis, _operators, _check = _bond_aligned_diatomic_qw_fixture(; bond_length = 1.4)
    expansion = coulomb_gaussian_expansion(doacc = false)
    bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(basis, expansion)
    parent_box = (
        1:length(basis.basis_x),
        1:length(basis.basis_y),
        1:length(basis.basis_z),
    )
    working_box = (2:8, 2:8, 2:12)
    midpoint = 0.0

    geometry = GaussletBases._nested_bond_aligned_diatomic_split_geometry(
        bundles,
        parent_box,
        working_box;
        bond_axis = :z,
        midpoint = midpoint,
        nside = 5,
        min_parallel_to_transverse_ratio = 0.4,
    )
    sliver_geometry = GaussletBases._nested_bond_aligned_diatomic_split_geometry(
        bundles,
        parent_box,
        working_box;
        bond_axis = :z,
        midpoint = midpoint,
        nside = 5,
        min_parallel_to_transverse_ratio = 0.75,
    )
    short_geometry = GaussletBases._nested_bond_aligned_diatomic_split_geometry(
        bundles,
        parent_box,
        (3:7, 3:7, 4:12);
        bond_axis = :z,
        midpoint = midpoint,
        nside = 5,
        min_parallel_to_transverse_ratio = 0.4,
    )
    split_basis = bond_aligned_homonuclear_qw_basis(
        family = :G10,
        bond_length = 1.4,
        core_spacing = 0.3,
        xmax_parallel = 8.0,
        xmax_transverse = 6.0,
        bond_axis = :z,
        nuclear_charge = 2.0,
    )
    split_diagnostics = bond_aligned_diatomic_nested_geometry_diagnostics(
        split_basis;
        nside = 5,
    )
    split_bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(split_basis, expansion)
    split_parent_box = (
        1:length(split_basis.basis_x),
        1:length(split_basis.basis_y),
        1:length(split_basis.basis_z),
    )
    split_geometry = GaussletBases._nested_bond_aligned_diatomic_split_geometry(
        split_bundles,
        split_parent_box,
        split_diagnostics.geometry.working_box;
        bond_axis = :z,
        midpoint = midpoint,
        nside = 5,
        min_parallel_to_transverse_ratio = 0.4,
    )

    @test geometry.count_eligible
    @test !geometry.unsplit_aspect_eligible
    @test geometry.shape_eligible
    @test !geometry.did_split
    @test geometry.split_index == 7
    @test geometry.working_box == working_box
    @test isnothing(geometry.shared_midpoint_box)
    @test geometry.child_boxes == [(2:8, 2:8, 2:6), (2:8, 2:8, 8:12)]
    @test all(widths[3] >= 0.4 * max(widths[1], widths[2]) for widths in geometry.child_physical_widths)

    @test split_geometry.count_eligible
    @test split_geometry.unsplit_aspect_eligible
    @test split_geometry.shape_eligible
    @test split_geometry.did_split
    @test !isnothing(split_geometry.shared_midpoint_box)
    @test length(split_geometry.child_boxes) == 2

    @test sliver_geometry.count_eligible
    @test !sliver_geometry.shape_eligible
    @test !sliver_geometry.did_split

    @test !short_geometry.count_eligible
    @test !short_geometry.did_split
end

@testset "Bond-aligned heteronuclear split geometry" begin
    centers_axis = collect(-5.0:1.0:5.0)
    @test GaussletBases._nested_diatomic_split_plane_index(
        centers_axis,
        1:11,
        0.0;
        prefer_midpoint_tie_side = :left,
    ) == 6
    @test GaussletBases._nested_diatomic_split_plane_index(
        centers_axis,
        1:11,
        0.0;
        prefer_midpoint_tie_side = :right,
    ) == 5

    basis = bond_aligned_heteronuclear_qw_basis(
        atoms = ("He", "H"),
        bond_length = 1.45,
        core_spacings = (0.25, 0.5),
        nuclear_charges = (2.0, 1.0),
        xmax_parallel = 6.0,
        xmax_transverse = 4.0,
        bond_axis = :z,
    )
    expansion = coulomb_gaussian_expansion(doacc = false)
    bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(basis, expansion)
    parent_box = (
        1:length(basis.basis_x),
        1:length(basis.basis_y),
        1:length(basis.basis_z),
    )
    working_box = (2:(length(basis.basis_x) - 1), 2:(length(basis.basis_y) - 1), 2:(length(basis.basis_z) - 1))
    midpoint = sum(GaussletBases._qwrg_axis_coordinate(nucleus, :z) for nucleus in basis.nuclei) / 2
    geometry = GaussletBases._nested_bond_aligned_diatomic_split_geometry(
        bundles,
        parent_box,
        working_box;
        bond_axis = :z,
        midpoint = midpoint,
        nside = 5,
        min_parallel_to_transverse_ratio = 0.4,
        use_midpoint_slab = false,
        prefer_midpoint_tie_side = :left,
    )

    @test geometry.count_eligible
    @test !geometry.unsplit_aspect_eligible
    @test geometry.shape_eligible
    @test !geometry.did_split
    @test geometry.shared_midpoint_box === nothing
    @test length(geometry.child_boxes) == 2
    @test sum(length(box[3]) for box in geometry.child_boxes) == length(working_box[3])
    @test all(widths[3] >= 0.4 * max(widths[1], widths[2]) for widths in geometry.child_physical_widths)
end

@testset "Bond-aligned diatomic nested fixed block" begin
    (
        basis,
        operators,
        check,
        expansion,
        source,
        fixed_block,
        parent_modes,
        _parent_ground,
        _projected,
        projected_vee,
        capture,
        projected_energy,
    ) = _bond_aligned_diatomic_nested_fixed_block_fixture(; bond_length = 1.4)

    @test source isa GaussletBases._CartesianNestedBondAlignedDiatomicSource3D
    @test fixed_block isa GaussletBases._NestedFixedBlock3D
    @test !source.geometry.did_split
    @test length(source.shared_shell_layers) >= 1
    @test length(source.child_sequences) == 1
    @test isnothing(source.geometry.shared_midpoint_box)
    @test isnothing(source.midpoint_slab_column_range)
    @test length(source.child_column_ranges) == 1
    @test length(source.child_column_ranges[1]) == size(source.child_sequences[1].coefficient_matrix, 2)
    @test size(fixed_block.overlap, 1) == size(source.sequence.coefficient_matrix, 2)
    @test source.sequence.working_box == (
        1:length(basis.basis_x),
        1:length(basis.basis_y),
        1:length(basis.basis_z),
    )
    @test length(source.sequence.support_indices) ==
        length(basis.basis_x) * length(basis.basis_y) * length(basis.basis_z)
    @test size(fixed_block.coefficient_matrix, 1) ==
        length(basis.basis_x) * length(basis.basis_y) * length(basis.basis_z)
    @test size(fixed_block.coefficient_matrix, 2) < size(fixed_block.coefficient_matrix, 1)
    @test source.sequence.packet.term_storage == :compact_production
    @test fixed_block.term_storage == :compact_production
    @test isnothing(fixed_block.gaussian_terms)
    @test isnothing(fixed_block.pair_terms)
    @test !isnothing(fixed_block.gaussian_sum)
    @test !isnothing(fixed_block.pair_sum)
    @test norm(fixed_block.overlap - I, Inf) < 1.0e-10
    @test all(isfinite, fixed_block.weights)
    @test minimum(fixed_block.weights) > 0.0
    @test all(isfinite, fixed_block.fixed_centers)
    @test capture > 0.998
    @test projected_energy < -1.2
    @test abs(projected_energy - parent_modes.values[1]) < 0.03
    @test 0.7 < projected_vee < 0.8
    @test abs(projected_vee - check.vee_expectation) < 5.0e-4
end

@testset "Bond-aligned diatomic compact nested fixed-block contract" begin
    (
        basis,
        _operators,
        _check,
        expansion,
        default_source,
        default_fixed_block,
        _parent_modes,
        _parent_ground,
        _projected,
        _projected_vee,
        _capture,
        _projected_energy,
    ) = _bond_aligned_diatomic_nested_fixed_block_fixture(; bond_length = 1.4)

    term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
    explicit_nested = bond_aligned_diatomic_nested_fixed_block(
        basis;
        expansion = expansion,
        term_coefficients = term_coefficients,
    )
    explicit_fixed_block = explicit_nested.fixed_block

    @test default_source.sequence.packet.term_storage == :compact_production
    @test default_fixed_block.term_storage == :compact_production
    @test isnothing(default_fixed_block.gaussian_terms)
    @test isnothing(default_fixed_block.pair_terms)
    @test !isnothing(default_fixed_block.gaussian_sum)
    @test !isnothing(default_fixed_block.pair_sum)

    @test explicit_fixed_block.term_storage == :compact_production
    @test isnothing(explicit_fixed_block.gaussian_terms)
    @test isnothing(explicit_fixed_block.pair_terms)
    @test !isnothing(explicit_fixed_block.gaussian_sum)
    @test !isnothing(explicit_fixed_block.pair_sum)

    @test explicit_nested.source.sequence.packet.term_storage == :compact_production

    @test norm(default_fixed_block.overlap - explicit_fixed_block.overlap, Inf) < 1.0e-12
    @test norm(default_fixed_block.coefficient_matrix - explicit_fixed_block.coefficient_matrix, Inf) < 1.0e-12
    @test norm(default_fixed_block.gaussian_sum - explicit_fixed_block.gaussian_sum, Inf) < 1.0e-12
    @test norm(default_fixed_block.pair_sum - explicit_fixed_block.pair_sum, Inf) < 1.0e-12

    @test GaussletBases._qwrg_fixed_block_one_body_matrix(default_fixed_block, expansion; Z = 1.0) ≈
        GaussletBases._qwrg_fixed_block_one_body_matrix(explicit_fixed_block, expansion; Z = 1.0) atol =
        1.0e-10 rtol = 1.0e-10
    @test GaussletBases._qwrg_fixed_block_interaction_matrix(default_fixed_block, expansion) ≈
        GaussletBases._qwrg_fixed_block_interaction_matrix(explicit_fixed_block, expansion) atol =
        1.0e-10 rtol = 1.0e-10

    @test_throws MethodError bond_aligned_diatomic_nested_fixed_source(
        basis;
        expansion = expansion,
        term_storage = :full_debug,
    )
    @test_throws MethodError bond_aligned_diatomic_nested_fixed_block(
        basis;
        expansion = expansion,
        term_storage = :full_debug,
    )
end

@testset "Bond-aligned diatomic nested QW consumer path" begin
    (
        _basis,
        parent_ops,
        parent_check,
        _expansion,
        _source,
        fixed_block,
        nested_ops,
        nested_check,
    ) = _bond_aligned_diatomic_nested_qw_fixture(; bond_length = 1.4)

    @test nested_ops isa GaussletBases.QiuWhiteResidualGaussianOperators
    @test nested_ops.basis === fixed_block
    @test nested_ops.gaussian_data === nothing
    @test nested_ops.interaction_treatment == :ggt_nearest
    @test nested_ops.residual_count == 0
    @test nested_ops.gausslet_count == size(fixed_block.overlap, 1)
    @test size(nested_ops.residual_centers) == (0, 3)
    @test size(nested_ops.residual_widths) == (0, 3)
    @test norm(nested_ops.overlap - I, Inf) < 1.0e-10
    @test norm(nested_ops.one_body_hamiltonian - transpose(nested_ops.one_body_hamiltonian), Inf) < 1.0e-12
    @test norm(nested_ops.interaction_matrix - transpose(nested_ops.interaction_matrix), Inf) < 1.0e-12
    @test all(orbital.kind == :nested_fixed for orbital in orbitals(nested_ops))
    @test all(isfinite, fixed_block.weights)
    @test minimum(fixed_block.weights) > 0.0
    @test nested_check.overlap_error < 1.0e-10
    @test nested_check.orbital_energy < -1.2
    @test 0.7 < nested_check.vee_expectation < 0.8
    @test abs(nested_check.orbital_energy - parent_check.orbital_energy) < 0.03
    @test abs(nested_check.vee_expectation - parent_check.vee_expectation) < 0.03
    @test parent_ops.residual_count == 0
end

@testset "Bond-aligned diatomic molecular supplement ordinary QW path" begin
    fixture = _bond_aligned_diatomic_hybrid_qw_fixture(; bond_length = 1.4)
    if fixture === nothing
        @test true
    else
        (
            _basis,
            parent_ops,
            parent_check,
            supplement,
            ordinary_ops,
            ordinary_check,
        ) = fixture

        @test supplement isa LegacyBondAlignedDiatomicGaussianSupplement
        @test supplement.atomic_source.lmax == 1
        @test ordinary_ops.gaussian_data === supplement
        @test ordinary_ops.interaction_treatment == :ggt_nearest
        @test ordinary_ops.residual_count > 0
        @test size(ordinary_ops.residual_centers, 1) == ordinary_ops.residual_count
        @test size(ordinary_ops.residual_widths, 1) == ordinary_ops.residual_count
        @test norm(ordinary_ops.overlap - I, Inf) < 1.0e-10
        @test ordinary_check.overlap_error < 1.0e-10
        @test ordinary_check.orbital_energy < -1.25
        @test 0.75 < ordinary_check.vee_expectation < 0.81
        @test abs(ordinary_check.orbital_energy - parent_check.orbital_energy) < 0.01
        @test abs(ordinary_check.vee_expectation - parent_check.vee_expectation) < 0.01
        @test parent_ops.residual_count == 0
    end
end

@testset "Bond-aligned heteronuclear molecular supplement ordinary QW path" begin
    (
        basis,
        parent_ops,
        parent_check,
        supplement,
        ordinary_ops,
        ordinary_check,
    ) = _bond_aligned_heteronuclear_hybrid_qw_fixture(; bond_length = 1.45)

    @test supplement isa LegacyBondAlignedHeteronuclearGaussianSupplement
    @test supplement.atomic_sources[1].atom == "He"
    @test supplement.atomic_sources[2].atom == "H"
    @test supplement.atomic_sources[1].basis_name == "cc-pVTZ"
    @test supplement.atomic_sources[2].basis_name == "cc-pVTZ"
    @test supplement.atomic_sources[1].lmax == 1
    @test supplement.atomic_sources[2].lmax == 1
    @test ordinary_ops.gaussian_data === supplement
    @test ordinary_ops.interaction_treatment == :ggt_nearest
    @test ordinary_ops.residual_count > 0
    @test size(ordinary_ops.residual_centers, 1) == ordinary_ops.residual_count
    @test size(ordinary_ops.residual_widths, 1) == ordinary_ops.residual_count
    @test norm(ordinary_ops.overlap - I, Inf) < 1.0e-10
    @test ordinary_check.overlap_error < 1.0e-10
    @test ordinary_check.orbital_energy < -2.0
    @test 0.8 < ordinary_check.vee_expectation < 1.5
    @test ordinary_check.orbital_energy < parent_check.orbital_energy
    @test ordinary_check.vee_expectation > parent_check.vee_expectation

    payload = bond_aligned_diatomic_geometry_payload(ordinary_ops)
    slice = bond_aligned_diatomic_plane_slice(
        payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-12,
    )
    @test payload.bond_axis == :z
    @test length(payload.nuclei) == 2
    @test count(point -> point.group_kind == :residual_gaussian, payload.points) == ordinary_ops.residual_count
    @test slice.selected_count > 0
    @test all(abs(point.y) <= slice.plane_tol for point in slice.points)
    @test parent_ops.residual_count == 0
end

@testset "Bond-aligned heteronuclear nested fixed block" begin
    (
        basis,
        parent_ops,
        parent_check,
        _supplement,
        ordinary_ops,
        ordinary_check,
        _expansion,
        source,
        fixed_block,
        _parent_modes,
        _parent_ground,
        _projected,
        projected_vee,
        capture,
        projected_energy,
    ) = _bond_aligned_heteronuclear_nested_fixed_block_fixture(; bond_length = 1.45)

    @test source isa GaussletBases._CartesianNestedBondAlignedDiatomicSource3D
    @test fixed_block isa GaussletBases._NestedFixedBlock3D
    @test !source.geometry.did_split
    @test source.geometry.shared_midpoint_box === nothing
    @test isnothing(source.midpoint_slab_column_range)
    @test length(source.shared_shell_layers) >= 1
    @test length(source.child_sequences) == 1
    @test length(source.child_column_ranges) == 1
    @test source.sequence.working_box == (
        1:length(basis.basis_x),
        1:length(basis.basis_y),
        1:length(basis.basis_z),
    )
    @test norm(fixed_block.overlap - I, Inf) < 1.0e-10
    @test all(isfinite, fixed_block.weights)
    @test minimum(fixed_block.weights) > 0.0
    @test all(isfinite, fixed_block.fixed_centers)
    @test parent_check.overlap_error < 1.0e-10
    @test ordinary_check.overlap_error < 1.0e-10
    @test capture > 0.99
    @test projected_energy < -2.0
    @test projected_vee > 0.8
    @test size(fixed_block.overlap, 1) < ordinary_ops.gausslet_count
    @test size(fixed_block.overlap, 1) < parent_ops.gausslet_count
    @test source.geometry.child_boxes[1][3] != source.geometry.child_boxes[2][3]
end

@testset "Bond-aligned heteronuclear nested QW consumer path" begin
    (
        _basis,
        parent_ops,
        parent_check,
        source,
        fixed_block,
        supplement,
        nested_ops,
        nested_check,
    ) = _bond_aligned_heteronuclear_nested_hybrid_qw_fixture(; bond_length = 1.45)

    @test nested_ops isa GaussletBases.QiuWhiteResidualGaussianOperators
    @test nested_ops.basis === fixed_block
    @test nested_ops.gaussian_data === supplement
    @test nested_ops.interaction_treatment == :ggt_nearest
    @test nested_ops.residual_count > 0
    @test nested_check.overlap_error < 1.0e-8
    @test nested_check.orbital_energy < -2.0
    @test nested_check.vee_expectation > 0.8
    @test abs(nested_check.orbital_energy - parent_check.orbital_energy) < 0.01
    @test nested_check.vee_expectation > parent_check.vee_expectation
    payload = bond_aligned_diatomic_geometry_payload(nested_ops, source)
    slice = bond_aligned_diatomic_plane_slice(
        payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-5,
    )
    @test count(point -> point.group_kind == :shared_midpoint_slab, payload.points) == 0
    @test count(point -> point.group_kind == :residual_gaussian, payload.points) == nested_ops.residual_count
    @test slice.selected_count > 0
    @test parent_ops.residual_count == 0
end

@testset "Bond-aligned diatomic molecular supplement nested QW path" begin
    fixture = _bond_aligned_diatomic_nested_hybrid_qw_fixture(; bond_length = 1.4)
    if fixture === nothing
        @test true
    else
        (
            _basis,
            parent_ops,
            parent_check,
            _source,
            fixed_block,
            supplement,
            nested_ops,
            nested_check,
        ) = fixture

        @test supplement isa LegacyBondAlignedDiatomicGaussianSupplement
        @test nested_ops.basis === fixed_block
        @test nested_ops.gaussian_data === supplement
        @test nested_ops.interaction_treatment == :ggt_nearest
        @test nested_ops.residual_count > 0
        @test size(nested_ops.residual_centers, 1) == nested_ops.residual_count
        @test size(nested_ops.residual_widths, 1) == nested_ops.residual_count
        @test norm(nested_ops.overlap - I, Inf) < 1.0e-10
        @test nested_check.overlap_error < 1.0e-10
        @test all(isfinite, fixed_block.weights)
        @test minimum(fixed_block.weights) > 0.0
        @test nested_check.orbital_energy < -1.25
        @test 0.75 < nested_check.vee_expectation < 0.81
        @test abs(nested_check.orbital_energy - parent_check.orbital_energy) < 0.03
        @test abs(nested_check.vee_expectation - parent_check.vee_expectation) < 0.01
        @test parent_ops.residual_count == 0
    end
end

@testset "Bond-aligned diatomic geometry payloads and plane slices" begin
    basis, ordinary_ops, _check = _bond_aligned_diatomic_qw_fixture(; bond_length = 1.4)
    basis_payload = bond_aligned_diatomic_geometry_payload(basis)
    ordinary_payload = bond_aligned_diatomic_geometry_payload(ordinary_ops)
    basis_slice = bond_aligned_diatomic_plane_slice(
        basis_payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-12,
    )

    @test basis_payload isa GaussletBases.BondAlignedDiatomicGeometryPayload3D
    @test length(basis_payload.nuclei) == 2
    @test basis_payload.bond_axis == :z
    @test length(basis_payload.points) ==
        length(basis.basis_x) * length(basis.basis_y) * length(basis.basis_z)
    @test Set(point.group_kind for point in basis_payload.points) == Set([:gausslet_product])
    @test length(basis_payload.box_outlines) == 1
    @test basis_payload.box_outlines[1].group_kind == :basis_box
    @test ordinary_payload.bond_axis == :z
    @test length(ordinary_payload.points) == ordinary_ops.gausslet_count
    @test Set(point.group_kind for point in ordinary_payload.points) == Set([:gausslet_product])
    @test basis_slice.plane_axis == :y
    @test basis_slice.plane_value == 0.0
    @test basis_slice.plane_tol == 1.0e-12
    @test basis_slice.total_count == length(basis_payload.points)
    expected_y_count = count(y -> abs(y) <= 1.0e-12, centers(basis.basis_y))
    @test basis_slice.selected_count == length(basis.basis_x) * expected_y_count * length(basis.basis_z)
    @test all(abs(point.y) <= basis_slice.plane_tol for point in basis_slice.points)
    @test length(basis_slice.nuclei) == 2

    (
        _nested_basis,
        _parent_ops,
        _parent_check,
        _expansion,
        source,
        fixed_block,
        _nested_ops,
        _nested_check,
    ) = _bond_aligned_diatomic_nested_qw_fixture(; bond_length = 1.4)
    nested_payload = bond_aligned_diatomic_geometry_payload(fixed_block, source)
    expected_nested_groups =
        source.geometry.did_split ?
        Set([:shared_shell_layer, :left_child, :shared_midpoint_slab, :right_child]) :
        Set([:shared_shell_layer, :shared_child])
    expected_nested_box_count =
        2 + length(source.geometry.child_boxes) + (isnothing(source.geometry.shared_midpoint_box) ? 0 : 1)

    @test nested_payload isa GaussletBases.BondAlignedDiatomicGeometryPayload3D
    @test length(nested_payload.points) == size(fixed_block.fixed_centers, 1)
    @test Set(point.group_kind for point in nested_payload.points) == expected_nested_groups
    @test length(nested_payload.box_outlines) == expected_nested_box_count
    @test nested_payload.box_outlines[1].group_kind == :parent_box
    @test nested_payload.box_outlines[2].group_kind == :working_box
    @test count(box -> box.group_kind == :child_box, nested_payload.box_outlines) == length(source.geometry.child_boxes)
    @test any(box -> box.group_kind == :shared_midpoint_slab_box, nested_payload.box_outlines) ==
        !isnothing(source.geometry.shared_midpoint_box)

    hybrid_fixture = _bond_aligned_diatomic_hybrid_qw_fixture(; bond_length = 1.4)
    @test hybrid_fixture !== nothing
    if hybrid_fixture !== nothing
        (
            _hybrid_basis,
            _hybrid_parent_ops,
            _hybrid_parent_check,
            _supplement,
            hybrid_ops,
            _hybrid_check,
        ) = hybrid_fixture
        hybrid_payload = bond_aligned_diatomic_geometry_payload(hybrid_ops)
        @test count(point -> point.group_kind == :residual_gaussian, hybrid_payload.points) == hybrid_ops.residual_count
    end

    nested_hybrid_fixture = _bond_aligned_diatomic_nested_hybrid_qw_fixture(; bond_length = 1.4)
    @test nested_hybrid_fixture !== nothing
    if nested_hybrid_fixture !== nothing
        (
            _basis2,
            _parent_ops2,
            _parent_check2,
            hybrid_source,
            _hybrid_fixed_block,
            _hybrid_supplement,
            hybrid_nested_ops,
            _hybrid_nested_check,
        ) = nested_hybrid_fixture
        hybrid_nested_payload = bond_aligned_diatomic_geometry_payload(hybrid_nested_ops, hybrid_source)
        expected_hybrid_nested_groups =
            hybrid_source.geometry.did_split ?
            Set([:shared_shell_layer, :left_child, :shared_midpoint_slab, :right_child, :residual_gaussian]) :
            Set([:shared_shell_layer, :shared_child, :residual_gaussian])
        @test count(point -> point.group_kind == :residual_gaussian, hybrid_nested_payload.points) ==
            hybrid_nested_ops.residual_count
        @test Set(point.group_kind for point in hybrid_nested_payload.points) ==
            expected_hybrid_nested_groups
    end
end

@testset "Bond-aligned diatomic raw source geometry and 3d export" begin
    nested_hybrid_fixture = _bond_aligned_diatomic_nested_hybrid_qw_fixture(; bond_length = 1.4)
    @test nested_hybrid_fixture !== nothing
    if nested_hybrid_fixture !== nothing
        (
            _basis,
            _parent_ops,
            _parent_check,
            source,
            _fixed_block,
            _supplement,
            hybrid_nested_ops,
            _hybrid_nested_check,
        ) = nested_hybrid_fixture

        fixed_payload = bond_aligned_diatomic_geometry_payload(hybrid_nested_ops, source)
        source_payload = bond_aligned_diatomic_source_geometry_payload(source)
        expected_source_groups =
            source.geometry.did_split ?
            Set([:shared_shell_region, :left_child_region, :shared_midpoint_slab_region, :right_child_region]) :
            Set([:shared_shell_region, :shared_child_region])
        expected_child_box_count = length(source.geometry.child_boxes)
        expected_shared_midpoint_box = !isnothing(source.geometry.shared_midpoint_box)
        expected_shared_shell_points = sum(length(layer.support_indices) for layer in source.shared_shell_layers)
        expected_left_child_points = source.geometry.did_split ? length(source.child_sequences[1].support_indices) : 0
        expected_right_child_points = source.geometry.did_split ? length(source.child_sequences[2].support_indices) : 0
        expected_shared_child_points = source.geometry.did_split ? 0 : length(source.child_sequences[1].support_indices)
        expected_midpoint_points =
            isnothing(source.geometry.shared_midpoint_box) ?
            0 :
            prod(length.(source.geometry.shared_midpoint_box))

        @test all(isnothing(sequence.support_states) for sequence in source.child_sequences)
        @test Set(point.group_kind for point in source_payload.points) == expected_source_groups
        @test length(source_payload.points) == prod(length.(source.geometry.parent_box))
        @test count(point -> point.group_kind == :shared_shell_region, source_payload.points) ==
            expected_shared_shell_points
        @test count(point -> point.group_kind == :left_child_region, source_payload.points) ==
            expected_left_child_points
        @test count(point -> point.group_kind == :right_child_region, source_payload.points) ==
            expected_right_child_points
        @test count(point -> point.group_kind == :shared_child_region, source_payload.points) ==
            expected_shared_child_points
        @test count(point -> point.group_kind == :shared_midpoint_slab_region, source_payload.points) ==
            expected_midpoint_points
        @test any(box.group_kind == :parent_box for box in source_payload.box_outlines)
        @test any(box.group_kind == :working_box for box in source_payload.box_outlines)
        @test count(box -> box.group_kind == :child_box, source_payload.box_outlines) == expected_child_box_count
        @test any(box -> box.group_kind == :shared_midpoint_slab_box, source_payload.box_outlines) ==
            expected_shared_midpoint_box

        fixed_slice = bond_aligned_diatomic_plane_slice(
            fixed_payload;
            plane_axis = :y,
            plane_value = 0.0,
            plane_tol = 1.0e-12,
        )
        debug_slice = bond_aligned_diatomic_plane_slice(
            fixed_payload;
            plane_axis = :y,
            plane_value = 0.0,
            plane_tol = 1.0e-5,
        )
        source_slice = bond_aligned_diatomic_plane_slice(
            source_payload;
            plane_axis = :y,
            plane_value = 0.0,
            plane_tol = 1.0e-5,
        )
        @test debug_slice.selected_count >= fixed_slice.selected_count
        @test source_slice.selected_count > debug_slice.selected_count

        mktemp() do path, io
            close(io)
            payload = write_bond_aligned_diatomic_points3d(path, source_payload)
            text = read(path, String)
            @test payload === source_payload
            @test occursin("# bond_axis = z", text)
            @test occursin("# point_count = $(length(source_payload.points))", text)
            @test occursin("# nucleus_count = 2", text)
            @test occursin("# columns = x y z role kind group_kind group_id label", text)
            @test occursin("# box label=parent_box", text)
            @test occursin("# box label=working_box", text)
            if source.geometry.did_split
                @test occursin("# box label=left_child_box", text)
                @test occursin("# box label=right_child_box", text)
                @test occursin("# box label=shared_midpoint_slab_box", text)
            else
                @test occursin("# box label=shared_child_box", text)
                @test !occursin("# box label=shared_midpoint_slab_box", text)
            end
            @test occursin("\tpoint\tsource_region\tshared_shell_region\t1\t", text)
            if source.geometry.did_split
                @test occursin("\tpoint\tsource_region\tleft_child_region\t1\t", text)
                @test occursin("\tpoint\tsource_region\tshared_midpoint_slab_region\t1\t", text)
                @test occursin("\tpoint\tsource_region\tright_child_region\t2\t", text)
            else
                @test occursin("\tpoint\tsource_region\tshared_child_region\t1\t", text)
            end
            @test occursin("\tnucleus\tnucleus\tnucleus\t1\tA", text)
            @test occursin("\tnucleus\tnucleus\tnucleus\t2\tB", text)
        end
    end
end

@testset "Bond-aligned diatomic doside / COMX trace diagnostics" begin
    (
        _basis,
        _parent_ops,
        _parent_check,
        _expansion,
        source,
        _fixed_block,
        _nested_ops,
        _nested_check,
    ) = _bond_aligned_diatomic_nested_qw_fixture(; bond_length = 1.4)

    traces = GaussletBases._bond_aligned_diatomic_doside_traces(
        source;
        symmetry_tol = 1.0e-8,
        zero_tol = 1.0e-8,
    )
    lost_center = filter(
        trace -> trace.symmetric_about_zero && !trace.contains_near_zero_center,
        traces,
    )
    expected_context_suffixes = Set([
        "face_xy/tangential_x",
        "face_xy/tangential_y",
        "face_xz/tangential_x",
        "face_xz/tangential_z",
        "face_yz/tangential_y",
        "face_yz/tangential_z",
        "edge_x/free_axis_x",
        "edge_y/free_axis_y",
        "edge_z/free_axis_z",
    ])

    @test length(traces) == 9 * length(source.shared_shell_layers)
    @test all(trace.group_kind == :shared_shell for trace in traces)
    @test Set(trace.layer_index for trace in traces) ==
        Set(1:length(source.shared_shell_layers))
    @test all(trace.symmetric_about_zero for trace in traces)
    @test isempty(lost_center)
    @test all(trace.contains_near_zero_center for trace in traces)
    for layer_index in 1:length(source.shared_shell_layers)
        layer_traces = filter(trace -> trace.layer_index == layer_index, traces)
        @test length(layer_traces) == 9
        @test Set(
            replace(
                trace.context_label,
                "shared_shell/layer_$(layer_index)/" => "",
            ) for trace in layer_traces
        ) == expected_context_suffixes
    end

    mktemp() do path, io
        close(io)
        written = GaussletBases._write_bond_aligned_diatomic_doside_trace(
            path,
            source;
            symmetry_tol = 1.0e-8,
            zero_tol = 1.0e-8,
        )
        text = read(path, String)
        @test length(written) == length(traces)
        @test occursin("# trace_count = $(length(traces))", text)
        @test occursin(
            "# note shared_child has no local side contractions; it remains a direct core block",
            text,
        )
        @test !occursin(
            "# note left_child has no local side contractions; it remains a direct core block",
            text,
        )
        @test !occursin(
            "# note right_child has no local side contractions; it remains a direct core block",
            text,
        )
        @test occursin("context_label = shared_shell/layer_1/face_xy/tangential_x", text)
        @test occursin(
            "context_label = shared_shell/layer_$(length(source.shared_shell_layers))/edge_z/free_axis_z",
            text,
        )
        @test occursin("parent_centers = [", text)
        @test occursin("localized_centers = [", text)
        @test occursin("contains_near_zero_center = true", text)
        @test occursin("even_retained_count = false", text)
    end
end

@testset "Bond-aligned diatomic shared-shell odd-retain experiment" begin
    (
        _basis,
        _parent_ops,
        _parent_check,
        _expansion,
        baseline_source,
        baseline_fixed_block,
        _baseline_parent_modes,
        _baseline_parent_ground,
        _baseline_projected,
        baseline_projected_vee,
        baseline_capture,
        baseline_projected_energy,
        _baseline_supplement,
        baseline_ops,
        baseline_check,
    ) = _bond_aligned_diatomic_nested_hybrid_qw_shared_shell_experiment_fixture(
        ;
        bond_length = 1.4,
        shared_shell_retain_xy = (4, 3),
        shared_shell_retain_xz = (4, 3),
        shared_shell_retain_yz = (4, 3),
    )
    (
        _basis2,
        _parent_ops2,
        _parent_check2,
        _expansion2,
        experiment_source,
        experiment_fixed_block,
        _experiment_parent_modes,
        _experiment_parent_ground,
        _experiment_projected,
        experiment_projected_vee,
        experiment_capture,
        experiment_projected_energy,
        _experiment_supplement,
        experiment_ops,
        experiment_check,
    ) = _bond_aligned_diatomic_nested_hybrid_qw_shared_shell_experiment_fixture(
        ;
        bond_length = 1.4,
    )

    baseline_traces = GaussletBases._bond_aligned_diatomic_doside_traces(baseline_source)
    experiment_traces = GaussletBases._bond_aligned_diatomic_doside_traces(experiment_source)
    baseline_trace_map = Dict(trace.context_label => trace for trace in baseline_traces)
    trace_map = Dict(trace.context_label => trace for trace in experiment_traces)
    targeted_contexts = (
        "shared_shell/layer_1/face_xy/tangential_x",
        "shared_shell/layer_1/face_xz/tangential_x",
        "shared_shell/layer_1/face_yz/tangential_y",
    )
    fixed_payload = bond_aligned_diatomic_geometry_payload(baseline_ops, baseline_source)
    experiment_payload = bond_aligned_diatomic_geometry_payload(experiment_ops, experiment_source)
    fixed_slice = bond_aligned_diatomic_plane_slice(
        fixed_payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-5,
    )
    experiment_slice = bond_aligned_diatomic_plane_slice(
        experiment_payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-5,
    )

    @test norm(experiment_fixed_block.overlap - I, Inf) < 1.0e-10
    @test all(isfinite, experiment_fixed_block.weights)
    @test minimum(experiment_fixed_block.weights) > 0.0
    @test size(experiment_fixed_block.overlap, 1) < size(baseline_fixed_block.overlap, 1)
    @test Set(baseline_trace_map[context].contains_near_zero_center for context in targeted_contexts) ==
        Set([true])
    @test Set(trace_map[context].contains_near_zero_center for context in targeted_contexts) ==
        Set([true])
    @test Set(baseline_trace_map[context].retained_count for context in targeted_contexts) ==
        Set([5])
    @test Set(trace_map[context].retained_count for context in targeted_contexts) ==
        Set([3])
    @test fixed_slice.selected_count > experiment_slice.selected_count > 0
    @test count(point -> point.group_kind == :shared_shell_layer, fixed_slice.points) >
        count(point -> point.group_kind == :shared_shell_layer, experiment_slice.points) > 0
    @test abs(experiment_projected_vee - baseline_projected_vee) < 2.0e-5
    @test abs(experiment_capture - baseline_capture) < 2.0e-5
    @test abs(experiment_projected_energy - baseline_projected_energy) < 2.0e-5
    @test abs(experiment_check.orbital_energy - baseline_check.orbital_energy) < 2.0e-5
    @test abs(experiment_check.vee_expectation - baseline_check.vee_expectation) < 2.0e-5
end

@testset "Bond-aligned diatomic shared-shell odd-retain confirmation at R=2.0" begin
    (
        _basis,
        _parent_ops,
        _parent_check,
        _expansion,
        baseline_source,
        baseline_fixed_block,
        _baseline_parent_modes,
        _baseline_parent_ground,
        _baseline_projected,
        baseline_projected_vee,
        baseline_capture,
        baseline_projected_energy,
        _baseline_supplement,
        baseline_ops,
        baseline_check,
    ) = _bond_aligned_diatomic_nested_hybrid_qw_shared_shell_experiment_fixture(
        ;
        bond_length = 2.0,
        shared_shell_retain_xy = (4, 3),
        shared_shell_retain_xz = (4, 3),
        shared_shell_retain_yz = (4, 3),
    )
    (
        _basis2,
        _parent_ops2,
        _parent_check2,
        _expansion2,
        experiment_source,
        experiment_fixed_block,
        _experiment_parent_modes,
        _experiment_parent_ground,
        _experiment_projected,
        experiment_projected_vee,
        experiment_capture,
        experiment_projected_energy,
        _experiment_supplement,
        experiment_ops,
        experiment_check,
    ) = _bond_aligned_diatomic_nested_hybrid_qw_shared_shell_experiment_fixture(
        ;
        bond_length = 2.0,
    )

    baseline_traces = GaussletBases._bond_aligned_diatomic_doside_traces(baseline_source)
    experiment_traces = GaussletBases._bond_aligned_diatomic_doside_traces(experiment_source)
    baseline_trace_map = Dict(trace.context_label => trace for trace in baseline_traces)
    trace_map = Dict(trace.context_label => trace for trace in experiment_traces)
    targeted_contexts = (
        "shared_shell/layer_1/face_xy/tangential_x",
        "shared_shell/layer_1/face_xz/tangential_x",
        "shared_shell/layer_1/face_yz/tangential_y",
    )
    baseline_payload = bond_aligned_diatomic_geometry_payload(baseline_ops, baseline_source)
    experiment_payload = bond_aligned_diatomic_geometry_payload(experiment_ops, experiment_source)
    baseline_slice = bond_aligned_diatomic_plane_slice(
        baseline_payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-5,
    )
    experiment_slice = bond_aligned_diatomic_plane_slice(
        experiment_payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-5,
    )

    @test norm(experiment_fixed_block.overlap - I, Inf) < 1.0e-10
    @test all(isfinite, experiment_fixed_block.weights)
    @test minimum(experiment_fixed_block.weights) > 0.0
    @test size(experiment_fixed_block.overlap, 1) < size(baseline_fixed_block.overlap, 1)
    @test Set(baseline_trace_map[context].contains_near_zero_center for context in targeted_contexts) ==
        Set([true])
    @test Set(trace_map[context].contains_near_zero_center for context in targeted_contexts) ==
        Set([true])
    @test Set(baseline_trace_map[context].retained_count for context in targeted_contexts) ==
        Set([5])
    @test Set(trace_map[context].retained_count for context in targeted_contexts) ==
        Set([3])
    @test baseline_slice.selected_count > experiment_slice.selected_count > 0
    @test count(point -> point.group_kind == :shared_shell_layer, baseline_slice.points) >
        count(point -> point.group_kind == :shared_shell_layer, experiment_slice.points) > 0
    @test abs(experiment_projected_vee - baseline_projected_vee) < 2.0e-5
    @test abs(experiment_capture - baseline_capture) < 2.0e-5
    @test abs(experiment_projected_energy - baseline_projected_energy) < 2.0e-5
    @test abs(experiment_check.orbital_energy - baseline_check.orbital_energy) < 2.0e-5
    @test abs(experiment_check.vee_expectation - baseline_check.vee_expectation) < 2.0e-5
end

@testset "Bond-aligned diatomic doside boundary correction on larger debug box" begin
    basis = bond_aligned_homonuclear_qw_basis(
        bond_length = 1.4,
        core_spacing = 0.5,
        xmax_parallel = 8,
        xmax_transverse = 5,
        bond_axis = :z,
    )
    source = bond_aligned_diatomic_nested_fixed_source(basis)
    traces = GaussletBases._bond_aligned_diatomic_doside_traces(source)
    trace_map = Dict(trace.context_label => trace for trace in traces)

    @test length(traces) == 27
    symmetric_traces = filter(trace -> trace.symmetric_about_zero, traces)
    @test !isempty(symmetric_traces)
    @test all(trace -> trace.retained_count == 3, symmetric_traces)
    @test all(trace -> trace.contains_near_zero_center, symmetric_traces)
    @test trace_map["left_child/layer_1/face_xy/tangential_x"].contains_near_zero_center
    @test trace_map["left_child/layer_1/face_yz/tangential_y"].contains_near_zero_center
    @test trace_map["right_child/layer_1/face_xy/tangential_x"].contains_near_zero_center
    @test trace_map["right_child/layer_1/face_yz/tangential_y"].contains_near_zero_center

    fixed_block = GaussletBases._nested_fixed_block(source)
    supplement = legacy_bond_aligned_diatomic_gaussian_supplement(
        "H",
        "cc-pVTZ",
        basis.nuclei;
        lmax = 1,
    )
    ops = ordinary_cartesian_qiu_white_operators(
        fixed_block,
        supplement;
        nuclear_charges = [1.0, 1.0],
        interaction_treatment = :ggt_nearest,
    )
    payload = bond_aligned_diatomic_geometry_payload(ops, source)
    slice = bond_aligned_diatomic_plane_slice(
        payload;
        plane_axis = :y,
        plane_value = 0.0,
        plane_tol = 1.0e-5,
    )

    @test slice.selected_count == 121
    @test count(point -> point.group_kind == :left_child, slice.points) == 44
    @test count(point -> point.group_kind == :shared_midpoint_slab, slice.points) == 9
    @test count(point -> point.group_kind == :right_child, slice.points) == 44
    @test count(point -> point.group_kind == :shared_shell_layer, slice.points) == 16
end

@testset "Bond-aligned diatomic plane projection export" begin
    hybrid_fixture = _bond_aligned_diatomic_hybrid_qw_fixture(; bond_length = 1.4)
    @test hybrid_fixture !== nothing
    if hybrid_fixture !== nothing
        (
            _basis,
            _parent_ops,
            _parent_check,
            _supplement,
            hybrid_ops,
            _hybrid_check,
        ) = hybrid_fixture
        payload = bond_aligned_diatomic_geometry_payload(hybrid_ops)
        mktemp() do path, io
            close(io)
            slice = write_bond_aligned_diatomic_plane_projection(
                path,
                payload;
                plane_axis = :y,
                plane_value = 0.0,
                plane_tol = 1.0e-12,
            )
            text = read(path, String)
            @test slice.plane_axis == :y
            @test slice.plane_value == 0.0
            @test slice.plane_tol == 1.0e-12
            @test occursin("# plane_axis = y", text)
            @test occursin("# plane_value = 0.0", text)
            @test occursin("# plane_tol = 1.0e-12", text)
            @test occursin("# bond_axis = z", text)
            @test occursin("# selected_count = $(slice.selected_count)", text)
            @test occursin("# total_count = $(slice.total_count)", text)
            @test occursin("# projection_axes = x z", text)
            @test occursin("role=point group_kind=gausslet_product group_id=1", text)
            selected_residual_count = count(
                point -> point.group_kind == :residual_gaussian,
                slice.points,
            )
            @test count(==('@'), text) == 1 + selected_residual_count + 2
            @test sum(occursin("role=point group_kind=residual_gaussian", line) for line in split(text, '\n')) ==
                selected_residual_count
            gausslet_pos = findfirst("role=point group_kind=gausslet_product", text)
            residual_pos = findfirst("role=point group_kind=residual_gaussian", text)
            nucleus_pos = findfirst("role=nucleus group_kind=nucleus", text)
            @test gausslet_pos !== nothing
            @test residual_pos !== nothing
            @test nucleus_pos !== nothing
            @test gausslet_pos < residual_pos < nucleus_pos
        end
    end

    nested_hybrid_fixture = _bond_aligned_diatomic_nested_hybrid_qw_fixture(; bond_length = 1.4)
    @test nested_hybrid_fixture !== nothing
    if nested_hybrid_fixture !== nothing
        (
            _basis2,
            _parent_ops2,
            _parent_check2,
            source,
            _fixed_block,
            _supplement2,
            hybrid_nested_ops,
            _nested_check2,
        ) = nested_hybrid_fixture
        payload = bond_aligned_diatomic_geometry_payload(hybrid_nested_ops, source)
        mktemp() do path, io
            close(io)
            slice = write_bond_aligned_diatomic_plane_projection(
                path,
                payload;
                plane_axis = :y,
                plane_value = 0.0,
                plane_tol = 1.0e-12,
            )
            text = read(path, String)
            @test slice.selected_count <= slice.total_count
            @test occursin("role=point group_kind=left_child group_id=1", text)
            @test occursin("role=point group_kind=shared_midpoint_slab group_id=1", text)
            @test occursin("role=point group_kind=right_child group_id=2", text)
            @test occursin("role=point group_kind=shared_shell_layer group_id=1", text)
            selected_residual_count = count(
                point -> point.group_kind == :residual_gaussian,
                slice.points,
            )
            @test sum(occursin("role=point group_kind=residual_gaussian", line) for line in split(text, '\n')) ==
                selected_residual_count
            left_pos = findfirst("role=point group_kind=left_child", text)
            slab_pos = findfirst("role=point group_kind=shared_midpoint_slab", text)
            right_pos = findfirst("role=point group_kind=right_child", text)
            shared_pos = findfirst("role=point group_kind=shared_shell_layer", text)
            residual_pos = findfirst("role=point group_kind=residual_gaussian", text)
            nucleus_pos = findfirst("role=nucleus group_kind=nucleus", text)
            @test left_pos !== nothing
            @test slab_pos !== nothing
            @test right_pos !== nothing
            @test shared_pos !== nothing
            @test residual_pos !== nothing
            @test nucleus_pos !== nothing
            @test left_pos < slab_pos < right_pos < shared_pos < residual_pos < nucleus_pos
            @test count(==('@'), text) == 4 + selected_residual_count + 2
        end
    end
end

end

if _test_group_enabled(:ordinary)
@testset "Ordinary Cartesian localized backend" begin
    expansion = _truncate_coulomb_expansion(coulomb_gaussian_expansion(doacc = false), 3)
    basis = build_basis(MappedUniformBasisSpec(:G10;
        count = 5,
        mapping = fit_asinh_mapping_for_strength(s = 0.5, npoints = 5, xmax = 6.0),
        reference_spacing = 1.0,
    ))

    proxy = ordinary_cartesian_ida_operators(
        basis;
        expansion = expansion,
        Z = 2.0,
        backend = :pgdg_experimental,
    )
    localized = ordinary_cartesian_ida_operators(
        basis;
        expansion = expansion,
        Z = 2.0,
        backend = :pgdg_localized_experimental,
    )
    (_, _, overlap_reference_localized, kinetic_reference_localized, gaussian_reference_localized) =
        _localized_numerical_reference_1d(basis, expansion.exponents)
    raw_localized = mapped_pgdg_localized(GaussletBases.mapped_pgdg_logfit_prototype(basis))
    oracle_localized = GaussletBases._mapped_ordinary_localized_oracle_operators(
        basis;
        exponents = expansion.exponents,
        center = 0.0,
    )
    raw_h1, _ = _cartesian_hydrogen_energy(
        overlap_matrix(raw_localized),
        kinetic_matrix(raw_localized),
        gaussian_factor_matrices(raw_localized; exponents = expansion.exponents, center = 0.0),
        expansion;
        Z = 2.0,
    )
    corrected_h1, _ = _cartesian_hydrogen_energy(
        localized.one_body_1d.overlap,
        localized.one_body_1d.kinetic,
        localized.one_body_1d.gaussian_factors,
        expansion;
        Z = 2.0,
    )
    oracle_h1, _ = _cartesian_hydrogen_energy(
        oracle_localized.overlap,
        oracle_localized.kinetic,
        oracle_localized.gaussian_factors,
        expansion;
        Z = 2.0,
    )
    reference_h1, _ = _cartesian_hydrogen_energy(
        overlap_reference_localized,
        kinetic_reference_localized,
        gaussian_reference_localized,
        expansion;
        Z = 2.0,
    )

    @test localized isa OrdinaryCartesianIDAOperators
    @test localized.backend == :pgdg_localized_experimental
    @test occursin("experimental=true", sprint(show, localized))
    @test norm(localized.one_body_1d.overlap - I, Inf) < 1.0e-10
    @test norm(localized.overlap_3d - I, Inf) < 1.0e-9
    @test norm(localized.one_body_1d.overlap - I, Inf) < norm(proxy.one_body_1d.overlap - I, Inf)
    @test norm(localized.overlap_3d - I, Inf) < norm(proxy.overlap_3d - I, Inf)
    @test norm(localized.one_body_1d.overlap - overlap_reference_localized, Inf) < 1.0e-10
    @test norm(localized.one_body_1d.kinetic - kinetic_reference_localized, Inf) <
          norm(kinetic_matrix(raw_localized) - kinetic_reference_localized, Inf)
    @test norm(corrected_h1 - reference_h1, Inf) < norm(raw_h1 - reference_h1, Inf)
    @test norm(oracle_h1 - reference_h1, Inf) < norm(corrected_h1 - reference_h1, Inf)
    @test localized.one_body_hamiltonian ≈ transpose(localized.one_body_hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test localized.interaction_matrix ≈ transpose(localized.interaction_matrix) atol = 1.0e-10 rtol = 1.0e-10
    @test minimum(diag(localized.interaction_matrix)) > 0.0
end

@testset "Legacy 1D hybrid mapped ordinary basis" begin
    (
        basis,
        core_gaussians,
        expansion,
        hybrid_reference,
        hybrid_analytic,
        one_body_reference,
        one_body_analytic,
        hard_reference,
        hard_analytic,
        energy_reference,
        energy_analytic,
        hard_energy_reference,
        hard_energy_analytic,
    ) = _quick_hybrid_mapped_ordinary_fixture()

    @test hybrid_reference isa HybridMappedOrdinaryBasis1D
    @test hybrid_analytic isa HybridMappedOrdinaryBasis1D
    @test hybrid_reference.backend == :numerical_reference
    @test hybrid_analytic.backend == :pgdg_localized_experimental
    @test length(core_gaussians) == 2
    @test length(hybrid_reference) == length(basis) + length(core_gaussians)
    @test length(hybrid_analytic) == length(basis) + length(core_gaussians)
    @test occursin("experimental=true", sprint(show, hybrid_analytic))
    @test !occursin("experimental=true", sprint(show, hybrid_reference))
    @test norm(one_body_reference.overlap - I, Inf) < 1.0e-10
    @test norm(one_body_analytic.overlap - I, Inf) < 1.0e-10
    @test one_body_analytic.kinetic ≈ transpose(one_body_analytic.kinetic) atol = 1.0e-10 rtol = 1.0e-10
    hybrid_factor_diff = maximum(
        norm(one_body_analytic.gaussian_factors[index] - one_body_reference.gaussian_factors[index], Inf)
        for index in eachindex(one_body_analytic.gaussian_factors)
    )
    hard_factor_diff = maximum(
        norm(hard_analytic.gaussian_factors[index] - hard_reference.gaussian_factors[index], Inf)
        for index in eachindex(hard_analytic.gaussian_factors)
    )
    @test hybrid_factor_diff < hard_factor_diff
    @test abs(energy_analytic - energy_reference) <
          abs(hard_energy_analytic - hard_energy_reference)
end

@testset "Legacy 1D hybrid ordinary mapped SHO smoke" begin
    (
        hybrid_analytic,
        centered_analytic,
        centered_hamiltonian,
    ) = _quick_ordinary_sho_smoke_fixture()

    @test occursin("experimental=true", sprint(show, hybrid_analytic))
    @test centered_hamiltonian.backend == :pgdg_localized_experimental
    @test centered_hamiltonian.overlap ≈ transpose(centered_hamiltonian.overlap) atol = 1.0e-10 rtol = 1.0e-10
    @test centered_hamiltonian.hamiltonian ≈ transpose(centered_hamiltonian.hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test norm(centered_hamiltonian.overlap - I, Inf) < 1.0e-10
    @test issorted(centered_analytic.eigenvalues)
    @test centered_analytic.eigenvalues[1] > 0.0
    @test abs(centered_analytic.eigenvalues[1] - centered_analytic.exact[1]) < 1.0e-2
    @test centered_analytic.kinetic_expectation > 0.0
    @test centered_analytic.displacement2_expectation > 0.0
end

if _RUN_SLOW_TESTS
    @testset "Ordinary mapped SHO spectra" begin
        (
            centered_reference,
            centered_analytic,
            shifted_reference,
            shifted_analytic,
            stress_reference,
            stress_analytic,
        ) = _slow_ordinary_sho_fixture()

        centered_diff = maximum(abs.(centered_reference.eigenvalues .- centered_analytic.eigenvalues))
        shifted_diff = maximum(abs.(shifted_reference.eigenvalues .- shifted_analytic.eigenvalues))
        stress_diff = maximum(abs.(stress_reference.eigenvalues .- stress_analytic.eigenvalues))

        @test centered_diff < 2.0e-4
        @test shifted_diff < 2.0e-4
        @test stress_diff > 1.0e-2
        @test maximum(abs.(centered_analytic.eigenvalues .- centered_analytic.exact)) <
              maximum(abs.(stress_analytic.eigenvalues .- stress_analytic.exact))
        @test maximum(
            abs.(
                abs.(centered_analytic.eigenvalues .- centered_analytic.exact) .-
                abs.(centered_reference.eigenvalues .- centered_reference.exact),
            ),
        ) < 5.0e-4
        @test maximum(
            abs.(
                abs.(shifted_analytic.eigenvalues .- shifted_analytic.exact) .-
                abs.(shifted_reference.eigenvalues .- shifted_reference.exact),
            ),
        ) < 5.0e-4
    end

    @testset "Mapped Cartesian hydrogen" begin
        basis, mapping_value, representation, overlap_1d, kinetic_1d, expansion, hamiltonian, energy =
            _quick_mapped_cartesian_hydrogen_fixture()
        (_, _, _, _, _, _, _, _, _, unmapped_energy) = _quick_cartesian_hydrogen_fixture()

        @test mapping_value isa AsinhMapping
        @test overlap_1d ≈ transpose(overlap_1d) atol = 1.0e-10 rtol = 1.0e-10
        @test kinetic_1d ≈ transpose(kinetic_1d) atol = 1.0e-10 rtol = 1.0e-10
        @test length(expansion) == 45
        @test size(hamiltonian) == (length(basis)^3, length(basis)^3)
        @test energy < unmapped_energy - 0.02
        @test energy < -0.46
        @test energy > -0.5
    end

    @testset "Mapped Coulomb expansion calibration" begin
        expansion_115 = coulomb_gaussian_expansion(doacc = true, maxu = 115.0)
        expansion_135 = coulomb_gaussian_expansion(doacc = true, maxu = 135.0)

        (_, _, _, _, _, energy_115, _, _, _, _) = _mapped_cartesian_hydrogen_comparison(expansion_115)
        (_, _, _, _, _, energy_135, _, _, _, _) = _mapped_cartesian_hydrogen_comparison(expansion_135)

        @test length(expansion_115) == 115
        @test length(expansion_135) == 135
        @test abs(energy_135 - energy_115) < 1.0e-8
    end

    @testset "Mapped PGDG hydrogen prototype" begin
        expansion = coulomb_gaussian_expansion(doacc = false)
        (
            basis,
            prototype,
            localized,
            refined,
            refined_localized,
            overlap_numeric,
            kinetic_numeric,
            hamiltonian_numeric,
            energy_numeric,
            overlap_pgdg,
            position_pgdg,
            kinetic_pgdg,
            hamiltonian_pgdg,
            energy_pgdg,
            overlap_localized,
            position_localized,
            kinetic_localized,
            hamiltonian_localized,
            energy_localized,
            overlap_refined,
            position_refined,
            kinetic_refined,
            hamiltonian_refined,
            energy_refined,
            overlap_refined_localized,
            position_refined_localized,
            kinetic_refined_localized,
            hamiltonian_refined_localized,
            energy_refined_localized,
        ) = _mapped_cartesian_hydrogen_comparison(expansion)

        representation_numeric = basis_representation(basis; operators = (:overlap, :position, :kinetic))
        transform_numeric, centers_numeric = _cleanup_comx_transform(
            representation_numeric.basis_matrices.overlap,
            representation_numeric.basis_matrices.position,
            integral_weights(basis),
        )
        overlap_numeric_localized = transform_numeric' * representation_numeric.basis_matrices.overlap * transform_numeric
        kinetic_numeric_localized = transform_numeric' * representation_numeric.basis_matrices.kinetic * transform_numeric
        position_numeric_localized = transform_numeric' * representation_numeric.basis_matrices.position * transform_numeric
        factor_numeric = gaussian_factor_matrix(basis; exponent = 0.35, center = 0.0)
        factor_pgdg = gaussian_factor_matrix(prototype; exponent = 0.35, center = 0.0)
        factor_refined = gaussian_factor_matrix(refined; exponent = 0.35, center = 0.0)
        factor_numeric_localized = transform_numeric' * factor_numeric * transform_numeric
        factor_localized = gaussian_factor_matrix(localized; exponent = 0.35, center = 0.0)
        factor_refined_localized = gaussian_factor_matrix(refined_localized; exponent = 0.35, center = 0.0)
        span_pre = _subspace_overlap_metric(basis, prototype)
        span_refined = _subspace_overlap_metric(basis, refined)
        projector_pre = _projector_difference_metric(basis, prototype)
        projector_refined = _projector_difference_metric(basis, refined)

        @test size(hamiltonian_numeric) == size(hamiltonian_pgdg)
        @test size(hamiltonian_numeric) == size(hamiltonian_localized)
        @test size(hamiltonian_numeric) == size(hamiltonian_refined)
        @test size(hamiltonian_numeric) == size(hamiltonian_refined_localized)
        @test overlap_pgdg ≈ transpose(overlap_pgdg) atol = 1.0e-10 rtol = 1.0e-10
        @test kinetic_pgdg ≈ transpose(kinetic_pgdg) atol = 1.0e-10 rtol = 1.0e-10
        @test overlap_localized ≈ I atol = 1.0e-10 rtol = 1.0e-10
        @test position_localized ≈ Diagonal(centers(localized)) atol = 1.0e-10 rtol = 1.0e-10
        @test overlap_numeric_localized ≈ I atol = 1.0e-10 rtol = 1.0e-10
        @test position_numeric_localized ≈ Diagonal(centers_numeric) atol = 1.0e-10 rtol = 1.0e-10
        @test overlap_refined ≈ transpose(overlap_refined) atol = 1.0e-10 rtol = 1.0e-10
        @test kinetic_refined ≈ transpose(kinetic_refined) atol = 1.0e-10 rtol = 1.0e-10
        @test overlap_refined_localized ≈ I atol = 1.0e-10 rtol = 1.0e-10
        @test position_refined_localized ≈ Diagonal(centers(refined_localized)) atol = 1.0e-10 rtol = 1.0e-10
        @test norm(kinetic_numeric_localized - kinetic_localized, Inf) < norm(kinetic_numeric - kinetic_pgdg, Inf)
        @test norm(factor_numeric_localized - factor_localized, Inf) < norm(factor_numeric - factor_pgdg, Inf)
        @test energy_numeric < -0.47
        @test energy_pgdg < -0.47
        @test energy_localized < -0.47
        @test energy_refined < -0.47
        @test energy_refined_localized < -0.47
        @test abs(energy_numeric - energy_pgdg) < 0.01
        @test abs(energy_localized - energy_pgdg) < 1.0e-10
        @test abs(energy_numeric - energy_localized) ≤ abs(energy_numeric - energy_pgdg) + 1.0e-10
        @test span_refined.min_sv > span_pre.min_sv
        @test projector_refined.frob < projector_pre.frob
        @test projector_refined.op < projector_pre.op
        @test abs(energy_numeric - energy_refined) < abs(energy_numeric - energy_pgdg)
        @test abs(energy_numeric - energy_refined_localized) ≤ abs(energy_numeric - energy_refined) + 1.0e-10
        @test norm(factor_numeric - factor_refined, Inf) ≤ norm(factor_numeric - factor_pgdg, Inf) + 0.02
        @test norm(factor_numeric_localized - factor_refined_localized, Inf) ≤ norm(factor_numeric_localized - factor_localized, Inf) + 0.02
    end

    @testset "Mapped ordinary backend hydrogen regimes" begin
        expansion = coulomb_gaussian_expansion(doacc = false)

        function backend_energy(s_value)
            basis = build_basis(MappedUniformBasisSpec(:G10;
                count = 5,
                mapping = fit_asinh_mapping_for_strength(s = s_value, npoints = 5, xmax = 6.0),
                reference_spacing = 1.0,
            ))
            reference = mapped_ordinary_one_body_operators(
                basis;
                exponents = expansion.exponents,
                backend = :numerical_reference,
            )
            analytic = mapped_ordinary_one_body_operators(
                basis;
                exponents = expansion.exponents,
                backend = :pgdg_experimental,
            )
            return (
                basis,
                mapped_cartesian_hydrogen_energy(reference, expansion; Z = 1.0),
                mapped_cartesian_hydrogen_energy(analytic, expansion; Z = 1.0),
            )
        end

        mild_basis, mild_energy_reference, mild_energy_analytic = backend_energy(0.5)
        moderate_basis, moderate_energy_reference, moderate_energy_analytic = backend_energy(1.0)
        stress_basis, stress_energy_reference, stress_energy_analytic = backend_energy(2.0)

        mild_diff = abs(mild_energy_reference - mild_energy_analytic)
        moderate_diff = abs(moderate_energy_reference - moderate_energy_analytic)
        stress_diff = abs(stress_energy_reference - stress_energy_analytic)

        @test mild_basis isa MappedUniformBasis
        @test moderate_basis isa MappedUniformBasis
        @test stress_basis isa MappedUniformBasis
        @test mild_diff < 3.0e-4
        @test moderate_diff < 1.0e-3
        @test stress_diff > moderate_diff
        @test stress_diff < 0.02
        @test stress_energy_analytic < -0.45
        @test stress_energy_reference < -0.45
    end

    @testset "Ordinary Cartesian IDA backend regimes" begin
        expansion = _truncate_coulomb_expansion(coulomb_gaussian_expansion(doacc = false), 3)

        function ida_difference(s_value)
            basis = build_basis(MappedUniformBasisSpec(:G10;
                count = 5,
                mapping = fit_asinh_mapping_for_strength(s = s_value, npoints = 5, xmax = 6.0),
                reference_spacing = 1.0,
            ))
            reference = ordinary_cartesian_ida_operators(
                basis;
                expansion = expansion,
                Z = 2.0,
                backend = :numerical_reference,
            )
            analytic = ordinary_cartesian_ida_operators(
                basis;
                expansion = expansion,
                Z = 2.0,
                backend = :pgdg_experimental,
            )
            return (
                norm(reference.one_body_hamiltonian - analytic.one_body_hamiltonian, Inf),
                norm(reference.interaction_matrix - analytic.interaction_matrix, Inf),
            )
        end

        mild_h1_diff, mild_vee_diff = ida_difference(0.5)
        stress_h1_diff, stress_vee_diff = ida_difference(2.0)

        @test mild_h1_diff < 0.05
        @test mild_vee_diff < 0.1
        @test stress_h1_diff > mild_h1_diff
        @test stress_vee_diff > mild_vee_diff
        @test stress_vee_diff < 1.0
    end

    @testset "Ordinary Cartesian IDA reference agreement" begin
        mild_basis, mild_expansion, mild_reference = _quick_ordinary_cartesian_ida_fixture(
            backend = :numerical_reference,
            mapped = true,
            s = 0.5,
        )
        (_, _, mild_analytic) = _quick_ordinary_cartesian_ida_fixture(
            backend = :pgdg_experimental,
            mapped = true,
            s = 0.5,
        )
        (_, _, identity_reference) = _quick_ordinary_cartesian_ida_fixture(
            backend = :numerical_reference,
            mapped = false,
        )
        (_, _, identity_analytic) = _quick_ordinary_cartesian_ida_fixture(
            backend = :pgdg_experimental,
            mapped = false,
        )

        @test mild_basis isa MappedUniformBasis
        @test mild_reference.backend == :numerical_reference
        @test !occursin("experimental=true", sprint(show, mild_reference))
        @test norm(identity_reference.one_body_hamiltonian - identity_analytic.one_body_hamiltonian, Inf) < 1.0e-6
        @test norm(identity_reference.interaction_matrix - identity_analytic.interaction_matrix, Inf) < 1.0e-6
        @test norm(mild_reference.one_body_hamiltonian - mild_analytic.one_body_hamiltonian, Inf) < 0.01
        @test norm(mild_reference.interaction_matrix - mild_analytic.interaction_matrix, Inf) < 1.0e-3
    end

    @testset "Ordinary Cartesian localized backend agreement" begin
        expansion = _truncate_coulomb_expansion(coulomb_gaussian_expansion(doacc = false), 3)

        function build_pair(s_value)
            basis = build_basis(MappedUniformBasisSpec(:G10;
                count = 5,
                mapping = fit_asinh_mapping_for_strength(s = s_value, npoints = 5, xmax = 6.0),
                reference_spacing = 1.0,
            ))
            reference = ordinary_cartesian_ida_operators(
                basis;
                expansion = expansion,
                Z = 2.0,
                backend = :numerical_reference,
            )
            localized = ordinary_cartesian_ida_operators(
                basis;
                expansion = expansion,
                Z = 2.0,
                backend = :pgdg_localized_experimental,
            )
            return reference, localized
        end

        mild_reference, mild_localized = build_pair(0.5)
        stress_reference, stress_localized = build_pair(2.0)

        mild_h1_diff = norm(mild_reference.one_body_hamiltonian - mild_localized.one_body_hamiltonian, Inf)
        mild_vee_diff = norm(mild_reference.interaction_matrix - mild_localized.interaction_matrix, Inf)
        stress_h1_diff = norm(stress_reference.one_body_hamiltonian - stress_localized.one_body_hamiltonian, Inf)
        stress_vee_diff = norm(stress_reference.interaction_matrix - stress_localized.interaction_matrix, Inf)

        @test norm(mild_localized.one_body_1d.overlap - I, Inf) < 1.0e-10
        @test norm(mild_localized.overlap_3d - I, Inf) < 1.0e-9
        @test mild_h1_diff < 0.02
        @test mild_vee_diff < 1.0e-2
        @test stress_h1_diff > mild_h1_diff
        @test stress_vee_diff > mild_vee_diff
        @test stress_h1_diff < 0.2
        @test stress_vee_diff < 0.2
    end
end

@testset "Cartesian hydrogen via Coulomb expansion" begin
    basis, representation, expansion, overlap_1d, kinetic_1d, gaussian_factors, overlap_3d, nuclear_3d, hamiltonian, energy =
        _quick_cartesian_hydrogen_fixture()

    @test length(basis) == 5
    @test size(overlap_1d) == (length(basis), length(basis))
    @test norm(overlap_1d - I, Inf) ≤ 1.0e-10
    @test size(overlap_3d) == (length(basis)^3, length(basis)^3)
    @test norm(overlap_3d - I, Inf) ≤ 1.0e-8
    @test hamiltonian ≈ transpose(hamiltonian) atol = 1.0e-10 rtol = 1.0e-10
    @test size(nuclear_3d) == size(hamiltonian)
    @test minimum(diag(nuclear_3d)) < 0.0
    @test -0.7 < energy < -0.3
end

end

if _test_group_enabled(:core)
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

end

if _test_group_enabled(:radial)
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

end

if _test_group_enabled(:angular)
@testset "Vendored angular sphere-point access" begin
    orders = sphere_point_set_orders()
    curated_orders = curated_sphere_point_set_orders()

    @test first(orders) == 10
    @test last(orders) == 580
    @test length(orders) == 122
    @test issorted(orders)
    @test 15 in orders
    @test 32 in orders
    @test 51 in orders
    @test 100 in orders
    @test 580 in orders
    @test_throws ArgumentError sphere_point_set(9)
    @test_throws ArgumentError sphere_point_set(101)

    set10 = sphere_point_set(10)
    set15_full = sphere_point_set(15)
    set100 = sphere_point_set(100)
    set580 = sphere_point_set(580)

    for set in (set10, set15_full, set100, set580)
        @test set isa SpherePointSet
        @test set.cardinality == set.order
        @test size(set.coordinates) == (set.cardinality, 3)
        @test all(isfinite, set.coordinates)
        @test set.nn_ratio >= 1.0
        @test set.provenance.source_tag == "optimized_sphere_points_full_vendor"
        @test set.provenance.source_project == "GaussletModules/Radial"
        @test occursin("SpherePoints.jld2", set.provenance.source_artifact)
        @test occursin("xyzsets", set.provenance.source_note)

        norms = sqrt.(sum(abs2, set.coordinates; dims = 2))
        @test maximum(abs.(norms .- 1.0)) < 1.0e-12
    end

    @test all(order in orders for order in curated_orders)
    @test set15_full.coordinates ≈ curated_sphere_point_set(15).coordinates atol = 0.0 rtol = 1.0e-14
    @test set15_full.nn_ratio ≈ curated_sphere_point_set(15).nn_ratio atol = 0.0 rtol = 1.0e-14
end

@testset "Curated angular sphere-point access" begin
    orders = curated_sphere_point_set_orders()

    @test orders == [15, 32, 51]
    @test_throws ArgumentError curated_sphere_point_set(14)

    set15 = curated_sphere_point_set(15)
    set32 = curated_sphere_point_set(32)
    set51 = curated_sphere_point_set(51)

    for set in (set15, set32, set51)
        @test set isa CuratedSpherePointSet
        @test set.cardinality == set.order
        @test size(set.coordinates) == (set.cardinality, 3)
        @test all(isfinite, set.coordinates)
        @test set.nn_ratio >= 1.0
        @test set.provenance.source_tag == "optimized_sphere_points_curated_subset"
        @test set.provenance.source_project == "GaussletModules/Radial"
        @test occursin("SpherePoints.jld2", set.provenance.source_artifact)
        @test occursin("xyzsets", set.provenance.source_note)

        norms = sqrt.(sum(abs2, set.coordinates; dims = 2))
        @test maximum(abs.(norms .- 1.0)) < 1.0e-12
    end

    @test set15.nn_ratio ≈ 1.1233641689852316 atol = 0.0 rtol = 1.0e-14
    @test set32.nn_ratio ≈ 1.0000000159730313 atol = 0.0 rtol = 1.0e-14
    @test set51.nn_ratio ≈ 1.0802896662822246 atol = 0.0 rtol = 1.0e-14
    @test set15.coordinates[1, :] ≈ [0.5449408412377406, -0.022962554521322145, 0.838160008972606] atol = 0.0 rtol = 1.0e-14
    @test set32.coordinates[1, :] ≈ [0.40002938831494556, -0.0010525966270966628, 0.9165017078678638] atol = 0.0 rtol = 1.0e-14
end

@testset "Explicit angular Fibonacci and optimization paths" begin
    fib10_a = fibonacci_sphere_point_set(10)
    fib10_b = fibonacci_sphere_point_set(10)
    fib15 = fibonacci_sphere_point_set(15)
    optimized10 = optimize_sphere_point_set(fib10_a; beta = 2.0, iters = 20, gtol = 1.0e-8)

    @test fib10_a isa SpherePointSet
    @test fib15 isa SpherePointSet
    @test fib10_a.order == fib10_a.cardinality == 10
    @test fib15.order == fib15.cardinality == 15
    @test fib10_a.coordinates ≈ fib10_b.coordinates atol = 0.0 rtol = 0.0
    @test fib10_a.nn_ratio == fib10_b.nn_ratio
    @test fib10_a.provenance.source_tag == "deterministic_fibonacci_seed"
    @test fib10_a.provenance.source_project == "GaussletBases"
    @test fib10_a.provenance.source_artifact == "in_memory_generation"
    @test occursin("Fibonacci sphere seed", fib10_a.provenance.source_note)
    @test occursin("no randomization", fib10_a.provenance.source_note)

    fib10_kappa = GaussletBases._sphere_point_kappa_from_beta(fib10_a.order, 2.0)
    fib10_initial_logdet =
        GaussletBases._sphere_point_logdet(
            GaussletBases._sphere_point_gaussian_gram(fib10_a.coordinates, fib10_kappa),
        )
    fib10_optimized_logdet =
        GaussletBases._sphere_point_logdet(
            GaussletBases._sphere_point_gaussian_gram(optimized10.coordinates, fib10_kappa),
        )

    @test optimized10 isa SpherePointSet
    @test optimized10.order == fib10_a.order
    @test optimized10.cardinality == fib10_a.cardinality
    @test optimized10.provenance.source_tag == "optimized_from_input_point_set"
    @test optimized10.provenance.source_project == "GaussletBases"
    @test optimized10.provenance.source_artifact == "in_memory_optimization"
    @test occursin("deterministic_fibonacci_seed", optimized10.provenance.source_note)
    @test occursin("beta=2.0", optimized10.provenance.source_note)
    @test occursin("iters=20", optimized10.provenance.source_note)
    @test occursin("gtol=1.0e-8", optimized10.provenance.source_note)
    @test maximum(abs.(sqrt.(sum(abs2, optimized10.coordinates; dims = 2)) .- 1.0)) < 1.0e-12
    @test optimized10.nn_ratio >= 1.0
    @test fib10_optimized_logdet ≥ fib10_initial_logdet - 1.0e-10
end

@testset "Shell-local injected angular basis" begin
    shell15 = _shell_local_injected_angular_fixture(15)
    shell32 = _shell_local_injected_angular_fixture(32)
    shell51 = _shell_local_injected_angular_fixture(51)

    @test shell15 isa ShellLocalInjectedAngularBasis
    @test shell32 isa ShellLocalInjectedAngularBasis
    @test shell51 isa ShellLocalInjectedAngularBasis

    @test shell15.l_inject == 1
    @test shell32.l_inject == 3
    @test shell51.l_inject == 4

    for shell in (shell15, shell32, shell51)
        diagnostics = shell_local_injected_angular_diagnostics(shell)
        expected_kinetic = Diagonal(diagnostics.expected_injected_kinetic_eigenvalues)

        @test shell.prototype_count == shell.point_set.order
        @test shell.final_count == shell.prototype_count
        @test shell.injected_count == (shell.l_inject + 1)^2
        @test shell.whitened_complement_count ≥ shell.final_count - shell.injected_count
        @test size(shell.prototype_overlap) == (shell.prototype_count, shell.prototype_count)
        @test size(shell.final_overlap) == (shell.final_count, shell.final_count)
        @test size(shell.final_kinetic) == (shell.final_count, shell.final_count)
        @test size(shell.injected_overlap) == (shell.injected_count, shell.final_count)
        @test sum(shell.shell_weights) ≈ 4 * pi atol = 1.0e-12 rtol = 1.0e-12
        @test minimum(shell.theta_nn) > 0.0
        @test minimum(shell.kappa) > 0.0
        @test diagnostics.overlap_error ≤ 1.0e-10
        @test diagnostics.injected_exactness_error ≤ 1.0e-10
        @test diagnostics.injected_kinetic_error ≤ 1.0e-9
        @test norm(shell.injected_kinetic - expected_kinetic, Inf) ≤ 1.0e-9
        @test issorted(diagnostics.expected_injected_kinetic_eigenvalues)
        @test sort(diagnostics.injected_kinetic_eigenvalues) ≈
              diagnostics.expected_injected_kinetic_eigenvalues atol = 1.0e-9 rtol = 1.0e-9
    end
end

@testset "Shell-local angular profiles" begin
    profile15 = _shell_local_angular_profile_fixture(15)
    profile15_again = shell_local_angular_profile(15)
    profile15_uncached =
        GaussletBases._build_shell_local_angular_profile_uncached(sphere_point_set(15))
    profile32 = _shell_local_angular_profile_fixture(32)
    overlap15_32 = _shell_local_angular_profile_overlap_fixture(15, 32)
    overlap15_32_again = adjacent_shell_local_angular_profile_overlap(profile15, profile32)

    @test profile15 isa ShellLocalAngularProfile
    @test profile15.key.order == 15
    @test profile15.key.point_set_source_tag == sphere_point_set(15).provenance.source_tag
    @test profile15.key.gauge_version == :v1_seed_order_dominant_positive
    @test profile15 === profile15_again
    @test build_shell_local_injected_angular_basis(15) === profile15.basis
    @test profile15.profile_id == profile15_uncached.profile_id
    @test profile15.labels == profile15_uncached.labels
    @test profile15.block_kinds == profile15_uncached.block_kinds
    @test profile15.diagnostics.grand_coefficients_checksum ==
          profile15_uncached.diagnostics.grand_coefficients_checksum

    expected_exact_labels = [
        "exact_l$(channel.l)_m$(channel.m)" for channel in profile15.basis.injected_channels.channel_data
    ]
    @test profile15.exact_labels == expected_exact_labels
    @test profile15.labels[1:profile15.basis.injected_count] == expected_exact_labels
    @test profile15.mixed_labels == ["mixed_$(i)" for i in 1:(profile15.basis.final_count - profile15.basis.injected_count)]
    @test profile15.block_kinds ==
          vcat(fill(:exact, profile15.basis.injected_count), fill(:mixed, profile15.basis.final_count - profile15.basis.injected_count))
    @test all(diag(profile15.basis.injected_overlap[:, 1:profile15.basis.injected_count]) .> 0.0)

    mixed_offset = profile15.basis.injected_count
    mixed_metadata = profile15.gauge_metadata
    @test mixed_metadata.mixed_orientation_strategy == :seed_order
    @test length(mixed_metadata.mixed_dominant_grand_indices) == profile15.basis.final_count - profile15.basis.injected_count
    @test length(mixed_metadata.mixed_signs) == profile15.basis.final_count - profile15.basis.injected_count
    for j in 1:(profile15.basis.final_count - profile15.basis.injected_count)
        idx = mixed_metadata.mixed_dominant_grand_indices[j]
        @test profile15.basis.grand_coefficients[idx, mixed_offset + j] > 0.0
    end

    @test overlap15_32 isa ShellLocalAngularProfileOverlap
    @test overlap15_32 === overlap15_32_again
    @test size(overlap15_32.overlap) == (profile15.basis.final_count, profile32.basis.final_count)
    @test overlap15_32.source_labels == profile15.labels
    @test overlap15_32.target_labels == profile32.labels
    @test overlap15_32.source_exact_count == profile15.basis.injected_count
    @test overlap15_32.target_exact_count == profile32.basis.injected_count
    @test overlap15_32.shell_independent
    @test overlap15_32.diagnostics.min_singular_value > 1.0e-8
    @test isfinite(overlap15_32.diagnostics.exact_block_inf_norm)
end

@testset "Atomic shell-local angular assembly" begin
    assembly = _atomic_shell_local_angular_fixture()
    diagnostics = atomic_shell_local_angular_diagnostics(assembly)
    scheduled_orders = assign_atomic_angular_shell_orders(
        [0.05, 0.5, 1.2, 3.0, 8.0];
        ord_min = 15,
        ord_max = 51,
        r_lo = 0.2,
        r_hi = 4.5,
        w_lo = 0.2,
        w_hi = 0.7,
    )
    scheduled_orders_full = assign_atomic_angular_shell_orders(
        [0.05, 0.5, 1.2, 3.0, 8.0];
        ord_min = 24,
        ord_max = 58,
        r_lo = 0.2,
        r_hi = 4.5,
        w_lo = 0.2,
        w_hi = 0.7,
    )

    @test assembly isa AtomicShellLocalInjectedAngularAssembly
    @test assembly.shell_orders == [15, 32, 51, 32]
    @test assembly.shell_dimensions == [15, 32, 51, 32]
    @test assembly.shell_offsets == [1, 16, 48, 99]
    @test assembly.shell_exact_lmax == [1, 3, 4, 3]
    @test assembly.profiles[1] === shell_local_angular_profile(15)
    @test assembly.profiles[2] === shell_local_angular_profile(32)
    @test assembly.profiles[2] === assembly.profiles[4]
    @test assembly.shells[2] === assembly.profiles[2].basis
    @test assembly.shells[2] === assembly.shells[4]
    @test size(assembly.overlap) == (130, 130)
    @test size(assembly.kinetic) == (130, 130)
    @test diagnostics.nshells == 4
    @test diagnostics.total_dim == 130
    @test diagnostics.max_shell_overlap_error ≤ 1.0e-10
    @test diagnostics.max_shell_injected_exactness_error ≤ 1.0e-10
    @test diagnostics.max_shell_injected_kinetic_error ≤ 1.0e-9
    @test diagnostics.max_diagonal_overlap_error ≤ 1.0e-10
    @test diagnostics.max_diagonal_injected_kinetic_error ≤ 1.0e-9
    @test diagnostics.max_pair_overlap_symmetry_error ≤ 1.0e-10
    @test diagnostics.max_pair_kinetic_symmetry_error ≤ 1.0e-10
    @test all(
        diagnostics.shell_interaction_lcap[i] ≥
        diagnostics.shell_interaction_lexpand[i] ≥
        diagnostics.shell_exact_lmax[i] for i in eachindex(diagnostics.shell_exact_lmax)
    )
    @test diagnostics.interaction_pair_plan ≥ diagnostics.interaction_exact_lower_bound
    @test diagnostics.max_shell_interaction_tail ≤ 1.0e-6

    @test scheduled_orders[1] == 15
    @test scheduled_orders[end] == 15
    @test all(order in sphere_point_set_orders() for order in scheduled_orders)
    @test any(order ∉ curated_sphere_point_set_orders() for order in scheduled_orders)
    @test maximum(scheduled_orders) ≥ 32
    @test all(order in sphere_point_set_orders() for order in scheduled_orders_full)
    @test any(order ∉ curated_sphere_point_set_orders() for order in scheduled_orders_full)
end

@testset "Atomic fixed-radial angular sequence" begin
    sequence = _paper_style_fixed_radial_angular_sequence_fixture()
    level10 = sequence.levels[1]
    level15 = sequence.levels[2]
    level32 = sequence.levels[3]
    sidecar10_15 = sequence.adjacent_overlaps[1]
    sidecar15_32 = sequence.adjacent_overlaps[2]
    sidecar10_32 = sequence.direct_overlaps[1]
    nr = length(sequence.shell_ids)

    @test sequence isa AtomicFixedRadialAngularSequence
    @test sequence.N_sph_values == [10, 15, 32]
    @test length(sequence.levels) == 3
    @test length(sequence.adjacent_overlaps) == 2
    @test length(sequence.direct_overlaps) == 1
    @test sequence.shell_ids == collect(1:nr)
    @test issorted(sequence.shell_centers_r)
    @test all(level -> level.radial_basis_id == sequence.radial_basis_id, sequence.levels)
    @test all(level -> level.shell_ids == sequence.shell_ids, sequence.levels)
    @test all(level -> level.shell_centers_r == sequence.shell_centers_r, sequence.levels)
    @test level10.profile === shell_local_angular_profile(10)
    @test level15.profile === shell_local_angular_profile(15)
    @test level32.profile === shell_local_angular_profile(32)
    @test all(profile -> profile === level10.profile, level10.payload.one_body.angular_assembly.profiles)
    @test all(profile -> profile === level15.profile, level15.payload.one_body.angular_assembly.profiles)
    @test all(profile -> profile === level32.profile, level32.payload.one_body.angular_assembly.profiles)
    @test level10.shell_dimensions == fill(level10.profile.basis.final_count, nr)
    @test level15.shell_dimensions == fill(level15.profile.basis.final_count, nr)
    @test level32.shell_dimensions == fill(level32.profile.basis.final_count, nr)

    level_payload = atomic_fixed_radial_angular_level_dense_payload(level10)
    @test level_payload.payload["H1"] ≈ level10.payload.hamiltonian atol = 0.0 rtol = 0.0
    @test level_payload.payload["Vee"] ≈ level10.payload.interaction atol = 0.0 rtol = 0.0
    @test level_payload.payload["shell_ids"] == sequence.shell_ids
    @test level_payload.payload["shell_centers_r"] == sequence.shell_centers_r
    @test level_payload.payload["shell_dimensions"] == level10.shell_dimensions
    @test level_payload.payload["within_shell_labels"] == level10.profile.labels
    @test level_payload.bridge_meta["sequence_id"] == sequence.sequence_id
    @test level_payload.bridge_meta["N_sph"] == 10
    @test level_payload.bridge_meta["angular_profile_id"] == level10.profile.profile_id
    @test level_payload.bridge_meta["gauge_version"] == string(level10.profile.key.gauge_version)

    overlap_payload = atomic_fixed_radial_angular_overlap_sidecar_payload(sidecar10_15)
    @test overlap_payload.payload["overlap"] ≈ sidecar10_15.overlap atol = 0.0 rtol = 0.0
    @test overlap_payload.payload["source_labels"] == level10.profile.labels
    @test overlap_payload.payload["target_labels"] == level15.profile.labels
    @test overlap_payload.bridge_meta["sequence_id"] == sequence.sequence_id
    @test overlap_payload.bridge_meta["source_N_sph"] == 10
    @test overlap_payload.bridge_meta["target_N_sph"] == 15
    @test overlap_payload.bridge_meta["pair_kind"] == "adjacent"
    @test overlap_payload.bridge_meta["shell_independent"]
    @test sidecar10_15.source_profile_id == level10.profile.profile_id
    @test sidecar10_15.target_profile_id == level15.profile.profile_id
    @test sidecar10_15.source_gauge_version == level10.profile.key.gauge_version
    @test sidecar10_15.target_gauge_version == level15.profile.key.gauge_version
    @test sidecar10_15.source_labels == level10.profile.labels
    @test sidecar10_15.target_labels == level15.profile.labels
    @test sidecar10_15.shell_independent
    @test sidecar10_15.pair_kind == :adjacent
    @test sidecar15_32.source_profile_id == level15.profile.profile_id
    @test sidecar15_32.target_profile_id == level32.profile.profile_id
    @test sidecar15_32.pair_kind == :adjacent
    @test sidecar10_32.source_N_sph == 10
    @test sidecar10_32.target_N_sph == 32
    @test sidecar10_32.source_profile_id == level10.profile.profile_id
    @test sidecar10_32.target_profile_id == level32.profile.profile_id
    @test sidecar10_32.source_gauge_version == level10.profile.key.gauge_version
    @test sidecar10_32.target_gauge_version == level32.profile.key.gauge_version
    @test sidecar10_32.source_labels == level10.profile.labels
    @test sidecar10_32.target_labels == level32.profile.labels
    @test sidecar10_32.shell_independent
    @test sidecar10_32.pair_kind == :direct

    legacy_payload = atomic_fixed_radial_legacy_dmrgatom_payload(level10)
    @test legacy_payload.payload["H1"] ≈ level10.payload.hamiltonian atol = 0.0 rtol = 0.0
    @test legacy_payload.payload["Vee"] ≈ level10.payload.interaction atol = 0.0 rtol = 0.0
    @test legacy_payload.payload["dims_per_shell"] == level10.shell_dimensions
    @test legacy_payload.meta_values["Z"] == 2
    @test legacy_payload.bridge_meta["angular_profile_id"] == level10.profile.profile_id
    @test legacy_payload.bridge_meta["gauge_version"] == string(level10.profile.key.gauge_version)
    @test legacy_payload.bridge_meta["basis_centers_are_representative"]
    @test legacy_payload.bridge_meta["basis_centers_center_policy_name"] ==
          "representative_dominant_prototype_direction_v1"
    dominant_indices =
        Int.(legacy_payload.bridge_meta["within_shell_dominant_prototype_indices"])
    @test length(dominant_indices) == level10.profile.basis.final_count
    @test all(index -> 1 <= index <= level10.profile.basis.prototype_count, dominant_indices)
    basis_centers = Matrix{Float64}(legacy_payload.payload["basis_centers"])
    @test size(basis_centers) ==
          (size(level10.payload.hamiltonian, 1), 3)
    point_coordinates = level10.profile.basis.point_set.coordinates
    shell_size = level10.profile.basis.final_count
    first_shell_centers = basis_centers[1:shell_size, :]
    expected_first_shell_centers =
        level10.shell_centers_r[1] .* point_coordinates[dominant_indices, :]
    @test first_shell_centers ≈ expected_first_shell_centers atol = 1.0e-12 rtol = 1.0e-12

    mktempdir() do dir
        level_path = joinpath(dir, "he_fixed_radial_level10.jld2")
        sidecar_path = joinpath(dir, "he_fixed_radial_10_15_overlap.jld2")
        direct_sidecar_path = joinpath(dir, "he_fixed_radial_10_32_overlap.jld2")
        legacy_path = joinpath(dir, "he_fixed_radial_level10.legacy_dmrgatom.jld2")
        @test write_atomic_fixed_radial_angular_level_jld2(level_path, level10) == level_path
        @test write_atomic_fixed_radial_angular_overlap_sidecar_jld2(sidecar_path, sidecar10_15) == sidecar_path
        @test write_atomic_fixed_radial_angular_overlap_sidecar_jld2(direct_sidecar_path, sidecar10_32) == direct_sidecar_path
        @test write_atomic_fixed_radial_legacy_dmrgatom_jld2(legacy_path, level10) == legacy_path
        jldopen(level_path, "r") do file
            @test String(file["bridge/format"]) == "angular_fixed_radial_dense_v1"
            @test Int(file["bridge/level_index"]) == 1
            @test Int(file["bridge/N_sph"]) == 10
            @test String(file["bridge/sequence_id"]) == sequence.sequence_id
            @test String(file["bridge/angular_profile_id"]) == level10.profile.profile_id
            @test Int.(file["shell_ids"]) == sequence.shell_ids
            @test Float64.(file["shell_centers_r"]) == sequence.shell_centers_r
            @test String.(file["within_shell_labels"]) == level10.profile.labels
        end
        jldopen(sidecar_path, "r") do file
            @test String(file["bridge/format"]) == "angular_fixed_radial_profile_overlap_v1"
            @test Int(file["bridge/source_N_sph"]) == 10
            @test Int(file["bridge/target_N_sph"]) == 15
            @test String(file["bridge/source_level_id"]) == level10.level_id
            @test String(file["bridge/target_level_id"]) == level15.level_id
            @test String(file["bridge/source_gauge_version"]) == string(level10.profile.key.gauge_version)
            @test String(file["bridge/target_gauge_version"]) == string(level15.profile.key.gauge_version)
            @test String(file["bridge/pair_kind"]) == "adjacent"
            @test Bool(file["bridge/shell_independent"])
            @test String.(file["source_labels"]) == level10.profile.labels
            @test String.(file["target_labels"]) == level15.profile.labels
            @test Matrix{Float64}(file["overlap"]) ≈ sidecar10_15.overlap atol = 0.0 rtol = 0.0
        end
        jldopen(direct_sidecar_path, "r") do file
            @test String(file["bridge/format"]) == "angular_fixed_radial_profile_overlap_v1"
            @test Int(file["bridge/source_N_sph"]) == 10
            @test Int(file["bridge/target_N_sph"]) == 32
            @test String(file["bridge/source_level_id"]) == level10.level_id
            @test String(file["bridge/target_level_id"]) == level32.level_id
            @test String(file["bridge/source_profile_id"]) == level10.profile.profile_id
            @test String(file["bridge/target_profile_id"]) == level32.profile.profile_id
            @test String(file["bridge/source_gauge_version"]) == string(level10.profile.key.gauge_version)
            @test String(file["bridge/target_gauge_version"]) == string(level32.profile.key.gauge_version)
            @test String(file["bridge/pair_kind"]) == "direct"
            @test Bool(file["bridge/shell_independent"])
            @test String.(file["source_labels"]) == level10.profile.labels
            @test String.(file["target_labels"]) == level32.profile.labels
            @test Matrix{Float64}(file["overlap"]) ≈ sidecar10_32.overlap atol = 0.0 rtol = 0.0
        end
        jldopen(legacy_path, "r") do file
            @test String(file["bridge/format"]) == "legacy_dmrgatom_dense_v1"
            @test Int(file["bridge/N_sph"]) == 10
            @test Int(file["meta/Z"]) == 2
            @test Int.(file["dims_per_shell"]) == level10.shell_dimensions
            @test Matrix{Float64}(file["H1"]) ≈ level10.payload.hamiltonian atol = 0.0 rtol = 0.0
            @test Matrix{Float64}(file["Vee"]) ≈ level10.payload.interaction atol = 0.0 rtol = 0.0
            @test String(file["bridge/basis_centers_center_policy_name"]) ==
                  "representative_dominant_prototype_direction_v1"
            @test Bool(file["bridge/basis_centers_are_representative"])
            @test Int.(file["bridge/within_shell_dominant_prototype_indices"]) ==
                  dominant_indices
            @test Matrix{Float64}(file["basis_centers"]) ≈ basis_centers atol = 0.0 rtol = 0.0
        end
    end
end

@testset "Atomic injected angular one-body benchmark" begin
    benchmark = _atomic_injected_angular_one_body_benchmark_fixture()
    diagnostics = atomic_injected_angular_one_body_diagnostics(benchmark)
    _, _, radial_ops, _, _, _ = _quick_radial_atomic_fixture()
    exact_atom = atomic_one_body_operators(radial_ops; lmax = diagnostics.exact_common_lmax)
    exact_benchmark_spectrum =
        GaussletBases._generalized_spectrum(benchmark.exact_hamiltonian, benchmark.exact_overlap)
    exact_atom_spectrum =
        GaussletBases._generalized_spectrum(exact_atom.hamiltonian, exact_atom.overlap)

    @test benchmark isa AtomicInjectedAngularOneBodyBenchmark
    @test benchmark.angular_assembly isa AtomicShellLocalInjectedAngularAssembly
    @test benchmark.angular_assembly.shell_radii == radial_ops.shell_centers_r
    @test benchmark.exact_common_lmax ≥ 1
    @test length(benchmark.exact_channels) == (benchmark.exact_common_lmax + 1)^2
    @test size(benchmark.overlap) == size(benchmark.hamiltonian)
    @test size(benchmark.overlap, 1) == sum(benchmark.angular_assembly.shell_dimensions)
    @test size(benchmark.exact_overlap) == size(exact_atom.overlap)
    @test size(benchmark.exact_hamiltonian) == size(exact_atom.hamiltonian)
    @test exact_benchmark_spectrum ≈ exact_atom_spectrum atol = 1.0e-10 rtol = 1.0e-10

    @test diagnostics.overlap_symmetry_error ≤ 1.0e-10
    @test diagnostics.hamiltonian_symmetry_error ≤ 1.0e-10
    @test diagnostics.min_overlap_eigenvalue > 1.0e-6
    @test diagnostics.projected_exact_overlap_error ≤ 1.0e-9
    @test diagnostics.projected_exact_hamiltonian_error ≤ 1.0e-8
    @test diagnostics.projected_exact_low_eigenvalue_count == 4
    @test diagnostics.projected_exact_low_eigenvalue_error ≤ 1.0e-8
    @test diagnostics.benchmark_ground_state_error ≤ 1.0e-8
    @test diagnostics.benchmark_ground_state_energy ≈ diagnostics.exact_ground_state_energy atol = 1.0e-8 rtol = 1.0e-8
end

@testset "Atomic injected angular Cartesian moments" begin
    rb, grid, radial_ops, benchmark, bundle, bundle_reconstructed =
        _paper_style_angular_cartesian_moments_fixture(15; Z = 2.0, lmax = 2)

    @test bundle isa AtomicInjectedAngularCartesianMomentBundle
    @test bundle.radial_operators === radial_ops
    @test bundle.angular_assembly === benchmark.angular_assembly
    @test bundle.shell_ranges == benchmark.shell_ranges
    @test bundle.S ≈ benchmark.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test bundle.radial_moment_source == :explicit_basis_and_grid
    @test bundle_reconstructed.radial_moment_source == :reconstructed_from_source_manifest
    @test bundle_reconstructed.S ≈ bundle.S atol = 0.0 rtol = 0.0
    @test bundle_reconstructed.X ≈ bundle.X atol = 1.0e-10 rtol = 1.0e-10
    @test bundle_reconstructed.Y ≈ bundle.Y atol = 1.0e-10 rtol = 1.0e-10
    @test bundle_reconstructed.Z ≈ bundle.Z atol = 1.0e-10 rtol = 1.0e-10
    @test bundle_reconstructed.X2 ≈ bundle.X2 atol = 1.0e-10 rtol = 1.0e-10
    @test bundle_reconstructed.Y2 ≈ bundle.Y2 atol = 1.0e-10 rtol = 1.0e-10
    @test bundle_reconstructed.Z2 ≈ bundle.Z2 atol = 1.0e-10 rtol = 1.0e-10
    @test bundle_reconstructed.R2 ≈ bundle.R2 atol = 1.0e-10 rtol = 1.0e-10

    for matrix in (bundle.S, bundle.X, bundle.Y, bundle.Z, bundle.X2, bundle.Y2, bundle.Z2, bundle.R2)
        @test all(isfinite, matrix)
        @test opnorm(matrix - transpose(matrix), Inf) ≤ 1.0e-8
    end

    @test opnorm(bundle.X, Inf) > 1.0e-8
    @test opnorm(bundle.Y, Inf) > 1.0e-8
    @test opnorm(bundle.Z, Inf) > 1.0e-8
    @test opnorm(bundle.X2, Inf) > 1.0e-8
    @test opnorm(bundle.Y2, Inf) > 1.0e-8
    @test opnorm(bundle.Z2, Inf) > 1.0e-8
    @test bundle.R2 ≈ bundle.X2 + bundle.Y2 + bundle.Z2 atol = 1.0e-8 rtol = 1.0e-8

    s_exact = build_atomic_injected_angular_one_body_benchmark(
        radial_ops;
        shell_orders = fill(10, length(radial_ops.shell_centers_r)),
        l_inject = 0,
    )
    s_bundle = build_atomic_injected_angular_cartesian_moments(
        s_exact;
        radial_basis = rb,
        radial_grid = grid,
    )
    exact_transform = s_exact.exact_transform
    projected_X = exact_transform * s_bundle.X * transpose(exact_transform)
    projected_Y = exact_transform * s_bundle.Y * transpose(exact_transform)
    projected_Z = exact_transform * s_bundle.Z * transpose(exact_transform)
    projected_X2 = exact_transform * s_bundle.X2 * transpose(exact_transform)
    projected_Y2 = exact_transform * s_bundle.Y2 * transpose(exact_transform)
    projected_Z2 = exact_transform * s_bundle.Z2 * transpose(exact_transform)
    projected_R2 = exact_transform * s_bundle.R2 * transpose(exact_transform)

    @test s_exact.exact_common_lmax == 0
    @test projected_X ≈ zeros(size(projected_X)) atol = 1.0e-10 rtol = 0.0
    @test projected_Y ≈ zeros(size(projected_Y)) atol = 1.0e-10 rtol = 0.0
    @test projected_Z ≈ zeros(size(projected_Z)) atol = 1.0e-10 rtol = 0.0
    @test projected_X2 ≈ projected_Y2 atol = 1.0e-9 rtol = 1.0e-9
    @test projected_Y2 ≈ projected_Z2 atol = 1.0e-9 rtol = 1.0e-9
    @test projected_R2 ≈ projected_X2 + projected_Y2 + projected_Z2 atol = 1.0e-9 rtol = 1.0e-9
end

@testset "Atomic injected angular HF-style benchmark" begin
    benchmark = _atomic_injected_angular_hf_style_benchmark_fixture()
    diagnostics = atomic_injected_angular_hf_style_diagnostics(benchmark)
    exact_ida = atomic_ida_operators(benchmark.one_body.radial_operators; lmax = benchmark.one_body.exact_common_lmax)
    exact_shell_major_interaction = atomic_ida_density_interaction_matrix(exact_ida; ordering = :shell_major)

    @test benchmark isa AtomicInjectedAngularHFStyleBenchmark
    @test benchmark.one_body isa AtomicInjectedAngularOneBodyBenchmark
    @test benchmark.exact_ida_reference isa AtomicIDAOperators
    @test size(benchmark.interaction) == size(benchmark.one_body.overlap)
    @test size(benchmark.exact_interaction) == size(benchmark.one_body.exact_overlap)
    @test benchmark.exact_interaction ≈ exact_shell_major_interaction atol = 1.0e-12 rtol = 1.0e-12

    @test diagnostics.interaction_symmetry_error ≤ 1.0e-10
    @test diagnostics.full_converged
    @test diagnostics.exact_converged
    @test diagnostics.full_residual ≤ 1.0e-8
    @test diagnostics.exact_residual ≤ 1.0e-8 || diagnostics.exact_iterations == 13
    @test diagnostics.full_electron_count_error ≤ 1.0e-10
    @test diagnostics.exact_electron_count_error ≤ 1.0e-10
    @test isfinite(diagnostics.full_energy)
    @test isfinite(diagnostics.exact_energy)
    @test abs(diagnostics.energy_difference_to_exact_reference) ≤ 1.0e-8
    @test diagnostics.ground_orbital_energy_error ≤ 1.0e-8

    hfdmrg = _local_hfdmrg_module()
    if hfdmrg !== nothing
        payload = build_atomic_injected_angular_hfdmrg_payload(benchmark)
        @test payload isa AtomicInjectedAngularHFDMRGHFAdapter
        @test payload.route == :dense_density_density
        @test payload.solver_mode == :restricted_closed_shell
        @test payload.nup == 1
        @test payload.ndn == 1
        hfdmrg_result = run_atomic_injected_angular_hfdmrg_hf(
            payload;
            hfmod = hfdmrg,
            nblockcenter = 1,
            maxiter = 40,
            cutoff = 1.0e-10,
            scf_cutoff = 1.0e-11,
            verbose = false,
        )
        @test hfdmrg_result.route == payload.route
        @test hfdmrg_result.solver_mode == :restricted_closed_shell
        @test hfdmrg_result.nblockcenter == 1
        @test hfdmrg_result.blocksize == min(size(payload.hamiltonian, 1), 64)
        @test abs(diagnostics.full_energy - hfdmrg_result.energy) ≤ 1.0e-8
    end
end

@testset "Atomic injected angular HFDMRG-facing HF adapter" begin
    adapter = _atomic_injected_angular_hfdmrg_hf_adapter_fixture()
    diagnostics = atomic_injected_angular_hfdmrg_hf_adapter_diagnostics(adapter)
    benchmark = _atomic_injected_angular_hf_style_benchmark_fixture()
    from_benchmark = build_atomic_injected_angular_hfdmrg_hf_adapter(benchmark)
    payload_from_benchmark = build_atomic_injected_angular_hfdmrg_payload(benchmark)
    from_small_ed =
        build_atomic_injected_angular_hfdmrg_hf_adapter(_atomic_injected_angular_small_ed_benchmark_fixture())
    benchmark_diagnostics = atomic_injected_angular_hfdmrg_hf_adapter_diagnostics(from_benchmark)
    open_shell_seeds =
        build_atomic_injected_angular_hfdmrg_hf_seeds(benchmark; nup = 2, ndn = 1)
    explicit_psiup0 = open_shell_seeds.psiup0[:, [2, 1]]
    explicit_psidn0 = open_shell_seeds.psidn0
    open_shell_adapter = build_atomic_injected_angular_hfdmrg_hf_adapter(
        benchmark;
        nup = 2,
        ndn = 1,
    )
    explicit_open_shell_adapter = build_atomic_injected_angular_hfdmrg_hf_adapter(
        benchmark;
        nup = 2,
        ndn = 1,
        psiup0 = explicit_psiup0,
        psidn0 = explicit_psidn0,
    )
    mixed_seed_adapter = build_atomic_injected_angular_hfdmrg_hf_adapter(
        benchmark;
        nup = 2,
        ndn = 1,
        psiup0 = explicit_psiup0,
    )

    @test adapter isa AtomicInjectedAngularHFDMRGHFAdapter
    @test adapter.one_body isa AtomicInjectedAngularOneBodyBenchmark
    @test isnothing(adapter.hf_style)
    @test adapter.route == :dense_density_density
    @test size(adapter.hamiltonian) == size(adapter.interaction)
    @test size(adapter.hamiltonian, 1) == diagnostics.basis_dim
    @test size(adapter.psiup0) == (diagnostics.basis_dim, diagnostics.nup)
    @test size(adapter.psidn0) == (diagnostics.basis_dim, diagnostics.ndn)
    @test diagnostics.nup == 1
    @test diagnostics.ndn == 1
    @test diagnostics.solver_mode == :restricted_closed_shell
    @test diagnostics.overlap_identity_error ≤ 2.0e-6
    @test diagnostics.hamiltonian_symmetry_error ≤ 1.0e-8
    @test diagnostics.interaction_symmetry_error ≤ 1.0e-10
    @test diagnostics.psiup0_orthogonality_error ≤ 1.0e-12
    @test diagnostics.psidn0_orthogonality_error ≤ 1.0e-12
    @test diagnostics.psiup0_source == :default_one_body_orbitals
    @test diagnostics.psidn0_source == :default_one_body_orbitals
    @test !diagnostics.has_benchmark_reference
    @test ismissing(diagnostics.benchmark_full_energy)
    @test ismissing(diagnostics.benchmark_exact_energy)
    @test from_benchmark.hf_style === benchmark
    @test payload_from_benchmark.hf_style === benchmark
    @test benchmark_diagnostics.has_benchmark_reference
    @test benchmark_diagnostics.solver_mode == :restricted_closed_shell
    @test isfinite(benchmark_diagnostics.benchmark_full_energy)
    @test isfinite(benchmark_diagnostics.benchmark_exact_energy)
    @test from_benchmark.route == adapter.route
    @test from_benchmark.nup == adapter.nup
    @test from_benchmark.ndn == adapter.ndn
    @test from_benchmark.hamiltonian ≈ adapter.hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test from_benchmark.interaction ≈ adapter.interaction atol = 1.0e-12 rtol = 1.0e-12
    @test from_benchmark.psiup0 ≈ adapter.psiup0 atol = 1.0e-12 rtol = 1.0e-12
    @test from_benchmark.psidn0 ≈ adapter.psidn0 atol = 1.0e-12 rtol = 1.0e-12
    @test payload_from_benchmark.route == from_benchmark.route
    @test payload_from_benchmark.solver_mode == from_benchmark.solver_mode
    @test payload_from_benchmark.hamiltonian ≈ from_benchmark.hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test payload_from_benchmark.interaction ≈ from_benchmark.interaction atol = 1.0e-12 rtol = 1.0e-12
    @test payload_from_benchmark.psiup0 ≈ from_benchmark.psiup0 atol = 1.0e-12 rtol = 1.0e-12
    @test payload_from_benchmark.psidn0 ≈ from_benchmark.psidn0 atol = 1.0e-12 rtol = 1.0e-12
    @test from_small_ed.route == adapter.route
    @test from_small_ed.nup == adapter.nup
    @test from_small_ed.ndn == adapter.ndn
    @test from_small_ed.hamiltonian ≈ adapter.hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test from_small_ed.interaction ≈ adapter.interaction atol = 1.0e-12 rtol = 1.0e-12
    @test from_small_ed.psiup0 ≈ adapter.psiup0 atol = 1.0e-12 rtol = 1.0e-12
    @test from_small_ed.psidn0 ≈ adapter.psidn0 atol = 1.0e-12 rtol = 1.0e-12
    @test size(open_shell_seeds.psiup0) == (diagnostics.basis_dim, 2)
    @test size(open_shell_seeds.psidn0) == (diagnostics.basis_dim, 1)
    @test open_shell_adapter.nup == 2
    @test open_shell_adapter.ndn == 1
    @test open_shell_adapter.solver_mode == :unrestricted
    @test open_shell_adapter.psiup0_source == :default_one_body_orbitals
    @test open_shell_adapter.psidn0_source == :default_one_body_orbitals
    @test opnorm(transpose(open_shell_adapter.psiup0) * open_shell_adapter.psiup0 - Matrix{Float64}(I, 2, 2), Inf) ≤ 1.0e-12
    @test opnorm(transpose(open_shell_adapter.psidn0) * open_shell_adapter.psidn0 - Matrix{Float64}(I, 1, 1), Inf) ≤ 1.0e-12
    @test explicit_open_shell_adapter.solver_mode == :unrestricted
    @test explicit_open_shell_adapter.psiup0_source == :explicit_seed
    @test explicit_open_shell_adapter.psidn0_source == :explicit_seed
    @test explicit_open_shell_adapter.psiup0 ≈ explicit_psiup0 atol = 1.0e-12 rtol = 1.0e-12
    @test explicit_open_shell_adapter.psidn0 ≈ explicit_psidn0 atol = 1.0e-12 rtol = 1.0e-12
    @test mixed_seed_adapter.solver_mode == :unrestricted
    @test mixed_seed_adapter.psiup0_source == :explicit_seed
    @test mixed_seed_adapter.psidn0_source == :default_one_body_orbitals
    @test mixed_seed_adapter.psiup0 ≈ explicit_psiup0 atol = 1.0e-12 rtol = 1.0e-12
    @test mixed_seed_adapter.psidn0 ≈ open_shell_seeds.psidn0 atol = 1.0e-12 rtol = 1.0e-12
end

@testset "Angular He legacy-trim HF payload anchors" begin
    hfdmrg = _local_hfdmrg_module()
    hfdmrg === nothing && return

    he_anchor_reference = -2.861679990485
    anchor_settings = (
        nblockcenter = 2,
        blocksize = 100,
        maxiter = 100,
        cutoff = 1.0e-8,
        scf_cutoff = 1.0e-9,
        verbose = false,
    )

    payload10 = _paper_style_angular_hfdmrg_payload_fixture(10)
    payload15 = _paper_style_angular_hfdmrg_payload_fixture(15)
    result10 = _solve_hfdmrg_from_payload_direct(payload10, hfdmrg; anchor_settings...)
    result15 = _solve_hfdmrg_from_payload_direct(payload15, hfdmrg; anchor_settings...)

    @test payload10.solver_mode == :restricted_closed_shell
    @test payload15.solver_mode == :restricted_closed_shell
    @test payload10.nup == 1
    @test payload10.ndn == 1
    @test payload15.nup == 1
    @test payload15.ndn == 1
    @test first(payload10.one_body.angular_assembly.shell_orders) == 10
    @test first(payload15.one_body.angular_assembly.shell_orders) == 15
    @test all(==(10), payload10.one_body.angular_assembly.shell_orders)
    @test all(==(15), payload15.one_body.angular_assembly.shell_orders)
    @test result10.energy ≈ he_anchor_reference atol = 5.0e-6 rtol = 1.0e-6
    @test result15.energy ≈ he_anchor_reference atol = 5.0e-6 rtol = 1.0e-6
    @test abs(result15.energy - he_anchor_reference) ≤ abs(result10.energy - he_anchor_reference)
end

@testset "Angular Be legacy-trim one-body anchors" begin
    for order in (10, 15)
        benchmark = _paper_style_angular_one_body_benchmark_fixture(order)
        full_spectrum = sort(real(eigvals(Hermitian(benchmark.hamiltonian), Hermitian(benchmark.overlap))))
        exact_spectrum =
            sort(real(eigvals(Hermitian(benchmark.exact_hamiltonian), Hermitian(benchmark.exact_overlap))))
        low_count = 8
        closed_shell_noninteracting = 2.0 * (full_spectrum[1] + full_spectrum[2])
        one_s_like_count = count(<(-4.0), full_spectrum[1:10])

        @test benchmark.exact_common_lmax == 1
        @test benchmark.angular_assembly.shell_orders == fill(order, length(benchmark.angular_assembly.shell_orders))
        @test all(lcap ≥ benchmark.exact_common_lmax + 4 for lcap in benchmark.angular_assembly.shell_kinetic_lcap)
        @test one_s_like_count == 1
        @test full_spectrum[1] < -7.0
        @test full_spectrum[2] > -3.0
        @test closed_shell_noninteracting > -22.0
        @test closed_shell_noninteracting < -18.0
        @test full_spectrum[1:low_count] ≈ exact_spectrum[1:low_count] atol = 1.0e-8 rtol = 1.0e-8
    end
end

@testset "Angular Ne legacy-trim one-body discriminator" begin
    for order in (10, 15, 32)
        benchmark = _paper_style_angular_one_body_benchmark_fixture(order; Z = 10.0, lmax = 6)
        full_spectrum = sort(real(eigvals(Hermitian(benchmark.hamiltonian), Hermitian(benchmark.overlap))))
        exact_spectrum =
            sort(real(eigvals(Hermitian(benchmark.exact_hamiltonian), Hermitian(benchmark.exact_overlap))))
        low_count = 8
        nocc = 5
        closed_shell_noninteracting = 2.0 * sum(full_spectrum[1:nocc])
        exact_noninteracting = 2.0 * sum(exact_spectrum[1:nocc])
        one_s_like_count = count(<(-20.0), full_spectrum[1:nocc])

        @test benchmark.angular_assembly.shell_orders == fill(order, length(benchmark.angular_assembly.shell_orders))
        if order in (10, 15)
            @test benchmark.exact_common_lmax == 1
        else
            @test benchmark.exact_common_lmax ≥ 3
        end
        @test one_s_like_count == 1
        @test full_spectrum[1] < -40.0
        @test full_spectrum[2] > -20.0
        @test closed_shell_noninteracting > -201.0
        @test closed_shell_noninteracting < -199.0
        @test closed_shell_noninteracting ≈ exact_noninteracting atol = 1.0e-5 rtol = 1.0e-7
        @test full_spectrum[1:low_count] ≈ exact_spectrum[1:low_count] atol = 1.0e-6 rtol = 1.0e-8
    end
end

@testset "Angular Ne legacy-trim interaction moment-span discriminator" begin
    for order in (10, 15)
        benchmark = _paper_style_angular_one_body_benchmark_fixture(order; Z = 10.0, lmax = 6)
        assembly = benchmark.angular_assembly
        bare_moment_lmax =
            maximum(maximum(keys(blocks)) for blocks in assembly.shell_moment_blocks)
        interaction_moment_lmax =
            maximum(maximum(keys(blocks)) for blocks in assembly.shell_interaction_moment_blocks)
        required_product_lmax = 2 * benchmark.exact_common_lmax
        current_interaction =
            GaussletBases._assemble_atomic_injected_angular_interaction(
                benchmark.radial_operators,
                assembly,
            )
        interaction_pair_plan = 2 * maximum(assembly.shell_interaction_lexpand)

        @test benchmark.exact_common_lmax == 1
        @test bare_moment_lmax == benchmark.exact_common_lmax
        @test bare_moment_lmax < required_product_lmax
        @test interaction_moment_lmax > bare_moment_lmax
        @test interaction_moment_lmax ≥ required_product_lmax
        @test interaction_pair_plan ≥ required_product_lmax
        @test opnorm(current_interaction - transpose(current_interaction), Inf) ≤ 1.0e-10
        @test minimum(diag(current_interaction)) > 0.0
    end
end

@testset "Angular Ne legacy-trim HF branch repair" begin
    benchmark10 = _paper_style_angular_hf_style_benchmark_fixture(10)
    benchmark15 = _paper_style_angular_hf_style_benchmark_fixture(15)
    benchmark32 = _paper_style_angular_hf_style_benchmark_fixture(32)
    energy10 = benchmark10.scf_result.energy
    energy15 = benchmark15.scf_result.energy
    energy32 = benchmark32.scf_result.energy

    @test energy10 < -128.4
    @test energy15 < -128.5
    @test abs(energy10 - energy32) < 0.1
    @test abs(energy15 - energy32) < 0.05
end

@testset "Atomic injected angular small-ED benchmark" begin
    benchmark = _atomic_injected_angular_small_ed_benchmark_fixture()
    diagnostics = atomic_injected_angular_small_ed_diagnostics(benchmark)
    exact_problem = benchmark.exact_reference_problem

    @test benchmark isa AtomicInjectedAngularSmallEDBenchmark
    @test benchmark.hf_style isa AtomicInjectedAngularHFStyleBenchmark
    @test benchmark.orbital_count == size(benchmark.hf_style.one_body.overlap, 1)
    @test benchmark.state_count == benchmark.orbital_count^2
    @test size(benchmark.orbital_overlap) == (benchmark.orbital_count, benchmark.orbital_count)
    @test size(benchmark.orbital_one_body) == size(benchmark.orbital_overlap)
    @test size(benchmark.orbital_interaction) == size(benchmark.orbital_overlap)
    @test exact_problem isa AtomicIDATwoElectronProblem

    @test diagnostics.orbital_overlap_symmetry_error ≤ 1.0e-10
    @test diagnostics.orbital_one_body_symmetry_error ≤ 1.0e-8
    @test diagnostics.orbital_interaction_symmetry_error ≤ 1.0e-10
    @test diagnostics.orbital_overlap_identity_error ≤ 2.0e-6
    @test diagnostics.state_overlap_identity_error_estimate ≤ 3.0e-6
    @test diagnostics.min_orbital_overlap_eigenvalue > 1.0e-6
    @test diagnostics.min_state_overlap_eigenvalue_estimate > 1.0e-6
    @test diagnostics.state_interaction_diagonal_min > 0.0
    @test diagnostics.state_interaction_diagonal_max > diagnostics.state_interaction_diagonal_min
    @test diagnostics.full_converged
    @test diagnostics.full_residual ≤ 1.0e-7
    @test diagnostics.exact_reference_energy ≈ 13.020668426715936 atol = 1.0e-8 rtol = 1.0e-8
    @test diagnostics.full_energy ≈ 12.97749161121589 atol = 1.0e-5 rtol = 1.0e-8
    @test diagnostics.energy_difference_to_exact_reference < -1.0e-3
end

@testset "Atomic Ylm one-body layer" begin
    rb, grid, radial_ops, channels, atom = _quick_radial_atomic_fixture()[1:5]

    @test length(channels) == 9
    @test channels[1] == YlmChannel(0, 0)
    @test channels[2] == YlmChannel(1, -1)
    @test channels[4] == YlmChannel(1, 1)
    @test channels[end] == YlmChannel(2, 2)
    @test size(atom.overlap) == (9 * length(rb), 9 * length(rb))
    @test size(atom.hamiltonian) == size(atom.overlap)

    for (i, channel) in enumerate(channels)
        block = channel_range(atom, i)
        direct_block = radial_ops.kinetic + radial_ops.nuclear + centrifugal(radial_ops, channel.l)
        @test block == ((i - 1) * length(rb) + 1):(i * length(rb))
        @test channel_range(atom, channel) == block
        @test channel_overlap(atom, i) ≈ radial_ops.overlap atol = 1.0e-12 rtol = 1.0e-12
        @test channel_overlap(atom, channel) ≈ radial_ops.overlap atol = 1.0e-12 rtol = 1.0e-12
        @test channel_hamiltonian(atom, i) ≈ direct_block atol = 1.0e-12 rtol = 1.0e-12
        @test channel_hamiltonian(atom, channel) ≈ direct_block atol = 1.0e-12 rtol = 1.0e-12

        for j in 1:length(channels)
            i == j && continue
            other = channel_range(atom, j)
            @test atom.overlap[block, other] == zeros(Float64, length(rb), length(rb))
            @test atom.hamiltonian[block, other] == zeros(Float64, length(rb), length(rb))
        end
    end
end

@testset "Gaunt table backend" begin
    table = GaussletBases.build_gaunt_table(2; Lmax = 4, atol = 1.0e-14, basis = :complex)
    rb, grid, radial_ops, channels, atom, ida = _quick_radial_atomic_fixture()

    @test table isa GaussletBases.GauntTable{Float64}
    @test GaussletBases.gaunt_lmax(table) == 2
    @test GaussletBases.gaunt_Lmax(table) == 4
    @test GaussletBases.gaunt_nnz(table) > 0
    @test GaussletBases.gaunt_hasblock(table, 0, 0, 0)
    @test !GaussletBases.gaunt_hasblock(table, 1, 0, 0)
    @test GaussletBases.gaunt_value(table, 0, 0, 0, 0, 0, 0) ≈ inv(sqrt(4 * pi)) atol = 1.0e-12 rtol = 1.0e-12

    for L in 0:GaussletBases.gaunt_Lmax(table)
        counted = 0
        for (l1, l2, entries) in GaussletBases.gaunt_each_block(table, L)
            @test GaussletBases.gaunt_legal_triple(l1, l2, L)
            previous = nothing
            for entry in entries
                counted += 1
                @test GaussletBases.gaunt_legal_ms(l1, entry.m1, l2, entry.m2, L, entry.M)
                @test entry.M == entry.m1 - entry.m2
                current = (entry.M, entry.m1, entry.m2)
                previous === nothing || @test previous <= current
                previous = current
            end
        end
        @test counted == GaussletBases.gaunt_nnz(table, L)

        expected_tensor = [
            GaussletBases.gaunt_value(
                table,
                L,
                channels[alpha].l,
                channels[alpha].m,
                channels[alphap].l,
                channels[alphap].m,
                M,
            )
            for alpha in 1:length(channels), alphap in 1:length(channels), M in -L:L
        ]
        @test gaunt_tensor(ida, L) ≈ expected_tensor atol = 1.0e-12 rtol = 1.0e-12
    end
end

@testset "Angular kernel sectorization" begin
    _, _, _, channels, _, ida = _quick_radial_atomic_fixture()
    sectors = ida.angular_sectors
    nchannels = length(channels)

    @test sectors.nchannels == nchannels
    @test length(sectors.pair_to_sector) == nchannels^2
    @test length(sectors.pair_to_local) == nchannels^2
    @test sum(length(sector.pair_indices) for sector in sectors.sectors) == nchannels^2

    for (sector_index, sector) in enumerate(sectors.sectors)
        @test issorted(sector.pair_indices)
        matrix_by_L = sectors.sector_matrices

        for local_index in eachindex(sector.pair_indices)
            pair_index = sector.pair_indices[local_index]
            alpha = sector.left_channel_indices[local_index]
            beta = sector.right_channel_indices[local_index]
            @test channels[alpha].m + channels[beta].m == sector.msum
            @test sectors.pair_to_sector[pair_index] == sector_index
            @test sectors.pair_to_local[pair_index] == local_index
        end

        @test length(matrix_by_L) == GaussletBases.gaunt_Lmax(ida.gaunt_table) + 1
        for level in matrix_by_L
            @test size(level[sector_index]) == (length(sector.pair_indices), length(sector.pair_indices))
        end
    end

    for L in 0:GaussletBases.gaunt_Lmax(ida.gaunt_table)
        expected_kernel = _direct_dense_angular_kernel(ida.gaunt_table, channels, L)
        @test angular_kernel(ida, L) ≈ expected_kernel atol = 1.0e-12 rtol = 1.0e-12

        for alpha in 1:nchannels, beta in 1:nchannels, alphap in 1:nchannels, betap in 1:nchannels
            value = GaussletBases._angular_kernel_sector_value(sectors, L, alpha, beta, alphap, betap)
            @test value ≈ expected_kernel[alpha, alphap, beta, betap] atol = 1.0e-12 rtol = 1.0e-12

            if channels[alpha].m + channels[beta].m != channels[alphap].m + channels[betap].m
                @test value == 0.0
            end
        end
    end
end

@testset "Hydrogen Ylm spectrum" begin
    rb, grid, radial_ops, atom = _quick_hydrogen_ylm_fixture()

    @test norm(atom.overlap - I, Inf) ≤ 1.0e-5
    spectrum = sort(real(eigen(Hermitian(atom.hamiltonian)).values))
    E0 = spectrum[1]

    l1_channels = [YlmChannel(1, m) for m in -1:1]
    l1_energies = [
        begin
            @test norm(channel_overlap(atom, channel) - I, Inf) ≤ 1.0e-5
            minimum(real(eigen(Hermitian(channel_hamiltonian(atom, channel))).values))
        end
        for channel in l1_channels
    ]
    l2_channels = [YlmChannel(2, m) for m in -2:2]
    l2_energies = [
        begin
            @test norm(channel_overlap(atom, channel) - I, Inf) ≤ 1.0e-5
            minimum(real(eigen(Hermitian(channel_hamiltonian(atom, channel))).values))
        end
        for channel in l2_channels
    ]

    @test abs(E0 + 0.5) ≤ 5.0e-5
    @test maximum(abs.(l1_energies .- l1_energies[1])) ≤ 1.0e-10
    @test maximum(abs.(l2_energies .- l2_energies[1])) ≤ 1.0e-10
    @test abs(l1_energies[1] + 0.125) ≤ 5.0e-3
    @test abs(l2_energies[1] + 1.0 / 18.0) ≤ 5.0e-3
end

end

if _test_group_enabled(:ida)
@testset "Atomic IDA ingredients" begin
    rb, grid, radial_ops, channels, atom, ida = _quick_radial_atomic_fixture()
    nchannels = length(channels)
    radial_dim = length(rb)

    @test ida isa AtomicIDAOperators
    @test length(orbitals(ida)) == nchannels * radial_dim
    @test size(radial_multipole(ida, 1)) == (radial_dim, radial_dim)
    @test size(gaunt_tensor(ida, 1)) == (nchannels, nchannels, 3)
    @test size(angular_kernel(ida, 1)) == (nchannels, nchannels, nchannels, nchannels)
    @test channel_overlap(ida, YlmChannel(1, 0)) ≈ channel_overlap(atom, YlmChannel(1, 0)) atol = 1.0e-12 rtol = 1.0e-12
    @test channel_hamiltonian(ida, YlmChannel(2, 1)) ≈ channel_hamiltonian(atom, YlmChannel(2, 1)) atol = 1.0e-12 rtol = 1.0e-12
    @test channel_range(ida, YlmChannel(1, -1)) == channel_range(atom, YlmChannel(1, -1))

    q0 = angular_kernel(ida, 0)
    for alpha in 1:nchannels, alphap in 1:nchannels, beta in 1:nchannels, betap in 1:nchannels
        expected = (alpha == alphap && beta == betap) ? 1.0 : 0.0
        @test q0[alpha, alphap, beta, betap] ≈ expected atol = 1.0e-12 rtol = 1.0e-12
    end

    @test gaunt_coefficient(ida, 0, 0, YlmChannel(0, 0), YlmChannel(0, 0)) ≈ inv(sqrt(4 * pi)) atol = 1.0e-12 rtol = 1.0e-12
    @test gaunt_coefficient(ida, 1, 0, YlmChannel(0, 0), YlmChannel(0, 0)) == 0.0
    @test gaunt_coefficient(ida, 1, 1, YlmChannel(1, 0), YlmChannel(1, 0)) == 0.0
    @test abs(gaunt_coefficient(ida, 1, 0, YlmChannel(1, 0), YlmChannel(0, 0))) > 1.0e-12

    orbital_data = orbitals(ida)
    @test orbital_data[1].index == 1
    @test orbital_data[1].channel == YlmChannel(0, 0)
    @test orbital_data[1].radial_index == 1
    @test orbital_data[radial_dim + 1].channel == YlmChannel(1, -1)
    @test orbital_data[radial_dim + 1].radial_index == 1
end

@testset "Atomic full-IDA dense export" begin
    _, _, _, channels, _, ida = _quick_radial_atomic_fixture()
    nchannels = length(channels)
    radial_dim = size(ida.radial_operators.overlap, 1)
    norbitals = length(orbitals(ida))
    shell_centers_r = Float64[Float64(value) for value in ida.radial_operators.shell_centers_r]
    expected_channel_l = Int[channel.l for channel in channels]
    expected_channel_m = Int[channel.m for channel in channels]
    package_version = string(Base.pkgversion(GaussletBases))
    perm = GaussletBases._atomic_shell_major_permutation(ida)
    expected_h1 = Matrix{Float64}(ida.one_body.hamiltonian[perm, perm])
    expected_vee = GaussletBases._ida_density_interaction_matrix(ida, orbitals(ida)[perm])
    expected_dims = fill(nchannels, radial_dim)
    expected_orders = collect(1:radial_dim)
    payload_data = fullida_dense_payload(
        ida;
        nelec = 2,
        meta = (example = "test_atomic_fullida_dense_export",),
    )

    mktempdir() do dir
        path = joinpath(dir, "atomic_fullida_dense_test.jld2")
        @test write_fullida_dense_jld2(
            path,
            ida;
            nelec = 2,
            meta = (example = "test_atomic_fullida_dense_export",),
        ) == path

        jldopen(path, "r") do file
            top_keys = Set(
                key isa AbstractVector ? join(string.(key), "/") : string(key) for key in keys(file)
            )
            @test "H1" in top_keys
            @test "Vee" in top_keys
            @test "dims_per_shell" in top_keys
            @test "orders" in top_keys
            @test "basis_centers" in top_keys
            @test "bridge" in top_keys
            @test "meta" in top_keys

            @test size(file["H1"]) == (norbitals, norbitals)
            @test size(file["Vee"]) == (norbitals, norbitals)
            @test size(file["basis_centers"]) == (norbitals, 3)
            @test Int.(file["dims_per_shell"]) == expected_dims
            @test Int.(file["orders"]) == expected_orders
            @test Matrix{Float64}(file["H1"]) ≈ expected_h1 atol = 0.0 rtol = 0.0
            @test Matrix{Float64}(file["Vee"]) ≈ expected_vee atol = 0.0 rtol = 0.0
            @test Matrix{Float64}(file["basis_centers"]) == zeros(norbitals, 3)

            @test String(file["bridge/format"]) == "fullida_dense_v1"
            @test Int(file["bridge/version"]) == 1
            @test String(file["bridge/site_type"]) == "Electron"
            @test String(file["bridge/interaction_model"]) == "density_density"
            @test String(file["bridge/model_detail"]) == "two_index_ida"
            @test String(file["bridge/source_branch"]) == "atomic_ida"
            @test String(file["bridge/onebody_key"]) == "H1"
            @test String(file["bridge/interaction_key"]) == "Vee"
            @test String(file["bridge/optional_interaction_key"]) == "Vps"
            @test !Bool(file["bridge/has_optional_interaction"])
            @test Int(file["bridge/norb"]) == norbitals
            @test Int(file["bridge/nelec"]) == 2
            @test Bool(file["bridge/has_nelec"])
            @test String(file["bridge/order/convention"]) == "shell_major_by_radial_index"
            @test String(file["bridge/order/within_shell"]) == "ylm_channel_order"
            @test Int.(file["bridge/order/dims_per_shell"]) == expected_dims
            @test Int.(file["bridge/order/shell_offsets"]) == collect(1:nchannels:(norbitals + 1))
            @test Int.(file["bridge/order/shell_index"]) == vcat((fill(shell, nchannels) for shell in 1:radial_dim)...)
            @test Int.(file["bridge/order/orders_per_shell"]) == expected_orders
            @test String(file["bridge/order/basis_centers_kind"]) == "origin_only_atomic_orbitals"
            @test Float64.(file["bridge/order/shell_centers_r"]) == shell_centers_r
            @test !any(isnan, Float64.(file["bridge/order/shell_centers_r"]))
            @test Float64(file["bridge/order/basis_radius"]) == maximum(shell_centers_r)
            @test Int.(file["bridge/order/permutation_from_in_memory"]) == perm
            @test file["bridge/order/shell_centers_r"] == payload_data.bridge_meta["order/shell_centers_r"]
            @test Matrix{Float64}(file["H1"]) == payload_data.payload["H1"]
            @test Matrix{Float64}(file["Vee"]) == payload_data.payload["Vee"]

            @test String(file["meta/producer"]) == "GaussletBases.write_fullida_dense_jld2"
            @test String(file["meta/producer_type"]) == "AtomicIDAOperators"
            @test Float64(file["meta/Z"]) == 2.0
            @test String(file["meta/source_branch"]) == "atomic_ida"
            @test String(file["meta/source_model"]) == "radial_atomic_ida"
            @test String(file["meta/interaction_model"]) == "density_density_ida"
            @test String(file["meta/interaction_detail"]) == "two_index_ida"
            @test Int(file["meta/nchannels"]) == nchannels
            @test Int(file["meta/nradial"]) == radial_dim
            @test String(file["meta/public_extent_kind"]) == "count"
            @test isnan(Float64(file["meta/public_rmax"]))
            @test Int(file["meta/public_count"]) == 6
            @test !Bool(file["meta/has_public_rmax"])
            @test Bool(file["meta/has_public_count"])
            @test String(file["meta/manifest/producer/package"]) == "GaussletBases"
            @test String(file["meta/manifest/producer/version"]) == package_version
            @test String(file["meta/manifest/producer/entrypoint"]) == "GaussletBases.write_fullida_dense_jld2"
            @test String(file["meta/manifest/producer/object_type"]) == "AtomicIDAOperators"
            @test String(file["meta/manifest/interaction/model"]) == "density_density_ida"
            @test String(file["meta/manifest/interaction/detail"]) == "two_index_ida"
            @test String(file["meta/manifest/source/branch"]) == "atomic_ida"
            @test String(file["meta/manifest/source/model"]) == "radial_atomic_ida"
            @test Float64(file["meta/manifest/source/atomic_charge"]) == 2.0
            @test String(file["meta/manifest/source/basis_spec_type"]) == "RadialBasisSpec"
            @test String(file["meta/manifest/source/basis_family"]) == "G10"
            @test String(file["meta/manifest/source/public_extent_kind"]) == "count"
            @test isnan(Float64(file["meta/manifest/source/public_rmax"]))
            @test Int(file["meta/manifest/source/public_count"]) == 6
            @test !Bool(file["meta/manifest/source/has_public_rmax"])
            @test Bool(file["meta/manifest/source/has_public_count"])
            @test Float64(file["meta/manifest/source/reference_spacing"]) == 1.0
            @test Int(file["meta/manifest/source/tails"]) == 3
            @test Int(file["meta/manifest/source/odd_even_kmax"]) == 2
            @test String(file["meta/manifest/source/supplement_kind"]) == "xgaussian"
            @test Int(file["meta/manifest/source/supplement_count"]) == 1
            @test Float64.(file["meta/manifest/source/supplement/xgaussian_alphas"]) == [0.2]
            @test String(file["meta/manifest/source/mapping/type"]) == "AsinhMapping"
            @test !Bool(file["meta/manifest/source/mapping/is_identity"])
            @test Float64(file["meta/manifest/source/mapping/a"]) == 1.0
            @test Float64(file["meta/manifest/source/mapping/c"]) == 0.15
            @test Float64(file["meta/manifest/source/mapping/s"]) == 0.15
            @test Float64(file["meta/manifest/source/mapping/tail_spacing"]) == 10.0
            @test Int(file["meta/manifest/source/radial_dimension"]) == radial_dim
            @test Int(file["meta/manifest/source/channel_count"]) == nchannels
            @test Int(file["meta/manifest/source/channel_lmax"]) == channels.lmax
            @test Int.(file["meta/manifest/source/channel_l"]) == expected_channel_l
            @test Int.(file["meta/manifest/source/channel_m"]) == expected_channel_m
            @test String(file["meta/manifest/source/channel_convention"]) == "ylm_channels_increasing_l_then_m"
            @test String(file["meta/example"]) == "test_atomic_fullida_dense_export"
        end
    end
end

@testset "Atomic sliced Hamiltonian export" begin
    _, _, _, channels, _, ida = _quick_radial_atomic_fixture()
    nchannels = length(channels)
    radial_dim = size(ida.radial_operators.overlap, 1)
    norbitals = length(orbitals(ida))
    shell_centers_r = Float64[Float64(value) for value in ida.radial_operators.shell_centers_r]
    expected_channel_l = Int[channel.l for channel in channels]
    expected_channel_m = Int[channel.m for channel in channels]
    package_version = string(Base.pkgversion(GaussletBases))
    orbital_perm, channel_perm = GaussletBases._atomic_sliced_permutation(ida)
    ordered_channels = ida.one_body.channels.channel_data[channel_perm]
    expected_dims = fill(nchannels, radial_dim)
    expected_offs = collect(1:nchannels:(norbitals + 1))
    expected_m = Int[channel.m for channel in ordered_channels]
    expected_l = Int[channel.l for channel in ordered_channels]
    expected_labels = String["r=1,l=$(channel.l),m=$(channel.m)" for channel in ordered_channels]
    components = GaussletBases._atomic_onebody_component_matrices(ida)
    expected_h1 = Matrix{Float64}(components.H1[orbital_perm, orbital_perm])
    expected_t = Matrix{Float64}(components.T[orbital_perm, orbital_perm])
    expected_vnuc = Matrix{Float64}(components.Vnuc[orbital_perm, orbital_perm])
    expected_vee = GaussletBases._ida_density_interaction_matrix(ida, orbitals(ida)[orbital_perm])
    payload_data = sliced_ham_payload(
        ida;
        nelec = 2,
        meta = (example = "test_atomic_sliced_export",),
    )

    mktempdir() do dir
        path = joinpath(dir, "atomic_sliced_ham_test.jld2")
        @test write_sliced_ham_jld2(
            path,
            ida;
            nelec = 2,
            meta = (example = "test_atomic_sliced_export",),
        ) == path

        jldopen(path, "r") do file
            top_keys = Set(
                key isa AbstractVector ? join(string.(key), "/") : string(key) for key in keys(file)
            )
            @test "layout" in top_keys
            @test "basis" in top_keys
            @test "ordering" in top_keys
            @test "onebody" in top_keys
            @test "twobody" in top_keys
            @test "meta" in top_keys

            @test Int(file["layout/nslices"]) == radial_dim
            @test Int.(file["layout/dims"]) == expected_dims
            @test Int.(file["layout/offs"]) == expected_offs
            @test Float64.(file["layout/slice_coord"]) == shell_centers_r
            @test !any(isnan, Float64.(file["layout/slice_coord"]))
            @test Int.(file["layout/slice_index"]) == collect(1:radial_dim)
            @test file["layout/slice_coord"] == payload_data.layout_values["slice_coord"]

            m_by_slice = [Int.(collect(v)) for v in file["basis/m_by_slice"]]
            l_by_slice = [Int.(collect(v)) for v in file["basis/l_by_slice"]]
            labels_by_slice = [String.(collect(v)) for v in file["basis/labels_by_slice"]]
            @test length(m_by_slice) == radial_dim
            @test length(l_by_slice) == radial_dim
            @test length(labels_by_slice) == radial_dim
            @test m_by_slice[1] == expected_m
            @test l_by_slice[1] == expected_l
            @test Int.(file["basis/m_flat"]) == vcat(fill(expected_m, radial_dim)...)
            @test Int.(file["basis/l_flat"]) == vcat(fill(expected_l, radial_dim)...)
            @test labels_by_slice[1] == expected_labels

            @test String(file["ordering/major"]) == "slice_major"
            @test String(file["ordering/slice_meaning"]) == "radial_shell"
            @test String(file["ordering/within_slice"]) == "l0_desc_mzigzag"
            @test occursin("slice-major by radial index", String(file["ordering/description"]))
            @test Int.(file["ordering/permutation_from_in_memory"]) == orbital_perm

            @test String(file["onebody/stored"]) == "coo"
            @test Bool(file["onebody/is_hermitian"])
            @test occursin("centrifugal", String(file["onebody/decomposition"]))
            @test String(file["twobody/stored"]) == "coo_all"
            @test String(file["twobody/convention"]) == "density_density_pairdiag_v1"
            @test String(file["twobody/symmetry"]) == "pair_diagonal_density_density"
            @test occursin("density-density interaction model", String(file["twobody/description"]))

            H1blocks = file["onebody/H1blocks"]
            Tblocks = file["onebody/Tblocks"]
            Vnucblocks = file["onebody/Vnucblocks"]
            Vblocks = file["twobody/Vblocks"]
            @test length(H1blocks) == radial_dim
            @test length(H1blocks[1]) == radial_dim
            @test length(Vblocks) == radial_dim
            @test length(Vblocks[1]) == radial_dim

            @test GaussletBases._coo_blocks_to_dense(H1blocks, expected_dims) ≈ expected_h1 atol = 0.0 rtol = 0.0
            @test GaussletBases._coo_blocks_to_dense(Tblocks, expected_dims) ≈ expected_t atol = 0.0 rtol = 0.0
            @test GaussletBases._coo_blocks_to_dense(Vnucblocks, expected_dims) ≈ expected_vnuc atol = 0.0 rtol = 0.0
            @test GaussletBases._pairdiag_blocks_to_density_matrix(Vblocks, expected_dims) ≈ expected_vee atol = 0.0 rtol = 0.0

            @test String(file["meta/format"]) == "atomic_ida_sliced_v1"
            @test String(file["meta/consumer_shape"]) == "slicedmrgutils.HamIO"
            @test Float64(file["meta/Z"]) == 2.0
            @test String(file["meta/producer"]) == "GaussletBases.write_sliced_ham_jld2"
            @test String(file["meta/producer_type"]) == "AtomicIDAOperators"
            @test String(file["meta/source_branch"]) == "atomic_ida"
            @test String(file["meta/source_model"]) == "radial_atomic_ida"
            @test String(file["meta/interaction_model"]) == "density_density_ida"
            @test String(file["meta/interaction_detail"]) == "two_index_ida_pair_diagonal"
            @test String(file["meta/twobody_encoding"]) == "pair_diagonal_density_density"
            @test String(file["meta/slice_kind"]) == "radial_shell"
            @test String(file["meta/slice_coord_kind"]) == "physical_radial_center"
            @test String(file["meta/slice_index_kind"]) == "radial_index"
            @test String(file["meta/orbital_ordering"]) == "slice_major_by_radial_index_then_l0_desc_mzigzag"
            @test Int(file["meta/nchannels"]) == nchannels
            @test Int(file["meta/nradial"]) == radial_dim
            @test String(file["meta/public_extent_kind"]) == "count"
            @test isnan(Float64(file["meta/public_rmax"]))
            @test Int(file["meta/public_count"]) == 6
            @test !Bool(file["meta/has_public_rmax"])
            @test Bool(file["meta/has_public_count"])
            @test Int(file["meta/norb"]) == norbitals
            @test Int(file["meta/nelec"]) == 2
            @test Bool(file["meta/has_nelec"])
            @test Int.(file["meta/permutation_from_in_memory"]) == orbital_perm
            @test String(file["meta/manifest/producer/package"]) == "GaussletBases"
            @test String(file["meta/manifest/producer/version"]) == package_version
            @test String(file["meta/manifest/producer/entrypoint"]) == "GaussletBases.write_sliced_ham_jld2"
            @test String(file["meta/manifest/producer/object_type"]) == "AtomicIDAOperators"
            @test String(file["meta/manifest/interaction/model"]) == "density_density_ida"
            @test String(file["meta/manifest/interaction/detail"]) == "two_index_ida_pair_diagonal"
            @test String(file["meta/manifest/source/branch"]) == "atomic_ida"
            @test String(file["meta/manifest/source/model"]) == "radial_atomic_ida"
            @test Float64(file["meta/manifest/source/atomic_charge"]) == 2.0
            @test String(file["meta/manifest/source/basis_spec_type"]) == "RadialBasisSpec"
            @test String(file["meta/manifest/source/basis_family"]) == "G10"
            @test String(file["meta/manifest/source/public_extent_kind"]) == "count"
            @test isnan(Float64(file["meta/manifest/source/public_rmax"]))
            @test Int(file["meta/manifest/source/public_count"]) == 6
            @test !Bool(file["meta/manifest/source/has_public_rmax"])
            @test Bool(file["meta/manifest/source/has_public_count"])
            @test Float64(file["meta/manifest/source/reference_spacing"]) == 1.0
            @test Int(file["meta/manifest/source/tails"]) == 3
            @test Int(file["meta/manifest/source/odd_even_kmax"]) == 2
            @test String(file["meta/manifest/source/supplement_kind"]) == "xgaussian"
            @test Int(file["meta/manifest/source/supplement_count"]) == 1
            @test Float64.(file["meta/manifest/source/supplement/xgaussian_alphas"]) == [0.2]
            @test String(file["meta/manifest/source/mapping/type"]) == "AsinhMapping"
            @test !Bool(file["meta/manifest/source/mapping/is_identity"])
            @test Float64(file["meta/manifest/source/mapping/a"]) == 1.0
            @test Float64(file["meta/manifest/source/mapping/c"]) == 0.15
            @test Float64(file["meta/manifest/source/mapping/s"]) == 0.15
            @test Float64(file["meta/manifest/source/mapping/tail_spacing"]) == 10.0
            @test Int(file["meta/manifest/source/radial_dimension"]) == radial_dim
            @test Int(file["meta/manifest/source/channel_count"]) == nchannels
            @test Int(file["meta/manifest/source/channel_lmax"]) == channels.lmax
            @test Int.(file["meta/manifest/source/channel_l"]) == expected_channel_l
            @test Int.(file["meta/manifest/source/channel_m"]) == expected_channel_m
            @test String(file["meta/manifest/source/channel_convention"]) == "ylm_channels_increasing_l_then_m"
            @test String(file["meta/example"]) == "test_atomic_sliced_export"
        end
    end
end

@testset "Atomic HamV6 compatibility export and interaction accessor" begin
    _, _, _, channels, _, ida = _quick_radial_atomic_fixture()
    nchannels = length(channels)
    radial_dim = size(ida.radial_operators.overlap, 1)
    norbitals = length(orbitals(ida))
    shell_centers_r = Float64[Float64(value) for value in ida.radial_operators.shell_centers_r]
    expected_dims = fill(nchannels, radial_dim)
    expected_offs = collect(1:nchannels:(norbitals + 1))
    package_version = string(Base.pkgversion(GaussletBases))
    shell_perm = GaussletBases._atomic_shell_major_permutation(ida)
    native_perm, _ = GaussletBases._atomic_sliced_permutation(ida)
    orbital_perm, channel_perm = GaussletBases._atomic_hamv6_permutation(ida)
    ordered_channels = ida.one_body.channels.channel_data[channel_perm]
    expected_m = Int[channel.m for channel in ordered_channels]
    expected_l = Int[channel.l for channel in ordered_channels]
    expected_labels = String["r=1,l=$(channel.l),m=$(channel.m)" for channel in ordered_channels]
    components = GaussletBases._atomic_onebody_component_matrices(ida)
    expected_h1 = Matrix{Float64}(components.H1[orbital_perm, orbital_perm])
    expected_t = Matrix{Float64}(components.T[orbital_perm, orbital_perm])
    expected_vnuc = Matrix{Float64}(components.Vnuc[orbital_perm, orbital_perm])
    expected_vee = GaussletBases._ida_density_interaction_matrix(ida, orbitals(ida)[orbital_perm])

    @test atomic_ida_density_interaction_matrix(ida) ≈
        GaussletBases._ida_density_interaction_matrix(ida, orbitals(ida)) atol = 0.0 rtol = 0.0
    @test atomic_ida_density_interaction_matrix(ida; ordering = :shell_major) ≈
        GaussletBases._ida_density_interaction_matrix(ida, orbitals(ida)[shell_perm]) atol = 0.0 rtol = 0.0
    @test atomic_ida_density_interaction_matrix(ida; ordering = :sliced_native) ≈
        GaussletBases._ida_density_interaction_matrix(ida, orbitals(ida)[native_perm]) atol = 0.0 rtol = 0.0
    @test atomic_ida_density_interaction_matrix(ida; ordering = :hamv6) ≈ expected_vee atol = 0.0 rtol = 0.0
    @test atomic_ida_density_interaction_matrix(ida; ordering = orbital_perm) ≈ expected_vee atol = 0.0 rtol = 0.0

    payload_data = atomic_hamv6_payload(
        ida;
        nelec = 2,
        meta = (example = "test_atomic_hamv6_export",),
    )

    mktempdir() do dir
        path = joinpath(dir, "atomic_hamv6_test.jld2")
        @test write_atomic_hamv6_jld2(
            path,
            ida;
            nelec = 2,
            meta = (example = "test_atomic_hamv6_export",),
        ) == path

        jldopen(path, "r") do file
            top_keys = Set(
                key isa AbstractVector ? join(string.(key), "/") : string(key) for key in keys(file)
            )
            @test "layout" in top_keys
            @test "basis" in top_keys
            @test "ordering" in top_keys
            @test "onebody" in top_keys
            @test "twobody" in top_keys
            @test "meta" in top_keys

            @test Int(file["layout/nslices"]) == radial_dim
            @test Int.(file["layout/dims"]) == expected_dims
            @test Int.(file["layout/offs"]) == expected_offs
            @test Float64.(file["layout/slice_coord"]) == shell_centers_r
            @test Int.(file["layout/slice_index"]) == collect(1:radial_dim)

            m_by_slice = [Int.(collect(v)) for v in file["basis/m_by_slice"]]
            l_by_slice = [Int.(collect(v)) for v in file["basis/l_by_slice"]]
            labels_by_slice = [String.(collect(v)) for v in file["basis/labels_by_slice"]]
            @test m_by_slice[1] == expected_m
            @test l_by_slice[1] == expected_l
            @test Int.(file["basis/m_flat"]) == vcat(fill(expected_m, radial_dim)...)
            @test Int.(file["basis/l_flat"]) == vcat(fill(expected_l, radial_dim)...)
            @test labels_by_slice[1] == expected_labels

            @test String(file["ordering/major"]) == "slice_major"
            @test String(file["ordering/slice_meaning"]) == "radial_shell"
            @test String(file["ordering/within_slice"]) == "mzigzag_then_l"
            @test occursin("for each m, l increasing", String(file["ordering/description"]))
            @test Int.(file["ordering/permutation_from_in_memory"]) == orbital_perm

            H1blocks = file["onebody/H1blocks"]
            Tblocks = file["onebody/Tblocks"]
            Vnucblocks = file["onebody/Vnucblocks"]
            Vblocks = file["twobody/Vblocks"]
            @test GaussletBases._coo_blocks_to_dense(H1blocks, expected_dims) ≈ expected_h1 atol = 0.0 rtol = 0.0
            @test GaussletBases._coo_blocks_to_dense(Tblocks, expected_dims) ≈ expected_t atol = 0.0 rtol = 0.0
            @test GaussletBases._coo_blocks_to_dense(Vnucblocks, expected_dims) ≈ expected_vnuc atol = 0.0 rtol = 0.0
            @test GaussletBases._pairdiag_blocks_to_density_matrix(Vblocks, expected_dims) ≈ expected_vee atol = 0.0 rtol = 0.0

            @test String(file["meta/format"]) == "atomic_hamv6_v1"
            @test String(file["meta/consumer_shape"]) == "slicedmrgutils.HamIO/HamV6"
            @test Float64(file["meta/Z"]) == 2.0
            @test String(file["meta/source_model"]) == "radial_atomic_ida"
            @test String(file["meta/interaction_model"]) == "density_density_ida"
            @test String(file["meta/interaction_detail"]) == "two_index_ida_pair_diagonal"
            @test String(file["meta/orbital_ordering"]) == "slice_major_by_radial_index_then_mzigzag_then_l"
            @test Int(file["meta/norb"]) == norbitals
            @test Int(file["meta/nelec"]) == 2
            @test String(file["meta/public_extent_kind"]) == "count"
            @test Int.(file["meta/permutation_from_in_memory"]) == orbital_perm
            @test String(file["meta/manifest/producer/package"]) == "GaussletBases"
            @test String(file["meta/manifest/producer/version"]) == package_version
            @test String(file["meta/manifest/producer/entrypoint"]) == "GaussletBases.write_atomic_hamv6_jld2"
            @test String(file["meta/example"]) == "test_atomic_hamv6_export"

            @test file["layout/slice_coord"] == payload_data.layout_values["slice_coord"]
            @test file["basis/m_flat"] == payload_data.basis_values["m_flat"]
            @test file["ordering/permutation_from_in_memory"] == payload_data.ordering_values["permutation_from_in_memory"]
            @test file["meta/Z"] == payload_data.meta_values["Z"]
        end
    end
end

@testset "Angular benchmark exact HamV6 bridge export" begin
    benchmark = _atomic_injected_angular_small_ed_benchmark_fixture()
    hf_payload = angular_benchmark_exact_hamv6_payload(
        benchmark.hf_style;
        nelec = 2,
        meta = (example = "test_angular_bridge_export",),
    )
    payload = angular_benchmark_exact_hamv6_payload(
        benchmark;
        nelec = 2,
        meta = (example = "test_angular_bridge_export",),
    )
    reference_payload = atomic_hamv6_payload(
        benchmark.hf_style.exact_ida_reference;
        nelec = 2,
        meta = (example = "test_angular_bridge_export",),
    )

    @test payload.layout_values == reference_payload.layout_values
    @test payload.basis_values == reference_payload.basis_values
    @test payload.ordering_values == reference_payload.ordering_values
    @test payload.onebody_values == reference_payload.onebody_values
    @test payload.twobody_values == reference_payload.twobody_values
    @test payload.layout_values == hf_payload.layout_values
    @test payload.basis_values == hf_payload.basis_values
    @test payload.ordering_values == hf_payload.ordering_values
    @test payload.onebody_values == hf_payload.onebody_values
    @test payload.twobody_values == hf_payload.twobody_values

    @test payload.meta_values["Z"] == 2.0
    @test payload.meta_values["consumer_shape"] == "slicedmrgutils.HamIO/HamV6"
    @test payload.meta_values["angular_bridge_kind"] == "exact_common_low_l_reference"
    @test payload.meta_values["angular_bridge_scope"] == "exact_common_low_l_reference_only"
    @test payload.meta_values["angular_bridge_consumer_language"] == "slicedmrgutils.HamIO/HamV6"
    @test payload.meta_values["angular_exact_common_lmax"] == benchmark.hf_style.one_body.exact_common_lmax
    @test payload.meta_values["angular_shell_orders"] == benchmark.hf_style.one_body.angular_assembly.shell_orders
    @test payload.meta_values["producer"] == "GaussletBases.write_angular_benchmark_exact_hamv6_jld2"
    @test payload.meta_values["producer_type"] == "AtomicInjectedAngularSmallEDBenchmark"
    @test hf_payload.meta_values["producer_type"] == "AtomicInjectedAngularHFStyleBenchmark"

    mktempdir() do dir
        path = joinpath(dir, "angular_exact_hamv6_bridge_test.jld2")
        @test write_angular_benchmark_exact_hamv6_jld2(
            path,
            benchmark;
            nelec = 2,
            meta = (example = "test_angular_bridge_export",),
        ) == path

        jldopen(path, "r") do file
            @test String(file["ordering/within_slice"]) == "mzigzag_then_l"
            @test String(file["meta/producer"]) == "GaussletBases.write_angular_benchmark_exact_hamv6_jld2"
            @test String(file["meta/producer_type"]) == "AtomicInjectedAngularSmallEDBenchmark"
            @test String(file["meta/angular_bridge_kind"]) == "exact_common_low_l_reference"
            @test String(file["meta/angular_bridge_scope"]) == "exact_common_low_l_reference_only"
            @test Int(file["meta/angular_exact_common_lmax"]) == benchmark.hf_style.one_body.exact_common_lmax
            @test Int.(file["meta/angular_shell_orders"]) == benchmark.hf_style.one_body.angular_assembly.shell_orders
            @test Float64(file["meta/Z"]) == 2.0
            @test file["layout/dims"] == payload.layout_values["dims"]
            @test file["basis/l_flat"] == payload.basis_values["l_flat"]
        end
    end
end

@testset "Atomic export source metadata supports rmax-based recipes" begin
    _, _, radial_ops, atom = _quick_hydrogen_ylm_fixture()
    ida = atomic_ida_operators(radial_ops; lmax = atom.channels.lmax)
    dense_payload = fullida_dense_payload(ida)
    sliced_payload = sliced_ham_payload(ida)

    @test String(dense_payload.meta_values["manifest/source/public_extent_kind"]) == "rmax"
    @test dense_payload.meta_values["Z"] == 1.0
    @test String(dense_payload.meta_values["source_model"]) == "radial_atomic_ida"
    @test String(dense_payload.meta_values["public_extent_kind"]) == "rmax"
    @test dense_payload.meta_values["public_rmax"] == 30.0
    @test dense_payload.meta_values["manifest/source/public_rmax"] == 30.0
    @test !Bool(dense_payload.meta_values["manifest/source/has_public_count"])
    @test Bool(dense_payload.meta_values["manifest/source/has_public_rmax"])
    @test String(dense_payload.meta_values["manifest/source/mapping/type"]) == "AsinhMapping"
    @test sliced_payload.meta_values["manifest/source/public_rmax"] == dense_payload.meta_values["manifest/source/public_rmax"]
    @test sliced_payload.meta_values["manifest/source/atomic_charge"] == dense_payload.meta_values["manifest/source/atomic_charge"]
    @test sliced_payload.meta_values["Z"] == dense_payload.meta_values["Z"]
    @test sliced_payload.meta_values["manifest/source/channel_l"] == dense_payload.meta_values["manifest/source/channel_l"]
    @test sliced_payload.meta_values["manifest/source/channel_m"] == dense_payload.meta_values["manifest/source/channel_m"]
end

@testset "Atomic IDA direct matrix" begin
    rb, grid, radial_ops, channels, atom, ida = _quick_radial_atomic_fixture()
    radial_dim = length(rb)
    nchannels = length(channels)
    norbitals = length(orbitals(ida))
    orbital_index(channel_index, radial_index) = (channel_index - 1) * radial_dim + radial_index

    density = zeros(Float64, norbitals, norbitals)
    density[orbital_index(1, 1), orbital_index(1, 1)] = 1.2
    density[orbital_index(3, 1), orbital_index(3, 1)] = 0.4
    density[orbital_index(1, 1), orbital_index(3, 1)] = 0.3
    density[orbital_index(3, 1), orbital_index(1, 1)] = 0.3
    density[orbital_index(7, 2), orbital_index(7, 2)] = 0.5
    density[orbital_index(7, 2), orbital_index(3, 2)] = -0.15
    density[orbital_index(3, 2), orbital_index(7, 2)] = -0.15
    density[orbital_index(4, 1), orbital_index(8, 3)] = 0.9
    density[orbital_index(8, 3), orbital_index(4, 1)] = 0.9

    direct = direct_matrix(ida, density)
    reference = _dense_direct_reference(ida, density)

    @test size(direct) == (norbitals, norbitals)
    @test direct ≈ reference atol = 1.0e-12 rtol = 1.0e-12
    @test direct ≈ transpose(direct) atol = 1.0e-12 rtol = 1.0e-12

    radial_projected_density = copy(density)
    for left_channel in 1:nchannels, right_channel in 1:nchannels, left_radial in 1:radial_dim, right_radial in 1:radial_dim
        left_radial == right_radial && continue
        radial_projected_density[orbital_index(left_channel, left_radial), orbital_index(right_channel, right_radial)] = 0.0
    end
    @test direct ≈ direct_matrix(ida, radial_projected_density) atol = 1.0e-12 rtol = 1.0e-12

    for left_channel in 1:nchannels, right_channel in 1:nchannels, left_radial in 1:radial_dim, right_radial in 1:radial_dim
        left_radial == right_radial && continue
        @test direct[orbital_index(left_channel, left_radial), orbital_index(right_channel, right_radial)] == 0.0
    end

    selection_density = zeros(Float64, norbitals, norbitals)
    selection_density[orbital_index(3, 1), orbital_index(1, 1)] = 1.0
    selection_density[orbital_index(1, 1), orbital_index(3, 1)] = 1.0
    selection_direct = direct_matrix(ida, selection_density)
    target_delta_m = channels[1].m - channels[3].m

    for left_channel in 1:nchannels, right_channel in 1:nchannels, radial_index in 1:radial_dim
        delta_m = channels[left_channel].m - channels[right_channel].m
        value = selection_direct[orbital_index(left_channel, radial_index), orbital_index(right_channel, radial_index)]
        if delta_m != target_delta_m && delta_m != -target_delta_m
            @test abs(value) ≤ 1.0e-12
        end
    end
end

@testset "Atomic IDA exchange matrix" begin
    rb, grid, radial_ops, channels, atom, ida = _quick_radial_atomic_fixture()
    radial_dim = length(rb)
    nchannels = length(channels)
    norbitals = length(orbitals(ida))
    orbital_index(channel_index, radial_index) = (channel_index - 1) * radial_dim + radial_index

    density = zeros(Float64, norbitals, norbitals)
    density[orbital_index(1, 1), orbital_index(1, 1)] = 0.9
    density[orbital_index(3, 1), orbital_index(1, 1)] = 0.25
    density[orbital_index(1, 1), orbital_index(3, 1)] = 0.25
    density[orbital_index(1, 1), orbital_index(3, 2)] = -0.2
    density[orbital_index(3, 2), orbital_index(1, 1)] = -0.2
    density[orbital_index(7, 3), orbital_index(4, 3)] = 0.35
    density[orbital_index(4, 3), orbital_index(7, 3)] = 0.35

    exchange = exchange_matrix(ida, density)
    reference = _dense_exchange_reference(ida, density)

    @test size(exchange) == (norbitals, norbitals)
    @test exchange ≈ reference atol = 1.0e-12 rtol = 1.0e-12
    @test exchange ≈ transpose(exchange) atol = 1.0e-12 rtol = 1.0e-12

    zeroed_density = copy(density)
    zeroed_density[orbital_index(1, 1), orbital_index(3, 2)] = 0.0
    zeroed_density[orbital_index(3, 2), orbital_index(1, 1)] = 0.0
    zeroed_exchange = exchange_matrix(ida, zeroed_density)
    @test exchange != zeroed_exchange

    selection_density = zeros(Float64, norbitals, norbitals)
    selection_density[orbital_index(3, 1), orbital_index(1, 2)] = 1.0
    selection_density[orbital_index(1, 2), orbital_index(3, 1)] = 1.0
    selection_exchange = exchange_matrix(ida, selection_density)
    target_msum = channels[3].m + channels[1].m

    for left_channel in 1:nchannels, right_channel in 1:nchannels, left_radial in 1:radial_dim, right_radial in 1:radial_dim
        value = selection_exchange[orbital_index(left_channel, left_radial), orbital_index(right_channel, right_radial)]
        if channels[left_channel].m + channels[right_channel].m != target_msum
            @test abs(value) ≤ 1.0e-12
        end
    end
end

@testset "Atomic IDA Fock matrix" begin
    rb, grid, radial_ops, channels, atom, ida = _quick_radial_atomic_fixture()
    radial_dim = length(rb)
    norbitals = length(orbitals(ida))
    orbital_index(channel_index, radial_index) = (channel_index - 1) * radial_dim + radial_index

    density = zeros(Float64, norbitals, norbitals)
    density[orbital_index(1, 1), orbital_index(1, 1)] = 0.9
    density[orbital_index(3, 1), orbital_index(1, 1)] = 0.25
    density[orbital_index(1, 1), orbital_index(3, 1)] = 0.25
    density[orbital_index(1, 1), orbital_index(3, 2)] = -0.2
    density[orbital_index(3, 2), orbital_index(1, 1)] = -0.2
    density[orbital_index(4, 2), orbital_index(4, 2)] = 0.4

    direct = direct_matrix(ida, density)
    exchange = exchange_matrix(ida, density)
    fock = fock_matrix(ida, density)
    reference = ida.one_body.hamiltonian + _dense_direct_reference(ida, density) - _dense_exchange_reference(ida, density)
    reference = 0.5 .* (reference .+ transpose(reference))

    @test size(fock) == size(ida.one_body.hamiltonian)
    @test fock ≈ ida.one_body.hamiltonian + direct - exchange atol = 1.0e-12 rtol = 1.0e-12
    @test fock ≈ reference atol = 1.0e-12 rtol = 1.0e-12
    @test fock ≈ transpose(fock) atol = 1.0e-12 rtol = 1.0e-12

    block = channel_range(ida, YlmChannel(0, 0))
    @test size(fock[block, block]) == (radial_dim, radial_dim)
end

@testset "Atomic IDA spin-aware Fock matrices" begin
    rb, grid, radial_ops, channels, atom, ida = _quick_radial_atomic_fixture()
    radial_dim = length(rb)
    norbitals = length(orbitals(ida))
    orbital_index(channel_index, radial_index) = (channel_index - 1) * radial_dim + radial_index

    density_alpha = zeros(Float64, norbitals, norbitals)
    density_beta = zeros(Float64, norbitals, norbitals)

    density_alpha[orbital_index(1, 1), orbital_index(1, 1)] = 0.8
    density_alpha[orbital_index(3, 1), orbital_index(1, 1)] = 0.2
    density_alpha[orbital_index(1, 1), orbital_index(3, 1)] = 0.2
    density_alpha[orbital_index(1, 1), orbital_index(3, 2)] = -0.15
    density_alpha[orbital_index(3, 2), orbital_index(1, 1)] = -0.15

    density_beta[orbital_index(4, 2), orbital_index(4, 2)] = 0.5
    density_beta[orbital_index(2, 1), orbital_index(2, 1)] = 0.3
    density_beta[orbital_index(2, 1), orbital_index(4, 2)] = 0.1
    density_beta[orbital_index(4, 2), orbital_index(2, 1)] = 0.1

    total_density = density_alpha + density_beta
    direct_total = direct_matrix(ida, total_density)
    exchange_alpha = exchange_matrix(ida, density_alpha)
    exchange_beta = exchange_matrix(ida, density_beta)

    fock_alpha = fock_matrix_alpha(ida, density_alpha, density_beta)
    fock_beta = fock_matrix_beta(ida, density_alpha, density_beta)

    @test fock_alpha ≈ ida.one_body.hamiltonian + direct_total - exchange_alpha atol = 1.0e-12 rtol = 1.0e-12
    @test fock_beta ≈ ida.one_body.hamiltonian + direct_total - exchange_beta atol = 1.0e-12 rtol = 1.0e-12
    @test fock_alpha ≈ transpose(fock_alpha) atol = 1.0e-12 rtol = 1.0e-12
    @test fock_beta ≈ transpose(fock_beta) atol = 1.0e-12 rtol = 1.0e-12

    beta_shift = copy(density_beta)
    beta_shift[orbital_index(1, 2), orbital_index(1, 2)] += 0.4
    beta_shift[orbital_index(1, 2), orbital_index(3, 2)] -= 0.1
    beta_shift[orbital_index(3, 2), orbital_index(1, 2)] -= 0.1
    fock_alpha_shift = fock_matrix_alpha(ida, density_alpha, beta_shift)
    @test fock_alpha_shift - fock_alpha ≈ direct_matrix(ida, beta_shift - density_beta) atol = 1.0e-12 rtol = 1.0e-12

    alpha_shift = copy(density_alpha)
    alpha_shift[orbital_index(1, 3), orbital_index(1, 3)] += 0.25
    alpha_shift[orbital_index(1, 3), orbital_index(3, 3)] += 0.05
    alpha_shift[orbital_index(3, 3), orbital_index(1, 3)] += 0.05
    fock_beta_shift = fock_matrix_beta(ida, alpha_shift, density_beta)
    @test fock_beta_shift - fock_beta ≈ direct_matrix(ida, alpha_shift - density_alpha) atol = 1.0e-12 rtol = 1.0e-12

    closed_shell_density = density_alpha
    closed_alpha = fock_matrix_alpha(ida, closed_shell_density, closed_shell_density)
    closed_beta = fock_matrix_beta(ida, closed_shell_density, closed_shell_density)
    @test closed_alpha ≈ closed_beta atol = 1.0e-12 rtol = 1.0e-12
    @test closed_alpha ≈ ida.one_body.hamiltonian + direct_matrix(ida, 2.0 .* closed_shell_density) - exchange_matrix(ida, closed_shell_density) atol = 1.0e-12 rtol = 1.0e-12

    spinless_helper = fock_matrix(ida, density_alpha)
    @test spinless_helper ≈ ida.one_body.hamiltonian + direct_matrix(ida, density_alpha) - exchange_matrix(ida, density_alpha) atol = 1.0e-12 rtol = 1.0e-12
end

@testset "Atomic IDA UHF" begin
    _rb, _grid, _radial_ops, ida, exact_problem, scf = _tiny_atomic_ida_uhf_fixture()
    exact_energy = ground_state_energy(exact_problem)
    norbitals = length(orbitals(ida))

    coeffs = [1.0 0.0; 0.0 1.0; 0.0 0.0]
    @test density_matrix(coeffs) ≈ [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.0] atol = 1.0e-12 rtol = 1.0e-12
    @test density_matrix(view(coeffs, :, 1)) ≈ [1.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0] atol = 1.0e-12 rtol = 1.0e-12

    initial_alpha = zeros(Float64, norbitals, norbitals)
    initial_beta = zeros(Float64, norbitals, norbitals)
    initial_alpha[1, 1] = 1.0
    initial_beta[1, 1] = 1.0
    step = uhf_step(ida, initial_alpha, initial_beta; nalpha = 1, nbeta = 1)

    @test step.fock_alpha ≈ fock_matrix_alpha(ida, initial_alpha, initial_beta) atol = 1.0e-12 rtol = 1.0e-12
    @test step.fock_beta ≈ fock_matrix_beta(ida, initial_alpha, initial_beta) atol = 1.0e-12 rtol = 1.0e-12
    @test step.density_alpha ≈ density_matrix(step.occupied_coefficients_alpha) atol = 1.0e-12 rtol = 1.0e-12
    @test step.density_beta ≈ density_matrix(step.occupied_coefficients_beta) atol = 1.0e-12 rtol = 1.0e-12
    @test step.energy ≈ uhf_energy(ida, step.density_alpha, step.density_beta) atol = 1.0e-12 rtol = 1.0e-12

    @test scf.converged
    @test 1 <= scf.iterations <= 80
    @test !isempty(scf.energies)
    @test length(scf.energies) == scf.iterations
    @test length(scf.residuals) == scf.iterations
    @test scf.fock_alpha ≈ transpose(scf.fock_alpha) atol = 1.0e-12 rtol = 1.0e-12
    @test scf.fock_beta ≈ transpose(scf.fock_beta) atol = 1.0e-12 rtol = 1.0e-12
    @test scf.fock_alpha ≈ fock_matrix_alpha(ida, scf.density_alpha, scf.density_beta) atol = 1.0e-12 rtol = 1.0e-12
    @test scf.fock_beta ≈ fock_matrix_beta(ida, scf.density_alpha, scf.density_beta) atol = 1.0e-12 rtol = 1.0e-12
    @test scf.energy ≈ uhf_energy(ida, scf.density_alpha, scf.density_beta) atol = 1.0e-12 rtol = 1.0e-12
    @test norm(scf.fock_alpha - scf.fock_beta, Inf) ≤ 1.0e-8
    @test norm(scf.density_alpha - scf.density_beta, Inf) ≤ 1.0e-8
    @test scf.energy >= exact_energy - 1.0e-9
    @test scf.energy < -2.7
end

@testset "Atomic IDA two-electron problem" begin
    rb, grid, radial_ops, ida, problem = _tiny_atomic_ida_two_electron_fixture()
    norbitals = length(orbitals(problem))
    states = two_electron_states(problem)

    @test problem isa AtomicIDATwoElectronProblem
    @test size(problem.overlap) == (norbitals^2, norbitals^2)
    @test size(problem.one_body) == size(problem.overlap)
    @test size(problem.two_body) == size(problem.overlap)
    @test size(problem.hamiltonian) == size(problem.overlap)
    @test length(states) == norbitals^2
    @test states[1].index == 1
    @test states[1].up_orbital == orbitals(problem)[1]
    @test states[1].down_orbital == orbitals(problem)[1]
    @test states[2].up_orbital == orbitals(problem)[1]
    @test states[2].down_orbital == orbitals(problem)[2]
    @test states[norbitals + 1].up_orbital == orbitals(problem)[2]
    @test states[norbitals + 1].down_orbital == orbitals(problem)[1]
    @test norm(problem.orbital_overlap - I, Inf) ≤ 1.0e-5
    @test norm(problem.overlap - I, Inf) ≤ 2.0e-5
    @test problem.hamiltonian ≈ transpose(problem.hamiltonian) atol = 1.0e-12 rtol = 1.0e-12
    @test problem.overlap ≈ transpose(problem.overlap) atol = 1.0e-12 rtol = 1.0e-12
    @test problem.two_body ≈ Diagonal(diag(problem.two_body)) atol = 1.0e-12 rtol = 1.0e-12

    coefficients = collect(range(1.0, length = length(states)))
    @test apply_hamiltonian(problem, coefficients) ≈ problem.hamiltonian * coefficients atol = 1.0e-12 rtol = 1.0e-12
    @test apply_overlap(problem, coefficients) ≈ problem.overlap * coefficients atol = 1.0e-12 rtol = 1.0e-12

    E0 = ground_state_energy(problem)
    @test isfinite(E0)
    @test -3.2 < E0 < -2.5
end

@testset "Atomic IDA Lanczos" begin
    _rb, _grid, _radial_ops, _ida, problem = _tiny_atomic_ida_lanczos_fixture()
    dense = ground_state_energy(problem)
    lanczos = lanczos_ground_state(problem; krylovdim = 200, maxiter = 200, tol = 1.0e-7)

    @test lanczos.converged
    @test lanczos.iterations <= 200
    @test lanczos.residual ≤ 1.0e-7
    @test abs(lanczos.value - dense) ≤ 1.0e-10
    @test abs(norm(lanczos.vector) - 1.0) ≤ 1.0e-10
    @test norm(problem.orbital_overlap - I, Inf) ≤ 5.0e-6
    @test norm(problem.overlap - I, Inf) ≤ 1.0e-5
    @test -2.95 < dense < -2.85
end

end

if _test_group_enabled(:misc)
@testset "REPL displays" begin
    (
        family,
        map,
        ub_spec,
        rb_spec,
        ub,
        rb,
        grid,
        ops,
        rep,
        channels,
        atom,
        ida,
        _tiny_rb,
        _tiny_grid,
        _tiny_radial_ops,
        tiny_ida,
        tiny_problem,
    ) = _quick_display_fixture()

    @test sprint(show, family) == "GaussletFamily(:G10)"
    @test occursin("AsinhMapping(", sprint(show, map))
    @test occursin("UniformBasisSpec(", sprint(show, ub_spec))
    @test occursin("RadialBasisSpec(", sprint(show, rb_spec))
    @test occursin("UniformBasis(length=3", sprint(show, ub))
    @test occursin("RadialBasis(length=6", sprint(show, rb))
    @test occursin("RadialQuadratureGrid(length=", sprint(show, grid))
    @test occursin("RadialAtomicOperators(size=(6, 6)", sprint(show, ops))
    @test occursin("BasisRepresentation1D(kind=:uniform", sprint(show, rep))
    @test occursin("YlmChannel(l=1, m=0)", sprint(show, YlmChannel(1, 0)))
    @test occursin("YlmChannelSet(lmax=2, nchannels=9)", sprint(show, channels))
    @test occursin("AtomicOneBodyOperators(nchannels=9", sprint(show, atom))
    @test occursin("AtomicOrbital(index=1, channel=YlmChannel(l=0, m=0), radial_index=1)", sprint(show, orbitals(ida)[1]))
    @test occursin("AtomicIDAOperators(nchannels=9", sprint(show, ida))
    @test occursin("AtomicIDATwoElectronState(index=1", sprint(show, two_electron_states(tiny_problem)[1]))
    @test occursin("AtomicIDATwoElectronProblem(norbitals=", sprint(show, tiny_problem))
end

if _RUN_SLOW_TESTS
    @testset "REPL displays (slow advanced objects)" begin
        (
            _family,
            _map,
            _ub_spec,
            hb_spec,
            _rb_spec,
            _ub,
            hb,
            _rb,
            _grid,
            _ops,
            _rep,
            partition,
            hierarchy,
            pgdg,
            augmented_pgdg,
            spec,
            global_layer,
            contracted_layer,
            _channels,
            _atom,
            _ida,
        ) = _slow_display_fixture()

        @test occursin("HalfLineBasisSpec(", sprint(show, hb_spec))
        @test occursin("HalfLineBasis(length=", sprint(show, hb))
        @test occursin("BasisPartition1D(nbasis=3, nboxes=3)", sprint(show, partition))
        @test occursin("BasisBox1D(index=1", sprint(show, boxes(partition)[1]))
        @test occursin("HierarchicalBasisPartition1D(nbasis=3, nboxes=5, nleaves=4)", sprint(show, hierarchy))
        @test occursin("HierarchicalBasisBox1D(index=4", sprint(show, boxes(hierarchy)[4]))
        @test occursin("LeafLocalPGDG1D(nleaves=4, nbasis=8, primitives_per_leaf=2, naugmented=0)", sprint(show, pgdg))
        @test occursin("LeafLocalPGDG1D(nleaves=4, nbasis=9, primitives_per_leaf=2, naugmented=1)", sprint(show, augmented_pgdg))
        @test occursin("LeafGaussianSpec1D(relative_position=0.5, width_scale=0.2)", sprint(show, spec))
        @test occursin("GlobalMappedPrimitiveLayer1D(nbasis=", sprint(show, global_layer))
        @test occursin("LeafBoxContractionLayer1D(nleaves=4", sprint(show, contracted_layer))
    end
end

end

if _test_group_enabled(:docs)
@testset "Documentation consistency" begin
    design = read(joinpath(_PROJECT_ROOT, "DESIGN.md"), String)
    readme = read(joinpath(_PROJECT_ROOT, "README.md"), String)
    atomic_setup = read(joinpath(_PROJECT_ROOT, "docs", "recommended_atomic_setup.md"), String)
    radial_workflow = read(joinpath(_PROJECT_ROOT, "docs", "first_radial_workflow.md"), String)
    docs_index = read(joinpath(_PROJECT_ROOT, "docs", "index.md"), String)
    docs_structure_plan = read(joinpath(_PROJECT_ROOT, "docs", "documentation_structure_plan.md"), String)
    docs_consolidation_note = read(joinpath(_PROJECT_ROOT, "docs", "documentation_consolidation_note.md"), String)
    documenter_transition_plan = read(joinpath(_PROJECT_ROOT, "docs", "documenter_transition_plan.md"), String)
    api_reference_plan = read(joinpath(_PROJECT_ROOT, "docs", "api_reference_first_pass.md"), String)
    docs_hardening_plan = read(joinpath(_PROJECT_ROOT, "docs", "documentation_hardening_plan.md"), String)
    docs_presentation_plan = read(joinpath(_PROJECT_ROOT, "docs", "documentation_presentation_plan.md"), String)
    docs_navigation_note = read(joinpath(_PROJECT_ROOT, "docs", "documentation_navigation_feel_note.md"), String)
    docs_deployment_note = read(joinpath(_PROJECT_ROOT, "docs", "documentation_deployment_note.md"), String)
    docs_public_polish_note = read(joinpath(_PROJECT_ROOT, "docs", "public_polish_note.md"), String)
    docs_public_story_note = read(joinpath(_PROJECT_ROOT, "docs", "public_story_status_note.md"), String)
    onboarding_plan = read(joinpath(_PROJECT_ROOT, "docs", "onboarding_frontdoor_plan.md"), String)
    docs_project = read(joinpath(_PROJECT_ROOT, "docs", "Project.toml"), String)
    docs_make = read(joinpath(_PROJECT_ROOT, "docs", "make.jl"), String)
    docs_workflow = read(joinpath(_PROJECT_ROOT, ".github", "workflows", "docs.yml"), String)
    docs_site_index = read(joinpath(_PROJECT_ROOT, "docs", "src", "index.md"), String)
    docs_site_manual = read(joinpath(_PROJECT_ROOT, "docs", "src", "manual", "index.md"), String)
    docs_site_examples_index = read(joinpath(_PROJECT_ROOT, "docs", "src", "examples", "index.md"), String)
    docs_site_radial = read(joinpath(_PROJECT_ROOT, "docs", "src", "tutorials", "first_radial_workflow.md"), String)
    docs_site_atomic = read(joinpath(_PROJECT_ROOT, "docs", "src", "explanations", "current_atomic_branch.md"), String)
    docs_site_angular_track = read(joinpath(_PROJECT_ROOT, "docs", "src", "explanations", "angular_research_track.md"), String)
    docs_site_ordinary = read(joinpath(_PROJECT_ROOT, "docs", "src", "explanations", "current_ordinary_branch.md"), String)
    docs_site_examples = read(joinpath(_PROJECT_ROOT, "docs", "src", "howto", "example_guide.md"), String)
    docs_site_developer = read(joinpath(_PROJECT_ROOT, "docs", "src", "developer", "index.md"), String)
    docs_site_developer_notes = read(joinpath(_PROJECT_ROOT, "docs", "src", "developer", "supporting_notes.md"), String)
    docs_site_reference_index = read(joinpath(_PROJECT_ROOT, "docs", "src", "reference", "index.md"), String)
    docs_site_reference_bases = read(joinpath(_PROJECT_ROOT, "docs", "src", "reference", "bases_and_mappings.md"), String)
    docs_site_reference_ops = read(joinpath(_PROJECT_ROOT, "docs", "src", "reference", "operators_and_diagnostics.md"), String)
    docs_site_reference_atomic = read(joinpath(_PROJECT_ROOT, "docs", "src", "reference", "atomic_and_ordinary.md"), String)
    docs_site_reference_export = read(joinpath(_PROJECT_ROOT, "docs", "src", "reference", "export.md"), String)
    current_atomic_branch = read(joinpath(_PROJECT_ROOT, "docs", "current_atomic_branch.md"), String)
    current_ordinary_branch = read(joinpath(_PROJECT_ROOT, "docs", "current_ordinary_branch.md"), String)
    atomic_mean_field_supporting = read(joinpath(_PROJECT_ROOT, "docs", "atomic_mean_field_supporting_notes.md"), String)
    ordinary_pgdg_supporting = read(joinpath(_PROJECT_ROOT, "docs", "ordinary_pgdg_supporting_notes.md"), String)
    architecture = read(joinpath(_PROJECT_ROOT, "docs", "architecture.md"), String)
    primitive_layer_note = read(joinpath(_PROJECT_ROOT, "docs", "intermediate_primitive_layer.md"), String)
    radial_primitive_note = read(joinpath(_PROJECT_ROOT, "docs", "radial_primitive_operator_layer.md"), String)
    ordinary_coulomb_note = read(joinpath(_PROJECT_ROOT, "docs", "ordinary_coulomb_expansion_path.md"), String)
    mapped_ordinary_note = read(joinpath(_PROJECT_ROOT, "docs", "mapped_ordinary_basis.md"), String)
    short_gaucoulomb_note = read(joinpath(_PROJECT_ROOT, "docs", "short_gaucoulomb_backend.md"), String)
    ordinary_pgdg_note = read(joinpath(_PROJECT_ROOT, "docs", "ordinary_pgdg_decision.md"), String)
    ordinary_pgdg_comx_note = read(joinpath(_PROJECT_ROOT, "docs", "ordinary_pgdg_comx.md"), String)
    ordinary_pgdg_refinement_note = read(joinpath(_PROJECT_ROOT, "docs", "ordinary_pgdg_proxy_refinement.md"), String)
    ordinary_pgdg_backend_note = read(joinpath(_PROJECT_ROOT, "docs", "ordinary_pgdg_backend_pivot.md"), String)
    atomic_ylm_note = read(joinpath(_PROJECT_ROOT, "docs", "atomic_ylm_layer.md"), String)
    atomic_ida_note = read(joinpath(_PROJECT_ROOT, "docs", "atomic_ida_layer.md"), String)
    atomic_two_electron_note = read(joinpath(_PROJECT_ROOT, "docs", "atomic_ida_two_electron.md"), String)
    atomic_direct_note = read(joinpath(_PROJECT_ROOT, "docs", "atomic_ida_direct.md"), String)
    atomic_exchange_note = read(joinpath(_PROJECT_ROOT, "docs", "atomic_ida_exchange.md"), String)
    atomic_fock_note = read(joinpath(_PROJECT_ROOT, "docs", "atomic_ida_fock.md"), String)
    atomic_spin_fock_note = read(joinpath(_PROJECT_ROOT, "docs", "atomic_ida_spin_fock.md"), String)
    atomic_uhf_note = read(joinpath(_PROJECT_ROOT, "docs", "atomic_ida_uhf.md"), String)
    atomic_export_note = read(joinpath(_PROJECT_ROOT, "docs", "hamiltonian_export_fullida_dense.md"), String)
    atomic_sliced_export_note = read(joinpath(_PROJECT_ROOT, "docs", "hamiltonian_export_sliced_blocks.md"), String)
    atomic_export_smoke_note = read(joinpath(_PROJECT_ROOT, "docs", "atomic_export_consumer_smoke.md"), String)
    atomic_export_active_note = read(joinpath(_PROJECT_ROOT, "docs", "atomic_export_active_space_validation.md"), String)
    gaunt_backend_note = read(joinpath(_PROJECT_ROOT, "docs", "gaunt_backend_note.md"), String)
    example_guide = read(joinpath(_PROJECT_ROOT, "docs", "example_guide.md"), String)
    ordinary_one_body_note = read(joinpath(_PROJECT_ROOT, "docs", "ordinary_pgdg_one_body_fidelity.md"), String)
    ordinary_hybrid_note = read(joinpath(_PROJECT_ROOT, "docs", "ordinary_pgdg_hybrid_regime.md"), String)
    ordinary_sho_note = read(joinpath(_PROJECT_ROOT, "docs", "ordinary_sho_spectral_test.md"), String)
    ordinary_vee_note = read(joinpath(_PROJECT_ROOT, "docs", "ordinary_cartesian_vee_validation.md"), String)
    ordinary_hybrid_vee_note = read(joinpath(_PROJECT_ROOT, "docs", "ordinary_cartesian_hybrid_vee_validation.md"), String)
    ordinary_hybrid_consolidation_note = read(joinpath(_PROJECT_ROOT, "docs", "ordinary_pgdg_hybrid_consolidation.md"), String)
    global_map_note = read(joinpath(_PROJECT_ROOT, "docs", "global_map_local_contraction.md"), String)
    leaf_pgdg_note = read(joinpath(_PROJECT_ROOT, "docs", "leaf_pgdg_1d.md"), String)
    global_contraction_note = read(joinpath(_PROJECT_ROOT, "docs", "global_mapped_leaf_contraction_1d.md"), String)
    terminology = read(joinpath(_PROJECT_ROOT, "docs", "terminology.md"), String)
    roadmap = read(joinpath(_PROJECT_ROOT, "ROADMAP.md"), String)
    status = read(joinpath(_PROJECT_ROOT, "STATUS.md"), String)

    @test !occursin("primitive_kinetic_matrix", design)
    @test !occursin("CombinedMapping", design)
    @test !occursin("ScaledMapping", design)
    @test !occursin("NoDiagonalApproximation", design)
    @test occursin("Gausslets are localized, orthogonal basis functions constructed from short", readme)
    @test occursin("They were developed to combine several properties", readme)
    @test occursin("mature **radial / atomic workflow**", readme)
    @test occursin("A first useful calculation", readme)
    @test occursin("Best first path through the repository", readme)
    @test occursin("examples/04_hydrogen_ground_state.jl", readme)
    @test occursin("examples/15_atomic_hydrogen_ylm.jl", readme)
    @test occursin("https://srwhite59.github.io/GaussletBases.jl/dev/", readme)
    @test occursin("https://img.shields.io/badge/docs-dev-blue.svg", readme)
    @test occursin("https://srwhite59.github.io/GaussletBases.jl/dev/manual/", readme)
    @test occursin("https://srwhite59.github.io/GaussletBases.jl/dev/reference/", readme)
    @test occursin("https://srwhite59.github.io/GaussletBases.jl/dev/developer/", readme)
    @test occursin("## Acknowledgments", readme)
    @test occursin("OpenAI Codex-style interactive", readme)
    @test occursin("author-driven", readme)
    @test occursin("AsinhMapping(c = c, s = s)", readme)
    @test occursin("ordinary Julia keyword arguments", readme)
    @test occursin("Lowest hydrogen energy", readme)
    @test occursin("diag = basis_diagnostics(rb)", readme)
    @test occursin("ops = atomic_operators(rb, grid; Z = Z, lmax = 2)", readme)
    @test !occursin("reference_spacing = 1.0", readme)
    @test !occursin("tails = 6", readme)
    @test !occursin("odd_even_kmax = 6", readme)
    @test !occursin("xgaussians = XGaussian[]", readme)
    @test occursin("Pkg.add(url = \"https://github.com/srwhite59/GaussletBases.jl\")", readme)
    @test occursin("s = 0.2", atomic_setup)
    @test occursin("odd_even_kmax = 6", atomic_setup)
    @test occursin("accuracy = :veryhigh", atomic_setup)
    @test occursin("examples/04_hydrogen_ground_state.jl", radial_workflow)
    @test occursin("Quickstart", radial_workflow)
    @test occursin("quadrature_rmax = 30.0", radial_workflow)
    @test occursin("ordinary Julia keyword arguments", radial_workflow)
    @test occursin("recommended_atomic_setup.md", radial_workflow)
    @test occursin("xgaussians", radial_workflow)
    @test occursin("example_guide.md", radial_workflow)
    @test occursin("if you are new", lowercase(docs_index))
    @test occursin("if you want radial atoms", lowercase(docs_index))
    @test occursin("if you want the ordinary cartesian branch", lowercase(docs_index))
    @test occursin("if you want advanced primitive and hierarchy work", lowercase(docs_index))
    @test occursin("if you want research and development notes", lowercase(docs_index))
    @test occursin("https://srwhite59.github.io/GaussletBases.jl/dev/", docs_index)
    @test occursin("atomic_mean_field_supporting_notes.md", docs_index)
    @test occursin("ordinary_pgdg_supporting_notes.md", docs_index)
    @test occursin("onboarding docs", lowercase(docs_structure_plan))
    @test occursin("current workflow docs", lowercase(docs_structure_plan))
    @test occursin("supporting notes", lowercase(docs_structure_plan))
    @test occursin("what duplication still remains", lowercase(docs_consolidation_note))
    @test occursin("atomic mean-field chain", lowercase(docs_consolidation_note))
    @test occursin("ordinary pgdg development chain", lowercase(docs_consolidation_note))
    @test occursin("documenter-based structure", lowercase(documenter_transition_plan))
    @test occursin("tutorials / getting started", lowercase(documenter_transition_plan))
    @test occursin("api reference from docstrings", lowercase(documenter_transition_plan))
    @test occursin("incremental", lowercase(documenter_transition_plan))
    @test occursin("curated api reference", lowercase(api_reference_plan))
    @test occursin("@docs", api_reference_plan)
    @test occursin("@autodocs", api_reference_plan)
    @test occursin("docstrings should answer questions", lowercase(api_reference_plan))
    @test occursin("docs build run automatically in ci", lowercase(docs_hardening_plan))
    @test occursin("small doctest slice", lowercase(docs_hardening_plan))
    @test occursin("checkdocs", lowercase(docs_hardening_plan))
    @test occursin("curated `@docs` pages", lowercase(docs_hardening_plan))
    @test occursin("user-facing julia package docs site", lowercase(docs_presentation_plan))
    @test occursin("home", lowercase(docs_presentation_plan))
    @test occursin("manual", lowercase(docs_presentation_plan))
    @test occursin("reference", lowercase(docs_presentation_plan))
    @test occursin("developer notes", lowercase(docs_presentation_plan))
    @test occursin("navigation feel", lowercase(docs_navigation_note))
    @test occursin("primary user documents", lowercase(docs_navigation_note))
    @test occursin("krylovkit", lowercase(docs_navigation_note))
    @test occursin("delivery surface", lowercase(docs_deployment_note))
    @test occursin("deploydocs", docs_deployment_note)
    @test occursin("github pages", lowercase(docs_deployment_note))
    @test occursin("public-facing polish pass", lowercase(docs_public_polish_note))
    @test occursin("examples more visible", lowercase(docs_public_polish_note))
    @test occursin("codex-style interactive coding assistance", lowercase(docs_public_polish_note))
    @test occursin("public opening of the", lowercase(docs_public_story_note))
    @test occursin("real second workflow", lowercase(docs_public_story_note))
    @test occursin("generalized eigensolve", lowercase(docs_public_story_note))
    @test occursin("front door", lowercase(onboarding_plan))
    @test occursin("what belongs in the readme", lowercase(onboarding_plan))
    @test occursin("what belongs in the radial quickstart", lowercase(onboarding_plan))
    @test occursin("what belongs only in the deeper setup note", lowercase(onboarding_plan))
    @test isdir(joinpath(_PROJECT_ROOT, "docs", "src"))
    @test occursin("Documenter", docs_project)
    @test occursin("makedocs", docs_make)
    @test occursin("doctest = true", docs_make)
    @test occursin("checkdocs = :none", docs_make)
    @test occursin("deploydocs(", docs_make)
    @test occursin("prettyurls = DOCS_CI", docs_make)
    @test occursin("canonical = DOCS_CI ?", docs_make)
    @test occursin("devbranch = \"main\"", docs_make)
    @test occursin("\"Manual\"", docs_make)
    @test occursin("\"Examples\"", docs_make)
    @test occursin("\"Reference\"", docs_make)
    @test occursin("\"Developer Notes\"", docs_make)
    @test occursin("hide(\"First radial workflow\"", docs_make)
    @test occursin("hide(\"Angular research track\"", docs_make)
    @test occursin("hide(\"Bases and mappings\"", docs_make)
    @test occursin("name: Docs", docs_workflow)
    @test occursin("julia --project=docs docs/make.jl", docs_workflow)
    @test occursin("Pkg.develop(PackageSpec(path=pwd()))", docs_workflow)
    @test occursin("contents: write", docs_workflow)
    @test occursin("GITHUB_TOKEN", docs_workflow)
    @test occursin("Start here", docs_site_index)
    @test occursin("Build the docs locally", docs_site_index)
    @test occursin("Primary documents", docs_site_index)
    @test occursin("Manual", docs_site_index)
    @test occursin("Examples", docs_site_index)
    @test occursin("Reference", docs_site_index)
    @test occursin("Developer Notes", docs_site_index)
    @test occursin("Manual, Algorithms, Examples, and Reference", docs_site_index)
    @test occursin("user-facing manual", lowercase(docs_site_manual))
    @test occursin("Who this manual is for", docs_site_manual)
    @test occursin("Recommended reading order", docs_site_manual)
    @test occursin("If you want more depth later", docs_site_manual)
    @test occursin("Developer Notes section is where", docs_site_manual)
    @test occursin("This page is the visible examples entry point", docs_site_examples_index)
    @test occursin("Start with these", docs_site_examples_index)
    @test occursin("Then choose a branch", docs_site_examples_index)
    @test occursin("Full curated guide", docs_site_examples_index)
    @test occursin("Bases and mappings", docs_site_radial)
    @test occursin("Operators and diagnostics", docs_site_radial)
    @test occursin("AsinhMapping(c = c, s = s)", docs_site_radial)
    @test occursin("Lowest hydrogen energy", docs_site_radial)
    @test occursin("current user-facing status read", lowercase(docs_site_atomic))
    @test occursin("density-density / IDA", docs_site_atomic)
    @test occursin("solver-facing export", lowercase(docs_site_atomic))
    @test occursin("developer/supporting material", lowercase(docs_site_atomic))
    @test occursin("current user-facing status read", lowercase(docs_site_ordinary))
    @test occursin("pgdg-style analytic route", lowercase(docs_site_ordinary))
    @test occursin("stress tests", lowercase(docs_site_ordinary))
    @test occursin("Core starting sequence", docs_site_examples)
    @test occursin("Ordinary Cartesian sequence", docs_site_examples)
    @test occursin("Developer Notes", docs_site_examples)
    @test occursin("lower-priority development, architecture, and history", lowercase(docs_site_developer))
    @test occursin("supporting note map", lowercase(docs_site_developer))
    @test occursin("flat", lowercase(docs_site_developer_notes))
    @test occursin("`docs/` tree", docs_site_developer_notes)
    @test occursin("jldoctest", docs_site_reference_bases)
    @test occursin("jldoctest", docs_site_reference_ops)
    @test occursin("jldoctest", docs_site_reference_atomic)
    @test occursin("First radial workflow", docs_site_reference_bases)
    @test occursin("Current atomic branch", docs_site_reference_atomic)
    @test occursin("Angular research track", docs_site_manual)
    @test occursin("active research track", lowercase(docs_site_angular_track))
    @test occursin("hooke remains important", lowercase(docs_site_angular_track))
    @test occursin("hf", lowercase(docs_site_angular_track))
    @test occursin("small ed", lowercase(docs_site_angular_track))
    @test occursin("hamio / hfdmrg-facing hf bridge", lowercase(docs_site_angular_track))
    @test occursin("sphere_point_set_orders", docs_site_angular_track)
    @test occursin("sphere_point_set(order)", docs_site_angular_track)
    @test occursin("curated_sphere_point_set", docs_site_angular_track)
    @test occursin("fibonacci_sphere_point_set", docs_site_angular_track)
    @test occursin("optimize_sphere_point_set", docs_site_angular_track)
    @test occursin("full vendored optimized", lowercase(docs_site_angular_track))
    @test occursin("sphere-point collection", lowercase(docs_site_angular_track))
    @test occursin("curated fixture subset", lowercase(docs_site_angular_track))
    @test occursin("explicit optional paths", lowercase(docs_site_angular_track))
    @test occursin("build_shell_local_injected_angular_basis", docs_site_angular_track)
    @test occursin("build_atomic_shell_local_angular_assembly", docs_site_angular_track)
    @test occursin("shell-to-atom angular assembly", lowercase(docs_site_angular_track))
    @test occursin("build_atomic_injected_angular_one_body_benchmark", docs_site_angular_track)
    @test occursin("one-electron angular benchmark", lowercase(docs_site_angular_track))
    @test occursin("build_atomic_injected_angular_hf_style_benchmark", docs_site_angular_track)
    @test occursin("hf-style benchmark", lowercase(docs_site_angular_track))
    @test occursin("build_atomic_injected_angular_hfdmrg_payload", docs_site_angular_track)
    @test occursin("build_atomic_injected_angular_hfdmrg_hf_adapter", docs_site_angular_track)
    @test occursin("build_atomic_injected_angular_hfdmrg_hf_seeds", docs_site_angular_track)
    @test occursin("run_atomic_injected_angular_hfdmrg_hf", docs_site_angular_track)
    @test occursin("in-memory hfdmrg-facing hf adapter", lowercase(docs_site_angular_track))
    @test occursin("primary stable handshake", lowercase(docs_site_angular_track))
    @test occursin("open-shell", lowercase(docs_site_angular_track))
    @test occursin("partially", lowercase(docs_site_angular_track))
    @test occursin("external paper-side drivers", lowercase(docs_site_angular_track))
    @test !occursin("tmp/work/angular_hfdmrg_payload_direct_scan.jl", docs_site_angular_track)
    @test occursin("build_atomic_injected_angular_small_ed_benchmark", docs_site_angular_track)
    @test occursin("small-ed benchmark", lowercase(docs_site_angular_track))
    @test occursin("write_angular_benchmark_exact_hamv6_jld2", docs_site_angular_track)
    @test occursin("exact common low-`l` reference", docs_site_angular_track)
    @test occursin("Angular research track", docs_site_atomic)
    @test occursin("shell-local injected angular basis construction and shell-to-atom assembly", lowercase(docs_site_atomic))
    @test occursin("experimental", lowercase(docs_site_atomic))
    @test !occursin("generalized eigen", lowercase(readme))
    @test !occursin("generalized eigen", lowercase(docs_site_index))
    @test !occursin("generalized eigen", lowercase(docs_site_manual))
    @test !occursin("generalized eigen", lowercase(docs_site_examples_index))
    @test occursin("first curated api-reference slice", lowercase(docs_site_reference_index))
    @test occursin("UniformBasisSpec", docs_site_reference_bases)
    @test occursin("AsinhMapping", docs_site_reference_bases)
    @test occursin("basis_diagnostics", docs_site_reference_ops)
    @test occursin("atomic_ida_operators", docs_site_reference_atomic)
    @test occursin("mapped_cartesian_hydrogen_energy", docs_site_reference_atomic)
    @test occursin("atomic_ida_density_interaction_matrix", docs_site_reference_export)
    @test occursin("atomic_hamv6_payload", docs_site_reference_export)
    @test occursin("angular_benchmark_exact_hamv6_payload", docs_site_reference_export)
    @test occursin("write_fullida_dense_jld2", docs_site_reference_export)
    @test occursin("write_atomic_hamv6_jld2", docs_site_reference_export)
    @test occursin("write_angular_benchmark_exact_hamv6_jld2", docs_site_reference_export)
    @test occursin("write_sliced_ham_jld2", docs_site_reference_export)
    @test occursin("build_atomic_injected_angular_hfdmrg_payload", docs_site_reference_export)
    @test occursin("HFDMRG.solve_hfdmrg", docs_site_reference_export)
    @test occursin("integral-diagonal approximation (IDA)", readme)
    @test occursin("integral-diagonal approximation (IDA)", docs_site_index)
    @test occursin("ordinary gausslets are the broad foundation", lowercase(architecture))
    @test occursin("radial gausslets are the current mature public-facing workflow", lowercase(architecture))
    @test occursin("first minimal atomic mean-field line", lowercase(architecture))
    @test occursin("it provides a small uhf kernel for the current atomic ida approximation", lowercase(architecture))
    @test occursin("global mapped layer plus local contraction", lowercase(architecture))
    @test occursin("where to navigate next", lowercase(architecture))
    @test occursin("docs/index.md", architecture)
    @test occursin("docs/current_atomic_branch.md", architecture)
    @test occursin("docs/current_ordinary_branch.md", architecture)
    @test occursin("PrimitiveSet1D", primitive_layer_note)
    @test occursin("BasisRepresentation1D", primitive_layer_note)
    @test occursin("HierarchicalBasisPartition1D", primitive_layer_note)
    @test occursin("GlobalMappedPrimitiveLayer1D", primitive_layer_note)
    @test occursin("LeafBoxContractionLayer1D", primitive_layer_note)
    @test occursin("LeafLocalPGDG1D", primitive_layer_note)
    @test occursin("the basis is not a black box", primitive_layer_note)
    @test occursin("one global mapped primitive layer", primitive_layer_note)
    @test occursin("primitive_set(rb)", radial_primitive_note)
    @test occursin("not to derive analytic formulas", lowercase(radial_primitive_note))
    @test occursin("prepares the way for ylm", lowercase(radial_primitive_note))
    @test occursin("coulomb expansion first", lowercase(ordinary_coulomb_note))
    @test occursin("shortgaucoulomb.jl", lowercase(ordinary_coulomb_note))
    @test occursin("3d grid", lowercase(ordinary_coulomb_note))
    @test occursin("mapped-`u` grid", ordinary_coulomb_note)
    @test occursin("MappedUniformBasisSpec", mapped_ordinary_note)
    @test occursin("fit_asinh_mapping_for_extent", mapped_ordinary_note)
    @test occursin("a general 3d tensor-product basis framework", lowercase(mapped_ordinary_note))
    @test occursin("current working asinh hydrogen study", lowercase(mapped_ordinary_note))
    @test occursin("not as a claim that", lowercase(mapped_ordinary_note))
    @test occursin("`asinhmapping` or the present coupled `c,s` tuning", lowercase(mapped_ordinary_note))
    @test occursin("source-of-truth", lowercase(short_gaucoulomb_note))
    @test occursin("almost verbatim", lowercase(short_gaucoulomb_note))
    @test occursin("coulombgaussianexpansion", lowercase(short_gaucoulomb_note))
    @test occursin("gaussian_factor_matrix", short_gaucoulomb_note)
    @test occursin("115", ordinary_pgdg_note)
    @test occursin("135", ordinary_pgdg_note)
    @test occursin("takeshi/yanai", lowercase(ordinary_pgdg_note))
    @test occursin("primitive/contraction-level analytic prototype", lowercase(ordinary_pgdg_note))
    @test occursin("comx/localization", lowercase(ordinary_pgdg_note))
    @test occursin("experimental", lowercase(ordinary_pgdg_comx_note))
    @test occursin("position-localization", lowercase(ordinary_pgdg_comx_note))
    @test occursin("post-comx analytic path", lowercase(ordinary_pgdg_comx_note))
    @test occursin("white–lindsey", lowercase(ordinary_pgdg_refinement_note))
    @test occursin("nearly identical span", lowercase(ordinary_pgdg_refinement_note))
    @test occursin("weighted log-quadratic gaussian fit", lowercase(ordinary_pgdg_refinement_note))
    @test occursin("hydrogen energy only as an end-to-end check", lowercase(ordinary_pgdg_refinement_note))
    @test occursin("numerical_reference", ordinary_pgdg_backend_note)
    @test occursin("pgdg_experimental", ordinary_pgdg_backend_note)
    @test occursin("mild-to-moderate distortion regime", lowercase(ordinary_pgdg_backend_note))
    @test occursin("validation route", lowercase(ordinary_pgdg_backend_note))
    @test occursin("YlmChannel", atomic_ylm_note)
    @test occursin("AtomicOneBodyOperators", atomic_ylm_note)
    @test occursin("hydrogen", lowercase(atomic_ylm_note))
    @test occursin("He / IDA", atomic_ylm_note)
    @test occursin("AtomicIDAOperators", atomic_ida_note)
    @test occursin("radial multipole", lowercase(atomic_ida_note))
    @test occursin("not solve the many-electron problem", lowercase(atomic_ida_note))
    @test occursin("channel-major", lowercase(atomic_ida_note))
    @test occursin("1 up, 1 down", atomic_two_electron_note)
    @test occursin("AtomicIDAOperators", atomic_two_electron_note)
    @test occursin("orthonormal", lowercase(atomic_two_electron_note))
    @test occursin("h \\otimes I + I \\otimes h", atomic_two_electron_note)
    @test occursin("ordinary Hermitian eigenproblem", atomic_two_electron_note)
    @test occursin("density-density", lowercase(atomic_two_electron_note))
    @test occursin("first genuinely interacting calculation", lowercase(atomic_two_electron_note))
    @test occursin("spatial one-particle density matrix", lowercase(atomic_direct_note))
    @test occursin("direct/hartree", lowercase(atomic_direct_note))
    @test occursin("radial-diagonal", lowercase(atomic_direct_note))
    @test occursin("dense `angular_kernel`", atomic_direct_note)
    @test occursin("exchange", lowercase(atomic_direct_note))
    @test occursin("exchange means", lowercase(atomic_exchange_note))
    @test occursin("radial-pair", lowercase(atomic_exchange_note))
    @test occursin("four-index coulomb contraction", lowercase(atomic_exchange_note))
    @test occursin("f = h + j - k", lowercase(atomic_exchange_note))
    @test occursin("f = h + j - k", lowercase(atomic_fock_note))
    @test occursin("not a full scf framework", lowercase(atomic_fock_note))
    @test occursin("choose occupations", lowercase(atomic_fock_note))
    @test occursin("uhf-style", lowercase(atomic_spin_fock_note))
    @test occursin("density_alpha", lowercase(atomic_spin_fock_note))
    @test occursin("density_beta", lowercase(atomic_spin_fock_note))
    @test occursin("uses the total density", lowercase(atomic_spin_fock_note))
    @test occursin("exchange term uses only the same-spin density", lowercase(atomic_spin_fock_note))
    @test occursin("spinless-model helper", lowercase(atomic_spin_fock_note))
    @test occursin("minimal uhf layer", lowercase(atomic_uhf_note))
    @test occursin("density_alpha", lowercase(atomic_uhf_note))
    @test occursin("density_beta", lowercase(atomic_uhf_note))
    @test occursin("orbital occupations", lowercase(atomic_uhf_note))
    @test occursin("uhf total energy", lowercase(atomic_uhf_note))
    @test occursin("fixed-point iteration", lowercase(atomic_uhf_note))
    @test occursin("fullida_dense_v1", atomic_export_note)
    @test occursin("dense bridge producer", lowercase(atomic_export_note))
    @test occursin("density-density", lowercase(atomic_export_note))
    @test occursin("two-index ida", lowercase(atomic_export_note))
    @test occursin("atomicidaoperators", lowercase(atomic_export_note))
    @test occursin("dense export came first", lowercase(atomic_sliced_export_note))
    @test occursin("layout/*", atomic_sliced_export_note)
    @test occursin("basis/*", atomic_sliced_export_note)
    @test occursin("ordering/*", atomic_sliced_export_note)
    @test occursin("onebody/*", atomic_sliced_export_note)
    @test occursin("twobody/*", atomic_sliced_export_note)
    @test occursin("pair-diagonal", lowercase(atomic_sliced_export_note))
    @test occursin("not a full four-index coulomb tensor", lowercase(atomic_sliced_export_note))
    @test occursin("two atomic export layers", lowercase(atomic_export_smoke_note))
    @test occursin("read_ham(...; validate=true)", atomic_export_smoke_note)
    @test occursin("compatible as written", lowercase(atomic_export_smoke_note))
    @test occursin("ordinary export", lowercase(atomic_export_smoke_note))
    @test occursin("schema compatibility is now established", lowercase(atomic_export_active_note))
    @test occursin("activespaceops.dense_h1", lowercase(atomic_export_active_note))
    @test occursin("build_dense_vblocks", lowercase(atomic_export_active_note))
    @test occursin("dense_h1", atomic_export_active_note)
    @test occursin("interaction_model = \"density_density_ida\"", atomic_export_active_note)
    @test occursin("defer ordinary export", lowercase(atomic_export_active_note))
    @test occursin("GauntTables", gaunt_backend_note)
    @test occursin("public atomic story should remain the same", lowercase(gaunt_backend_note))
    @test occursin("src/atomic_ida.jl", gaunt_backend_note)
    @test occursin("13_global_leaf_contraction.jl", example_guide)
    @test occursin("Core starting sequence", example_guide)
    @test occursin("Radial and atomic sequence", example_guide)
    @test occursin("Ordinary Cartesian sequence", example_guide)
    @test occursin("Primitive-layer and contraction sequence", example_guide)
    @test occursin("Hierarchy and current corrected direction", example_guide)
    @test occursin("Prototype side branch", example_guide)
    @test occursin("15_atomic_hydrogen_ylm.jl", example_guide)
    @test occursin("23_cartesian_hydrogen_coulomb_expansion.jl", example_guide)
    @test occursin("24_mapped_cartesian_hydrogen.jl", example_guide)
    @test occursin("25_mapped_cartesian_hydrogen_backends.jl", example_guide)
    @test occursin("33_ordinary_cartesian_1s2_vee.jl", example_guide)
    @test occursin("38_qiu_white_reference_vee.jl", example_guide)
    @test occursin("26_ordinary_cartesian_ida.jl", example_guide)
    @test occursin("27_ordinary_cartesian_ida_localized_backends.jl", example_guide)
    @test occursin("28_ordinary_one_body_fidelity.jl", example_guide)
    @test occursin("29_hybrid_mapped_cartesian_hydrogen.jl", example_guide)
    @test occursin("30_ordinary_sho_spectra.jl", example_guide)
    @test occursin("legacy/internal experimental regressions", lowercase(example_guide))
    @test occursin("16_atomic_ida_ingredients.jl", example_guide)
    @test occursin("19_atomic_ida_direct.jl", example_guide)
    @test occursin("20_atomic_ida_exchange.jl", example_guide)
    @test occursin("21_atomic_ida_fock.jl", example_guide)
    @test occursin("22_atomic_ida_uhf.jl", example_guide)
    @test occursin("31_atomic_fullida_dense_export.jl", example_guide)
    @test occursin("32_atomic_sliced_export.jl", example_guide)
    @test occursin("17_atomic_ida_two_electron.jl", example_guide)
    @test occursin("18_atomic_ida_two_electron_lanczos.jl", example_guide)
    @test occursin("docs/index.md", example_guide)
    @test occursin("docs/current_atomic_branch.md", example_guide)
    @test occursin("docs/current_ordinary_branch.md", example_guide)
    @test occursin("localized pgdg route", lowercase(ordinary_one_body_note))
    @test occursin("remaining issue is `H1`", ordinary_one_body_note)
    @test occursin("kinetic", lowercase(ordinary_one_body_note))
    @test occursin("aligned-kinetic route", ordinary_one_body_note)
    @test occursin("make the ordinary cartesian `vee` layer explicit", lowercase(ordinary_vee_note))
    @test occursin("doubly occupied noninteracting `1s` state", ordinary_vee_note)
    @test occursin("(5 / 8) z", lowercase(ordinary_vee_note))
    @test occursin("1.25 eh", lowercase(ordinary_vee_note))
    @test occursin("friendlier hybrid/core-supported regime", lowercase(ordinary_hybrid_vee_note))
    @test occursin("backbone-core gaussian contributions", lowercase(ordinary_hybrid_vee_note))
    @test occursin("there is **not yet** a separate residual-gaussian transfer", lowercase(replace(ordinary_hybrid_vee_note, "‑" => "-")))
    @test occursin("1.25 eh", lowercase(ordinary_hybrid_vee_note))
    @test occursin("radial branch should stay numerical", lowercase(ordinary_hybrid_note))
    @test occursin("white-lindsey", lowercase(replace(ordinary_hybrid_note, "–" => "-")))
    @test occursin("core gaussian", lowercase(ordinary_hybrid_note))
    @test occursin("stress-test regime", lowercase(ordinary_hybrid_note))
    @test occursin("low-energy sho", lowercase(ordinary_sho_note))
    @test occursin("low-momentum / low-order smooth completeness", lowercase(ordinary_sho_note))
    @test occursin("nearly identical span/subspace", lowercase(ordinary_sho_note))
    @test occursin("friendlier hybrid/core-supported regime", lowercase(ordinary_sho_note))
    @test occursin("ordinary pgdg backend is good enough", lowercase(ordinary_hybrid_consolidation_note))
    @test occursin("legacy `basissets` file", lowercase(ordinary_hybrid_consolidation_note))
    @test occursin("readbasis.jl", lowercase(ordinary_hybrid_consolidation_note))
    @test occursin("puregaussiangausslet.jl", lowercase(ordinary_hybrid_consolidation_note))
    @test occursin("old/basisgen.jl", lowercase(ordinary_hybrid_consolidation_note))
    @test occursin("mapping heuristics should be read as provisional", lowercase(ordinary_hybrid_consolidation_note))
    @test occursin("radial branch should stay numerical", lowercase(ordinary_hybrid_consolidation_note))
    @test occursin("cc-pvdz", lowercase(ordinary_hybrid_consolidation_note))
    @test occursin("cc-pvtznc", lowercase(ordinary_hybrid_consolidation_note))
    @test occursin("prototype", example_guide)
    @test occursin("supporting notes", lowercase(example_guide))
    @test occursin("what the atomic branch is today", lowercase(current_atomic_branch))
    @test occursin("atomic_ylm_layer.md", current_atomic_branch)
    @test occursin("atomic_ida_layer.md", current_atomic_branch)
    @test occursin("atomic_ida_uhf.md", current_atomic_branch)
    @test occursin("hamiltonian_export_fullida_dense.md", current_atomic_branch)
    @test occursin("hamiltonian_export_sliced_blocks.md", current_atomic_branch)
    @test occursin("atomic_export_consumer_smoke.md", current_atomic_branch)
    @test occursin("atomic_export_active_space_validation.md", current_atomic_branch)
    @test occursin("31_atomic_fullida_dense_export.jl", current_atomic_branch)
    @test occursin("32_atomic_sliced_export.jl", current_atomic_branch)
    @test occursin("supporting notes for the atomic line", lowercase(current_atomic_branch))
    @test occursin("atomic_mean_field_supporting_notes.md", current_atomic_branch)
    @test occursin("broad general atomic hf workflow", lowercase(current_atomic_branch))
    @test occursin("what the ordinary branch is today", lowercase(current_ordinary_branch))
    @test occursin("ordinary_cartesian_vee_validation.md", current_ordinary_branch)
    @test occursin("ordinary_cartesian_hybrid_vee_validation.md", current_ordinary_branch)
    @test occursin("ordinary_pgdg_hybrid_regime.md", current_ordinary_branch)
    @test occursin("ordinary_sho_spectral_test.md", current_ordinary_branch)
    @test occursin("ordinary_pgdg_hybrid_consolidation.md", current_ordinary_branch)
    @test occursin("supporting notes for the ordinary line", lowercase(current_ordinary_branch))
    @test occursin("ordinary_pgdg_supporting_notes.md", current_ordinary_branch)
    @test occursin("legacy/internal supporting notes", lowercase(current_ordinary_branch))
    @test occursin("asinhmapping is the current working map", lowercase(replace(current_ordinary_branch, "`" => "")))
    @test occursin("fixed `a = 1/(2Z)` with `s` solved from `count` and `xmax`", current_ordinary_branch)
    @test occursin("real second workflow", lowercase(current_ordinary_branch))
    @test occursin("recommended supporting-note order", lowercase(atomic_mean_field_supporting))
    @test occursin("atomic_ida_direct.md", atomic_mean_field_supporting)
    @test occursin("atomic_ida_spin_fock.md", atomic_mean_field_supporting)
    @test occursin("mergeable later", lowercase(atomic_mean_field_supporting))
    @test occursin("recommended supporting-note order", lowercase(ordinary_pgdg_supporting))
    @test occursin("ordinary_pgdg_decision.md", ordinary_pgdg_supporting)
    @test occursin("ordinary_pgdg_backend_pivot.md", ordinary_pgdg_supporting)
    @test occursin("asinhmapping", lowercase(ordinary_pgdg_supporting))
    @test startswith(atomic_direct_note, "> **Status:** supporting note.")
    @test startswith(atomic_exchange_note, "> **Status:** supporting note.")
    @test startswith(atomic_fock_note, "> **Status:** supporting note.")
    @test startswith(atomic_spin_fock_note, "> **Status:** supporting note.")
    @test startswith(ordinary_pgdg_note, "> **Status:** supporting development note.")
    @test startswith(ordinary_pgdg_comx_note, "> **Status:** supporting development note.")
    @test startswith(ordinary_pgdg_refinement_note, "> **Status:** supporting development note.")
    @test startswith(ordinary_pgdg_backend_note, "> **Status:** supporting development note.")
    @test startswith(global_map_note, "> **Note for new users:**")
    @test startswith(global_contraction_note, "> **Note for new users:**")
    @test startswith(leaf_pgdg_note, "> **Note for new users:**")
    @test occursin("LeafLocalPGDG1D", leaf_pgdg_note)
    @test occursin("1D only", leaf_pgdg_note)
    @test occursin("hierarchy-driven local basis generation", leaf_pgdg_note)
    @test occursin("no historical nested-driver port", leaf_pgdg_note)
    @test occursin("LeafGaussianSpec1D", leaf_pgdg_note)
    @test occursin("augment_leaf_pgdg", leaf_pgdg_note)
    @test occursin("GlobalMappedPrimitiveLayer1D", global_contraction_note)
    @test occursin("LeafBoxContractionLayer1D", global_contraction_note)
    @test occursin("one common basis over the whole region", global_contraction_note)
    @test occursin("the uncontracted global mapped basis is usable directly", lowercase(global_contraction_note))
    @test occursin("the basis is not the quadrature grid", terminology)
    @test occursin("Primitive set", terminology)
    @test occursin("contract_primitive_matrix", terminology)
    @test occursin("An exact non-diagonal radial electron-electron layer", roadmap)
    @test occursin("The first actual He / IDA-style solve", roadmap)
    @test occursin("geometry-aware grouping", roadmap)
    @test occursin("three directions are now real in the code", lowercase(roadmap))
    @test occursin("real but newer ordinary / cartesian workflow", lowercase(roadmap))
    @test !occursin("Gausslets.jl", readme)
    @test occursin("three clear layers in the public story", lowercase(status))
    @test occursin("mature extension of that workflow: the atomic angular line", lowercase(status))
    @test occursin("mature extension of that workflow: minimal atomic mean-field", lowercase(status))
    @test occursin("newer public-facing workflow: ordinary cartesian", lowercase(status))
    @test occursin("a broad general hf/rhf/uhf workflow", lowercase(status))
    @test occursin("prototype line inside the advanced research layer", lowercase(status))
    @test occursin("one global mapped primitive layer", status)
end

end

if _test_group_enabled(:examples)
@testset "Example scripts (quick smoke subset)" begin
    @test _run_example_script("01_first_gausslet.jl")
end

if _RUN_SLOW_TESTS
    @testset "Example scripts" begin
        @test _run_example_script("01_first_gausslet.jl")
        @test _run_example_script("02_radial_basis.jl")
        @test _run_example_script("03_radial_operators.jl")
        @test _run_example_script("04_hydrogen_ground_state.jl")
        @test _run_example_script("23_cartesian_hydrogen_coulomb_expansion.jl")
        @test _run_example_script("24_mapped_cartesian_hydrogen.jl")
        @test _run_example_script("25_mapped_cartesian_hydrogen_backends.jl")
        @test _run_example_script("05_primitive_sets.jl")
        @test _run_example_script("06_basis_contraction.jl")
        @test _run_example_script("07_position_contraction.jl")
        @test _run_example_script("08_basis_representation.jl")
        @test _run_example_script("09_basis_partition.jl")
        @test _run_example_script("10_hierarchical_partition.jl")
        @test _run_example_script("11_leaf_pgdg.jl")
        @test _run_example_script("12_leaf_pgdg_augmentation.jl")
        @test _run_example_script("13_global_leaf_contraction.jl")
        @test _run_example_script("14_radial_primitive_operators.jl")
        @test _run_example_script("15_atomic_hydrogen_ylm.jl")
        @test _run_example_script("16_atomic_ida_ingredients.jl")
        @test _run_example_script("19_atomic_ida_direct.jl")
        @test _run_example_script("20_atomic_ida_exchange.jl")
        @test _run_example_script("21_atomic_ida_fock.jl")
        @test _run_example_script("22_atomic_ida_uhf.jl")
        @test _run_example_script("17_atomic_ida_two_electron.jl")
        @test _run_example_script("18_atomic_ida_two_electron_lanczos.jl")
    end

    @testset "README example slice" begin
        Z = 2.0
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
        diag = basis_diagnostics(rb)
        grid = radial_quadrature(rb)
        ops = atomic_operators(rb, grid; Z = Z, lmax = 2)

        @test diag.overlap_error < 1.0e-5
        @test diag.D < 1.0e-3
        @test length(quadrature_points(grid)) > 0
        @test quadrature_weights(grid)[1] > 0.0
        @test size(ops.overlap) == (length(rb), length(rb))
        @test size(centrifugal(ops, 2)) == (length(rb), length(rb))
        @test size(multipole(ops, 1)) == (length(rb), length(rb))
    end
end

end
