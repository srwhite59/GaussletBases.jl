using GaussletBases
using LinearAlgebra
using Test

const CRD = GaussletBases.CartesianReferenceDensity
const CRG = GaussletBases.CartesianResidualGaussians
const CFR = GaussletBases.CartesianFinalBasisRealization
const CGR = GaussletBases.CartesianGaussianRawBlocks

struct RepresentedHartreeProxyLayer
    primitives::PrimitiveSet1D
    coefficients::Matrix{Float64}
end
struct RepresentedHartreeBundles
    axes::NamedTuple
end
import GaussletBases: primitive_set, stencil_matrix, _nested_axis_lengths, _nested_axis_pgdg
primitive_set(layer::RepresentedHartreeProxyLayer) = layer.primitives
stencil_matrix(layer::RepresentedHartreeProxyLayer) = layer.coefficients
_nested_axis_lengths(::RepresentedHartreeBundles) = (1, 1, 1)
_nested_axis_pgdg(bundle::RepresentedHartreeBundles, axis::Symbol) = bundle.axes[axis]

function oracle_field(pair_matrix, density)
    n = size(density, 1)
    field = zeros(Float64, n, n)
    for q in 1:n, p in 1:n, s in 1:n, r in 1:n
        pq = gaussian_coulomb_pair_index(p, q, n)
        rs = gaussian_coulomb_pair_index(r, s, n)
        field[p, q] += density[r, s] * pair_matrix[pq, rs]
    end
    return 0.5 .* (field .+ transpose(field))
end

function represented_fixture(nA)
    alpha = 0.7
    primitive = Gaussian(center = 0.0, width = inv(sqrt(2alpha)))
    normalization = GaussletBases.GaussianAnalyticIntegrals.
        polynomial_gaussian_shell_prefactor(alpha, 0)
    layer = RepresentedHartreeProxyLayer(
        PrimitiveSet1D([primitive]), reshape([normalization], 1, 1))
    proxy = (; x = layer, y = layer, z = layer, ncart = 1)
    axis = (; overlap = ones(1, 1), base_layer = layer)
    bundles = RepresentedHartreeBundles((; x = axis, y = axis, z = axis))
    block = CFR.CartesianTerminalBasisBlock(
        :identity, [1], [(1, 1, 1)], nothing, 1:1)
    basis = CFR.CartesianTerminalBasisRealization([block], 1, 0.0)
    centers = nA == 1 ? [0.65] : collect(range(-0.8, 0.8; length = nA))
    orbitals = [CartesianGaussianShellOrbitalRepresentation3D(
        "a$(index)", (0, 0, 0), (center, 0.15 * index, -0.1 * index),
        [0.45 + 0.12 * index], [1.0],
        :axiswise_normalized_cartesian_gaussian) for (index, center) in pairs(centers)]
    supplement = CartesianGaussianShellSupplementRepresentation3D(
        :synthetic, orbitals, (; fixture = :represented_hartree))
    donor = CFR._r3a_qw_supplement(supplement)
    parent_overlap = CGR.gaussian_non_nuclear_overlap_blocks(proxy, donor).ga.overlap
    X = CFR._r3a_project_parent_ga(basis, parent_overlap)
    S_AA = GaussletBases._cartesian_supplement_cross_overlap(supplement, supplement)
    residual = CRG.build_residual_gaussian_basis(1, X, S_AA,
        [orbital.label for orbital in orbitals], [orbital.center for orbital in orbitals],
        collect(1:nA); residual_occupation_cutoff = 1.0e-12)
    final_dimension = 1 + residual.residual_dimension
    alpha_final = reshape(normalize(collect(1.0:final_dimension)), :, 1)
    beta_final = reshape(normalize(reverse(collect(1.0:final_dimension))), :, 1)
    density = CRD.represented_mixed_density(basis, bundles, supplement, residual,
        alpha_final, [0.7], beta_final, [0.4]; proxy)
    expansion = CoulombGaussianExpansion([0.65, 0.25], [0.35, 1.1];
        del = 0.0, s = 0.0, c = 0.0, maxu = 0.0)
    potential = CRD.direct_represented_hartree_potential(
        density, expansion; proxy)
    field = CRD.represented_hartree_field(potential)
    g_orbital = CartesianGaussianShellOrbitalRepresentation3D(
        "g", (0, 0, 0), (0.0, 0.0, 0.0), [alpha], [1.0],
        :axiswise_normalized_cartesian_gaussian)
    mixed_orbitals = vcat([g_orbital], orbitals)
    pair_matrix = gaussian_coulomb_pair_matrix(
        mixed_orbitals; expansion, max_orbitals = length(mixed_orbitals))
    mixed_alpha = vcat(density.C_G_alpha, density.C_A_alpha)
    mixed_beta = vcat(density.C_G_beta, density.C_A_beta)
    P = 0.7 .* (mixed_alpha * transpose(mixed_alpha)) +
        0.4 .* (mixed_beta * transpose(mixed_beta))
    expected_raw = oracle_field(pair_matrix, P)
    expected_native = CRG.transform_augmented_operator(
        expected_raw[1:1, 1:1], expected_raw[1:1, 2:end],
        expected_raw[2:end, 2:end], residual)
    return (; density, potential, field, pair_matrix, P, expected_raw, expected_native,
        proxy, supplement, residual, expansion)
end

function represented_state_field(fixture, coefficients)
    density = CRD.represented_mixed_density(fixture.density.basis, fixture.density.bundles,
        fixture.supplement, fixture.residual, reshape(coefficients, :, 1), [1.0],
        zeros(length(coefficients), 0), Float64[]; proxy = fixture.proxy)
    return CRD.represented_hartree_field(CRD.direct_represented_hartree_potential(
        density, fixture.expansion; proxy = fixture.proxy))
end

@testset "represented mixed-density direct Hartree" begin
    small = represented_fixture(1)
    large = represented_fixture(2)
    for fixture in (small, large)
        density, field = fixture.density, fixture.field
        nA = length(fixture.supplement.orbitals)
        @test density.diagnostics.alpha_charge ≈ 0.7 atol = 1.0e-12
        @test density.diagnostics.beta_charge ≈ 0.4 atol = 1.0e-12
        @test density.diagnostics.total_charge ≈ 1.1 atol = 1.0e-12
        @test density.diagnostics.alpha_gram_error <= 1.0e-12
        @test density.diagnostics.beta_gram_error <= 1.0e-12
        @test density.diagnostics.residual_cross_error <= 1.0e-12
        @test density.diagnostics.residual_metric_error <= 1.0e-12
        @test field.GG ≈ fixture.expected_raw[1:1, 1:1] atol = 2.0e-11
        @test field.GA ≈ fixture.expected_raw[1:1, 2:(nA + 1)] atol = 2.0e-11
        @test field.AA ≈ fixture.expected_raw[2:(nA + 1), 2:(nA + 1)] atol = 2.0e-11
        @test field.native ≈ fixture.expected_native atol = 2.0e-11
        @test norm(field.native - transpose(field.native), Inf) <= 1.0e-12
        @test field.no_half_E0 ≈ sum(fixture.P .* fixture.expected_raw) atol = 2.0e-11
        @test field.hartree_energy == 0.5 * field.no_half_E0
        mixed_only = copy(fixture.P)
        mixed_only[diagind(mixed_only)] .= 0.0
        @test norm(oracle_field(fixture.pair_matrix, mixed_only), Inf) > 1.0e-4
        n = size(field.native, 1)
        delta = Diagonal(vcat([0.2, -0.2], zeros(max(0, n - 2)))) |> Matrix
        mixed_transform = [ones(1, 1) fixture.residual.T_G;
            zeros(nA, 1) fixture.residual.T_A]
        delta_mixed = mixed_transform * delta * transpose(mixed_transform)
        energy(P) = 0.5 * sum(P .* oracle_field(fixture.pair_matrix, P))
        epsilon = 1.0e-5
        derivative = (energy(fixture.P + epsilon * delta_mixed) -
            energy(fixture.P - epsilon * delta_mixed)) / (2epsilon)
        @test derivative ≈ sum(delta .* field.native) atol = 2.0e-10
        @test tr(delta) == 0.0
    end
    sresource = small.potential.action.diagnostics
    lresource = large.potential.action.diagnostics
    @test sresource.retained_one_dimensional_cache_entries == lresource.retained_one_dimensional_cache_entries == 0
    @test lresource.source_component_count > sresource.source_component_count
    @test lresource.source_block_pair_count == div(
        lresource.source_component_count * (lresource.source_component_count + 1), 2)
    @test all(>(0), lresource.primitive_scratch_entries)
    @test lresource.peak_additional_workspace_bytes <
        8 * lresource.source_block_pair_count^2 + lresource.output_bytes
    @test_throws DimensionMismatch CRD.represented_mixed_density(
        small.density.basis, small.density.bundles, small.supplement, small.residual,
        zeros(1, 1), [1.0], small.density.final_beta,
        small.density.occupations_beta; proxy = small.proxy)
    altered = copy(small.density.final_alpha); altered[1, 1] *= 1.1
    @test_throws ArgumentError CRD.represented_mixed_density(
        small.density.basis, small.density.bundles, small.supplement, small.residual,
        altered, [0.7], small.density.final_beta,
        small.density.occupations_beta; proxy = small.proxy)
    @test_throws ArgumentError CRD.represented_mixed_density(small.density.basis,
        small.density.bundles, small.supplement, small.residual, small.density.final_alpha,
        [-0.1], small.density.final_beta, small.density.occupations_beta; proxy = small.proxy)
    saved = small.residual.T_G[1, 1]; small.residual.T_G[1, 1] += 1.0e-3
    @test_throws ArgumentError CRD.represented_mixed_density(small.density.basis,
        small.density.bundles, small.supplement, small.residual, small.density.final_alpha,
        [0.7], small.density.final_beta, [0.4]; proxy = small.proxy)
    small.residual.T_G[1, 1] = saved
    saved = small.expansion.coefficients[1]; small.expansion.coefficients[1] += 1.0e-3
    @test_throws ArgumentError CRD.represented_hartree_field(small.potential)
    small.expansion.coefficients[1] = saved
    saved = small.density.C_G_alpha[1, 1]; small.density.C_G_alpha[1, 1] += 1.0e-3
    @test_throws ArgumentError CRD.represented_hartree_field(small.potential)
    small.density.C_G_alpha[1, 1] = saved
    saved = small.supplement.orbitals[1].coefficients[1]; small.supplement.orbitals[1].coefficients[1] *= 1.1
    @test_throws ArgumentError CRD.represented_mixed_density(small.density.basis,
        small.density.bundles, small.supplement, small.residual, small.density.final_alpha,
        [0.7], small.density.final_beta, [0.4]; proxy = small.proxy)
    small.supplement.orbitals[1].coefficients[1] = saved
    transform = [ones(1, 1) small.residual.T_G; zeros(1, 1) small.residual.T_A]
    aa_field = represented_state_field(small, transform \ [0.0, 1.0])
    donor = CFR._r3a_qw_supplement(small.supplement)
    atomic_gg = CGR.atomic_reference_hartree_gg_block(small.density.basis,
        small.density.bundles, small.supplement, reshape([1.0], 1, 1); expansion = small.expansion)
    atomic_ga_aa = CGR.atomic_reference_hartree_ga_aa_blocks(small.proxy, donor,
        small.supplement, reshape([1.0], 1, 1); expansion = small.expansion)
    @test aa_field.GG ≈ atomic_gg.GG atol = 2.0e-11
    @test aa_field.GA ≈ atomic_ga_aa.GA atol = 2.0e-11
    @test aa_field.AA ≈ atomic_ga_aa.AA atol = 2.0e-11
    @test norm(represented_state_field(small, [1.0, 0.0]).GG, Inf) > 1.0e-4
end
