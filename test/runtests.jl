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

include(joinpath(@__DIR__, "support", "radial_fixtures.jl"))

if _test_group_enabled(:radial)
    include(joinpath(@__DIR__, "radial", "runtests.jl"))
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
    include(joinpath(@__DIR__, "core", "runtests.jl"))
end

if _test_group_enabled(:nested)
    include(joinpath(@__DIR__, "nested", "runtests.jl"))
end

if _test_group_enabled(:ordinary)
    include(joinpath(@__DIR__, "ordinary", "runtests.jl"))
end

if _test_group_enabled(:diatomic)
    include(joinpath(@__DIR__, "diatomic", "runtests.jl"))
end




if _test_group_enabled(:angular)
    include(joinpath(@__DIR__, "angular", "runtests.jl"))
end

if _test_group_enabled(:ida)
    include(joinpath(@__DIR__, "ida", "runtests.jl"))
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
    @test oc
