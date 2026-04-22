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

    full_ops_localized = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :by_center,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
    )
    atom_a_ops_localized = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = [1.0, 0.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
    )
    atom_b_ops_localized = ordinary_cartesian_qiu_white_operators(
        basis;
        nuclear_charges = [0.0, 1.0],
        nuclear_term_storage = :total_only,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
    )

    @test full_ops_localized.gausslet_backend == :pgdg_localized_experimental
    @test full_ops_localized.nuclear_term_storage == :by_center
    @test !isnothing(full_ops_localized.kinetic_one_body)
    @test !isnothing(full_ops_localized.nuclear_one_body_by_center)
    @test length(full_ops_localized.nuclear_one_body_by_center) == 2
    @test assembled_one_body_hamiltonian(full_ops_localized) ≈
          full_ops_localized.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test assembled_one_body_hamiltonian(full_ops_localized; nuclear_charges = [1.0, 0.0]) ≈
          atom_a_ops_localized.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
    @test assembled_one_body_hamiltonian(full_ops_localized; nuclear_charges = [0.0, 1.0]) ≈
          atom_b_ops_localized.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10

    (
        _nested_basis,
        _nested_parent_ops,
        _nested_parent_check,
        expansion,
        _nested_source,
        fixed_block,
        _parent_modes,
        _parent_ground,
        _projected,
        _projected_vee,
        _capture,
        _projected_energy,
    ) = _bond_aligned_diatomic_nested_fixed_block_fixture(; bond_length = 1.4)

    full_nested_ops_localized = ordinary_cartesian_qiu_white_operators(
        fixed_block;
        nuclear_charges = [1.0, 1.0],
        nuclear_term_storage = :by_center,
        expansion = expansion,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
    )
    atom_a_nested_ops_localized = ordinary_cartesian_qiu_white_operators(
        fixed_block;
        nuclear_charges = [1.0, 0.0],
        nuclear_term_storage = :total_only,
        expansion = expansion,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
    )
    atom_b_nested_ops_localized = ordinary_cartesian_qiu_white_operators(
        fixed_block;
        nuclear_charges = [0.0, 1.0],
        nuclear_term_storage = :total_only,
        expansion = expansion,
        interaction_treatment = :ggt_nearest,
        gausslet_backend = :pgdg_localized_experimental,
    )

    @test full_nested_ops_localized.gausslet_backend == :pgdg_localized_experimental
    @test full_nested_ops_localized.nuclear_term_storage == :by_center
    @test !isnothing(full_nested_ops_localized.kinetic_one_body)
    @test !isnothing(full_nested_ops_localized.nuclear_one_body_by_center)
    @test length(full_nested_ops_localized.nuclear_one_body_by_center) == 2
    localized_nested_bundles = GaussletBases._qwrg_bond_aligned_axis_bundles(
        fixed_block.parent_basis,
        expansion;
        gausslet_backend = :pgdg_localized_experimental,
    )
    localized_nested_factorized_basis =
        GaussletBases._nested_factorized_parent_basis(fixed_block)
    direct_contracted_nested_nuclear =
        GaussletBases._qwrg_bond_aligned_direct_contracted_nuclear_one_body_by_center(
            fixed_block.parent_basis,
            localized_nested_factorized_basis,
            localized_nested_bundles.bundle_x,
            localized_nested_bundles.bundle_y,
            localized_nested_bundles.bundle_z,
            expansion,
        )
    reference_nested_parent_nuclear = GaussletBases._qwrg_diatomic_nuclear_one_body_by_center(
        fixed_block.parent_basis,
        localized_nested_bundles.bundle_x,
        localized_nested_bundles.bundle_y,
        localized_nested_bundles.bundle_z,
        expansion,
    )
    reference_nested_contracted_nuclear = [
        GaussletBases._qwrg_contract_parent_symmetric_matrix(
            fixed_block.coefficient_matrix,
            matrix,
        ) for matrix in reference_nested_parent_nuclear
    ]
    @test length(direct_contracted_nested_nuclear) == length(reference_nested_contracted_nuclear)
    for nucleus_index in eachindex(direct_contracted_nested_nuclear, reference_nested_contracted_nuclear)
        @test direct_contracted_nested_nuclear[nucleus_index] ≈
              reference_nested_contracted_nuclear[nucleus_index] atol = 1.0e-10 rtol = 1.0e-10
    end
    @test assembled_one_body_hamiltonian(full_nested_ops_localized) ≈
          full_nested_ops_localized.one_body_hamiltonian atol = 1.0e-12 rtol = 1.0e-12
    @test assembled_one_body_hamiltonian(full_nested_ops_localized; nuclear_charges = [1.0, 0.0]) ≈
          atom_a_nested_ops_localized.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
    @test assembled_one_body_hamiltonian(full_nested_ops_localized; nuclear_charges = [0.0, 1.0]) ≈
          atom_b_nested_ops_localized.one_body_hamiltonian atol = 1.0e-10 rtol = 1.0e-10
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
    @test length(nested_payload.shell_provenance) == length(source.shared_shell_layers)
    @test nested_payload.shell_provenance[1].source_box == source.shared_shell_layers[1].provenance.source_box
    @test nested_payload.shell_provenance[1].next_inner_box ==
        source.shared_shell_layers[1].provenance.next_inner_box
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
        first_shell = source.shared_shell_layers[1].provenance
        first_payload_shell = source_payload.shell_provenance[1]
        first_shell_source_dims = string(
            length(first_shell.source_box[1]), "x",
            length(first_shell.source_box[2]), "x",
            length(first_shell.source_box[3]),
        )
        first_shell_next_inner_dims = string(
            length(first_shell.next_inner_box[1]), "x",
            length(first_shell.next_inner_box[2]), "x",
            length(first_shell.next_inner_box[3]),
        )

        @test all(isnothing(sequence.support_states) for sequence in source.child_sequences)
        @test Set(point.group_kind for point in source_payload.points) == expected_source_groups
        @test length(source_payload.points) == prod(length.(source.geometry.parent_box))
        @test length(source_payload.shell_provenance) == length(source.shared_shell_layers)
        @test first_shell.source_box == source.geometry.parent_box
        @test first_shell.source_point_count ==
            prod(length.(first_shell.source_box)) - prod(length.(first_shell.next_inner_box))
        @test first_shell.retained_fixed_count == length(source.sequence.layer_column_ranges[1])
        @test first_payload_shell.source_box == first_shell.source_box
        @test first_payload_shell.next_inner_box == first_shell.next_inner_box
        @test first_payload_shell.source_point_count == first_shell.source_point_count
        @test first_payload_shell.retained_fixed_count == first_shell.retained_fixed_count
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
        @test source_slice.shell_provenance == source_payload.shell_provenance

        mktemp() do path, io
            close(io)
            payload = write_bond_aligned_diatomic_points3d(path, source_payload)
            text = read(path, String)
            @test payload === source_payload
            @test occursin("# bond_axis = z", text)
            @test occursin("# point_count = $(length(source_payload.points))", text)
            @test occursin("# nucleus_count = 2", text)
            @test occursin("# columns = x y z role kind group_kind group_id label", text)
            @test occursin("# shell label=shared_shell_1", text)
            @test occursin("source_box=$(first_shell_source_dims)", text)
            @test occursin("source_points=$(first_shell.source_point_count)", text)
            @test occursin("retained_fixed_count=$(first_shell.retained_fixed_count)", text)
            @test occursin("next_inner_box=$(first_shell_next_inner_dims)", text)
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

        mktemp() do path, io
            close(io)
            slice = write_bond_aligned_diatomic_plane_projection(
                path,
                source_payload;
                plane_axis = :y,
                plane_value = 0.0,
                plane_tol = 1.0e-5,
            )
            text = read(path, String)
            @test slice.selected_count == source_slice.selected_count
            @test occursin("# shell label=shared_shell_1", text)
            @test occursin("source_box=$(first_shell_source_dims)", text)
            @test occursin("next_inner_box=$(first_shell_next_inner_dims)", text)
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
        @test occursin("# shared_shell layer=1 source_box=", text)
        @test occursin(" source_points=", text)
        @test occursin(" retained_fixed_count=", text)
        @test occursin(" next_inner_box=", text)
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
