using Test
using Logging
using LinearAlgebra
using JLD2
using SparseArrays

using GaussletBases

const _PROJECT_ROOT = dirname(@__DIR__)
const _RUN_SLOW_TESTS = get(ENV, "GAUSSLETBASES_SLOW_TESTS", "0") == "1"
const _FIXTURE_CACHE = Dict{Symbol,Any}()
const _TEST_GROUP_ENV = strip(get(ENV, "GAUSSLETBASES_TEST_GROUPS", "all"))
const _AVAILABLE_TEST_GROUPS = (
    :radial,
    :core,
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
                     angular_kernels[level_index][left_channel, left_source, right_source, right_channel] *
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

if _test_group_enabled(:core)
    include(joinpath(@__DIR__, "core", "runtests.jl"))
end

if _test_group_enabled(:angular)
    include(joinpath(@__DIR__, "angular", "runtests.jl"))
end

if _test_group_enabled(:ida)
    include(joinpath(@__DIR__, "ida", "runtests.jl"))
end


if _test_group_enabled(:misc)
    include(joinpath(@__DIR__, "misc", "runtests.jl"))
end

if _test_group_enabled(:docs)
    include(joinpath(@__DIR__, "docs", "runtests.jl"))
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
