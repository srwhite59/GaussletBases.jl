using Test
using LinearAlgebra

using GaussletBases

const _PROJECT_ROOT = dirname(@__DIR__)
const _RUN_SLOW_TESTS = get(ENV, "GAUSSLETBASES_SLOW_TESTS", "0") == "1"
const _FIXTURE_CACHE = Dict{Symbol,Any}()

_cached_fixture(key::Symbol, builder::Function) = get!(_FIXTURE_CACHE, key) do
    builder()
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

if _RUN_SLOW_TESTS
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

@testset "Documentation consistency" begin
    design = read(joinpath(_PROJECT_ROOT, "DESIGN.md"), String)
    readme = read(joinpath(_PROJECT_ROOT, "README.md"), String)
    atomic_setup = read(joinpath(_PROJECT_ROOT, "docs", "recommended_atomic_setup.md"), String)
    radial_workflow = read(joinpath(_PROJECT_ROOT, "docs", "first_radial_workflow.md"), String)
    architecture = read(joinpath(_PROJECT_ROOT, "docs", "architecture.md"), String)
    primitive_layer_note = read(joinpath(_PROJECT_ROOT, "docs", "intermediate_primitive_layer.md"), String)
    radial_primitive_note = read(joinpath(_PROJECT_ROOT, "docs", "radial_primitive_operator_layer.md"), String)
    ordinary_coulomb_note = read(joinpath(_PROJECT_ROOT, "docs", "ordinary_coulomb_expansion_path.md"), String)
    mapped_ordinary_note = read(joinpath(_PROJECT_ROOT, "docs", "mapped_ordinary_basis.md"), String)
    short_gaucoulomb_note = read(joinpath(_PROJECT_ROOT, "docs", "short_gaucoulomb_backend.md"), String)
    ordinary_pgdg_note = read(joinpath(_PROJECT_ROOT, "docs", "ordinary_pgdg_decision.md"), String)
    ordinary_pgdg_comx_note = read(joinpath(_PROJECT_ROOT, "docs", "ordinary_pgdg_comx.md"), String)
    ordinary_pgdg_refinement_note = read(joinpath(_PROJECT_ROOT, "docs", "ordinary_pgdg_proxy_refinement.md"), String)
    atomic_ylm_note = read(joinpath(_PROJECT_ROOT, "docs", "atomic_ylm_layer.md"), String)
    atomic_ida_note = read(joinpath(_PROJECT_ROOT, "docs", "atomic_ida_layer.md"), String)
    atomic_two_electron_note = read(joinpath(_PROJECT_ROOT, "docs", "atomic_ida_two_electron.md"), String)
    atomic_direct_note = read(joinpath(_PROJECT_ROOT, "docs", "atomic_ida_direct.md"), String)
    atomic_exchange_note = read(joinpath(_PROJECT_ROOT, "docs", "atomic_ida_exchange.md"), String)
    atomic_fock_note = read(joinpath(_PROJECT_ROOT, "docs", "atomic_ida_fock.md"), String)
    atomic_spin_fock_note = read(joinpath(_PROJECT_ROOT, "docs", "atomic_ida_spin_fock.md"), String)
    atomic_uhf_note = read(joinpath(_PROJECT_ROOT, "docs", "atomic_ida_uhf.md"), String)
    gaunt_backend_note = read(joinpath(_PROJECT_ROOT, "docs", "gaunt_backend_note.md"), String)
    example_guide = read(joinpath(_PROJECT_ROOT, "docs", "example_guide.md"), String)
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
    @test occursin("atom-centered radial calculations", readme)
    @test occursin("primitive layers, contraction, partitions, and hierarchy", readme)
    @test occursin("Best first path through the repository", readme)
    @test occursin("examples/04_hydrogen_ground_state.jl", readme)
    @test occursin("examples/15_atomic_hydrogen_ylm.jl", readme)
    @test occursin("13_global_leaf_contraction.jl", readme)
    @test occursin("recommended_atomic_setup.md", readme)
    @test occursin("atomic_ylm_layer.md", readme)
    @test occursin("atomic_ida_uhf.md", readme)
    @test occursin("mapped_ordinary_basis.md", readme)
    @test occursin("first_radial_workflow.md", readme)
    @test occursin("example_guide.md", readme)
    @test occursin("architecture.md", readme)
    @test occursin("intermediate_primitive_layer.md", readme)
    @test occursin("terminology.md", readme)
    @test occursin("diag = basis_diagnostics(rb)", readme)
    @test occursin("ops = atomic_operators(rb, grid; Z = Z, lmax = 2)", readme)
    @test occursin("Pkg.add(url = \"https://github.com/srwhite59/GaussletBases.jl\")", readme)
    @test occursin("s = 0.2", atomic_setup)
    @test occursin("odd_even_kmax = 6", atomic_setup)
    @test occursin("accuracy = :veryhigh", atomic_setup)
    @test occursin("examples/04_hydrogen_ground_state.jl", radial_workflow)
    @test occursin("Quickstart", radial_workflow)
    @test occursin("quadrature_rmax = 30.0", radial_workflow)
    @test occursin("example_guide.md", radial_workflow)
    @test occursin("ordinary gausslets are the broad foundation", lowercase(architecture))
    @test occursin("radial gausslets are the current mature public-facing workflow", lowercase(architecture))
    @test occursin("first minimal atomic mean-field line", lowercase(architecture))
    @test occursin("it provides a small uhf kernel for the current atomic ida approximation", lowercase(architecture))
    @test occursin("global mapped layer plus local contraction", lowercase(architecture))
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
    @test occursin("GauntTables", gaunt_backend_note)
    @test occursin("public atomic story should remain the same", lowercase(gaunt_backend_note))
    @test occursin("src/atomic_ida.jl", gaunt_backend_note)
    @test occursin("13_global_leaf_contraction.jl", example_guide)
    @test occursin("15_atomic_hydrogen_ylm.jl", example_guide)
    @test occursin("23_cartesian_hydrogen_coulomb_expansion.jl", example_guide)
    @test occursin("24_mapped_cartesian_hydrogen.jl", example_guide)
    @test occursin("16_atomic_ida_ingredients.jl", example_guide)
    @test occursin("19_atomic_ida_direct.jl", example_guide)
    @test occursin("20_atomic_ida_exchange.jl", example_guide)
    @test occursin("21_atomic_ida_fock.jl", example_guide)
    @test occursin("22_atomic_ida_uhf.jl", example_guide)
    @test occursin("17_atomic_ida_two_electron.jl", example_guide)
    @test occursin("18_atomic_ida_two_electron_lanczos.jl", example_guide)
    @test occursin("docs/atomic_ida_direct.md", example_guide)
    @test occursin("docs/atomic_ida_exchange.md", example_guide)
    @test occursin("docs/atomic_ida_fock.md", example_guide)
    @test occursin("docs/atomic_ida_spin_fock.md", example_guide)
    @test occursin("docs/atomic_ida_uhf.md", example_guide)
    @test occursin("docs/ordinary_coulomb_expansion_path.md", example_guide)
    @test occursin("docs/mapped_ordinary_basis.md", example_guide)
    @test occursin("small atomic ida / hf line", lowercase(example_guide))
    @test occursin("prototype", example_guide)
    @test occursin("radial atomic work", lowercase(example_guide))
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
    @test !occursin("Gausslets.jl", readme)
    @test occursin("Established public-facing path: radial calculations", status)
    @test occursin("First atomic angular step", status)
    @test occursin("First minimal atomic mean-field line", status)
    @test occursin("ordinary cartesian hydrogen branch", lowercase(status))
    @test occursin("a broad general hf/rhf/uhf workflow", lowercase(status))
    @test occursin("Prototype line", status)
    @test occursin("one global mapped primitive layer", status)
end

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
