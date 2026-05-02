function _experimental_high_order_identity_basis(
    count::Int;
    reference_spacing::Real = 1.0,
)
    return build_basis(MappedUniformBasisSpec(:G10;
        count = count,
        mapping = IdentityMapping(),
        reference_spacing = Float64(reference_spacing),
    ))
end

function _experimental_high_order_distorted_he_basis(
    count::Int;
    Z::Real = 2.0,
    d::Real = 0.2,
    tail_spacing::Real = 10.0,
    reference_spacing::Real = 1.0,
)
    return build_basis(MappedUniformBasisSpec(:G10;
        count = count,
        mapping = white_lindsey_atomic_mapping(
            Z = Z,
            d = d,
            tail_spacing = tail_spacing,
        ),
        reference_spacing = Float64(reference_spacing),
    ))
end

function _experimental_high_order_he_singlet_case(
    sides::AbstractVector{<:Integer};
    parent_side::Int = maximum(sides),
    backend::Symbol = :numerical_reference,
    reference_spacing::Real = 1.0,
    expansion::CoulombGaussianExpansion = coulomb_gaussian_expansion(doacc = false),
    Z::Real = 2.0,
    krylovdim::Int = 32,
    maxiter::Int = 32,
    tol::Real = 1.0e-8,
)
    basis = _experimental_high_order_identity_basis(
        parent_side;
        reference_spacing = reference_spacing,
    )
    stack = GaussletBases._experimental_high_order_doside_stack_3d(
        basis;
        backend = backend,
        sides = sides,
    )
    data = GaussletBases._experimental_high_order_doside_he_singlet_data(
        stack;
        expansion = expansion,
        Z = Z,
        krylovdim = krylovdim,
        maxiter = maxiter,
        tol = tol,
    )
    return (basis = basis, stack = stack, data = data)
end

@testset "Experimental high-order physical 1D polynomial blocks stay sane" begin
    cases = (
        (_experimental_high_order_identity_basis(11), 5, 5),
        (_experimental_high_order_distorted_he_basis(11), 7, 5),
        (_experimental_high_order_distorted_he_basis(11), 11, 5),
        (_experimental_high_order_distorted_he_basis(11), 10, 6),
    )

    for (basis, side, doside) in cases
        axis_data = GaussletBases._experimental_high_order_axis_data_1d(
            basis;
            backend = :numerical_reference,
        )
        physical = GaussletBases._experimental_high_order_physical_block_1d(
            axis_data,
            side;
            doside = doside,
        )
        current = GaussletBases._experimental_high_order_block_1d(
            axis_data,
            side;
            doside = doside,
        )

        block = physical.block
        diagnostics = physical.diagnostics
        parent_overlap = axis_data.overlap
        physical_coefficients = Matrix{Float64}(block.coefficients)
        current_coefficients = Matrix{Float64}(current.coefficients)

        @test size(block.local_coefficients, 2) == doside
        @test diagnostics.precleanup_overlap_spectrum.kept_rank == doside
        @test diagnostics.precleanup_overlap_spectrum.minimum_eigenvalue > 1.0e-12
        @test diagnostics.overlap_error < 1.0e-8
        @test all(isfinite, diagnostics.parent_polynomial_projection_errors)
        @test all(isfinite, diagnostics.parent_polynomial_capture_fractions)
        @test first(diagnostics.parent_polynomial_projection_errors) < 1.0e-12
        @test diagnostics.max_parent_polynomial_projection_error < 1.0e-12
        @test diagnostics.min_parent_polynomial_capture_fraction > 1.0 - 1.0e-12
        @test issorted(block.localized_centers)

        current_in_physical = current_coefficients - physical_coefficients * (
            transpose(physical_coefficients) * parent_overlap * current_coefficients
        )
        physical_in_current = physical_coefficients - current_coefficients * (
            transpose(current_coefficients) * parent_overlap * physical_coefficients
        )

        @test norm(transpose(current_in_physical) * parent_overlap * current_in_physical, Inf) < 1.0e-8
        @test norm(transpose(physical_in_current) * parent_overlap * physical_in_current, Inf) < 1.0e-8
    end
end

@testset "Experimental high-order physical 3D full blocks stay sane" begin
    cases = (
        (_experimental_high_order_identity_basis(11), 5, 5),
        (_experimental_high_order_distorted_he_basis(11), 7, 5),
        (_experimental_high_order_distorted_he_basis(11), 10, 6),
    )

    for (basis, side, doside) in cases
        axis_data = GaussletBases._experimental_high_order_axis_data_1d(
            basis;
            backend = :numerical_reference,
        )
        physical = GaussletBases._experimental_high_order_physical_full_block_3d(
            axis_data,
            side;
            doside = doside,
        )
        diagnostics = physical.diagnostics

        @test diagnostics.full_block_dimension == doside^3
        @test diagnostics.full_block_overlap_error < 1.0e-8
        @test diagnostics.full_block_overlap_spectrum.kept_rank == doside^3
        @test diagnostics.full_block_overlap_spectrum.minimum_eigenvalue > 1.0e-12
        @test diagnostics.current_in_physical_residual < 1.0e-8
        @test diagnostics.physical_in_current_residual < 1.0e-8
    end
end

@testset "Experimental high-order physical 3D shell selection stays sane" begin
    cases = (
        (_experimental_high_order_identity_basis(11), 5, 5, (faces = 54, edges = 36, corners = 8)),
        (_experimental_high_order_distorted_he_basis(11), 7, 5, (faces = 54, edges = 36, corners = 8)),
        (_experimental_high_order_distorted_he_basis(11), 10, 6, (faces = 96, edges = 48, corners = 8)),
    )

    for (basis, side, doside, expected_shell_kinds) in cases
        axis_data = GaussletBases._experimental_high_order_axis_data_1d(
            basis;
            backend = :numerical_reference,
        )
        physical = GaussletBases._experimental_high_order_physical_shell_3d(
            axis_data,
            side;
            doside = doside,
        )
        diagnostics = physical.diagnostics

        @test diagnostics.shell_dimension == GaussletBases._experimental_high_order_expected_shell_dimension(doside)
        @test diagnostics.expected_shell_dimension == GaussletBases._experimental_high_order_expected_shell_dimension(doside)
        @test diagnostics.shell_overlap_error < 1.0e-8
        @test diagnostics.shell_overlap_spectrum.kept_rank == diagnostics.shell_dimension
        @test diagnostics.shell_overlap_spectrum.minimum_eigenvalue > 1.0e-12
        @test diagnostics.shell_kind_counts == expected_shell_kinds
        @test diagnostics.current_shell_in_physical_residual < 1.0e-8
        @test diagnostics.physical_shell_in_current_residual < 1.0e-8
    end
end

@testset "Experimental high-order reduced one-body contraction stays sane" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    identity_data = GaussletBases._experimental_high_order_physical_reduced_one_body_data(
        _experimental_high_order_identity_basis(11),
        5;
        doside = 5,
        backend = :numerical_reference,
        expansion = expansion,
        Z = 2.0,
        direct_comparison = :always,
    )
    identity_diagnostics = identity_data.diagnostics

    @test identity_diagnostics.direct_comparison_performed
    @test identity_diagnostics.full_overlap_error < 1.0e-8
    @test identity_diagnostics.full_hamiltonian_error < 1.0e-8
    @test identity_diagnostics.shell_overlap_error < 1.0e-8
    @test identity_diagnostics.shell_hamiltonian_error < 1.0e-8
    @test identity_diagnostics.full_summary.overlap_error < 1.0e-8
    @test identity_diagnostics.shell_summary.overlap_error < 1.0e-8
    @test abs(identity_diagnostics.full_summary.ground_energy - identity_diagnostics.direct_full_summary.ground_energy) < 1.0e-8
    @test abs(identity_diagnostics.shell_summary.ground_energy - identity_diagnostics.direct_shell_summary.ground_energy) < 1.0e-8

    distorted_cases = (
        (_experimental_high_order_distorted_he_basis(5), 5, 5),
    )

    for (basis, side, doside) in distorted_cases
        data = GaussletBases._experimental_high_order_physical_reduced_one_body_data(
            basis,
            side;
            doside = doside,
            backend = :numerical_reference,
            expansion = expansion,
            Z = 2.0,
            direct_comparison = :never,
        )
        diagnostics = data.diagnostics

        @test !diagnostics.direct_comparison_performed
        @test isnothing(data.parent_data.parent_overlap)
        @test isnothing(data.parent_data.parent_hamiltonian)
        @test isnothing(data.direct_full)
        @test isnothing(data.direct_shell)
        @test isnan(diagnostics.full_overlap_error)
        @test isnan(diagnostics.full_hamiltonian_error)
        @test isnan(diagnostics.shell_overlap_error)
        @test isnan(diagnostics.shell_hamiltonian_error)
        @test diagnostics.full_summary.overlap_error < 1.0e-8
        @test diagnostics.shell_summary.overlap_error < 1.0e-8
        @test diagnostics.full_summary.overlap_spectrum.minimum_eigenvalue > 1.0e-12
        @test diagnostics.shell_summary.overlap_spectrum.minimum_eigenvalue > 1.0e-12
        @test isfinite(diagnostics.full_summary.ground_energy)
        @test isfinite(diagnostics.shell_summary.ground_energy)
        @test diagnostics.full_summary.ground_energy <= diagnostics.shell_summary.ground_energy + 1.0e-10
        @test isnothing(diagnostics.direct_full_summary)
        @test isnothing(diagnostics.direct_shell_summary)
    end
end

@testset "Experimental high-order axis one-body cache reuses repeated requests" begin
    basis = _experimental_high_order_identity_basis(5)
    axis_data = GaussletBases._experimental_high_order_axis_data_1d(
        basis;
        backend = :numerical_reference,
    )
    expansion = coulomb_gaussian_expansion(doacc = false)

    first_one_body = GaussletBases._experimental_high_order_axis_one_body_1d(
        axis_data;
        exponents = expansion.exponents,
        center = 0.0,
    )
    second_one_body = GaussletBases._experimental_high_order_axis_one_body_1d(
        axis_data;
        exponents = expansion.exponents,
        center = 0.0,
    )

    @test first_one_body === second_one_body
    @test length(axis_data.one_body_cache) == 1

    parent_data = GaussletBases._experimental_high_order_parent_one_body_data(
        basis;
        axis_data = axis_data,
        backend = :numerical_reference,
        expansion = expansion,
        Z = 2.0,
        include_parent_projection_data = false,
    )
    @test parent_data.axis_data === axis_data

    reused_stack = GaussletBases._experimental_high_order_doside_stack_3d(
        basis;
        axis_data = axis_data,
        backend = :numerical_reference,
        doside = 5,
        sides = [5],
    )
    fresh_stack = GaussletBases._experimental_high_order_doside_stack_3d(
        basis;
        backend = :numerical_reference,
        doside = 5,
        sides = [5],
    )
    @test reused_stack.coefficient_matrix ≈ fresh_stack.coefficient_matrix atol = 1.0e-12 rtol = 1.0e-12
end

@testset "Experimental high-order PGDG axis one-body consumption matches ordinary backend" begin
    basis = _experimental_high_order_distorted_he_basis(5)
    axis_data = GaussletBases._experimental_high_order_axis_data_1d(
        basis;
        backend = :pgdg_localized_experimental,
    )
    expansion = coulomb_gaussian_expansion(doacc = false)

    ordinary = mapped_ordinary_one_body_operators(
        basis;
        exponents = expansion.exponents,
        center = 0.0,
        backend = :pgdg_localized_experimental,
    )
    cached = GaussletBases._experimental_high_order_axis_one_body_1d(
        axis_data;
        exponents = expansion.exponents,
        center = 0.0,
    )

    @test cached.backend == :pgdg_localized_experimental
    @test cached.overlap ≈ ordinary.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test cached.kinetic ≈ ordinary.kinetic atol = 1.0e-12 rtol = 1.0e-12
    @test length(cached.gaussian_factors) == length(ordinary.gaussian_factors)
    for (cached_factor, ordinary_factor) in zip(cached.gaussian_factors, ordinary.gaussian_factors)
        @test cached_factor ≈ ordinary_factor atol = 1.0e-12 rtol = 1.0e-12
    end

    block = GaussletBases._experimental_high_order_physical_block_1d(
        axis_data,
        5;
        doside = 5,
    ).block
    cached_reduced = GaussletBases._experimental_high_order_contract_one_body_1d(
        cached,
        block.coefficients,
    )
    ordinary_reduced = GaussletBases._experimental_high_order_contract_one_body_1d(
        ordinary,
        block.coefficients,
    )

    @test cached_reduced.overlap ≈ ordinary_reduced.overlap atol = 1.0e-12 rtol = 1.0e-12
    @test cached_reduced.kinetic ≈ ordinary_reduced.kinetic atol = 1.0e-12 rtol = 1.0e-12
    @test length(cached_reduced.gaussian_factors) == length(ordinary_reduced.gaussian_factors)
    for (cached_factor, ordinary_factor) in zip(cached_reduced.gaussian_factors, ordinary_reduced.gaussian_factors)
        @test cached_factor ≈ ordinary_factor atol = 1.0e-12 rtol = 1.0e-12
    end
end

@testset "Experimental high-order PGDG reduced one-body matches dense direct projection on distorted case" begin
    basis = _experimental_high_order_distorted_he_basis(5)
    expansion = coulomb_gaussian_expansion(doacc = false)
    data = GaussletBases._experimental_high_order_physical_reduced_one_body_data(
        basis,
        5;
        doside = 5,
        backend = :pgdg_localized_experimental,
        expansion = expansion,
        Z = 2.0,
        direct_comparison = :always,
    )
    diagnostics = data.diagnostics

    @test diagnostics.direct_comparison_performed
    @test diagnostics.full_overlap_error < 1.0e-8
    @test diagnostics.full_hamiltonian_error < 1.0e-8
    @test diagnostics.shell_overlap_error < 1.0e-8
    @test diagnostics.shell_hamiltonian_error < 1.0e-8
    @test !isnothing(data.parent_data.parent_overlap)
    @test !isnothing(data.parent_data.parent_hamiltonian)
    @test !isnothing(data.direct_full)
    @test !isnothing(data.direct_shell)
    @test abs(diagnostics.full_summary.ground_energy - diagnostics.direct_full_summary.ground_energy) < 1.0e-8
    @test abs(diagnostics.shell_summary.ground_energy - diagnostics.direct_shell_summary.ground_energy) < 1.0e-8
end

@testset "Experimental high-order doside stack core dimensions" begin
    for (sides, expected_dimension) in (
        ([5], 125),
        ([5, 7], 223),
        ([5, 7, 9, 11], 419),
    )
        basis = _experimental_high_order_identity_basis(maximum(sides))
        stack = GaussletBases._experimental_high_order_doside_stack_3d(
            basis;
            backend = :numerical_reference,
            sides = sides,
        )
        @test stack isa GaussletBases.ExperimentalHighOrderDosideStack3D
        @test stack.backend == :numerical_reference
        @test stack.doside == 5
        @test stack.sides == Int[sides...]
        @test stack.parent_side == maximum(sides)
        @test size(stack.coefficient_matrix, 2) == expected_dimension
        @test stack.diagnostics.stack_dimension == expected_dimension
    end
end

@testset "Experimental high-order doside shell-kind counts" begin
    basis = _experimental_high_order_identity_basis(11)
    stack = GaussletBases._experimental_high_order_doside_stack_3d(
        basis;
        backend = :numerical_reference,
        sides = [5, 7, 9, 11],
    )

    @test length(stack.shell_layers) == 3
    @test stack.diagnostics.shell_dimensions == [98, 98, 98]
    for shell in stack.shell_layers
        @test length(shell.shell_labels) == 98
        @test shell.shell_kind_counts == (faces = 54, edges = 36, corners = 8)
    end
end

@testset "Experimental high-order even-doside stack dimensions and shell counts" begin
    for (parent_side, doside, sides, expected_dimension, expected_shell_dimension, expected_shell_kinds) in (
        (11, 4, [4, 6, 8, 10], 232, 56, (faces = 24, edges = 24, corners = 8)),
        (11, 6, [6, 8, 10], 520, 152, (faces = 96, edges = 48, corners = 8)),
    )
        basis = _experimental_high_order_identity_basis(parent_side)
        stack = GaussletBases._experimental_high_order_doside_stack_3d(
            basis;
            backend = :numerical_reference,
            doside = doside,
            sides = sides,
        )

        @test stack.doside == doside
        @test stack.sides == Int[sides...]
        @test size(stack.coefficient_matrix, 2) == expected_dimension
        @test stack.diagnostics.stack_dimension == expected_dimension
        @test all(==(expected_shell_dimension), stack.diagnostics.shell_dimensions)
        @test all(==(expected_shell_kinds), [shell.shell_kind_counts for shell in stack.shell_layers])
    end
end

@testset "Experimental high-order even-doside matched count ladders" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    basis = _experimental_high_order_distorted_he_basis(11)
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
    )

    lower4 = GaussletBases._experimental_high_order_lower_route_data(
        bundle;
        expansion = expansion,
        outer_sides = [4, 6, 8, 10],
        comparator_nside = 4,
    )
    lower6 = GaussletBases._experimental_high_order_lower_route_data(
        bundle;
        expansion = expansion,
        outer_sides = [6, 8, 10],
        comparator_nside = 6,
    )

    high4_counts = Int[]
    for sides in ([4], [4, 6], [4, 6, 8], [4, 6, 8, 10])
        stack = GaussletBases._experimental_high_order_doside_stack_3d(
            basis;
            backend = :numerical_reference,
            doside = 4,
            sides = sides,
        )
        push!(high4_counts, size(stack.coefficient_matrix, 2))
    end

    high6_counts = Int[]
    for sides in ([6], [6, 8], [6, 8, 10])
        stack = GaussletBases._experimental_high_order_doside_stack_3d(
            basis;
            backend = :numerical_reference,
            doside = 6,
            sides = sides,
        )
        push!(high6_counts, size(stack.coefficient_matrix, 2))
    end

    @test [row.function_count for row in lower4] == high4_counts == [64, 120, 176, 232]
    @test [row.function_count for row in lower6] == high6_counts == [216, 368, 520]
end

@testset "Experimental high-order even-doside distorted He+ data stay sane" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    basis = _experimental_high_order_distorted_he_basis(11)

    for (doside, sides) in (
        (4, [4, 6, 8, 10]),
        (6, [6, 8, 10]),
    )
        stack = GaussletBases._experimental_high_order_doside_stack_3d(
            basis;
            backend = :numerical_reference,
            doside = doside,
            sides = sides,
        )
        data = GaussletBases._experimental_high_order_doside_heplus_data(
            stack;
            expansion = expansion,
            Z = 2.0,
        )

        @test stack.diagnostics.parent_mapping_family == :white_lindsey_atomic_he_d0p2
        @test stack.diagnostics.overlap_error < 1.0e-8
        @test stack.diagnostics.overlap_spectrum.minimum_eigenvalue > 1.0 - 1.0e-8
        @test stack.diagnostics.overlap_spectrum.maximum_eigenvalue < 1.0 + 1.0e-8
        @test isfinite(data.ground_energy)
        @test data.overlap_error < 1.0e-8
        @test data.diagnostics.projected_overlap_spectrum.minimum_eigenvalue > 1.0 - 1.0e-8
        @test data.diagnostics.projected_overlap_spectrum.maximum_eigenvalue < 1.0 + 1.0e-8
    end
end

@testset "Experimental high-order doside stack orthogonality and span" begin
    for sides in ([5], [5, 7], [5, 7, 9, 11])
        basis = _experimental_high_order_identity_basis(maximum(sides))
        axis_data = GaussletBases._experimental_high_order_axis_data_1d(
            basis;
            backend = :numerical_reference,
        )
        parent_overlap = GaussletBases._experimental_high_order_parent_overlap_3d(axis_data)
        stack = GaussletBases._experimental_high_order_doside_stack_3d(
            basis;
            backend = :numerical_reference,
            sides = sides,
        )
        coefficients = Matrix{Float64}(stack.coefficient_matrix)
        @test norm(transpose(coefficients) * parent_overlap * coefficients - I, Inf) < 1.0e-8
        @test stack.diagnostics.overlap_error < 1.0e-8
        @test stack.diagnostics.contracted_weights_finite
        @test all(isfinite, stack.contracted_weights)
    end

    basis = _experimental_high_order_identity_basis(11)
    axis_data = GaussletBases._experimental_high_order_axis_data_1d(
        basis;
        backend = :numerical_reference,
    )
    parent_overlap = GaussletBases._experimental_high_order_parent_overlap_3d(axis_data)
    stack = GaussletBases._experimental_high_order_doside_stack_3d(
        basis;
        backend = :numerical_reference,
        sides = [5, 7, 9, 11],
    )
    union_coefficients = GaussletBases._experimental_high_order_full_block_union_coefficients(
        axis_data,
        [5, 7, 9, 11];
        doside = 5,
    )
    stack_coefficients = Matrix{Float64}(stack.coefficient_matrix)
    union_matrix = Matrix{Float64}(union_coefficients)
    residual = union_matrix - stack_coefficients * (transpose(stack_coefficients) * parent_overlap * union_matrix)
    residual_metric = transpose(residual) * parent_overlap * residual
    union_metric = transpose(union_matrix) * parent_overlap * union_matrix
    @test count(>(1.0e-8), svdvals(union_metric)) == 419
    @test norm(residual_metric, Inf) < 1.0e-8
end

@testset "Experimental high-order doside projected overlap and He+ ladder" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    energies = Float64[]

    for sides in ([5], [5, 7], [5, 7, 9], [5, 7, 9, 11])
        basis = _experimental_high_order_identity_basis(maximum(sides))
        stack = GaussletBases._experimental_high_order_doside_stack_3d(
            basis;
            backend = :numerical_reference,
            sides = sides,
        )
        data = GaussletBases._experimental_high_order_doside_heplus_data(
            stack;
            expansion = expansion,
            Z = 2.0,
        )
        @test norm(data.projected_overlap - I, Inf) < 1.0e-8
        @test data.overlap_error < 1.0e-8
        @test isfinite(data.ground_energy)
        @test all(isfinite, data.orbital_energies)
        push!(energies, data.ground_energy)
    end

    @test all(diff(energies) .<= 1.0e-8)
end

@testset "Experimental high-order doside He+ matches orthonormalized full-union reference" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    basis = _experimental_high_order_identity_basis(11)
    axis_data = GaussletBases._experimental_high_order_axis_data_1d(
        basis;
        backend = :numerical_reference,
    )
    stack = GaussletBases._experimental_high_order_doside_stack_3d(
        basis;
        backend = :numerical_reference,
        sides = [5, 7, 9, 11],
    )
    union_coefficients = GaussletBases._experimental_high_order_orthonormalized_full_block_union_coefficients(
        axis_data,
        [5, 7, 9, 11];
        doside = 5,
    )
    stack_data = GaussletBases._experimental_high_order_doside_heplus_data(
        stack;
        expansion = expansion,
        Z = 2.0,
    )
    union_data = GaussletBases._experimental_high_order_doside_heplus_data(
        basis,
        union_coefficients;
        backend = :numerical_reference,
        expansion = expansion,
        Z = 2.0,
    )

    @test norm(union_data.projected_overlap - I, Inf) < 1.0e-8
    @test abs(stack_data.ground_energy - union_data.ground_energy) < 1.0e-8
end

@testset "Experimental high-order doside He singlet action symmetry and Hermiticity" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    basis = _experimental_high_order_identity_basis(7)
    stack = GaussletBases._experimental_high_order_doside_stack_3d(
        basis;
        backend = :numerical_reference,
        sides = [5, 7],
    )
    problem = GaussletBases._experimental_high_order_he_singlet_problem(
        stack;
        expansion = expansion,
        Z = 2.0,
    )
    n = size(problem.projected_hamiltonian, 1)
    seed_a = reshape(collect(1.0:(n * n)), n, n)
    seed_b = reshape(collect((n * n + 1.0):(2n * n)), n, n)
    A = GaussletBases._symmetrize_ida_matrix(seed_a)
    B = GaussletBases._symmetrize_ida_matrix(seed_b)
    A ./= norm(A)
    B ./= norm(B)

    HA = GaussletBases._experimental_high_order_he_singlet_action(problem, A)
    HB = GaussletBases._experimental_high_order_he_singlet_action(problem, B)

    @test norm(HA - transpose(HA), Inf) < 1.0e-10
    @test norm(HB - transpose(HB), Inf) < 1.0e-10
    @test all(isfinite, HA)
    @test all(isfinite, HB)
    @test abs(
        GaussletBases._experimental_high_order_frobenius_dot(A, HB) -
        GaussletBases._experimental_high_order_frobenius_dot(HA, B)
    ) < 1.0e-8
end

@testset "Experimental high-order doside He singlet ladder" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    energies = Float64[]

    for sides in ([5], [5, 7], [5, 7, 9], [5, 7, 9, 11])
        basis = _experimental_high_order_identity_basis(maximum(sides))
        stack = GaussletBases._experimental_high_order_doside_stack_3d(
            basis;
            backend = :numerical_reference,
            sides = sides,
        )
        data = GaussletBases._experimental_high_order_doside_he_singlet_data(
            stack;
            expansion = expansion,
            Z = 2.0,
            krylovdim = 32,
            maxiter = 32,
            tol = 1.0e-8,
        )
        @test isfinite(data.ground_energy)
        @test isfinite(data.residual)
        @test data.residual < 1.0e-5
        @test all(isfinite, data.ground_matrix)
        @test norm(data.ground_matrix - transpose(data.ground_matrix), Inf) < 1.0e-10
        push!(energies, data.ground_energy)
    end

    @test all(diff(energies) .<= 1.0e-5)
end

@testset "Experimental high-order doside He singlet matches orthonormalized full-union reference" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    basis = _experimental_high_order_identity_basis(11)
    axis_data = GaussletBases._experimental_high_order_axis_data_1d(
        basis;
        backend = :numerical_reference,
    )
    stack = GaussletBases._experimental_high_order_doside_stack_3d(
        basis;
        backend = :numerical_reference,
        sides = [5, 7, 9, 11],
    )
    union_coefficients = GaussletBases._experimental_high_order_orthonormalized_full_block_union_coefficients(
        axis_data,
        [5, 7, 9, 11];
        doside = 5,
    )
    stack_data = GaussletBases._experimental_high_order_doside_he_singlet_data(
        stack;
        expansion = expansion,
        Z = 2.0,
        krylovdim = 32,
        maxiter = 32,
        tol = 1.0e-8,
    )
    union_data = GaussletBases._experimental_high_order_doside_he_singlet_data(
        basis,
        union_coefficients;
        backend = :numerical_reference,
        expansion = expansion,
        Z = 2.0,
        krylovdim = 32,
        maxiter = 32,
        tol = 1.0e-8,
    )

    @test abs(stack_data.ground_energy - union_data.ground_energy) < 1.0e-6
end

@testset "Experimental high-order doside He singlet spacing sensitivity and conditioning" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    ladder = ([5], [5, 7], [5, 7, 9], [5, 7, 9, 11])

    for spacing in (0.9, 1.0, 1.1)
        energies = Float64[]
        for sides in ladder
            case = _experimental_high_order_he_singlet_case(
                sides;
                reference_spacing = spacing,
                expansion = expansion,
            )
            stack = case.stack
            data = case.data

            @test isfinite(data.ground_energy)
            @test isfinite(data.residual)
            @test stack.diagnostics.overlap_error < 1.0e-8
            @test stack.diagnostics.overlap_spectrum.minimum_eigenvalue > 1.0 - 1.0e-8
            @test stack.diagnostics.overlap_spectrum.maximum_eigenvalue < 1.0 + 1.0e-8
            @test stack.diagnostics.overlap_spectrum.condition_number < 1.0 + 1.0e-8
            @test data.diagnostics.projected_overlap_spectrum.minimum_eigenvalue > 1.0 - 1.0e-8
            @test data.diagnostics.projected_overlap_spectrum.maximum_eigenvalue < 1.0 + 1.0e-8
            @test data.diagnostics.projected_overlap_spectrum.condition_number < 1.0 + 1.0e-8
            @test all(spectrum -> spectrum.kept_rank == 98, stack.diagnostics.shell_cleanup_spectra)
            @test all(spectrum -> spectrum.minimum_eigenvalue > 1.0e-3, stack.diagnostics.shell_cleanup_spectra)
            push!(energies, data.ground_energy)
        end

        @test all(diff(energies) .<= 1.0e-5)
    end
end

@testset "Experimental high-order doside He singlet parent-box sensitivity" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    energies = Float64[]

    for parent_side in (11, 13)
        case = _experimental_high_order_he_singlet_case(
            [5, 7, 9, 11];
            parent_side = parent_side,
            expansion = expansion,
        )
        stack = case.stack
        data = case.data

        @test stack.parent_side == parent_side
        @test stack.parent_side == length(case.basis)
        @test stack.diagnostics.parent_padding == parent_side - 11
        @test stack.diagnostics.overlap_error < 1.0e-8
        @test stack.diagnostics.overlap_spectrum.minimum_eigenvalue > 1.0 - 1.0e-8
        @test stack.diagnostics.overlap_spectrum.maximum_eigenvalue < 1.0 + 1.0e-8
        @test isfinite(data.ground_energy)
        push!(energies, data.ground_energy)
    end

    @test abs(energies[2] - energies[1]) < 1.0e-8
end

@testset "Experimental high-order doside He singlet shell participation diagnostics" begin
    case = _experimental_high_order_he_singlet_case([5, 7, 9, 11])
    stack = case.stack
    data = case.data
    participation = data.diagnostics.shell_participation

    @test participation.block_labels == stack.block_labels
    @test length(participation.occupation_fractions) == length(stack.block_column_ranges)
    @test length(participation.peak_orbital_fractions) == length(stack.block_column_ranges)
    @test all(value -> isfinite(value) && value >= -1.0e-12, participation.occupation_fractions)
    @test all(value -> isfinite(value) && value >= -1.0e-12, participation.peak_orbital_fractions)
    @test abs(sum(participation.occupation_fractions) - 1.0) < 1.0e-10
    @test participation.total_core_fraction > 0.95
    @test participation.total_outer_shell_fraction > 0.0
    @test participation.total_outer_shell_fraction < 0.05
    @test participation.largest_outer_shell_fraction < 0.05
    @test participation.largest_outer_orbital_fraction < 0.01
end

@testset "Experimental high-order doside He singlet moment-risk audit" begin
    case = _experimental_high_order_he_singlet_case([5, 7, 9, 11])
    stack = case.stack
    data = case.data
    moment_risk = stack.diagnostics.moment_risk
    audit = data.diagnostics.moment_risk_audit
    participation = data.diagnostics.shell_participation

    @test length(moment_risk.column_diagnostics) == size(stack.coefficient_matrix, 2)
    @test moment_risk.outer_column_count == 3 * 98
    @test moment_risk.worst_outer !== nothing
    @test length(moment_risk.top_outer) == 10
    @test all(diagnostic -> isfinite(diagnostic.center_drift), moment_risk.column_diagnostics)
    @test all(diagnostic -> isfinite(diagnostic.max_normalized_mu3), moment_risk.column_diagnostics)
    @test all(diagnostic -> isfinite(diagnostic.max_normalized_mu4), moment_risk.column_diagnostics)
    @test moment_risk.worst_outer.risk_score == max(
        moment_risk.worst_outer.max_normalized_mu3,
        moment_risk.worst_outer.max_normalized_mu4,
    )

    @test audit !== nothing
    @test audit.outer_column_count == moment_risk.outer_column_count
    @test audit.worst_outer_direction !== nothing
    @test audit.worst_outer_direction.block_label != :side5_full
    @test isfinite(audit.worst_outer_direction.risk_score)
    @test isfinite(audit.worst_outer_direction.state_weight)
    @test length(audit.top_k_weights) == 3
    @test audit.top_k_weights[1].k == 1
    @test audit.top_k_weights[2].k == 5
    @test audit.top_k_weights[3].k == 10
    @test audit.top_k_weights[1].total_state_weight < 1.0e-4
    @test audit.top_k_weights[2].total_state_weight < 5.0e-3
    @test audit.top_k_weights[3].total_state_weight < 0.75 * participation.total_outer_shell_fraction
    @test length(audit.shell_block_weights) == length(stack.block_labels)
    @test audit.shell_block_weights[1].block_label == :side5_full
    @test abs(sum(weight.total_state_weight for weight in audit.shell_block_weights) - 1.0) < 1.0e-10
end

@testset "Experimental high-order doside distorted-parent He+ benchmark surface" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    basis = _experimental_high_order_distorted_he_basis(15)
    benchmark = GaussletBases._experimental_high_order_distorted_parent_heplus_benchmark(
        basis;
        backend = :numerical_reference,
        expansion = expansion,
        Z = 2.0,
    )
    high_outer = benchmark.high_order_rows[end]

    @test benchmark.mapping_family == :white_lindsey_atomic_he_d0p2
    @test benchmark.parent_side == 15
    @test benchmark.parent_reference_dimension == 15^3
    @test isfinite(benchmark.parent_reference_energy)
    @test benchmark.lower_order_comparator == :one_center_atomic_legacy_profile_fixed_block
    @test benchmark.comparator_nside == 5
    @test isfinite(high_outer.energy)
    @test high_outer.overlap_error < 1.0e-8
    @test high_outer.support_overlap_error < 1.0e-8
    @test high_outer.support_overlap_minimum_eigenvalue > 1.0 - 1.0e-8
    @test high_outer.support_overlap_maximum_eigenvalue < 1.0 + 1.0e-8
    @test length(benchmark.lower_order_transfer_rows) == 4
    @test length(benchmark.high_order_transfer_rows) == 4
    @test all(row -> row.admitted, benchmark.lower_order_transfer_rows)
    @test all(row -> row.admitted, benchmark.high_order_transfer_rows)
    @test all(row -> row.used_for_reuse, benchmark.lower_order_transfer_rows)
    @test all(row -> row.used_for_reuse, benchmark.high_order_transfer_rows)
    @test length(benchmark.lower_order_rows) == 4
    @test length(benchmark.high_order_rows) == 4
    @test length(benchmark.comparison_rows) == 4
    @test [row.function_count for row in benchmark.lower_order_rows] == [125, 223, 321, 419]
    @test [row.function_count for row in benchmark.high_order_rows] == [125, 223, 321, 419]
    @test [row.outer_side for row in benchmark.comparison_rows] == [5, 7, 9, 11]
    @test [row.function_count for row in benchmark.comparison_rows] == [125, 223, 321, 419]
    @test all(row -> isfinite(row.energy), benchmark.lower_order_rows)
    @test all(row -> isfinite(row.error), benchmark.lower_order_rows)
    @test all(row -> isfinite(row.energy), benchmark.high_order_rows)
    @test all(row -> isfinite(row.error), benchmark.high_order_rows)
    @test all(row -> row.overlap_error < 1.0e-8, benchmark.lower_order_rows)
    @test all(row -> row.overlap_error < 1.0e-8, benchmark.high_order_rows)
    @test all(row -> row.support_overlap_error < 1.0e-8, benchmark.lower_order_rows)
    @test all(row -> row.support_overlap_error < 1.0e-8, benchmark.high_order_rows)
    @test all(row -> isfinite(row.lower_order_error), benchmark.comparison_rows)
    @test all(row -> isfinite(row.high_order_error), benchmark.comparison_rows)
end

@testset "Experimental high-order doside He singlet PGDG smoke" begin
    case = _experimental_high_order_he_singlet_case(
        [5, 7];
        backend = :pgdg_localized_experimental,
        krylovdim = 16,
        maxiter = 16,
        tol = 1.0e-6,
    )
    data = case.data

    @test isfinite(data.ground_energy)
    @test isfinite(data.residual)
    @test all(isfinite, data.ground_matrix)
    @test norm(data.ground_matrix - transpose(data.ground_matrix), Inf) < 1.0e-8
end
