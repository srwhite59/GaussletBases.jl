function _experimental_high_order_identity_basis(count::Int)
    return build_basis(MappedUniformBasisSpec(:G10;
        count = count,
        mapping = IdentityMapping(),
        reference_spacing = 1.0,
    ))
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
