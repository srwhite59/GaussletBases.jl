@testset "Cartesian nested shell interface" begin
    function _fixed_a_multi_face_basis(count::Int; a::Float64 = 0.25, xmax::Float64 = 10.0, tail_spacing::Float64 = 10.0)
        endpoint = (count - 1) / 2
        s = asinh(xmax / a) / (endpoint - xmax / tail_spacing)
        basis = build_basis(MappedUniformBasisSpec(:G10;
            count,
            mapping = AsinhMapping(; a, s, tail_spacing),
            reference_spacing = 1.0,
        ))
        return basis, s
    end

    expansion = coulomb_gaussian_expansion(doacc = false)
    basis, s = _fixed_a_multi_face_basis(13)
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    interval = 2:(length(basis) - 1)
    term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
    shell = GaussletBases._nested_rectangular_shell(
        bundle,
        interval,
        interval,
        interval;
        retain_xy = (4, 3),
        retain_xz = (4, 3),
        retain_yz = (4, 3),
        term_coefficients,
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
    @test !hasproperty(packet, :gaussian_terms)
    @test !hasproperty(packet, :pair_terms)
    @test !hasproperty(packet, :term_storage)
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
