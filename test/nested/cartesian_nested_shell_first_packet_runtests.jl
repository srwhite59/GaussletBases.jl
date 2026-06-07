# Integration/slow test. Do not include in default nested runner.

@testset "Cartesian nested shell first packet" begin
    function _fixed_a_nested_shell_basis(count::Int; a::Float64 = 0.25, xmax::Float64 = 10.0, tail_spacing::Float64 = 10.0)
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
    basis, s = _fixed_a_nested_shell_basis(13)
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    pgdg = bundle.pgdg_intermediate
    term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
    interval = 2:(length(basis) - 1)
    shell = GaussletBases._nested_xy_shell_pair(
        bundle,
        interval,
        interval;
        retain_x = 4,
        retain_y = 3,
        term_coefficients,
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
    @test !hasproperty(packet, :gaussian_terms)
    @test !hasproperty(packet, :pair_terms)
    @test !hasproperty(packet, :term_storage)
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
