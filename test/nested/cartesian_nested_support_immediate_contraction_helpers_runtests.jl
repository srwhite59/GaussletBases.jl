# Integration/slow test. Do not include in default nested runner.

@testset "Cartesian nested support immediate contraction helpers" begin
    expansion = coulomb_gaussian_expansion(doacc = false)
    basis = build_basis(MappedUniformBasisSpec(:G10;
        count = 13,
        mapping = AsinhMapping(a = 0.25, s = asinh(10.0 / 0.25) / (6.0 - 1.0), tail_spacing = 10.0),
        reference_spacing = 1.0,
    ))
    bundle = GaussletBases._mapped_ordinary_gausslet_1d_bundle(
        basis;
        exponents = expansion.exponents,
        backend = :numerical_reference,
        refinement_levels = 0,
    )
    pgdg = bundle.pgdg_intermediate
    term_coefficients = Float64[Float64(value) for value in expansion.coefficients]
    interval = 2:(length(basis) - 1)
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
