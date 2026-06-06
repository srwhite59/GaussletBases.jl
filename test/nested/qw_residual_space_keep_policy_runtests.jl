@testset "QW residual-space keep policy is near-null-only and stabilized" begin
    # Literal residual-overlap spectrum observed on the anchored one-center
    # Ne legacy-profile case:
    # parent side = 29, working box = 2:28, nside = 7, supplement lmax = 1.
    residual_overlap_eigenvalues = Float64[
        6.486197469e-08,
        3.165964397e-06,
        3.165964398e-06,
        3.165964398e-06,
        5.681904400e-06,
        1.681965647e-05,
        3.337404514e-05,
        5.805312472e-05,
        5.805312472e-05,
        5.805312472e-05,
        7.256172691e-05,
        1.406818079e-04,
        1.406818079e-04,
        1.406818079e-04,
        1.927015773e-04,
        1.927015773e-04,
        1.927015773e-04,
        1.995510583e-04,
        4.261857498e-04,
        4.261857498e-04,
        4.261857498e-04,
        5.359433116e-04,
        1.945893481e-03,
        1.945893481e-03,
        1.945893481e-03,
    ]
    gausslet_overlap = Matrix{Float64}(I, 1, 1)
    overlap_ga = zeros(Float64, 1, length(residual_overlap_eigenvalues))
    overlap_aa = Matrix(Diagonal(residual_overlap_eigenvalues))
    near_null_diagnostics = diagnose_qwrg_residual_space(
        gausslet_overlap,
        overlap_ga,
        overlap_aa;
        keep_policy = :near_null_only,
        keep_abs_tol = GaussletBases._qwrg_atomic_residual_keep_tol(),
        accept_tol = GaussletBases._qwrg_atomic_residual_accept_tol(),
    )
    near_null_data = GaussletBases._qwrg_residual_space(
        gausslet_overlap,
        overlap_ga,
        overlap_aa;
        keep_policy = :near_null_only,
        keep_abs_tol = GaussletBases._qwrg_atomic_residual_keep_tol(),
        accept_tol = GaussletBases._qwrg_atomic_residual_accept_tol(),
    )
    legacy_alias_diagnostics = diagnose_qwrg_residual_space(
        gausslet_overlap,
        overlap_ga,
        overlap_aa;
        keep_policy = :legacy_profile,
        keep_abs_tol = GaussletBases._qwrg_atomic_residual_keep_tol(),
        accept_tol = GaussletBases._qwrg_atomic_residual_accept_tol(),
    )

    @test near_null_diagnostics.keep_policy == :near_null_only
    @test near_null_diagnostics.gaussian_count == 25
    @test near_null_diagnostics.supplement_numerical_rank == 25
    @test near_null_diagnostics.residual_numerical_rank == 25
    @test near_null_diagnostics.kept_count == 24
    @test near_null_diagnostics.discarded_count == 1
    @test near_null_diagnostics.keep_tol ≈ 1.0e-7 atol = 0.0 rtol = 0.0
    @test near_null_diagnostics.accept_tol ≈ 1.0e-7 atol = 0.0 rtol = 0.0
    @test near_null_diagnostics.kept_block_stabilization_null_tol ≈ 1.0e-12 atol = 1.0e-15 rtol = 0.0
    @test near_null_diagnostics.kept_block_stabilization_correction_passes >= 1
    @test near_null_diagnostics.kept_block_stabilization_clipped_count == 0
    @test near_null_diagnostics.kept_block_stabilization_dropped_count == 0
    @test near_null_diagnostics.kept_block_pre_stabilization_overlap_error < 1.0e-12
    @test near_null_diagnostics.kept_block_post_stabilization_overlap_error < 1.0e-12
    @test near_null_diagnostics.kept_block_pre_stabilization_symmetry_defect < 1.0e-12
    @test near_null_diagnostics.kept_block_post_stabilization_symmetry_defect < 1.0e-12
    @test near_null_diagnostics.kept_block_pre_stabilization_negative_count == 0
    @test near_null_diagnostics.kept_block_post_stabilization_negative_count == 0
    @test near_null_diagnostics.kept_block_pre_stabilization_near_null_count == 0
    @test near_null_diagnostics.kept_block_post_stabilization_near_null_count == 0
    @test norm(near_null_data.final_overlap - I, Inf) < 1.0e-10
    @test legacy_alias_diagnostics.keep_policy == :near_null_only
    @test legacy_alias_diagnostics.kept_count == near_null_diagnostics.kept_count
    @test legacy_alias_diagnostics.keep_tol == near_null_diagnostics.keep_tol
    @test legacy_alias_diagnostics.kept_block_post_stabilization_overlap_error ==
        near_null_diagnostics.kept_block_post_stabilization_overlap_error

    nsynthetic = 69
    synthetic_raw_overlap = Matrix{Float64}(I, nsynthetic, nsynthetic)
    synthetic_coefficients = Matrix{Float64}(I, nsynthetic, nsynthetic)
    @inbounds for i in 1:nsynthetic, j in 1:nsynthetic
        synthetic_coefficients[i, j] += 8.0e-9 * sin(Float64(i + 2 * j))
    end
    synthetic_stabilization = GaussletBases._qwrg_stabilize_residual_coefficients(
        synthetic_raw_overlap,
        synthetic_coefficients,
    )
    @test synthetic_stabilization.pre_error > 1.0e-8
    @test synthetic_stabilization.post_error < 1.0e-10
    @test synthetic_stabilization.post_symmetry_defect < 1.0e-12
    @test synthetic_stabilization.pre_negative_count == 0
    @test synthetic_stabilization.post_negative_count == 0
    @test synthetic_stabilization.dropped_count == 0
    @test synthetic_stabilization.correction_passes >= 1
end
