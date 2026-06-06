@testset "One-center atomic ns=9 legacy-profile residual stabilization closes center-extraction failure" begin
    if !_RUN_SLOW_TESTS
        @test true
    else
        data = _one_center_atomic_ns9_legacy_profile_qw_fixture()
        @test data.residual_data.diagnostics.kept_count == 24
        @test data.residual_data.diagnostics.keep_policy == :near_null_only
        @test data.residual_data.diagnostics.keep_tol ≈ 1.0e-7 atol = 0.0 rtol = 0.0
        @test data.residual_data.diagnostics.accept_tol ≈ 1.0e-7 atol = 0.0 rtol = 0.0
        @test data.residual_data.diagnostics.kept_block_pre_stabilization_overlap_error > 0.0
        @test data.residual_data.diagnostics.kept_block_post_stabilization_overlap_error <=
            data.residual_data.diagnostics.kept_block_pre_stabilization_overlap_error
        @test data.residual_data.diagnostics.kept_block_post_stabilization_overlap_error < 1.0e-9
        @test data.residual_data.diagnostics.kept_block_post_stabilization_symmetry_defect < 1.0e-9
        @test data.residual_data.diagnostics.kept_block_stabilization_dropped_count == 0
        @test norm(data.residual_data.final_overlap - I, Inf) < 1.0e-7
        @test data.operators.residual_count == 24
        @test norm(data.operators.overlap - I, Inf) < 1.0e-7
        check = GaussletBases.ordinary_cartesian_1s2_check(data.operators)
        @test isfinite(check.orbital_energy)
        @test check.overlap_error < 1.0e-7
    end
end
