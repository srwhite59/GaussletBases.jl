@testset "Atomic residual keep policy rejects relative_case_scale on public QW routes" begin
    if !_legacy_basisfile_available()
        @test true
    else
        source_basis_qw, _legacy_qw, _ordinary_l0, _ordinary_l0_check = _qiu_white_full_nearest_fixture()
        supplement = legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 1)
        err = try
            ordinary_cartesian_qiu_white_operators(
                source_basis_qw,
                supplement;
                expansion = coulomb_gaussian_expansion(doacc = false),
                Z = 2.0,
                interaction_treatment = :ggt_nearest,
                residual_keep_policy = :relative_case_scale,
            )
            nothing
        catch caught
            caught
        end
        @test err isa ArgumentError
        @test occursin(":near_null_only", sprint(showerror, err))
        @test occursin(":legacy_profile", sprint(showerror, err))
    end
end
