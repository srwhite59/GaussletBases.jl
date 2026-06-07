# Integration/slow test. Do not include in default nested runner.

@testset "Atomic hybrid He orbital transfer remains stable across same-parent different-final-contraction change" begin
    fixture = _atomic_hybrid_he_same_parent_stress_fixture()

    @test fixture.source_ops isa OrdinaryCartesianOperators3D
    @test fixture.target_ops isa OrdinaryCartesianOperators3D
    @test fixture.source_rep.metadata.parent_kind == :cartesian_plus_supplement_raw
    @test fixture.target_rep.metadata.parent_kind == :cartesian_plus_supplement_raw
    @test fixture.source_working_box == 2:6
    @test fixture.target_working_box == (1:7, 1:7, 1:7)
    @test length(fixture.source_ops.orbital_data) == 134
    @test length(fixture.target_ops.orbital_data) == 232
    @test fixture.transfer.diagnostics.transfer_path == :hybrid_mixed_raw_cross_overlap_transfer
    @test fixture.transfer.diagnostics.transferred_residual_inf < 1.0e-10

    @test fixture.source_observables.metric_norm_error < 1.0e-12
    @test fixture.target_observables.metric_norm_error < 1.0e-12
    @test fixture.aligned_transferred_observables.metric_norm_error < 1.0e-12

    source_self_overlap = cross_overlap(fixture.source_rep, fixture.source_rep)
    target_self_overlap = cross_overlap(fixture.target_rep, fixture.target_rep)
    cross_overlap_source_target = cross_overlap(fixture.source_rep, fixture.target_rep)
    @test source_self_overlap ≈ fixture.source_ops.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test target_self_overlap ≈ fixture.target_ops.overlap atol = 1.0e-10 rtol = 1.0e-10
    @test maximum(abs, cross_overlap_source_target) <= 1.0 + 1.0e-10
    @test maximum(svdvals(cross_overlap_source_target)) <= 1.0 + 1.0e-10

    @test fixture.target_observables.total < fixture.source_observables.total
    @test fixture.aligned_overlap_to_target > 0.999995
    @test abs(
        fixture.aligned_transferred_observables.one_body - fixture.target_observables.one_body,
    ) < 1.0e-4
    @test abs(
        fixture.aligned_transferred_observables.vee - fixture.target_observables.vee,
    ) < 2.0e-4
    @test abs(
        fixture.aligned_transferred_observables.total - fixture.target_observables.total,
    ) < 5.0e-4

    mktempdir() do dir
        source_path = joinpath(dir, "he_source_hybrid.jld2")
        target_path = joinpath(dir, "he_target_hybrid.jld2")

        write_cartesian_basis_bundle_jld2(source_path, fixture.source_ops)
        write_cartesian_basis_bundle_jld2(target_path, fixture.target_ops)

        source_bundle = read_cartesian_basis_bundle(source_path)
        target_bundle = read_cartesian_basis_bundle(target_path)
        @test hasproperty(source_bundle.basis.parent_data, :exact_cartesian_supplement_overlap)
        @test hasproperty(source_bundle.basis.parent_data, :exact_supplement_overlap)
        @test hasproperty(target_bundle.basis.parent_data, :exact_cartesian_supplement_overlap)
        @test hasproperty(target_bundle.basis.parent_data, :exact_supplement_overlap)

        disk_source_self = cross_overlap(source_path, source_path)
        disk_target_self = cross_overlap(target_path, target_path)
        disk_cross = cross_overlap(source_path, target_path)
        @test disk_source_self ≈ fixture.source_ops.overlap atol = 1.0e-10 rtol = 1.0e-10
        @test disk_target_self ≈ fixture.target_ops.overlap atol = 1.0e-10 rtol = 1.0e-10
        @test disk_source_self ≈ source_self_overlap atol = 1.0e-10 rtol = 1.0e-10
        @test disk_target_self ≈ target_self_overlap atol = 1.0e-10 rtol = 1.0e-10
        @test disk_cross ≈ cross_overlap_source_target atol = 1.0e-10 rtol = 1.0e-10
        @test maximum(abs, disk_cross) <= 1.0 + 1.0e-10
        @test maximum(svdvals(disk_cross)) <= 1.0 + 1.0e-10

        disk_transfer = transfer_orbitals(
            fixture.source_observables.orbital,
            source_path,
            target_path,
        )

        @test disk_transfer.coefficients ≈ fixture.transfer.coefficients atol = 1.0e-10 rtol = 1.0e-10
        @test disk_transfer.diagnostics.transfer_path ==
            fixture.transfer.diagnostics.transfer_path
    end
end
