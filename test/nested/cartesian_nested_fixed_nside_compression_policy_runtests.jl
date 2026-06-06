@testset "Cartesian nested fixed-nside compression policy" begin
    (
        basis,
        bundle,
        shell1,
        shell2,
        shell3,
        grow_sequence,
        shrinking_sequence,
        fixed_grow,
        fixed_shrink,
        legacy,
        baseline,
        grow_ops,
        shrink_ops,
        baseline_check,
        grow_check,
        shrink_check,
    ) = _nested_qiu_white_nside_sequence_fixture()

    @test GaussletBases._nested_shrunk_interval(4:14, 0; nside = 5) == 4:14
    @test GaussletBases._nested_shrunk_interval(4:14, 1; nside = 5) == 5:13
    @test GaussletBases._nested_shrunk_interval(4:14, 2; nside = 5) == 6:12
    @test GaussletBases._nested_shrunk_interval(4:14, 3; nside = 5) == 7:11
    @test GaussletBases._nested_shrunk_interval(4:14, 4; nside = 5) == 7:11

    @test shrinking_sequence isa GaussletBases._CartesianNestedShellSequence3D
    @test length(shrinking_sequence.shell_layers) == 3
    @test shrinking_sequence.shell_layers[1] === shell1
    @test shrinking_sequence.shell_layers[2] === shell2
    @test shrinking_sequence.shell_layers[3] === shell3
    @test length(grow_sequence.core_indices) == 5^3
    @test length(shrinking_sequence.core_indices) == 5^3
    @test first(shrinking_sequence.core_column_range) == 1
    @test last(shrinking_sequence.core_column_range) == 3^3
    @test isempty(intersect(shrinking_sequence.core_indices, shell1.support_indices))
    @test isempty(intersect(shrinking_sequence.core_indices, shell2.support_indices))
    @test isempty(intersect(shrinking_sequence.core_indices, shell3.support_indices))
    @test isempty(intersect(shell1.support_indices, shell2.support_indices))
    @test isempty(intersect(shell1.support_indices, shell3.support_indices))
    @test isempty(intersect(shell2.support_indices, shell3.support_indices))

    @test fixed_grow isa GaussletBases._NestedFixedBlock3D
    @test fixed_shrink isa GaussletBases._NestedFixedBlock3D
    @test grow_ops.gausslet_count == 287
    @test shrink_ops.gausslet_count == 189
    @test baseline.gausslet_count == 17^3
    @test norm(fixed_grow.overlap - I, Inf) < 1.0e-10
    @test norm(fixed_shrink.overlap - I, Inf) < 1.0e-10
    @test grow_check.overlap_error < 1.0e-10
    @test shrink_check.overlap_error < 1.0e-10
    @test isfinite(grow_check.orbital_energy)
    @test isfinite(grow_check.vee_expectation)
    @test isfinite(shrink_check.orbital_energy)
    @test isfinite(shrink_check.vee_expectation)
    @test grow_check.orbital_energy < 0.0
    @test shrink_check.orbital_energy < 0.0
    @test grow_check.vee_expectation > 0.0
    @test shrink_check.vee_expectation > 0.0
    @test shrink_ops.gausslet_count < grow_ops.gausslet_count
end
