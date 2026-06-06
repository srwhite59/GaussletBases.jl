@testset "Cartesian nested shell sequence fixed-block" begin
    (
        basis,
        bundle,
        shell1,
        shell2,
        shell_plus_core,
        shell_sequence,
        fixed_shell_plus_core,
        fixed_sequence,
        legacy,
        baseline,
        shell_plus_core_ops,
        shell_sequence_ops,
        baseline_check,
        shell_plus_core_check,
        shell_sequence_check,
    ) = _nested_qiu_white_shell_sequence_fixture()

    @test shell_sequence isa GaussletBases._CartesianNestedShellSequence3D
    @test length(shell_sequence.shell_layers) == 2
    @test shell_sequence.shell_layers[1] === shell1
    @test shell_sequence.shell_layers[2] === shell2
    @test first(shell_sequence.core_column_range) == 1
    @test last(shell_sequence.core_column_range) == length(shell_sequence.core_indices)
    @test length(shell_sequence.layer_column_ranges) == 2
    @test length(shell_sequence.core_indices) == 11^3
    @test isempty(intersect(shell_sequence.core_indices, shell1.support_indices))
    @test isempty(intersect(shell_sequence.core_indices, shell2.support_indices))
    @test isempty(intersect(shell1.support_indices, shell2.support_indices))

    @test fixed_sequence isa GaussletBases._NestedFixedBlock3D
    @test fixed_sequence.parent_basis === basis
    @test fixed_sequence.shell === shell_sequence
    @test shell_plus_core_ops.gausslet_count == 1385
    @test shell_sequence_ops.gausslet_count == 1439
    @test baseline.gausslet_count == 17^3
    @test norm(fixed_shell_plus_core.overlap - I, Inf) < 1.0e-10
    @test norm(fixed_sequence.overlap - I, Inf) < 1.0e-10
    @test shell_plus_core_check.overlap_error < 1.0e-10
    @test shell_sequence_check.overlap_error < 1.0e-10
    @test shell_plus_core_check.orbital_energy < 0.0
    @test shell_sequence_check.orbital_energy < 0.0
    @test shell_plus_core_check.vee_expectation > 0.0
    @test shell_sequence_check.vee_expectation > 0.0
    @test abs(shell_plus_core_check.orbital_energy - baseline_check.orbital_energy) < 1.0e-4
    @test abs(shell_plus_core_check.vee_expectation - baseline_check.vee_expectation) < 1.0e-4
    @test abs(shell_sequence_check.orbital_energy - baseline_check.orbital_energy) < 1.0e-4
    @test abs(shell_sequence_check.vee_expectation - baseline_check.vee_expectation) < 1.0e-4
    @test abs(shell_sequence_check.orbital_energy - shell_plus_core_check.orbital_energy) < 1.0e-4
    @test abs(shell_sequence_check.vee_expectation - shell_plus_core_check.vee_expectation) < 1.0e-4
end
