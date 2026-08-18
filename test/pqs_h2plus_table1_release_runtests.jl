using GaussletBases
using Test

@testset "PQS Table I matched H2+ comparison" begin
    reference = -0.6026342144949465
    comparison = pqs_h2plus_comparison(
        independent_reference_total_energy_hartree = reference)
    rows = comparison.rows

    @test comparison isa PQSH2PlusComparison
    @test comparison.independent_reference_total_energy_hartree == reference
    @test comparison.parent_axis_counts == (21, 21, 29)
    @test comparison.parent_dimension == 12789
    @test comparison.direct_core_columns == 275
    @test comparison.complete_shell_count == 8
    @test comparison.complete_shell_columns == 960
    @test comparison.slab_count == 2
    @test comparison.slab_columns == 50
    @test getfield.(rows, :basis) == (:parent, :pqs, :white_lindsey)
    @test getfield.(rows, :dimension) == (12789, 1285, 1285)

    expected_capture = (1.0, 0.9999709356338458, 0.999956210164623)
    expected_electronic = (
        -1.1022405729510154,
        -1.1019722712680333,
        -1.1018754077372366,
    )
    expected_total = (
        -0.6022405729510154,
        -0.6019722712680333,
        -0.6018754077372366,
    )
    expected_contraction = (0.0, 0.0002683016829820861, 0.0003651652137788)
    expected_error = (
        0.000393641543931067,
        0.0006619432269131531,
        0.0007588067577099,
    )
    @test collect(getfield.(rows, :capture)) ≈
        collect(expected_capture) atol = 2.0e-8 rtol = 0
    @test collect(getfield.(rows, :electronic_energy_hartree)) ≈
        collect(expected_electronic) atol = 1.0e-7 rtol = 0
    @test getfield.(rows, :nuclear_repulsion_energy_hartree) == (0.5, 0.5, 0.5)
    @test collect(getfield.(rows, :total_energy_hartree)) ≈
        collect(expected_total) atol = 1.0e-7 rtol = 0
    @test collect(getfield.(rows, :contraction_error_hartree)) ≈
        collect(expected_contraction) atol = 1.0e-7 rtol = 0
    @test collect(getfield.(rows, :total_error_hartree)) ≈
        collect(expected_error) atol = 1.0e-7 rtol = 0
    @test_throws ArgumentError pqs_h2plus_comparison(
        independent_reference_total_energy_hartree = NaN)
end
