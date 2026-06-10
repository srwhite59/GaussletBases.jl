# Runtime role: post-CPB WL/Cartesian gausslet-only H2+ acceptance path.
#
# This is a bounded scientific acceptance check for the post-CPB gausslet-only
# White-Lindsey path. It assembles overlap, kinetic, and by-center
# electron-nuclear local blocks through CPB provider machinery, then solves the
# one-electron H2+ electronic problem at R = 2.0 bohr. It does not call
# ordinary_cartesian_ida_operators, hand-build direct Cartesian product
# matrices, add GTO supplements, PQS retained transforms, route/global
# refactors, or broad metadata assertions.

using Test

include("cartesian_wl_cpb_acceptance_support.jl")

const _H2PLUS_R2_EXACT_ELECTRONIC_ENERGY = -1.1026342144949465
const _H2PLUS_R2_EXACT_TOTAL_ENERGY = -0.6026342144949465
const _H2PLUS_OLD_DIRECT_ELECTRONIC_BASELINE = -1.0654839328172023
const _H2PLUS_OLD_DIRECT_TOTAL_BASELINE = -0.5654839328172023

@testset "post-CPB WL gausslet-only H2+ scientific acceptance" begin
    R = 2.0
    proton_repulsion = 1 / R
    result = _wl_cpb_acceptance_result(
        system = :h2plus,
        centers = ((0.0, 0.0, -R / 2), (0.0, 0.0, R / 2)),
        charges = (1.0, 1.0),
        axis_counts = (x = 15, y = 15, z = 17),
        old_direct_electronic_baseline = _H2PLUS_OLD_DIRECT_ELECTRONIC_BASELINE,
        physical_electronic_reference = _H2PLUS_R2_EXACT_ELECTRONIC_ENERGY,
        physical_total_reference = _H2PLUS_R2_EXACT_TOTAL_ENERGY,
        proton_repulsion = proton_repulsion,
    )
    report = result.report
    println("post-CPB WL gausslet-only H2+ acceptance report: ", report)

    @test report.route == :post_cpb_wl_gausslet_only
    @test report.system == :h2plus
    @test report.q == 5
    @test report.n_s == report.q
    @test report.core_spacing == 0.15
    @test report.axis_counts == (x = 15, y = 15, z = 17)
    @test report.parent_support_size == 3825
    @test report.basis_dimension == 3825
    @test report.center_count == 2
    @test report.centers == ((0.0, 0.0, -1.0), (0.0, 0.0, 1.0))
    @test report.nuclear_charges == (1.0, 1.0)
    @test report.cpb_local_operator_blocks_materialized
    @test report.cpb_overlap_block_status === :materialized_cpb_overlap_operator_block
    @test report.cpb_kinetic_block_status === :materialized_cpb_kinetic_operator_block
    @test report.cpb_electron_nuclear_block_statuses == (
        :materialized_cpb_electron_nuclear_by_center_local_block,
        :materialized_cpb_electron_nuclear_by_center_local_block,
    )
    @test report.cpb_local_electron_nuclear_by_center_blocks_materialized
    @test !report.direct_cartesian_product_assembly_used
    @test !report.ordinary_cartesian_ida_operators_used
    @test report.nuclear_charge_application_stage ===
          :wl_acceptance_hamiltonian_assembly
    @test report.solve_kind == :generalized
    @test size(result.overlap_matrix) == (3825, 3825)
    @test size(result.hamiltonian) == (3825, 3825)
    @test isfinite(report.lowest_electronic_energy)
    @test isfinite(report.born_oppenheimer_total_energy)
    @test report.proton_proton_repulsion == 0.5
    @test report.lowest_electronic_energy > _H2PLUS_R2_EXACT_ELECTRONIC_ENERGY
    @test report.born_oppenheimer_total_energy > _H2PLUS_R2_EXACT_TOTAL_ENERGY
    @test report.lowest_electronic_energy < -1.09
    @test report.born_oppenheimer_total_energy < -0.59
    @test report.distance_from_physical_electronic_reference > 0.0
    @test report.distance_from_physical_electronic_reference < 0.02
    @test report.distance_from_physical_total_reference > 0.0
    @test report.distance_from_physical_total_reference < 0.02
    @test abs(report.distance_from_old_direct_electronic_baseline) < 0.04
    @test abs(
        report.born_oppenheimer_total_energy - _H2PLUS_OLD_DIRECT_TOTAL_BASELINE,
    ) < 0.04
end
