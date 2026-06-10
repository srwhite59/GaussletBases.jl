# Runtime role: post-CPB WL/Cartesian gausslet-only H atom acceptance path.
#
# This is a bounded scientific acceptance check for the post-CPB gausslet-only
# White-Lindsey path. It assembles overlap, kinetic, and electron-nuclear
# by-center local blocks through CPB provider machinery, then solves the
# one-electron H atom problem. It does not call ordinary_cartesian_ida_operators,
# add GTO supplements, PQS retained transforms, CPB/GTO bundle consumption,
# route/global refactors, or broad metadata assertions.

using Test

include("cartesian_wl_cpb_acceptance_support.jl")

const _H_EXACT_ENERGY = -0.5
const _H_OLD_DIRECT_SMALL_FITTED_ENERGY = -0.4706400351534759
const _H_OLD_DIRECT_COARSE_DISTORTED_ENERGY = -0.4966106635473884

@testset "post-CPB WL gausslet-only H atom scientific acceptance" begin
    result = _wl_cpb_acceptance_result(
        system = :h_atom,
        centers = ((0.0, 0.0, 0.0),),
        charges = (1.0,),
        axis_counts = (x = 15, y = 15, z = 15),
        old_direct_electronic_baseline = _H_OLD_DIRECT_COARSE_DISTORTED_ENERGY,
        physical_electronic_reference = _H_EXACT_ENERGY,
    )
    report = result.report
    println("post-CPB WL gausslet-only H atom acceptance report: ", report)

    @test report.route == :post_cpb_wl_gausslet_only
    @test report.system == :h_atom
    @test report.q == 5
    @test report.n_s == report.q
    @test report.core_spacing == 0.15
    @test report.axis_counts == (x = 15, y = 15, z = 15)
    @test report.parent_support_size == 3375
    @test report.basis_dimension == 3375
    @test report.center_count == 1
    @test report.nuclear_charges == (1.0,)
    @test report.cpb_local_operator_blocks_materialized
    @test report.cpb_overlap_block_status === :materialized_cpb_overlap_operator_block
    @test report.cpb_kinetic_block_status === :materialized_cpb_kinetic_operator_block
    @test report.cpb_electron_nuclear_block_statuses ==
          (:materialized_cpb_electron_nuclear_by_center_local_block,)
    @test report.cpb_local_electron_nuclear_by_center_blocks_materialized
    @test !report.direct_cartesian_product_assembly_used
    @test !report.ordinary_cartesian_ida_operators_used
    @test report.nuclear_charge_application_stage ===
          :wl_acceptance_hamiltonian_assembly
    @test report.solve_kind == :generalized
    @test size(result.overlap_matrix) == (3375, 3375)
    @test size(result.hamiltonian) == (3375, 3375)
    @test isfinite(report.lowest_electronic_energy)
    @test report.lowest_electronic_energy > _H_EXACT_ENERGY
    @test report.lowest_electronic_energy < -0.48
    @test report.distance_from_physical_electronic_reference > 0.0
    @test report.distance_from_physical_electronic_reference < 0.03
    @test abs(report.distance_from_old_direct_electronic_baseline) < 0.03
    @test abs(
        report.lowest_electronic_energy - _H_OLD_DIRECT_SMALL_FITTED_ENERGY,
    ) < 0.03
end
