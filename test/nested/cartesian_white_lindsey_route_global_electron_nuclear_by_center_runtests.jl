using Test
using GaussletBases

include("cartesian_white_lindsey_adapter_fixture_helpers.jl")

const WLRouteGlobalENCPBM = GaussletBases.CartesianPairBlockMaterialization

function _wl_route_global_en_parent_axis_bundle()
    doside_source_1d = _lw_adapter_doside_source_1d()
    return (;
        x = doside_source_1d,
        y = doside_source_1d,
        z = doside_source_1d,
    )
end

function _wl_route_global_en_centers()
    return (
        (;
            center_key = :proton_a,
            center_index = 1,
            nuclear_charge = 1.0,
            location = (0.0, 0.0, 0.0),
        ),
    )
end

@testset "route-global decomposed WL electron-nuclear by-center adapter" begin
    seed_report = GaussletBases._white_lindsey_low_order_materialized_seed_report()
    expansion = coulomb_gaussian_expansion(doacc = false)
    center_records = _wl_route_global_en_centers()

    by_center = WLRouteGlobalENCPBM.route_global_electron_nuclear_by_center_matrices(
        seed_report;
        parent_axis_counts = (7, 7, 7),
        parent_axis_bundle_object = _wl_route_global_en_parent_axis_bundle(),
        coulomb_expansion = expansion,
        center_records,
    )

    @test by_center.status ==
          :materialized_route_global_electron_nuclear_by_center_matrix_set
    @test isnothing(by_center.blocker)
    @test by_center.inventory_source_kind ==
          :white_lindsey_low_order_materialized_seed_ranges
    @test by_center.retained_dimension == 223
    @test by_center.retained_dimension_status ==
          :available_from_decomposed_wl_unit_column_ranges
    @test by_center.unit_count == 27
    @test by_center.pair_count == 378
    @test by_center.center_count == 1
    @test by_center.by_center_matrix_count == 1
    @test by_center.center_indices == (1,)
    @test by_center.center_keys == (:proton_a,)
    @test by_center.by_center
    @test !by_center.centers_summed
    @test by_center.nuclear_charge_recorded
    @test !by_center.nuclear_charge_applied
    @test by_center.route_global_by_center_matrices_materialized
    @test by_center.global_electron_nuclear_by_center_matrices_materialized
    @test by_center.global_one_body_term_matrix_materialized
    @test !by_center.hamiltonian_data_materialized
    @test !by_center.full_parent_window_cpb_used
    @test !by_center.direct_cartesian_product_assembly_used
    @test !by_center.ordinary_cartesian_ida_operators_used
    @test !by_center.route_driver_wiring
    @test !by_center.exports_materialized

    for result in by_center.matrix_results
        @test result.status ==
              :materialized_route_global_electron_nuclear_by_center_matrix
        @test result.retained_dimension == 223
        @test result.unit_count == 27
        @test result.pair_count == 378
        @test result.local_pair_block_count == 378
        @test result.placeable_record_count == 378
        @test result.global_electron_nuclear_by_center_matrix_materialized
        @test result.global_one_body_term_matrix_materialized
        @test result.by_center
        @test !result.centers_summed
        @test result.nuclear_charge_recorded
        @test !result.nuclear_charge_applied
        @test result.charge_application_stage ==
              :acceptance_or_hamiltonian_assembly
        @test size(result.matrix) == (223, 223)
        @test all(isfinite, result.matrix)
        @test !result.full_parent_window_cpb_used
        @test !result.direct_cartesian_product_assembly_used
        @test !result.ordinary_cartesian_ida_operators_used
        @test !result.hamiltonian_data_materialized
        @test !result.route_driver_wiring
    end
end
