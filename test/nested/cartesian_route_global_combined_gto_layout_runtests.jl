# Runtime role: route-global combined gausslet+GTO layout contract.
#
# This test covers metadata-only layout arithmetic for appending a small GTO
# supplement sector to a decomposed WL retained gausslet sector. It does not
# assemble overlap, Hamiltonian, route-driver data, exports, or artifacts.

using Test
using GaussletBases

const CombinedGTOLayoutCPBM = GaussletBases.CartesianPairBlockMaterialization

function _combined_gto_layout_inventory()
    return (;
        status = :available_white_lindsey_decomposed_unit_pair_inventory,
        blocker = nothing,
        source_kind = :test_decomposed_wl_unit_pair_inventory,
        unit_count = 2,
        pair_count = 3,
        retained_dimension = 5,
        retained_dimension_status =
            :available_from_decomposed_wl_unit_column_ranges,
    )
end

function _combined_gto_layout_supplement()
    return basis_representation(
        legacy_atomic_gaussian_supplement("H", "cc-pVTZ"; lmax = 0),
    )
end

@testset "route-global combined gausslet+GTO basis layout" begin
    inventory = _combined_gto_layout_inventory()
    supplement = _combined_gto_layout_supplement()
    layout =
        CombinedGTOLayoutCPBM.route_global_combined_gto_basis_layout(
            inventory,
            supplement,
        )

    @test layout.object_kind == :route_global_combined_gto_basis_layout
    @test layout.status == :available_route_global_combined_gto_basis_layout
    @test isnothing(layout.blocker)
    @test layout.layout_kind ==
          :decomposed_wl_gausslet_plus_gto_supplement
    @test layout.gausslet_retained_range == 1:5
    @test layout.gausslet_retained_dimension == 5
    @test layout.gausslet_retained_dimension_source ==
          :decomposed_wl_unit_pair_inventory_retained_dimension
    @test layout.decomposed_unit_count == 2
    @test layout.decomposed_unit_pair_count == 3
    @test layout.gto_supplement_range == 6:8
    @test layout.gto_supplement_orbital_count == 3
    @test layout.gto_supplement_orbital_labels == ("s1", "s2", "s3")
    @test layout.total_combined_dimension == 8
    @test layout.combined_basis_dimension == 8
    @test layout.block_layout_keys ==
          (:gausslet_gausslet, :gausslet_gto, :gto_gausslet, :gto_gto)
    @test layout.gausslet_gausslet_block_layout.row_range == 1:5
    @test layout.gausslet_gausslet_block_layout.column_range == 1:5
    @test layout.gausslet_gto_block_layout.row_range == 1:5
    @test layout.gausslet_gto_block_layout.column_range == 6:8
    @test layout.gto_gausslet_block_layout.row_range == 6:8
    @test layout.gto_gausslet_block_layout.column_range == 1:5
    @test layout.gto_gto_block_layout.row_range == 6:8
    @test layout.gto_gto_block_layout.column_range == 6:8
    @test layout.mixed_cpb_gto_blocks_orientation ==
          :gausslet_rows_by_gto_columns
    @test layout.gto_gto_blocks_orientation == :gto_rows_by_gto_columns
    @test layout.by_center_nuclear_blocks_separated
    @test layout.nuclear_charge_application_stage ==
          :future_combined_hamiltonian_assembly
    @test layout.combined_overlap_layout_available
    @test layout.combined_hamiltonian_layout_available
    @test !layout.combined_overlap_matrix_materialized
    @test !layout.combined_hamiltonian_matrix_materialized
    @test !layout.route_global_combined_matrix_materialized
    @test !layout.gto_route_global_blocks_materialized
    @test !layout.hamiltonian_assembly
    @test !layout.hamiltonian_data_materialized
    @test !layout.full_parent_window_cpb_used
    @test !layout.direct_cartesian_product_assembly_used
    @test !layout.ordinary_cartesian_ida_operators_used
    @test !layout.pqs_transforms_materialized
    @test !layout.exports_or_artifacts
    @test !hasproperty(layout, :overlap_matrix)
    @test !hasproperty(layout, :hamiltonian_matrix)
    @test !hasproperty(layout, :dense_block)
    @test !hasproperty(layout, :mixed_blocks)
    @test !hasproperty(layout, :gto_blocks)
end

@testset "blocked combined GTO layout paths" begin
    supplement = _combined_gto_layout_supplement()
    missing_inventory =
        CombinedGTOLayoutCPBM.route_global_combined_gto_basis_layout(
            (; status = :blocked_white_lindsey_decomposed_unit_pair_inventory),
            supplement,
        )
    @test missing_inventory.status ==
          :blocked_route_global_combined_gto_basis_layout
    @test missing_inventory.blocker ==
          :missing_decomposed_wl_unit_pair_inventory
    @test missing_inventory.block_layout == :unavailable

    missing_dimension =
        CombinedGTOLayoutCPBM.route_global_combined_gto_basis_layout(
            (;
                status = :available_white_lindsey_decomposed_unit_pair_inventory,
                retained_dimension = nothing,
            ),
            supplement,
        )
    @test missing_dimension.status ==
          :blocked_route_global_combined_gto_basis_layout
    @test missing_dimension.blocker == :missing_decomposed_wl_retained_dimension

    empty_supplement =
        (; supplement_kind = :test_empty_supplement, orbitals = ())
    empty_layout =
        CombinedGTOLayoutCPBM.route_global_combined_gto_basis_layout(
            _combined_gto_layout_inventory(),
            empty_supplement,
        )
    @test empty_layout.status ==
          :blocked_route_global_combined_gto_basis_layout
    @test empty_layout.blocker == :empty_gto_supplement_orbitals
end
