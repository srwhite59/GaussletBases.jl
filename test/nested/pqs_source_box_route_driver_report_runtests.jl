using Test
using GaussletBases

const _ROUTE_DRIVER_STANDARD_UNIT_INVENTORY_KEYS = (
    :unit_count,
    :unit_keys,
    :retained_unit_kinds,
    :source_families,
    :source_dimensions,
    :retained_counts,
    :retained_dimension,
    :retained_counts_materialized,
    :retained_ranges_materialized,
    :pair_count,
    :pair_family_counts,
    :pair_families,
    :output_representations,
)
const _PQS_ROUTE_DRIVER_TEST_ROUTE_KIND = :be2_cartesian_nesting_route_driver_spine

function _pqs_route_driver_report_for_test(; route_family = :pqs_source_box)
    return GaussletBases._pqs_source_box_route_driver_dry_run(
        ;
        route_family,
        route_kind = _PQS_ROUTE_DRIVER_TEST_ROUTE_KIND,
        atom_symbols = ("Be", "Be"),
        nuclear_charges = (4, 4),
        atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
        radius = 15.0,
        parent_axis_counts = (x = 9, y = 7, z = 9),
        map_backend = :pgdg_localized_experimental,
        q = 5,
        n_s = 5,
        reference_spacing = 1.0,
        tail_spacing = 10.0,
        q_to_core_spacing_rule = :standard_pqs_ns_equals_q,
        core_spacing = nothing,
        probe_parent_axis_construction = false,
        parent_axis_probe_backend = :pgdg_localized_experimental,
        parent_axis_probe_family = :G10,
        probe_raw_product_box_plans = false,
        raw_product_box_probe_backend = :pgdg_localized_experimental,
        route_shape = (:pqs_left, :product, :pqs_right),
        product_body_rule = :centered_single_z_slab,
        pqs_retained_rule = :boundary_comx_product_mode_selection,
        product_retained_rule = :product_doside_retained_unit,
        terms = (
            :overlap,
            :position_x,
            :position_y,
            :position_z,
            :x2_x,
            :x2_y,
            :x2_z,
            :kinetic,
        ),
        pair_factor_normalization = :density_normalized,
        support_dense_direct_allowed = false,
        reference_only_authorities = (:support_row_oracle, :dense_parent_projection),
    )
end

function _pqs_route_driver_check_standard_unit_inventory(
    report;
    retained_dimension,
    unit_keys,
    pair_count,
    pair_family_counts,
    retained_counts_materialized,
    retained_ranges_materialized,
    pair_families,
    output_representations,
)
    @test report.retained_dimension == retained_dimension
    @test report.recipe_metadata.route_kind == _PQS_ROUTE_DRIVER_TEST_ROUTE_KIND
    @test Tuple(unit.unit_key for unit in report.retained_units) == unit_keys
    @test length(report.pair_entries) == pair_count
    @test report.pair_family_counts == pair_family_counts

    @test hasproperty(report, :standard_unit_inventory)
    summary = report.standard_unit_inventory
    @test Tuple(keys(summary)) == _ROUTE_DRIVER_STANDARD_UNIT_INVENTORY_KEYS
    @test summary.unit_count == length(unit_keys)
    @test summary.unit_keys == unit_keys
    @test summary.retained_unit_kinds ==
          Tuple(unit.retained_unit_kind for unit in report.retained_units)
    @test summary.source_families ==
          Tuple(unit.source_family for unit in report.retained_units)
    @test summary.source_dimensions == report.source_dimensions
    @test summary.retained_counts == report.retained_counts
    @test summary.retained_dimension == retained_dimension
    @test summary.retained_counts_materialized == retained_counts_materialized
    @test summary.retained_ranges_materialized == retained_ranges_materialized
    @test summary.pair_count == pair_count
    @test summary.pair_family_counts == pair_family_counts
    @test summary.pair_families == pair_families
    @test summary.output_representations == output_representations
    return nothing
end

function _pqs_route_driver_check_report_output_sections(report)
    tmpdir = mktempdir()
    try
        textfile = joinpath(tmpdir, "report.txt")
        open(textfile, "w") do io
            redirect_stdout(io) do
                GaussletBases._pqs_source_box_route_driver_print_details(report)
            end
        end
        text = read(textfile, String)
        @test occursin("[standard_unit_inventory]", text)
        @test occursin("[retained_units]", text)
        @test occursin("[pair_inventory]", text)

        tsvfile = joinpath(tmpdir, "report.tsv")
        stdout_file = joinpath(tmpdir, "save_stdout.txt")
        open(stdout_file, "w") do io
            redirect_stdout(io) do
                GaussletBases._pqs_source_box_route_driver_save(
                    report;
                    save_artifact = false,
                    save_tsv = true,
                    outfile = joinpath(tmpdir, "unused.jld2"),
                    tsvfile,
                )
            end
        end
        tsv = read(tsvfile, String)
        @test occursin("standard_unit_inventory\tunit_count", tsv)
        @test occursin("retained_unit\t", tsv)
        @test occursin("pair_entry\t", tsv)
    finally
        rm(tmpdir; recursive = true, force = true)
    end
    @test !isfile("pqs_source_box_route_driver_report.jld2")
    @test !isfile("pqs_source_box_route_driver_report.tsv")
    return nothing
end

@testset "Route-driver standard unit inventory report" begin
    pqs_report = _pqs_route_driver_report_for_test()
    _pqs_route_driver_check_standard_unit_inventory(
        pqs_report;
        retained_dimension = 221,
        unit_keys = (:pqs_left, :pqs_right, :product),
        pair_count = 6,
        pair_family_counts =
            (pqs_pqs = 3, pqs_product = 2, product_pqs = 0, product_product = 1),
        retained_counts_materialized = true,
        retained_ranges_materialized = true,
        pair_families = (:pqs_pqs, :pqs_product, :product_product),
        output_representations = (:retained_two_index_density_density,),
    )
    _pqs_route_driver_check_report_output_sections(pqs_report)

    white_lindsey_report =
        _pqs_route_driver_report_for_test(route_family = :white_lindsey_low_order)
    _pqs_route_driver_check_standard_unit_inventory(
        white_lindsey_report;
        retained_dimension = nothing,
        unit_keys = (:low_order_units,),
        pair_count = 1,
        pair_family_counts = (
            pqs_pqs = 0,
            pqs_product = 0,
            product_pqs = 0,
            product_product = 0,
            white_lindsey_low_order = 1,
        ),
        retained_counts_materialized = false,
        retained_ranges_materialized = false,
        pair_families = (:white_lindsey_low_order,),
        output_representations = (:low_order_nested_cartesian_basis,),
    )
    _pqs_route_driver_check_report_output_sections(white_lindsey_report)
end
