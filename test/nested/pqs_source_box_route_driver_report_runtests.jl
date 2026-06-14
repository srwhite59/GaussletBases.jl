# Integration/slow test. Do not include in default nested runner.

using Test
using GaussletBases
using JLD2

const _PQS_ROUTE_DRIVER_TEST_ROUTE_KIND = :be2_cartesian_nesting_route_driver_spine

function _pqs_route_driver_report_for_test(;
    route_family = :pqs_source_box,
    route_kind = _PQS_ROUTE_DRIVER_TEST_ROUTE_KIND,
    atom_symbols = ("Be", "Be"),
    nuclear_charges = (4, 4),
    atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0)),
    parent_axis_counts = (x = 9, y = 7, z = 9),
)
    return GaussletBases._pqs_source_box_route_driver_dry_run(
        ;
        route_family,
        route_kind,
        atom_symbols,
        nuclear_charges,
        atom_locations,
        radius = 15.0,
        parent_axis_counts,
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

function _pqs_route_driver_check_live_report_smoke(report)
    @test report.retained_dimension == 221
    @test report.recipe_metadata.route_kind == _PQS_ROUTE_DRIVER_TEST_ROUTE_KIND
    @test Tuple(unit.unit_key for unit in report.retained_units) ==
          (:pqs_left, :pqs_right, :product)
    @test length(report.pair_entries) == 6
    @test report.parent_contract.system_classification == :bond_aligned_diatomic
    return nothing
end

function _pqs_route_driver_check_report_output_sections(report)
    tmpdir = mktempdir()
    try
        jld2file = joinpath(tmpdir, "report.jld2")
        open(joinpath(tmpdir, "jld2_stdout.txt"), "w") do io
            redirect_stdout(io) do
                GaussletBases._pqs_source_box_route_driver_save(
                    report;
                    save_artifact = true,
                    save_tsv = false,
                    outfile = jld2file,
                    tsvfile = joinpath(tmpdir, "unused.tsv"),
                )
            end
        end
        jldopen(jld2file, "r") do file
            top_keys = Set(
                key isa AbstractVector ? join(string.(key), "/") : string(key) for key in keys(file)
            )
            @test "report" in top_keys
            @test !("materialization" in top_keys)
        end

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
        @test occursin("parent_contract\tsystem_classification", tsv)
        @test occursin("standard_unit_inventory\tunit_count", tsv)
    finally
        rm(tmpdir; recursive = true, force = true)
    end
    @test !isfile("cartesian_ham_builder_report.jld2")
    @test !isfile("cartesian_ham_builder_report.tsv")
    return nothing
end

@testset "Route-driver standard unit inventory report" begin
    pqs_report = _pqs_route_driver_report_for_test()
    _pqs_route_driver_check_live_report_smoke(pqs_report)
    _pqs_route_driver_check_report_output_sections(pqs_report)
end
