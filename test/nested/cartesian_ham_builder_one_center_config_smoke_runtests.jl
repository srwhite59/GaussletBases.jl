using Test
using JLD2

const _CARTESIAN_HAM_BUILDER_DRIVER =
    normpath(joinpath(@__DIR__, "..", "..", "bin", "cartesian_ham_builder.jl"))

function _cartesian_ham_builder_with_args(f, args::Vector{String})
    saved_args = copy(ARGS)
    empty!(ARGS)
    append!(ARGS, args)
    try
        return f()
    finally
        empty!(ARGS)
        append!(ARGS, saved_args)
    end
end

function _write_one_center_be_driver_config(
    configfile;
    reportfile,
    tsvfile,
    basisfile,
    hamfile,
)
    write(
        configfile,
        """
route_family = :white_lindsey_low_order
route_kind = :one_center_be_white_lindsey_driver_config_smoke
atom_symbols = ("Be",)
nuclear_charges = (4,)
atom_locations = ((0.0, 0.0, 0.0),)
radius = 15.0
parent_axis_counts = (x = 7, y = 7, z = 7)
map_backend = :pgdg_localized_experimental

q = 5
n_s = 5
reference_spacing = 1.0
tail_spacing = 10.0
q_to_core_spacing_rule = :standard_pqs_ns_equals_q
core_spacing = 0.15
parent_axis_probe_backend = :pgdg_localized_experimental
raw_product_box_probe_backend = :pgdg_localized_experimental

materializer_backend = :pgdg_localized_experimental
materializer_nside = 5
materialize_route = true
probe_route_configured_one_center_materializer = false
save_artifact = true
save_tsv = true
save_basis_artifact = true
save_ham_artifact = true
outfile = $(repr(reportfile))
tsvfile = $(repr(tsvfile))
basisfile = $(repr(basisfile))
hamfile = $(repr(hamfile))
white_lindsey_Z = 4.0
white_lindsey_expansion = GaussletBases.coulomb_gaussian_expansion(doacc = false)
""",
    )
    return configfile
end

function _jld2_top_keys(file)
    return Set(
        key isa AbstractVector ? join(string.(key), "/") : string(key) for key in keys(file)
    )
end

@testset "cartesian_ham_builder one-center PGDG config smoke" begin
    mktempdir() do dir
        configfile = joinpath(dir, "one_center_be_driver_config.jl")
        reportfile = joinpath(dir, "one_center_be_report.jld2")
        tsvfile = joinpath(dir, "one_center_be_report.tsv")
        basisfile = joinpath(dir, "one_center_be_basis.jld2")
        hamfile = joinpath(dir, "one_center_be_ham.jld2")
        stdoutfile = joinpath(dir, "one_center_be_driver_stdout.txt")
        _write_one_center_be_driver_config(
            configfile;
            reportfile,
            tsvfile,
            basisfile,
            hamfile,
        )

        open(stdoutfile, "w") do io
            redirect_stdout(io) do
                _cartesian_ham_builder_with_args([configfile]) do
                    include(_CARTESIAN_HAM_BUILDER_DRIVER)
                end
            end
        end

        @test isfile(reportfile)
        @test isfile(tsvfile)
        @test isfile(basisfile)
        @test isfile(hamfile)

        jldopen(reportfile, "r") do file
            top_keys = _jld2_top_keys(file)
            @test "report" in top_keys
            @test "materialization" in top_keys
            materialization = file["materialization"]
            @test materialization.status ==
                  :materialized_route_configured_one_center_report_available
            @test materialization.route_configured_shellization_consumed
            @test materialization.shellization_source ==
                  :route_configured_one_center_low_order
            @test materialization.route_configured_materializer_backend_requested ==
                  :pgdg_localized_experimental
            @test materialization.route_configured_materializer_backend_consumed ==
                  :pgdg_localized_experimental
            @test materialization.route_configured_materializer_d_requested == 0.15
            @test materialization.route_configured_materializer_d_consumed == 0.15
            @test materialization.route_configured_materializer_nside_requested == 5
            @test materialization.route_configured_materializer_nside_consumed == 5
            @test materialization.retained_dimension == 223
            @test materialization.ham_artifact_written
            @test materialization.basis_artifact_written
        end

        jldopen(basisfile, "r") do file
            top_keys = _jld2_top_keys(file)
            @test "basis" in top_keys
            @test "meta" in top_keys
            @test !("ham" in top_keys)
            @test String(file["meta/route_configured_materializer_backend_requested"]) ==
                  "pgdg_localized_experimental"
            @test String(file["meta/route_configured_materializer_backend_consumed"]) ==
                  "pgdg_localized_experimental"
            @test file["meta/route_configured_materializer_d_requested"] == 0.15
            @test file["meta/route_configured_materializer_d_consumed"] == 0.15
            @test file["meta/route_configured_materializer_nside_requested"] == 5
            @test file["meta/route_configured_materializer_nside_consumed"] == 5
            @test Bool(file["meta/route_configured_shellization_consumed"])
            @test String(file["meta/shellization_source"]) ==
                  "route_configured_one_center_low_order"
        end

        jldopen(hamfile, "r") do file
            top_keys = _jld2_top_keys(file)
            @test "basis" in top_keys
            @test "ham" in top_keys
            @test "meta" in top_keys
            @test size(file["ham/overlap"]) == (223, 223)
            @test size(file["ham/one_body_hamiltonian"]) == (223, 223)
            @test size(file["ham/interaction_matrix"]) == (223, 223)
            @test String(file["meta/route_configured_materializer_backend_requested"]) ==
                  "pgdg_localized_experimental"
            @test String(file["meta/route_configured_materializer_backend_consumed"]) ==
                  "pgdg_localized_experimental"
            @test file["meta/route_configured_materializer_d_requested"] == 0.15
            @test file["meta/route_configured_materializer_d_consumed"] == 0.15
            @test file["meta/route_configured_materializer_nside_requested"] == 5
            @test file["meta/route_configured_materializer_nside_consumed"] == 5
            @test Bool(file["meta/route_configured_shellization_consumed"])
            @test String(file["meta/shellization_source"]) ==
                  "route_configured_one_center_low_order"
        end

        tsv = read(tsvfile, String)
        @test occursin(
            "route_materialization\troute_configured_materializer_backend_requested\t:pgdg_localized_experimental",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_materializer_backend_consumed\t:pgdg_localized_experimental",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_materializer_d_requested\t0.15",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_materializer_d_consumed\t0.15",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_materializer_nside_requested\t5",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_materializer_nside_consumed\t5",
            tsv,
        )
        @test occursin(
            "route_materialization\troute_configured_shellization_consumed\ttrue",
            tsv,
        )
    end
end
