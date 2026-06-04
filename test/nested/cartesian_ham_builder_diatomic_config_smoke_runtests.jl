using Test
using JLD2

const _CARTESIAN_HAM_BUILDER_DIATOMIC_DRIVER =
    normpath(joinpath(@__DIR__, "..", "..", "bin", "cartesian_ham_builder.jl"))

function _cartesian_ham_builder_diatomic_with_args(f, args::Vector{String})
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

function _write_diatomic_be2_driver_config(
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
route_kind = :be2_cartesian_diatomic_driver_config_smoke
atom_symbols = ("Be", "Be")
nuclear_charges = (4, 4)
atom_locations = ((-2.0, 0.0, 0.0), (2.0, 0.0, 0.0))
radius = 15.0
parent_axis_counts = (x = 9, y = 7, z = 9)
map_backend = :pgdg_localized_experimental

q = 5
n_s = 5
reference_spacing = 1.0
tail_spacing = 10.0
q_to_core_spacing_rule = :standard_pqs_ns_equals_q
core_spacing = 0.15
parent_axis_probe_backend = :pgdg_localized_experimental
parent_axis_probe_family = :G10
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
white_lindsey_expansion = GaussletBases.coulomb_gaussian_expansion(doacc = false)
""",
    )
    return configfile
end

function _diatomic_jld2_top_keys(file)
    return Set(
        key isa AbstractVector ? join(string.(key), "/") : string(key) for key in keys(file)
    )
end

@testset "cartesian_ham_builder diatomic PGDG config smoke" begin
    mktempdir() do dir
        configfile = joinpath(dir, "be2_driver_config.jl")
        reportfile = joinpath(dir, "be2_report.jld2")
        tsvfile = joinpath(dir, "be2_report.tsv")
        basisfile = joinpath(dir, "be2_basis.jld2")
        hamfile = joinpath(dir, "be2_ham.jld2")
        stdoutfile = joinpath(dir, "be2_driver_stdout.txt")
        _write_diatomic_be2_driver_config(
            configfile;
            reportfile,
            tsvfile,
            basisfile,
            hamfile,
        )

        open(stdoutfile, "w") do io
            redirect_stdout(io) do
                _cartesian_ham_builder_diatomic_with_args([configfile]) do
                    include(_CARTESIAN_HAM_BUILDER_DIATOMIC_DRIVER)
                end
            end
        end

        @test isfile(reportfile)
        @test isfile(tsvfile)
        @test isfile(basisfile)
        @test isfile(hamfile)

        jldopen(reportfile, "r") do file
            top_keys = _diatomic_jld2_top_keys(file)
            @test "report" in top_keys
            @test "materialization" in top_keys
            materialization = file["materialization"]
            @test materialization.status ==
                  :materialized_route_configured_diatomic_shellization_available
            @test materialization.retained_dimension == 1271
            @test materialization.route_configured_shellization_consumed
            @test materialization.shellization_source ==
                  :route_configured_bond_aligned_diatomic_source
            @test !materialization.route_configured_diatomic_seed_fallback
            @test materialization.route_configured_materializer_backend_requested ==
                  :pgdg_localized_experimental
            @test materialization.route_configured_materializer_backend_consumed ==
                  :pgdg_localized_experimental
            @test materialization.route_configured_materializer_d_requested == 0.15
            @test materialization.route_configured_materializer_d_consumed == 0.15
            @test materialization.route_configured_materializer_nside_requested == 5
            @test materialization.route_configured_materializer_nside_consumed == 5
            @test materialization.basis_artifact_written
            @test materialization.ham_artifact_written
            @test materialization.ham_bundle_export_status ==
                  :available_route_configured_diatomic_ham_bundle_payload
            @test materialization.ham_artifact_status ==
                  :written_route_configured_diatomic_ham_bundle
        end

        jldopen(basisfile, "r") do file
            top_keys = _diatomic_jld2_top_keys(file)
            @test "basis" in top_keys
            @test "meta" in top_keys
            @test !("ham" in top_keys)
            @test String(file["basis/format"]) == "cartesian_basis_bundle_v1"
            @test String(file["basis/basis_kind"]) == "nested_fixed_block"
            @test file["basis/final_dimension"] == 1271
            @test length(file["basis/final_integral_weights"]) == 1271
            @test all(isfinite, file["basis/final_integral_weights"])
            @test minimum(file["basis/final_integral_weights"]) > 0.0
            @test length(file["basis/basis_labels"]) == 1271
            @test size(file["basis/basis_centers"]) == (1271, 3)
            @test Bool(file["meta/has_ham"]) == false
            @test String(file["meta/materialized_report_kind"]) ==
                  "cartesian_shellization_route_bond_aligned_diatomic_materialization"
            @test String(file["meta/shellization_source"]) ==
                  "route_configured_bond_aligned_diatomic_source"
            @test Bool(file["meta/route_configured_shellization_consumed"])
            @test String(file["meta/export_status"]) == "basis_only"
            @test String(file["meta/basis_export_status"]) ==
                  "supported_route_configured_diatomic_basis_only_fixed_block"
            @test String(file["meta/ham_export_status"]) ==
                  "artifact_local_basis_only_no_ham_payload"
            @test Bool(file["meta/ham_export_blocker/is_nothing"])
            @test Bool(file["meta/companion_ham_artifact_requested"])
            @test String(file["meta/companion_ham_artifact_status"]) ==
                  "companion_route_configured_diatomic_ham_artifact_ready"
            @test String(file["meta/companion_ham_export_status"]) ==
                  "available_route_configured_diatomic_ham_bundle_payload"
            @test Bool(file["meta/companion_ham_export_blocker/is_nothing"])
        end

        jldopen(hamfile, "r") do file
            top_keys = _diatomic_jld2_top_keys(file)
            @test "basis" in top_keys
            @test "ham" in top_keys
            @test "meta" in top_keys
            @test String(file["basis/format"]) == "cartesian_basis_bundle_v1"
            @test String(file["basis/basis_kind"]) == "nested_fixed_block"
            @test String(file["ham/format"]) == "cartesian_hamiltonian_bundle_v1"
            @test String(file["ham/model_kind"]) == "ordinary_cartesian_operators"
            @test String(file["ham/interaction_treatment"]) == "ggt_nearest"
            @test String(file["ham/gausslet_backend"]) ==
                  "pgdg_localized_experimental"
            @test file["basis/final_dimension"] == 1271
            @test size(file["ham/overlap"]) == (1271, 1271)
            @test size(file["ham/one_body_hamiltonian"]) == (1271, 1271)
            @test size(file["ham/interaction_matrix"]) == (1271, 1271)
            @test length(file["ham/basis_integral_weights"]) == 1271
            @test file["ham/basis_integral_weights"] ==
                  file["basis/final_integral_weights"]
            @test length(file["ham/orbital_labels"]) == 1271
            @test size(file["ham/basis_centers"]) == (1271, 3)
            @test file["ham/default_nuclear_charges"] == [4.0, 4.0]
            @test String(file["ham/nuclear_term_storage"]) == "by_center"
            @test file["ham/nuclear_one_body_by_center/count"] == 2
            @test Bool(file["meta/has_ham"])
            @test String(file["meta/export_status"]) == "basis_and_ham"
            @test String(file["meta/shellization_source"]) ==
                  "route_configured_bond_aligned_diatomic_source"
            @test Bool(file["meta/route_configured_shellization_consumed"])
            @test String(file["meta/ham_preflight_status"]) ==
                  "available_route_configured_diatomic_ham_adapter"
            @test String(file["meta/ham_operator_payload_status"]) ==
                  "available_route_configured_diatomic_operator_payload"
            @test String(file["meta/ham_interaction_status"]) ==
                  "available_route_configured_diatomic_density_density_interaction_matrix"
            @test String(file["meta/ham_export_status"]) ==
                  "available_route_configured_diatomic_ham_bundle_payload"
            @test Bool(file["meta/ham_export_blocker/is_nothing"])
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
            "route_materialization\tham_artifact_status\t:written_route_configured_diatomic_ham_bundle",
            tsv,
        )
    end
end
