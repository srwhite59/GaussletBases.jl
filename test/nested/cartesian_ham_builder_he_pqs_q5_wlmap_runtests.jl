using Test
using JLD2

const _HE_PQS_DRIVER = normpath(joinpath(@__DIR__, "..", "..", "bin", "cartesian_ham_builder.jl"))
const _HE_PQS_INPUT = normpath(joinpath(@__DIR__, "..", "driver_inputs", "he_pqs_q5_wlmap.jl"))

@testset "cartesian_ham_builder He PQS q5 WL-map artifact" begin
    mktempdir() do dir
        outfile = joinpath(dir, "he_pqs_q5_wlmap.jld2")
        tsvfile = joinpath(dir, "he_pqs_q5_wlmap.tsv")
        saved_args = copy(ARGS)
        empty!(ARGS)
        append!(ARGS, [_HE_PQS_INPUT, "outfile=$(repr(outfile))", "tsvfile=$(repr(tsvfile))"])
        try
            include(_HE_PQS_DRIVER)
        finally
            empty!(ARGS)
            append!(ARGS, saved_args)
        end

        @test isfile(outfile)
        jldopen(outfile, "r") do file
            @test file["basis/final_dimension"] == 419
            @test file["basis/shell_layer_count"] == 3
            @test Tuple(file["basis/retained_per_shell"]) == (98, 98, 98)
            @test file["basis/final_overlap_identity_error"] < 1.0e-10
            @test file["physics/h1_lowest"] ≈ -1.991334820314074 atol = 1.0e-10 rtol = 0.0
            @test file["physics/h1_j_self_coulomb"] ≈ 1.2420423900074902 atol = 1.0e-10 rtol = 0.0
            @test file["density_interaction/density_gauge"] ===
                  :pre_final_localized_positive_weight
            @test file["comparison/delta_h1"] ≈ 9.649649361120893e-6 atol = 1.0e-10 rtol = 0.0
            @test file["comparison/delta_h1_j"] ≈ -4.997485057112172e-6 atol = 1.0e-10 rtol = 0.0
        end
    end
end
