include("h2_pqs_q5_independent_source_box_r4_supplement_preflight.jl")

materialize_route = true
save_basis_artifact = false
save_ham_artifact = true
residual_gto_provider_blocks = :one_body_and_density_provider
hamfile = joinpath(
    tempdir(),
    "h2_pqs_q5_independent_source_box_r4_gto_ham.jld2",
)
