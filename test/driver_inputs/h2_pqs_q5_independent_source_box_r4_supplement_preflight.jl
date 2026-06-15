include("h2_pqs_q5_independent_source_box_r4_h1_j.jl")

supplement_policy = :mwg_residual_gto
comparison_ready = false
comparison_blocker = :independent_pqs_supplement_preflight_only
physics_endpoint_ready = false
physics_endpoint_blocker = :missing_provider_gto_supplement_blocks
run_final_basis = true
run_h1 = true
run_h1_j = true
run_private_rhf = false
