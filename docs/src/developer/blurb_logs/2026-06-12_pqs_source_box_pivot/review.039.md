Review 039: accepted.

The pass promoted the direct-retained PQS final-basis H1 seam into a durable
workflow gate and shrank old helper-vocabulary integration coverage around
`_pqs_current_route_safe_term_matrices(...)` and
`_pqs_current_route_safe_term_authority_comparison(...)`.

Manager validation:

```text
julia --project=. -e 't = @elapsed include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); println("elapsed_s=", t)'
julia --project=. -e 'Meta.parseall(read("test/nested/bond_aligned_diatomic_high_order_recipe_opt_in_source_construction_integration_runtests.jl", String)); println("integration parse ok")'
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
```

The new H1 gate passed with 29 checks in about 5 seconds. The edited slow
integration file was parse-checked rather than run through the full broad
harness.

Caveat: the new H1 gate is 375 lines, mostly fixture setup. It is accepted
because it is a real workflow/physics-style gate and it replaces old
helper-vocabulary pressure, but do not use it as precedent for adding another
large test. The next PQS work should use `tmp/work` probes unless a permanent
gate clearly replaces more old coverage.

Next step: run a physically interpretable one-center PQS H1 diagnostic through
the direct-retained final-basis path. Do not start IDA or RHF until that one-body
physics target is understood.

-- repo-manager@macmini
