Review 222 - accepted

Pass 222 retired the old H2 221-dimensional source-box diagnostic scaffold and
audited the next GTO/MWG supplement seam.

Accepted deletion:

```text
deleted test/driver_inputs/h2_pqs_q5_source_box_diagnostic_r4.jl
deleted test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl
```

The 221 route implementation itself was not touched. The reference audit found
no live source, `bin`, default-runner, or current physics workflow references;
remaining hits are historical handoff/log text only.

Line budget:

```text
0 added
68 deleted
net -68
```

Validation:

```text
rg -n "h2_pqs_q5_source_box_diagnostic_r4" test bin src docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot
  passed; only historical log hits remain

julia --project=. -e 'using GaussletBases; println("load ok")'
  passed

git diff --check
git diff --cached --check
  passed
```

Committed and pushed:

```text
cfe470db Retire H2 source-box diagnostic scaffold
```

Supplement audit result:

The current usable supplement machinery is around the ordinary WL/QW operator
path, legacy supplement constructors, MWG component kernels, and the
route-global combined-GTO CPBM surfaces. The accepted H2 physical PQS endpoint
already has the common physical support plan, final PQS basis, H1, H1-J/density
interaction, and private RHF diagnostic. It does not yet have provider-level
GTO supplement blocks, combined raw GTO matrices, raw moment matrices, a
residual MWG representation, or a route-owned bridge from the physical PQS
final-basis payload into the combined-GTO supplement surfaces.

The smallest future implementation seam should be a private H2 physical
supplement preflight/payload boundary:

```text
common physical support/intermediate gausslet plan
+ retained_transform_kind = :pqs
+ supplement_policy = :mwg_residual_gto
```

Before adding that seam, clean up stale component-smoke/CR2 sidecar machinery
if it is no longer live. That older path mixes source-box diagnostic vocabulary
with final-residual MWG sidecar concepts and could pull the next supplement work
toward an obsolete route-shadow contract.

-- repo-manager@macmini
