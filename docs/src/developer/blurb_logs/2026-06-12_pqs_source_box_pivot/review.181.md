Accepted pass 181.

The corrected 419-dimensional fixed-q He basis now has a focused H1-J /
pre-final density-interaction smoke:

```text
final dimension = 419
H1 energy = -1.9866819751748936
H1 orbital reconstruction error = 2.0872192862952943e-14
self-Coulomb = 1.2261626003119184
density gauge = pre_final_localized_positive_weight
density-interaction status = materialized_pqs_complete_core_shell_pre_final_density_interaction
```

The density inputs come from
`CartesianContractedParentMetrics._pqs_source_box_ida_factor_provenance(...)`,
specifically `axis_weights` and `raw_axis_pair_factor_terms`, with interaction
path `:ida_gausslet_source_box`. The test does not use retained diagnostic
weights, density-normalized pair terms as authority, signed-final-weight
division, raw-no-division, or old fixed-block oracle data.

Validation passed:

```text
julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test
```

The source/test line budget is net-negative by 128 lines:

```text
36 added, 164 deleted
```

The deleted Be2 fingerprint was standalone artifact/route scaffold from the
paused and explicitly non-comparable Be2/WL/PQS artifact line. It was not in
the default nested runner and had no live source/test caller.

Next step should be a local, ignored RHF/HF diagnostic probe over this same
419-dimensional H1-J payload to get the first comparable scalar against the old
WL Fig.8-style He result. Do not wire RHF into the route or add a tracked
solver test yet.

-- repo-manager@macmini
