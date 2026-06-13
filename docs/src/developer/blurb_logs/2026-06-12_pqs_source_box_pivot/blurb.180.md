Pass 180 - restore a one-center H1 gate on the corrected 419-dimensional fixed-q PQS basis.

Pass 179 fixed the source/final-basis inventory. The route-owned fixed-q one-center PQS source plan now gives:

```text
parent side = 11
core support = 125
three shell layers
raw_source_dims per shell = (5,5,5)
retained per shell = 98
final retained dimension = 419
```

This pass should put H1 back on that corrected basis. The current
`test/nested/pqs_direct_retained_final_h1_runtests.jl` is now an inventory gate
despite its name. That was acceptable for pass 179, but it should not remain
the state of the repo.

Physics target:

```text
one-center He q=5/n_s=5 fixed-q PQS basis
H1 only
no H1-J
no RHF/SCF/DIIS
no DMRG/HFDMRG
no AHGBS residual/supplement layer
no Be2/Cr2 artifact work
no exports/artifacts/public API
```

Implementation scope:

- Prefer editing only:
  - `test/nested/pqs_direct_retained_final_h1_runtests.jl`
  - possibly deleting stale standalone test scaffold:
    `test/nested/pqs_source_box_route_driver_complete_core_shell_ham_payload_runtests.jl`

- You probably should not need source changes. If H1 construction fails on the
  419 basis because a source helper still assumes the old compact 223 route,
  write `.agent_handoffs/ATTENTION.md` with the exact blocker and stop rather
  than broadening source code casually.

Expected H1 construction:

- Use the pass-179 fixed-q fixture.
- Keep the 419 source/final-basis assertions.
- Add the smallest H1 construction using:

  ```julia
  GaussletBases.pqs_multilayer_complete_core_shell_h1_payload(
      fixture.plan;
      final_basis = fixture.final_basis,
      coulomb_expansion = fixture.expansion,
      center_records = ((;
          center_key = :origin,
          center_index = 1,
          location = (0.0, 0.0, 0.0),
          charge = 2.0,
      ),),
      gaussian_factor_terms_by_center =
          fixture.bundle.pgdg_intermediate.gaussian_factor_terms,
      metadata = (; fixture = :pqs_fixed_q_he_h1_gate),
  )
  ```

  Adjust names to match local style; do not duplicate old explicit-box or fixed-block oracle code.

The test should assert only the live H1 contract:

```text
h1 payload materialized
final H1 matrix finite
final H1 matrix symmetric within tolerance
ordinary symmetric solve
lowest H1 energy finite
lowest H1 energy negative and plausibly He-like
H1 final dimension == 419
no generalized-overlap solve
no H1-J/RHF/export/artifact materialization
```

Do not assert a tight reviewed physics energy yet unless a stable value is observed and justified in the response. This is a restored H1 smoke on the corrected basis, not the final Fig.8 comparison.

Line-budget rule:

```text
git diff --numstat -- src test
sum(deleted) > sum(added)
```

Good deletion candidate:

- `test/nested/pqs_source_box_route_driver_complete_core_shell_ham_payload_runtests.jl`

Why it is a candidate:

- it is a standalone old compact `223` complete-core/shell Ham payload scaffold;
- `rg` shows it is not included by `test/nested/runtests.jl`;
- it preserves the previous compact route fixture rather than the current 419 physics-facing atom target.

Before deleting it, verify quickly that it has no default-runner include and no live role beyond old compact scaffold. If that is true, delete it outright. Do not move it to another default lane. If you find a live caller or a distinct scientific contract, report that and find a directly related smaller deletion instead.

Validation:

```text
julia --project=. test/nested/pqs_direct_retained_final_h1_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test
```

If H1 construction makes the focused test much slower, report the test time and the apparent cost center. Do not run a broad suite.

Reporting requirements:

- exact files changed/deleted;
- observed final dimension and H1 lowest energy;
- matrix finiteness/symmetry result;
- whether old compact 223 Ham scaffold was deleted and why;
- source/test line additions, deletions, and net;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Write the result to `.agent_handoffs/response.180.md` and copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.180.md
```

-- repo-manager@macmini
