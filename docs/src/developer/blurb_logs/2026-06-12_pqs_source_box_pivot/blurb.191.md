Pass 191 - audit the driver-owned H2 PQS target and harvest one safe stale-test deletion.

Purpose:

He q=5/n_s=5 is now driver-owned for H1/H1-J and optional private RHF. The
next physics target should be H2, but do not implement H2 yet. First define the
driver-owned H2 target clearly enough that the implementation does not repeat
the earlier Be2 mistake of comparing unlike bases/routes.

This pass also keeps the line-reduction discipline active. It should delete one
clearly stale development scaffold test/probe if, and only if, it has no live
runner include or source caller.

Task type:

Audit plus deletion. No H2 implementation.

H2 context to recover:

The old WL/QW H2 reference the user has cited is approximately:

```text
H2 at R = 4.0
WL/default route final dimension around 481
WL/endcap-panel route final dimension around 461
HF totals around:
  -0.910938264352
  -0.910977315003
```

Do not treat these numbers as authoritative until you find the exact source in
the repo or old references. The audit should recover the exact fixture,
geometry, route, basis dimensions, and values from tracked tests, docs, reports,
or old fixtures.

Read/inspect:

```text
bin/cartesian_ham_builder.jl
test/driver_inputs/he_pqs_q5_wlmap.jl
test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_runtests.jl
test/nested/cartesian_ham_builder_he_pqs_q5_wlmap_rhf_runtests.jl
test/ordinary/runtests.jl
test/diatomic/runtests.jl
test/nested/*diatomic*
test/nested/*endcap*
test/nested/*ham_builder*
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/
```

Use `rg` first. Do not run broad tests.

Audit questions:

1. What is the exact old WL/QW H2 reference?

   Report:

   ```text
   geometry
   nuclear charges
   mapping/backend
   q/n_s or nside
   route kind
   retained/final dimension
   H1 value if available
   J/self-Coulomb if available
   RHF/HF total if available
   source file or report line
   ```

2. What should the first H2 PQS driver target be?

   It should be driver-owned, visible-stage, and comparable to the WL reference
   in parent mapping/shellification policy. Be explicit about whether this is:

   ```text
   default WL/QW H2 route
   endcap-panel WL/QW H2 route
   PQS source-box analog of one of those
   ```

   Do not overinterpret `q = n_s`. State the retained-basis/shellification
   equivalence requirement.

3. What driver input file would be needed?

   Candidate future file:

   ```text
   test/driver_inputs/h2_pqs_q5_wlmap.jl
   ```

   Report the exact variables it would set, but do not create it in this pass.

4. What artifact keys should the H2 driver test read?

   Keep the pattern from He:

   ```text
   basis/*
   physics/*
   density_interaction/*
   private_rhf/* if requested
   comparison/*
   ```

   Do not request final-basis `S` as downstream working data. It can remain a
   diagnostic identity-error scalar only.

5. What stale tests/probes are safe to delete before H2 implementation?

   Find one development scaffold test/probe whose only role is old route
   transition pressure. Requirements for deletion:

   - no include in `test/nested/runtests.jl`;
   - no include in `test/nested/integration_runtests.jl`;
   - no include in top-level or group runners;
   - no source caller;
   - not one of the new explicit He endpoint tests;
   - not a scientific acceptance/reference test.

   A quick manager inventory showed many explicit/manual nested tests, so do
   not delete by filename pattern alone. Verify the specific file.

   If you cannot identify one safe stale test/probe, write
   `.agent_handoffs/ATTENTION.md` explaining the blocker and stop.

Allowed deletion:

- Delete exactly one clearly stale development scaffold test/probe, plus any
  dead include if one exists.

Forbidden deletion:

- Do not delete the He driver endpoint tests.
- Do not delete WL He/H atom acceptance tests.
- Do not delete H2/WL reference tests or fixtures.
- Do not delete a file only because it is explicit/manual; explicit physics
  endpoints are allowed.

Trust boundary:

- No H2 implementation.
- No new driver input files.
- No new source code.
- No new tests.
- No RHF changes.
- No Be2/Cr2, HFDMRG, DMRG, exports, public solver behavior, ECP, or
  Qiu-White correction work.
- Do not request interactive approval or sandbox escalation. If approval is
  genuinely required, write `.agent_handoffs/ATTENTION.md` and stop.

Line-budget rule:

This pass should be net-negative for:

```text
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Expected:

```text
added: 0
deleted: >0
net: negative
```

The audit response and curated log do not count toward source/test/bin budget.

Validation:

```text
rg -n "<deleted_file_basename>|<deleted_helper_or_testset_name>" test src bin
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test bin tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
git status --short --branch
```

Do not run multi-minute endpoint tests in this audit/deletion pass.

Report:

- exact old WL/QW H2 reference source and values found;
- recommended first H2 PQS driver target;
- future H2 driver input shape;
- future H2 artifact/check shape;
- deleted stale test/probe and why deletion is safe;
- line-budget arithmetic;
- validation results;
- deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

Write the result to `.agent_handoffs/response.191.md` and copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.191.md
```

-- repo-manager@macmini
