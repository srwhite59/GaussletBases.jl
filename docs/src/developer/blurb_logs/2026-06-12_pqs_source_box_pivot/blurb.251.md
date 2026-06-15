# Pass 251 blurb - clarify independent H2 PQS driver input taxonomy

Role: repo-doer.

Read before starting:

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/pqs_manager_running_log.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.250.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.250.md`

## Current state

Independent H2 PQS now has source plan, final basis, H1, and H1-J diagnostic
materialization under focused driver overrides:

```text
source_plan_status = :available_pqs_diatomic_physical_gausslet_core_shell_source_plan
final_basis_status = :available_pqs_physical_gausslet_final_basis
final_dimension = 471
h1_status = :materialized_pqs_physical_gausslet_h1_solve
h1_j_status = :materialized_pqs_physical_gausslet_h1_j_payload
h1_j_self_coulomb = 0.4569117646737236
fake_pqs_enabled = false
source_backed_fixed_source_oracle_used = false
```

The checked-in independent input is still a target/readiness input:

- `test/driver_inputs/h2_pqs_q5_independent_source_box_r4.jl`

It defaults all physics flags off:

```julia
run_final_basis = false
run_h1 = false
run_h1_j = false
run_private_rhf = false
```

That was correct during readiness work, but it is now easy for future agents to
misread. The endpoint manifest is also stale: it still has a generic planned
independent-H2 row, and it references older H2 diagnostic input/test names that
are not present in `test/driver_inputs` or `test/nested`.

## Task

Do a taxonomy/input cleanup pass only. Do not change route/source behavior.

Required:

1. Keep the existing independent input as the no-physics readiness input.
   You may rename it only if you update all references cleanly; otherwise leave
   the filename stable and clarify its role in the manifest.

2. Add explicit small driver input variants for the already-demonstrated
   independent H2 PQS stages:

   ```text
   test/driver_inputs/h2_pqs_q5_independent_source_box_r4_final_basis.jl
   test/driver_inputs/h2_pqs_q5_independent_source_box_r4_h1.jl
   test/driver_inputs/h2_pqs_q5_independent_source_box_r4_h1_j.jl
   ```

   Prefer tiny include/override files if that works with Julia include semantics.
   For example, a variant should reuse the base independent input, then only set
   the relevant flags:

   ```julia
   include("h2_pqs_q5_independent_source_box_r4.jl")
   run_final_basis = true
   run_h1 = false
   run_h1_j = false
   run_private_rhf = false
   ```

   H1 and H1-J variants should similarly set only the needed flags. Do not
   request private RHF in any of these files.

3. Update `docs/src/developer/cartesian_driver_endpoint_manifest.md`.

   The manifest should now distinguish:

   - independent H2 PQS target/readiness, endpoint not ready;
   - independent H2 PQS final-basis diagnostic, dimension 471, not solver-ready;
   - independent H2 PQS H1 diagnostic, dimension 471, not solver-ready;
   - independent H2 PQS H1-J density diagnostic, dimension 471, not solver-ready;
   - fake-PQS source-backed WL/QW reproduction remains fake and separate.

   Remove or correct stale rows that point to nonexistent H2 source-box
   diagnostic input/test files, unless those files actually exist and are still
   live.

4. If there is a lightweight test or manifest assertion that references the old
   planned row, update it. Do not add a new slow endpoint test in this pass.

## Strict exclusions

Do not modify:

- H2 route/source-plan/final-basis/H1/H1-J implementation;
- fake-PQS implementation;
- RHF/private RHF behavior;
- supplements, CR2, exports, or public API readiness.

Do not run slow H2 driver validations just to prove the variant files. Passes
248-250 already validated the slow route seams.

## Line budget

Keep scoped `src + test + bin` net-negative if you touch any source/test/bin
file. For docs-only plus driver-input-file additions, report both:

```text
scoped src+test+bin
docs+test/driver_inputs manifest/input taxonomy
```

Because adding explicit tiny input variants may be a small positive outside the
usual `src + test + bin` scope, offset it where practical by deleting stale
manifest rows and stale staged-test mirror assertions. Do not weaken active
route checks.

If this cannot be kept modest and clearly simplifying, write
`.agent_handoffs/ATTENTION.md` and stop.

## Validation

Run only lightweight validation:

- `git diff --check`;
- a parse/include smoke for the new driver input variants if cheap;
- package load only if you touched Julia source or if the include smoke needs it.

Do not run:

- full H2 driver route checks;
- H1/H1-J/RHF/supplement/CR2 tests;
- full nested suite;
- broad staged low-order policy gates;
- `test/nested/cartesian_pair_stage_low_order_policy_runtests.jl`.

## Report

In `response.251.md`, report:

- new/renamed input files;
- manifest rows changed;
- confirmation fake-PQS remains separate;
- confirmation no source/route behavior changed;
- validation commands;
- line budget, including any input/docs positive cost and any deletion offset.

Do not commit. Leave the worktree ready for manager review.

-- repo-manager@macmini
