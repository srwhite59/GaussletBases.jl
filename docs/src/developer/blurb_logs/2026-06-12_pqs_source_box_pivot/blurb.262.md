Pass 262 - pay down support-partition exception by shrinking terminal assembly flat mirrors

Context:
- Current HEAD should include
  `4aa03d8e Materialize independent H2 PQS support partition`.
- Pass 261 was accepted as a deliberate implementation exception:
  scoped `src + test + bin` was `+435 / -0`.
- The new support-partition payload is valuable, but the next cleanup-capable
  pass should pay down stale test/assertion pressure.
- A mature target remains in
  `test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl`: the
  terminal-shellification section still locks many flat assembly/summary mirror
  fields and deferred-status aliases. Similar old-flat report/assembly mirror
  assertions have already been deleted in earlier passes.

Task:
Do a tests-only shrink of
`test/nested/cartesian_assembly_stage_low_order_policy_runtests.jl`.

Target deletion surface:
- In the terminal policy block, remove stale flat-field mirror assertions that
  primarily preserve the old terminal-shellification report/assembly key cloud,
  especially duplicated checks across `terminal_assembly` and
  `terminal_summary` for:
  - `terminal_shellification_*_available` mirrors;
  - deferred pair/materialization status aliases;
  - duplicated matrix/materialization false flags;
  - duplicated pair inventory source/count/family helper mirrors;
  - exact object identity mirrors for scaffold/unit inventories where a compact
    status/count check is enough.

Preserve compact live checks:
- The terminal policy is selected and clearly not active source authority.
- The route remains summary/deferred only.
- Region/unit inventory still has the key physical facts:
  - region count;
  - midpoint slab count;
  - distorted product box count;
  - unit inventory available or equivalent compact count check;
  - owned supports/shellification regions are not CPBs.
- The assembled stage still carries the major stage objects:
  `shells`, `units`, `transforms`, `pairs`.

Strict exclusions:
- Do not edit `src/`.
- Do not edit driver inputs, docs other than the response if needed, or running
  log.
- Do not delete the whole test file.
- Do not change default or atom-growth policy coverage except if a tiny adjacent
  stale mirror assertion is obviously redundant.
- Do not run broad stale low-order integration gates as validation.
- Do not touch independent H2 PQS support partition/provider code.

Validation:
- `git diff --check`.
- Run a lightweight syntax/package-load validation.
- If you run the edited test file, time it and report elapsed seconds. If it is
  slow or fails on unrelated stale terminal-shellification exact-field checks,
  do not chase broad failures; report them and keep the deletion scoped.

Line budget:
- Scoped `src + test + bin` must be net-negative.
- Aim for at least `-50` lines if safe, but do not remove the compact live
  terminal-policy checks listed above just to hit a number.

Report:
- Lines deleted from the test file.
- What stale assertion families were removed.
- What compact terminal-policy checks remain.
- Validation command(s).
- Scoped line count for `src + test + bin`.
- Deletion/shrinkage result:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
