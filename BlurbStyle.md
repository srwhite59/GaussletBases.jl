# Blurb Style Notes

These notes describe how repo-manager should write doer blurbs for
`GaussletBases`. They supplement `AGENTS.md`, `JuliaStyle.md`, and
`docs/src/developer/manager_doer_collaboration_contract_2026-05-26.md`.

The goal is to make each blurb executable without requiring the doer to
rediscover the manager's searches, reasoning, or intended code boundary.

## Default Structure

Use the standard doer blurb shape:

1. Purpose
2. Why now
3. Current state
4. Exact task
5. Trust boundary
6. Artifacts to use
7. What to record
8. Decision rule
9. Report back

Do not force every short blurb into numbered headings, but make sure the same
information is present when the task is technical or risky.

## Include What The Manager Already Knows

If the manager has already searched the repo, name the files, functions,
fixtures, commits, artifacts, or docs that matter. Do not expect the doer to
repeat the same search and arrive at the same interpretation.

Named surfaces are starting points, not authority. If a named file, function,
fixture, or artifact is missing or inconsistent with the live repo, the doer
should stop and report the mismatch instead of silently substituting a different
path.

Prefer:

```text
Use the old nested/fixed-block WL route:
  test/runtests.jl:
    _nested_qiu_white_nearest_fixture
    _nested_qiu_white_shell_sequence_fixture
  ordinary_cartesian_qiu_white_operators(fixed_block, ...)
  _qwrg_fixed_block_interaction_matrix(...)

Do not use the flat/full-product helpers as the nested comparator:
  _qwrg_diatomic_overlap_matrix
  _qwrg_diatomic_kinetic_matrix
  _qwrg_diatomic_nuclear_one_body_by_center
  _qwrg_diatomic_interaction_matrix
```

Avoid:

```text
Compare to the old WL code.
```

The short version hides a distinction the manager already knows. The longer
version saves a wasted pass.

## State Exclusions Explicitly

Every blurb near a confusing boundary should include a "do not use" or
"not in this pass" list. This is especially important when two code paths have
similar names or both look plausible.

Examples:

- Do not use full-parent CPB fallback.
- Do not use ordinary Cartesian IDA operators as acceptance paths.
- Do not use `_qwrg_diatomic_*` flat/full-product helpers as old nested WL.
- Do not treat `pgdg_intermediate.pair_factor_terms` as final retained IDA
  authority.
- Do not add public route wiring, exports, artifacts, PQS changes, or GTO paths
  unless the blurb explicitly requests them.

Exclusions are not filler. They protect the contract.

## Mark The Task Type

Say whether the pass is:

- implementation;
- audit/read-only;
- timing/profiling probe;
- docs-only;
- cleanup/deletion;
- review before commit;
- physics acceptance work.

If it is a probe, say that it is not a test and not production code. If it is
docs-only, say that no production code or tests should change. If it is a
cleanup pass, state the deletion or shrinkage target.

## Tie Work To A Live Target

Each blurb should name the physics target, active numerical kernel, or deletion
target that justifies the work.

Good targets:

- decomposed WL He H1/J/RHF;
- H/H2+ final-basis GTO acceptance;
- retained IDA weight boundary;
- shellification-backed WL inventory performance;
- deletion of an obsolete full-window acceptance path.

Weak targets:

- make architecture more complete;
- add more coverage;
- preserve old behavior;
- might be useful later.

If the target is weak, stop and ask whether the pass should happen.

## Use The Manager Running Log For Cartesian/PQS Work

For Cartesian/PQS work, repo-manager should read
`docs/src/developer/pqs_manager_running_log.md` before drafting a blurb.

Each blurb should be consistent with the long-term goals, anti-goals, and
current medium-term goals in that log. If a proposed pass changes the
medium-term goals, the blurb should say so explicitly and the manager should
update the running log after acceptance.

Running-log updates should be compact but substantive for accepted strategic
passes: usually 100-250 words or a short bullet list, not a two-sentence tick.
They should preserve goal advancement, provenance/authority interpretation,
guardrails, and the next blocker. Purely mechanical passes may use a 1-3
sentence tick only when the entry explicitly says there was no strategic change.

Every 5 accepted Cartesian/PQS passes, repo-manager should add a medium-term
goal checkpoint. Every 10-20 accepted passes, or after a major correction,
repo-manager should add a strategic compression entry. Blurbs that initiate a
checkpoint or compression pass should say so explicitly.

Do not ask the doer to rediscover the strategic state from old pass logs when
the running log already records it.

## Respect Carrying Cost

Every line of source, test, documentation, compatibility glue, metadata, and
adapter code has carrying cost. A blurb that asks for new artifacts should say
what live contract, clarity improvement, duplication reduction, performance
gain, or current workflow justifies that cost.

When two viable approaches have similar value, prefer the one that leaves fewer
stale surfaces. For cleanup and retirement tasks, ask what old code or tests
become removable rather than only what new coverage should be added.

## Give Decision Rules

Doers should know when to continue, stop, or report a blocker.

Examples:

```text
If side-13 H1 improves but J remains far from 1.25, stop and report the IDA
convention evidence before RHF interpretation.
```

```text
If the old nested route cannot be matched to the fixture without broad
construction work, stop with the exact missing object instead of falling back to
the flat route.
```

```text
If the optimized kernel differs from the existing decomposed matrix by more
than 1e-12, do not update baselines; report the mismatch and artifact path.
```

## Require Useful Reporting

Tell the doer exactly what to report. For performance-sensitive work, ask for:

- fixture and dimensions;
- cold and warm timings;
- phase timings;
- physics values that prove the same computation is being timed;
- artifact paths;
- whether the work is ready, needs optimization, or is prototype-only.

For cleanup work, ask for:

- deleted;
- simplified;
- quarantined;
- not deleted because;
- exact remaining caller or blocker.

## Keep Physics Tests Primary

Physics tests expose convention bugs that metadata tests miss. When a pass
touches numerical meaning, include the simplest relevant physics check before a
larger endpoint.

Examples:

- He+ H1 orbital energy before He RHF;
- hydrogenic 1s self-Coulomb `5Z/8` before RHF interpretation;
- H/H2+ gausslet-only before GTO supplement;
- final-basis ordinary solve before treating a raw generalized solve as
  accepted.

Do not ask for broad helper-name or metadata-field tests unless those fields
are the active contract.

## Mention Existing Templates

For long manager/doer loops, use the file-baton template rather than writing
loop files from memory:

- `docs/src/developer/file_baton_loop_template.md`
- `.agent_handoffs/_template/`, when present locally

For ordinary chat blurbs, use this file plus the standard structure in:

- `docs/src/developer/manager_doer_collaboration_contract_2026-05-26.md`

## Formatting For Paste

When the user asks for a blurb to paste:

- introduce it with `Here is the blurb:`;
- put the blurb text after that line;
- end the blurb with `--` on a line by itself;
- put any manager commentary outside the blurb, after the delimiter or before
  the introduction;
- do not use block quotes or extra quoting that makes pasting annoying.

The delimiter matters because it makes clear where the pasted blurb ends.

## Example Skeleton

```text
Here is the blurb:

Purpose:
...

Why now:
...

Current state:
...

Use these surfaces:
...

Do not use:
...

Exact task:
...

Decision rule:
...

Record:
...

Report back:
...

--
```

The skeleton is optional. The information is not.
