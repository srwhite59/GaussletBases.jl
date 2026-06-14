# Pass 230 - Independent H2 PQS Source-Box Recovery Audit

## Purpose

Restart the loop on the real PQS/H2 problem after the fake-PQS quarantine.

The H2 463 source-backed route is now explicitly fake-PQS: it is a useful
driver golden regression for reproducing old WL/QW fixed-source data, but it is
not evidence for an independent PQS retained transform. Do not build on it as
if it were real PQS.

The next target is an independent H2 PQS route:

```text
H2 R = 4.0, q = n_s = 5
fake_pqs/enabled = false
source_backed_fixed_source_oracle_used = false
retained_transform_authority = :pqs_source_box_construction
```

## Why Now

The current running log says the active medium-term goals are:

- MT1: keep fake-PQS quarantined;
- MT2: recover independent H2 PQS;
- MT3: preserve common physical support vocabulary;
- MT4: delay MWG/GTO supplement work until retained-transform authority is clear;
- MT5: keep cleanup pressure;
- MT6: classify old Cartesian flat paths instead of promoting them by accident.

The next pass should clarify where independent PQS construction can start. It
should not resume supplement provider blocks or mutate the fake route.

## Task Type

Read-only audit. No source edits, no tests, no driver input changes, no tracked
docs except the required response copy.

## Required Reading

Read these first:

```text
AGENTS.md
BlurbStyle.md
docs/src/developer/pqs_manager_running_log.md
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/summary.md
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.229.md
```

Then inspect the relevant source surfaces. Start with:

```text
src/pqs_multilayer_shell_region_plan.jl
src/pqs_multilayer_shell_source_plan.jl
src/pqs_multilayer_support_one_body.jl
src/pqs_multilayer_support_density.jl
src/pqs_source_box_diatomic_complete_core_shell.jl
src/CartesianRouteCore.jl
src/CartesianShellification.jl
src/CartesianTerminalLowering.jl
src/CartesianRetainedUnits.jl
src/CartesianPairOperatorPlans.jl
src/CartesianPairBlockMaterialization.jl
```

If any named file is absent or has clearly moved, report the exact live
replacement path you found. Do not silently substitute a different conceptual
route.

## Questions To Answer

Answer these concretely, with file/function references:

1. What independent PQS source-box construction already exists and is still
   usable?
2. Which parts of the old 221-dimensional H2 diagnostic route are reusable, and
   which are diagnostic-only?
3. How should `atom_contact_core` be constructed without importing WL/QW
   fixed-source coefficients?
4. Can `shared_shell_1` and `shared_shell_2` be built from filled PQS source
   CPBs and PQS retained rules rather than old fixed-source transforms?
5. What route kind and driver input name should represent the real H2 PQS
   target?
6. What exact blocker prevents independent source-plan materialization?

Be especially careful around `atom_contact_core`. If the old retained count
`251` is not PQS-generated, say so. Do not normalize it as a PQS result merely
because the fake-PQS route reproduced the WL/QW 463 endpoint.

## Trust Boundary

Allowed:

- read source/tests/docs/logs;
- run `rg`, `git grep`, `git show`, `sed`, `git status`, and similar cheap
  inspection commands;
- inspect tracked artifacts or ignored `tmp/work` scripts if they clarify
  route history.

Forbidden in this pass:

- source edits;
- test edits;
- driver input edits;
- new artifacts;
- H1, H1-J, RHF, supplement, CR2, export, or public API work;
- fake-PQS route mutation except reading it as a guardrail;
- importing fake-PQS or WL/QW coefficients into the independent target;
- broad test runs.

## Decision Rules

If a route surface is ambiguous, stop at the ambiguity and report it. Do not
create a bridge or compatibility adapter.

If the only way to reach a 463-dimensional H2 endpoint is to reuse the old
source-backed WL/QW retained transform, say that the independent PQS route is
still blocked.

If you find a plausible independent source-plan seam, describe the smallest
next implementation pass. It should preferably be target/readiness only, not
final basis or H1.

## Report Back

Write:

```text
.agent_handoffs/response.230.md
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.230.md
```

Your response should include:

- concise audit verdict;
- concrete source surfaces inspected;
- answers to the six questions above;
- recommended next pass;
- carrying-cost/deletion note, even if no files were changed;
- explicit statement that no source/test/driver files were edited.

Sign:

```text
-- repo-doer@macmini
```

After writing `response.230.md`, continue polling for `blurb.231.md`,
`ATTENTION.md`, or `STOP.md`.
