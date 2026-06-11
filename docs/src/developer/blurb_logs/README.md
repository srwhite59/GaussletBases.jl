# Blurb Logs

This directory holds curated, tracked records of important manager/doer
instructions and results. It is for auditability: future agents, GitHub
connector reviews, and archive reviews should be able to inspect the actual
task blurbs, doer responses, and manager reviews without relying on chat paste
history.

Use this directory for selected nontrivial work only. It is not a raw transcript
dump and it is not a replacement for durable design docs, retirement ledgers,
numerical-contract notes, or archive reports.

## Relationship To `.agent_handoffs`

`.agent_handoffs/` is the local live baton/polling mechanism. It may contain
`.tmp` files, runner state, stale polling artifacts, local paths, and
machine-specific coordination details. It should remain local and mostly
untracked.

`docs/src/developer/blurb_logs/` is the tracked curated record. It should contain
only reviewed, useful records that are worth preserving in the repository.

Do not move `.agent_handoffs` contents here mechanically. Curate first, exclude
local clutter, and preserve only material that helps future review.

## File Shape

For a logged work line, create a directory such as:

```text
docs/src/developer/blurb_logs/YYYY-MM-DD_short_topic/
```

Useful files are:

```text
summary.md
blurb.NNN.md
response.NNN.md
review.NNN.md
```

`summary.md` should explain why the line exists, what live contract it protects,
and where milestone interpretation lives in the normal docs.

Use `blurb.NNN.md`, `response.NNN.md`, and `review.NNN.md` only for curated
records. They do not need to include every small pass.

## What To Include

Include records when they preserve an important boundary, such as:

- a physics target that controls implementation scope;
- a live numerical convention such as IDA weight placement;
- a performance/scaling target that prevents a wrong algorithmic shape;
- a cleanup/deletion target where drift would preserve obsolete code;
- a manager correction that future agents are likely to forget.

Each implementation blurb should include or explicitly waive:

- purpose and live target;
- current state;
- exact task;
- trust boundary and exclusions;
- artifacts or code surfaces to inspect;
- decision rules;
- validation;
- reporting requirements;
- deletion/shrinkage reporting;
- test policy.

Named files, functions, fixtures, and artifacts are starting points, not
authority. If they are missing or inconsistent with the live repo, the doer
should stop and report the mismatch instead of silently substituting a different
path.

## What Not To Include

Do not include:

- `.tmp` files;
- raw huge stdout;
- local runner state;
- stale polling loops;
- machine-only clutter;
- large artifacts;
- private or sensitive material;
- generated runtime environments;
- material that belongs in a normal design doc or numerical contract instead.

Large timing tables or probe outputs should generally stay under `tmp/work/`
or another artifact location and be summarized here only when they matter.

## Required Deletion And Shrinkage Section

Every nontrivial implementation blurb should include this section, or explicitly
say why it does not apply:

```text
Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.
```

This is required because additive work can hide conceptual drift. A pass that
adds source, tests, metadata, adapters, or docs should also account for what it
made simpler, obsolete, or still unresolved.

## Test Policy

Do not add tests by default.

Add a test only if it:

- protects a live contract;
- replaces or shrinks older coverage;
- catches a non-obvious bug not covered by an endpoint or compact
  module-contract test.

Prefer `tmp/work` probes for exploratory audits. Long-term tests should
preferentially be scientific/workflow acceptance gates plus compact
module-contract tests. Avoid helper-vocabulary tests, broad metadata tests, and
tests that preserve obsolete paths unless the blurb names them as temporary
oracle/reference coverage.

## Carrying Cost

Source, tests, docs, metadata, compatibility glue, adapters, artifacts, and
blurb-log entries all have carrying cost. New artifacts should earn that cost
by protecting a live contract, improving clarity, reducing duplication,
improving performance, or enabling a current workflow.

When two approaches have similar value, prefer the one that leaves fewer stale
surfaces. For cleanup and retirement tasks, ask what old code or tests become
removable rather than only what new coverage should be added.

## Response Expectations

A doer response should report:

- what changed;
- what was validated;
- what did not change;
- whether the work stayed inside the trust boundary;
- blockers or next exact missing pieces;
- deletion/shrinkage status;
- test additions and why they were justified.

If no deletion happened, say why no existing surface was made obsolete.

## Review Expectations

A manager review should inspect actual files and artifacts, not only the doer
prose. It should record:

- whether the change matched the blurb;
- validation run;
- any corrections made before commit;
- commit/push status when applicable;
- the next live target or stop condition.

## Milestone Interpretation

Milestone interpretation belongs in normal developer docs, numerical-contract
notes, retirement ledgers, or archive reports. The blurb log can point to those
docs, but it should not be the only place where a scientific result, route
contract, or retirement decision is explained.

## Related Guidance

Use these repo-local references with this directory:

- `AGENTS.md`
- `JuliaStyle.md`
- `BlurbStyle.md`
- `docs/src/developer/manager_doer_collaboration_contract_2026-05-26.md`
- `docs/src/developer/file_baton_loop_template.md`
- `docs/src/developer/cartesian_route_retirement_ledger.md`

## Lightweight Templates

### `blurb.NNN.md`

```text
Purpose:

Why now:

Current state:

Exact task:

Trust boundary:

Artifacts/code surfaces:

Decision rules:

Validation:

Deletion/shrinkage report required:

Report back:
```

### `response.NNN.md`

```text
What changed:

Validation:

Deletion/shrinkage:

Tests:

Blocked or next:

Final status:
```

### `review.NNN.md`

```text
Review result:

Corrections made:

Validation:

Commit/push:

Next target:
```
