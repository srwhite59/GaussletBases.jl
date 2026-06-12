# Root Live File-Baton Loop

This document defines the no-`Go` manager/doer loop used when the manager and
doer should keep working from files until the user stops the loop or a real
design/blocker needs discussion.

This is a variant of the older file-baton loop. The key difference is that the
live baton files live directly in `.agent_handoffs/`; do not create a new
`.agent_handoffs/<run_id>/` directory for this mode.

## Purpose

Use this mode when:

- the manager should keep polling for doer responses without waiting for chat
  prompts;
- the doer should keep polling for the next blurb after each response;
- the work should continue until `STOP.md`, `ATTENTION.md`, user interruption,
  or a design decision that needs discussion;
- tracked audit records should still be preserved under
  `docs/src/developer/blurb_logs/`.

## Live Files

The live, local, disposable files are:

```text
.agent_handoffs/RUN.md
.agent_handoffs/state.md
.agent_handoffs/DOER_STARTUP.md
.agent_handoffs/DOER_STARTUP_PASTE.md
.agent_handoffs/MANAGER_STARTUP.md
.agent_handoffs/RESPONSE_TEMPLATE.md
.agent_handoffs/blurb.NNN.md
.agent_handoffs/response.NNN.md
.agent_handoffs/review.NNN.md
.agent_handoffs/ATTENTION.md
.agent_handoffs/STOP.md
```

Write in-progress files as `*.tmp` and atomically rename them to the final
filename when complete.

The durable, curated copies belong in:

```text
docs/src/developer/blurb_logs/<topic>/
```

The tracked log is for auditability. The live `.agent_handoffs/` files are the
polling mechanism and may contain local state.

## State File

`state.md` is the live source of truth for the loop state. Keep it short:

```text
status: waiting_for_doer | waiting_for_manager | attention | stopped
current_pass: 001
next_expected_file: blurb.001.md | response.001.md
tracked_log_dir: docs/src/developer/blurb_logs/<topic>
active_target: short target description
last_committed_hash: <hash>
continue_policy: continue until STOP, ATTENTION, user stop, or design blocker
```

## Manager Loop

The manager should:

1. Read `AGENTS.md`, this document, `.agent_handoffs/RUN.md`, and
   `.agent_handoffs/state.md` after reentry or compaction.
2. Publish each blurb as `.agent_handoffs/blurb.NNN.md.tmp`, then rename to
   `.agent_handoffs/blurb.NNN.md`.
3. Copy the curated blurb to the tracked log when it is important enough to
   preserve.
4. Set `state.md` to `waiting_for_doer` and `next_expected_file:
   response.NNN.md`.
5. Poll for `.agent_handoffs/response.NNN.md` or `.agent_handoffs/ATTENTION.md`.
6. When a response appears, inspect `git status`, read the response, inspect
   actual diffs/artifacts, validate at risk-appropriate scope, commit/push if
   appropriate, write `review.NNN.md`, copy the curated review to the tracked
   log, and publish the next blurb.
7. Continue polling. Do not wait for a chat `Go`.

If the next step requires design discussion, write `ATTENTION.md` and stop the
polling loop. If the user stops the loop or the work reaches a clean pause
point, write `STOP.md`.

## Doer Loop

The doer should:

1. Read `AGENTS.md`, `.agent_handoffs/RUN.md`,
   `.agent_handoffs/state.md`, and `.agent_handoffs/DOER_STARTUP.md`.
2. Poll for the exact `blurb.NNN.md` named by `state.md`.
3. Ignore `*.tmp` files.
4. Before each pass, check for `STOP.md` and `ATTENTION.md`.
5. Execute only the current blurb. Do not self-assign follow-up work.
6. Write `.agent_handoffs/response.NNN.md.tmp`, then atomically rename it to
   `.agent_handoffs/response.NNN.md`.
7. Also write the curated response copy to the tracked log if the blurb asks
   for it.
8. Continue polling for `blurb.NNN+1.md`, `STOP.md`, or `ATTENTION.md`.

Publishing one response is not a stop condition.

## Polling

Default polling:

- poll for baton files about once per minute;
- do not use process listings or environment diagnostics unless the current
  pass truly needs them;
- do not request escalation just to poll files;
- after a long wait, continue polling unless `STOP.md`, `ATTENTION.md`, or a
  user instruction changes the loop state.

If no expected baton file appears for about one hour, stop the unattended loop
and report the last known state. Do not silently convert a missing response
into manager-does-the-work unless the user explicitly asks for recovery.

## Unattended Escalation And Validation

In unattended baton mode, the doer must not request UI escalation. If a command
needs permission or sandbox escape, write `.agent_handoffs/ATTENTION.md` with
the exact command, reason, and blocker, then stop. The manager poll treats
`ATTENTION.md` as a hard stop for user/design review.

Validation commands may run longer than the one-minute file polling interval.
If validation progress is unclear for a long time, prefer a status artifact or
response checkpoint over indefinite waiting. A user may set a shorter temporary
limit for a specific run, but that does not change the standing baton rule.

## Commit And Review Rules

Only the manager commits and pushes unless a blurb explicitly delegates that.

Before commit, the manager should:

- inspect `git status --short --branch`;
- inspect the diff and relevant artifacts;
- run risk-appropriate validation;
- stage only intended files;
- commit and push;
- verify the branch is clean.

## Boundaries

The root live baton does not relax the repo policies:

- physics targets drive implementation permission;
- deletion/shrinkage reporting is required for nontrivial implementation;
- do not add tests by default;
- prefer `tmp/work` probes for exploratory audits;
- do not preserve old route authority through compatibility glue unless there
  is a live caller and the blurb says so.

## Related Documents

- `AGENTS.md`
- `BlurbStyle.md`
- `docs/src/developer/file_baton_loop_template.md`
- `docs/src/developer/blurb_logs/README.md`
- `docs/src/developer/manager_doer_collaboration_contract_2026-05-26.md`
