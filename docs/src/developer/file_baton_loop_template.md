# File-Baton Loop Template

Use this template when setting up a manager/doer loop. Do not hand-write loop
files from memory.

The local working copy may contain an ignored helper directory at
`.agent_handoffs/_template/`. If available, copy it first:

```text
cp -R .agent_handoffs/_template .agent_handoffs/<run_id>
```

Then replace placeholders in the copied files:

```text
<run_id>
<manager>
<doer>
<max_passes>
<last_committed_hash>
<purpose>
<boundaries>
<current_boundary>
```

## Required Files

Each loop should contain:

```text
RUN.md
state.md
DOER_STARTUP.md
DOER_STARTUP_PASTE.md
RESPONSE_TEMPLATE.md
blurb.001.md
```

Longer or resumed loops may also contain:

```text
MANAGER_STARTUP.md
review.NNN.md
response.NNN.md
ATTENTION.md
STOP.md
```

## Required Doer Startup Paste

The startup paste given to the doer must include this rule directly. Do not rely
on `state.md` alone:

```text
You are <doer> in the continuous file-baton loop for this GaussletBases run:
.agent_handoffs/<run_id>/

Read:
- AGENTS.md
- JuliaStyle.md, if this pass edits Julia code
- .agent_handoffs/<run_id>/RUN.md
- .agent_handoffs/<run_id>/state.md
- .agent_handoffs/<run_id>/DOER_STARTUP.md

Then enter the polling loop described there. Act only on the next complete
blurb.NNN.md named by state.md. Ignore *.tmp files. Before each pass, check for
ATTENTION.md and STOP.md.

After completing pass NNN, write response.NNN.md.tmp and atomically rename it to
response.NNN.md.

Important: publishing response.NNN.md does not end the loop. Do not stop after
one pass. Keep polling for blurb.NNN+1.md, STOP.md, ATTENTION.md, max_passes,
or a real blocker. If the next expected baton file does not appear immediately,
poll for up to one hour before reporting that the manager did not advance the
baton.

During polling, do not request escalation for auxiliary status checks such as
process listings, Dropbox/file-provider diagnostics, or environmental hints. If
such a check would need approval, skip it and keep polling the baton files.
Request escalation only for work that is actually required by the current
implementation or validation pass; if approval is not available and the pass
cannot continue, report the concrete blocker in ATTENTION.md.

Do not self-assign tasks, do not edit manager-owned files, do not commit or
push unless the blurb explicitly says so, preserve unrelated dirty files, and
sign each response:
-- <doer>
```

## Stop Conditions

The doer should stop only when one of these appears:

- `STOP.md`
- `ATTENTION.md`
- max pass count reached after the current response is written
- a real blocker
- explicit user or manager instruction to stop

Writing one response is not a stop condition.

## Manager Checklist

After each response, the manager should:

1. Inspect `git status --short`.
2. Read `response.NNN.md`.
3. Inspect the actual diff.
4. Run risk-appropriate validation.
5. Commit/push only intended files, if appropriate.
6. Write `review.NNN.md`.
7. Update `state.md`.
8. Publish either `blurb.NNN+1.md`, `STOP.md`, or `ATTENTION.md`.

If the next pass is not ready, write `STOP.md` or `ATTENTION.md`; do not leave
the doer guessing.
