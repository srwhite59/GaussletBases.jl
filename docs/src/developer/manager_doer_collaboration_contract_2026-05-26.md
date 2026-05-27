# Manager-Doer Collaboration Contract

## Purpose

Record the working pattern that has been effective when a manager agent directs
doer-agent implementation under loose user supervision.

This note is meant to be reusable by future manager/doer pairs. It is about
communication, handbacks, validation judgment, and user-facing summaries. It is
not specific to one code path.

## Roles

The user sets priorities, corrects direction, and decides when a lane is worth
continuing. The user is supervising loosely, not managing every implementation
choice.

The manager owns direction and integration:

- reads doer handbacks
- checks the actual repository state
- decides whether the work stayed inside the requested boundary
- chooses appropriate validation
- commits and pushes when the change is ready and the tree can be staged safely
- gives the next doer blurb unless user attention is needed first

The doer owns execution:

- implements the scoped task
- preserves unrelated dirty files
- reports what changed, what was validated, and what remains out of scope
- does not quietly widen the trust boundary

Default handback signatures are part of the role contract:

- manager handbacks use `-- repo-manager@<host>`
- implementation doer handbacks use `-- repo-doer@<host>`

Use a more specific role name only when the startup prompt, role identity file,
or user explicitly gives one. A terminal nickname such as "repo-doer2" does not
by itself require a different signature if the role is still the implementation
doer.

The signature goes on messages that hand control back to the user, not on
progress updates while the agent is still continuing work. It must be the final
line.

The manager should not turn every handback into a full technical lecture for
the user. The manager is translating doer execution into a clear view of where
the work stands.

## Standard cadence

When the user pastes a doer handback, the expected manager flow is:

1. Read the handback and inspect the repo state.
2. Review the actual diff or artifact, not just the prose report.
3. Run validation sized to the risk of the change.
4. Commit and push the intended changes when they are ready, if that is within
   the current operating agreement.
5. Give a short user-facing summary.
6. Immediately give the next doer blurb unless something needs the user's
   attention.

The short summary should usually be about ten lines or less. It should answer:

- what changed
- why it matters
- what confidence we have
- what remains next

It should not lead with file lists, diff shape, or git mechanics unless those
are the point of the handback.

## Summary style

The user prefers summaries that are plain and status-oriented.

Good summary style:

- explain the work in user-facing language
- state where the line now stands
- say why the step matters
- mention validation briefly
- keep git details secondary
- do not assume the user knows the code details

Less useful summary style:

- too many file names
- too much diff language
- too many internal helper names
- long lists of implementation details before the bottom line
- treating every handback like a code-review report

An effective summary shape is:

> We now have a safe wrapper path around the existing construction code. It
> does not replace the working builders; it records what was requested, runs the
> existing builder, and checks that the result matches the intended route. The
> next step is to extend that wrapper to the next route family.

That style gives the user a clear management view without requiring them to
reconstruct the implementation from internal names.

## Standard doer blurb format

When issuing the next task to a doer, use the standard blurb structure:

1. Purpose
2. Why now
3. Current state
4. Exact task
5. Trust boundary
6. Artifacts to use
7. What to record
8. Decision rule
9. Report back

The blurb should be concrete enough for a doer to execute without guessing, but
it should keep the trust boundary explicit. The trust boundary is especially
important when plumbing work is near numerical kernels or scientific policy.

## Commit and push convention

The user generally delegates commit/push judgment to the manager.

The manager may commit and push directly when:

- the change has been reviewed
- validation is appropriate for the risk
- only intended files are staged
- unrelated dirty files are preserved
- the change is not waiting on a user decision

If the situation is ambiguous, the manager should either hold the commit or make
the next blurb include explicit commit commands for the doer.

## Validation judgment

Do not run broad slow tests automatically just because a doer reported them.
Tests in this repo can take a long time.

The manager should choose validation based on the changed surface:

- metadata-only or audit-only changes can usually be checked with parse/load,
  whitespace, and targeted smoke tests
- changes touching shared numerical behavior need broader regression coverage
- changes that add a heavy production path need performance and allocation
  evidence, not just tiny correctness tests

If the manager intentionally does not rerun a broad test group, say that
clearly and treat the doer's reported full pass as supporting evidence rather
than as locally reproduced evidence.

## Trust-boundary reminders

The manager should keep repeating the important boundary in each doer blurb.
Examples from the current repo style:

- no new Hamiltonian kernels unless explicitly requested
- no backend default changes unless explicitly requested
- no geometry-policy changes unless explicitly requested
- no dense parent matrices or metric packets in metadata/audit plumbing
- preserve PGDG analytic-integral routes
- numerical quadrature only on explicit reference or diagnostic paths
- preserve unrelated dirty files

These reminders prevent a doer from turning a plumbing step into a numerical
policy change.

## When to stop and ask the user

Most handbacks should flow directly into review, commit, and the next blurb.
Pause for the user when:

- the doer widened the trust boundary
- validation fails in a way that changes the decision
- there is a scientific or policy choice rather than an engineering cleanup
- the repo state makes safe staging unclear
- the next step depends on priority rather than sequence

Otherwise, keep momentum. The manager should not ask for routine confirmation
after every successful small handback.

## Reusable bottom line

The best manager behavior in this pattern is:

- translate doer work into a clear status picture
- guard the trust boundary
- validate enough, not maximally
- integrate clean work without unnecessary delay
- give the next executable blurb
- keep the user oriented at the level of direction and confidence
