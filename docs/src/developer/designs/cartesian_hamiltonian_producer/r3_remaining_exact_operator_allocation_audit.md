# R3 Remaining Exact-Operator Allocation Audit

Status: approved measurement-only authority. This is not production source
authority and does not extend `HP-R3GG-FN-01`.

## Decision

After `954c86cd`, the terminal final-basis `G-G` product-workspace lane is
crossed. The Cr2 q4 exact augmented-operator wrapper is now measured at about
`5.8389s / 4605.517 MiB`. The remaining allocation is reported mainly outside
the nine terminal product buffers, including unit-nuclear `U_GG` Gaussian-sum
work and route/raw-block setup.

That evidence is not yet specific enough to approve another source lane. The
next authority step is a narrow measurement-only audit that separates the
remaining allocation buckets and identifies the exact owner/file/function for a
possible later amendment.

Do not send implementation work into unit-nuclear `U_GG`, route/raw-block
setup, or any continuation of `HP-R3GG-FN-01` until this audit has produced a
specific target and a separate docs-only source amendment approves it.

## Approved ID

- `HP-R3REM-AUDIT-01` - measurement-only audit of remaining R3/RG exact
  augmented-operator allocation after terminal `G-G` product workspace reuse.

`HP-R3REM-AUDIT-01` approves only ignored measurement/probe work. It does not
authorize production source edits, committed tests, new helper files, public
API, artifacts, metadata, report/status/payload fields, route cleanup, or Cr2
facade/artifact workflow.

## Measurement Scope

The audit should use Cr2 q4 as the main scale proxy, with H2/Be2 sanity only if
needed to keep replay code honest.

Separate timing/allocation for:

- total exact augmented-operator wrapper;
- neutral non-nuclear `G-A`/`A-A` raw-block construction;
- neutral nuclear `G-A`/`A-A` raw-block construction;
- terminal `G-G` kinetic/moment product assemblies, confirming they remain
  crossed after `HP-R3GG-FN-01`;
- unit-nuclear `U_GG` Gaussian-sum construction by center;
- transformation of exact augmented nuclear blocks, if separable from
  `U_GG` construction;
- route/stage setup and raw-block setup not attributable to the numerical
  kernels above;
- audit/replay overhead distinct from production helper work.

The audit may create ignored `tmp/work` probes. It must not commit those probes
or use them as new production APIs.

## Required Output

The audit handoff must classify the remaining allocation into exact buckets and
recommend one of:

- approve a future unit-nuclear `U_GG` Gaussian-sum source lane;
- approve a future route/raw-block setup source lane;
- accept the current cost for now;
- run another narrower measurement because attribution is still ambiguous.

If a source lane is recommended, the handoff must name the proposed owner,
files, functions, forbidden surfaces, validation gates, line budget, expected
deletion/simplification, and failure rule. Without those details, no source
handoff is approved.

## Forbidden

This audit does not approve:

- edits to `src`, `test`, `tools`, or `bin`;
- extension of `HP-R3GG-FN-01` into unit-nuclear Gaussian sums;
- unit-nuclear `U_GG` source optimization;
- route/stage setup cleanup or raw-block setup cleanup;
- final-basis `G-G` product-matrix changes;
- `G-A`/`A-A` raw-block changes;
- residual Gaussian selection, orientation, transform, MWG, or IDA changes;
- Qiu-White semantic changes;
- parent construction or terminal basis realization changes;
- persistent caches/workspaces, metadata, report/status/payload fields,
  artifact schema changes, public API/export, Cr2 facade support, or Cr2
  artifact workflow.

## Validation For This Docs Pass

- `git diff --check`;
- focused search confirming `HP-R3REM-AUDIT-01` is measurement-only and not
  source authority;
- confirmation that no `src`, `test`, `tools`, or `bin` files changed.
