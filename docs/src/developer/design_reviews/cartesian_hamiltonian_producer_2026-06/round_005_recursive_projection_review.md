# Round 005 Recursive Projection Freeze Review

Reviewer: repo-manager@atlas recommendation, pasted by user

Date: 2026-06-20

## Verdict

Do not bind the Slice A IDs at `2325b846` yet. Run one final, tightly scoped,
uncommitted recursive-projection spike, then freeze without another broad review
round if it passes.

## Main Findings

- The previous spike used shell-local approximations for previous PQS blocks,
  not recursively projected coefficients and effective supports.
- Later direct records were not explicitly reported against earlier PQS blocks.
- The projection algorithm should avoid numerical support pollution:
  if `norm(C_previous' * S * X, Inf) <= projection_atol`, do not subtract and do
  not enlarge effective support.
- Projection and audit should compute block cross actions incrementally:
  `C_left' * S_lr * C_right`, with at most one support-pair workspace live.

## Required Design Corrections

1. Split sign ownership:
   - `HP-FN-00` performs projection plus shell-local Lowdin only.
   - `HP-FN-01` derives weights, validates direct weights, sign-canonicalizes
     completed block columns, and constructs `CartesianTerminalBasisBlock`.
2. Specify the overlap-action contract explicitly.
3. Remove per-object/helper line targets and keep only the Slice A budget:
   target `150`, redesign threshold `225`, net source decrease required.
4. Require the reviewed Cr2 fixture to produce a real terminal basis; distorted
   COMX rejection is allowed only when the typed transform inventory actually
   contains distorted COMX.

## Final Spike Questions

1. Does true recursive projection preserve locality?
2. Do later direct records remain orthogonal to previous direct/PQS blocks?
3. Does the projection threshold prevent roundoff-only support pollution?
4. Is pairwise/incremental overlap fast enough for one-center, H2, and Cr2?
