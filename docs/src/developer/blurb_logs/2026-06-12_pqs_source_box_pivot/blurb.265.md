Pass 265 - audit remaining PQS blurb-log compression windows

Context:
- Current HEAD should include
  `e6402d5f Retire legacy source-box fixture builders`.
- Passes 263 and 264 removed 6,980 scoped `src + test + bin` lines from old
  current-route metadata export and legacy source-box fixture-builder
  scaffolding.
- A separate remaining bloat surface is tracked documentation history:
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot`.
- Passes 1-100 have already been consolidated into summary form, but raw
  blurb/response/review logs after pass 100 still account for a large tracked
  documentation surface.
- Preliminary manager-side inventory found roughly:
  - passes 101-150: 150 raw files, about 12.8k lines;
  - passes 151-200: 145 raw files, about 16.2k lines;
  - passes 201-229: 82 raw files, about 9.9k lines;
  - passes 230-current: 103 raw files, about 10.3k lines.
  Please verify these counts rather than assuming them.
- The recent active tail, especially roughly passes 230-current, is still
  operationally useful because it records fake-PQS quarantine, independent H2
  PQS recovery, final-basis/H1/H1-J/RHF diagnostics, supplement preflight, and
  provider-block support-partition staging.

Task:
Do a read-only docs-history compression audit for the PQS blurb-log directory.

Directory to audit:
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot`

Required inventory:
- Count tracked `blurb.NNN.md`, `response.NNN.md`, and `review.NNN.md` files
  by pass ranges:
  - `101-150`
  - `151-200`
  - `201-229`
  - `230-current`
- Compute line counts by range and file kind.
- Identify the latest/current raw logs that should remain intact for active
  operations.
- Identify any current docs, state files, or recent blurbs that explicitly
  depend on raw logs from old ranges.

Compression proposal:
- Propose summary files to create, for example:
  - `pqs_passes_101_150_summary.md`
  - `pqs_passes_151_200_summary.md`
  - `pqs_passes_201_229_fake_pqs_leadup_summary.md`
  or better names/windows if the audit suggests them.
- List exact raw files likely safe to delete after those summaries exist.
- Separate files that are safe-to-delete-after-summary from files that should
  remain full for now.

Guardrails for summaries:
- Preserve durable decisions, route-authority changes, accepted dimensions and
  counts, blocker transitions, validation landmarks, deletion/paydown events,
  and what future agents must not repeat.
- Especially preserve:
  - fake-PQS H2 463 demotion to source-backed WL/QW reproduction only;
  - fake-PQS retained counts `(251, 98, 114)` are not independent PQS
    authority;
  - independent H2 PQS retained counts `(275, 98, 98)` and final dimension
    `471`;
  - independent route must use `fake_pqs/enabled = false` and
    `retained_transform_authority = :pqs_source_box_construction`;
  - gausslet-only results must not be compared to supplemented WL/QW
    references;
  - provider-block work requires explicit support tiling/row ownership before
    CPB provider calls.
- Assume active tail `230-current` should remain full unless the audit gives a
  concrete reason otherwise.

Strict exclusions:
- Do not edit files.
- Do not delete logs.
- Do not create summaries in this pass.
- Do not touch source, tests, driver inputs, support-partition code, provider
  blocks, supplements, CR2/export, or public API.
- Do not chase untracked `.agent_handoffs` as an archival source of truth; the
  target is tracked docs history.

Validation:
- No Julia validation is needed.
- If no tracked files are edited, report validation as not applicable.
- If you use a helper command that writes nothing, include enough command detail
  for manager to reproduce the inventory.

Report:
- Counts of tracked files by pass range and kind.
- Line counts by pass range and kind.
- Current/latest raw logs to keep intact.
- Proposed summary files to create.
- Exact raw files likely safe to delete after each summary.
- Risks and guardrails.
- Deletion/shrinkage planning card:
  - deletion candidate:
  - why obsolete:
  - live source callers / docs dependencies:
  - expected line savings:
  - risk class: green / yellow / red:
  - minimum validation for the later deletion pass:
  - what not to run:
  - replacement/current authority:

-- repo-manager@macmini
