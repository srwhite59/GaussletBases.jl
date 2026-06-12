# MWG / Decomposed WL Cleanup Blurb Log

This log records curated manager/doer exchanges for the current decomposed
White-Lindsey, MWG, and scientific-acceptance line. It exists to make the
actual operating instructions auditable and to prevent drift back toward
additive helper/test work without deletion or shrinkage accountability.

This directory is not itself an implementation milestone. Milestone
interpretation remains in the normal developer docs, especially
`docs/src/developer/numerical_contracts.md`,
`docs/src/developer/cartesian_route_retirement_ledger.md`, and related archive
reports.

## Current Theme

The active line is driven by scientific/workflow acceptance pressure:

- decomposed White-Lindsey H/H2+/He physics checks;
- corrected IDA density-weight boundary;
- factorized decomposed WL operator construction for larger He boxes;
- residual GTO / MWG representation needed before He + GTO RHF can be valid.

The useful long-term tests should be scientific/workflow acceptance gates plus
compact module-contract tests. They should not be broad helper-vocabulary tests
or metadata inventories that only preserve development scaffolding.

## Why This Log Exists

Recent passes showed that required sections, especially deletion/shrinkage
reporting, can disappear after compaction or under task pressure. This log is a
tracked reminder that nontrivial implementation work must report whether it
deleted, simplified, or replaced old surfaces, and why any new test or helper
earns its carrying cost.

## Current Contract

Current manager blurbs in this line should include:

- a live physics or workflow target;
- explicit exclusions and trust boundary;
- exact known code surfaces when the manager has already searched;
- decision rules for stopping versus continuing;
- validation sized to risk;
- deletion/shrinkage reporting;
- test justification.

The test policy is:

- do not add tests by default;
- add tests only for live contracts, replacement/shrinkage of older coverage, or
  non-obvious bugs not covered by endpoint/module-contract tests;
- prefer `tmp/work` probes for exploratory audits;
- keep long-term tests focused on scientific/workflow acceptance and compact
  module contracts.

## Known Technical Holes And Lessons

An earlier scientific-test pressure point was honest decomposed WL H/H2+
acceptance: the route needed decomposed overlap, kinetic, and nuclear
by-center placement rather than full-window CPB or direct Cartesian fallbacks.
That kind of hole is exactly what the blurb log should preserve, because it
shows why scientific tests are better long-term guards than helper-vocabulary
coverage.

The current live He + GTO question is later in the stack:

- final-basis residual GTO directions are now represented as MWG/effective
  Gaussian data;
- raw GTO density-density is still not accepted as the final electron-electron
  operator;
- residual MWG density-density now matches the old nested fixed-block QW/MWG
  oracle on the side13 He + GTO fixture;
- the Fig. 8 AHGBS-9 S-only `n_s = 5`, `d = 0.3` reproduction now matches the
  447-function structure and is within about `0.558 mHa` of the plotted energy,
  which is close enough to stop that reproduction audit for now;
- the `n_s = 7` He RHF probe is now in the microhartree range against Fig. 8,
  so the next scientific target is diatomic H2 at the restricted closed-shell
  HF level.

## Current Cleanup Rule

Every nontrivial implementation blurb should ask what old code, test, metadata,
or compatibility path became unnecessary. If nothing was removed, the response
should say why no existing surface was made obsolete and identify any remaining
stale or duplicate surface to retire next.
