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
  closing the immediate atomic He accuracy check;
- the H2 `R = 4.0` old nested/QW restricted HF probe reproduced documented
  S+P molecular rows to roundoff, so the old diatomic route is trusted enough
  as an oracle for this line;
- the next scientific target is Be atom S+P GTO residual behavior: old
  nested/QW should be used as oracle, while the live question is how far the
  newer decomposed/final-basis route can honestly go and what precise blocker
  remains before valid Be RHF;
- the first Be S+P attempt correctly stopped because `/Users/srw/BasisSets`
  was absent; the authorized retry should use
  `/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/BasisSets`, which
  contains Be `cc-pV5Z`;
- the authorized Be S+P old nested/QW oracle is now clean at q/ns `5 / 5`,
  with fixed dimension `615`, `21` residual S+P supplement directions, final
  dimension `636`, and RHF total `-14.574514244574694`;
- the next replacement-oriented blocker is
  `:missing_driver_owned_decomposed_be_sp_fixture_wiring`: removing it should
  turn H/H2+ fixture-local GTO wiring into a reusable decomposed
  atom+supplement seam and make old nested/QW less necessary as route
  authority;
- the private decomposed atom+supplement seam now materializes the Be S+P
  final-basis one-electron route with final dimension `636`; the next blocker
  is not old-oracle data, but phase-attributed final density-density/RHF
  materialization through that seam;
- the Be S+P final-basis RHF route now matches the old nested/QW oracle to
  about `3.2e-14 Ha`; the measured cost center is `mixed_gto_blocks`, about
  `188.6` seconds out of a `357.4` second probe;
- after hoisting reusable GTO/GTO self blocks, the Be S+P probe still matches
  the old oracle to about `3.2e-14 Ha`, while total time is about `342.5`
  seconds and `mixed_gto_blocks` is about `177.2` seconds; subphase timing
  shows the true remaining bottleneck is per-unit mixed CPB/GTO local block
  construction, about `168.1` seconds over `131` retained units;
- the active one-center atomic mixed-GTO route now has a factorized/projected
  retained-basis path, removing the per-unit provider-local hot loop for Be
  S+P; the Be probe still matches the old nested/QW oracle to about
  `5.5e-14 Ha`, total time is about `174.3` seconds, and `mixed_gto_blocks`
  is about `9.18` seconds;
- warm/cold attribution after the mixed-GTO replacement shows the remaining
  long cold Be timings are mostly compilation: cold route is about `169.9`
  seconds, but same-process warm route is about `0.523` seconds and warm total
  is about `2.40` seconds, dominated by RHF rather than route construction;
- a small synthetic atom+GTO precompile workload now covers the final-basis
  route shape without user-local `BasisSets`; it costs about `49s` to rebuild
  package precompile and improves the Be S+P fresh-process route from about
  `169.9s` to `159.0s`, so the next cold-latency target is specialization
  pressure rather than a larger precompile fixture;
- the first specialization audit found that hot decomposed WL inventory results
  still encode retained-unit count in the type via tuple-valued `unit_keys` and
  `unit_summaries`; the side-7 precompile workload compiles a 27-unit result,
  while side-15 Be needs a 131-unit result. Small-inventory pair summaries also
  compile a tuple-heavy shape that production side-15 omits.

## Current Cleanup Rule

Every nontrivial implementation blurb should ask what old code, test, metadata,
or compatibility path became unnecessary. If nothing was removed, the response
should say why no existing surface was made obsolete and identify any remaining
stale or duplicate surface to retire next.
