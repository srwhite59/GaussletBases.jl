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

The He + GTO and Be S+P work later in the stack established three durable
lessons. First, final-basis residual GTO directions should be represented as
MWG/effective Gaussian data; raw GTO density-density must not be accepted as the
final electron-electron operator. Second, old nested/QW can be trusted as an
oracle when it reproduces known He/H2 rows or Be S+P values, but it should not
be treated as the target architecture. Third, the newer decomposed/final-basis
route became scientifically credible only after matching the old Be S+P oracle
to roundoff and replacing per-unit mixed CPB/GTO local construction with a
factorized retained-basis path.

The final state of this loop is performance debt, not physics uncertainty. Be
S+P with `q/ns = 5 / 5` has fixed gausslet dimension `615`, `21` residual S+P
supplement directions, final dimension `636`, and RHF total
`-14.574514244574694` from the old oracle. The decomposed/final-basis route
matches that oracle to roundoff; warm route construction is subsecond, but cold
fresh-process latency still reflects compilation and specialization pressure,
especially moment-capable mixed GTO payloads and density prerequisites.

## Current Cleanup Rule

Every nontrivial implementation blurb should ask what old code, test, metadata,
or compatibility path became unnecessary. If nothing was removed, the response
should say why no existing surface was made obsolete and identify any remaining
stale or duplicate surface to retire next.

## Consolidated Pass Summary: 001-025

Passes 001-025 moved the MWG/decomposed-WL line from residual-GTO plumbing into
a validated Be S+P atom route, then paused after compile-attribution work. The
early passes added combined GTO residual moment matrices and final-basis
residual MWG density-density blocks while preserving the trust boundary: raw
GTO density-density was not accepted as the final electron-electron operator.
Side13 He + GTO RHF improved over gausslet-only He but sat about `2.7 mHa`
below the exact He HF reference, so it stayed diagnostic until matched old
nested/QW MWG comparison showed roundoff agreement with the decomposed route.
That established the issue as Hamiltonian/fixture semantics, not a decomposed
route bug.

The loop then used old nested/QW as oracle rather than route authority. The
Fig. 8 AHGBS-9 He audit recovered the `n_s = 5`, `d = 0.3`, 447-function
structure but remained `0.558 mHa` below the plotted energy, while the `n_s = 7`
He points reached microhartree-to-tens-of-microhartree agreement, enough to stop
the atomic-He reproduction audit. H2 `R = 4.0` old nested/QW restricted HF
reproduced documented S+P molecular rows to roundoff, making the old diatomic
route trustworthy as an oracle. The first Be S+P attempt correctly blocked on a
missing `/Users/srw/BasisSets`; the authorized GaussletModules basis retry
produced the old nested/QW oracle with final dimension `636` and RHF total
`-14.574514244574694`.

Passes 012-015 removed the key replacement blocker,
`:missing_driver_owned_decomposed_be_sp_fixture_wiring`, by adding a private
decomposed atom+GTO final-basis seam. The one-electron Be S+P route materialized
with retained gausslet dimension `615`, `21` supplement directions, final
dimension `636`, and clean final overlap/Hamiltonian checks. The full
final-basis density-density/RHF route then matched the old nested/QW oracle to
roundoff, with RHF total near `-14.57451424457464`. Phase timing identified
`mixed_gto_blocks` as the real cost center; hoisting GTO/GTO self blocks helped
only modestly, and the decisive fix was a factorized retained-basis mixed-GTO
path for the one-center atomic case, reducing `mixed_gto_blocks` from about
`177s` to about `9.18s` while preserving oracle agreement.

Passes 016-022 separated algorithmic runtime from cold compilation. Warm/cold
attribution showed the post-mixed-GTO Be route was about `170s` cold but only
about `0.52s` warm, with total warm time dominated by RHF rather than operator
construction. A small synthetic atom+GTO precompile workload was added, but its
benefit was limited, so the loop turned to specialization shape instead of
enlarging precompile. The hot decomposed-WL inventory result was changed from
tuple-sized summaries to vector-backed `unit_keys` / `unit_summaries` and
compact `pair_summaries`, cutting Be fresh-process route time from about `159s`
to about `31s`. Narrow matrix-set compute objects then replaced report-shaped
residual moment and one-electron staging, improving residual-moment cold phase
from about `6.91s` to `1.07s` and cold route time to about `25.5s`, with
physics still pinned to the Be oracle.

Passes 023-025 deliberately stopped short of speculative refactoring. Direct
one-electron helper timing did not justify rewriting inner factorized helpers;
full-route attribution showed timing sinks and metadata were not the cause, and
density-density-enabled work accounted for about `13.5s` of the remaining cold
route gap. The final density result constructor and old residual-MWG kernel
were not bottlenecks; remaining pressure was traced to density prerequisites,
especially moment-capable mixed GTO block payload shape, residual moment
prerequisites, and gausslet density construction. The recommended state after
pass 025 was to record this as known GTO compile debt, avoid further polish
unless it blocks the next physics target, and pivot to the PQS source-box-first
plan with the baton-loop mechanics in place.

The detailed per-pass files for passes 001-025 were removed after this
consolidation to keep the curated log from accumulating obsolete scaffolding.
