Purpose:

Build an old nested/QW MWG oracle comparison for the side13 He + GTO
final-basis diagnostic, and use it to decide whether the below-HF result is an
expected MWG/IDA approximation effect or a convention mismatch.

Why now:

The focused convention audit found no clear bug in the new decomposed
final-basis RHF path. The result still cannot become acceptance because it is
below the He HF reference by about `2.7 mHa`. The next useful comparison is not
another acceptance test; it is a matched or intentionally matched old QW/MWG
oracle. White-Lindsey Fig. 8 is also relevant: the manager recollection is that
for `ns = 5`, `d = 0.3`, 447 basis functions, and a very accurate S-only GTO
supplement, Fig. 8 shows an error of a few tenths of a mH, with larger error
for smaller `d`.

Current state:

- new decomposed side13 He + GTO diagnostic:
  - `d = 0.1`, `s = 1.0`, side count `13`
  - final dimension `520`
  - RHF total `-2.864394277839975`
  - below He HF reference by about `2.7 mHa`
- side13 gausslet-only RHF total: `-2.8364979997009137`
- no repo-local Fig. 8 numeric table has been found yet
- old QW/MWG reference code exists and should be used only as oracle/reference,
  not as new route authority

Named starting points:

- `docs/src/algorithms/qiu_white_residual_gaussian_route.md`
- `docs/ordinary_cartesian_qiu_white_reference.md`
- `src/ordinary_qw_operator_assembly.jl`
- `ordinary_cartesian_qiu_white_operators(...)`
- `_qwrg_mwg_interaction_components`
- `_qwrg_final_residual_mwg_component_blocks`
- `test/ordinary/runtests.jl`, especially:
  - nested fixed-block QW/MWG coverage around the existing
    `ordinary_cartesian_qiu_white_operators(fixed_block, supplement; ...)`
    calls;
  - the "Legacy He s MWG residual interaction" block;
  - the "Qiu-White MWG pair-density normalization matches IDA convention"
    block.

Named surfaces are starting points, not authority. If any named surface is
missing or inconsistent with the live repo, stop and report the actual live
surface before implementing.

Exact task:

Create a `tmp/work` developer probe that compares the new decomposed side13
He + GTO final-basis diagnostic against the old nested/QW MWG reference path.

The probe should try, in this order:

1. A matched side13 oracle:
   - same He center and `Z = 2`;
   - same or intentionally matched side13 gausslet shell/fixed-block geometry;
   - same or intentionally matched S-only GTO supplement;
   - `interaction_treatment = :mwg`;
   - ordinary/QW final-basis RHF or the same scalar/RHF diagnostic used by the
     new probe.

2. A paper-like Fig. 8 oracle if the repo can construct it cleanly:
   - `ns = 5`;
   - `d = 0.3`;
   - approximately 447 basis functions if that follows from the available
     constructor;
   - very accurate S-only GTO supplement if available from existing basis
     plumbing;
   - compare the resulting error with the manager's Fig. 8 lead.

If either fixture cannot be constructed cleanly from existing public/internal
surfaces, do not build a new framework. Report the exact blocker and the
closest available intentionally matched fixture.

Compare and record:

- fixture parameters and basis dimensions;
- residual counts;
- final overlap identity error;
- one-electron energy contribution;
- electron-electron contribution;
- RHF total or the closest existing old-QW scalar/RHF diagnostic;
- difference between old QW/MWG and new decomposed MWG for the matched fixture;
- whether both are below/above the exact He HF reference;
- whether the Fig. 8 trend about smaller `d` worsening error is consistent with
  the available data;
- timing at a coarse level only.

Decision rules:

- If old QW/MWG agrees with the new decomposed result for the matched fixture,
  report that the below-HF behavior is likely MWG/IDA approximation behavior,
  not a new decomposed-route bug. Do not add an acceptance test yet unless the
  manager explicitly asks.
- If old QW/MWG and new decomposed MWG disagree materially, localize the
  difference to one of:
  - final-basis one-electron matrix;
  - residual centers/widths;
  - `V_gg`;
  - `V_gR`;
  - `V_RR`;
  - RHF convention.
- If the Fig. 8 source or fixture is not available, say so clearly and do not
  treat the recollected numbers as a citation.

Trust boundary:

Use `tmp/work` probes first. Do not add a new acceptance test. Do not change
production code unless the probe exposes a narrow, clear bug. Do not use old
fixed-block matrices as the new route authority. Do not introduce raw GTO
density-density final operators, full-parent CPB fallback, direct Cartesian
product fallback, ordinary Cartesian IDA fallback, PQS, exports, or a
generalized final-basis solve.

Validation:

- run the probe;
- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`;
- if production code changes are made for a clear bug, run the most local
  relevant existing test and explain why.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- artifact paths;
- whether a matched old QW/MWG oracle was constructed;
- whether a Fig. 8 source or reproducible proxy was found;
- old-vs-new energy/component comparison;
- whether the below-HF side13 result now looks expected, erroneous, or still
  unresolved;
- validation run;
- deletion/shrinkage report.
