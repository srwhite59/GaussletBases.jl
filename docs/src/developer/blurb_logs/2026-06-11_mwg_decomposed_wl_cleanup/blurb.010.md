Purpose:

Start the Be atom S+P GTO physics probe. Use the existing old nested/QW atomic
route as the trusted oracle, then identify how far the newer decomposed /
final-basis route can honestly go with atomic S+P residual GTOs.

Why now:

The H2 `R = 4.0` old nested/QW diagnostic reproduced documented S+P molecular
HF rows to roundoff, so we do not need to keep proving the old diatomic route.
The active new-route supplement work has mostly proven S-only fixtures. Be is a
better next target because it is one-center, four-electron, closed-shell, and
can test atomic S+P GTO supplement/residual behavior without adding molecular
owner and geometry complications.

Known repo context:

- `docs/cr_atomic_cartesian_support_note.md` says the active atomic
  ordinary-QW Cartesian route supports `lmax = 1` and higher; the older
  one-dimensional hybrid path is still S-only and should not be used.
- `docs/src/developer/be_atom_supplement_source_facts_checkpoint_2026-06-02.md`
  records Be atom cc-pV5Z source/facts and MWG residual-component probes.
- Those Be probes use `legacy_atomic_gaussian_supplement("Be", "cc-pV5Z";
  basisfile = "/Users/srw/BasisSets")` with a machine-local basis file. Keep
  this as a developer probe, not a tracked test dependency.
- The checkpoint used `lmax = 5`; this pass should start with `lmax = 1`
  because the live question is S+P support, not high-l completeness.

Exact task:

1. Build a `tmp/work` Be atom S+P old nested/QW oracle probe.
   - atom `Be`, `Z = 4.0`, all-electron
   - closed-shell RHF with four electrons
   - `legacy_atomic_gaussian_supplement("Be", "cc-pV5Z"; lmax = 1,
     basisfile = "/Users/srw/BasisSets")`
   - use the active atomic ordinary-QW / nested-QW path, not the older
     one-dimensional hybrid branch
   - use `interaction_treatment = :mwg`
   - use `residual_keep_policy = :near_null_only`
   - start with q=5 or the closest q=5 checkpoint-compatible setup; q=6 may be
     run only if q=5 is clean and cheap enough

2. Record the oracle facts:
   - fixed dimension
   - raw S+P supplement orbital count
   - residual count and discarded count
   - final dimension
   - final overlap identity error
   - RHF one-electron, electron-electron, and total energies
   - density trace / electron count
   - residual centers and widths finite/positive
   - timing for fixed-block construction, operator construction, and RHF solve

3. Then audit the current decomposed/final-basis Be S+P path as far as it can
   honestly go.
   - Do not fake electron-electron or RHF.
   - If one-body final-basis S+P materializes, compare dimensions and one-body
     facts to the old QW oracle where meaningful.
   - If residual MWG electron-electron or final-basis interaction assembly is
     missing for S+P, report the exact blocker.
   - The likely acceptable outcome is a precise blocker, not a forced RHF.

Trust boundary:

Use `tmp/work` probes only. Do not add an acceptance test. Do not add public
APIs, exports, driver behavior, production route defaults, PQS paths, ECP
paths, raw GTO density-density final operators, generalized final-basis solves,
or full-parent/direct Cartesian fallbacks. Do not silently substitute a
repo-vendored H/He basis if `/Users/srw/BasisSets` is missing Be cc-pV5Z; stop
and report the missing basisfile blocker.

Decision rules:

- If the old nested/QW S+P oracle cannot be built, stop and localize whether
  the blocker is basisfile availability, fixture construction, supplement
  loading, QW operators, MWG interaction, or RHF convention.
- If the old oracle builds but the new route cannot proceed, report the exact
  new-route blocker and do not patch around it with a fallback.
- If the new route reaches a one-body final-basis result, compare it against
  the oracle at the level the code actually supports.
- Do not run q=6 until q=5 is clean and the q=6 run is justified by a specific
  question.

Artifacts:

- `tmp/work/be_atom_sp_qw_oracle_probe.jl`
- `tmp/work/be_atom_sp_qw_oracle_summary.txt`
- optional TSV if multiple q values or route variants are run

Validation:

- run the probe;
- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`;
- if production code changes are made for a tiny unambiguous bug, run the most
  local relevant existing test and explain why.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- artifact paths;
- exact old oracle fixture and whether `/Users/srw/BasisSets` was used;
- S+P supplement orbital count and residual count;
- RHF energy components and total;
- timing;
- how far the new decomposed/final-basis Be S+P path got;
- exact blocker if the new path cannot yet run valid RHF;
- validation run;
- deletion/shrinkage report.
