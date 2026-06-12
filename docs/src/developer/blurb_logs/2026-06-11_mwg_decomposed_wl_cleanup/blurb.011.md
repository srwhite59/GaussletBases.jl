Purpose:

Retry the Be atom S+P GTO physics probe using the specific Be `cc-pV5Z` basis
source that is actually present on this machine.

Why now:

The previous pass correctly stopped because `/Users/srw/BasisSets` was absent.
That was a path mismatch, not absence of Be `cc-pV5Z` data. Manager review
found Be `cc-pV5Z` in the GaussletModules basis files. This follow-up
authorizes one exact alternate source so the old nested/QW oracle can be built
without asking doer to guess.

Authorized basisfile:

Use:

```julia
basisfile = "/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/BasisSets"
```

This file contains `#BASIS SET: Be  cc-pV5Z`. Do not use `/Users/srw/BasisSets`
in this retry. Do not search for or substitute another basis source unless this
exact authorized file is missing, in which case stop and report.

Exact task:

1. Rerun the Be atom S+P old nested/QW oracle probe.
   - atom `Be`, `Z = 4.0`, all-electron
   - closed-shell RHF with four electrons
   - `legacy_atomic_gaussian_supplement("Be", "cc-pV5Z"; lmax = 1,
     basisfile = "/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/BasisSets")`
   - use the active atomic ordinary-QW / nested-QW path, not the older
     one-dimensional hybrid branch
   - use `interaction_treatment = :mwg`
   - use `residual_keep_policy = :near_null_only`
   - start with q=5 or the closest q=5 checkpoint-compatible setup; q=6 may be
     run only if q=5 is clean and the q=6 question is explicitly useful

2. Record the old oracle facts:
   - fixed dimension
   - raw S+P supplement orbital count
   - residual count and discarded count
   - final dimension
   - final overlap identity error
   - RHF one-electron, electron-electron, and total energies
   - density trace / electron count
   - residual centers and widths finite/positive
   - timing for fixed-block construction, operator construction, and RHF solve

3. Only after the old oracle is clean, audit the current decomposed/final-basis
   Be S+P path as far as it can honestly go.
   - Do not fake electron-electron or RHF.
   - If one-body final-basis S+P materializes, compare dimensions and one-body
     facts to the old QW oracle where meaningful.
   - If residual MWG electron-electron or final-basis interaction assembly is
     missing for S+P, report the exact blocker.
   - A precise blocker is an acceptable outcome.

Trust boundary:

Use `tmp/work` probes only. Do not add an acceptance test. Do not add public
APIs, exports, driver behavior, production route defaults, PQS paths, ECP
paths, raw GTO density-density final operators, generalized final-basis solves,
or full-parent/direct Cartesian fallbacks. Do not broaden this into high-l
Be, Be2, Cr, or H2 work.

Decision rules:

- If the authorized GaussletModules basisfile is missing or cannot load Be
  `cc-pV5Z`, stop and report the exact basisfile blocker.
- If the old nested/QW S+P oracle cannot be built, stop and localize whether
  the blocker is fixture construction, supplement loading, QW operators, MWG
  interaction, or RHF convention.
- If the old oracle builds but the new route cannot proceed, report the exact
  new-route blocker and do not patch around it with a fallback.
- Do not run q=6 until q=5 is clean and q=6 answers a named question.

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
- exact basisfile used;
- S+P supplement orbital count and residual count;
- RHF energy components and total;
- timing;
- how far the new decomposed/final-basis Be S+P path got;
- exact blocker if the new path cannot yet run valid RHF;
- validation run;
- deletion/shrinkage report.
