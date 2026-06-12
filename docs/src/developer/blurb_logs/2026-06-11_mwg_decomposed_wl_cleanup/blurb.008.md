Purpose:

Move on from the `n_s = 5` Fig. 8 reproduction audit and run the first
`n_s = 7` He HF diagnostic. This should test accuracy at a higher
White-Lindsey shell order while staying in the same He RHF physics family.

Why now:

The Fig. 8 `n_s = 5`, `d = 0.3`, AHGBS-9 S-only reproduction is structurally
good enough for now:

- AHGBS-9 S-only loads from the external GaussletModules `BasisSets` file.
- The old nested/QW MWG fixture reaches the expected 447-function structure:
  `419` gausslet functions plus `28` residual S directions.
- The reproduction energy is `-2.862102144533723` versus the Fig. 8 target
  `-2.861543784624258`, a difference of about `-0.558 mHa`.
- The manager considers this close enough to stop the `n_s = 5` reproduction
  audit for now.

The previous side13 cc-pVTZ diagnostic was about `2.7 mHa` below the He HF
reference, which felt large. The AHGBS-9/Fig. 8 reproduction gets the
comparison into the sub-mH regime. The next accuracy target should therefore be
`n_s = 7`, still He RHF, not more `n_s = 5` convention archaeology.

Archive target:

- data table:
  `/Users/srw/Library/CloudStorage/Dropbox/chatarchive/references/papers/canonical/White and Lindsey - 2023 - Nested Gausslet Basis Sets_fig8_heuhf_data.tsv`
- provenance note:
  `/Users/srw/Library/CloudStorage/Dropbox/chatarchive/references/papers/canonical/White and Lindsey - 2023 - Nested Gausslet Basis Sets_fig8_heuhf_provenance.md`

Useful Fig. 8 `n_s = 7` table points:

- `d = 0.10`: energy `-2.8616756776234777`, signed error
  `+4.317988761215e-06 Ha`
- `d = 0.15`: energy `-2.8616793497179782`, signed error
  `+6.458942607424e-07 Ha`
- `d = 0.20`: energy `-2.8616833781331921`, signed error
  `-3.382520953110e-06 Ha`

Exact task:

Create a `tmp/work` developer probe for `n_s = 7` He RHF using the old
nested/QW MWG path first. This is an accuracy diagnostic and paper-target
comparison, not a new acceptance test.

Use:

- He, `Z = 2`
- AHGBS-9 S-only supplement, `lmax = 0`
- `restrictedHF = true`
- `doside = n_s = 7`
- the same external basis source used in the successful `n_s = 5` reproduction:
  `/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/BasisSets`
- `interaction_treatment = :mwg`
- ordinary final-basis solve; no generalized final solve

Try at least:

1. `n_s = 7`, `d = 0.15`, because Fig. 8 reports the smallest positive error
   in the `n_s = 7` family there.
2. If runtime is reasonable, also run `d = 0.10` and `d = 0.20` to confirm the
   local trend around the minimum.

For each point, record:

- constructor and mapping parameters used;
- physical endpoints and any legacy box convention notes;
- gausslet count;
- residual S-GTO count;
- final dimension;
- final overlap identity error;
- H1 lowest orbital energy;
- IDA self-Coulomb diagnostic;
- RHF one-electron energy;
- RHF electron-electron energy;
- RHF total energy;
- error relative to the Fig. 8 table row for that `n_s,d`;
- error relative to the He HF reference
  `-2.861679995612238878775544`;
- iterations and coarse timing.

Decision rules:

- If the `n_s = 7`, `d = 0.15` point is within a few mH of the Fig. 8 table,
  report it and continue no further unless runtime is cheap.
- If it is in the sub-mH regime, run the adjacent `d = 0.10` and `d = 0.20`
  points if practical.
- If the mismatch is unexpectedly large, stop after the first point and report
  whether the issue appears to be basis loading, dimension mismatch,
  mapping/box setup, or RHF/MWG convention.
- Do not try to fix mapping conventions in this pass unless there is a tiny,
  unambiguous bug.

Trust boundary:

Use `tmp/work` probes only. Do not add an acceptance test. Do not add public
driver/export paths. Do not change production code unless an unambiguous tiny
bug blocks the diagnostic. Do not silently replace AHGBS-9 with cc-pVTZ. Do not
introduce raw GTO density-density final operators, full-parent CPB fallback,
direct Cartesian product fallback, ordinary Cartesian IDA fallback, PQS, or a
generalized final-basis solve.

Artifacts:

- `tmp/work/fig8_he_ns7_rhf_probe.jl`
- `tmp/work/fig8_he_ns7_rhf_summary.txt`
- `tmp/work/fig8_he_ns7_rhf.tsv` if multiple points are run

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
- points run;
- final dimensions and residual counts;
- energies and errors vs Fig. 8 and He HF;
- timing;
- whether `n_s = 7` looks accurate enough to become the next comparison target;
- validation run;
- deletion/shrinkage report.
