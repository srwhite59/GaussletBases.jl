Purpose:

Audit the remaining mismatch between the repo's 447-function AHGBS-9
old-QW/MWG reproduction and the White-Lindsey Fig. 8 `n_s = 5`, `d = 0.3`
legacy result.

Why now:

The previous pass solved the basis-count problem: AHGBS-9 S-only loads from the
external GaussletModules basis file and the old nested/QW MWG fixture has final
dimension `447`, matching the Fig. 8/`He.5.3` count. The energy still does not
match. The repo reproduction gives `-2.862102144533723`, while the Fig. 8 /
`He.5.3` target is `-2.861543784624258`, a difference of about
`-0.558 mHa`. This must stay an audit, not an acceptance test.

Current known facts:

- archive Fig. 8 target table:
  `/Users/srw/Library/CloudStorage/Dropbox/chatarchive/references/papers/canonical/White and Lindsey - 2023 - Nested Gausslet Basis Sets_fig8_heuhf_data.tsv`
- provenance note:
  `/Users/srw/Library/CloudStorage/Dropbox/chatarchive/references/papers/canonical/White and Lindsey - 2023 - Nested Gausslet Basis Sets_fig8_heuhf_provenance.md`
- detailed legacy log:
  `/Users/srw/Library/CloudStorage/Dropbox/Nestedruns/He/He.5.3`
- target:
  - `n_s = 5`
  - `d = 0.3`
  - `AHGBS-9`
  - `lmaxadd = 0`
  - `restrictedHF = true`
  - final dimension `447`
  - `uhfen = -2.861543784624258`
- repo reproduction:
  - final dimension `447`
  - gausslet/residual counts `419 / 28`
  - H1 lowest `-1.9999998632985623`
  - IDA self-Coulomb `1.2496940228276845`
  - RHF total `-2.862102144533723`

Legacy log facts to inspect:

- `(xmin, xmax, xbasradius, ymax) = (-8.0, 8.0, 8.0, 8.0)`
- `(doside, dwidth, jflat, wi, gscalefac, polylim, lmaxadd, restrictedHF) =
  (5, 10.0, 39//2, 6.0, 1.4142135623730951, 9, 0, true)`
- `doInvsqrt = true`
- `(rangeg, nlet) = (-8:8, 17)`
- `norm(O1 - eye(Ntot)) = 4.847810690997331e-11`
- `Ntot = 447`
- energy iterations:
  - iteration 1: `-2.861533876725111`
  - iteration 2: `-2.8615437773499623`
  - iteration 3: `-2.861543784624258`
  - iteration 4: `-2.86154378464922`
- `uhfen = -2.861543784624258`

Exact task:

Create a `tmp/work` audit probe/report that localizes the remaining mismatch.
Do not add an acceptance test.

Audit in this order:

1. Mapping/box convention:
   - reconstruct, as closely as possible, the legacy `He.5.3` one-dimensional
     mapping from the log/source data: `xmin/xmax`, `rangeg`, `nlet`,
     `dwidth`, `jflat`, `wi`, `gscalefac`, `doInvsqrt`;
   - compare the 1D coordinate grid/backbone coordinates to the repo
     `white_lindsey_atomic_mapping(Z = 2, d = 0.3, tail_spacing = 10.0)`
     constructor;
   - record coordinate endpoint differences and whether any existing repo
     constructor can match the legacy coordinates without new framework work.

2. Operator/component comparison:
   - compare the repo reproduction's one-body lowest orbital and
     `ordinary_cartesian_1s2_check` self-Coulomb with the legacy log values if
     available;
   - isolate whether the mismatch is mostly in one-body, density-density
     self-Coulomb, or the RHF update/energy formula;
   - if useful, compare the first iteration energy using the initial one-body
     orbital against the legacy iteration-1 energy.

3. RHF convention:
   - inspect the legacy log and any available legacy source for how
     `restrictedHF = true` computes `uhfen`;
   - compare with the compact probe helper's closed-shell formula
     `2 tr(rho h) + direct - exchange`;
   - determine whether the repo helper's twenty-iteration convergence vs the
     legacy four-iteration report is a harmless damping/update difference or an
     energy-convention mismatch.

4. Basis/source handling:
   - confirm again that AHGBS-9 S-only source is the external
     GaussletModules `BasisSets` block and not the vendored repo basis file;
   - do not vendor or copy the basis file in this pass;
   - report whether a future accepted Fig. 8 test would need a curated small
     basis extract or an explicit external-data requirement.

Decision rules:

- If the mismatch is explainable by mapping/box convention and an existing repo
  constructor can match the legacy coordinate grid, run that matched probe and
  report the new energy.
- If the mismatch is explainable by RHF energy/update convention, report the
  exact convention and whether the repo helper should be corrected or merely
  kept probe-local.
- If neither is resolved, report the smallest next source artifact needed
  before further coding.
- Do not change production code unless a narrow bug is unambiguous.

Trust boundary:

Use `tmp/work` probes first. Do not add an acceptance test. Do not add public
driver/export paths. Do not introduce raw GTO density-density final operators,
full-parent CPB fallback, direct Cartesian product fallback, ordinary Cartesian
IDA fallback, PQS, or generalized final-basis solve. Do not silently replace
AHGBS-9 with cc-pVTZ.

Artifacts:

- `tmp/work/fig8_he_rhf_legacy_convention_audit.jl`
- `tmp/work/fig8_he_rhf_legacy_convention_audit_summary.txt`
- optional TSV if multiple mapping/RHF variants are compared

Validation:

- run the audit probe;
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
- mapping/box comparison result;
- one-body vs density-density vs RHF-convention attribution;
- whether an existing repo constructor can reproduce the legacy coordinates;
- whether the `0.558 mHa` mismatch is resolved, reduced, or still open;
- validation run;
- deletion/shrinkage report.
