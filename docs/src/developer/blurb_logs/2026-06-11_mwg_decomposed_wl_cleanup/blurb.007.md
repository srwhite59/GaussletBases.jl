Purpose:

Use the located White-Lindsey Fig. 8 He RHF data to define the next comparison
target. The immediate task is to audit whether the repo can reproduce the
paper-like `n_s = 5`, `d = 0.3`, AHGBS-9 S-only setup, not to add a new
acceptance test.

Why now:

The old nested fixed-block QW/MWG oracle matches the current decomposed side13
He + GTO diagnostic to roundoff, so the decomposed route is not the local
culprit. The missing piece is the external paper target. The archive now has a
curated Fig. 8 extract and provenance note, and those show that the previous
`d = 0.3` proxy was not a real Fig. 8 reproduction: it used the repo's current
cc-pVTZ S-only supplement path and produced final dimension `422`, while the
paper setup used AHGBS-9 S-only and labels the `n_s = 5` curve near `447`
basis functions.

Archive authority:

- data table:
  `/Users/srw/Library/CloudStorage/Dropbox/chatarchive/references/papers/canonical/White and Lindsey - 2023 - Nested Gausslet Basis Sets_fig8_heuhf_data.tsv`
- provenance note:
  `/Users/srw/Library/CloudStorage/Dropbox/chatarchive/references/papers/canonical/White and Lindsey - 2023 - Nested Gausslet Basis Sets_fig8_heuhf_provenance.md`

Important Fig. 8 setup from the provenance note:

- He RHF
- `Z = 2`, one atom
- `basisaddname = "AHGBS-9"`
- `lmaxadd = 0`
- `restrictedHF = true`
- `doside = n_s`
- `corespacing = d`
- `dwidth = 10.0`
- `gscalefac = sqrt(2.0)`
- reference energy:
  `-2.861679995612238878775544 Ha`
- manuscript caption says box size `R_b = 7 bohr`, while legacy logs mention
  `xmin = -8.0`, `xmax = 8.0`, and `basradius = 8.0`; record this convention
  mismatch instead of hiding it.

Key table targets:

- `n_s = 5`, `d = 0.3`:
  - energy `-2.8615437846242582`
  - signed error `+1.362109879808e-04 Ha`
  - about `0.136 mHa`
- `n_s = 5`, `d = 0.1`:
  - energy `-2.8610408147225441`
  - signed error `+6.391808896948e-04 Ha`
  - about `0.639 mHa`
- best plotted point:
  - `n_s = 9`, `d = 0.1`
  - signed error about `+1.849157982292e-08 Ha`

Exact task:

Create a `tmp/work` audit probe/report that answers whether the repo can build
a Fig. 8-like He RHF fixture from existing constructors.

Do this in order:

1. Confirm whether an AHGBS-9 S-only He supplement can be loaded through the
   existing basis plumbing, likely via `legacy_atomic_gaussian_supplement(...)`
   or the current named-basis loader.

2. If AHGBS-9 is available, build the closest old nested/QW MWG Fig. 8
   reproduction for:
   - `n_s = 5`
   - `d = 0.3`
   - `lmax = 0`
   - AHGBS-9 S-only supplement
   - `restrictedHF = true`

3. Resolve and report the mapping/box convention used by the repo constructor:
   - how `corespacing = d`, `dwidth = 10.0`, and `gscalefac = sqrt(2.0)` map
     to current `MappedUniformBasisSpec` / `AsinhMapping` parameters;
   - whether the repo fixture corresponds to caption `R_b = 7` or legacy
     `xmin/xmax = +/-8` behavior;
   - the parent side count and final basis dimension.

4. Compare the produced energy to the Fig. 8 table target
   `-2.8615437846242582` for `n_s = 5`, `d = 0.3`.

5. If AHGBS-9 is not available, do not substitute cc-pVTZ silently. Report the
   exact missing loader/source blocker and keep the prior cc-pVTZ side13 result
   labeled as a diagnostic, not a Fig. 8 proxy.

6. If the final dimension is not near the Fig. 8 label `447`, report why:
   - supplement count differs;
   - shell/fixed-block construction differs;
   - residual keep policy differs;
   - mapping/box convention differs;
   - or source basis is unavailable.

Trust boundary:

Use `tmp/work` probes first. Do not add an acceptance test. Do not change
production code unless a tiny, clear loader bug blocks an otherwise available
AHGBS-9 basis. Do not introduce raw GTO density-density final operators,
full-parent CPB fallback, direct Cartesian product fallback, ordinary Cartesian
IDA fallback, PQS, exports, or a generalized final-basis solve.

Do not use cc-pVTZ as a stand-in for AHGBS-9 when claiming Fig. 8 agreement.
cc-pVTZ may be mentioned only as the previous diagnostic supplement.

Artifacts:

- `tmp/work/fig8_he_rhf_target_reproduction_probe.jl`
- `tmp/work/fig8_he_rhf_target_reproduction_summary.txt`
- optional TSV if multiple attempted fixtures are compared

Validation:

- run the probe;
- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`;
- if production code changes are made for a clear loader bug, run the most
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
- whether AHGBS-9 S-only was available;
- exact constructor and mapping parameters used;
- final dimension and residual count;
- energy and error relative to the Fig. 8 `n_s = 5`, `d = 0.3` target;
- whether the repo can currently reproduce the Fig. 8 point;
- if not, exact blocker;
- validation run;
- deletion/shrinkage report.
