Purpose:

Move from atomic He to the first diatomic H2 HF diagnostic. This should be a
developer-only probe at the restricted closed-shell HF level, using the
existing old nested/QW H2 route as the first reference path.

Why now:

The `n_s = 7` He RHF probe is accurate enough to close the current atomic
accuracy check:

- `n_s = 7`, `d = 0.10` old nested/QW MWG RHF total:
  `-2.861673961528321`
- error vs Fig. 8 row: `+1.716095156645281e-6 Ha`
- error vs He HF reference: `+6.034083917860755e-6 Ha`

That is good enough to move to a new physics target. The next useful target is
diatomic H2 at the HF level, because it exercises two-center geometry,
bond-aligned nested construction, molecular nuclear terms, residual ownership
on two centers, and the restricted closed-shell density-density HF machinery in
a setting more demanding than He.

Existing H2 reference line:

The repo already has an old-standard H2 chemistry reproduction note:

- `docs/src/developer/high_order_endcap_panel_h2_chemistry_reproduction_2026-05-16.md`
- `docs/src/developer/high_order_mainline_import_readiness_2026-05-15.md`

The key first fixture there is:

- H2 at `R = 4.0` bohr
- `core_spacing = 0.5`
- `xmax_parallel = 6.0`
- `xmax_transverse = 4.0`
- H/cc-pVTZ molecular Gaussian supplement with `lmax = 1` (S/P)
- `nside = 5`
- finite IDA/QW density-density Hamiltonian
- restricted closed-shell HF

Known reference HF totals from the reproduction note:

- default complete rectangular route:
  - final dimension `481`
  - residual count `18`
  - HF total `-0.910938264352`
- endcap/panel `q = 4`, `L = 4` route:
  - final dimension `461`
  - residual count `18`
  - HF total `-0.910977315003`

Exact task:

Create a `tmp/work` developer probe for H2 restricted closed-shell HF at
`R = 4.0` using the old nested/QW route first. Do not add an acceptance test.

Use existing constructors where possible:

- `bond_aligned_homonuclear_qw_basis(...)`
- `bond_aligned_diatomic_nested_fixed_block(...)`
- `legacy_bond_aligned_diatomic_gaussian_supplement("H", "cc-pVTZ", nuclei; lmax = 1)`
- `ordinary_cartesian_qiu_white_operators(...)`

Start with the old-standard/default complete rectangular route. If it is cheap
and clean, also run the endcap/panel `q = 4`, `L = 4` route using the existing
opt-in shared-shell policy. Do not implement new high-order machinery.

Interaction route:

The historical H2 chemistry rows used `interaction_treatment = :ggt_nearest`.
Use that first so the diagnostic compares to the documented old-standard HF
totals. If `:mwg` is cheap to run as an additional comparison, it may be
reported as a secondary diagnostic, but do not let that expand the task.

Record for each route:

- route label and constructor options;
- bond length and nuclear coordinates;
- fixed block size;
- final dimension;
- residual count;
- final overlap identity error;
- residual owner set;
- interaction treatment;
- one-electron and electron-electron HF contributions if available;
- restricted HF total;
- difference from the documented HF total;
- HF convergence iterations;
- coarse timing for fixed-block construction, operator construction, HF solve,
  and total.

Decision rules:

- If the default complete rectangular route reproduces the documented HF total
  to a small numerical tolerance, report success and stop unless the endcap
  route is cheap.
- If the default route does not reproduce the documented HF total, stop and
  localize whether the difference is fixture construction, supplement, QW
  operator interaction treatment, or HF solver convention.
- If the endcap route runs, compare it with the documented `-0.910977315003`
  HF total but keep it secondary.
- Do not run ED in this pass. This is HF only.

Trust boundary:

Use `tmp/work` probes only. Do not add an acceptance test. Do not add public
driver/export paths. Do not change production code unless an unambiguous tiny
bug blocks the diagnostic. Do not introduce decomposed-route changes, raw GTO
density-density final operators, full-parent CPB fallback, direct Cartesian
product fallback, ordinary Cartesian IDA fallback, PQS, or generalized
final-basis solves.

Artifacts:

- `tmp/work/h2_r4_qw_hf_probe.jl`
- `tmp/work/h2_r4_qw_hf_summary.txt`
- optional TSV if multiple route variants are run

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
- routes run;
- HF totals and differences from documented references;
- dimensions/residual counts;
- timing;
- whether H2 R = 4 HF is a clean next physics target;
- validation run;
- deletion/shrinkage report.
