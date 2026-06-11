Purpose:

Run the first side13 He + GTO final-basis RHF diagnostic using the
now-materialized final-basis density-density matrix.

Why now:

Pass 002 materialized final-basis density-density:

```text
V_final =
[ V_gg   V_gR
  V_Rg   V_RR ]
```

for side13 He + H cc-pVTZ `lmax = 0`. The route is now ready for a diagnostic
RHF solve in the final orthonormal basis. This should be a probe first, not a
new acceptance test.

Current state:

- gausslet-only side13 He RHF total is about `-2.8364979997009137`.
- gausslet-only best fixed-`ns = 5` exploratory side17 point is about
  `-2.858531351214`, but side17 is not a routine gate.
- side13 He + GTO final-basis one-electron machinery exists.
- side13 He + GTO final-basis density-density now materializes.
- raw GTO density-density is not accepted as final operator data.

Exact task:

Create and run a `tmp/work` developer probe for side13 He + GTO final-basis RHF.

The probe should:

1. Build the side13 decomposed WL He fixture:
   - `Z = 2`
   - `d = 0.1`
   - `s = 1.0`
   - side count 13
   - shellification-backed decomposed WL inventory

2. Build the He cc-pVTZ `lmax = 0` GTO supplement.

3. Build the existing combined one-electron matrices and final-basis projection.

4. Build the final-basis density-density matrix using:
   - existing decomposed WL `V_gg`;
   - residual MWG `V_gR` / `V_RR`;
   - no raw GTO density-density final operator.

5. Run restricted closed-shell RHF in the final orthonormal basis using:
   - final one-electron Hamiltonian;
   - final density-density matrix;
   - ordinary symmetric final-basis solve.

6. Compare against:
   - side13 gausslet-only RHF total;
   - He HF reference near `-2.861679995612234`;
   - exact He reference only as context, not as the HF target.

Trust boundary:

Do not add an acceptance test in this pass.

Do not add public driver/export/artifact paths.

Do not use generalized-overlap final solve as the accepted path.

Do not use raw GTO/GTO density-density as final operator data.

Do not use full-parent CPB fallback, flat `_qwrg_diatomic_*` fallback, old
fixed-block matrix authority, ordinary Cartesian IDA fallback, PQS, or side17
long-running gates.

Implementation preference:

Use a `tmp/work` probe first. If a tiny private helper is needed to avoid
duplicating an existing final-basis solve utility, keep it scoped and report why
the probe could not reuse existing code. Do not add production code unless the
missing helper is clearly part of the live final-basis RHF path.

Record in the artifact:

- fixture parameters and endpoints;
- retained gausslet dimension;
- raw supplement count;
- retained residual count;
- final dimension;
- final overlap identity error;
- final one-electron Hamiltonian shape;
- final density-density shape;
- `V_gR` and `V_RR` shapes;
- RHF convergence status;
- iterations;
- final one-electron energy;
- final electron-electron energy;
- final RHF total;
- improvement over side13 gausslet-only RHF;
- error from He HF reference;
- timing split:
  - inventory;
  - one-electron combined/final-basis build;
  - density-density build;
  - RHF solve;
  - total elapsed;
- anti-fallback flags.

Artifacts:

- `tmp/work/side13_he_gto_final_basis_rhf_probe.jl`
- `tmp/work/side13_he_gto_final_basis_rhf_summary.txt`

Validation:

- run the probe;
- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`.

Decision rules:

- If RHF improves over side13 gausslet-only and remains above the HF reference,
  report it as a plausible diagnostic.
- If RHF drops below the HF reference, stop and audit the final-basis
  density-density or RHF convention before making any acceptance claim.
- If final overlap is not identity within tolerance, stop; do not use a
  generalized final solve.
- If timing is very slow, report the phase split before optimizing.
- If any path falls back to raw GTO density-density, full-parent CPB, direct
  Cartesian product, or ordinary IDA, stop and report.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- artifact paths;
- final RHF energy and comparison to side13 gausslet-only / HF reference;
- whether the result is suitable for a future acceptance fixture;
- validation run;
- deletion/shrinkage report.
