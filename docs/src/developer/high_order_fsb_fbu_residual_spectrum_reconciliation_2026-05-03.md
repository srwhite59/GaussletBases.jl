# High-order full-shell basis / full-block union residual-spectrum reconciliation

Date: 2026-05-03

This note supersedes the overly strong conclusion in commit `e768d3c`.

That earlier pass correctly showed that the current **full-shell basis**
(FSB) exactly spans the current reduced transformed-block target. But it did
**not** test completeness against the true distorted local cube. The present
audit separates those targets explicitly and shows that the older distorted
residual-sector result is **not** contradicted.

## Summary

Three distinct target notions must be kept separate:

1. `reduced_physical_transformed_block`
   - the current physical-coordinate transformed block built from the retained
     `doside` one-dimensional span
   - for `doside = 5`, side `11`, this has `125` raw columns

2. `debug_u_transformed_block`
   - the older compatibility/debug transformed block
   - in the bounded cases tested here, it has the same reduced target size and
     behaves the same as the physical route when judged against the true local
     cube

3. `true_local_distorted_cube`
   - the actual centered local cube in the parent basis
   - for side `11`, this has `1331` raw columns

Main conclusion:

- `e768d3c` proved exactness only for the reduced transformed-block target
- the new rank/null audit shows that the **full-block union** (FBU) and FSB are
  exact for that reduced target, but still leave a large residual sector when
  tested against the true local distorted cube
- so the older distorted residual-sector story survives

## Why this audit was needed

The older note

- [Residual Gausslet Hybrid Report](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/docs/src/developer/residual_gausslet_hybrid_report_2026-04-27.md)

recorded nonzero distorted residual spectra. For example, in the `11^3`,
`n_s = 5` case, it reported:

- side `9`: trace `2.061642e-03`, largest eigenvalue `1.245477e-04`
- side `11`: trace `2.356384e-03`, largest eigenvalue `1.696710e-04`

The recent He+ FSB-versus-FBU sweep instead gave roundoff agreement. The risk
was that these two results were answering different target questions.

That is exactly what the present audit confirms.

## Construction audited here

Parent / mapping family:

- `IdentityMapping()` control
- `AsinhMapping(c = 0.2, s = scale * sqrt(0.4), tail_spacing = 10.0)`

Primary bounded case:

- `count = 11`, `doside = 5`, sides `5, 7, 9, 11`
- scales:
  - `identity`
  - `s/s0 = 1.0`

Spot check:

- `count = 13`, `doside = 5`, sides `5, 7, 9, 11, 13`
- `s/s0 = 1.0`

The current intended physical route was used for the working basis:

- physical full blocks from
  `_experimental_high_order_physical_full_block_3d(...)`
- physical shells accumulated into the FSB by metric projection and Lowdin
  cleanup

The debug route was retained only as a control:

- `_experimental_high_order_tensor_shell_3d(...)`

## What was measured

For each route, side, and target kind, the audit reported:

- raw target column count
- target Gram rank and spectrum before cleanup
- FBU union raw column count
- FBU union metric rank before cleanup
- FBU cleaned dimension
- FSB cleaned dimension
- residual trace / largest eigenvalue / counts above `1e-8`, `1e-10`,
  `1e-12` after projecting the target against the accumulated FSB span

This explicitly checks both:

- whether the target itself is reduced or full-rank
- whether FBU cleanup is discarding many redundant or null directions

## Driver and artifacts

Driver:

- [tmp/work/high_order_fsb_fbu_residual_spectrum_reconciliation.jl](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/tmp/work/high_order_fsb_fbu_residual_spectrum_reconciliation.jl)

Artifacts from the audit pass:

- [TSV table](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/tmp/work/high_order_fsb_fbu_residual_spectrum_reconciliation_2026-05-03_175030.tsv)
- [text summary](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/tmp/work/high_order_fsb_fbu_residual_spectrum_reconciliation_2026-05-03_175030.txt)

The older reduced-target-only artifacts from `e768d3c` remain useful as a
control:

- [TSV table](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/tmp/work/high_order_fsb_fbu_residual_spectrum_reconciliation_2026-05-03_173141.tsv)
- [text summary](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/tmp/work/high_order_fsb_fbu_residual_spectrum_reconciliation_2026-05-03_173141.txt)

## Results

### 1. The “physical full block” at side 11 is reduced, not a true cube

For `count = 11`, `doside = 5`, side `11`:

- `reduced_physical_transformed_block`: `125` raw columns
- `debug_u_transformed_block`: `125` raw columns
- `true_local_distorted_cube`: `1331` raw columns

So the current physical “full block” is the full tensor product of the
retained transformed one-dimensional block, not the true `11^3` local cube.

This is the key reason `e768d3c` was not yet a valid reconciliation.

### 2. FBU cleanup really is dropping many redundant directions

For `count = 11`, sides `5:2:11`:

- FBU raw union column count: `500`
- FBU union metric rank before cleanup: `419`
- FBU cleaned dimension: `419`
- FSB cleaned dimension: `419`

For `count = 13`, sides `5:2:13`:

- FBU raw union column count: `625`
- FBU union metric rank before cleanup: `517`
- FBU cleaned dimension: `517`
- FSB cleaned dimension: `517`

So the raw FBU union carries many redundant directions before Lowdin cleanup.
That was not checked in `e768d3c`, and it matters.

### 3. Reduced transformed-block targets are still exact

For both routes, both identity and distorted cases:

- the residual against `reduced_physical_transformed_block` is numerical zero
- the residual against `debug_u_transformed_block` is numerical zero

Example, `count = 11`, `doside = 5`, side `11`, `s/s0 = 1.0`:

- reduced physical target:
  - trace `8.715e-27`
  - largest eigenvalue `2.660e-27`
  - counts above `1e-8`, `1e-10`, `1e-12`: `0, 0, 0`

This confirms the narrow statement from `e768d3c`:

- the current FSB exactly spans the current reduced transformed-block target

### 4. The true local distorted cube leaves a large residual sector

This is the decisive audit result.

For `count = 11`, `doside = 5`, side `11`, both identity and `s/s0 = 1.0`,
and for both `physical_x` and `debug_u` routes:

- true local cube raw columns: `1331`
- true local cube metric rank: `1331`
- residual trace after projection against accumulated FSB: `9.120e+02`
- largest residual eigenvalue: `1.000e+00`
- counts above `1e-8`, `1e-10`, `1e-12`: `912, 912, 912`

For the `count = 13`, side `13` spot check:

- true local cube raw columns: `2197`
- true local cube metric rank: `2197`
- residual trace: `1.680e+03`
- largest residual eigenvalue: `1.000e+00`
- counts above thresholds: `1680, 1680, 1680`

So the old “nonzero residual sector under distortion” story is not contradicted
by the current route. Once the target is the true local cube rather than the
reduced transformed block, a large residual sector remains.

### 5. Physical-x and debug-U do not separate in this audit

In the bounded cases tested here:

- `physical_x` and `debug_u` are both exact on the reduced transformed-block
  target
- and both show the same large residual on the true local cube target

So the current discrepancy is not “physical route versus debug route.” It is
“reduced transformed-block target versus true local-cube target.”

## Reconciliation

The tension is now resolved as follows:

1. The older residual report was about completeness relative to a fuller local
   cube notion, and it found a real distorted residual sector.
2. The recent same-backend He+ FSB-versus-FBU sweep, and the earlier
   `e768d3c` reconciliation, were effectively testing only the current reduced
   transformed-block target.
3. On that reduced target, FSB and FBU are exact to roundoff.
4. On the true local distorted cube target, a large residual sector remains.

But this audit does **not** reproduce the older distortion-defect homotopy.
Here, the true-local-cube residual is already large at identity. So the old
story of a small zero-at-identity residual that then grew with distortion is
still historically unresolved.

So `e768d3c` should be read narrowly and is superseded by this audit:

- it established exactness only for the current reduced transformed-block
  target
- it did **not** establish true distorted local-cube completeness
- and it did not determine which historical target/construction produced the
  older zero-at-identity homotopy

## Consequence for the recent He+ agreement

The recent He+ FSB/FBU roundoff agreement is still explainable, but only in the
reduced-target sense:

- the current FSB and FBU agree because they span the same reduced transformed
  block union
- that says nothing by itself about completeness against the true local
  distorted cube

So the He+ agreement is not a contradiction of the older residual report. It is
answering a narrower question.

## Performance note

This was a bounded scratch/docs audit with reused parent setup per case. The
timings were acceptable for `count = 11` and the `count = 13` spot check, so no
broader optimization work was needed here.

## Bottom line

The corrected conclusion is:

- the current FSB exactly reproduces the current reduced transformed-block
  target
- the side `11` physical full block in that reduced sense has `125` columns,
  not `1331`
- FBU cleanup drops many redundant directions before final cleanup
- when the target is the true local distorted cube, a large nonzero residual
  sector remains
- therefore the older distorted residual-sector result is **not** contradicted
  by the recent reduced-target FSB/FBU studies
- but the old small-at-identity homotopy is still not reconstructed here, so
  the remaining open question is which historical target/construction produced
  that behavior
