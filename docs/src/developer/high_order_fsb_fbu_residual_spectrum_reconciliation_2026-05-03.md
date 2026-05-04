# High-order full-shell basis / full-block union residual-spectrum reconciliation

Date: 2026-05-03

This note reconciles the tension between:

- older residual reports that showed a real distorted residual sector, and
- the recent bounded same-backend He+ study in which the **full-shell basis**
  (FSB) and the **full-block union** (FBU) agreed to roundoff.

The result of the present pass is clear:

- on the current physical-coordinate polynomial route, the old distorted
  residual sector is **not** reproduced
- shell-only and full-block-before-add residuals are large, as expected
- but once the current physical shell is added, the full-block-after-add
  residual collapses to numerical zero at every tested side

So, for the current physical route, the later He+ FSB/FBU roundoff agreement is
fully explainable: the final FSB and FBU really are the same subspace to
numerical precision in the bounded cases tested here.

## Why this pass was needed

The older note

- [Residual Gausslet Hybrid Report](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/docs/src/developer/residual_gausslet_hybrid_report_2026-04-27.md)

recorded nonzero distorted residual spectra. For example, in the `11^3`,
`n_s = 5` case, the distorted-shell hybrid reported:

- side `9`: trace `2.061642e-03`, largest eigenvalue `1.245477e-04`
- side `11`: trace `2.356384e-03`, largest eigenvalue `1.696710e-04`

and the homotopy table showed monotone growth of the residual from the identity
limit to the fully distorted case.

That older evidence was strong enough that a direct reconciliation was needed.

## Construction used here

This pass used the **current physical-coordinate polynomial route** throughout.

Parent / mapping family:

- `IdentityMapping()` control
- `AsinhMapping(c = 0.2, s = scale * sqrt(0.4), tail_spacing = 10.0)`

Primary bounded case:

- `count = 11`, `doside = 5`, sides `5,7,9,11`
- scales:
  - `identity`
  - `s/s0 = 1.0, 0.8, 0.6, 0.4, 0.2`

Spot-check case:

- `count = 13`, `doside = 5`, sides `5,7,9,11,13`
- `s/s0 = 1.0`

### FBU construction

For each side, the physical full block was built from:

- `_experimental_high_order_physical_full_block_3d(...)`

The final FBU was the overlap-metric-cleaned union of those physical full-block
coefficient matrices.

### FSB construction

The final FSB was built as:

1. the first physical full block
2. then, for each larger side, the physical shell coefficients from that side
3. each added after projection against the accumulated span and Lowdin cleanup

So this is the current intended physical-`x` shell/full construction, not the
older compatibility/debug route.

## Residual targets that were distinguished

To avoid collapsing different notions of “residual”, the pass measured three
targets at each side:

1. `full_block_before_add`
   - residual of the current side’s physical full block against the accumulated
     FSB span **before** adding the current shell

2. `shell_only_before_add`
   - residual of the current side’s physical shell coefficients against the
     accumulated FSB span **before** adding the current shell

3. `full_block_after_add`
   - residual of the current side’s physical full block against the accumulated
     FSB span **after** adding the current shell

The last one is the most relevant comparison to the older distorted residual
story, because it asks whether the current structured shell leaves any leftover
full-block sector after being added.

For each residual matrix `Y`, the Gram matrix

`G = Y' S_parent Y`

was analyzed, and the following were reported:

- `trace(G)`
- largest eigenvalue
- counts above `1e-8`, `1e-10`, `1e-12`
- target dimension
- accumulated FSB dimension before/after the shell addition

The final He+ metrics from the recent bounded study were also carried along:

- `E_FSB - E_FBU`
- FBU ground-state capture deficiency by FSB
- final FSB/FBU cross-overlap error

## Driver and artifacts

Driver:

- [tmp/work/high_order_fsb_fbu_residual_spectrum_reconciliation.jl](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/tmp/work/high_order_fsb_fbu_residual_spectrum_reconciliation.jl)

Saved artifacts:

- [TSV table](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/tmp/work/high_order_fsb_fbu_residual_spectrum_reconciliation_2026-05-03_173141.tsv)
- [text summary](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/tmp/work/high_order_fsb_fbu_residual_spectrum_reconciliation_2026-05-03_173141.txt)

## Main results

### 1. The current shell does carry the new directions before it is added

For the primary `count=11`, `doside=5` case, at every distorted scale:

- `shell_only_before_add` had trace very close to `98`
- the largest eigenvalue was essentially `1`
- there were `98` eigenvalues above all three thresholds

So the new shell is not redundant before it is added. It carries a full
rank-`98` new sector, just as one would expect.

Likewise:

- `full_block_before_add` had trace around `93` to `98` for outer sides
- and also had `98` eigenvalues above the thresholds

So before the shell is added, the current full block also has a substantial
missing sector relative to the accumulated FSB span.

### 2. After shell addition, the current full block residual disappears

This is the decisive result.

For every tested side and every tested scale:

- `full_block_after_add` had trace around `1e-27` to `1e-26`
- largest eigenvalue around `1e-30` to `1e-27`
- zero eigenvalues above `1e-12`

Examples for `count=11`, `doside=5`, `s/s0 = 1.0`:

- side `7`, `full_block_after_add`
  - trace `1.662e-26`
  - largest eigenvalue `2.824e-27`
  - counts above `1e-8`, `1e-10`, `1e-12`: `0, 0, 0`

- side `9`, `full_block_after_add`
  - trace `1.522e-26`
  - largest eigenvalue `3.192e-27`
  - counts above thresholds: `0, 0, 0`

- side `11`, `full_block_after_add`
  - trace `8.715e-27`
  - largest eigenvalue `2.660e-27`
  - counts above thresholds: `0, 0, 0`

This same qualitative pattern held for:

- identity
- all distorted homotopy points
- the `count=13`, `s/s0=1.0` spot check

### 3. Final FSB/FBU agreement stays at roundoff level

Because the post-add full-block residual vanishes, the final FSB/FBU agreement
from the earlier bounded He+ note is no longer mysterious.

For the `count=11` homotopy:

- max `|E_FSB - E_FBU| = 2.60e-14`
- max FBU-state capture deficiency `= 1.22e-15`
- final cross-overlap errors stayed around `1e-13` to `1e-12`

For the `count=13`, `s/s0=1.0` spot check:

- `|E_FSB - E_FBU| = 9.10e-15`
- capture deficiency `= 1.55e-15`
- final cross-overlap error `= 5.05e-13`

So the one-electron agreement is explained by the basis result itself:

- after each shell is added, the current physical FSB already spans the current
  physical full block
- therefore the final FSB and FBU are the same subspace to numerical precision

## Reconciliation with the older residual report

The older residual report showed a genuine distorted residual sector.

The present pass shows that the **current physical-coordinate polynomial route**
does not reproduce that behavior.

So the reconciliation is:

1. The older report was measuring a real residual sector for the construction
   used at that time.
2. The current physical route, as now implemented, does not leave a measurable
   post-shell full-block residual in the bounded cases tested here.
3. Therefore the current He+ FSB/FBU roundoff agreement is not surprising; it
   follows directly from the current basis algebra.

What remains uncertain is historical, not numerical:

- exactly which change between the older distorted-shell residual study and the
  current physical route removed that post-shell residual sector

But the numerical reconciliation for the current code is now straightforward:

- the old distorted residual sector is **not** reproduced on the current
  physical route
- the current physical shell additions close the current physical full blocks
  exactly to numerical resolution

## Performance note

The pass stayed bounded as intended.

Measured total per-case times were:

- `count=11` homotopy points: about `1.19 s` to `2.42 s`
- `count=13`, `s/s0=1.0` spot check: about `2.81 s`

The code reused:

- basis setup
- axis data
- parent overlap / one-body package
- physical full blocks per case/scale

so the pass remained cheap enough to run as a correctness diagnostic.

## Bottom line

For the current physical-coordinate high-order route:

- the shell-only and full-block-before-add residuals are real and large
- but the full-block-after-add residual is numerically zero
- therefore the final FSB and FBU agree to roundoff in the bounded cases tested

So the recent He+ FSB/FBU roundoff agreement is still explainable, and it is
explained by the current residual-spectrum result itself rather than by a
transfer-metric artifact.
