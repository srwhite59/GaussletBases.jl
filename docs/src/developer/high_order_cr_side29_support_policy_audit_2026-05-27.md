# High-Order Cr Side-29 Support Policy Audit

## Status

This is a repo-side occupied-capture and support-policy audit only. It does
not claim Hamiltonian, two-electron, same-density, or energy validation.

The source evidence is CR2's existing side-29 PySCF occupied-orbital handoff:

- report:
  `~/Dropbox/codexhome/work/cr2/reports/cr_ccecp_side29_parent_occupied_handoff_2026-05-25.md`
- artifact directory:
  `~/Dropbox/codexhome/work/cr2/runs/2026-05-25_cr_ccecp_side29_parent_occupied_handoff`
- repo-side cumulative support scan:
  `tmp/work/high_order_cr_side29_cumulative_capture_report.txt`
- repo-side failed cropped high-order gate:
  `tmp/work/high_order_cr_side29_h1_capture_gate.txt`

The side-29 raw parent coefficient artifacts are present in the CR2 handoff
directory. This note does not rerun PySCF or recompute those coefficients.

## Fixture Identity

The audited handoff is:

- family: `G10`
- nuclear charge label: `Z = 24`
- mapping spacing: `d = 0.0125`
- `dZ = 0.3`
- target radius: `10.0`
- parent side: `29`
- parent dimension: `24389`
- downstream backend label: `:pgdg_localized_experimental`

The existing repo-side scan verified parent-grid ordering and reported
`numerical_reference_used=false`. The older side-27 `Z_eff = 14` fixture is a
different target and is not part of this audit unless a matching CR2 occupied
handoff is provided.

## Support Scan

The cumulative centered-box scan shows where the PySCF occupied norm lives in
the side-29 parent representation:

| support side | support count | alpha capture | beta capture | worst alpha orbital | worst beta orbital |
|---:|---:|---:|---:|---|---|
| 19 | 6859 | `0.8573543538` | `0.9948604749` | `alpha_mo9_col10_target10`: `0.1671875456` | `beta_mo3_col4_target14`: `0.9936608993` |
| 21 | 9261 | `0.9515023713` | `0.9999696891` | `alpha_mo9_col10_target10`: `0.5951237606` | `beta_mo2_col3_target13`: `0.9999608623` |
| 23 | 12167 | `0.9917274664` | `0.9999999218` | `alpha_mo9_col10_target10`: `0.9212791136` | `beta_mo0_col1_target11`: `0.9999998932` |
| 25 | 15625 | `0.9995071782` | `0.9999999881` | `alpha_mo9_col10_target10`: `0.9951302270` | `beta_mo0_col1_target11`: `0.9999999762` |
| 27 | 19683 | `0.9999916374` | `0.9999999977` | `alpha_mo9_col10_target10`: `0.9999164369` | `beta_mo0_col1_target11`: `0.9999999938` |
| 29 | 24389 | `1.0000000000` | `1.0000000000` | all captured in the full parent | all captured in the full parent |

The same report identified side `25` as the first centered box above `0.99`
alpha capture and side `27` as the first above `0.999` alpha capture. Beta
capture reaches the same thresholds earlier.

## Interpretation

The side-29 parent representation was adequate for the CR2 occupied-orbital
handoff. The full parent captures the reported raw occupied columns as:

- alpha: `9.9994606594` captured out of `9.9994606594`
- beta: `3.9996329181` captured out of `3.9996329181`

The catastrophic cropped high-order captures were therefore a support-policy
failure, not evidence that the side-29 parent box itself was inadequate. The
failed cropped rows were:

| policy | retained dim | covered support | alpha capture | beta capture |
|---|---:|---:|---:|---:|
| `baseline_centered_complete_doside4_sides4_6_8` | 176 | 512 | `0.0000783340` | `0.0000369714` |
| `terminal_core4_side6_full_side8_z_side_panels` | 144 | reported terminal/cropped support | `0.0000281008` | `0.0000144174` |

Those policies stopped near the center and omitted the outer support required
by the alpha occupied space, especially `alpha_mo9_col10_target10`.

The current default high-order doside policy is aligned with the support
evidence: omitted `sides` means full-parent coverage, so a side-29,
`doside = 5` request uses the ladder `[5, 7, ..., 29]`. Explicit cropped
`sides` remain a diagnostic or research request and should not be interpreted
as the default Cr support policy.

## Remaining Risk

The existing CR2 evidence supports the new full-parent default support policy.
It does not validate a high-order Hamiltonian or energy route.

The remaining high-order risk is now retained span and rank within the
full-domain contracted parent, plus any future diatomic recipe design. Support
coverage should still be audited for explicit cropped experiments, but it is
not the blocker for the side-29 parent representation when the full-parent
default ladder is used.
