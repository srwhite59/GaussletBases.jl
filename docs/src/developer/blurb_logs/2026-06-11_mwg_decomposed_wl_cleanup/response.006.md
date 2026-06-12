Implemented the audit as a developer-only oracle probe.

Artifacts:

- `tmp/work/side13_he_gto_old_qw_mwg_oracle_probe.jl`
- `tmp/work/side13_he_gto_old_qw_mwg_oracle_summary.txt`
- `tmp/work/side13_he_gto_old_qw_mwg_oracle.tsv`

Matched side13 old nested fixed-block QW/MWG oracle:

- fixture: He, `d = 0.1`, `s = 1.0`, side count `13`
- endpoints: `(-8.565228460168399, 8.565228460168399)`
- final dimension: `520`
- gausslet/residual counts: `517 / 3`
- final overlap identity error: `5.318190332559425e-12`
- H1 lowest energy: `-1.9995342256997053`
- IDA self-Coulomb J: `1.2441722148113061`
- RHF one-electron energy: `-3.8890817574267578`
- RHF electron-electron energy: `1.0246874795868057`
- RHF total: `-2.864394277839952`
- error from He HF reference: `-0.002714282227717746`
- RHF iterations: `19`

Comparison to current decomposed side13 diagnostic from
`tmp/work/side13_he_gto_final_basis_rhf_summary.txt`:

- decomposed final dimension: `520`
- decomposed retained residual count: `3`
- decomposed RHF one-electron energy: `-3.889081757426776`
- decomposed RHF electron-electron energy: `1.0246874795868013`
- decomposed RHF total: `-2.864394277839975`
- old-minus-decomposed total: `2.3092638912203256e-14`
- old-minus-decomposed one-electron: `1.8207657603852567e-14`
- old-minus-decomposed electron-electron: `4.440892098500626e-15`

Decision:

- The matched old nested fixed-block QW/MWG route agrees with the current decomposed side13 He + GTO final-basis diagnostic to roundoff.
- The below-HF result is therefore not localized by this oracle to a new decomposed-route convention mismatch.
- This supports treating the current side13 He + GTO result as a MWG/IDA approximation behavior or broader model/fixture issue, not as an immediate route-specific bug.

Paper-like proxy:

- I did not find a repo-local Fig. 8 numeric table to validate against.
- I ran only a clearly labeled constructible old-QW proxy: `d = 0.3`, `s = 1.0`, side count `11`.
- It produced final dimension `422`, not the cited paper-like `~447` functions, so it is not a paper-value validation.
- Proxy RHF total: `-2.8633283682267505`
- Proxy error from He HF reference: `-0.0016483726145164468`

Old surfaces inspected/used:

- `one_center_atomic_full_parent_fixed_block(...)`
- `ordinary_cartesian_qiu_white_operators(fixed_block, supplement; interaction_treatment = :mwg)`
- `ordinary_cartesian_1s2_check(...)`
- the existing ordinary QW MWG tests around legacy He S residual interactions and pair-density normalization.

Validation:

- Ran `julia --project=. tmp/work/side13_he_gto_old_qw_mwg_oracle_probe.jl`
- Ran `julia --project=. -e 'using GaussletBases; println("load ok")'`
- Ran `git diff --check`

Deletion/shrinkage report:

- No production code, tests, docs, metadata, or compatibility paths became obsolete from this audit-only pass.
- Nothing was deleted or simplified because the task was explicitly an oracle comparison and the old nested fixed-block QW/MWG route remains useful as a reference comparator.
- No new tracked test was added. The probe is developer-only under `tmp/work` and does not replace or expand routine coverage.
- Remaining stale/duplicate surface to retire next: none identified by this pass. The actionable follow-up is interpretive: decide whether to keep He + GTO side13 diagnostic-only or next audit the MWG/IDA approximation itself against external/paper data if an authoritative table is located.

-- repo-doer@macmini
