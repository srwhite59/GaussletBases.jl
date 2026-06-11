What changed:

- Added an audit artifact under `tmp/work`:
  - `tmp/work/side13_he_gto_final_basis_rhf_convention_audit.txt`
- No production code changed.
- No tests were added.
- No acceptance claim was made.

Audit result:

- No clear implementation bug was found in this focused pass.
- The probe-local RHF formula matches the gausslet-only He RHF convention:
  - one doubly occupied spatial orbital;
  - density trace near `1`;
  - electron count near `2`;
  - Fock form `h + 2J - K`;
  - energy form `2 tr(rho h) + direct - exchange`;
  - scalar one-orbital convention `E = 2h + V`.
- The final-basis one-electron path is the intended ordinary final-basis path:
  - final overlap identity error `2.4938717757549966e-11`;
  - solve kind `ordinary_symmetric`;
  - no generalized final-basis solve;
  - no final self-overlap used as working metric data.
- The final-basis density-density path is using the intended pieces:
  - `V_gg` from corrected decomposed WL density-density;
  - `V_gR` and `V_RR` from residual MWG effective centers/widths;
  - no raw GTO/GTO density-density accepted as final operator data;
  - final density-density shape `(520, 520)`;
  - `V_gR` shape `(517, 3)`;
  - `V_RR` shape `(3, 3)`;
  - symmetry error `0.0`;
  - finite entries.

Component comparison:

- Side13 gausslet-only RHF:
  - one-electron energy `-3.8457802974697652`
  - electron-electron energy `1.0092822977688516`
  - total `-2.8364979997009137`
- Side13 He + GTO final-basis RHF:
  - one-electron energy `-3.889081757426776`
  - electron-electron energy `1.0246874795868013`
  - total `-2.864394277839975`
- Adding the residual GTO supplement:
  - lowers the one-electron energy by about `0.0433014599570106`
  - raises the electron-electron energy by about `0.0154051818179497`
  - improves the total by about `0.0278962781390613`
- Relative to the He HF reference `-2.861679995612234`, the diagnostic total is lower by about `0.00271428222774084`.

Interpretation:

- The below-HF result should not be accepted against the all-electron HF reference.
- The most likely interpretation from this audit is that the current residual MWG/IDA Hamiltonian is approximate and not guaranteed variational relative to exact all-electron HF.
- The result must remain diagnostic until compared with an old nested/QW MWG oracle for the same or intentionally matched fixture, or with a located White-Lindsey Fig. 8 source table.

White-Lindsey Fig. 8/source audit:

- I did not find a repo-local paper PDF or source table for White-Lindsey Fig. 8 in the searched docs/src/test/tmp/data paths.
- `docs/src/algorithms/qiu_white_residual_gaussian_route.md` documents the Qiu-White MWG algorithm, but does not contain the Fig. 8 numeric table.
- Local notes such as `docs/atomic_core_compression_vs_distortion_tradeoff.md` contain related atomic hybrid diagnostics, but are not the requested Fig. 8 source.
- The manager recollection about `ns = 5`, `d = 0.3`, 447 basis functions, and S-only supplement remains a lead, not a citation.

Validation:

- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - pass: `load ok`
- `git diff --check`
  - pass
- I did not rerun the full side13 He + GTO RHF probe because this was an audit-only pass and no source code changed. The audit used the existing fresh probe artifact from response 004.

Deletion/shrinkage:

- No old production code, test, metadata, or compatibility path became unnecessary from this audit-only pass.
- Nothing was deleted because no clear bug or obsolete path was identified.
- No new test was added; the audit artifact remains under `tmp/work`.
- Remaining stale/duplicate surfaces to retire next:
  - the probe-local RHF helper should be deleted or replaced by a small shared final-basis RHF utility if this route becomes accepted;
  - the residual MWG final-basis density-density path needs an old-MWG/internal-Hamiltonian oracle before acceptance;
  - if a White-Lindsey Fig. 8 source table is located later, the current uncited recollection should be replaced by that source.

Trust boundary:

- No raw GTO density-density fallback was added.
- No full-parent CPB, flat `_qwrg_diatomic_*` authority, direct Cartesian product, ordinary Cartesian IDA, PQS, or generalized final solve was used.
- No public driver/export/artifact path was added.

-- repo-doer@macmini
