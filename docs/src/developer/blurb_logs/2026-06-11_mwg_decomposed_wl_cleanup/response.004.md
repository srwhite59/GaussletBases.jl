What changed:

- Added and ran a developer probe:
  - `tmp/work/side13_he_gto_final_basis_rhf_probe.jl`
- The probe writes:
  - `tmp/work/side13_he_gto_final_basis_rhf_summary.txt`
- No production code was changed.
- No acceptance test was added.

Probe result:

- Fixture:
  - `Z = 2`
  - `q/ns = 5/5`
  - `d = 0.1`
  - `s = 1.0`
  - side count `13`
  - endpoints `(-8.565228460168399, 8.565228460168399)`
- Decomposed WL inventory:
  - retained gausslet dimension `517`
  - unit count `105`
  - pair count `5565`
- GTO supplement:
  - He cc-pVTZ `lmax = 0`
  - raw supplement count `3`
  - retained residual count `3`
  - dropped supplement count `0`
  - final dimension `520`
- Final-basis diagnostics:
  - final overlap identity error `2.4938717757549966e-11`
  - final overlap rank `520`
  - final overlap condition estimate `1.000000000026454`
  - final one-electron Hamiltonian shape `(520, 520)`
  - final density-density shape `(520, 520)`
  - `V_gR` shape `(517, 3)`
  - `V_RR` shape `(3, 3)`
  - final density-density symmetry error `0.0`
  - final density-density finite `true`
- RHF:
  - status `converged_final_basis_he_gto_restricted_hartree_fock`
  - iterations `19`
  - solve kind `ordinary_symmetric`
  - one-electron energy `-3.889081757426776`
  - electron-electron energy `1.0246874795868013`
  - total RHF energy `-2.864394277839975`
  - density trace `1.0000000000000007`
  - electron count `2.0000000000000013`

Comparison:

- Side13 gausslet-only RHF total: `-2.8364979997009137`
- He + GTO improvement over side13 gausslet-only: `0.02789627813906126`
- He HF reference: `-2.861679995612234`
- Error from He HF reference: `-0.0027142822277408385`

Decision:

- The diagnostic improves substantially over side13 gausslet-only.
- It is not suitable as a future acceptance fixture yet because it is below the He HF reference by about `2.7 mHa`.
- Per the blurb decision rule, this should stop at diagnostic status and the next pass should audit the final-basis density-density or RHF convention before any acceptance claim.

Timing:

- parent axis setup: `0.053445542 s`
- inventory: `2.145188041 s`
- one-electron combined/final-basis build: `93.916680916 s`
- residual moment and MWG representation: `24.80738875 s`
- density-density build: `20.676718417 s`
- RHF solve: `2.099251541 s`
- total elapsed: `144.885561958 s`

Anti-fallback checks:

- raw GTO density-density used as final operator: `false`
- generalized final-basis solve: `false`
- full-parent CPB used: `false`
- direct Cartesian product assembly used: `false`
- ordinary Cartesian IDA operators used: `false`
- PQS transforms materialized: `false`
- exports/artifacts: `false`

Validation:

- `julia --project=. tmp/work/side13_he_gto_final_basis_rhf_probe.jl`
  - pass; summary artifact written
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - pass: `load ok`
- `git diff --check`
  - pass

Deletion/shrinkage:

- No old production code, test, metadata, or compatibility path became unnecessary from this probe-only pass.
- Nothing was deleted because the pass only exercised the newly materialized final-basis density-density path and found a physics/convention stop condition.
- No new test was added; the probe is exploratory and remains under `tmp/work`.
- Remaining stale/duplicate surfaces to retire next:
  - if the convention audit confirms the final-basis RHF path, the probe-local RHF helper should be replaced by a small shared final-basis RHF utility or deleted with the probe;
  - if the convention audit finds an error, the residual MWG density-density/RHF convention surface should be corrected before any acceptance test is added.

-- repo-doer@macmini
