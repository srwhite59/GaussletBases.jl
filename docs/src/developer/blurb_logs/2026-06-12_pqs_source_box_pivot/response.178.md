Pass 178 response - WL-matched multi-shell He q=5 PQS target audit

Current git status at audit close:

```text
## main...origin/main
?? docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/blurb.178.md
```

Files changed:
- No source, test, tracked generator, or artifact files were edited.
- This response and the curated response copy are the only files written for the pass.

Read-only inspection commands used:
- `sed` on the requested developer docs and blurb responses.
- `sed` / `rg` / `wc -l` on the requested source and test surfaces.
- `git status --short --branch`.
- No Julia command was run.

Old WL He construction inventory:

1. Current compact decomposed WL q/ns = 5/5 readiness path
   - File: `test/nested/cartesian_wl_gausslet_he_atom_acceptance_runtests.jl`
   - Construction: one-center He, `Z = 2`, q/ns = 5/5, `d = 0.2`, `s = sqrt(dZ) = 0.6324555320336759`, `tail_spacing = 10.0`, parent side count 7.
   - Geometry/inventory: mapped parent axis count `(7, 7, 7)`, direct core retained range `1:125`, one shell retained range `126:223`, retained dimension 223, 27 decomposed units, 378 upper-triangular unit pairs.
   - Route status: current decomposed WL/shellification-backed route, not old nested/QW; no full-parent CPB, direct Cartesian, or ordinary Cartesian IDA fallback.
   - Recorded values: compact H1 convention diagnostic in `numerical_contracts.md` records He+ H1 about `-1.878770102537269` and self-Coulomb about `1.4202542835594492`; the current tracked test asserts closed-shell RHF total `-2.3944175346639884` with electron-electron contribution `1.3106054775285387`.
   - Role: route readiness / convention acceptance for the decomposed WL plumbing, not the high-accuracy Fig. 8 target.

2. Larger decomposed WL gausslet-only side13 path
   - Construction: `AsinhMapping(c = 0.1, s = 1.0, tail_spacing = 10.0)`, parent side 13, endpoints about `+/-8.565228460168399`.
   - Inventory: retained dimension 517, 105 retained units, 5,565 unit pairs, source kind `:cartesian_shellification_retained_unit_pair_plan`.
   - Recorded values: H1 `-1.9748150892830352`, IDA self-Coulomb `1.2158294767735702`, RHF total `-2.8364979997009137`.
   - Role: current decomposed WL exploratory larger-box path, not a permanent acceptance gate.

3. Best fixed-ns = 5 decomposed WL gausslet-only exploratory point
   - Construction: `d = 0.075`, `s = 0.75`, side count 17.
   - Inventory: retained dimension 713, 157 retained units, 12,403 unit pairs.
   - Recorded values: H1 error about `0.002753`, self-Coulomb error about `0.003389`, RHF total about `-2.858531351214`.
   - Role: best current gausslet-only exploratory point, but too heavy/not reviewed as a routine gate.

4. Fig. 8 old nested/QW MWG n_s = 5, d = 0.3 target
   - Construction: He RHF with AHGBS-9 S-only supplement, `doside = n_s = 5`, `corespacing = d = 0.3`, `dwidth/tail_spacing = 10.0`, `gscalefac = sqrt(2.0)`.
   - Repo constructor used: `MappedUniformBasisSpec(:G10; count = 11, mapping = white_lindsey_atomic_mapping(Z = 2, d = 0.3, tail_spacing = 10.0), reference_spacing = 1.0)`.
   - Mapping note: repo constructor gives `AsinhMapping(a=0.38729833462074165, s=0.7745966692414834, tail_spacing=10.0)`, `a*s = 0.3`, endpoints about `+/-5.892850307983052`. Legacy logs mention `basradius = 8.0` and backbone coordinates about `+/-5.470267`, so exact paper mapping/log convention is not identical.
   - Inventory: old fixed-block gausslet dimension 419, AHGBS-9 S-only residual count 28, final dimension 447.
   - Reconstructed shell interpretation: fixed `5^3 = 125` direct core plus three fixed-q boundary shell sectors, `3 * 98`, gives `125 + 294 = 419`; AHGBS residuals are a separate supplement layer, not part of the gausslet-only PQS target.
   - Recorded values: H1 lowest `-1.9999998632985623`, IDA self-Coulomb `1.2496940228276845`, RHF one-electron `-3.887511967391441`, electron-electron `1.0254098228577182`, RHF total `-2.862102144533723`.
   - Accuracy: about `-0.558 mHa` from the Fig. 8 `n_s=5,d=0.3` row.
   - Route status: old nested/QW MWG oracle/reproduction, not current decomposed WL route authority.

5. Fig. 8 old nested/QW MWG n_s = 7, d = 0.10 target
   - Construction: AHGBS-9 S-only, `n_s = 7`, `d = 0.10`, tail/dwidth 10.0.
   - Repo constructor: count 21, endpoints about `+/-7.117581897720067`.
   - Inventory: gausslet count 1869, residual count 28, final dimension 1897. Fixed-q shell interpretation is `7^3 = 343` core plus seven `218`-retained shells: `343 + 7*218 = 1869`.
   - Recorded values: H1 `-1.9999998641649084`, J `1.2500024042126645`, RHF total `-2.861673961528321`.
   - Accuracy: `+1.716095156645281e-6 Ha` versus the Fig. 8 row and `+6.034083917860755e-6 Ha` versus the He HF reference.
   - Route status: old nested/QW MWG oracle/reproduction; accurate enough to close the immediate atomic He accuracy check.

Correct first PQS q=5 atom target:

The first serious PQS atom target should be the gausslet-only fixed-q multi-shell He construction, not the old q-ladder and not the compact one-shell 223 gate:

```text
He, Z = 2
q/ns = 5/5
WL-style fixed local 5^3 core
multiple fixed-q=5 surrounding shells
start with the n_s=5,d=0.3 Fig. 8 geometry family
expected gausslet-only retained dimension: 419 = 125 + 3*98
AHGBS-9 S-only residuals: separate later comparison layer, +28 -> 447
```

The current compact `5^3 core + one 98 shell = 223` target is still useful as a route seam/readiness gate, but it is not the desired physics target. The desired first PQS comparison is fixed 5^3 core plus the multiple q=5 shells needed to match the old WL n_s=5 inventory.

Can current PQS express that target?

Not exactly.

Current PQS can express:
- the compact one-shell `223` H1/Ham diagnostic path in `test/nested/pqs_direct_retained_final_h1_runtests.jl` and `test/nested/pqs_source_box_route_driver_complete_core_shell_ham_payload_runtests.jl`;
- route-owned multi-layer planning through `pqs_multilayer_shell_region_plan(...)` and `pqs_multilayer_shell_source_plan(...)`;
- final basis, H1, H1/J diagnostics, and private RHF diagnostic helpers through `src/pqs_multilayer_complete_core_shell_h1.jl` and `src/pqs_multilayer_complete_core_shell_rhf.jl`;
- side13/core7 multi-layer PQS probes with final dimension 1549 and RHF `-2.8372556463894707`.

Current PQS does not yet express the WL fixed-q shell inventory. In `src/pqs_multilayer_shell_source_plan.jl`, `_pqs_multilayer_realize_shell_source_plan(...)` sets each layer's `q` and `raw_source_dims` from `length.(inner_box)`. Therefore a side11/core5 three-layer plan would grow layer source sizes as q=5, q=7, q=9, with retained shell counts `98 + 218 + 386`, not the WL fixed q=5 counts `98 + 98 + 98`. That would produce a different PQS basis from the old WL 419-gausslet target.

Smallest next implementation seam:

Add a fixed-source-shape path to the existing one-center multi-layer PQS source plan. The smallest pass should teach the shellification/lowering-backed `pqs_multilayer_shell_source_plan(bundles, region_plan; ...)` path to consume each layer's lowering/contract `source_mode_shape` when present, rather than always deriving `q` from the layer inner-box side. The explicit-box bridge can keep its current growing-q behavior for compatibility/probe use. A focused tmp/work probe should build the He `n_s=5,d=0.3` style parent with a fixed 5^3 core and three q=5 shell layers, assert source-plan/final-basis dimensions `125 + 3*98 = 419`, then run H1/J only and compare directionally against the old nested/QW values. Do not run RHF or add AHGBS residuals in that first pass.

Deletion/shrinkage candidates for the next pass:
- Simplify `test/nested/pqs_direct_retained_final_h1_runtests.jl` by removing the explicit-box bridge duplicate (`explicit_box_plan`, `explicit_box_final_basis`, and the equality assertions against the active shellification/lowering plan) once the fixed-q route-owned source plan is validated.
- Keep the 223 test only as a compact seam smoke, or replace it later with one compact fixed-q multi-shell H1 smoke if the new target becomes the reviewed gate.
- Do not add a new long-term test for every old probe. If a fixed-q multi-shell H1 gate is promoted, shrink older q-ladder/explicit-box/oracle pressure instead.
- `test/nested/pqs_source_box_route_driver_complete_core_shell_ham_payload_runtests.jl` is a later shrink candidate because it asserts many private Ham-payload labels for the 223 diagnostic path; it should not be the main physics guard once a WL-matched atom endpoint exists.

Stale assumptions corrected:
- The Be2/Cr2 artifact work should wait; atom-first He physics is the cleaner comparison target.
- `final_dimension == 223` is a route seam, not a physics acceptance claim.
- The old PQS q-ladder grew q/core size with one shell; it is not the desired WL-style construction.
- Existing PQS multi-layer support is real but currently grows q per layer; it is not yet a fixed-q WL-matched shell inventory.
- AHGBS-9 S-only residuals belong to a later supplemented comparison layer; the first PQS target should be gausslet-only.

Blocker:
- No scientific policy choice blocks the next engineering pass. The first blocker is concrete: fixed-q source-mode shape is not yet consumed by the route-owned multi-layer PQS source plan.

-- repo-doer@macmini
