Pass 102 response

No-edit RHF contract design audit completed.

Files changed:

- none before this response/log

Read-only inspection used:

- `AGENTS.md` / baton startup docs already re-read after compaction
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/successor_handoff_2026-06-12_pqs_source_box_pivot.md`
- `/Users/srw/Dropbox/chatarchive/handoff/software_packets/gausslet-methods-fundamentals.md`
- `src/pqs_multilayer_complete_core_shell_h1.jl`
- `src/pqs_multilayer_support_density.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`
- recent pass logs/reviews around 097-101, plus older accepted RHF convention
  notes surfaced by search
- `git status --short --branch`

Proposed RHF contract

- RHF should be a separate private complete core/shell PQS contract, not an
  extension of `pqs_multilayer_complete_core_shell_h1_j_payload(...)`.
- Working name for the first compact object/payload boundary:
  `_PQSCompleteCoreShellRHFDiagnosticPayload`, or the module-level equivalent
  `pqs_multilayer_complete_core_shell_rhf_payload(...)`.
- Category should be explicit: private diagnostic/prototype, not production
  route, not public API, not export/artifact behavior.
- It should consume the already route-owned complete core/shell objects:
  source plan, final basis, H1 payload/final H1 matrix, route density inputs,
  Coulomb expansion, electron count / occupation policy, and SCF controls.
- It should return a compact status/blocker/summary plus a private SCF object,
  not a cloud of report aliases.

Module/seam ownership

- Do not put Fock/SCF construction inside the H1/J helper. H1/J should remain a
  one-orbital diagnostic over the accepted density gauge.
- Do not put the SCF loop in `src/pqs_multilayer_support_density.jl`. That file
  should stay limited to support weights and raw support pair numerator inputs.
- Preferred implementation seam: a new private module/file such as
  `src/pqs_multilayer_complete_core_shell_rhf.jl`, included near the current
  H1 helper. It can consume:
  - `pqs_multilayer_complete_core_shell_h1_payload(...)`;
  - `pqs_multilayer_support_weights(...)`;
  - `pqs_multilayer_support_pair_raw_numerator_matrix(...)`;
  - `CartesianFinalBasisRealization.pqs_complete_core_shell_pre_final_density_interaction(...)`.
- The route driver should later consume only the compact RHF payload, probably
  through `_PQSCompleteCoreShellDiagnosticRoutePayload` or a sibling private
  route payload. It should not own Fock algebra.

Density/orbital/occupation conventions

- Density gauge: use the accepted localized pre-final positive-weight density
  gauge for electron-electron contractions.
- Relationship to H1/J: the current H1/J diagnostic uses the same density
  interaction to evaluate one occupied H1 orbital. RHF should reuse the gauge
  convention, but not reuse the H1/J helper as the SCF owner.
- Final basis remains the orthonormal one-body/SCF orbital basis.
- Final-to-prefinal map must be explicit:

```text
c_prefinal = combined_lowdin_cleanup * c_final
```

- Fock construction should keep the same shape documented by earlier accepted
  probes:

```text
G_pre   = 2 * diag(V_pre * n_pre) - rho_pre .* V_pre
F_final = H_final + L' * G_pre * L
```

  where `L` is the final-to-prefinal map from the materialized density
  interaction (`final_to_pre_final_coefficients` / combined Lowdin cleanup).
- Electron count must be an explicit input before code. The safe first contract
  is closed-shell RHF only:
  - `electron_count::Int`
  - even electron count required
  - `nocc = electron_count ÷ 2`
  - spatial orbital occupancy `2`
  - no fractional occupation, spin polarization, open shell, unrestricted HF,
    or charge inference unless separately approved.
- For neutral one-center probes, nuclear charge can motivate the electron count,
  but the RHF payload should not silently infer it as route behavior.

Required RHF inputs

- available `:pqs_multilayer_shell_source_plan`
- available complete core/shell final basis
- materialized complete core/shell H1 payload/final H1 matrix
- route-owned support weights and raw pair numerator terms
- Coulomb expansion coefficients matching the raw pair terms
- materialized pre-final density interaction with positive pre-final weights
- explicit closed-shell electron count / occupation policy
- SCF controls: max iterations, energy tolerance, density/Fock residual
  tolerance, damping/mixing policy if any, initial guess policy, symmetry
  tolerance

Required RHF outputs

- status/blocker and compact summary
- final-basis molecular/orbital coefficients and orbital energies
- occupied-orbital count, occupation vector, electron-count trace checks
- final density matrix summary and pre-final mapped density summary
- Fock matrix symmetry/residual/convergence diagnostics
- SCF iteration count, convergence flag, final residuals, and optional compact
  per-iteration history
- energy components:
  - one-electron electronic energy
  - electron-electron energy
  - electronic total
  - nuclear repulsion term, explicitly zero for one-center when applicable
  - total molecular energy

Convergence/SCF object before route adoption

- A private SCF object should exist before any route adoption, e.g.
  `_PQSCompleteCoreShellRHFSCFResult`.
- It should carry convergence status, blocker, iteration count, tolerances,
  initial guess source, final orbital energies, final coefficients, compact
  density/Fock residuals, energy components, and non-promotion flags.
- Route adoption should require this object to be materialized and summarized;
  the driver should not infer convergence from a single energy scalar.

Forbidden shortcuts

- Do not use signed final-weight division as the density contract.
- Do not use raw no-division density as the density contract.
- Do not use density-normalized pair terms as route authority.
- Do not treat retained diagnostic weights as final IDA/quadrature weights.
- Do not use fixed-block/shell/support-row pair matrices as active route
  authority; they remain oracle/debug.
- Do not make H1/J materialization imply RHF readiness.
- Do not add GTO, IDA/MWG, exports, artifacts, public API, fixture promotion,
  or production route behavior in the first RHF pass.
- Do not hide electron-count or occupation conventions behind defaults.

What remains diagnostic/private

- H1/J payload and report aliases remain private diagnostics.
- Pre-final positive weights remain a private density-interaction gauge for
  this complete core/shell lane, not final IDA weights.
- RHF, when first added, should also be private diagnostic/prototype until a
  separate fixture/performance/science review promotes it.

Smallest safe first implementation pass

- Stop before implementation unless manager accepts the above contract.
- If approved, the smallest first pass should add only a private module-level
  blocked/available RHF input contract object. It should validate and summarize
  required inputs and occupation policy, but not run SCF yet.
- A slightly larger but still bounded follow-up would add a one-step Fock
  builder over an externally supplied final density matrix, returning
  `F_final` plus energy-component diagnostics. Full SCF iteration should come
  after that object boundary exists.

Git status:

```text
## main...origin/main
```

Deletion/shrinkage forecast:

- deleted: none in this no-edit audit.
- simplified: future RHF work can avoid adding RHF fields to the H1/J helper
  and can keep the driver helper to compact payload consumption.
- quarantined: current H1/J and pre-final density gauge remain diagnostic/private;
  older RHF ladder probes remain historical `tmp/work` evidence, not tests.
- not deleted because: no implementation was assigned; H1/J report aliases and
  density helpers are still active diagnostic seams.
- exact remaining caller/blocker: no RHF caller exists yet. The blocker is an
  approved private RHF input/SCF payload boundary with explicit electron-count,
  occupation, density-gauge, and convergence conventions.

-- repo-doer@macmini
