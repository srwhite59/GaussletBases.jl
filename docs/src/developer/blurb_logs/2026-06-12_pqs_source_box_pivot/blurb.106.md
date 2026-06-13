Pass 106 - one-step RHF density contraction seam, private only

Baseline:

- Current pushed HEAD should include `9d6a0da1 Add PQS RHF input contract`.
- Pass 105 added `_pqs_multilayer_complete_core_shell_rhf_input_contract(...)`
  in `src/pqs_multilayer_complete_core_shell_rhf.jl`.
- Governing boundaries:
  - `docs/src/developer/pqs_source_box_operator_framework.md`
  - `docs/src/developer/pqs_source_box_fixture_policy.md`
  - `docs/src/developer/successor_handoff_2026-06-12_pqs_source_box_pivot.md`
  - `docs/src/developer/manager_doer_collaboration_contract_2026-05-26.md`

Task:

Add the next private seam toward RHF: a one-step closed-shell
Fock/energy-style contraction over an externally supplied final-basis density.

This is still not SCF. The goal is to make the density contraction convention
explicit and testable before any iteration loop exists.

Important decision rule:

- First inspect the existing complete core/shell density helpers, especially:
  - `pqs_complete_core_shell_pre_final_density_interaction(...)`
  - `pqs_complete_core_shell_pre_final_orbital_self_coulomb(...)`
  - `_pqs_complete_core_shell_restricted_one_orbital_interaction_energy(...)`
  in `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`.
- If the existing two-index pre-final density interaction gives a clear
  closed-shell density contraction, implement the smallest private helper.
- If the contraction is not unambiguous, do not invent a physics convention.
  Report a no-edit blocker and the exact ambiguity.

Preferred implementation surface if unambiguous:

- Stay in `src/pqs_multilayer_complete_core_shell_rhf.jl`.
- Add a private helper with a name like:
  `_pqs_multilayer_complete_core_shell_rhf_one_step_payload(...)`
  or a similarly narrow name.
- Inputs should be explicit and compact:
  - available RHF input contract from pass 105;
  - materialized density interaction, or the H1/J payload containing it;
  - externally supplied final-basis density matrix;
  - optional metadata.
- Validate:
  - RHF input contract status is available;
  - density interaction status is materialized
    `:materialized_pqs_complete_core_shell_pre_final_density_interaction`;
  - final density is a finite square matrix with dimension equal to the final
    dimension in the contract;
  - final density is symmetric within a small tolerance;
  - electron trace matches `electron_count` within tolerance, or report a
    precise blocker;
  - fixture role is carried through unchanged.

Expected payload shape:

- `object_kind`
- `status`
- `blocker`
- `missing_inputs`
- `final_density`
- one-step matrix if materialized, e.g. `fock_matrix` or
  `effective_fock_matrix`
- compact energy components if materialized, e.g. `one_body_energy`,
  `two_body_energy`, `total_energy`
- `summary`
- `metadata`

The summary must stay compact. Do not duplicate large matrices in summary or
metadata.

Required nonclaims:

- `scf_materialized = false`
- `rhf_converged = false`
- `driver_route_materialized = false`
- `exports_materialized = false`
- `artifacts_materialized = false`
- no public API/export
- no route-driver/report fields

Tests:

- Add one small focused synthetic test file if code is added.
- Use tiny matrices and synthetic NamedTuple inputs; do not run the heavy
  source-box dry-run.
- Cover:
  - available one-step contraction for a closed-shell two-electron density;
  - missing/unavailable RHF input contract blocker;
  - nonsymmetric density blocker;
  - electron-trace mismatch blocker.
- Keep numerical assertions simple and derived from the tiny fixture.

Exclusions:

- Do not add an SCF loop.
- Do not solve RHF.
- Do not infer electron count from nuclei or charges.
- Do not wire the route driver.
- Do not add driver options or report aliases.
- Do not touch GTO, IDA/MWG, exports, artifacts, or production route behavior.
- Do not promote compact route-smoke fixtures to physics acceptance.
- Do not request interactive escalation/approval; if a required command or
  design decision truly needs approval, stop by writing `ATTENTION.md`.

Validation:

- If code is added:
  - `julia --project=. test/nested/<new focused test>.jl`
  - `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - `git diff --check`
- If no code is added due to ambiguity:
  - `git status --short --branch`
  - report the ambiguity and exact recommended next decision.

Report back:

- Whether the contraction was unambiguous or blocked.
- Helper/object names if added.
- Status/blocker labels.
- Validation commands/results.
- Git status.
- Deletion/shrinkage report:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
