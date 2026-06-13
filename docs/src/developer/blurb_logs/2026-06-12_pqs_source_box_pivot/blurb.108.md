Pass 108 - private RHF initial-density payload, no SCF

Baseline:

- Current pushed HEAD should include `0ae8190f Cross-check PQS RHF one-step convention`.
- Private RHF surfaces currently live in
  `src/pqs_multilayer_complete_core_shell_rhf.jl`.
- Existing private helpers:
  - `_pqs_multilayer_complete_core_shell_rhf_input_contract(...)`
  - `_pqs_multilayer_complete_core_shell_rhf_one_step_payload(...)`

Task:

Add a private closed-shell initial-density payload for the future SCF loop.

Goal:

Given an available RHF input contract and a materialized H1 payload, build the
ordinary H1-Aufbau spin-summed final-basis density:

```text
H1 matrix -> symmetric eigensolve -> occupy lowest nocc orbitals with occupancy 2
P0 = 2 * C_occ * C_occ'
```

This is only the initial density object. It must not call the one-step Fock
helper, run SCF, compute RHF convergence, wire the route driver, or add report
fields.

Preferred implementation surface:

- Stay in `src/pqs_multilayer_complete_core_shell_rhf.jl`.
- Add a private helper with a name like:
  `_pqs_multilayer_complete_core_shell_rhf_initial_density_payload(...)`.
- Inputs:
  - available RHF input contract;
  - materialized H1 payload, or an H1/J/diagnostic payload from which the H1
    payload can already be extracted by existing private helper logic;
  - optional metadata.
- Reuse existing private property/extraction helpers where possible.

Validation and behavior:

- Validate input contract status/kind.
- Validate H1 payload status/kind and finite square H1 matrix.
- Validate H1 matrix dimension equals the input contract final dimension.
- Validate `nocc <= final_dimension`.
- Symmetrize the H1 matrix before diagonalization, as existing H1 solve code
  does.
- Produce:
  - occupied orbital coefficient matrix;
  - spin-summed final density matrix;
  - eigenvalues, or at least compact occupied eigenvalue summary;
  - electron trace;
  - compact summary and metadata.

Expected status/blocker labels:

- available/materialized:
  `:materialized_pqs_multilayer_complete_core_shell_rhf_initial_density_payload`
- blocked:
  `:blocked_pqs_multilayer_complete_core_shell_rhf_initial_density_payload`
- blockers such as:
  - `:missing_rhf_input_contract`
  - `:missing_h1_payload`
  - `:h1_dimension_mismatch`
  - `:nonfinite_h1_matrix`
  - `:insufficient_final_dimension_for_occupation`

Required nonclaims:

- `initial_density_source = :h1_aufbau`
- `scf_materialized = false`
- `rhf_converged = false`
- `rhf_energy_materialized = false`
- `driver_route_materialized = false`
- `exports_materialized = false`
- `artifacts_materialized = false`

Tests:

- Add one focused synthetic test file or extend the existing RHF one-step test
  only if that is clearly cleaner.
- Prefer a new focused file if the initial-density contract would otherwise
  clutter the one-step convention test.
- Use tiny synthetic matrices.
- Cover:
  - available closed-shell two-electron density from diagonal H1;
  - density trace equals electron count;
  - density is symmetric and idempotent in the closed-shell convention
    (`P^2 ≈ 2P` for occupancy 2);
  - missing input contract blocker;
  - insufficient dimension for occupation blocker if easy.

Exclusions:

- Do not add SCF.
- Do not call `_pqs_multilayer_complete_core_shell_rhf_one_step_payload`.
- Do not compute RHF energy.
- Do not route-wire or add report aliases/options.
- Do not infer electron count from nuclei.
- Do not touch GTO, IDA/MWG, exports, artifacts, or production route behavior.
- Do not run the heavy source-box dry-run.
- Do not request interactive escalation/approval; if approval is truly needed,
  write `ATTENTION.md` and stop.

Validation:

- `julia --project=. test/nested/<focused initial density test>.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Report back:

- Helper/object names.
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
