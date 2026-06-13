Pass 105 response

Private RHF input-contract object/helper added. No Fock matrix, SCF iteration,
RHF energy, driver wiring, report fields, exports, artifacts, or fixture
promotion were added.

Files changed:

- `src/GaussletBases.jl`
- `src/pqs_multilayer_complete_core_shell_rhf.jl`
- `test/nested/pqs_multilayer_complete_core_shell_rhf_input_contract_runtests.jl`
- `test/nested/runtests.jl`

Helper/object name and fields:

- Added private helper:
  `_pqs_multilayer_complete_core_shell_rhf_input_contract(...)`
- Returned compact object fields:
  - `object_kind`
  - `status`
  - `blocker`
  - `missing_inputs`
  - `electron_count`
  - `occupation`
  - `fixture_role`
  - `summary`
  - `metadata`
- Available status:
  - `:available_pqs_multilayer_complete_core_shell_rhf_input_contract`
- Blocked status:
  - `:blocked_pqs_multilayer_complete_core_shell_rhf_input_contract`

Validated input contract:

- source plan must be an available `:pqs_multilayer_shell_source_plan`;
- final basis must be an available `:pqs_complete_core_shell_final_basis`;
- H1 payload must be a materialized
  `:pqs_multilayer_complete_core_shell_h1_payload` with a final H1 matrix;
- density inputs must be `:available_complete_core_shell_density_inputs` with
  axis weights and raw pair factor terms;
- Coulomb expansion must be present, and term count compatibility is summarized
  when density-input term count is available;
- `electron_count` is explicit only;
- closed-shell RHF policy only:
  - positive integer electron count;
  - even electron count;
  - `nocc = electron_count ÷ 2`;
  - occupancy `2`;
- fixture role must be explicit and currently one of:
  - `:route_smoke`
  - `:physics_endpoint`

Blockers covered:

- `:missing_source_plan`
- `:missing_final_basis`
- `:missing_h1_payload`
- `:missing_density_inputs`
- `:missing_coulomb_expansion`
- `:missing_electron_count`
- `:unsupported_open_shell_rhf_input`
- `:missing_fixture_role`
- also guarded: `:invalid_electron_count`,
  `:unsupported_fixture_role`, and
  `:coulomb_expansion_term_count_mismatch`

Tests added:

- Added
  `test/nested/pqs_multilayer_complete_core_shell_rhf_input_contract_runtests.jl`.
- The test uses synthetic compact NamedTuple inputs, not the heavy source-box
  dry-run.
- It catches:
  - available closed-shell route-smoke contract;
  - missing electron count blocked at `:missing_electron_count`;
  - odd electron count blocked at `:unsupported_open_shell_rhf_input`;
  - missing fixture role blocked at `:missing_fixture_role`;
  - missing H1 payload blocked at `:missing_h1_payload`.
- Registered the focused test in `test/nested/runtests.jl` because it is a fast
  synthetic contract check.

Validation commands/results:

- `julia --project=. test/nested/pqs_multilayer_complete_core_shell_rhf_input_contract_runtests.jl`
  - passed: `30/30`
  - package precompilation took about 56 s; the test itself reported `0.2s`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - passed: `load ok`
- `git diff --check`
  - passed

Git status:

```text
## main...origin/main
 M src/GaussletBases.jl
 M test/nested/runtests.jl
?? src/pqs_multilayer_complete_core_shell_rhf.jl
?? test/nested/pqs_multilayer_complete_core_shell_rhf_input_contract_runtests.jl
```

Deletion/shrinkage report:

- deleted: none.
- simplified: RHF now has a separate private input-contract boundary, reducing
  pressure to grow H1/J helpers or driver report aliases into RHF/SCF surfaces.
- quarantined: RHF remains private diagnostic/prototype; compact fixture roles
  remain explicit `:route_smoke` or `:physics_endpoint`; H1/J remains
  diagnostic/private.
- not deleted because: this was the first additive boundary pass; no existing
  RHF implementation or report surface existed to delete.
- exact remaining caller/blocker: no production caller yet. The remaining
  blocker before Fock/SCF is an approved private Fock-builder/SCF-result seam
  consuming this input contract, still with explicit fixture role and electron
  count.

-- repo-doer@macmini
