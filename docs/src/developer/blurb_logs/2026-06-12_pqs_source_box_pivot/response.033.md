Pass 033 completed.

Scope:

- Audit-only pass for the proposed `CartesianFinalBasisRealization` extraction.
- No source code, tests, exports, or route behavior were changed.

Files changed:

- `.agent_handoffs/response.033.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.033.md`

Files inspected:

- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
- `src/cartesian_pair_block_materialization/pqs_source_shell_final_basis.jl`
- `src/cartesian_pair_block_materialization/CartesianPairBlockMaterialization.jl`
- `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/*`

Candidate classification:

```text
pqs_source_shell_realization_final_basis
    move_now
    Reason: owns shell/support realization, Lowdin cleanup, and final overlap
    diagnostics. It depends on `CartesianRawProductSources` plus LinearAlgebra,
    not on CPBM internals.

pqs_source_shell_final_one_body_from_boundary_matrix
    move_now
    Reason: pure retained-boundary operator -> final-basis transfer. Input is
    structural `NamedTuple` retained-source one-body output. This is exactly
    the proposed final-basis operator-transfer seam.

pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block
    needs_manager_decision
    Reason: conceptually belongs with final-basis operator transfer, but the
    current signature is concretely typed as `PairBlockMaterializationResult`.
    Moving it before CPBM would require either a small structural-input
    signature change or moving/sharing that CPBM result type. Keeping the exact
    concrete type while also having CPBM import/reexport the new module creates
    an avoidable module-cycle problem.

pqs_source_shell_projected_one_body_matrix
    oracle_helper_move_with_module
    Reason: useful shell-support oracle/comparison helper. It should move with
    the module only if labeled oracle/reference, not production operator path.

pqs_source_shell_final_one_electron_hamiltonian
    leave_in_CPBM_for_now
    Reason: narrow final one-body Hamiltonian assembly is close to the
    final-basis seam, but the near-term plan explicitly says it may stay in
    CPBM for the first extraction. Leaving it avoids expanding the first module
    move beyond realization/operator transfer.
```

Proposed module layout:

```text
src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl
src/cartesian_final_basis_realization/pqs_source_shell_final_basis.jl
```

Proposed include order:

- Include `CartesianFinalBasisRealization` in `src/GaussletBases.jl` after
  `CartesianRawProductSources` and before `CartesianPairBlockMaterialization`.
- The new module should use `..CartesianRawProductSources` and `LinearAlgebra`.
- Avoid depending on `CartesianPairBlockMaterialization` from the new module.

CPBM compatibility choice:

- Preserve existing internal compatibility with
  `GaussletBases.CartesianPairBlockMaterialization.pqs_source_shell_*` by
  adding CPBM aliases/reexports to the new module functions.
- Suggested CPBM shape:

```julia
using ..CartesianFinalBasisRealization
const CFBR = CartesianFinalBasisRealization
const pqs_source_shell_realization_final_basis =
    CFBR.pqs_source_shell_realization_final_basis
```

- Repeat for moved functions and keep CPBM exports stable for now.
- Do not add broad compatibility adapters; use direct aliases only.

Tests and caller impact:

- Active source callers are all definitions/exports in CPBM plus the existing
  CPBM contract section in
  `test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
  around the final-basis block.
- Historical references in `docs/src/developer/blurb_logs/...` are log history
  and should not be rewritten.
- No test caller must change if CPBM aliases are kept.
- A later cleanup should move or shrink the final-basis section out of the large
  CPBM contract file into a compact FinalBasisRealization module-contract test.

CPBM contract-test sections that can move/shrink:

- The section starting at the synthetic `final_basis_source_box` fixture through
  final Hamiltonian assembly is now mixed: final-basis realization/operator
  transfer plus CPBM retained-source inputs.
- After extraction, keep only a CPBM alias/load smoke if needed, and move the
  actual final-basis checks to a compact module-contract file.
- Keep the shell-support projected helper checks as oracle-only coverage, not
  as production path pressure.

Validation target for eventual mechanical extraction:

- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- A compact FinalBasisRealization contract test covering:
  - `q=5/L=5` source dims/count `5 x 5 x 5 / 125`;
  - retained boundary count `98`;
  - final retained count `98`;
  - final overlap identity;
  - final boundary overlap/kinetic transfer;
  - by-center nuclear transfer if the type-boundary decision is resolved;
  - no generalized final solve, IDA, RHF, driver, export, or artifact claim.
- A small CPBM compatibility check only if aliases are introduced.
- `git diff --check`.

Validation run for this audit:

- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.
- No Julia tests were run because this was audit-only and no executable code was edited.

Deletion/shrinkage report:

- This audit identifies CPBM ownership that can be removed next pass:
  shell-realization final-basis construction, Lowdin final overlap diagnostics,
  retained-boundary one-body final transfer, and oracle shell-support projection.
- The existing CPBM contract-test final-basis section can move/shrink instead
  of being duplicated. It should become a compact FinalBasisRealization contract
  test plus, at most, a CPBM alias smoke.
- `pqs_source_shell_projected_one_body_matrix` should remain oracle-only in the
  new module.
- No code was deleted in this audit because the blurb explicitly prohibited code
  movement.
- Remaining stale/duplicate surface for the next pass: CPBM still owns the
  final-basis realization file and docs comments describe it as a CPBM file.
- Direct retained-boundary overlap/kinetic kernels remain a later pass after
  the module extraction.

-- repo-doer@macmini
