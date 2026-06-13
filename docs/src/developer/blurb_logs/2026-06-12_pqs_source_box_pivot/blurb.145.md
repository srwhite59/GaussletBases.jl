Pass 145 - audit diatomic source-plan adapter contract

Purpose:

Do a no-edit audit of the remaining Be2/PQS blocker:

```text
:missing_pqs_multilayer_shell_source_plan_adapter_contract
```

The question is not whether the new diatomic source-realization payload has
useful facts. It does. The question is whether those facts can honestly be
adapted to the existing one-center `:pqs_multilayer_shell_source_plan` contract,
or whether Be2/PQS needs a distinct private source-plan object and a later
consumer seam.

Read first:

- `src/pqs_multilayer_complete_core_shell_h1.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `test/nested/pqs_direct_retained_final_h1_runtests.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/fixture_role_policy_pqs_h1_h1j_rhf.md`

Audit questions:

1. What exact source-plan object fields/properties are consumed by:
   - `pqs_multilayer_complete_core_shell_final_basis(...)`
   - `pqs_multilayer_complete_core_shell_h1_payload(...)`
   - complete-core/shell density/H1-J helpers, if relevant

2. Which of those requirements are already represented by the new
   `_PQSDiatomicCompleteCoreShellSourceRealizationPayload`?

3. Which requirements are missing as real structured data, not scalar aliases?
   Pay particular attention to:
   - shell/source coefficient matrices versus shape metadata
   - support row order versus retained column order
   - precleanup retained ranges versus route retained ranges
   - final-to-pre-final transform expectations
   - product/core/body placement
   - support-one-body and support-density helper expectations
   - object kind and convention labels

4. Is it semantically correct for the Be2/PQS realization to return
   `object_kind = :pqs_multilayer_shell_source_plan`?

   If yes, specify the minimal structured adapter and exact facts it must
   materialize.

   If no, specify the smallest new private object/consumer seam needed instead.

5. Identify the smallest safe implementation pass after this audit:
   - adapter to existing source-plan contract, or
   - new diatomic source-plan object, or
   - blocked payload with a sharper missing fact.

Trust boundary:

- No code edits.
- No commits.
- No route-driver behavior changes.
- No source-plan/final-basis/H1/H1-J/Ham materialization.
- No RHF/SCF/DIIS work.
- No WL payload, public API, exports, artifacts, hfdmrg, or CR2 execution.
- Do not promote shell/support-row contraction or raw product-box probes to
  route authority.
- Do not reinterpret retained diagnostic/self-integral weights as
  IDA/quadrature weights.
- Do not ask for interactive approval during unattended baton work. If approval
  would be required, write `.agent_handoffs/ATTENTION.md` and stop.

Validation:

- Read-only audit only.
- `git status --short --branch`
- No Julia commands are required unless you need a short local introspection
  probe. If you do run one, report why and keep it local/ignored.

Report back:

- Existing source-plan consumer contract, with file/line references.
- Which required fields are satisfied by the diatomic source-realization payload.
- Which required fields remain missing.
- Clear decision: existing adapter is semantically valid, or new private object
  is required.
- Recommended next implementation pass.
- Git status.
- Deletion/shrinkage forecast:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
