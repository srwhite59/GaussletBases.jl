Pass 145 review - accepted

Response reviewed:
- `.agent_handoffs/response.145.md`

Summary:
- No files were edited and no Julia commands were run.
- The audit identified the exact consumer contract for the existing one-center
  `:pqs_multilayer_shell_source_plan`.
- It also confirmed that the new diatomic source-realization payload should not
  pretend to be that one-center object kind.

Key decision:
- Be2/PQS should get its own private source-plan object:
  `:pqs_diatomic_complete_core_shell_source_plan`.
- The existing one-center shell source-plan remains valid for the compact H1
  fixture and should not be retrofitted to cover the diatomic product/body plus
  left/right PQS route.

Accepted rationale:
- The old source-plan object encodes one-center nested shell semantics:
  rectangular core/outer boxes, shell records, layer count, and coverage checks.
- The Be2/PQS source realization is a different route shape:
  product core/body sector plus separated left/right PQS source sectors.
- Satisfying old guard clauses by fabricating old fields would preserve syntax
  while breaking route semantics.

Next step:
- Add the private diatomic source-plan object and materializer from the
  source-realization/raw-box payload.
- Stop before final-basis/H1/H1-J/Ham materialization.
- The next blocker should become a clear missing diatomic final-basis/consumer
  seam, not the old `:pqs_multilayer_shell_source_plan_adapter_contract`.

-- repo-manager@macmini
