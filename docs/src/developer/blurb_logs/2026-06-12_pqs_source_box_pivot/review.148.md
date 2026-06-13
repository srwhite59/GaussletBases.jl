Pass 148 review - accepted

Response reviewed:
- `.agent_handoffs/response.148.md`

Summary:
- No files were edited and no Julia commands were run.
- The audit confirms that the private Be2/PQS source plan and final basis have
  the numerical ingredients needed for H1, but the old one-center support/H1
  helpers correctly hard-gate on `:pqs_multilayer_shell_source_plan`.

Accepted decision:
- Do not broaden the old one-center source-plan guard as the first move.
- Add a private diatomic H1 payload that reuses lower matrix/transfer primitives
  while preserving the distinct `:pqs_diatomic_complete_core_shell_source_plan`
  contract.

Key constraints for the next pass:
- H1 can be materialized as a private diagnostic route object.
- H1-J, density-density, RHF, WL, public API, exports, artifacts, hfdmrg, and
  CR2 remain out of scope.
- The old one-center source-plan/H1 helpers remain untouched unless a tiny local
  extraction is clearly needed and behavior-preserving.

Next step:
- Implement `_PQSDiatomicCompleteCoreShellH1Payload` and a private helper that
  consumes the diatomic source plan plus final basis, builds support kinetic and
  electron-nuclear matrices, transfers them to the final basis, and solves H1.
- Keep the full Ham payload blocked after H1, with a sharper missing H1-J or
  electron-electron consumer blocker.

-- repo-manager@macmini
