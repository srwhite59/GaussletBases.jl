Pass 152 review - accepted

Response reviewed:
- `.agent_handoffs/response.152.md`

Summary:
- No files were edited and no Julia commands were run.
- The audit recommends a private factor-first Hamiltonian handoff payload, not
  a public export/HamV6 adapter and not a dense `Vee` shortcut.

Accepted decision:
- Add `_PQSDiatomicCompleteCoreShellHamiltonianHandoffPayload`.
- The handoff should carry references to the private source plan, final basis,
  H1 payload, and Ham-input/density interaction payload.
- It should expose compact inspect metadata for ordering, gauge, raw-pair
  convention, nuclear metadata, electron-count/spin status, and nonclaims.

Important boundary:
- This handoff can make the Be2/PQS constructor output inspectable.
- It still must not claim CR2/HFDMRG readiness, public export readiness, HamV6
  compatibility, RHF, or a completed WL/PQS comparison.

Next step:
- Implement the private handoff payload and move readiness to a sharper
  consumer/export blocker such as
  `:missing_diatomic_hamiltonian_consumer_contract`.

-- repo-manager@macmini
