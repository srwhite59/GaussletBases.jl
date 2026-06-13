Pass 150 review - accepted

Response reviewed:
- `.agent_handoffs/response.150.md`

Summary:
- No files were edited and no Julia commands were run.
- The audit separates two different post-H1 concepts:
  - H1-J: a private one-orbital self-Coulomb diagnostic scalar.
  - Ham input / electron-electron payload: reusable data needed for a downstream
    Hamiltonian comparison.

Accepted decision:
- Do not implement private Be2/PQS H1-J first.
- Implement a private Be2/PQS Ham-input/electron-electron payload first.

Reasoning:
- H1-J would mostly consume the same density interaction needed by a future
  Hamiltonian handoff, then report a scalar diagnostic.
- The medium-term Be2 WL/PQS comparison needs H1 plus an electron-electron
  representation and ordering/convention metadata, not just the H1-J scalar.
- A Ham-input payload can keep the useful density-interaction object route-owned
  without claiming public export, HamV6 compatibility, RHF readiness, or CR2
  readiness.

Next step:
- Add `_PQSDiatomicCompleteCoreShellHamInputPayload`.
- Build only the private support-density provenance, support weights, support
  raw pair numerator, pre-final density interaction, and compact metadata.
- Do not build H1-J, full Ham payload, public API, exports, hfdmrg, CR2, or RHF.

-- repo-manager@macmini
