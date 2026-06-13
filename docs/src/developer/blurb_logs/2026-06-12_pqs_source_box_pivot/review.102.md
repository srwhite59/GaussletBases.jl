Pass 102 review

Accepted no-edit RHF contract design.

Key accepted points:

- RHF must be a separate private complete core/shell PQS contract, not an
  extension of `pqs_multilayer_complete_core_shell_h1_j_payload(...)`.
- The first RHF object should be private diagnostic/prototype only.
- Electron count and closed-shell occupation must be explicit inputs.
- Fock/SCF construction should not live in `pqs_multilayer_support_density.jl`
  or the H1/J helper.
- The accepted pre-final positive-weight density gauge can inform RHF, but H1/J
  materialization does not imply RHF readiness.

Do not implement RHF yet. The next pass should define the fixture/science rule
boundary so the current compact H1/J fixture does not become an accepted physics
gate by inertia.

-- repo-manager@macmini
