Pass 113 manager review

Accepted as a local ignored probe.

The compact real PQS route-owned diagnostic fixture reaches the private RHF
input contract, initial density, initial one-step payload, and SCF payload with
explicit `electron_count = 4` and `fixture_role = :route_smoke`. The private SCF
payload does not converge within 25 iterations without damping/mixing:

- status: `:blocked_pqs_multilayer_complete_core_shell_rhf_scf_payload`
- blocker: `:scf_not_converged`
- last total energy: `-10.032171604675229`
- last one-body energy: `-14.21915637345403`
- last two-body energy: `4.1869847687788013`
- last density change: `4.4593642533563571e-05`
- last energy change: `1.6560067273019285e-07`

This is useful route-smoke diagnostic evidence, not physics acceptance. No route
wiring or tracked code/test/doc change was made.

Decision:

- Do not add damping/mixing yet.
- Do not wire RHF into the route driver yet.
- Next pass should collect a compact iteration trace from the same ignored
  local probe path, so we can distinguish slow monotone convergence from
  oscillation or instability before choosing a damping/mixing seam.

Validation/status:

- Local ignored probe ran in `104.527284583` seconds.
- Reported git status clean and even with origin/main, apart from ignored
  `tmp/work` artifacts.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: none; this was evidence gathering.
- quarantined: ignored probe artifacts remain local diagnostic evidence.
- not deleted because: the summary is useful for the next SCF behavior decision.
- exact remaining caller/blocker: compact real PQS private SCF is blocked by
  `:scf_not_converged`; next need is iteration-trace evidence before any
  damping/mixing implementation.

-- repo-manager@macmini
