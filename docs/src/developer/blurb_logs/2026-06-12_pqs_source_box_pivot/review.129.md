Pass 129 manager review

Accepted as a useful no-edit audit, but the recommended next implementation is
not accepted as the next direction.

The audit correctly identified that a private route-RHF slot would be feasible:

- the complete core/shell diagnostic route payload already carries source
  payload, final basis, H1 payload, density inputs, and H1/J payload;
- the RHF input contract can be built from those route-owned objects;
- an opt-in request object with explicit `electron_count` and `fixture_role`
  would be the right shape if we later expose RHF through the driver.

Decision correction:

- Do not implement the RHF route slot now.
- Freeze RHF feature work except for correctness fixes to the existing private
  diagnostic helper.
- RHF has served its current purpose as a Hamiltonian validator: it forced the
  density convention explicit, exercised the H1/J data, and exposed the SCF
  convergence-gating issue.
- The next priority is the Hamiltonian constructor payload boundary, not a more
  capable RHF surface.

Boundary:

- Serious HF/DMRG belongs to `codexhome/work/hfdmrg`.
- Downstream Cr2 testing belongs to the CR2 agent once GaussletBases exposes a
  clean Hamiltonian line.
- GaussletBases should focus on building clean final-basis/Hamiltonian payloads
  with explicit conventions and handoff data.

Next pass:

- No-edit Hamiltonian constructor payload audit.
- Identify what the PQS route should produce as a route-owned Ham payload:
  final basis metadata, H1 matrix, electron-electron representation, center and
  nuclear metadata, density/pair-factor convention labels, and downstream
  nonclaim/export boundaries.
- Do not implement export, public API, CR2 handoff, or route-default behavior
  in the audit.

Validation/status:

- Pass 129 made no edits and reported clean status.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: RHF route adoption is deferred; the next direction is a Ham
  constructor boundary.
- quarantined: RHF remains private route-smoke/control diagnostic only.
- not deleted because: the existing private RHF helper remains useful as a
  validator for Hamiltonian construction.
- exact remaining caller/blocker: no route-owned Hamiltonian constructor payload
  is defined yet; existing H1/J/RHF diagnostics are not the product surface.

-- repo-manager@macmini
