Pass 130 manager review

Accepted as the right no-edit audit and the correct pivot away from RHF feature
work.

The audit identified a clean private Hamiltonian constructor boundary:

- consume the existing route-owned complete core/shell diagnostic payload;
- carry final basis, H1 payload/final Hamiltonian, density inputs, H1/J payload,
  density interaction, Coulomb/source provenance, compact summaries, convention
  labels, and nonclaim flags;
- do not derive authority from report aliases;
- do not add public API, export/artifact writing, CR2 handoff, IDA/MWG
  promotion, RHF product surface, or serious-HF claims.

Decision:

- Implement the private Ham payload seam next.
- Do not implement RHF route adoption.
- Do not expose report fields in the first pass.
- Test at the assembly/private payload level so we can inspect object-carrying
  fields without expanding `cartesian_report(...)`.

Target:

- Add `_pqs_source_box_route_driver_complete_core_shell_ham_payload(...)`.
- Add a compact `complete_core_shell_ham_payload` slot to
  `_PQSCompleteCoreShellDiagnosticRoutePayload`.
- Keep heavy objects inside the payload and expose compact summary/convention
  data there only.
- Add one focused route-smoke test that builds the one-center PQS route through
  `cartesian_assembly(...)` and inspects the private payload.

Validation/status:

- Pass 130 made no edits and reported clean status.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: downstream Ham discussion can now target one private payload seam
  instead of RHF route wiring.
- quarantined: RHF remains private validator only; report/export/public
  surfaces remain out of this pass.
- not deleted because: this was an audit; compatibility report/export surfaces
  remain untouched.
- exact remaining caller/blocker: the route owns enough H1/J objects, but no
  private Ham payload yet carries them as the Hamiltonian-constructor boundary.

-- repo-manager@macmini
