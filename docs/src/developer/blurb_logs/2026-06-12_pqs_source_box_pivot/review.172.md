Manager review for pass 172

Accepted.

The artifact now declares the HF-facing conventions CR2 asked for:

- PQS final basis uses `:orthonormal_identity_by_contract`, stores no overlap
  matrix, and keeps identity defect at `0.0`.
- WL final basis uses `:orthonormal_identity_with_diagnostic_matrix`, keeps its
  stored overlap as diagnostic route data, and reports identity defect around
  `4e-14`.
- Both routes expose the common
  `two_body/interaction_matrix_representation_kind =
  :final_basis_density_density_matrix`.
- Both routes carry
  `:missing_reviewed_density_density_hf_fock_energy_convention`, and
  solver/export readiness remains false.

The line budget is satisfied across tracked `src`, `test`, and generator code:
40 added, 194 deleted, net -154. Deleting the standalone reviewed-overlap
placement-plan metadata test is acceptable because related placement behavior
is still exercised by the broader overlap placement/facts and route tests.

Next GaussletBases work should be generator performance/size reduction before
another CR2 handoff.

-- repo-manager@macmini
