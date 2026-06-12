Review 082: needs correction before commit.

The helper is functionally on target and validation passed, but the returned
payload is broader than the requested compact H1 assembly product. It currently
returns consumed/intermediate objects:

- `source_plan = plan`;
- `final_basis`;
- `support_kinetic`;
- `support_nuclear_by_center`.

This is exactly the kind of broad route-report shape we have been trying to
avoid. The H1 payload should expose the final-level H1 assembly products and a
compact summary, not all consumed inputs and support-space intermediates.

Keep:

- `final_kinetic`;
- `final_nuclear_by_center`;
- `final_hamiltonian`;
- `h1`;
- `summary`;
- compact metadata/nonclaims.

Remove from the returned payload unless there is a concrete live caller:

- `source_plan`;
- `final_basis`;
- `support_kinetic`;
- `support_nuclear_by_center`.

The axis-layer/origin-factor nuclear convention check may stay in the test as a
separate small support-helper comparison. It should not require broadening the
H1 payload. If the test needs the explicit support nuclear matrix for that
comparison, build it in the convention-check section rather than returning it
from the payload.

Validation observed before correction:

- focused H1 gate passed, 46/46;
- load check passed;
- `git diff --check` passed.

Do not commit pass 082 until the payload return shape is trimmed.

-- repo-manager@macmini
