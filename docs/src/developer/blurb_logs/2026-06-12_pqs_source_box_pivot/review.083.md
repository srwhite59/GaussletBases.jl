Review 083: accepted.

The H1 payload return shape is now appropriately narrow. It exposes:

- `final_kinetic`;
- `final_nuclear_by_center`;
- `final_hamiltonian`;
- `h1`;
- `summary`;
- compact metadata/nonclaims.

The broad consumed/intermediate fields are gone from the returned payload:

- `source_plan`;
- `final_basis`;
- `support_kinetic`;
- `support_nuclear_by_center`.

The axis-layer/origin-factor nuclear convention check remains separate from
the payload and builds its support matrices in the test's convention-check
section. That keeps the active H1 payload final-level while preserving the
focused support-factor convention check.

Manager validation:

- `julia --project=. -e 't = @elapsed include("test/nested/pqs_direct_retained_final_h1_runtests.jl"); println("elapsed_s=", t)'`
  passed, 46/46, elapsed about 6.59s.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage:

- manual H1 assembly moved out of the test into a module-owned H1 seam;
- broad payload intermediates were removed before acceptance;
- support-space helper usage remains internal to the H1 payload or isolated in
  the small nuclear factor convention check.

-- repo-manager@macmini
