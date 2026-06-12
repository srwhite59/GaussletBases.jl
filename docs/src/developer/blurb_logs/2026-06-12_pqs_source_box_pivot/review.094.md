Accepted with a carry-forward caveat.

The implementation adds the intended route-owned region/source-plan seam:
`cartesian_assembly(...)` now builds a private complete core/shell source-plan
payload from the selected terminal shellification/lowering plans and the
parent axis bundle, then passes the region plan, source plan, and Coulomb
expansion into the H1/J diagnostic helper. The H1/J missing-input list now
correctly drops region plan, source plan, and Coulomb expansion.

The active path remains blocked on:

- `pqs_multilayer_complete_core_shell_final_basis`;
- `pqs_multilayer_complete_core_shell_h1_payload`;
- axis weights;
- raw pair numerator terms.

The explicit-box bridge was not used as authority, and no broad tests or new
permanent tests were added. Load and diff checks passed.

Caveat:

The source-plan payload is currently a local assembly-phase object. Because
H1/J is still blocked, the blocked H1/J payload does not retain dense region or
source-plan objects in the final assembly/report result. That is good for
report and compile pressure, but the next pass should be careful: either use
the local payload immediately to build final basis/H1 inside the same assembly
seam, or introduce only a compact non-dense carrier/summary. Do not add more
report fields just to prove the source plan exists.

Second caveat:

The helper currently obtains `coulomb_gaussian_expansion(doacc = false)` inside
the assembly helper. This matches the current diagnostic path, but longer term
the Coulomb expansion should be a named route/parent input rather than an
implicit assembly choice before fixture policy is promoted.

Validation accepted:

- doer's compact dry-run smoke showing region/source/Coulomb no longer missing;
- `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- `git diff --check`.

Deletion/shrinkage review:

- Nothing became removable yet; this pass fills an input seam but still leaves
  final-basis/H1/J blocked.
- No tests were added, which is correct.
- The H1 gate was not expanded.
- Fixture-local complete core/shell H1 construction remains an oracle/probe
  surface until the driver owns final-basis and H1 payload construction.

Next pass should target final-basis/H1 payload construction from the local
source-plan payload, while preserving the test diet and avoiding RHF.

-- repo-manager@macmini
