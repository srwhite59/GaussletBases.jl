Purpose:

Audit the PQS electron-nuclear source-factor seam before implementation.

Why now:

The retained PQS source-mode route now has one-unit overlap and kinetic
matrices. A meaningful H/He+ one-electron diagnostic next needs nuclear
attraction:

```text
V_nuc = -Z * sum_t c_t Gx_t * Gy_t * Gz_t
```

This is the first non-safe PQS one-body term in the pivot. It should reuse the
repo's existing Gaussian Coulomb expansion and 1D factor conventions. Do not
guess the convention from the overlap/kinetic helper shape.

Exact task:

Do a read-only audit plus optional `tmp/work` probe. Do not change production
source or tests in this pass unless the audit finds a trivial documentation
correction.

Find and report:

- existing Coulomb/Gaussian expansion helpers used by WL electron-nuclear and
  old nested/QW code;
- how those helpers represent coefficients, exponents, centers, charges, and
  signs;
- whether source-space 1D Gaussian factor matrices already exist for PQS-like
  source modes, or what minimal helper would be needed;
- where the PQS source-space electron-nuclear block should live if implemented;
- how retained source-mode contraction should apply after raw source-space
  nuclear block construction;
- exact nonclaims to preserve: no shell realization, no Lowdin, no IDA, no
  Hamiltonian/driver adoption;
- first physics diagnostic after implementation, likely one-center H or He+
  H1, not RHF.

Code surfaces to inspect:

```text
src/cartesian_pair_block_materialization/white_lindsey_electron_nuclear.jl
src/cartesian_pair_block_materialization/pqs_source_safe_terms.jl
src/cartesian_pair_block_materialization/pqs_source_one_body.jl
src/cartesian_nested_faces.jl
src/cartesian_nested_diatomic.jl
docs/src/developer/pqs_source_box_operator_framework.md
docs/src/developer/numerical_contracts.md
```

Use `rg` to find relevant helpers rather than reading whole large files.

Trust boundary:

No implementation of electron-nuclear, density-density, IDA, Hamiltonian,
driver routes, shell realization, Lowdin cleanup, exports, artifacts,
full-parent fallback, or direct Cartesian fallback in this pass.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path appears likely to become
  unnecessary once PQS retained source-mode nuclear exists;
- what, if anything, was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete by this
  audit;
- whether any new artifact was added and why it earned its carrying cost;
- any stale or duplicate surfaces to retire after implementation.

Validation:

- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- no Julia tests unless production code or tests are changed.

Report back:

- write `.agent_handoffs/response.006.md.tmp`, then atomically rename to
  `.agent_handoffs/response.006.md`;
- also write the curated copy to
  `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.006.md`;
- include inspected surfaces;
- include exact recommended implementation seam;
- include blockers or convention risks;
- include validation run;
- include deletion/shrinkage report;
- sign `-- repo-doer@macmini`.

-- repo-manager@macmini
