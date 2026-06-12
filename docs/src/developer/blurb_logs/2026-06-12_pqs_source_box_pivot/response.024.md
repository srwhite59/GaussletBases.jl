Status: complete.

Created and ran the audit probe:

- `tmp/work/pqs_retained_nuclear_boundary_audit_probe.jl`
- `tmp/work/pqs_retained_nuclear_boundary_audit_probe_output.txt`

Result:

The existing CPBM retained-source centered electron-nuclear by-center block matches an independent shell-support boundary oracle on the real projected-q-shell fixture.

Fixture:

- `current_box = (1:5, 1:5, 1:5)`
- `inner_box = (2:4, 2:4, 2:4)`
- `q = 5`, `L = 5`
- source modes: `(5, 5, 5)`, count `125`
- retained/source boundary count: `98`
- shell support count from the old oracle plan: `98`

Origin center:

- center index: `1`
- location: `(0.0, 0.0, 0.0)`
- recorded charge: `2.0`
- term: `:retained_source_electron_nuclear_by_center`
- block shape: `(98, 98)`
- retained block vs `P' * V_shell_support * P` max error: `1.3877787807814457e-16`
- retained block symmetry error: `3.4694469519536142e-18`
- boundary oracle symmetry error: `6.9388939039072284e-18`
- nuclear charge recorded: `true`
- nuclear charge applied: `false`
- centers summed: `false`
- uncharged by-center convention: `true`

Off-origin center:

- center index: `2`
- location: `(0.15, -0.2, 0.25)`
- recorded charge: `3.0`
- term: `:retained_source_electron_nuclear_by_center`
- block shape: `(98, 98)`
- retained block vs `P' * V_shell_support * P` max error: `1.3877787807814457e-16`
- retained block symmetry error: `1.0408340855860843e-17`
- boundary oracle symmetry error: `1.3877787807814457e-17`
- nuclear charge recorded: `true`
- nuclear charge applied: `false`
- centers summed: `false`
- uncharged by-center convention: `true`

Nonclaims verified by the probe:

- `_pqs_current_route_safe_term_matrices(...)` was not called.
- no old fixed-block matrix authority was used.
- no H1 solve, Hamiltonian assembly, IDA, density-density, RHF, driver route, exports, or artifacts were materialized.

Answer to pass-024 questions:

- Yes, retained-source electron-nuclear by-center equals the shell-projected boundary nuclear operator on the real fixture to roundoff.
- Yes, it preserves charge-recorded/not-applied and centers-not-summed metadata.
- The off-origin center is also clean, which checks the axis-centering path.

Recommended pass-025 implementation target:

Extend the CPBM shell-final one-body boundary materialization path for the retained-source by-center nuclear term. The narrow next step is to allow `pqs_source_shell_final_one_body_from_boundary_matrix(...)` or a by-center-specific sibling to consume `:retained_source_electron_nuclear_by_center`, preserve center metadata, keep nuclear charge unapplied, and keep centers separated. Do not assemble H1 or sum/apply charges until a separate Hamiltonian assembly pass.

Validation:

- `julia --project=. tmp/work/pqs_retained_nuclear_boundary_audit_probe.jl` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.
- `git status --short` showed no tracked changes before writing this response.

Deletion/shrinkage report:

- No production code, permanent tests, metadata contracts, or compatibility paths were changed in this probe-only pass.
- Nothing was deleted because the audit only establishes whether the existing retained-source nuclear block can become an input to the shell-final boundary materialization path.
- The temporary `tmp/work` probe earned its carrying cost by replacing speculation about the nuclear boundary source with a concrete shell-support oracle comparison on the real fixture.
- No new permanent test replaces or shrinks older coverage yet.
- Remaining stale/duplicate surfaces to retire next: the retained-source/PQS H1 scratch probes and old current-route safe-term oracle vocabulary should be revisited after nuclear final-boundary materialization and explicit Hamiltonian-stage charge application exist.

-- repo-doer@macmini
