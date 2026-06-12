Pass 040 complete.

Summary:
- Updated the ignored probe `tmp/work/pqs_final_basis_h1_probe.jl` to record one-center physical H1 diagnostics.
- Wrote the requested ignored artifact:
  - `tmp/work/pqs_one_center_physical_h1_probe_summary.txt`
- No tracked source, test, or docs files were changed in this pass.

Fixture:
- `fixture_kind = :one_center_pqs_shell_boundary_h1`
- `source_box_scope = :boundary_shell_source_modes_only`
- `source_mode_dims = (5, 5, 5)`
- `source_mode_count = 125`
- `q = 5`
- `L = 5`
- `nuclear_charge = 1.0`
- physical extent is known:
  - axis centers: `(-7.999999999999998, -1.9955875976247728, 0.0, 1.9955875976247728, 7.999999999999998)`
  - axis min/max: about `-8.0 / 8.0`
- boundary/final retained count:
  - `98 / 98`

Direct-retained path flags:
- retained overlap direct-boundary product used: `true`
- retained kinetic direct-boundary product used: `true`
- retained centered by-center nuclear direct-boundary product used: `true`
- active retained raw source one-body blocks materialized:
  - overlap: `false`
  - kinetic: `false`
  - centered by-center nuclear: `false`
- `_pqs_current_route_safe_term_matrices(...)` called: `false`
- ordinary symmetric eigensolve diagnostic used: `true`
- generalized overlap solve used: `false`

H1 result:
- exact hydrogenic H reference for `Z = 1`: `-0.5`
- direct final-basis H1 lowest eigenvalue: `-0.08171962129085239`
- exact-reference error: `0.41828037870914758`
- Interpretation: this is physically interpretable but not a good H acceptance fixture, because it is only the PQS boundary-shell source-mode final basis, not a complete one-center parent/final basis.

Shell-support oracle comparison:
- shell-support oracle H1 lowest eigenvalue: `-0.081719621290852015`
- H1 oracle delta: `3.7470027081099033e-16`
- final Hamiltonian versus shell-support oracle max error: `1.3322676295501878e-15`
- final overlap identity error: `2.2360255762214364e-14`
- Hamiltonian shape: `(98, 98)`
- Hamiltonian finite: `true`
- Hamiltonian symmetry error: `2.7755575615628914e-17`

Coarse timings:
- total elapsed: `3.002011s`
- parent axis setup: `0.053449s`
- shell layer: `1.170188s`
- shell realization plan: `0.034579s`
- final basis: `0.075054s`
- direct retained overlap block: `0.063004s`
- direct retained kinetic block: `0.076742s`
- direct retained centered nuclear block: `0.404883s`
- final overlap transfer: `0.017613s`
- final kinetic transfer: `0.000054s`
- final nuclear transfer: `0.022547s`
- final Hamiltonian assembly: `0.105352s`
- shell-support oracle construction: `0.068457s`
- eigensolve: `0.000566s`

Validation:
- `julia --project=. tmp/work/pqs_final_basis_h1_probe.jl` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage report:
- No old coverage or code was removed in this pass.
- This was probe-only because the current physical diagnostic is boundary-shell scoped and not yet a physical acceptance gate.
- The immediate next blocker/target before IDA/RHF is not the retained one-body route: it is defining the final IDA density-weight source for the PQS final basis and applying those weights at the final density-interaction boundary after projection/Lowdin.
- A separate physical acceptance target should use a complete one-center PQS/final basis rather than only the retained boundary shell if H/He+ energy quality is meant to be judged.

-- repo-doer@macmini
