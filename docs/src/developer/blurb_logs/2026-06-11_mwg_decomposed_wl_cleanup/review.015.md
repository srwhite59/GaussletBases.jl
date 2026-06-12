Review result:

Accepted as the intended replacement of the measured mixed-GTO hot path for the
active one-center atomic Be S+P route.

Main result:

- status: `:materialized_final_basis_be_sp_rhf_probe`
- route status:
  `:materialized_decomposed_atom_gto_final_basis_density_density_route`
- final dimension: `636`
- retained gausslet / supplement counts: `615 / 21`
- final overlap rank: `636`
- final overlap identity error: `1.4233342578506758e-10`
- final density-density symmetry error: `0.0`
- RHF iterations: `27`
- new RHF total: `-14.574514244574639`
- old nested/QW oracle total: `-14.574514244574694`
- RHF total delta from old oracle: `5.5067062021407764e-14`

Timing result:

The factorized retained-basis mixed-GTO path removed the measured per-unit
provider/local hot loop from the active Be route.

```text
previous route elapsed              339.913994709s
new route elapsed                   171.710878583s
previous total elapsed              342.53537725s
new total elapsed                   174.253094s
previous mixed_gto_blocks           177.196382542s
new mixed_gto_blocks                9.177082084s
previous per-unit provider local    168.066350504s
new factorized_projection_total     0.338807417s
```

The remaining top-level time is now outside the mixed-GTO hot path:

```text
residual_moment_matrices      35.501277s
electron_nuclear_by_center    28.503277125s
overlap                       20.965129041s
gausslet_density_density      19.720307167s
kinetic                       13.06034275s
mixed_gto_blocks              9.177082084s
final_basis_density_density   8.165179s
```

Interpretation:

The old slow path was not numerical quadrature over a real-space grid. It was
repeated analytic Gaussian/GTO local-block construction through CPB-local
provider machinery. The new route uses the existing decomposed WL factorized
retained-basis sidecar and old-QW-shaped parent-axis/GTO cross tables, then
projects directly into retained gausslet rows. It is currently restricted to
the one-center origin atomic `UnitPairIndexTable` case, with unsupported cases
falling back to the existing provider/local path.

Deletion/shrinkage review:

The active Be S+P hot path no longer uses the `131`-unit per-unit mixed
provider/local block materialization loop. That loop remains as compatibility
and reference fallback for non-`UnitPairIndexTable`, multi-center, or off-origin
cases. No tests were added, which is appropriate here because the Be probe is
the live route validation and another metadata test would duplicate it.

Remaining stale/duplicate surface:

The fallback provider/local mixed-GTO route should become explicit
reference/compatibility machinery if the factorized path is generalized beyond
one-center origin atomic cases. Before further algorithm work, the remaining
`174` second Be probe should be split into cold compilation versus warm runtime,
because the current top-level phases are likely a mix of first-call compilation
and real construction cost.

Validation reviewed:

- `julia --project=. tmp/work/be_atom_sp_decomposed_final_basis_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- artifact: `tmp/work/be_atom_sp_decomposed_final_basis_summary.txt`
- artifact: `tmp/work/be_atom_sp_decomposed_final_basis_phase_timings.tsv`
- artifact: `tmp/work/be_atom_sp_mixed_gto_subphase_timings.tsv`

Next target:

Run a warm/cold timing attribution pass for the full Be S+P route after the
mixed-GTO replacement. Do not start a new kernel rewrite until the remaining
top-level phases are separated into cold compilation, cached warm runtime, and
actual algorithmic cost.
