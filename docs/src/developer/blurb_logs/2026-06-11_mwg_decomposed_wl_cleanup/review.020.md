Accepted as a measurement-only pass.

The inventory-summary cleanup did what it was supposed to do. The targeted
hot inventory report fields no longer encode the retained-unit count in their
concrete types:

```text
unit_keys        Vector{Symbol}
unit_summaries   Vector{NamedTuple{...}}
pair_summaries   compact count/status NamedTuple
```

The fresh-process Be S+P route timing moved from about 159s cold route time to
about 31s, while warm route time stayed essentially unchanged at about 0.49s.
That does not prove every compile problem is solved, but it is strong evidence
that the tuple-shaped inventory summaries were real specialization debt, not
just cosmetic report clutter.

The physics check remains clean:

```text
Be S+P RHF total        -14.574514244574639
old nested/QW oracle    -14.574514244574694
delta                    5.51e-14 Ha
final dimension          636
retained gausslets       615
units / pairs            131 / 8646
```

No source, tests, or docs changed in this pass beyond the response log. No
tests were added, and the doer used the existing `tmp/work` probes as requested.

Remaining cold pressure is now clearly elsewhere. The largest measured cold
phase is `residual_moment_matrices`: about 6.9s cold and about 0.003s warm.
The atom+GTO route currently builds that phase by calling six separate
decomposed WL moment builders, then passes six report-shaped results into
`route_global_combined_gto_residual_moment_matrices(...)`.

Next target should be narrow: reduce residual-moment route shape/compile
pressure and keep compute payloads separate from audit reports. This is not a
license to replace every `NamedTuple` with a struct. Use structs only for stable
concepts that are passed repeatedly as real compute objects; keep compact
summaries as summaries.

Validation reported by doer:

- `julia --project=. tmp/work/atom_gto_specialization_shape_audit.jl`
- `julia --project=. tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_timing_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

-- repo-manager@macmini
