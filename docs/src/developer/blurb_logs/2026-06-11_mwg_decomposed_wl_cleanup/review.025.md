Accepted as a density-route attribution closeout.

The probe gives enough evidence to pause the GTO compile-optimization loop
after recording the remaining debt. The final density result constructor is not
the bottleneck:

```text
density_result_only cold      0.050661875s
density_result_only warm      0.000007334s
final density result type     38 properties, type-name length 1761
```

The old residual-MWG kernel itself is also not the main cold source:

```text
residual_mwg_kernel_only cold 0.264626459s
residual_mwg_kernel_only warm 0.013451375s
```

The remaining measured pressure is earlier in the density prerequisites:

```text
residual_moment_matrices_only cold 7.178579542s
mixed_moment_blocks_only warm      0.368781584s
mixed moment block type length     21712
gausslet_density_only cold         2.308185875s
```

So the next implementation target, if we return to this, is not a broad
atom+GTO route-result rewrite. It is the moment-capable mixed GTO block payload
and residual moment prerequisite path. The full atom+GTO route result is still
large and awkward, but this probe does not make it the next measured bottleneck.

No source or test changes were made. No tests were added. The new artifacts are
ignored `tmp/work` probes, which is appropriate for this audit stage.

Recommended pause point:

- Record this as known GTO compile debt.
- Do not continue polishing unless it blocks the next physics target.
- Pivot to the PQS source-box-first plan after setting up the baton-loop
  mechanics.

Validation reported by doer:

- `julia --project=. tmp/work/be_atom_sp_density_route_compile_attribution_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

-- repo-manager@macmini
