Purpose:

Measure the effect of the inventory-summary shape cleanup before doing another
struct or precompile pass.

Why now:

The hot decomposed WL inventory result no longer encodes inventory size in
`unit_keys`, `unit_summaries`, or `pair_summaries`:

```text
unit_keys       Vector{Symbol}
unit_summaries  Vector{NamedTuple{...}}
pair_summaries  compact count/status NamedTuple
```

This removes the concrete side-7 versus side-15 type mismatch that the
specialization audit found. But the behavior validation run was not a clean
fresh-process timing comparison. Do not infer compile-time improvement until
it is measured.

Exact task:

1. Rerun the atom+GTO specialization-shape audit:

   ```text
   tmp/work/atom_gto_specialization_shape_audit.jl
   ```

   Confirm that `unit_keys`, `unit_summaries`, and `pair_summaries` no longer
   differ by inventory-size-specialized tuple types. Report any remaining
   shape differences.

2. Rerun the Be S+P warm/cold timing probe:

   ```text
   tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_timing_probe.jl
   ```

   Compare against the previous post-precompile numbers:

   ```text
   cold route elapsed      159.006798084s
   cold total elapsed      161.9370205s
   warm route elapsed        0.498273209s
   warm total elapsed        2.486322293s
   ```

3. Decide the next target from evidence:
   - If cold route time materially improves, record the result and identify
     the next largest remaining specialization source.
   - If cold route time does not materially improve, do not revert the cleanup;
     it still removed bad type shape. Instead rank the remaining candidates:
     `pair_keys` small/large shape split, parent-axis metadata,
     giant atom+GTO route result, factorized sidecar, residual moment bundles,
     or timing closures.
   - If warm route time regresses, stop and report before further changes.

4. Do not implement a broad refactor in this pass. A tiny local cleanup is
   allowed only if the measurement directly identifies it and the behavior
   validation remains cheap.

Trust boundary:

No public APIs, exports, route defaults, new acceptance tests by default, PQS,
ECP, Be2, high-l Be, H2, Cr, full-parent CPB fallback, direct Cartesian
fallback, ordinary Cartesian IDA fallback, raw GTO final density-density, or
generalized final-basis solve.

Test policy:

Do not add tests. Use `tmp/work` probes and existing focused tests only.

Validation:

- run the specialization-shape audit;
- run the Be S+P warm/cold timing probe;
- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`;
- if any source changes are made, rerun the focused inventory test or other
  directly affected cheap test.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- files changed, if any;
- before/after specialization-shape audit table;
- before/after Be warm/cold timing table;
- whether the inventory cleanup improved cold route latency;
- ranked next compile-pressure candidates;
- validation run;
- deletion/shrinkage report.
