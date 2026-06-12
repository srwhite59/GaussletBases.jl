Purpose:

Attribute the remaining Be atom+GTO cold compile pressure at the full-route
level before making more source changes.

Why now:

The one-electron matrix-set subphase audit did not justify an inner helper
rewrite. Direct helper timings were much smaller than the full-route cold phase:

```text
direct cold full one-electron matrix set      0.2198505s
full-route cold one-electron matrix set       4.365419666s
```

The remaining cold cost is therefore likely caused by the full atom+GTO route
context rather than simple local matrix construction.

Exact task:

Create or update an ignored `tmp/work` developer probe that runs the current
Be S+P atom+GTO route in controlled variants and reports cold/warm timings.
Do not change production source unless the probe identifies one tiny,
well-supported cleanup.

Measure at least these variants:

1. Current route shape:
   - phase timing sinks enabled;
   - mixed GTO subphase timing enabled;
   - `build_density_density = true`;
   - current metadata shape.

2. Timing-sink-disabled route:
   - `phase_timings = nothing`;
   - `mixed_gto_subphase_timings = nothing`;
   - otherwise same fixture and `build_density_density = true`.

3. One-electron-only route:
   - `build_density_density = false`;
   - timing sinks enabled and disabled if practical.

4. Minimal metadata route:
   - use the same physical fixture, but pass the smallest metadata value that
     still exercises the route;
   - compare with the current metadata shape.

5. If practical without source changes, separately time route-result
   construction or summarize its concrete type length/property count for
   one-electron-only and full density-density results.

For each variant, report:

```text
cold route elapsed
warm route elapsed
top cold phases when phase timings are enabled
RHF energy or one-electron route status, depending on variant
fallback/final-basis flags when available
```

Interpretation required:

- If disabling timing sinks materially reduces cold time, identify whether the
  next cleanup should move timing wrappers out of the hot route or reduce
  closure specialization while preserving useful TimeG/coarse timing.
- If `build_density_density = false` removes most cold cost, focus next on
  final density-density / residual MWG route shape.
- If minimal metadata materially reduces cold time, focus next on a
  route-owned parent-axis context or compact metadata summaries.
- If route-result type/property size appears dominant, propose a lean compute
  result plus separate audit summary, but do not implement it in this pass.
- If none of these variants explains the cold cost, rank likely old-QW/GTO
  analytic kernels or precompile targets using the evidence.

Trust boundary:

Measurement first. No public API/export/default changes. No PQS, ECP, Be2, H2,
high-l Be, driver defaults, acceptance fixtures, raw GTO final density-density,
generalized final solve, full-parent CPB fallback, direct Cartesian fallback,
or ordinary Cartesian IDA fallback.

Test policy:

Do not add tests. Use ignored `tmp/work` probes. If no production source
changes are made, no Julia test is required beyond a load check.

Validation:

- run the new/updated full-route attribution probe;
- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`;
- if production source changes are made, run the smallest directly affected
  focused test and the Be warm/cold probe.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- files changed, if any;
- probe artifact paths;
- variant timing table;
- interpretation and ranked next target;
- whether any source change was made;
- validation run;
- deletion/shrinkage report.

-- repo-manager@macmini
