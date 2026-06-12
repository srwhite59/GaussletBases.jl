Accepted as an audit-only pass.

The result is useful because it stopped a likely over-optimization. The direct
subphase probe did not reproduce the full route's cold
`decomposed_wl_one_electron_matrix_set` cost:

```text
direct cold factorized basis lookup/extract   0.038492417s
direct cold overlap helper                    0.047759s
direct cold kinetic helper                    0.000734333s
direct cold electron-nuclear helper           0.359268209s
direct cold full matrix set                   0.2198505s

full-route cold one-electron matrix set       4.365419666s
```

That means the obvious inner helper rewrite is not justified by current
evidence. The remaining cold cost is more likely route-context specialization,
TimeG/timing wrapper shape, metadata/report payload shape, or first-call
compilation caused by the full atom+GTO call graph.

The rerun route timing stayed stable:

```text
cold route elapsed   24.996941042s
warm route elapsed    0.502428875s
Be S+P RHF total     -14.574514244574639
old nested/QW oracle -14.574514244574694
```

No source or test changes were made. No tests were added. The new artifact is
only an ignored `tmp/work` timing probe, which is the right surface for this
kind of exploratory audit.

Next target:

Do not rewrite the inner one-electron helpers yet. The next pass should run a
full-route attribution probe with controlled variants:

- phase timing sinks enabled versus disabled;
- `build_density_density = false` versus full density-density route;
- minimal metadata versus the current metadata shape;
- one-electron-only route result construction versus full report payload;
- current fresh-process cold route and same-process warm route.

The goal is to identify whether the remaining cold pressure is real numerical
kernel compilation, instrumentation/closure pressure, giant route result shape,
metadata shape, or density-density/final-basis stages.

-- repo-manager@macmini
