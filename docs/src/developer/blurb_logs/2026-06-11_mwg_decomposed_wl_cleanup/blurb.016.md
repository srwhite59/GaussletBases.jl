Purpose:

Attribute the remaining Be S+P decomposed/final-basis runtime after the
factorized mixed-GTO replacement. Separate cold compilation from warm runtime
before starting another algorithmic rewrite.

Why now:

The mixed-GTO replacement succeeded:

- Be S+P RHF total: `-14.574514244574639`
- old nested/QW oracle total: `-14.574514244574694`
- delta: `5.5e-14 Ha`
- total elapsed: about `342.5s -> 174.3s`
- `mixed_gto_blocks`: about `177.2s -> 9.18s`
- factorized mixed-GTO projection subphase: about `0.34s`

The previous bottleneck is gone. The remaining top-level timings from the cold
probe are now:

```text
residual_moment_matrices      35.501277s
electron_nuclear_by_center    28.503277125s
overlap                       20.965129041s
gausslet_density_density      19.720307167s
kinetic                       13.06034275s
mixed_gto_blocks              9.177082084s
final_basis_density_density   8.165179s
```

These numbers are probably a mixture of first-call compilation and real
runtime. Do not optimize from cold timings alone.

Exact task:

1. Build or update a `tmp/work` warm/cold Be S+P timing probe for the existing
   private decomposed atom+supplement seam. Use the same fixture:
   - Be, `Z = 4`, q/ns `5 / 5`;
   - Be `cc-pV5Z`, `lmax = 1`;
   - authorized GaussletModules `BasisSets` path;
   - density-density/RHF enabled.

2. Run at least:
   - one cold process timing, using the current phase timing table;
   - one same-process warmup followed by one or more measured warm route calls.

3. Preserve durable artifacts under `tmp/work`, including:
   - a summary text file;
   - a phase timing TSV for cold and warm timings;
   - any subphase TSV if a remaining top-level phase needs drilldown.

4. Report which phases are mostly compilation and which remain warm runtime
   costs. Do not infer algorithmic bottlenecks from cold-only data.

5. If one warm phase remains clearly dominant, inspect that code surface enough
   to name the next replacement or precompile target. Do not implement the
   optimization in this pass unless it is a tiny obvious cleanup that does not
   change route shape.

6. If compilation dominates the remaining cold cost, propose the narrowest
   production precompile workload extension. Do not add it in this pass unless
   the evidence is clear and the workload is small.

Trust boundary:

Prefer `tmp/work` probes and documentation/log updates. Do not add public APIs,
exports, route defaults, tests by default, PQS, ECP, Be2, high-l Be, generalized
final-basis solves, raw GTO final density-density, full-parent CPB fallback,
direct Cartesian fallback, or ordinary Cartesian IDA fallback.

Test policy:

Do not add tests for this timing audit. If source code changes are limited to
tiny instrumentation, validate with the Be probe, `using GaussletBases`, and
`git diff --check`. If no source changes are needed, say so.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- files changed;
- cold versus warm timing table;
- which phases are compilation-heavy;
- which phases remain warm runtime costs;
- whether a precompile workload extension is justified;
- whether an algorithmic optimization target is justified;
- validation run;
- deletion/shrinkage report.
