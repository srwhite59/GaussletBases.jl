Purpose:

Reduce the Be S+P decomposed/final-basis runtime by attacking the measured
dominant phase:

```julia
route_global_mixed_gto_blocks_from_decomposed_units(...)
```

Do not guess. First attribute time inside `mixed_gto_blocks`, then remove the
clearest repeated work.

Why now:

The Be S+P final-basis RHF path now works and matches the old nested/QW oracle:

- new RHF total `-14.574514244574662`
- old oracle total `-14.574514244574694`
- delta `3.2e-14 Ha`

The full run took about `357.4` seconds. The dominant phase is
`mixed_gto_blocks` at about `188.6` seconds. RHF is only about `2.6` seconds,
and final density-density assembly is about `8.4` seconds, so the next target
is not the solver or final assembly.

Known likely repeated work:

`route_global_mixed_gto_blocks_from_decomposed_units(...)` loops over retained
units. In the current code, the per-unit path can rebuild provider-level
CPB/GTO local blocks and GTO/GTO self blocks repeatedly. GTO/GTO self overlap,
kinetic, moment, and nuclear blocks depend on the supplement and centers, not
on the retained unit, so they should not be recomputed for every unit.

Exact task:

1. Add narrow subphase attribution inside the mixed-GTO route-global phase.
   Attribute at least:
   - per-unit provider/local block construction;
   - unit coefficient/support coefficient construction;
   - retained contraction;
   - row placement / coverage bookkeeping;
   - GTO/GTO self-block construction or reuse.

2. Rerun the Be S+P probe enough to identify the dominant mixed-GTO subphase.
   Use the existing fixture:
   - Be, `Z = 4`, q/ns `5 / 5`;
   - Be `cc-pV5Z`, `lmax = 1`;
   - authorized GaussletModules `BasisSets` path;
   - `build_density_density = true`.

3. If the subphase timing confirms repeated GTO/GTO self-block construction,
   hoist/cache those self blocks once per route-global mixed-GTO call and
   reuse them across retained units.

4. If the dominant cost is instead mixed CPB/GTO row construction, do not make
   a broad rewrite in this pass. Report the measured subphase and propose the
   smallest next computational replacement, such as a factorized mixed-GTO
   projection path, with enough detail for review.

5. Preserve all physics:
   - final dimension `636`;
   - RHF total matching old oracle within tight numerical tolerance;
   - no fallback flags changing from false to true;
   - raw GTO density-density still not accepted as final operator data;
   - no generalized final-basis solve.

Trust boundary:

Keep changes private/internal. Do not add public APIs, exports, route defaults,
PQS, ECP, high-l Be, Be2, Cr, H2, full-parent CPB fallback, direct Cartesian
fallback, ordinary Cartesian IDA fallback, raw GTO final density-density, or
generalized final-basis solves.

Do not add tests by default. Use `tmp/work` probes for timing/validation. If a
small production change is made, validate it through the Be probe and the load
check; run an existing local test only if the changed surface has a cheap
focused test.

Decision rules:

- If a narrow cache/hoist reduces mixed-GTO time without changing results, do
  it and report before/after timings.
- If subphase timing shows a different bottleneck, stop after recording it and
  publish the next proposed optimization target.
- If the optimization would require a broad framework or public route change,
  stop and report the smaller design first.
- If the optimization makes a compatibility field or local bundle path
  unnecessary for the active route, delete or shrink it in the same pass.

Artifacts:

- update `tmp/work/be_atom_sp_decomposed_final_basis_probe.jl`;
- update `tmp/work/be_atom_sp_decomposed_final_basis_summary.txt`;
- update or add a TSV for mixed-GTO subphase timings.

Validation:

- run the Be S+P probe;
- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`;
- if a cheap focused existing test covers the changed mixed-GTO path, run it;
  otherwise explain why the Be probe is the direct validation.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- files changed;
- mixed-GTO subphase timing table;
- before/after Be total time and mixed-GTO time, if optimized;
- RHF energy comparison to old oracle;
- whether repeated GTO/GTO self blocks were removed from the hot loop;
- any remaining dominant subphase;
- validation run;
- deletion/shrinkage report.
