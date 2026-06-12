Purpose:

Audit cold compile specialization pressure in the atom+GTO final-basis route.
Do not add a larger precompile fixture yet.

Why now:

The first atom+GTO precompile workload is small and repo-local, but its benefit
is limited:

```text
package precompile after edit      48956 ms
cached fresh using GaussletBases   0.654013292s

Be S+P cold route before           169.936205334s
Be S+P cold route after            159.006798084s
Be S+P cold total before           172.476275959s
Be S+P cold total after            161.9370205s

Be S+P warm route after              0.498273209s
Be S+P warm total after              2.486322293s
```

The synthetic side-7 workload does not fully precompile the side-15 Be route
specializations. Making the precompile fixture bigger may just move
acceptance-scale work into package precompile. First identify what is causing
side/shape-specific compilation.

Exact task:

1. Audit the atom+GTO final-basis route for specialization-heavy object shapes.
   Focus on:
   - route result summaries and metadata `NamedTuple`s;
   - retained-unit and factorized retained-basis sidecars;
   - unit/pair inventory shapes;
   - large tuples or value-specialized fields;
   - helper signatures that specialize on full staged object types when only a
     small summary or matrix is needed.

2. As part of the audit, identify repeated field clouds that are really stable
   module-owned concepts and should eventually become small structs. Do not
   start from "replace NamedTuples with structs everywhere." Start from nouns
   already moving through the route. Likely candidates include:
   - `AxisTriplet{T}` for recurring Cartesian `x/y/z` bundles;
   - `ParentAxisContext3D` for parent axis counts, bundles, and 1D operator
     factors;
   - `DecomposedWLFactorizedBasis3D` for the retained WL basis represented in
     parent-product coordinates;
   - `ResidualMomentMatrices` for position and x2 axis triples;
   - `AtomGTOFinalBasisMatrices` as the lean compute product consumed by RHF.

   The intended direction is:

   ```text
   hot compute path:
       compact structs for stable numerical concepts

   audit/report path:
       compact NamedTuple summaries built for probes, docs, and debugging
   ```

   Avoid creating one giant route struct that simply mirrors the current giant
   nested `NamedTuple`. That would improve neither conceptual clarity nor cold
   compile pressure.

3. Use lightweight evidence, not guesswork. Good options include:
   - static inspection with `rg` and focused file reads;
   - a small `tmp/work` trace-compile or method-instance probe if useful;
   - comparing side-7 synthetic workload shape versus side-15 Be shape;
   - timing a small alternate synthetic fixture only if it answers a specific
     specialization question.

4. Do not implement a broad refactor in this pass. If you find a tiny obvious
   cleanup, it must be local, behavior-preserving, and justified by the audit.
   Otherwise report the exact specialization source and the smallest proposed
   fix.

5. Do not add tests by default. This is a compile/specialization audit; use
   `tmp/work` probes for evidence.

6. Do not broaden `src/precompile_workloads.jl` to side-15 or Be-like
   acceptance scale in this pass. A larger workload is allowed only after the
   audit shows it is the right answer and earns carrying cost.

Trust boundary:

No public APIs, exports, route defaults, new acceptance tests, PQS, ECP, Be2,
high-l Be, H2, Cr, full-parent CPB fallback, direct Cartesian fallback,
ordinary Cartesian IDA fallback, raw GTO final density-density, or generalized
final-basis solve.

Validation:

- If no source changes: run `git diff --check` and
  `julia --project=. -e 'using GaussletBases; println("load ok")'`.
- If source changes: also rerun the relevant small probe or Be timing probe
  needed to prove behavior is unchanged.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- files inspected;
- any `tmp/work` probes/artifacts created;
- the top suspected specialization sources, ranked;
- the stable concepts that should be structs, if any, and which current
  `NamedTuple`/field-cloud surfaces they would replace;
- whether any tiny cleanup was made;
- whether a future precompile workload should change, and why;
- validation run;
- deletion/shrinkage report.
