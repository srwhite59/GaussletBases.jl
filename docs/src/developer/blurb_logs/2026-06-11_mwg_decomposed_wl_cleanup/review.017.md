Review result:

Accepted as a limited precompile-workload extension. It is small, repo-local,
and exercises the intended atom+GTO final-basis route shape without turning the
Be acceptance fixture into package precompile work. The measured benefit is
modest, so this should not be followed by simply making the precompile fixture
larger.

What changed:

- `include("precompile_workloads.jl")` moved later in `src/GaussletBases.jl`
  so the atom+GTO supplement representation and old QW/GTO cross-table helpers
  are available.
- `src/precompile_workloads.jl` now has a tiny synthetic one-S Gaussian
  supplement fixture at the origin.
- The workload calls the existing private seam:
  `_white_lindsey_decomposed_atom_gto_final_basis_route(...)`.
- It does not solve RHF and does not add acceptance assertions, public APIs,
  route defaults, exports, artifacts, fallbacks, or tests.

Route surfaces exercised:

- shellification-backed decomposed WL inventory;
- route-global overlap, kinetic, and electron-nuclear by center;
- factorized mixed gausslet/GTO blocks;
- combined one-electron matrix assembly;
- final-basis projection;
- residual moment matrix construction;
- residual MWG representation;
- gausslet density-density;
- final-basis density-density.

Measured effect:

```text
package precompile after edit      48956 ms
cached fresh using GaussletBases   0.654013292s

Be S+P cold route before           169.936205334s
Be S+P cold route after            159.006798084s
Be S+P cold total before           172.476275959s
Be S+P cold total after            161.9370205s

Be S+P warm route before             0.523052584s
Be S+P warm route after              0.498273209s
```

Physics remained fixed:

- post-workload warm RHF total: `-14.574514244574639`
- old nested/QW oracle total: `-14.574514244574694`
- delta from old oracle: `5.5067062021407764e-14 Ha`
- final dimension: `636`

Interpretation:

The workload earns limited carrying cost because it is narrow and independent
of user-local `BasisSets`, but it only reduces Be cold route latency by about
`6%`. The side-7 synthetic fixture does not fully precompile the side-15 Be
specializations. A larger precompile fixture would risk moving acceptance-scale
work into package precompile. The next target should be specialization pressure
in the route objects and staged summaries, not a bigger workload.

Deletion/shrinkage review:

No old production code, tests, metadata, or compatibility path became
unnecessary. No tests were added. The generated probe remains under `tmp/work`.
Nothing should be deleted based on this pass alone.

Validation reviewed:

- `julia --project=. tmp/work/atom_gto_precompile_synthetic_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `julia --project=. -e 't = @elapsed using GaussletBases; println("load_elapsed_s=", t); println("load ok")'`
- `julia --project=. tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_timing_probe.jl`
- `git diff --check`

Next target:

Audit cold compile specialization pressure in the atom+GTO final-basis route.
Look for large or value-specialized staged object shapes, route summaries,
retained-unit sidecars, and `NamedTuple`/tuple-heavy paths that force new
compilation for side-15 Be after the side-7 synthetic workload. Do not add a
larger precompile fixture until the specialization source is understood.
