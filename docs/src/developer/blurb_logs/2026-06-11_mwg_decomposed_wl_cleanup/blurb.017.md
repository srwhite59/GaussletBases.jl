Purpose:

Add a narrow production precompile workload for the atom+GTO final-basis route
shape, using the warm/cold Be S+P timing audit as justification.

Why now:

After the factorized mixed-GTO replacement, the Be S+P route is no longer slow
when warm:

```text
label          route elapsed     RHF elapsed       total elapsed
cold           169.936205334s    2.540070625s     172.476275959s
warm measured    0.523052584s    1.880134542s       2.403187126s
```

The physics is unchanged:

- warm RHF total: `-14.574514244574639`
- old nested/QW oracle total: `-14.574514244574694`
- delta: `5.5e-14 Ha`

The remaining long cold phase timings are mostly compilation. Do not start
another algorithmic rewrite from the cold table.

Current precompile state:

`src/precompile_workloads.jl` only precompiles the small decomposed WL one-body
route: overlap, kinetic, and electron-nuclear by-center. It does not cover the
atom+GTO final-basis route shape, factorized mixed-GTO blocks, residual moment
matrices, gausslet density-density, or final-basis density-density.

Exact task:

1. Inspect `src/precompile_workloads.jl` and the private atom+GTO final-basis
   seam:

   ```julia
   _white_lindsey_decomposed_atom_gto_final_basis_route(...)
   ```

2. Add the smallest precompile workload that exercises the existing
   atom+GTO final-basis route surfaces responsible for the cold Be latency:
   - shellification-backed decomposed WL overlap/kinetic/electron-nuclear for
     the atom+GTO route shape;
   - factorized mixed-GTO block construction;
   - residual moment matrix construction;
   - gausslet density-density;
   - final-basis density-density.

3. Do not depend on the user-local GaussletModules `BasisSets` path. Use or add
   a tiny repo-local/synthetic supplement fixture that exercises the same code
   paths. It does not need to be chemically complete; it is a compile workload,
   not a physics acceptance test.

4. Keep the workload narrow:
   - no RHF solve;
   - no scientific assertions;
   - no public API, export, route default, or artifact path;
   - no PQS, ECP, high-l Be, Be2, H2, Cr, full-parent CPB fallback, direct
     Cartesian fallback, ordinary Cartesian IDA fallback, raw GTO final
     density-density, or generalized final-basis solve.

5. Measure before/after enough to decide whether the workload earns its
   carrying cost:
   - package precompile cost after the workload edit;
   - fresh-process `using GaussletBases` after cache is built;
   - fresh-process Be S+P route or warm/cold probe timing after the workload,
     if practical;
   - same-process warm timing should remain essentially unchanged.

Decision rules:

- If a small synthetic supplement can compile the route shape cleanly and
  materially reduces fresh-process Be route latency, keep it.
- If the workload becomes large, depends on user-local data, or starts encoding
  acceptance behavior, stop and report the smaller design instead.
- If precompile cost is too high relative to cold-latency benefit, revert or
  leave the workload out and report the evidence.
- If this makes old ad hoc timing/probe code unnecessary, delete or quarantine
  it in the same pass; otherwise explain why nothing became removable.

Validation:

- run `julia --project=. -e 'using GaussletBases; println("load ok")'`;
- run `git diff --check`;
- run the Be S+P timing probe or route probe if practical after precompile;
- no tests by default unless a source edit needs a cheap focused test and that
  test protects a live contract.

Deletion/shrinkage report required:

- what old code, test, metadata, or compatibility path became unnecessary;
- what was deleted or simplified;
- if nothing was deleted, why no existing surface was made obsolete;
- whether any new test replaces/shrinks older coverage or is genuinely new
  live-contract coverage;
- any remaining stale or duplicate surfaces to retire next.

Report back:

- files changed;
- exact synthetic/repo-local supplement fixture used;
- route surfaces exercised by the workload;
- package precompile cost;
- fresh-process load and Be route timing before/after if measured;
- whether the workload earns carrying cost;
- validation run;
- deletion/shrinkage report.
