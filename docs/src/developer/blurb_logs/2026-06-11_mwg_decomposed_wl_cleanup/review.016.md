Review result:

Accepted as the necessary warm/cold attribution pass after the mixed-GTO
algorithmic bottleneck was removed. The result changes the next target from
"rewrite another route phase" to "compile the route shape intentionally."

Main timing result:

```text
label          route elapsed     RHF elapsed       total elapsed
cold           169.936205334s    2.540070625s     172.476275959s
warmup           0.549007292s    1.890630792s       2.439638084s
warm measured    0.523052584s    1.880134542s       2.403187126s
```

The Be S+P physics stayed fixed:

- cold RHF total: `-14.574514244574639`
- warm RHF total: `-14.574514244574639`
- old nested/QW oracle total: `-14.574514244574694`
- warm delta from oracle: `5.5067062021407764e-14 Ha`
- final dimension: `636`
- retained gausslet dimension: `615`
- units / pairs: `131 / 8646`

Cold versus warm phase attribution:

```text
phase                         cold              warm
residual_moment_matrices      36.155260708s     0.002820542s
electron_nuclear_by_center    28.065223s        0.049385625s
overlap                       20.711751125s     0.005225084s
gausslet_density_density      19.2412885s       0.004252083s
kinetic                       12.651555542s     0.090058958s
mixed_gto_blocks               9.199361667s     0.216573708s
final_basis_density_density    8.028131375s     0.014046792s
```

Interpretation:

The remaining long cold route phases are almost entirely compilation. The warm
route construction is about half a second. The measured warm total is dominated
by the probe-local RHF solve at about `1.88` seconds, not by operator
construction. Among warm route phases, mixed GTO remains the largest at about
`0.217` seconds, but that is no longer large enough to justify another
algorithmic rewrite before precompile work.

Deletion/shrinkage review:

No production code, tests, metadata, or compatibility paths became unnecessary
in this timing-only pass. No test was added; the timing probe lives under
`tmp/work`, which is appropriate for exploratory attribution. The next
constructive task is a precompile workload extension, not deletion.

Validation reviewed:

- `julia --project=. tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_timing_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- artifact:
  `tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_timing_summary.txt`
- artifact:
  `tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_phase_timings.tsv`
- artifact:
  `tmp/work/be_atom_sp_decomposed_final_basis_warm_cold_mixed_gto_subphase_timings.tsv`

Next target:

Add a narrow production precompile workload for the atom+GTO final-basis route
shape. It should use a repo-local or synthetic supplement fixture rather than
the user-local GaussletModules `BasisSets` path, should not solve RHF or encode
acceptance logic, and should be measured for package precompile cost versus
fresh-process route latency.
