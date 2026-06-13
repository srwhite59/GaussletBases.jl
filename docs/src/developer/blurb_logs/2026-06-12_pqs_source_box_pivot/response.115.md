Pass 115 response - compact PQS RHF damping sweep

Scope:
- Local ignored probe only.
- No tracked source/test/doc edits.
- No production damping/mixing, route wiring, report fields, public API, exports, artifacts, GTO, IDA/MWG, or fixture promotion.

Local artifacts:
- `tmp/work/pqs_compact_route_owned_rhf_damping_sweep.jl`
- `tmp/work/pqs_compact_route_owned_rhf_damping_sweep.tsv`
- `tmp/work/pqs_compact_route_owned_rhf_damping_sweep_summary.txt`

Command:
- `julia --project=. tmp/work/pqs_compact_route_owned_rhf_damping_sweep.jl`

Elapsed:
- `109.971126833` seconds, measured with Julia-level `@elapsed` inside the probe.
- Cold route build/runtime exceeded 60 seconds as expected and was reported before launch.

Sweep table:

```text
alpha  converged  iterations  last_total_energy    last_density_change      last_energy_change       density_trace        monotone_energy  density_2cycle  energy_2cycle
1.0    false      100         -10.032185109282882  4.6790049683054225e-05  1.9429320197161815e-07  4                    true             false           false
0.75   false      100         -10.032180416340417  3.4635898012613797e-05  1.3860344694194282e-07  4.0000000000000009   true             false           false
0.5    false      100         -10.032175952536253  2.2727970654667917e-05  8.7695255501785141e-08  4                    true             false           false
0.25   false      100         -10.032171718404282  1.1155951330749092e-05  4.1518049442856864e-08  4.0000000000000009   true             false           false
0.1    false      100         -10.032149287715004  5.3185693829616021e-06  2.2378837272185592e-06  4.0000000000000027   true             false           false
```

Best alpha:
- No alpha converged to `density_atol = 1.0e-8`, `energy_atol = 1.0e-10` within 100 iterations.
- `alpha = 0.1` gave the smallest final density change, `5.3185693829616021e-06`, but its final energy change was worse than the mid-range alphas.
- `alpha = 0.25` gave the best balanced residuals among this sweep, but still missed the density tolerance by about three orders of magnitude.

Diagnosis:
- Simple density damping improves the density residual monotonically as alpha decreases, but does not solve convergence.
- Total energy is monotone decreasing for all alpha values.
- No obvious density or energy two-cycle signal was detected.
- Density trace stayed at 4 within roundoff, so the local damping rule did not materially drift electron count.

Recommended next pass:
- Do not add production/private damping controls yet.
- Run a residual/Fock convention audit next:
  - compare density residual definitions before and after mixing;
  - check whether the Fock build should use spin-summed density, orbital density, or an additional factor convention in the Coulomb/exchange term;
  - inspect whether energy bookkeeping is variational for the same density used to build the Fock.
- If that audit finds no convention issue, then a longer local probe with smaller alpha or acceleration can be considered, but damping alone is not yet a justified implementation seam.

Validation/status:
- `git status --short --branch`

```text
## main...origin/main
```

Ignored artifact status:

```text
!! tmp/work/pqs_compact_route_owned_rhf_damping_sweep.jl
!! tmp/work/pqs_compact_route_owned_rhf_damping_sweep.tsv
!! tmp/work/pqs_compact_route_owned_rhf_damping_sweep_summary.txt
```

Deletion/shrinkage report:
- deleted: none
- simplified: none; this was a local damping sweep only.
- quarantined: RHF SCF remains private route-smoke diagnostic output.
- not deleted because: ignored sweep artifacts are useful local evidence for the residual/convention decision.
- exact remaining caller/blocker: compact route-owned private RHF SCF remains blocked by nonconvergence; simple density damping did not converge, so the remaining blocker is a residual/Fock convention audit before adding any SCF-control seam.

-- repo-doer@macmini
