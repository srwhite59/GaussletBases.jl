Accepted.

Pass 027 proved the real projected-q-shell final-basis H1 path as a probe:

```text
raw source plan + PQS retained rule
-> final shell basis
-> final overlap/kinetic
-> final by-center nuclear
-> Hamiltonian-stage charge application/summing
-> ordinary symmetric eigensolve
```

Key results:

```text
final overlap identity error = 2.2360255762214364e-14
final Hamiltonian vs shell-support oracle error = 1.3322676295501878e-15
H1 lowest eigenvalue = -0.08171962129085239
oracle H1 lowest eigenvalue = -0.081719621290852015
H1 delta = 3.7470027081099033e-16
```

The probe used an ordinary symmetric eigensolve, did not use generalized
overlap solve logic, did not call `_pqs_current_route_safe_term_matrices(...)`,
and did not claim IDA, density-density, RHF, driver, export, or artifact
behavior.

Manager validation:

- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- reviewed `tmp/work/pqs_final_basis_h1_probe_output.txt`

Next target:

Do not immediately add a broad permanent test. First audit and shrink stale PQS
H1/readiness/oracle surfaces that the explicit final-basis path has superseded.
If a permanent gate is added later, it should be compact and physics-first.

-- repo-manager@macmini
