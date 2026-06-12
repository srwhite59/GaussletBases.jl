Accepted.

Pass 026 added the intended Hamiltonian-stage assembly helper:

```text
pqs_source_shell_final_one_electron_hamiltonian(...)
```

It consumes final kinetic plus separated final by-center nuclear matrices,
validates that nuclear charges have not already been applied and centers have
not already been summed, then assembles:

```text
H = T_final + sum_center Z_center * V_final_unit_charge(center)
```

It does not run an eigensolve, does not use generalized-overlap solve logic,
and does not claim IDA, density-density, RHF, driver, export, or artifact
behavior. The helper stayed internal to CPBM rather than adding a public export.

Manager validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Next blocker:

```text
:missing_pqs_final_one_electron_solve
```

The next pass should be a real projected-q-shell H1 probe, not a permanent test:
build final overlap/kinetic/nuclear/H through the new seams, verify the final
overlap is identity, run an ordinary symmetric eigensolve, and compare the final
Hamiltonian/energy to a shell-support oracle.

-- repo-manager@macmini
