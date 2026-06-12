Accepted.

Pass 025 added the intended by-center final-basis nuclear helper:

```text
pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block(...)
```

It consumes the live `PairBlockMaterializationResult` for
`:retained_source_electron_nuclear_by_center`, validates retained
PQS-source-mode block space, preserves center metadata, and computes:

```text
V_final(center) = L' * V_boundary(center) * L
```

The helper keeps the nuclear convention intact:

```text
nuclear charge recorded = true
nuclear charge applied = false
centers summed = false
uncharged by-center convention = true
```

and does not claim H1, Hamiltonian assembly, IDA, density-density, RHF, driver,
export, or artifact behavior.

Manager validation:

- `julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Next blocker:

```text
:missing_pqs_final_one_electron_hamiltonian_assembly
```

The next pass should add Hamiltonian-stage assembly only: consume final kinetic
and separated final by-center nuclear matrices, apply recorded nuclear charges,
sum centers, and produce the final one-electron Hamiltonian matrix. Do not run
an eigensolve yet.

-- repo-manager@macmini
