Review 035: accepted.

The pass did the intended test ownership cleanup. The moved
`CartesianFinalBasisRealization` behavior now has compact module-contract
coverage, and the large CPBM contract file no longer carries the detailed
checks for those moved functions.

Independent manager validation:

```text
julia --project=. test/nested/cartesian_final_basis_realization_contract_runtests.jl
julia --project=. test/nested/cartesian_pair_block_materialization_contract_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
```

The new CFBR test passed with 31 checks in about 0.4s. The focused CPBM
contract still passed in about 13s.

The CPBM contract shrinkage is real: the moved final-basis checks were removed,
the CPBM section keeps only a tiny alias smoke plus CPBM-owned nuclear and
Hamiltonian checks, and the new test replaces rather than duplicates the moved
coverage.

Next pass can move to algorithm shape: direct retained-boundary overlap and
kinetic construction, compared against the existing raw-source-block selector
path as oracle.

-- repo-manager@macmini
