Pass 045 response
=================

Status: implemented the focused 223-function complete core/shell H1 operator
transfer and ordinary symmetric solve.

Result
------

The complete one-center route now runs mechanically through:

```text
route-owned complete core/surrounding-shell final basis
-> route-owned product-factor support operators
-> final-basis one-body transfer
-> H = K + Z * V_center
-> ordinary symmetric H1 solve
```

New/changed CFBR functions in
`src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`:

- `pqs_complete_core_shell_final_one_body_matrix`
- `pqs_complete_core_shell_final_one_electron_hamiltonian`
- `pqs_complete_core_shell_final_h1_solve`
- private validators/metadata helpers for those functions

No package/public export was added.

Probe artifact
--------------

Updated ignored H1 probe:

```text
tmp/work/pqs_complete_core_shell_h1_probe.jl
tmp/work/pqs_complete_core_shell_h1_probe_summary.txt
```

Key output:

```text
status	materialized_pqs_complete_core_shell_final_h1_solve
blocker	nothing
final_dimension	223
core_range	1:125
shell_range	126:223
final_overlap_identity_error	7.0978733509910171e-14
transferred_overlap_identity_error	7.0887740122316245e-14
hamiltonian_symmetry_error	3.9968028886505635e-15
hamiltonian_finite	true
h1_energy	-2.0638461028784776e+00
h1_error_vs_minus_half	-1.5638461028784776e+00
solve_kind	ordinary_symmetric
route_owned_product_operators_used	true
old_fixed_block_matrix_authority_used	false
current_route_safe_term_matrices_used	false
generalized_overlap_solve_used	false
nuclear_charge_applied_stage	hamiltonian_assembly
nuclear_centers_summed_stage	hamiltonian_assembly
ida_data_materialized	false
density_density_materialized	false
```

Important numerical caveat:

The route now runs, but the H / Z=1 energy is not a physics acceptance result.
The observed value, `-2.0638461028784776`, is far below the exact hydrogen value
`-0.5`. I did not hide this with a generalized solve or fallback route. The
next pass should audit the product-factor one-body convention for the combined
support operator, especially the nuclear/operator weighting convention, before
promoting this path to scientific acceptance.

Oracle comparison:

`oracle_delta_vs_fixed_block` is recorded as
`not_available_mapping_mismatch`. The old fixed-block packet remains
oracle-only, and I did not use it as active matrix authority.

Files edited
------------

Tracked source:

- `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`

Ignored probe artifacts:

- `tmp/work/pqs_complete_core_shell_h1_probe.jl`
- `tmp/work/pqs_complete_core_shell_h1_probe_summary.txt`

Validation
----------

Ran:

```text
julia --project=. tmp/work/pqs_complete_core_shell_h1_probe.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
```

Results:

- focused H1 probe completed and materialized the 223-function H1 solve;
- load check passed with `load ok`;
- `git diff --check` passed.

Deletion/shrinkage report
-------------------------

- The old fixed-block packet is less necessary as a route comparator for
  mechanical H1 materialization; it remains oracle-only and was not consumed as
  active authority.
- No source/test path was deleted because this pass only added the missing
  final-basis one-body transfer and H1 seam.
- No permanent test was added. The validation remains an ignored focused probe,
  so routine test burden did not grow.
- No compatibility shim, generalized solve, IDA, density-density, RHF, GTO,
  driver, export, or artifact path was added.
- Remaining stale/duplicate surfaces to retire next:
  - the 98-function boundary-shell H1 path should stay mechanical/nonphysical;
  - old fixed-block H1 comparison should remain oracle-only;
  - the next blocker is scientific/operator convention validation for the
    route-owned product-factor H1 value, not basis independence or mechanical
    operator placement.

-- repo-doer@macmini
