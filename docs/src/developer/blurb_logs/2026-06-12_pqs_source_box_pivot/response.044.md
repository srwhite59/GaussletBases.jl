Pass 044 response
=================

Status: implemented the route-owned complete core/surrounding-shell
final-basis realization. No H1 operator placement was added.

Result
------

The independent `223`-function final basis now materializes through a compact
CFBR helper:

```julia
GaussletBases.CartesianFinalBasisRealization.pqs_complete_core_shell_final_basis
```

Object kind:

```julia
:pqs_complete_core_shell_final_basis
```

The focused probe uses the intended geometry:

```text
current_box = (1:7, 1:7, 1:7)
inner_box = (2:6, 2:6, 2:6)
raw_source_dims = (5, 5, 5)
direct core support count = 125
surrounding shell support count = 218
shell retained count = 98
total final dimension = 223
```

Key probe numbers:

```text
status	available_pqs_complete_core_shell_final_basis
blocker	nothing
core_support_count	125
shell_support_count	218
shell_retained_count	98
final_retained_count	223
core_range	1:125
shell_range	126:223
core_shell_support_disjoint	true
pre_final_overlap_rank	223
pre_final_overlap_expected_dimension	223
pre_final_overlap_eigenvalue_min	9.9999999999970413e-01
pre_final_overlap_eigenvalue_max	1.0000000000002984e+00
shell_final_overlap_identity_error	9.2592600253738055e-14
final_overlap_identity_error	7.0978733509910171e-14
final_basis_materialized	true
old_fixed_block_matrix_authority_used	false
current_route_safe_term_matrices_used	false
generalized_overlap_solve_used	false
```

This resolves pass 043's rank-deficient basis attempt. The successful path uses
an independent surrounding shell sector from `current_box = (1:7)^3` minus
`inner_box = (2:6)^3`, instead of treating the 98 shell directions as extra
columns inside the same `5 x 5 x 5` source-mode space.

Files edited
------------

Tracked source files:

- `src/cartesian_final_basis_realization/CartesianFinalBasisRealization.jl`
  - includes the new final-basis helper file;
  - no package/public export was added.

- `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`
  - adds `pqs_complete_core_shell_final_basis`;
  - builds the combined direct-core plus surrounding-shell pre-final overlap;
  - applies a combined Lowdin cleanup;
  - returns final overlap diagnostics and final coefficients;
  - explicitly does not materialize one-body operators, H1, IDA, RHF, driver
    wiring, exports, or artifacts.

Ignored probe artifacts:

- `tmp/work/pqs_complete_core_shell_final_basis_probe.jl`
- `tmp/work/pqs_complete_core_shell_final_basis_probe_summary.txt`

Validation
----------

Ran:

```text
julia --project=. tmp/work/pqs_complete_core_shell_final_basis_probe.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --no-index --check /dev/null src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl
```

Results:

- focused final-basis probe passed and wrote the available `223`-basis summary;
- load check passed with `load ok`;
- `git diff --check` passed;
- no-index whitespace check on the new source file produced no diagnostics.

Deletion/shrinkage report
-------------------------

- The old White-Lindsey low-order fixed-block packet is less necessary for
  final-basis realization. It remains useful as an oracle for geometry/counts,
  but not as active final-basis authority.
- No old source/test path was deleted because this pass only fills the missing
  final-basis realization seam; operator placement and H1 are still absent.
- No new permanent test was added. The validation is an ignored focused probe,
  matching the blurb's request and avoiding routine test burden.
- No compatibility shim was added.
- Remaining stale/duplicate surfaces to retire next:
  - the 98-function boundary-shell H1 path should remain mechanical/nonphysical;
  - the private `white_lindsey_low_order` seed/fixed-block path should stay
    oracle-only;
  - pass 045 can now add route-owned core-core, core-shell, and shell-shell
    overlap/kinetic/by-center nuclear placement against the new final-basis
    realization.

-- repo-doer@macmini
