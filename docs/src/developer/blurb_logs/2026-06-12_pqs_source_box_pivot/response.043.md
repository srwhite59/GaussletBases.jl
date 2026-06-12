Pass 043 response
=================

Status: blocked on the smaller combined-basis realization issue.

I did not add production source code. I initially tested a compact CPBM helper
shape, but removed it before handoff because it would not enable the requested
route and would add source carrying cost without a live consumer. The focused
probe now records the exact failure without changing production code.

Result
------

The attempted active basis

```text
125 raw/core source modes + 98 shell Lowdin-cleaned boundary modes
```

does not form a 223-dimensional independent final basis with the currently
available shell realization data. Its pre-final combined overlap has shape
`(223, 223)` but rank `125`.

Exact blocker:

```julia
:missing_route_owned_shell_projection_lowdin_for_combined_basis
```

Rank diagnostic:

```julia
:combined_core_shell_overlap_rank_deficient
```

This is more precise than pass 042's broader operator-block blocker. The
core-shell/shell-shell direct product block loops are not the first blocker:
the available shell Lowdin data describes boundary modes inside the same
5x5x5 source-mode space, so combining all 125 source modes with those 98
shell-final directions duplicates the source space. A route-owned complete
core-plus-surrounding-shell realization must define the shell as an independent
surrounding-shell final sector, not as extra directions in the same 125-mode
source space.

Probe artifact
--------------

Updated ignored developer probe:

```text
tmp/work/pqs_complete_core_shell_h1_probe.jl
tmp/work/pqs_complete_core_shell_h1_probe_summary.txt
```

Key output:

```text
status	blocked
blocker	missing_route_owned_shell_projection_lowdin_for_combined_basis
exact_rank_blocker	combined_core_shell_overlap_rank_deficient
expected_core_retained_count	125
expected_shell_retained_count	98
expected_total_retained_dimension	223
core_range	1:125
shell_range	126:223
core_source_mode_count	125
boundary_source_mode_count	98
shell_final_retained_count	98
pre_final_overlap_shape	(223, 223)
pre_final_overlap_rank	125
pre_final_overlap_expected_dimension	223
pre_final_overlap_symmetry_error	6.3108872417680944e-30
pre_final_overlap_eigenvalue_min	-2.2527994226297727e-15
pre_final_overlap_eigenvalue_max	2.0000000000002682e+00
final_overlap_identity_error	not_materialized
h1_energy	not_materialized
solve_kind	not_materialized
old_fixed_block_matrix_authority_used	false
current_route_safe_term_matrices_used	false
generalized_overlap_solve_used	false
full_parent_cpb_fallback_used	false
shell_support_projection_active_path_used	false
```

Files edited
------------

Tracked source/test files changed: none.

Ignored probe artifacts updated under `tmp/work/`:

- `tmp/work/pqs_complete_core_shell_h1_probe.jl`
- `tmp/work/pqs_complete_core_shell_h1_probe_summary.txt`

No production object/function was added or changed.

Validation
----------

Ran:

```text
julia --project=. tmp/work/pqs_complete_core_shell_h1_probe.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
```

Results:

- focused probe completed and wrote the blocked summary artifact;
- load check passed with `load ok`;
- `git diff --check` passed;
- final `git status --short --branch` was clean.

Deletion/shrinkage report
-------------------------

- No production code, test, metadata, or compatibility path became unnecessary
  because the 223-dimensional route did not materialize.
- I deleted the temporary production helper attempt before handoff. That avoided
  adding a blocked CPBM source path with no live consumer.
- No existing source/test surface was simplified because the pass ended at a
  design blocker, not a successful route replacement.
- No new permanent test was added. The probe is ignored under `tmp/work`, so it
  is exploratory/live-audit coverage rather than routine test burden.
- Remaining stale/duplicate surfaces to retire next:
  - the 98-function boundary-shell H1 path should remain explicitly mechanical
    and nonphysical;
  - the private White-Lindsey low-order 223 fixed-block packet should remain
    oracle-only;
  - the next implementation should create the route-owned complete
    core-plus-surrounding-shell final-basis realization before adding combined
    H1 operator placement.

Recommended next implementation target
--------------------------------------

Build or expose a route-owned one-center complete-shell final-basis realization
whose sectors are genuinely independent:

```text
inner 5x5x5 core direct sector
+ surrounding shell sector
-> combined final overlap / Lowdin cleanup
```

Only after that object is available should the route add core-core,
core-shell, and shell-shell retained overlap/kinetic/by-center nuclear block
placement and the ordinary symmetric H1 solve.

-- repo-doer@macmini
