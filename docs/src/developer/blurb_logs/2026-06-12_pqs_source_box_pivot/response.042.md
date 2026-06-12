Pass 042 response
=================

Status: blocked on a smaller exact route-owned implementation gap.

I did not implement the 223-function H1 solve in this pass. The live repo has a
223-mode complete one-center core-plus-shell object, but it is the private
White-Lindsey low-order seed/fixed-block packet:

- `_white_lindsey_low_order_materialized_seed_fixture`
- `_white_lindsey_low_order_materialized_seed_route_units`
- `_white_lindsey_low_order_materialized_seed_operator_inventory`

That object has the desired shape:

- core retained range: `1:125`
- shell retained range: `126:223`
- total retained dimension: `223`
- overlap identity error in the seed fixed block: `3.3306690738754696e-15`

But its operator source is `:nested_fixed_block` and its inventory status is
`:private_development_seed`. Using it as the active H1 path would make the old
fixed-block packet the route authority, which the blurb explicitly forbids.

Smaller exact blocker:

```julia
:missing_route_owned_combined_core_shell_retained_operator_blocks
```

More specifically, the current PQS final-basis seam only handles the
98-function retained boundary shell path:

- `pqs_source_shell_realization_final_basis`
- `pqs_source_shell_final_one_body_from_boundary_matrix`
- `pqs_source_shell_final_electron_nuclear_by_center_from_boundary_block`
- `pqs_source_shell_final_one_electron_hamiltonian`

It does not yet have a route-owned combined core/shell final-basis realization
or direct retained operator placement for:

- core-core blocks;
- core-shell blocks;
- shell-shell blocks;

for the required terms:

- overlap;
- kinetic;
- electron-nuclear by-center.

Probe artifact
--------------

Added ignored developer probe:

```text
tmp/work/pqs_complete_core_shell_h1_probe.jl
tmp/work/pqs_complete_core_shell_h1_probe_summary.txt
```

Key probe output:

```text
status	blocked
blocker	:missing_route_owned_combined_core_shell_retained_operator_blocks
smaller_blocker_context	:complete_223_geometry_exists_only_as_private_wl_seed_oracle
expected_total_retained_dimension	223
seed_retained_dimension	223
seed_core_range	1:125
seed_shell_range	126:223
seed_inventory_route_family	white_lindsey_low_order
seed_inventory_status	private_development_seed
seed_operator_source	nested_fixed_block
requested_active_route_runs	false
h1_solve_materialized	false
old_fixed_block_matrix_authority_used	false
current_route_safe_term_matrices_used	false
missing_combined_basis_object	:pqs_combined_core_shell_final_basis_realization
missing_block_families	(:core_core, :core_shell, :shell_shell)
missing_terms	(:overlap, :kinetic, :electron_nuclear_by_center)
```

The seed inventory did not expose a clean physical endpoint summary through the
route-owned PQS surface. I left that as a nonclaim in the artifact rather than
deriving physical extent from private seed internals.

Files changed
-------------

Tracked source/test files changed: none.

Ignored probe artifacts added under `tmp/work/`:

- `tmp/work/pqs_complete_core_shell_h1_probe.jl`
- `tmp/work/pqs_complete_core_shell_h1_probe_summary.txt`

Response files written:

- `.agent_handoffs/response.042.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.042.md`

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
- `git diff --check` passed.

Deletion/shrinkage report
-------------------------

- No CPBM/CFBR/source helper was deleted or simplified because no new
  route-owned combined core/shell path was implemented.
- No old oracle-only path became unnecessary. The private White-Lindsey
  low-order seed remains useful as an oracle/reference, but it must not be
  promoted to active PQS H1 authority.
- No new permanent test was added. The only new executable artifact is the
  ignored focused probe under `tmp/work`, so it does not increase routine test
  burden.
- Remaining stale/duplicate surface to retire next: the boundary-shell-only
  98-function H1 path is still a useful focused mechanical contract but should
  stay explicitly nonphysical until a combined core/shell direct-retained route
  exists. The private seed/fixed-block 223-mode packet should remain oracle-only
  until route-owned core-core/core-shell/shell-shell operator placement exists.

Recommended next implementation blurb
-------------------------------------

Add a narrow `pqs_combined_core_shell_final_basis_realization` object and direct
retained one-body block placement for core-core, core-shell, and shell-shell
sectors. It should consume route-owned complete one-center source/final-basis
data and materialize overlap, kinetic, and uncharged by-center nuclear blocks in
the 223-dimensional final basis without using `_pqs_current_route_safe_term_matrices`
or the private nested fixed-block packet as authority.

-- repo-doer@macmini
