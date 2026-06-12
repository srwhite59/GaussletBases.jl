Pass 059 complete: PQS multi-layer shell source seam audit.

Files edited:

- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`

Ignored artifacts written:

- `tmp/work/pqs_multilayer_shell_source_seam_audit.jl`
- `tmp/work/pqs_multilayer_shell_source_seam_audit_summary.txt`

Exact one-cell-boundary code sites:

- `src/cartesian_nested_faces.jl:_nested_projected_q_shell_validate_boxes`
  - Enforces `first(inner) == first(current) + 1` and `last(inner) == last(current) - 1` on every axis.
  - Also requires source-mode dimensions to match `q/L` metadata and fit inside the current physical support box.
- `src/cartesian_nested_faces.jl:_nested_projected_q_shell_boundary_rectangular_pieces`
  - Rechecks the same one-cell boundary before decomposing a descriptor into six slab pieces.
- The data shape is therefore a single one-cell raw boundary shell around a strict inner box. A single descriptor cannot represent `(1:13)^3` around `(4:10)^3` directly.

Repeated one-cell composition viability:

- The prototype built three legal one-cell descriptors for the WL-aligned side13 geometry:
  - shell1: current `(3:11)^3`, inner `(4:10)^3`, raw dims `(7,7,7)`, support 386, retained 218
  - shell2: current `(2:12)^3`, inner `(3:11)^3`, raw dims `(9,9,9)`, support 602, retained 386
  - shell3: current `(1:13)^3`, inner `(2:12)^3`, raw dims `(11,11,11)`, support 866, retained 602
- Shell supports were disjoint:
  - combined shell duplicate count: 0
  - core/shell duplicate count: 0
- Core `(4:10)^3` plus the three shell supports covered the full side13 parent:
  - core support count: 343
  - combined shell support count: 1854
  - combined support count: 2197
  - parent support count: 2197
  - combined support covers parent: true
- Each shell realization was well-conditioned at the shell isometry level:
  - shell1 isometry error about `1.17e-14`
  - shell2 isometry error about `2.59e-14`
  - shell3 isometry error about `2.58e-14`

Final-basis realization impact:

- `pqs_complete_core_shell_final_basis(...)` does not need immediate multi-shell API generalization for the first side13 probe.
- The prototype block-concatenated the three shell isometries into one collapsed shell-sector coefficient matrix with shape `(1854, 1206)`.
- Passing the core plus collapsed shell sector to the existing helper succeeded:
  - final basis status: `available_pqs_complete_core_shell_final_basis`
  - final dimension: 1549
  - shell final retained count: 1206
  - final overlap identity error: about `5.51e-13`
- This proves the first implementation can keep the current final-basis helper and add a route-owned multi-layer shell source/plan in front of it.

Recommended smallest implementation path:

1. Add a route-owned multi-layer PQS shell source plan that accepts the WL-aligned parent/core geometry.
2. Build repeated legal one-cell shell descriptors per layer, from inner core outward:
   - `(3:11)^3/(4:10)^3`
   - `(2:12)^3/(3:11)^3`
   - `(1:13)^3/(2:12)^3`
3. Validate support disjointness and parent coverage.
4. Concatenate shell support rows and block-diagonal shell isometry matrices into a collapsed shell sector.
5. Feed the collapsed shell sector to `pqs_complete_core_shell_final_basis(...)`.
6. Only after that seam is route-owned, add final-basis/H1 for the side13 fixture. Do not run RHF until H1 is sensible.

Oracle/reference guidance:

- Useful references:
  - existing `_nested_projected_q_shell_layer`
  - `_nested_projected_q_shell_staged_unit_descriptor`
  - `CartesianContractedParentMetrics._pqs_shell_realization_plan`
  - existing complete core/shell final-basis transfer helpers
- Must not become route authority:
  - old fixed-block matrices
  - q-only one-shell ladder probes as physical gates
  - full-parent support/fixed-block shortcuts

Prototype result and timing:

- `tmp/work/pqs_multilayer_shell_source_seam_audit.jl` passed.
- Geometry:
  - parent count 13
  - `AsinhMapping(c=0.1,s=1.0,tail_spacing=10.0)`
  - coordinate endpoints `(-8.565228460168399, 8.565228460168399)`
  - central spacing about `0.11574195942963642`
- Timing:
  - shell descriptor build total: about `2.136491` s
  - combined final-basis build: about `1.886935` s
- No H1 or RHF was run.

Docs update:

- Added a concise note to the PQS near-term plan recording that repeated one-cell shell descriptors can be collapsed into the existing complete core/shell final-basis helper, so the next implementation target is the multi-layer source plan.

Validation:

- `julia --project=. tmp/work/pqs_multilayer_shell_source_seam_audit.jl` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage report:

- No production code or permanent tests were added.
- Nothing was deleted because this was an audit/prototype pass; the active missing piece is a new route-owned source-plan seam, not a dead implementation.
- The failed single-descriptor multi-shell attempt from pass 058 is now less relevant: repeated one-cell descriptors are the viable source shape.
- The new artifact is a temporary developer prototype, not long-term coverage.
- Remaining stale/duplicate surfaces to retire next: the one-shell q ladder should stay exploratory/reference-only for physics; fixed-block matrices should remain oracle/timing references only, not authority for multi-layer PQS.

-- repo-doer@macmini
