Pass 058 complete: PQS WL-aligned He fixture design audit.

Files edited:

- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`

Ignored artifacts written:

- `tmp/work/pqs_wl_aligned_fixture_design_audit.jl`
- `tmp/work/pqs_wl_aligned_fixture_design_audit_summary.txt`

WL fixture controls found:

- Tiny active/readiness He fixture: q/ns = 5/5, Z = 2, `d = 0.2`, `s = sqrt(dZ) = 0.6324555320336759`, `tail_spacing = 10.0`, `parent_side_count = 7`, physical endpoints about `+/-0.9666` bohr. This remains a compact one-shell bridge/oracle-style acceptance fixture.
- Larger WL He side13 probe: `AsinhMapping(c = 0.1, s = 1.0, tail_spacing = 10.0)`, parent side 13, endpoints about `+/-8.565`, retained dimension 517, 105 retained units, 5,565 unit pairs.
- Best fixed-ns exploratory gausslet-only He point: `d = 0.075`, `s = 0.75`, side 17, retained dimension 713, RHF total about `-2.858531351214`; exploratory, not a routine acceptance gate.
- GTO supplement enters only after the gausslet-only final-basis route is valid. H/H2+ have final-basis residual GTO acceptance; He+GTO still needs consistent final-basis electron-electron work.

PQS fixture controls found:

- Current complete core/shell q ladder uses `current_box = (1:N)^3`, `inner_box = (2:N-1)^3`, `raw_source_dims = (N-2)^3`.
- That means every q5/q7/q9/q11 point has exactly one surrounding shell layer.
- The q ladder uses `AsinhMapping(a = 0.25, s = chosen_for_R8, tail_spacing = 10.0)`, so it keeps box radius near 8 while varying q/source/core size. It is route/scaling evidence, not WL-style physical convergence.
- The PQS physical fixture knobs are parent axis count, `AsinhMapping(a,s,tail_spacing)` or equivalent `c/s/tail`, `current_box`, `inner_box`, raw source dimensions, `q/L/selected_q`, surrounding shell-layer count, and final dimension.

Multi-shell support result:

- The cheap side13/core7/three-shell final-basis smoke was blocked before final-basis realization.
- Attempted fixture:
  - parent count 13
  - `AsinhMapping(c = 0.1, s = 1.0, tail_spacing = 10.0)`
  - `current_box = (1:13)^3`
  - `inner_box = (4:10)^3`
  - `raw_source_dims = (7,7,7)`
  - `q = 7`
  - intended three surrounding shell layers
- Exact blocker:
  - `ArgumentError: projected q-shell requires a one-cell raw boundary around a strict inner box`
- Source audit found this guard in `src/cartesian_nested_faces.jl`.
- Conclusion: the current complete core/shell PQS construction does not yet support more than one surrounding shell by only changing `current_box` / `inner_box` / raw source dimensions. It needs a multi-layer PQS shell/source producer; extending the q-only ladder is the wrong next move.

Recommended first WL-aligned PQS fixture:

- After the multi-layer shell/source blocker is fixed, try the side13-style final-basis/H1 fixture:
  - parent count 13
  - `AsinhMapping(c = 0.1, s = 1.0, tail_spacing = 10.0)`
  - `current_box = (1:13)^3`
  - `inner_box = (4:10)^3`
  - `raw_source_dims = (7,7,7)`
  - `q = 7`
  - three surrounding shell layers
- This keeps the adequate R >= 8 box and WL-reviewed d/s controls while adding real shell depth. Run final-basis/H1 first; do not interpret RHF until that seam is sane.

Docs update:

- Added a concise note to `pqs_near_term_final_basis_realization_plan.md` saying future PQS physical tests should be WL-aligned fixture probes, not more q-only scaling, and recording the current one-cell-boundary blocker.

Validation:

- `julia --project=. tmp/work/pqs_wl_aligned_fixture_design_audit.jl` passed and wrote the summary artifact.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage report:

- No production code or permanent tests were added.
- Nothing was deleted because this was a read-only/probe-only audit and the blocker is an active missing producer, not a dead surface.
- The q9/q11 q-only ladder became less plausible as future gate material; it remains exploratory route/scaling evidence only.
- The new probe is an ignored developer artifact, not new long-term coverage.
- Remaining stale/duplicate surface to retire next: avoid treating one-cell complete core/shell q-scaling probes as physical convergence fixtures; add a true multi-layer shell/source route before any PQS He gate promotion.

-- repo-doer@macmini
