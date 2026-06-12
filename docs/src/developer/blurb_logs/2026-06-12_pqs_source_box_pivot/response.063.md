Pass 063 complete: side13 PQS core/shell ladder probe.

Files edited:

- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`

Ignored artifacts written:

- `tmp/work/pqs_multilayer_shell_side13_core_ladder_probe.jl`
- `tmp/work/pqs_multilayer_shell_side13_core_ladder_probe_summary.txt`

Probe setup:

- parent count: 13
- mapping: `AsinhMapping(c=0.1,s=1.0,tail_spacing=10.0)`
- outer box: `(1:13)^3`
- He HF reference: `-2.861679995612239`
- No full-parent CPB fallback, ordinary Cartesian IDA fallback, old fixed-block pair authority, GTO, driver route, export, or artifact path was used.

Ladder results:

| core side | shell layers | final dim | H1 Z=2 | J | RHF total | HF error | total time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 5 | 4 | 1429 | `-1.9747770050414415` | `1.2155987566118234` | `-2.8366493020529360` | `+0.0250306935593030` | `11.804098s` |
| 7 | 3 | 1549 | `-1.9755618232013417` | `1.2169264388860319` | `-2.8372556463894707` | `+0.0244243492227683` | `9.831373s` |
| 9 | 2 | 1717 | `-1.9755847315579729` | `1.2169096340532488` | `-2.8373268318655693` | `+0.0243531637466696` | `11.561454s` |
| 11 | 1 | 1933 | `-1.9755865712480807` | `1.2169106627468296` | `-2.8373293921422116` | `+0.0243506034700274` | `14.145297s` |

Additional diagnostics:

- All four points converged in 8 RHF iterations.
- Final overlap identity errors stayed small:
  - core side 5: `1.5734514385856713e-12`
  - core side 7: `5.5061337433937041e-13`
  - core side 9: `2.1616129025625597e-12`
  - core side 11: `4.5461377717881390e-12`
- Pre-final density weights were positive for every point.
- Pre-final pair matrices were finite and symmetric, with symmetry errors around `4e-15` to `8e-15`.
- Electron counts stayed at 2 to roundoff.

Timing split:

- Core side 5:
  - plan/final basis `2.994805s`
  - one-body/H1 `2.107032s`
  - density interaction `1.258491s`
  - RHF `5.443771s`
- Core side 7:
  - plan/final basis `1.032836s`
  - one-body/H1 `1.738500s`
  - density interaction `0.936585s`
  - RHF `6.123452s`
- Core side 9:
  - plan/final basis `1.105776s`
  - one-body/H1 `1.908234s`
  - density interaction `0.878349s`
  - RHF `7.669096s`
- Core side 11:
  - plan/final basis `1.174396s`
  - one-body/H1 `2.296156s`
  - density interaction `0.908820s`
  - RHF `9.765925s`
- Total all points: `47.342222s`

Interpretation:

- The ladder is monotonic but nearly plateaued after core side 7.
- The side13/core7 result is not a special artifact of that direct-core split.
- At fixed side13 parent/mapping, repartitioning core size and shell depth does not close the remaining roughly `0.02435` Hartree gap to the He HF reference.
- The next fixture question is the physical parent/mapping rule, not more core/shell repartitioning on the same side13 parent.
- This remains probe-only; no fixture was promoted as an acceptance gate.

Docs update:

- Added a concise result table/interpretation to the PQS near-term plan.
- The note explicitly says no new `Z,d,s,ns` fixture rule was adopted.

Validation:

- `julia --project=. tmp/work/pqs_multilayer_shell_side13_core_ladder_probe.jl` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage report:

- No production code or permanent tests were added in this pass.
- Nothing was deleted because this was a probe-only fixture-control audit; it did not obsolete an active implementation.
- The core side 7 split is less special after this pass; it remains one reference point in a shallow side13 ladder rather than a gate.
- The new ladder script is an ignored developer artifact, not long-term coverage.
- Remaining stale/duplicate surfaces to retire next: avoid promoting same-parent core/shell repartitioning probes as physical convergence evidence; keep fixed-block and q-only one-shell surfaces as oracle/scaling references only.

-- repo-doer@macmini
