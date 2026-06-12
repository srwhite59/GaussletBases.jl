Pass 061 complete: multi-layer PQS side13 pre-RHF H1/J diagnostic.

Files edited:

- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`

Ignored artifacts written:

- `tmp/work/pqs_multilayer_shell_side13_j_probe.jl`
- `tmp/work/pqs_multilayer_shell_side13_j_probe_summary.txt`

H1 and J diagnostic numbers:

- Fixture:
  - parent count: 13
  - mapping: `AsinhMapping(c=0.1,s=1.0,tail_spacing=10.0)`
  - endpoints: `(-8.565228460168399, 8.565228460168399)`
  - core: `(4:10)^3`
  - outer box: `(1:13)^3`
  - shell layers: 3
  - final dimension: 1549
  - final overlap identity error: `5.5061337433937041e-13`
- Z = 2 H1:
  - lowest energy: `-1.9755618232013417`
  - error vs hydrogenic `-2.0`: `0.02443817679865834`
- H1 self-Coulomb:
  - `J = 1.2169264388860319`
  - hydrogenic reference `5Z/8 = 1.25`
  - error vs reference: `-0.03307356111396809`

Density-interaction weight and symmetry diagnostics:

- Support raw pair numerator:
  - shape: `(2197, 2197)`
  - symmetry error: `0.0`
- Support weights:
  - count: 2197
  - min/max: `0.032401911987067523` / `9.138107632725481`
  - positive count: 2197
- Pre-final density interaction:
  - status: `materialized_pqs_complete_core_shell_pre_final_density_interaction`
  - density gauge: `pre_final_localized_positive_weight`
  - pre-final weight count: 1549
  - pre-final weight min/max: `0.032401911987067523` / `9.1381076327257`
  - positive / negative / near-zero counts: `1549 / 0 / 0`
  - pre-final pair matrix shape: `(1549, 1549)`
  - finite: true
  - symmetry error: `4.218847493575595e-15`
  - pre-final weight division applied: true
  - signed final weight division used: false
  - raw no-division used: false
  - fixed-block pair data authority used: false

Comparison to WL side13:

- WL side13 H1 reference: `-1.9748150892830352`
- PQS/WL H1 delta: `-0.0007467339183064858`
- WL side13 J reference: `1.2158294767735702`
- PQS/WL J delta: `0.0010969621124616857`

Decision:

- RHF is a reasonable next probe.
- Rationale: H1 and J are coherent against the WL side13 diagnostics, support/pre-final weights are positive, the pre-final pair matrix is finite and symmetric, and the probe used the intended pre-final density-interaction seam without fixed-block pair authority.
- This is still not an acceptance gate and does not settle the provisional PQS spacing rule.

Docs update:

- Added a concise non-acceptance status note with the H1/J numbers and the conclusion that side13 PQS RHF is reasonable as the next probe.

Validation:

- `julia --project=. tmp/work/pqs_multilayer_shell_side13_j_probe.jl` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage report:

- No production code or permanent tests were added in this pass.
- Nothing was deleted because this was a probe/status pass; it validated the current density seam rather than replacing an old implementation.
- The old fixed-block pair data surface became less relevant for this path because the multi-layer pre-final density interaction materialized without it.
- The new probe is an ignored developer artifact, not long-term test coverage.
- Remaining stale/duplicate surfaces to retire next: keep q-only one-shell PQS probes as route/scaling evidence only, and keep fixed-block pair data as oracle/reference rather than active authority.

-- repo-doer@macmini
