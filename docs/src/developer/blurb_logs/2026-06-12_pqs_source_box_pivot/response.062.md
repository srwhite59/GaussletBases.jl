Pass 062 complete: side13 multi-layer PQS He RHF probe.

Files edited:

- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`

Ignored artifacts written:

- `tmp/work/pqs_multilayer_shell_side13_rhf_probe.jl`
- `tmp/work/pqs_multilayer_shell_side13_rhf_probe_summary.txt`

RHF result:

- Fixture:
  - parent count: 13
  - mapping: `AsinhMapping(c=0.1,s=1.0,tail_spacing=10.0)`
  - endpoints: `(-8.565228460168399, 8.565228460168399)`
  - core: `(4:10)^3`
  - outer box: `(1:13)^3`
  - shell layers: 3
  - final dimension: 1549
  - final overlap identity error: `5.5061337433937041e-13`
- Pre-RHF checks carried through:
  - Z = 2 H1: `-1.9755618232013417`
  - H1 self-Coulomb J: `1.2169264388860319`
  - density interaction: `materialized_pqs_complete_core_shell_pre_final_density_interaction`
  - pre-final pair matrix: `(1549,1549)`, finite, symmetry error `4.218847493575595e-15`
  - pre-final weights: all positive
  - signed final weight division used: false
  - raw no-division used: false
  - fixed-block pair data authority used: false
- RHF:
  - converged: true
  - iterations: 8
  - one-electron energy: `-3.8469458612019864`
  - electron-electron energy: `1.0096902148125155`
  - total energy: `-2.8372556463894707`
  - final density trace: `1.0`
  - pre-final density trace: `1.0000000000000038`
  - electron count: `2.0`
  - pre-final electron count proxy: `2.0000000000000075`
  - Fock symmetry error: `0.0`

Timing:

- plan/final basis: `3.184488` s
- one-body/H1: `2.253455` s
- density interaction: `1.189725` s
- RHF solve: `6.461152` s
- total: `13.088820` s

Comparison:

- He HF reference: `-2.861679995612239`
- Error vs He HF reference: `+0.02442434922276826` Hartree
- WL side13 RHF reference: `-2.8364979997009137`
- Delta vs WL side13: `-0.0007576466885570454` Hartree

Convention status:

- The final/pre-final RHF convention remains accepted for this probe:
  - final one-electron Hamiltonian in the orthonormal final basis;
  - pre-final positive-weight density interaction;
  - `c_prefinal = combined_lowdin_cleanup * c_final`;
  - no signed-final-weight density;
  - no raw no-division density;
  - no fixed-block pair authority.

Decision:

- The side13 multi-layer PQS RHF result is physically sane and close to the WL side13 probe.
- It remains probe-only/non-acceptance.
- It suggests a future compact gate candidate only after the PQS fixture rule is reviewed and settled.

Docs update:

- Added a concise non-acceptance note with the RHF result, comparison to WL side13, and fixture-rule caveat.

Validation:

- `julia --project=. tmp/work/pqs_multilayer_shell_side13_rhf_probe.jl` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.
- `git diff --check` passed.

Deletion/shrinkage report:

- No production code or permanent tests were added in this pass.
- Nothing was deleted because this was a probe/status pass, not a cleanup or replacement pass.
- The old fixed-block pair data surface became less necessary as active authority because the side13 RHF probe used the route-owned multi-layer source plan and pre-final density interaction.
- The new RHF script is an ignored developer artifact, not permanent coverage.
- Remaining stale/duplicate surfaces to retire next: keep q-only one-shell probes and fixed-block pair data as exploratory/oracle references, not acceptance route authority.

-- repo-doer@macmini
