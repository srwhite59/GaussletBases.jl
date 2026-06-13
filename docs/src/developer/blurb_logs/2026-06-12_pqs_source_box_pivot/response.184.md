Pass 184 complete.

Ignored local probe files:
- `tmp/work/pqs_he_419_wl_parent_mapping_probe.jl`
- `tmp/work/pqs_he_419_wl_parent_mapping_probe_summary.txt`
- Confirmed ignored by `.gitignore:12:tmp/`.

Probe command:

```text
julia --project=. tmp/work/pqs_he_419_wl_parent_mapping_probe.jl
```

Parent mapping:
- mapping: `AsinhMapping(a=0.38729833462074165, s=0.7745966692414834, tail_spacing=10.0)`
- `asinh_a = 0.38729833462074165`
- `asinh_s = 0.7745966692414834`
- `asinh_c = 0.3`
- physical endpoints: `(-5.892850307983052, 5.892850307983052)`
- constructor matched pass-183 WL:
  `MappedUniformBasisSpec(:G10; count = 11, mapping = white_lindsey_atomic_mapping(Z = 2.0, d = 0.3, tail_spacing = 10.0), reference_spacing = 1.0)`

Source/final-basis inventory:
- source plan status: `available_pqs_multilayer_shell_source_plan`
- region plan status: `available_pqs_multilayer_shell_region_plan`
- source kind:
  `shellification_backed_repeated_one_cell_projected_q_shell_layers`
- shell layer count: `3`
- source dims per shell: `((5, 5, 5), (5, 5, 5), (5, 5, 5))`
- retained per shell: `(98, 98, 98)`
- fixed source mode shape used: `true`
- core support count: `125`
- shell support count: `1206`
- shell final retained count: `294`
- final dimension: `419`
- final overlap identity error: `3.465561171367426e-13`

WL-parent PQS diagnostics:
- H1 status: `materialized_pqs_multilayer_complete_core_shell_h1_payload`
- H1 lowest energy: `-1.991334820314074`
- H1 matrix finite: `true`
- H1 symmetry error: `1.7763568394002505e-14`
- H1-J status:
  `materialized_pqs_multilayer_complete_core_shell_h1_j_payload`
- H1-J self-Coulomb: `1.2420423900074902`
- density gauge: `pre_final_localized_positive_weight`
- density-interaction status:
  `materialized_pqs_complete_core_shell_pre_final_density_interaction`
- private RHF status:
  `materialized_pqs_multilayer_complete_core_shell_rhf_scf_payload`
- private RHF blocker: `nothing`
- private RHF converged: `true`
- private RHF iterations: `7`
- private RHF total energy: `-2.850817886618113`
- private RHF one-body energy: `-3.8713594674599054`
- private RHF two-body energy: `1.0205415808417926`
- private RHF density trace: `2.0000000000000004`
- private RHF commutator residual: `2.8556726532069643e-9`
- private RHF idempotency error: `3.469446951953614e-17`
- elapsed time: `7.6675995` seconds

Deltas versus pass-183 WL 419:
- WL-parent PQS RHF total - WL RHF total:
  `-1.4383600322798173e-5`
- WL-parent PQS H1 - WL H1:
  `9.649649361120893e-6`
- WL-parent PQS H1-J self-Coulomb - WL H1 self-Coulomb:
  `-4.997485057112172e-6`
- WL-parent PQS RHF one-body - WL RHF one-body:
  `4.944121432171045e-5`
- WL-parent PQS RHF two-body - WL RHF two-body:
  `-6.382481464473067e-5`

Deltas versus current custom-mapping PQS 419:
- WL-parent PQS RHF total - custom-mapping PQS RHF total:
  `-0.0014307641877961963`
- WL-parent PQS H1 - custom-mapping PQS H1:
  `-0.004652845139180295`
- WL-parent PQS H1-J self-Coulomb - custom-mapping PQS H1-J self-Coulomb:
  `0.0158797896955718`
- WL-parent PQS RHF one-body - custom-mapping PQS RHF one-body:
  `-0.0037948736054369547`
- WL-parent PQS RHF two-body - custom-mapping PQS RHF two-body:
  `0.0023641094176405364`

Interpretation:
- The WL-parent PQS result moved from the pass-183 raw PQS/WL RHF gap of
  about `+0.00141638 Ha` to `-1.43836e-5 Ha`.
- The H1 gap moved from about `+0.00466249 Ha` to `+9.65e-6 Ha`.
- The H1-J self-Coulomb gap moved from about `-0.0158848 Ha` to `-4.997e-6 Ha`.
- This strongly indicates the apparent pass-183 PQS/WL gap was mostly a parent
  lattice/mapping mismatch, not evidence of a large PQS contraction/source
  representation error.
- At this scale, the remaining diagnostic difference is small; the largest
  remaining RHF component difference is the two-body term at about
  `-6.38e-5 Ha`, partly offset by the one-body term at about `+4.94e-5 Ha`.
- I did not tune PQS or WL and did not promote private RHF beyond a Hamiltonian
  diagnostic.

Validation:
- `julia --project=. tmp/work/pqs_he_419_wl_parent_mapping_probe.jl`
  - passed and wrote
    `tmp/work/pqs_he_419_wl_parent_mapping_probe_summary.txt`
- `julia --project=. -e 'pqs_one=-3.8713594674599054; wl_one=-3.871408908674227; custom_one=-3.8675645938544685; println("delta_wl_rhf_one_body=", pqs_one-wl_one); println("delta_custom_pqs_rhf_one_body=", pqs_one-custom_one)'`
  - printed the one-body deltas above
- `git status --short --branch`
  - before writing response files: `## main...origin/main`
- `git diff --numstat -- src test tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - empty
- `git check-ignore -v tmp/work/pqs_he_419_wl_parent_mapping_probe.jl tmp/work/pqs_he_419_wl_parent_mapping_probe_summary.txt`
  - both ignored by `.gitignore:12:tmp/`

Tracked source/test/generator changes:
- none

Deletion/shrinkage report:
- deleted: none
- simplified: none; this was an ignored local diagnostic probe
- quarantined: WL-parent PQS comparison remains an ignored local probe and
  summary, not a tracked test or route behavior
- not deleted because: ignored probe/summary are local evidence for the mapping
  mismatch decision
- exact remaining caller/blocker: no tracked caller added; manager decision
  remains whether to align the tracked focused He fixture parent mapping or
  keep this as a diagnostic-only probe

-- repo-doer@macmini
