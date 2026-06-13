Pass 182 complete.

Tracked edits:
- None to `src/`, `test/`, docs, generator, or artifacts other than response
  handoff files.

Ignored local probe files:
- `tmp/work/pqs_he_419_rhf_diagnostic_probe.jl`
- `tmp/work/pqs_he_419_rhf_diagnostic_probe_summary.txt`
- Confirmed ignored by `.gitignore:12:tmp/`.

Probe command:

```text
julia --project=. tmp/work/pqs_he_419_rhf_diagnostic_probe.jl
```

Probe setup:
- Included `test/nested/pqs_direct_retained_final_h1_runtests.jl` to reuse
  `_pqs_h1_complete_fixture()`.
- Rebuilt the corrected fixed-q one-center He `419`-dimensional final basis.
- Built the same H1 payload as pass 180.
- Built density inputs with
  `GaussletBases._pqs_source_box_route_driver_complete_core_shell_density_inputs(...)`.
- Built the same H1-J / pre-final density-interaction payload as pass 181.
- Built the private RHF input contract with:
  - `electron_count = 2`
  - `fixture_role = :route_smoke`
- Ran private diagnostic RHF SCF with:
  - `mixing_kind = :fock_diis`
  - `max_iterations = 50`
  - `residual_atol = 1.0e-8`
  - `density_atol = 1.0e-8`
  - `energy_atol = 1.0e-10`

Probe results:
- input contract status:
  `available_pqs_multilayer_complete_core_shell_rhf_input_contract`
- input contract blocker: `nothing`
- SCF status: `materialized_pqs_multilayer_complete_core_shell_rhf_scf_payload`
- SCF blocker: `nothing`
- converged: `true`
- iterations: `7`
- final total energy: `-2.8493871224303167`
- one-body energy: `-3.8675645938544685`
- two-body energy: `1.018177471424152`
- density trace: `1.9999999999999993`
- commutator residual: `3.2524884209517158e-9`
- idempotency error: `4.163336342344337e-17`
- H1 energy from input path: `-1.9866819751748936`
- H1-J self-Coulomb from input path: `1.2261626003119184`
- density inputs status: `available_complete_core_shell_density_inputs`
- density-interaction status:
  `materialized_pqs_complete_core_shell_pre_final_density_interaction`
- density gauge: `pre_final_localized_positive_weight`
- DIIS used count: `6`
- DIIS fallback count: `0`
- private diagnostic only: `true`
- probe elapsed time reported by Julia `@elapsed`: `3.637304334` seconds

Validation:
- `julia --project=. tmp/work/pqs_he_419_rhf_diagnostic_probe.jl`
  - passed; included focused H1/H1-J test passed `51/51` in `6.3s`
  - probe SCF converged as reported above
- `git status --short --branch`
  - before writing response files: `## main...origin/main`
- `git diff --numstat -- src test`
  - empty
- `git check-ignore -v tmp/work/pqs_he_419_rhf_diagnostic_probe.jl tmp/work/pqs_he_419_rhf_diagnostic_probe_summary.txt`
  - both ignored by `.gitignore:12:tmp/`

Interpretation:
- This is private diagnostic/probe behavior only.
- I am not calling this a physics endpoint.
- I did not compare it as accepted against White--Lindsey. Any later comparison
  should state explicitly that the old WL reference is Fig.8-style supplemented
  WL, while this PQS value is gausslet-only/no AHGBS residual.

Line-budget / diff check:
- `git diff --numstat -- src test` is empty, as required for this no-tracked-
  source/test probe pass.

Deletion/shrinkage report:
- deleted: none
- simplified: none; this was an ignored local diagnostic probe only
- quarantined: RHF remains private diagnostic/probe behavior in ignored
  `tmp/work`, with no route-driver adoption
- not deleted because: ignored probe and summary remain as local evidence for
  the pass-182 scalar result
- exact remaining caller/blocker: no tracked caller added; next decision is
  whether manager wants to promote any RHF scalar check or keep it diagnostic

-- repo-doer@macmini
