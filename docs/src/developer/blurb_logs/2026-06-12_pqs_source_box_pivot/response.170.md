Pass 170 implementation complete.

Files changed:
- `tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`

Line-budget:
- `src`/`test` line-budget was not triggered.
- `git diff --numstat -- src test` was empty.

WL source path actually used:
- In-process route-configured diatomic atom-growth White-Lindsey route.
- Route inputs mirror the accepted Be2 config shape: `route_family = :white_lindsey_low_order`, Be2 charges/locations, `parent_axis_counts = (x = 9, y = 7, z = 9)`, `q = 5`, `n_s = 5`, `materializer_backend = :pgdg_localized_experimental`, `materializer_nside = 5`, `low_order_shellization_policy = :atom_growth_complete_rectangular`, and `route_configured_diatomic_ham_interaction_treatment = :ggt_nearest`.
- The generator reads the atom-growth `basis_adapter` and `ham_adapter` directly from the route materialization. It does not promote the old WL seed/one-center/debug paths.

Generated output paths:
- `tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_inspection_bundle.jld2`
- `tmp/work/be2_wl_pqs_cr2_inspection_artifact/be2_wl_pqs_handoff_fingerprint.tsv`
- `git status --ignored --short tmp/work/be2_wl_pqs_cr2_inspection_artifact` shows both generated outputs as ignored.
- Bundle size observed: `241M`; TSV size observed: `884B`.

Readback summary:
- `producer/repo_commit = 2faa1e8189d887595396f7676264b47f8c21f34c`
- `producer/dirty = true` because the generator edit is local and uncommitted.
- PQS route:
  - status: `available_diatomic_complete_core_shell_hamiltonian_handoff_payload`
  - CR2 read-only inspector ready: `true`
  - CR2 solver ready: `false`
- WL route:
  - status: `available_route_configured_diatomic_atom_growth_ham_payload`
  - final dimension: `2287`
  - overlap shape: `(2287, 2287)`
  - one-body/H1 shape: `(2287, 2287)`
  - two-body shape: `(2287, 2287)`
  - two-body representation kind: `final_basis_density_density_matrix`
  - interaction treatment: `ggt_nearest`
  - CR2 read-only inspector ready: `true`
  - CR2 solver ready: `false`
  - CR2 export ready: `false`

Unavailable placeholders and blockers:
- supplement/residual-GTO: `unavailable`
- Qiu-White atom-local HF inputs: `unavailable`
- correction/EGOI/stationary/cusp metadata: `unavailable`
- MWG/IDA route-configured diatomic ham support: `pending_route_configured_diatomic_mwg_operator_support`
- solver/HFDMRG/RHF/DMRG/HamV6 export readiness: false
- WL CR2 handoff blocker after read-only population: `missing_hfdmrg_density_density_contract`

Validation commands/results:
- `julia --project=. tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl`
  - First run reached the WL operator build and failed on a mapper field access (`shellization_authority` lives on the atom-growth materialized report). Fixed in the generator.
  - Rerun passed and printed:
    - `pqs_status=available_diatomic_complete_core_shell_hamiltonian_handoff_payload`
    - `cr2_read_only_inspector_ready=true`
    - `cr2_solver_ready=false`
    - `white_lindsey_status=available_route_configured_diatomic_atom_growth_ham_payload`
- `julia --project=. -e 'using JLD2; ... readback ...'`
  - Passed; summary values listed above.
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
  - Passed: `load ok`
- `git diff --check`
  - Passed.
- `git diff --numstat -- src test`
  - Empty.
- `git status --ignored --short tmp/work/be2_wl_pqs_cr2_inspection_artifact`
  - Only the generator is modified; generated JLD2/TSV outputs are ignored.

Final `git status --short --branch`:
```text
## main...origin/main
 M tmp/work/be2_wl_pqs_cr2_inspection_artifact/generate_be2_wl_pqs_cr2_inspection_artifact.jl
```

Deletion/shrinkage report:
deleted: none
simplified: none in `src`/`test`; the generator now uses the compact atom-growth route/adapter seam instead of adapting seed/one-center WL paths
quarantined: old WL seed/one-center/debug paths remain out of the artifact; supplement/residual-GTO, QW atom-local HF, correction/EGOI/stationary/cusp, MWG/IDA, and solver/export surfaces remain explicitly unavailable or pending
not deleted because: this pass was generator-only and the existing source route-driver/test surfaces remain the live authority for WL route construction
exact remaining caller/blocker: downstream solver/export handoff is still blocked by `:missing_hfdmrg_density_density_contract`; MWG/IDA WL ham support remains blocked by `:pending_route_configured_diatomic_mwg_operator_support`

-- repo-doer@macmini
