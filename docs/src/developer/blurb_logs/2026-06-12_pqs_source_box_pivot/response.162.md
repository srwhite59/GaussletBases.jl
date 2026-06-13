Pass 162 response - choose next downstream Hamiltonian contract

Role: repo-doer@macmini

Task type: no-edit design/decision audit. I edited no files, ran no Julia
tests, and ran no CR2/HFDMRG commands.

Recommended next contract option:

Choose option 2 first: a CR2 read-only inspector contract, built on the current
private Be2/PQS handoff and thin consumer contract.

Why:

- The medium-term goal is CR2 comparing Be2 WL and PQS. A read-only CR2
  inspector directly supports that comparison by exposing dimensions, H1
  spectrum/readiness, nuclear metadata, electron/spin convention, two-body
  representation labels, and explicit downstream blockers.
- The current handoff is already enough for inspection. It is not enough to
  claim HFDMRG density-density solver readiness because final-space `V::N x N`
  and Fock/energy conventions are not reviewed.
- This avoids adding another payload layer. The implementation should add a
  compact inspectability summary/contract view from the existing handoff and
  consumer, not a new route object.

Options rejected/deferred:

- HFDMRG density-density `H,V` contract: defer. It is the likely first solver
  contract, but it requires a reviewed final-space two-body convention before
  setting `hfdmrg_density_density_ready = true`. The current
  `:pre_final_density_interaction` handoff is inspectable, not solver-ready.
- HamV6/sliced-integral/export contract: defer. No `V6`, `Vblocks`, dense
  four-index object, export format, or artifact behavior exists.
- More cleanup only: do not make it the sole next target, but the next
  implementation must pay for any CR2 inspector edit by deleting remaining
  duplicated handoff/nonclaim surfaces under the line-negative rule.

Exact implementation seam for the next pass:

- Stay in `src/pqs_source_box_diatomic_complete_core_shell.jl`.
- Do not add a new payload type.
- Add or derive a compact CR2 read-only inspection view from the existing
  `diatomic_complete_core_shell_hamiltonian_consumer_contract_payload` and its
  `source_handoff`.
- Expected compact readiness vocabulary:
  - `cr2_read_only_inspector_ready = true`
  - `cr2_ready = false`
  - `hfdmrg_density_density_ready = false`
  - `hfdmrg_sliced_ready = false`
  - `hamv6_export_ready = false`
  - `exports_materialized = false`
  - `artifacts_materialized = false`
  - `public_api = false`
- Keep overall Hamiltonian readiness blocked on:
  `:missing_hfdmrg_density_density_contract`.
- Add `:cr2_read_only_inspection_contract` to available objects only if this
  can be done without creating a new payload cloud.
- Keep `:cr2_handoff_format` or a more precise
  `:cr2_solver_handoff_format` blocked for actual CR2 bundle/export behavior.

Line-budget strategy for the next implementation:

- The next source/test pass must be net-line-negative across `src` + `test`.
- Pay for the CR2 inspector by deleting remaining duplicated nonclaim flags in
  the handoff layer. Current repeated handoff false flags are still visible in:
  - conventions: `src/pqs_source_box_diatomic_complete_core_shell.jl:3117-3120`
  - handoff summary: `3170-3178`
  - handoff metadata: `3193-3199`
- Reuse `_pqs_source_box_route_driver_diatomic_downstream_ham_nonclaims()` for
  those handoff summaries where names match, or delete handoff-level duplicate
  flags that are not tested and are now represented by consumer/readiness.
- Test changes should replace copied-field assertions with one compact
  assertion that the CR2 read-only inspector contract is available while CR2
  solver/export readiness remains false.
- Proposed budget target: at least 25 deleted lines, no more than 20 added
  lines, net negative by at least 5 lines. If the CR2 view needs more code than
  that, stop for manager exception or write `ATTENTION.md`.

Physics/numerical convention that must be verified before saying ready:

- For CR2 read-only inspection: only label and compare conventions; do not
  claim solver readiness. The inspector must state that the two-body object is
  `:pre_final_density_interaction`, with density gauge
  `:pre_final_localized_positive_weight` and raw pair convention
  `:raw_numerator`.
- Before HFDMRG readiness: verify whether the pre-final density-interaction
  representation can be transformed or interpreted as HFDMRG's final-space
  density-density `V::N x N`, including direct/exchange and energy conventions.
- Before CR2 production readiness: choose an actual CR2 handoff format and
  seed-orbital/export contract.

Expected focused validation for the next pass:

```text
julia --project=. test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
git diff --numstat -- src test
```

No CR2/HFDMRG runs should be needed for the first read-only inspector contract.

Git status:

```text
git status --short --branch
## main...origin/main
```

Deletion/shrinkage report:

- deleted: none in this audit
- simplified: none in this audit
- quarantined: none in this audit
- not deleted because: audit only
- exact remaining caller/blocker:
  active readiness remains blocked on
  `:missing_hfdmrg_density_density_contract`; next pass should add only a
  read-only CR2 inspection contract/view if it can be paid for line-negatively,
  and must keep CR2 solver/export, HFDMRG, HamV6, dense `Vee`, final-space `V`,
  `V6`/`Vblocks`, H1/J, and RHF/SCF unavailable

-- repo-doer@macmini
