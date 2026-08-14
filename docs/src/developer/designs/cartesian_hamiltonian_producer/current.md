# Current Cartesian Hamiltonian Producer Status

This page owns live implementation status, active work, blockers, and next
steps. It does not grant authority or restate subsystem contracts.

Start with [README](README.md), this page, and [invariants](invariants.md), then
read only the assigned generated registry entry and linked canonical contract.
`authority.toml` remains the record-level source. Silence here is neutral; an
explicit conflict fails closed. History, reviews, reports, and the manager log
are evidence rather than startup authority.

The broad documentation reorganization and machine-authority cutover are
complete. Remaining documentation work is maintenance attached to concrete
source/contract findings, not another migration campaign.

## Implemented Facilities

| Area | Current implemented state | Canonical contract |
| --- | --- | --- |
| Terminal basis | Disjoint owned terminal supports, support-local realization, structural cross-block overlap, blockwise exact operators | [Terminal basis and base assembly](terminal_basis_and_base_assembly.md), [common shell decomposition](common_terminal_shell_decomposition.md) |
| Diatomic shell corrections | Common physical shells, angular-z extension, neutral face products, compact PQS/WL thin slabs, matched aspect-aware PQS/WL shared-shell source modes, and private semantic per-shell source-q refinement/coarsening | [Common shell decomposition](common_terminal_shell_decomposition.md), [matched aspect source modes](pqs_complete_shell_aspect_source_modes.md), [semantic source-q overrides](pqs_semantic_shell_q_overrides.md) |
| Base producer | Exported facade plus blockwise exact H1, localized IDA, and direct `CartesianIDAHamiltonian` construction for the implemented atom/diatomic scope | [Terminal basis and base assembly](terminal_basis_and_base_assembly.md), [R1 base producer](r1_public_base_producer.md), [one-center atoms](r1_one_center_base_atoms.md) |
| Ordinary artifacts | Unchanged matrix payload and reader, plus facade-written native-order labels, recipe truth, and the implemented subset of construction-native source provenance | [Artifact manifest](cartesian_hamiltonian_artifact_manifest.md) |
| Composition and driver | PQS/WL, base/supplemented composition through shared producer boundaries; canonical artifact-producing driver; terminal inventory and due-diligence reporting | [Composition contract](nesting_supplement_composition_plan.md), [driver workflow](cartesian_driver_usability_workflow.md), [due diligence](terminal_shellification_due_diligence.md) |
| Mapping and source span | Expert `s_factor` and opt-in mapped-COMX source spans with provenance; defaults remain unchanged | [Mapping s_factor](pqs_mapping_s_factor.md), [mapped COMX](mapped_comx_source_span.md) |
| Coulomb policy | One producer-wide expansion reaches parent/PGDG, base IDA, residual-GTO, and MWG. Compact45 and high135 are implemented | [Coulomb accuracy](coulomb_accuracy_policy.md) |
| Residual Gaussians | Owner-local residual selection, one final merge, exact augmented one-body operators, final-basis MWG/IDA, current `1e-6` production cutoff, and opt-in numerical-complete `[G,R_num]` additive composition | [Residual Gaussian domain](residual_gaussian_domain_module.md), [numerical-complete basis](numerical_complete_residual_basis.md), [orthogonality/cutoff](residual_gaussian_orthogonality_robustness.md) |
| Parent-backed functions | Root-exported expert PRF composition is accepted with exact descriptor-to-PRF source binding. Internal injection, exact one-body, numerical-complete residual, separated interaction, screened-Hartree delegation, and parent-IDA mechanics remain its numerical owners | [Parent residual functions](parent_residual_functions.md), [parent-backed injected composition](parent_backed_injected_composition.md) |
| Direct-G injection | Default-off in-memory compatibility path; ordinary behavior is invariant and enabled artifacts remain unsupported | [Direct-G injection](residual_gaussian_injection_hybrid.md) |
| Protected-localized basis | Compact-main protected replacement, exact localized one-body matrices, inherited-site `Vee_L`; direct `C' V C` is rejected | [Protected-localized basis](protected_localized_basis.md) |
| Protected persistence | Opt-in protected Hamiltonian artifacts with native locality metadata, plus same-parent ladder bundles and exact cross-overlap transfer | [Protected artifact](protected_localized_artifact.md), [protected ladder](protected_localized_ladder.md) |
| Reference infrastructure | Converged atomic HF packets, neutral reference-Hartree numerics, in-memory screened direct-Hartree correction, and additive protected molecular references | [Atomic packets](atomic_hf_reference_packets.md), [reference Hartree numerics](reference_hartree_numerics.md), [screened Hartree](screened_hartree_correction_assembly.md), [additive references](protected_additive_reference_correction.md) |
| Representation transfer | External-GTO import by cross overlap and standalone protected native `S_LG` sidecars; source self-overlap remains validation-only | [External GTO import](external_gto_orbital_import.md) |

These are bounded internal/repo facilities. “Implemented” does not imply a
public export, solver workflow, broad molecular support, or production Cr2
claim.

## Active And Pending Work

| Lane | State | Exact next boundary |
| --- | --- | --- |
| `HP-PQS-ASPECTSHELL-*` | Implemented/completed maintenance; eligible PQS/WL shells share `(ns,ns,L)` and equal aggregate dimensions | Preserve parent/PQS parity and fail on missing shape, axis-count, aggregate-column, or Gram inconsistencies |
| `HP-PQS-PAPER-H2-DRV-FN-01` | Five-row fixed-state measurement replayed with matched bare/supplemented dimensions | Rerun the external same-density Coulomb oracle against the new WL fingerprints before interpreting method accuracy |
| `HP-R1-ESECTOR-*` | Explicit charged sectors implemented; maintenance | Preserve exact basis/operator independence from `nup`/`ndn`, positive nonzero sectors, neutral compatibility, supplemented parity, and charged artifact readback |
| `HP-PQS-PRF-CONSUMER-*` | Implemented/completed maintenance | Preserve exact descriptor-to-PRF source binding, explicit consumer targets, additive/fixed-span parity, and category-owned unscreened Hamiltonian assembly |
| `HP-PQS-COULOMB-ACCURACY-*` | Standard60 and canonical-driver exposure approved, not implemented | Add the fixed audited K60 resolver and fingerprint provenance; accept compact/standard/high in facade and driver without changing the compact default |
| `HP-REP-MIXDENS-HARTREE-*` | Direct represented-density Hartree producer implemented; bounded repository validation passed | Run the separate Cr2 complete-field, external-sector, energy, derivative, and resource acceptance; no fitted route or production endpoint claim |
| `HP-SLICE-HCHAIN-*` | Six-export producer implemented; validation completed; maintenance | Preserve the O(1), zero-allocation structural-bandwidth query, exact represented off-band zeros, compact operator physics, and existing consumer contract; no further export, field, solver, or dense matrix |
| `HP-RG-PROTECT-EGOI-*` | Measurement completed; specialized retained-GTO helper/test deferred with no execution grant | Preserve the archived experiment only; any new helper requires docs-only reactivation tied to a current physics target |
| `HP-RG-SPECTRAL-AUDIT-01` | Measurement-only | Characterize the surviving low residual-sector mode; no pruning or spectral guard is approved |
| `HP-RHO0-XPAIR-AUDIT-01` | Deferred measurement question | Exchange/direct pairing may be revisited on H/Be/Be2 only; it is not a current blocker or source lane |
| Existing execution IDs | Post-cutover conformance remediation | Reconcile the bounded discrepancies recorded in the [2026-07-12 execution audit](reviews/execution_conformance_audit_2026-07-12.md), beginning with fail-fast correctness and misleading completed-test claims |

Other approved IDs may remain actionable even when not listed here. Read the
assigned registry entry. Completed retirement, superseded, rejected, and
historical audit IDs are not active work.

## Current Physics Target

The clean `R=2` full-parent/PQS/White-Lindsey replay now uses one parent,
common physical shells, `275` core functions, `50` slab functions, and `960`
complete-shell functions in each terminal family. Bare and supplemented
dimensions are matched at `1285/1303`. The parent and both PQS rows are
exactly unchanged from the pre-correction evidence; WL rows were rebuilt from
the corrected per-axis inner counts.

The fixed-state interaction extension remains implemented at clean commit
`7d2b6dc61`. It holds the common parent H2+ orbital fixed after separate
projection into the parent, bare, and supplemented spaces; uses matrix-free
parent high135 IDA, terminal/site IDA, and the established residual MWG
blocks; and reports the density-density energy decomposition without RHF/SCF
or artifacts. The matched repo replay is complete, but external oracle
interpretation must use the corrected WL state fingerprints before a method
accuracy claim. The `padding=20.0`, `tail_spacing=2.0` combination, other H2
endpoints, supplement or cutoff scans, capture at other campaign points, `q`
ladders, and geometry curves remain forbidden. This changes no canonical
driver, public facade, producer default, artifact schema family, or solver API.

## Current Blockers And Follow-Ups

1. **Matched paper-oracle interpretation.** Repository construction and the
   clean five-row replay are complete. The external same-density Coulomb
   oracle must be rerun against the corrected WL fingerprints before method
   accuracy or curve evidence is interpreted.
2. **Standard Coulomb implementation.** The analytic K60 preset and artifact
   fingerprint are approved but not source-backed. The controlled Cr2
   screened comparison stays `:high` and must not be changed mid-comparison.
3. **Parent-backed consumer studies.** The compact in-memory expert PRF API is
   accepted with exact descriptor-to-PRF source binding. Hooke and other
   consumers must still construct and justify their own parent-space targets,
   beginning with the proposed one-center Be `1s/2s` study.
   Transition-density exchange, exact PRF-GTO
   interactions, automatic selection, artifacts, and endpoint acceptance
   remain unapproved.
4. **Deferred retained-GTO EGOI.** Generic matrix-level EGOI remains
   implemented. The specialized retained-GTO adapter is archived rather than
   active; there is no source/test blocker or execution lane.
5. **Residual spectral interpretation.** Tightening the RG cutoff removed
   marginal residuals but did not remove the measured low two-owner mode.
   Injection and cutoff changes are not substitutes for a separately approved
   safety policy.
6. **Protected atoms and counterpoise.** One-center protected compactness,
   separated kinetic/unit-nuclear persistence, and counterpoise sidecars remain
   separate future designs.
7. **Screened-reference exactness policy.** Current source uses determinant
   orbitals for `P0/q0`, the density fit for `E0`, and the fitted potential for
   approximate `J0`, with consistency error reported. Determinant-exact
   `J0/E0` would be a new scientific amendment, not a conformance repair; do
   not interpret a run under the other convention without deciding this first.
8. **Complete molecular represented Hartree.** REQ-084 correctly stopped before
   molecular-full Cr2 interpretation. The direct internal producer is approved,
   but source and its complete-field certificate are pending. The first
   implementation preflight identified the absent nonmaterializing contracted
   pair-product action; result wrappers alone are incomplete scaffolding. An
   AO-projected field, occupied action, random probes, IDA transition-product
   substitute, or auxiliary fit cannot close this blocker.

Durable numerical and workflow guardrails live in
[invariants](invariants.md); the task-specific contract map lives in
[README](README.md). Endpoint work also follows the operational-facts and
terminal due-diligence rules in `AGENTS.md`.
