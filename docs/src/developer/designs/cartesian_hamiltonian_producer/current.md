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
| Diatomic shell corrections | Angular-z extension, neutral face products, compact PQS/WL thin slabs, aspect-aware PQS shared complete-shell `(q,q,L)` source modes, and private semantic per-shell source-q variation; refinement is implemented and coarsening is approved pending | [Common shell decomposition](common_terminal_shell_decomposition.md), [aspect-aware source modes](pqs_complete_shell_aspect_source_modes.md), [semantic source-q overrides](pqs_semantic_shell_q_overrides.md) |
| Base producer | Exported facade plus blockwise exact H1, localized IDA, and direct `CartesianIDAHamiltonian` construction for the implemented atom/diatomic scope | [Terminal basis and base assembly](terminal_basis_and_base_assembly.md), [R1 base producer](r1_public_base_producer.md), [one-center atoms](r1_one_center_base_atoms.md) |
| Ordinary artifacts | Unchanged matrix payload and reader, plus facade-written native-order labels, recipe truth, and the implemented subset of construction-native source provenance | [Artifact manifest](cartesian_hamiltonian_artifact_manifest.md) |
| Composition and driver | PQS/WL, base/supplemented composition through shared producer boundaries; canonical artifact-producing driver; terminal inventory and due-diligence reporting | [Composition contract](nesting_supplement_composition_plan.md), [driver workflow](cartesian_driver_usability_workflow.md), [due diligence](terminal_shellification_due_diligence.md) |
| Mapping and source span | Expert `s_factor` and opt-in mapped-COMX source spans with provenance; defaults remain unchanged | [Mapping s_factor](pqs_mapping_s_factor.md), [mapped COMX](mapped_comx_source_span.md) |
| Coulomb policy | One producer-wide expansion reaches parent/PGDG, base IDA, residual-GTO, and MWG. Compact45 and high135 are implemented | [Coulomb accuracy](coulomb_accuracy_policy.md) |
| Residual Gaussians | Owner-local residual selection, one final merge, exact augmented one-body operators, final-basis MWG/IDA, current `1e-6` production cutoff, and opt-in numerical-complete `[G,R_num]` additive composition | [Residual Gaussian domain](residual_gaussian_domain_module.md), [numerical-complete basis](numerical_complete_residual_basis.md), [orthogonality/cutoff](residual_gaussian_orthogonality_robustness.md) |
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
| `HP-PQS-COULOMB-ACCURACY-*` | Standard60 and canonical-driver exposure approved, not implemented | Add the fixed audited K60 resolver and fingerprint provenance; accept compact/standard/high in facade and driver without changing the compact default |
| `HP-PQS-SHELLQ-OVERRIDE-*` | Refinement implemented; coarsening amendment approved, not implemented | Reuse the existing post-shellification override path; accept non-Boolean integer `source_q >= 3` except redundant `source_q == route_q`, then validate route-7 to 6/5 without changing parent support or ownership |
| `HP-RG-PROTECT-EGOI-*` | Measurement completed; retained-GTO helper/test approved pending | Implement only retained original `s1+s2`, local symmetric products, `M2`, and exact-zero disallowed `DeltaV`; the uncommitted `hamiltonian_corrections.jl` WIP is not accepted authority |
| `HP-RG-SPECTRAL-AUDIT-01` | Measurement-only | Characterize the surviving low residual-sector mode; no pruning or spectral guard is approved |
| `HP-RHO0-XPAIR-AUDIT-01` | Deferred measurement question | Exchange/direct pairing may be revisited on H/Be/Be2 only; it is not a current blocker or source lane |
| `HP-RETIRE-CARRIED-SPACE-*` | Retirement approved, deletion pending | Pass 409 removed the QW/high-order cluster; now delete its orphaned 266-line internal carried-space adapter and sole include without replacement |
| Existing execution IDs | Post-cutover conformance remediation | Reconcile the bounded discrepancies recorded in the [2026-07-12 execution audit](reviews/execution_conformance_audit_2026-07-12.md), beginning with fail-fast correctness and misleading completed-test claims |

Other approved IDs may remain actionable even when not listed here. Read the
assigned registry entry. Completed retirement, superseded, rejected, and
historical audit IDs are not active work.

## Current Physics Target

The current consumer-facing target is the authorized controlled Cr2
screened-Hartree off/on fixed-density comparison under the currently approved
fitted-cloud reference convention. The numerical-complete H2/Be2 source and
validation gates have passed:

- one high-accuracy `M=[G,R_num]` Hamiltonian basis that preserves `G`;
- one imported external-GTO occupied start represented by the source-backed
  cross-overlap infrastructure;
- identical starting state for screened and unscreened runs;
- consumer-owned solver continuation, spin/occupation diagnostics, and
  interpretation.

This is a consumer measurement, not permission for a Cr2-specific producer
branch, committed Cr2 endpoint, solver API, corrected artifact, or production
energy claim.

## Current Blockers And Follow-Ups

1. **Standard Coulomb implementation.** The analytic K60 preset and artifact
   fingerprint are approved but not source-backed. The controlled Cr2
   screened comparison stays `:high` and must not be changed mid-comparison.
2. **Retained-GTO EGOI helper.** The protected target/mask convention is
   accepted, but committed source does not yet implement the helper.
3. **Residual spectral interpretation.** Tightening the RG cutoff removed
   marginal residuals but did not remove the measured low two-owner mode.
   Injection and cutoff changes are not substitutes for a separately approved
   safety policy.
4. **Protected atoms and counterpoise.** One-center protected compactness,
   separated kinetic/unit-nuclear persistence, and counterpoise sidecars remain
   separate future designs.
5. **Screened-reference exactness policy.** Current source uses determinant
   orbitals for `P0/q0`, the density fit for `E0`, and the fitted potential for
   approximate `J0`, with consistency error reported. Determinant-exact
   `J0/E0` would be a new scientific amendment, not a conformance repair; do
   not interpret a run under the other convention without deciding this first.

Durable numerical and workflow guardrails live in
[invariants](invariants.md); the task-specific contract map lives in
[README](README.md). Endpoint work also follows the operational-facts and
terminal due-diligence rules in `AGENTS.md`.
