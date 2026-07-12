# Current Cartesian Hamiltonian Producer Status

This page owns a concise snapshot of implemented facilities, active work,
blockers, and next steps. It does not grant source authority or restate
subsystem mathematics.

For a task, authority requires agreement among:

1. the source-ID whitelist in `AGENTS.md`;
2. the assigned entry in `registry.md`;
3. the canonical subsystem contract linked by that entry;
4. any explicit boundary on this page.

Silence here is neutral. A conflict is not: fail closed and request a
docs-only reconciliation.

Normal producer startup is:

- [README](README.md) for orientation and task-specific navigation;
- this page for live status;
- [invariants](invariants.md) for cross-subsystem guardrails;
- only the assigned registry entry and linked subsystem contract.

Before numerical implementation, also consult
`docs/src/developer/algorithm_implementation_index.md`. Endpoint work must
follow the operational-facts and terminal due-diligence rules in `AGENTS.md`.
History, reviews, reports, implementation ledgers, and the manager running log
are evidence, not normal startup authority.

## Implemented Facilities

| Area | Current implemented state | Canonical contract |
| --- | --- | --- |
| Terminal basis | Disjoint owned terminal supports, support-local realization, structural cross-block overlap, blockwise exact operators | [Terminal basis and base assembly](terminal_basis_and_base_assembly.md), [common shell decomposition](common_terminal_shell_decomposition.md) |
| Diatomic shell corrections | Angular-z extension, neutral face products, compact PQS/WL thin slabs, aspect-aware PQS shared complete-shell `(q,q,L)` source modes, and private semantic per-shell source-q refinement | [Common shell decomposition](common_terminal_shell_decomposition.md), [aspect-aware source modes](pqs_complete_shell_aspect_source_modes.md), [semantic source-q overrides](pqs_semantic_shell_q_overrides.md) |
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
| `HP-RG-PROTECT-EGOI-*` | Measurement completed; retained-GTO helper/test approved pending | Implement only retained original `s1+s2`, local symmetric products, `M2`, and exact-zero disallowed `DeltaV`; the uncommitted `hamiltonian_corrections.jl` WIP is not accepted authority |
| `HP-RG-SPECTRAL-AUDIT-01` | Measurement-only | Characterize the surviving low residual-sector mode; no pruning or spectral guard is approved |
| `HP-RHO0-XPAIR-AUDIT-01` | Deferred measurement question | Exchange/direct pairing may be revisited on H/Be/Be2 only; it is not a current blocker or source lane |
| Documentation architecture | The reviewed candidate discrepancies, one follow-up live dependency, and pre-rehearsal tooling gaps are reconciled in a transition-bound candidate that still derives 150 execution IDs | Repeat the independent semantic/tooling rehearsal; no generated authority view or cutover is approved |

Other approved IDs may remain actionable even when not listed here. Read the
assigned registry entry. Completed retirement, superseded, rejected, and
historical audit IDs are not active work.

## Current Physics Target

The current consumer-facing target is the authorized controlled Cr2
screened-Hartree off/on fixed-density comparison. The numerical-complete
H2/Be2 source and validation gates have passed:

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
5. **Authority cutover remains unapproved.** The no-go findings in the durable
   [review result](reviews/authority_transition_rehearsal_2026-07-12.md) are
   reconciled: owned paths no longer absorb dependencies, missing paths and
   dependency records are restored, full-document whitelist context is
   checked, generated links are conservative, standalone previews are warned,
   and manifests bind candidate, prose, checker, snapshot, and Git inputs. The
   transition snapshot still says semantic parity requires manual review and
   the reconciled candidate is pending another independent rehearsal.

## Non-Negotiable Boundaries

- The terminal due-diligence report must be inspected before interpreting any
  endpoint, residual, injection, screened-Hartree, EGOI, Be2, or Cr2 result.
- Final working bases are orthonormal. Normal representation transfer uses
  only `S_BA = <B|A>` and `C_B = S_BA*C_A`; self-overlaps are diagnostics.
- Exact one-body operators transform through the actual represented basis.
  A two-index IDA/MWG density interaction is not transformed by `C' V C`.
- Residual directions are selected owner-locally and merged once. Rejected
  broad protected directions never become MWG residual channels.
- `Vnuc_G` remains Galerkin in screened Hartree. Represented occupied
  determinants define `P0/q0`; density fits define `E0`; fitted potentials
  are approximate fast `J0` evaluators.
- One Hamiltonian construction uses one resolved Coulomb expansion. Atomic
  packet RHF/density/potential roles remain separately qualified exceptions.
- No element-specific defaults, automatic mapping/Coulomb tuning, ECP,
  generalized-overlap final workflow, public solver, or Cr2 production path is
  implied by existing internal facilities.

## Canonical Reading

Read only the relevant contract:

- [Terminal and producer invariants](invariants.md)
- [Terminal basis and base assembly](terminal_basis_and_base_assembly.md)
- [Residual Gaussian domain](residual_gaussian_domain_module.md)
- [Numerical-complete residual basis](numerical_complete_residual_basis.md)
- [Protected-localized basis](protected_localized_basis.md)
- [Protected-localized artifacts](protected_localized_artifact.md)
- [Retained-GTO EGOI](retained_gto_egoi.md)
- [Protected ladders](protected_localized_ladder.md)
- [Atomic HF reference packets](atomic_hf_reference_packets.md)
- [Screened Hartree formalism](screened_hartree_residual_density.md)
- [External GTO import](external_gto_orbital_import.md)
- [Coulomb accuracy policy](coulomb_accuracy_policy.md)
- [Semantic per-shell PQS source-q overrides](pqs_semantic_shell_q_overrides.md)

Use [registry.md](registry.md) only for the assigned ID. Use
[README.md](README.md) for the broader task-specific document map.
