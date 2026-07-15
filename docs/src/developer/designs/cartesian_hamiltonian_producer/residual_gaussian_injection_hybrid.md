# Default-Off Direct-G Residual Injection

Status: implemented internal compatibility facility under
`HP-RG-INJECT-FN-01`; disabled by default and not the current protected-main
construction target.

This page is the canonical contract for the historical direct-`G` injection
path. The registry owns permission, lifecycle, and exact source surfaces. This
page owns the numerical construction, default-off behavior, failure rules, and
its boundary with the current protected-localized architecture.

## Purpose

Ordinary residual Gaussian (RG) construction projects an original supplement
direction `a` out of an orthonormal terminal gausslet space `G`. If its
residual occupation is `lambda`, the normalized complement behaves like:

```text
r = (I - P_G) a / sqrt(lambda)
```

This is fragile as `lambda -> 0`: discarding the direction loses exact
supplement-span information, while retaining it creates a normalized tiny
complement that can behave like a ghost residual function.

Direct-`G` injection supplies a third fate. A near-gausslet supplement mode
replaces its represented direction inside `G` rather than becoming an RG/MWG
channel. The implementation is retained as an internal, default-off
compatibility and measurement surface. Current protected work instead replaces
directions over the compact main space `M = [G, R_compact]`; see
[Protected-localized basis convention](protected_localized_basis.md).

## Policy And Thresholds

The internal option is:

```text
residual_injection_cutoff = lambda_inj
```

Its contract is:

```text
lambda_inj <= 0    injection disabled
lambda_inj > 0     classify near-G modes for direct-G replacement
```

When enabled, `lambda_inj` must be at least the active
`residual_occupation_cutoff`. This avoids an accidental band that is neither
injected nor retained. No public/default policy is approved.

Keep these numerical decisions distinct:

```text
candidate_overlap_cutoff       supplement metric rank cleanup
residual_injection_cutoff      near-G injection classification
injected_overlap_cutoff        duplicate global injection cleanup
residual_occupation_cutoff     true RG retention after replacement
identity_atol                  final residual identity validation
```

Candidate metric rank is not residual occupation. None of these thresholds is
a kinetic or one-body spectral guard.

## Candidate Classification

For each physical owner, form the owner-local supplement overlap `S_AA` and
rank-clean it with the established absolute/relative Gram rule. In the cleaned
orthonormal candidate basis `A_tilde`, define:

```text
C = G' S A_tilde
M_res = I - C'C
M_res v_i = lambda_i v_i
y_i = A_tilde v_i
```

Classification applies to the principal modes `y_i`, not to raw contracted
GTO columns:

- `lambda_i <= lambda_inj`: provisional injected mode;
- `lambda_i > lambda_inj`: provisional true-residual candidate.

Candidate Gram cleanup and residual occupation answer different questions and
must not share a cutoff.

Injection is global because injected modes do not become owner-local MWG
channels. True residual selection remains owner-local:

1. classify stable candidate principal modes separately for each owner;
2. concatenate provisional injected modes;
3. globally orthonormalize and remove duplicate injected directions using the
   separate injected-overlap rank rule;
4. construct one global injected target `Y`;
5. residualize each owner's remaining candidates against the replacement
   sector;
6. apply owner-local true-RG occupation selection;
7. perform the ordinary final inter-owner residual merge.

Global injection cleanup is not permission for global residual selection.

## Replacement Geometry

For the global orthonormal injected target `Y`, form:

```text
B = G' S Y
```

Construct an orthonormal complement `Q_perp` in the coordinate space of
`G`:

```text
B' Q_perp = 0
Q_perp' Q_perp = I
F = [Y, G Q_perp]
```

`F` has the same dimension as `G`. Injection is replacement, not append.
Required guards are:

- `dim(Y) < dim(G)`;
- `B` has full injected rank and acceptable conditioning;
- `F' S F` is identity to tolerance;
- true residuals are orthogonal to `F`, not merely to the old `G`;
- duplicate cross-owner injected modes are removed explicitly;
- final `F-R` and `R-R` metric checks pass.

If the injected target is too large, rank deficient, or unstable under global
cleanup, stop. Do not silently discard target directions, append them, or turn
them into MWG residuals.

This compatibility helper constructs the replacement span. That is only the
first part of the full angular-style localized injection construction. A caller
claiming angular-style localization must additionally project the old localized
`G` seeds into `F` and symmetrically Lowdin-orthogonalize those projected seeds.
Do not treat `[Y, G Q_perp]` itself as evidence that localized representatives
were recovered. The current protected-localized architecture owns that complete
relocalization step.

## Operators And Interaction

With injection enabled:

- exact kinetic, moment, and per-center unit-nuclear operators are transformed
  into the in-memory `[F, R]` basis;
- remaining true RGs retain the existing exact augmented one-body transform;
- the injected base sector inherits the original gausslet-sector IDA
  interaction convention;
- residual MWG/IDA descriptors and channels apply only to true residuals.

The inherited IDA treatment is an explicit two-body approximation. Injected
directions do not receive residual-GTO/MWG descriptors. Exact one-body
transformation does not authorize a two-index interaction congruence or a new
density convention.

The current supplemented artifact writer rejects
`residual_injection_cutoff > 0`. Enabled injection is in-memory only; this
contract does not authorize artifact provenance or persistence.

## Current Implementation Boundary

The implementation is internal and default-off. With
`residual_injection_cutoff <= 0`, ordinary RG selection, exact augmented
operators, MWG/IDA interaction, endpoint values, and artifact behavior must
remain unchanged within roundoff.

Approved implementation owners are:

- `src/cartesian_residual_gaussians/residual_basis.jl`;
- `src/cartesian_residual_gaussians/augmented_operators.jl`;
- `src/cartesian_residual_gaussians/mwg_interaction.jl`;
- narrow compatibility plumbing in
  `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl`.

Compact implementation facts may include injection and rank thresholds,
injected dimension/owner counts, `B/Q_perp` transform authority, and the
existing true-residual `T_G/T_A` blocks. Do not add discarded spectra,
inventory payloads, status fields, or report frameworks to construction
objects.

There is no current production caller that should turn this path on as the
protected construction. The protected-localized builder uses the separate
compact-main geometry and operator contract.

## Validation And Failure Behavior

The source-backed lane was validated with:

- default-off H2 residual-GTO/MWG endpoint parity;
- ignored default-off Cr/Cr2 replay against ordinary production RG counts;
- ignored enabled replay of injected/true-RG counts, `B` rank/condition,
  `F' S F`, `F' S R`, `R' S R`, finite/symmetric exact one-body
  matrices, and low `K_RR/H1_RR` spectra;
- package load and `git diff --check`.

Stop rather than broaden the implementation if enabled injection requires:

- artifact/schema/provenance changes;
- driver/public API/default changes;
- source outside the registered owners;
- terminal-basis, raw-block, route, or shellification rewrites;
- persistent dense `nG x nG` workspace beyond the approved compact
  transform;
- global residual selection, automatic spectral pruning, solver/HF workflow,
  or Cr2-specific production behavior.

No committed large fixture or endpoint is part of this compatibility lane.

## Historical Measurement

`HP-RG-INJECT-AUDIT-01` is completed historical measurement authority. The
first Cr2 reconstruction found that direct-`G` injection was numerically
well-conditioned but did not remove the low two-owner residual sector:

```text
lambda_inj        injected   true RG count   min K_RR   min H1_RR
0                 0          132             0.428594   -7.349209
1.0e-4            38         100             0.445326   -7.061948
```

At `lambda_inj = 1.0e-4`, the reported `B` condition was about `1.002`
and the implicit `F' S F` error about `3.9e-12`. The audit's zero point was
not exactly production-equivalent, so these values are evidence, not an
endpoint baseline.

The durable interpretation is:

- direct injection addresses the singular-complement construction problem;
- it is not by itself a demonstrated Cr2 residual-sector fix;
- cutoff-only RG remains an unsuitable long-term construction principle;
- residual spectral safety is a separate diagnostic/design question.

Detailed chronology remains in the manager running log and repository history.

## Successor Contracts

Do not use this direct-`G` compatibility path as a substitute for:

- [Occupied-first injection geometry](occupied_first_injection.md), which owns
  mandatory occupied protection and optional supplement capture;
- [Protected-localized basis convention](protected_localized_basis.md), which
  owns replacement over `M = [G, R_compact]`, exact localized one-body
  operators, and inherited-site `Vee_L`;
- [Protected-localized artifact contract](protected_localized_artifact.md);
- [Retained-GTO local-product EGOI](retained_gto_egoi.md);
- [Protected-localized ladder bundles](protected_localized_ladder.md);
- [Protected additive atomic reference correction](protected_additive_reference_correction.md);
- [Screened Hartree residual-density formalism](screened_hartree_residual_density.md).

Rejected broad protected candidates never become MWG residuals. Direct
`C' V C` interaction rotation remains invalid under the protected-localized
contract.

## Do Not Confuse

- Candidate overlap rank is not residual occupation.
- Residual occupation is not integral weight or numerical rank.
- Injected directions are not residual Gaussians.
- Injection is replacement, not append.
- Direct-`G` injection is not compact-main protected injection.
- Global injected-subspace cleanup is not global residual selection.
- Exact one-body transformation is not an interaction rotation.
- An enabled in-memory path is not artifact or production-workflow authority.
