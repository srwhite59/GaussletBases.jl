# Occupied-First Injection Geometry

Status: implemented internal subsystem under
`HP-RG-OCC-FIRST-INJECT-FN-01` and
`HP-RG-OCC-FIRST-INJECT-TEST-01`.

This page is the canonical contract for occupied-first global injection
geometry. The registry owns permission, lifecycle, and file surfaces; this
page owns the mandatory occupied subspace, capture mathematics, diagnostics,
failure behavior, and current consumer boundary.

## Purpose And Scope

The helper makes an explicitly supplied Hartree-Fock (HF) occupied subspace
mandatory before optional supplement injection is selected. This prevents the
reference determinant from depending on hard-coded orbital labels or an
optional capture cutoff.

The implemented operation is:

```text
occupied HF subspace: mandatory
remaining supplement: selected by represented-span capture
weak-capture directions: reported and rejected
rejected directions: never converted into MWG residual channels
```

This is internal geometry and selection logic. It does not construct a
Hamiltonian, write an artifact, or run HF.

## Inputs And Provenance

`occupied_first_injection_geometry(...)` consumes:

- orthonormal base/final gausslet dimension `nG`;
- mixed overlap `X_GA = <G|A>`;
- supplement metric `S_AA`;
- occupied coefficients `Y_occ` in supplement coordinates;
- an optional-injection capture cutoff;
- optional labels, owners, channels, and provenance for diagnostics.

`Y_occ` must come from an explicit atomic reference packet or external-GTO
import contract. The supplying owner must validate basis order, overlap, and
fingerprint before calling the geometry helper. The helper accepts the
already-resolved coefficient block and carries supplied provenance; it does
not infer occupied orbitals from labels.

Require:

```text
Y_occ' * S_AA * Y_occ = I
0 < size(Y_occ, 2) < nG
```

All supplied matrices must be finite. Occupations remain reference provenance,
not geometry input. The optional cutoff applies only to nonoccupied supplement
directions.

## Physical Mixed-Overlap Validation

The mixed overlap must describe a physical projection of the supplement into
an orthonormal base span. Define the complement metric:

```text
R_AA = S_AA - X_GA' * X_GA
```

The helper derives a private numerical `capture_tol` from the active overlap
tolerances and requires:

```text
lambda_min(R_AA) >= -capture_tol
```

A materially negative complement means `X_GA` and `S_AA` do not describe one
consistent overlap geometry. It is a hard failure, not a reason to clamp an
eigenvalue or weaken the injection cutoff.

## Pre-Inclusion Occupied Capture

Before making the occupied subspace mandatory, report its capture by the
original base span:

```text
s_occ = svdvals(X_GA * Y_occ)
```

The implementation records `s_occ` and
`occupied_base_capture_min = minimum(abs2.(s_occ))`. This is the diagnostic
that reveals whether mandatory inclusion will make a small or large change.
Poor pre-inclusion capture is reported; it is not confused with recovery after
the subspace has been forcibly included.

## Mandatory Inclusion And Recovery

The occupied supplement directions form the mandatory injected block. The
corresponding replacement base sector is the occupied target plus the part of
the old base span orthogonal to it:

```text
F_occ = Y_occ + (G intersect Y_occ_perp)
```

The helper obtains a stable supplement-metric orthonormal mandatory block and
its projection into `G`, then constructs the base-coordinate complement. The
mandatory block cannot be removed by the optional cutoff.

After inclusion, compute separate recovery singular values from the overlap of
the mandatory block with the original `Y_occ`. Require roundoff recovery and
report:

```text
occupied_recovery_after_mandatory_inclusion_singular_values
occupied_recovery_after_mandatory_inclusion_loss
```

Pre-inclusion capture and post-inclusion recovery are different objects. Do
not restore the deleted `weakest_occupied_capture` alias, which hid that
distinction by reporting a quantity near one by construction.

## Optional Supplement Capture

Rank-clean and orthonormalize the full supplement metric, remove the mandatory
occupied subspace, and evaluate the remaining supplement directions against
the represented occupied-first span. Diagonalize both the full and mandatory-
complement capture operators.

Every raw capture eigenvalue must satisfy:

```text
-capture_tol <= lambda <= 1 + capture_tol
```

The implementation records unclamped full/complement ranges. A material value
outside the projector range is a physical-contract failure. Tolerance-sized
roundoff may be normalized for presentation only; it must not hide malformed
geometry.

Select optional directions by:

```text
keep lambda >= optional_capture_cutoff
reject lambda < optional_capture_cutoff
```

The current default cutoff is `0.99`. It is an optional supplement-selection
policy, never protection for `Y_occ`. Kept and rejected directions are
reported with channel/owner labels when supplied. Rejected weak-capture
directions do not become residual-GTO or multiwavelet-Gaussian (MWG) channels.

## Outputs And Diagnostics

The helper returns mandatory, optional, and combined supplement coefficient
blocks, their base projections, capture spectra, kept/rejected indices, and a
compact diagnostic record. Required diagnostics include:

- supplied provenance and occupied metric error;
- pre-inclusion occupied singular values and squared minimum capture;
- post-inclusion recovery singular values and loss;
- `capture_tol` and complement-metric minimum eigenvalue;
- unclamped full and complement capture ranges;
- full and optional-complement capture spectra;
- cutoff, kept/rejected counts, and weakest kept/strongest rejected capture;
- final and mandatory ranks;
- kept/rejected direction summaries;
- explicit confirmation that rejected directions were not made MWG residuals.

## Failure Behavior

Stop when dimensions disagree, data are nonfinite, `Y_occ` is empty or not
orthonormal in `S_AA`, the mandatory rank reaches the base dimension, the
mixed-overlap complement is materially negative, capture eigenvalues leave
their physical range, mandatory recovery is not roundoff accurate, or the
combined injected projection is rank deficient.

Do not repair failure by dropping an occupied direction, selecting it only by
cutoff, relabeling orbitals, flooring capture eigenvalues, converting a rejected
direction into MWG, or silently changing packet/import provenance.

## Current Consumer Boundary

`occupied_first_injection_geometry(...)` is implemented and tested, but it has
no protected-localized builder caller. The current protected builder uses
staged protected-original geometry over:

```text
M = [G, R_compact]
```

The occupied-first helper operates on the base `G` geometry and is not a
direct call-site substitute for that staged construction. Combining mandatory
occupied protection with one shared `R_compact`, placed atomic packets, and
additive molecular `P0/J0/E0` belongs to
`HP-RG-PROTECT-ADDREF-*` and
[Protected additive atomic reference correction](protected_additive_reference_correction.md).

## Ownership And Tests

Implementation and test ownership is recorded in the compact registry entry.
The implementation lives in
`src/cartesian_residual_gaussians/residual_basis.jl` and may consume
already-owned packet/import coefficient data without changing those owners.

Validation uses:

- `test/misc/runtests.jl` for the tiny pre/post and malformed-capture contract;
- `test/nested/cartesian_occupied_first_injection_runtests.jl` for real Be/Ne
  cc-pV5Z, `lmax = 1`, `ns = 5` packet-driven PQS overlaps and terminal due
  diligence.

The historical measurement audit and implementation evidence remain in the
manager running log, especially Passes 323, 324, 344, and 345.

## Dependencies

- [Atomic HF reference packets](atomic_hf_reference_packets.md) for packet
  determinant identity and convergence;
- [External GTO orbital import](external_gto_orbital_import.md) for imported
  occupied coefficient identity;
- [Residual Gaussian domain module](residual_gaussian_domain_module.md) for
  general residual/injection ownership and MWG boundaries.

## Exclusions

This contract does not approve:

- direct substitution into staged protected-original geometry;
- protected additive-reference implementation beyond its separate authority;
- screened-Hartree correction changes, EGOI, exchange, or row-gauge rho0/P0;
- shell-local injection or fake-RDM hierarchy changes;
- public driver/API/defaults or automatic physics policy;
- artifact schema, writer, reader, or provenance changes;
- solver workflow, endpoint claims, or Cr/Cr2 production work.
