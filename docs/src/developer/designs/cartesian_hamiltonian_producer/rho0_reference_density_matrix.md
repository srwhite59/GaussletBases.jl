# Rho0 And Reference-Density Correction History

Status: compact historical account. The row-gauge, fixed-`P0`, mixed-seam,
approximate-Fock, old-anchor, and corrected-Hamiltonian audits are completed or
superseded evidence, not active source lanes.

## Current Ownership

The durable contracts now live in separate pages:

- [Reference Hartree numerics](reference_hartree_numerics.md) owns exact
  neutral `GG/GA/AA` Hartree blocks and protected fixed/localized transforms;
- [Screened Hartree residual density](screened_hartree_residual_density.md)
  owns the live `Delta_J0/C` physics;
- [Screened Hartree correction assembly](screened_hartree_correction_assembly.md)
  owns its implemented API;
- [Protected additive reference correction](protected_additive_reference_correction.md)
  owns additive molecular `P0/J0/E0` composition.

This page grants no correction source authority.

## Durable Historical Result

The original row-gauge line was rejected because

```text
(J*w)/w
diag(J)
u_direct
```

are different objects. Row action, a matrix diagonal, and the derivative of a
density-dependent Hartree approximation cannot substitute for one another.
`HP-RG-RHO0-GAL-AUDIT-01` is therefore completed historical evidence.

The successor audits used a fixed represented density matrix `P0` and the
affine matching equations

```text
Delta_F0_sigma = F_exact0_sigma[P0] - F_app0_sigma[P0]

C0 =
    E_exact0[P0]
  - E_app0[P0]
  - sum_sigma Tr(P0_sigma * Delta_F0_sigma).
```

Those audits established useful representation and exact-Hartree facts, but
did not produce the live correction policy. The exact mixed-Hartree seam was
successfully extracted into the neutral implemented infrastructure governed by
[reference Hartree numerics](reference_hartree_numerics.md).

Lifecycle summary:

- `HP-RHO0-REFDENS-AUDIT-01`: completed/superseded fixed-`P0` measurement;
- `HP-RHO0-REFDENS-MIXH-AUDIT-01`: completed mixed-Hartree seam audit;
- `HP-RHO0-FAPP-AUDIT-01`: completed approximate-Fock seam audit;
- `HP-RHO0-CORR-AUDIT-01`: completed/superseded corrected-Hamiltonian audit;
- `HP-RHO0-ANCHOR-FN-01` / `HP-RHO0-ANCHOR-TEST-01`: superseded, with no
  authority for new work;
- `HP-RHO0-FAPP-FN-01` / `HP-RHO0-FAPP-TEST-01`: implemented, caller-free,
  and dormant retirement candidates;
- `HP-RHO0-JANCHOR-FN-01` / `HP-RHO0-JANCHOR-TEST-01`: source-backed but
  superseded in use by
  `src/cartesian_reference_density/screened_hartree_correction.jl`;
  they are dormant retirement candidates.

The dormant helpers remain source facts until a separate deletion pass. This
page grants them no new caller or source authority.

Historical evidence is preserved in manager running-log Passes 271-297 and in
Git history, including source commits `8a7fcc00a`, `fce74302c`, and
`16834e974`.

## Deferred XPAIR Question

`HP-RHO0-XPAIR-AUDIT-01` remains an approved but deferred measurement
question. It is not a current blocker, source lane, or instruction to resume
the old correction program.

If explicitly reactivated by design-manager, it is limited to ignored
H/Be/Be2 measurements comparing direct-Hartree correction and inherited
approximate exchange-like contributions, with exact or supplement-space
exchange diagnostics where feasible. It does not authorize tracked source,
exact exchange implementation, corrected production Hamiltonians, artifacts,
public workflow, solver integration, Cr/Cr2, or use of a direct-only
diagnostic as physical Hamiltonian behavior.

## Unapproved Candidate IDs

`HP-RHO0-REFDENS-FN-01` and `HP-RHO0-REFDENS-ERI-01` remain unapproved
historical planning names. They are not active-status vocabulary, do not
belong in the AGENTS source whitelist, and grant no future source or numerical
authority.

## Historical Navigation

The full decision trail remains available without duplicating it here:

- manager running-log Passes 271-280: row-gauge retirement, fixed-`P0`, and
  mixed-seam audits;
- Passes 283-288: neutral exact Hartree implementation;
- Passes 289-297: approximate-Fock, anchor, corrected-Hamiltonian, and XPAIR
  investigation;
- Git history for the authority and implementation commits.

The former long implementation memo is retained only as a
[historical stub](rho0_reference_density_implementation_plan.md).
