# Historical Rho0 Reference-Density Implementation Plan

Status: superseded historical stub. This file grants no active source,
measurement, artifact, solver, or workflow authority.

The original plan guided the transition from row-gauge diagnostics to fixed
`P0`, identified the missing exact mixed-Hartree seam, and staged the `GG`,
`GA/AA`, and protected exact-side implementation. That numerical
infrastructure is complete and now belongs to
[reference Hartree numerics](reference_hartree_numerics.md).

The correction-policy experiments that followed are completed or superseded.
The live correction contracts are:

- [Screened Hartree residual density](screened_hartree_residual_density.md);
- [Screened Hartree correction assembly](screened_hartree_correction_assembly.md);
- [Protected additive reference correction](protected_additive_reference_correction.md).

The compact historical fixed-`P0` account and deferred
`HP-RHO0-XPAIR-AUDIT-01` question remain in
[rho0 and reference-density correction history](rho0_reference_density_matrix.md).
`HP-RHO0-REFDENS-FN-01` and `HP-RHO0-REFDENS-ERI-01` remain unapproved
planning names.

Historical detail is preserved in manager running-log Passes 271-297 and Git
history. In particular, Passes 280-282 record review of this plan, while
commits `efaee93f6`, `daac231d0`, and `40a6f7e99` record the implemented
neutral Hartree slices.
