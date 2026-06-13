Pass 115 manager review

Accepted as a local ignored damping sweep.

Simple density damping improved the density residual but did not converge the
compact route-owned private RHF SCF probe within 100 iterations for any tested
alpha:

- alpha 1.0: density change `4.679e-5`
- alpha 0.75: density change `3.464e-5`
- alpha 0.5: density change `2.273e-5`
- alpha 0.25: density change `1.116e-5`
- alpha 0.1: density change `5.319e-6`

Total energy was monotone decreasing and there was no obvious two-cycle, but
the density tolerance remains missed by orders of magnitude. This is not enough
to justify adding production/private damping controls.

Decision:

- Do not add damping/mixing to production code yet.
- Do not wire RHF into the route driver.
- Next pass should be a residual/Fock convention audit, focused on density
  convention, factor-of-two choices, and whether the energy/Fock residual is
  being measured for the same density used to build the Fock.

Validation/status:

- Local ignored sweep elapsed: `109.971126833` seconds.
- Reported git status clean and even with origin/main, apart from ignored local
  `tmp/work` artifacts.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: none; this was evidence gathering.
- quarantined: ignored sweep artifacts remain local diagnostic evidence.
- not deleted because: sweep evidence is needed for the residual/convention
  decision.
- exact remaining caller/blocker: compact route-owned private RHF SCF remains
  blocked by nonconvergence; simple damping did not resolve it.

-- repo-manager@macmini
