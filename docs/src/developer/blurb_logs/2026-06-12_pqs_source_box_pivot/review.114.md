Pass 114 manager review

Accepted as a local ignored trace probe.

The trace confirms that the compact route-owned private RHF SCF loop is not
diverging and is not showing an obvious two-cycle. Total energy decreases
monotonically over 50 iterations, but the density residual bottoms near
`4.409e-5` and then drifts upward to `4.544e-5`, so the undamped fixed-point
update is not converging to the current tolerance.

Decision:

- Do not add route wiring.
- Do not add production damping/mixing yet.
- Next pass should run a local ignored controlled damping/mixing sweep on the
  same compact route-smoke fixture. If damping works locally, then design a
  small private SCF-control seam; if it does not, audit the residual/Fock
  convention before changing production code.

Validation/status:

- Local ignored trace probe elapsed: `96.844649625` seconds.
- Reported git status clean and even with origin/main, apart from ignored local
  `tmp/work` artifacts.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: none; this was evidence gathering.
- quarantined: ignored trace artifacts remain local diagnostic evidence.
- not deleted because: trace evidence is useful for the damping decision.
- exact remaining caller/blocker: compact route-owned private RHF SCF remains
  blocked by `:scf_not_converged` without damping/mixing.

-- repo-manager@macmini
