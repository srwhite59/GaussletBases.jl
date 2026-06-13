Pass 162 review - downstream contract decision

Accepted. The audit recommends the right next step: a CR2 read-only inspection
view, not HFDMRG solver readiness and not HamV6/export.

Reasoning:

- The medium-term target is a CR2-side Be2 WL/PQS comparison.
- The current Be2/PQS handoff is already enough for inspection: dimensions,
  H1 matrix/spectrum, nuclear metadata, electron/spin convention, density
  interaction representation, and explicit blockers.
- It is not enough for HFDMRG density-density readiness because the final-space
  `H,V` solver/Fock/energy convention has not been reviewed.
- It is not enough for HamV6/sliced-integral/export readiness because no
  `V6`, `Vblocks`, dense four-index object, export format, or artifact behavior
  exists.

Decision:

Pass 163 should add only a compact CR2 read-only inspection view from the
existing consumer/handoff. No new payload type. No solver/export claim. The
readiness blocker should remain:

```text
:missing_hfdmrg_density_density_contract
```

and CR2 production/export should remain blocked through a solver/export format
label such as `:cr2_solver_handoff_format`.

Line budget:

Pass 163 must remain net-line-negative across tracked `src` + `test` files.
The implementation can pay for the small CR2 view by removing remaining
duplicated handoff-level nonclaim flags and copied test assertions.

Deletion/shrinkage accounting:

- deleted: none by pass 162
- simplified: next downstream target narrowed to CR2 read-only inspection
- quarantined: HFDMRG solver readiness, HamV6/export, dense `Vee`, final-space
  `V`, `V6`/`Vblocks`, CR2 production format, H1/J, and RHF/SCF remain out of
  scope
- not deleted because: audit only
- exact remaining caller/blocker: active readiness remains blocked on
  `:missing_hfdmrg_density_density_contract`

-- repo-manager@macmini
