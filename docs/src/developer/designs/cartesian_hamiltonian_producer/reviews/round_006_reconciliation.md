# Round 006 Reconciliation

Status: final recursive spike reconciled; Slice A ready for freeze decision

Repo-manager reconciled the final uncommitted recursive-projection spike into
`docs/src/developer/cartesian_hamiltonian_producer_design.md`. The design
remains unapproved and is not implementation authority until Slice A IDs are
explicitly frozen and bound.

## Accepted Changes

- Set the default projection threshold to `projection_atol = 1.0e-12`.
- Explained that `projection_atol` is a roundoff-subtraction threshold, distinct
  from final cross-overlap acceptance.
- Added a Slice A support-pair workspace cap of `64 MiB` unless amended.
- Required tiling or streaming for local overlap actions that exceed the cap.

## Evidence From Spike

- One-center, H2, and Cr2 used the same terminal realization inputs once
  terminal records were reached.
- True recursive accepted blocks were used.
- All projection residuals were below `1e-12` and skipped.
- Effective supports remained unchanged.
- Shell ranks remained `98`.
- Final cross overlaps were small:
  - one-center: `8.073e-16`
  - H2: `3.095e-15`
  - Cr2: `5.463e-14`
- Later Cr2 direct records remained orthogonal to prior blocks at about
  `2.4e-15 .. 5.3e-15`.

## Remaining Implementation Caveats

- One-center public staging still reaches terminal records through an old
  route-shape skeleton. Slice A must connect one-center terminal records to the
  typed terminal contract without freezing that skeleton as the public contract.
- Cr2's largest local dense workspace in the spike was `175.928 MiB`; production
  Slice A must tile or stream that local action to satisfy the `64 MiB` cap.

## Freeze Implication

No further broad design review is recommended. The next docs-only step is to
freeze and bind the Slice A IDs if the user/manager accepts this evidence:

```text
HP-OBJ-01
HP-OBJ-02
HP-FILE-01
HP-FN-00
HP-FN-01
HP-FN-02
HP-WIRE-01
```

`HP-FN-03`, `HP-FN-04`, and `HP-FN-05` remain future candidates.
