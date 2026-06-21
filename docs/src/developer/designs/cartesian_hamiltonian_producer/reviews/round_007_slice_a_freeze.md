# Round 007 Slice A Freeze

Status: Slice A frozen and bound

Date: 2026-06-20

Repo-manager froze Slice A after the final recursive-projection spike and user
approval to proceed.

## Approved Slice A IDs

```text
HP-OBJ-01
HP-OBJ-02
HP-FILE-01
HP-FN-00
HP-FN-01
HP-FN-02
HP-WIRE-01
```

These IDs are now implementation authority through
`docs/src/developer/cartesian_hamiltonian_producer_design.md` and the binding
rule in `AGENTS.md`.

## Still Candidate-Only

```text
HP-FN-03
HP-FN-04
HP-FN-05
```

These do not authorize one-body assembly, IDA assembly, Hamiltonian artifact
production, or driver simplification.

## Freeze Boundary

Slice A authorizes terminal-basis realization only:

- one-center atomic, contact-core H2, and separated Cr2 terminal records use
  the same terminal-basis entry point;
- projection uses recursive accepted blocks with `projection_atol = 1.0e-12`;
- local overlap/projection work respects the `64 MiB` support-pair workspace
  cap through tiling or streaming;
- the current terminal source-realization preflight path is deleted by the
  implementation.

Unlisted production surfaces require a prior docs-only amendment.
