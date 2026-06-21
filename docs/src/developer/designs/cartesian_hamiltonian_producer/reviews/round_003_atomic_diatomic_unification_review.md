# Round 003 Atomic/Diatomic Unification Review

Reviewer: ChatGPT-Pro milestone recommendation, pasted by user

Date: 2026-06-20

## Verdict

Settle atomic/diatomic unification in the design before Slice A implementation.
Do not expand Slice A to the full atomic Hamiltonian, but require one-center
terminal-basis realization through the same generic entry point as H2 and Cr2.

## Main Findings

- One-center PQS already reaches the same terminal shellification and lowering
  spine as bond-aligned diatomics.
- Atomic and diatomic geometry producers may differ, but everything after
  terminal support, retained, and transform records should be shared.
- If Slice A is implemented only under the diatomic target payload, it can
  silently acquire bond-axis, two-center, role-name, or contact-core assumptions.
- The current `CartesianTerminalBasisBlock` and
  `CartesianTerminalBasisRealization` candidates are suitable because they carry
  no atom-count, route-kind, bond-axis, or role-specific fields.

## Required Design Changes

1. Add an atomic/diatomic unification invariant.
2. Require Slice A validation on one-center atomic, contact-core H2, and
   separated Cr2 terminal records.
3. Add a wiring design ID that makes the generic terminal-basis stage consume
   typed terminal records and forbids dispatch on system classification, atom
   count, route kind, bond axis, or terminal role names.
4. Preserve the old atomic H1 path only as a migration oracle until the generic
   one-body slice reproduces the reviewed atomic endpoint.
5. Freeze later Slice B/C/D contracts so they consume
   `CartesianTerminalBasisRealization` without atomic/diatomic branches.

## Deferred Scope

The atomic Hamiltonian, hydrogen energy endpoint, one-body operators, IDA, and
artifact production remain future slices. Slice A only realizes and validates
the terminal basis.
