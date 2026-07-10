# Complete-Core-Shell RHF Retirement

Status: completed by `28e9b2c84`. This is a historical deletion record;
`HP-RETIRE-CCS-RHF-FN-01` and `HP-RETIRE-CCS-RHF-TEST-01` are closed and no
longer authorize source or test work.

## Decision

The old complete-core-shell RHF payload stack was retired. The canonical
Cartesian producer is the staged driver and `CartesianIDAHamiltonian` artifact
path; atom/diatomic, base/supplemented, and PQS/WL consumers do not use the old
RHF payload workflow.

## Removed Surface

Commit `28e9b2c84` removed:

```text
src/pqs_multilayer_complete_core_shell_rhf.jl
```

and its include from `src/GaussletBases.jl`, deleting 1,879 source lines. The
removed file carried route-era input-contract, SCF-payload, one-step-payload,
and blocked-status vocabulary. The accepted pre-deletion audit found no live
`src`, `bin`, `test`, or `tool` caller outside the file and root include.

## Current Guardrail

Current focused scans find no live reference to the removed stack. Do not
restore it, add an adapter or replacement payload, or infer authority from the
historical source/test plan retained in git history. A new docs-only amendment
is required before any related source work.
