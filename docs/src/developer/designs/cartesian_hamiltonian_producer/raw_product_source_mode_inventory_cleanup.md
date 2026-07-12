# Raw Product Source-Mode Inventory Cleanup

Status: implemented by `34cf8f106` under `HP-RAW-SRCMODE-FN-01` and
`HP-RAW-SRCMODE-TEST-01`.

This page is a historical pointer, not current source authority. The pass moved
the variable-length `RawProductBoxPlan` mode and column inventories to vectors
while preserving fixed `NTuple{3,Int}` mode coordinates/dimensions,
deterministic source-mode order, column association, retained-rule behavior,
and manifest source provenance.

The current contract is
[Route/stage metadata](route_stage_metadata_contract.md). Registry entries own
lifecycle, maintenance permission, and source/test surfaces. The cleanup did
not change pair/source-box algorithms, terminal realization, numerical values,
public inputs, driver behavior, or artifact schemas.
