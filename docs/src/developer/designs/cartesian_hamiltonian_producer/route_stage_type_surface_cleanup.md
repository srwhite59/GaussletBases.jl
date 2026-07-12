# Route/Stage Type-Surface Cleanup

Status: implemented by `118a639bf` under `HP-ROUTE-STAGE-TYPE-FN-01` and
`HP-ROUTE-STAGE-TYPE-TEST-01`.

This page is a historical pointer, not current source authority. The pass
vectorized route-sized shellification/lowering summaries, removed duplicate
shellification-owned lowering-contract mirrors, and narrowed stage returns
without changing their ordered facts or behavior. Canonical terminal-lowering
plans remain owned by `CartesianTerminalLowering`.

The current contract is
[Route/stage metadata](route_stage_metadata_contract.md). Registry entries own
lifecycle, maintenance permission, and source/test surfaces. The cleanup did
not change route or shellification semantics, numerical values, public or
driver contracts, artifacts, or materialization policy.
