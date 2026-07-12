# Contract-Plan Vector Cleanup

Status: implemented by `5938ddbc5` under `HP-CONTRACT-VEC-FN-01` and
`HP-CONTRACT-VEC-TEST-01`.

This page is a historical pointer, not current source authority. The pass moved
variable-length terminal-lowering available/selected contracts and retained-unit
transform contracts to vector-backed plan storage while preserving accessor
facts, selection semantics, summaries, and deterministic order.

The current contract is
[Route/stage metadata](route_stage_metadata_contract.md). Registry entries own
lifecycle, maintenance permission, and source/test surfaces. The cleanup did
not change per-contract `source_cpbs`, fixed coordinate/product-box values,
route or numerical behavior, public inputs, driver behavior, or artifacts.
