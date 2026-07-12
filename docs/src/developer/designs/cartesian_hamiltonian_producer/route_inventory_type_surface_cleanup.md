# Route Inventory Type-Surface Cleanup

Status: implemented by `c985723c7` under `HP-ROUTE-INV-FN-01` and
`HP-ROUTE-INV-TEST-01`.

This page is a historical pointer, not current source authority. The pass
replaced runtime-keyed retained-unit and pair-family `NamedTuple` inventories
with ordered vector-backed rows and label lookup, leaving labels as data rather
than concrete type parameters.

The current contract is
[Route/stage metadata](route_stage_metadata_contract.md). Registry entries own
lifecycle, maintenance permission, and source/test surfaces. The cleanup did
not change recipe, shellification, lowering, numerical, public, driver, or
artifact behavior.
