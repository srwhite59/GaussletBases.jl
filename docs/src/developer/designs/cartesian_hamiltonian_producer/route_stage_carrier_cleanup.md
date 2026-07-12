# Route/Stage Carrier Cleanup

Status: implemented by `8c3df2ad9` under
`HP-ROUTE-STAGE-CARRIER-FN-01` and `HP-ROUTE-STAGE-CARRIER-TEST-01`.

This page is a historical pointer, not current source authority. The pass
stopped carrying the full route skeleton and duplicated shellification,
lowering, support, retained-rule, and transform-plan payloads through later
stage signatures. Variable inventories became vectors or compact ordered rows;
later stages retain only their live plans, realizations, and summaries.

The current contract is
[Route/stage metadata](route_stage_metadata_contract.md). Registry entries own
lifecycle, maintenance permission, and source/test surfaces. Stopping the
wholesale carry did not retire route-skeleton construction, active
complete-core-shell support, pair/assembly/report stages, or tools, and did not
change route, shellification, numerical, public, driver, or artifact behavior.
