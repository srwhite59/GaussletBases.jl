# Cartesian Execution Conformance Audit - 2026-07-12

Status: read-only static audit of exact committed baseline `c8c1a4911`.
This is review evidence, not authority and not numerical certification.

## Scope And Method

The audit covered all `150` execution records derived from schema-v3
`authority.toml`, grouped across all `40` canonical contract combinations.
Eight disjoint reviews compared each record with its canonical contract,
owned committed paths, function signatures and defaults, tests, and live
callers. The unrelated uncommitted `hamiltonian_corrections.jl` work was
excluded by reading committed content.

Classifications mean:

- `MATCHED`: static implementation, authority, and contract agree;
- `DOCUMENTED_GAP`: partial or pending behavior is described accurately;
- `NUMERICAL_VALIDATION_REQUIRED`: static construction agrees, but the
  claimed endpoint or physical gate is not currently reproduced by a
  committed positive test;
- `DISCREPANCY`: source, test, authority metadata, or canonical wording does
  not agree and requires a bounded reconciliation.

No Julia endpoint or numerical test was run. `MATCHED` therefore does not
certify scientific accuracy beyond the existing committed evidence.

## Result

| Classification | Records |
| --- | ---: |
| `MATCHED` | 107 |
| `DOCUMENTED_GAP` | 11 |
| `NUMERICAL_VALIDATION_REQUIRED` | 8 |
| `DISCREPANCY` | 24 |
| **Total** | **150** |

Machine structure remained internally synchronized: `232` total records,
`44` hashed documents, and `150` execution IDs. All owned paths and generated
views passed the permanent checker. The audit separately found one dependency
cycle, described below.

This Pass 398 immediately removes that cycle and extends the checker to reject
future cycles or execution dependencies on closed/no-grant records. The other
findings remain remediation candidates after this review record lands.

## Highest-Priority Findings

1. `HP-PQS-SHELLQ-OVERRIDE-TEST-01` is marked completed, but the committed
   test covers H2 base geometry and rejection behavior rather than the required
   padded-Be2 packet capture, `J0/E0`, correction, and due-diligence gate.
2. White-Lindsey terminal realization is numerically compact, but live
   lowering/inventory records still encode retired identity/trivial-embedding
   semantics in `cartesian_terminal_lowering/`, `cartesian_retained_units/`,
   and `cartesian_retained_unit_transform_contracts/`.
3. The authority graph contained
   `HP-R3U-FN-01 -> HP-R3U-ZDI-FN-01 -> HP-R3U-FN-01`. The correct contract
   layering is only `HP-R3U-ZDI-FN-01 -> HP-R3U-FN-01`; the generic facade is
   the prerequisite of its z-axis scope extension.
4. `HP-RG-PROTECT-ONEBODY-FN-01` transforms and symmetrizes the fixed-sector
   matrices but does not produce the finite, dimension, symmetry, and geometry
   diagnostics required by its authority and canonical contract.
5. Protected ladder readback omits some written manifest facts, and restart
   trace loss is computed from occupied-column count instead of the reported
   source trace.
6. Direct packet helpers do not all reject unconverged packets, and the packet
   writer does not validate every required ordinary-fit consistency field.
7. The screened-Hartree correction test labels a determinant-density exact
   fixture as `:density_fit`; its algebra check is valid, but the source label
   is not evidence for the fitted-density path.

## Other Reconciliation Findings

- `HP-FN-03`: the one-body kernel does not reject nonfinite Gaussian
  coefficients at its owned boundary.
- `HP-R1-ATOM-WIRE-01`: a caller-free atom-only materialization helper cluster
  remains despite the shared-stage-only contract.
- `HP-DRV-SHELLDD-FN-01`: in-memory `warning_flags` and `warning_summary`
  shapes are reversed relative to the canonical contract.
- `HP-ROUTE-STAGE-TYPE-FN-01`: a tracked inspection generator still calls
  retired route-stage signatures, although active package callers are current.
- `HP-COMP-SUPPATOM-TEST-01`, `HP-COMP-SUPPWL-TEST-01`, and
  `HP-COMP-WLNS-TEST-01`: named positive behaviors are absent from their owned
  committed test files.
- `HP-R3U-ZDI-TEST-01`: the canonical page overstates committed Be2 coverage;
  authority and source currently own H2 assertions only.
- `HP-COMP-THINSLAB-FN-01` and `HP-COMP-THINSLAB-META-FN-01`: stale identity
  slab kinds and counters remain in compatibility/inventory code even though
  the active numerical realizer uses compact face products.
- `HP-PQS-MAP-SFACTOR-TEST-01`: no valid nonunit factor or explicit-`1.0`
  parity case exists in the owned test.
- `HP-MCOMX-OBJ-01`: a private helper accepts arbitrary positive `lambda`
  despite the fixed `0.5` contract; the live route still uses `0.5`.
- `HP-RG-INJECT-FN-01`: enabled direct injection lacks the required final
  identity check and material-negative-Gram rejection. The path is default-off
  and has no production caller.
- `HP-RG-NUMCOMP-FN-01`: one paragraph still calls semantic shell-q override
  approved-pending although it is implemented.
- `HP-RG-OCC-FIRST-INJECT-FN-01`: the contract overstates returned individual
  base projections; source returns only the combined projection.
- `HP-RG-PROTECT-INJECT-FN-01`: authority says geometry returns `Q_perp/F`,
  while the operator owner currently reconstructs them from geometry output.
- `HP-PQS-ATOMREF-PACKET-FN-01`: direct packet consumption validation is
  incomplete even though current screened/additive callers validate first.

## Scientific Policy Finding

Current source, tests, canonical contracts, and Passes 327, 328, and 360 agree
on the implemented packet policy:

```text
occupied determinant -> P0 and q0
density fit          -> compressed reference cloud and E0
potential fit        -> approximate fast J0
Tr(P0*J0_fit)-E0_fit -> reported consistency error
```

The neutral mixed-Hartree kernels remain available as determinant-density or
density-fit oracle infrastructure. Changing the production interpretation so
that determinant density is authoritative for exact `J0` and `E0` would be a
new scientific contract amendment, not a correction of source/document drift.
The controlled Cr2 screened comparison should not be interpreted under a
different convention until that choice is explicit.

## Complete Disposition Matrix

### Matched - 107

- Terminal/R1: `HP-FILE-01`, `HP-FN-00`, `HP-FN-01`, `HP-FN-02`,
  `HP-FN-04`, `HP-FN-05`, `HP-OBJ-01`, `HP-OBJ-02`, `HP-WIRE-01`,
  `HP-R1-ART-01`, `HP-R1-CORE-FN-01`, `HP-R1-FILE-01`, `HP-R1-FN-01`,
  `HP-R1-TEST-01`, `HP-R1-WIRE-01`, `HP-COMP-ATOMBOX-FN-01`,
  `HP-R1-ATOM-FN-01`, `HP-R1-ATOM-TEST-01`, `HP-COMP-NSCORE-FN-01`.
- Driver/artifact/route: `HP-DRV-FILE-01`, `HP-DRV-FN-01`,
  `HP-DRV-INV-FN-01`, `HP-DRV-NEST-FN-01`, `HP-DRV-NEST-WIRE-01`,
  `HP-DRV-STAGE-FN-01`, `HP-DRV-STAGE-WIRE-01`, `HP-DRV-ATOM-CLEAN-01`,
  `HP-DRV-ATOM-FN-01`, `HP-DRV-ATOM-WIRE-01`, `HP-HAM-MANIFEST-FN-01`,
  `HP-NEST-ART-FN-01`, `HP-R3-ART-01`, `HP-CONTRACT-VEC-FN-01`,
  `HP-RAW-SRCMODE-FN-01`, `HP-ROUTE-INV-FN-01`,
  `HP-ROUTE-STAGE-CARRIER-FN-01`.
- Composition/WL/R3U: `HP-COMP-BASEDIAT-FN-01`,
  `HP-COMP-BASEDIAT-TEST-01`, `HP-COMP-NS-FN-01`, `HP-COMP-NS-TEST-01`,
  `HP-COMP-WLNS-FN-01`, `HP-ROUTE-RECIPE-FN-01`, `HP-WLTERM-FILE-01`,
  `HP-WLTERM-FN-01`, `HP-WLTERM-TEST-01`, `HP-WLTERM-WIRE-01`,
  `HP-R3U-FILE-01`, `HP-R3U-TEST-01`, `HP-R3U-WIRE-01`,
  `HP-R3U-ZDI-FN-01`, `HP-R3U-ZDI-WIRE-01`.
- Shell/source/mapping: `HP-COMP-ANGBOX-FN-01`,
  `HP-COMP-FACEPROD-FN-01`, `HP-COMP-SHELLGEOM-DIAT-FN-01`,
  `HP-COMP-SHELLGEOM-FN-01`, `HP-PQS-ASPECTSHELL-FN-01`,
  `HP-MCOMX-DRV-FN-01`, `HP-MCOMX-FILE-01`, `HP-MCOMX-TERM-FN-01`,
  `HP-MCOMX-WIRE-01`, `HP-PQS-MAP-SFACTOR-FN-01`.
- Residual/R3: `HP-R3-FN-01`, `HP-R3-FN-02`, `HP-RG-FILE-01`,
  `HP-RG-OBJ-01`, `HP-RG-FN-01`, `HP-RG-FN-02`, `HP-RG-FN-03`,
  `HP-RG-FN-04`, `HP-RG-WIRE-01`, `HP-RG-TEST-01`, `HP-R3-OBJ-01`,
  `HP-R3-FN-03`, `HP-R3-TEST-01`, `HP-RG-CUTOFF-FN-02`,
  `HP-RG-CUTOFF-TEST-02`, `HP-RG-ORTHO-FN-01`, `HP-RG-ORTHO-TEST-01`.
- Protected/transfer: `HP-RG-OCC-FIRST-INJECT-TEST-01`,
  `HP-RG-PROTECT-ART-FN-01`, `HP-RG-PROTECT-ARTLOC-FN-01`,
  `HP-REP-XGTO-IMPORT-FN-01`, `HP-REP-XGTO-IMPORT-TEST-01`,
  `HP-REP-XGTO-PROTECT-SIDECAR-FN-01`,
  `HP-REP-XGTO-PROTECT-SIDECAR-TEST-01`.
- Reference: `HP-PQS-SCREEN-HARTREE-CORR-FN-01`,
  `HP-RHO0-MIXH-FEXACT-FN-01`, `HP-RHO0-MIXH-GAAA-FN-01`,
  `HP-RHO0-MIXH-GAAA-TEST-01`, `HP-RHO0-MIXH-GG-FN-01`,
  `HP-RHO0-MIXH-GG-TEST-01`.
- Raw blocks/Coulomb reuse: `HP-CGRB-FILE-01`, `HP-CGRB-FN-01`,
  `HP-CGRB-FN-02`, `HP-CGRB-TEST-01`, `HP-CGRB-WIRE-01`,
  `HP-CGRB-NN-FILE-01`, `HP-CGRB-NN-FN-01`, `HP-CGRB-NN-TEST-01`,
  `HP-CGRB-NN-WIRE-01`, `HP-R3GG-FN-01`, `HP-R3GG-TEST-01`,
  `HP-R3UN-FN-01`, `HP-R3UN-TEST-01`, `HP-R3BASE-DRV-WIRE-01`,
  `HP-R3BASE-FN-01`, `HP-R3BASE-TEST-01`.

### Documented Gaps - 11

`HP-COMP-NSCORE-TEST-01`, `HP-HAM-MANIFEST-SRC-FN-01`,
`HP-HAM-MANIFEST-SRC-TEST-01`, `HP-HAM-MANIFEST-TEST-01`,
`HP-NEST-ART-TEST-01`, `HP-RG-PROTECT-EGOI-FN-01`,
`HP-RG-PROTECT-EGOI-TEST-01`, `HP-PQS-ATOMREF-PACKET-TEST-01`,
`HP-RG-PROTECT-ADDREF-TEST-01`, `HP-PQS-COULOMB-ACCURACY-FN-01`,
`HP-PQS-COULOMB-ACCURACY-TEST-01`.

### Numerical Validation Required - 8

`HP-COMP-ATOMBOX-TEST-01`, `HP-COMP-SUPPATOM-FN-01`,
`HP-COMP-SUPPWL-FN-01`, `HP-WLDIAT-PARITY-FN-01`,
`HP-PQS-SHELLQ-OVERRIDE-FN-01`, `HP-MCOMX-FN-01`,
`HP-RG-NUMCOMP-TEST-01`, `HP-RG-PROTECT-ADDREF-FN-01`.

### Discrepancies - 24

`HP-FN-03`, `HP-R1-ATOM-WIRE-01`, `HP-DRV-SHELLDD-FN-01`,
`HP-ROUTE-STAGE-TYPE-FN-01`, `HP-COMP-SUPPATOM-TEST-01`,
`HP-COMP-SUPPWL-TEST-01`, `HP-COMP-WLNS-TEST-01`,
`HP-COMP-WLDIAT-FN-01`, `HP-WLDIAT-COMPACT-FN-01`, `HP-R3U-FN-01`,
`HP-R3U-ZDI-TEST-01`, `HP-COMP-THINSLAB-FN-01`,
`HP-COMP-THINSLAB-META-FN-01`, `HP-PQS-SHELLQ-OVERRIDE-TEST-01`,
`HP-MCOMX-OBJ-01`, `HP-PQS-MAP-SFACTOR-TEST-01`, `HP-RG-INJECT-FN-01`,
`HP-RG-NUMCOMP-FN-01`, `HP-RG-OCC-FIRST-INJECT-FN-01`,
`HP-RG-PROTECT-INJECT-FN-01`, `HP-RG-PROTECT-ONEBODY-FN-01`,
`HP-RG-PROTECT-LADDER-BUNDLE-FN-01`, `HP-PQS-ATOMREF-PACKET-FN-01`,
`HP-PQS-SCREEN-HARTREE-CORR-TEST-01`.

## Recommended Remediation Order

1. **Authority integrity:** remove the R3U reverse dependency and enforce
   acyclic, active-prerequisite dependency semantics.
2. **Fail-fast correctness:** packet convergence/write validation, nonfinite
   one-body coefficients, direct-injection Gram/identity checks, and protected
   one-body diagnostics.
3. **Consumer data correctness:** ladder readback/trace-loss and due-diligence
   warning-shape reconciliation.
4. **Stale path removal:** WL/thin-slab identity inventory, dead atom-only
   materialization helpers, obsolete tracked generator calls, and the latent
   mapped-COMX lambda knob.
5. **Validation completion:** shell-q padded Be2, `s_factor` nonunit parity,
   numerical-complete/additive positive consumers, mapped-COMX, and
   supplemented atom/WL gates. Add only tests that protect live behavior.
6. **Wording-only corrections:** occupied-first returned fields,
   numerical-complete shell-q status, and ZDI H2/Be2 test scope.

Each remediation should use existing IDs and canonical owners where the
authority already covers the fix. A separate docs-only amendment is required
if a fix needs a new source path, persistent shape, public surface, or changed
scientific convention.
