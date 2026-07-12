# Authority Transition Rehearsal - 2026-07-12

Status: completed read-only review; authority cutover is **not ready**.

This is historical review evidence. It does not grant work and does not replace
the prose [registry](../registry.md), marked `AGENTS.md` whitelist, current
status, invariants, or canonical subsystem contracts.

## Baseline And Inputs

- repository HEAD: `32119419315d2e6a873856fe952e2e3d90b2afea`;
- candidate SHA-256:
  `0028eaddbad0317f193508d59ec77d40e8443cd1d8b08cf46a0bfc11e8d1d1b0`;
- candidate checker SHA-256:
  `5125246d1b0312207ba5712ba0e88658975e0bb6bfcf28bd90aba9a58e82f4fd`;
- candidate records: `232`;
- document inventory: `44`;
- candidate-derived execution IDs: `150`;
- generated registry preview: `3,412` lines;
- authoritative prose registry: `4,844` lines.

Two independent external renders were byte-identical. All committed shadow,
candidate, and transition checks passed before review. Those checks proved
determinism and current structural parity; they did not prove complete semantic
parity.

## Review Coverage

Three independent reviewers covered disjoint candidate ranges:

- records 1-80: `80/80` reviewed;
- records 81-160: `80/80` reviewed;
- records 161-232: `72/72` reviewed.

A fourth reviewer covered transition tooling, generated views, CI integration,
and the full `232`-record/`150`-execution-ID inventory without duplicating the
record-level semantic review.

The execution set itself matched exactly. The no-go result comes from
record-level ownership/dependency drift and transition-tooling limitations.

## Confirmed Candidate Reconciliation

The next bounded metadata pass must reconcile these records against the current
prose registry and canonical contracts:

| ID | Required reconciliation |
| --- | --- |
| `HP-DRV-STAGE-FN-01` | Add the three existing operator-factoring source owners named by the registry. |
| `HP-MCOMX-FILE-01` | Add the existing PQS source-axis owner. |
| `HP-OBJ-03` | Classify the live rejection boundary as canonical, not history. |
| `HP-PQS-ATOMREF-PACKET-FN-01` | Bind the neutral mixed-Hartree owner as a dependency; do not silently make that helper packet-owned. |
| `HP-PQS-SCREEN-HARTREE-CORR-FN-01` | Add module/include wiring and the bounded packet-helper source surface named by the canonical contract. |
| `HP-R3BASE-FN-01` | Remove the protected-ladder consumer from owned paths; its authority is separate. |
| `HP-R3U-FILE-01` | Restore the test surface and exact nested validation path. |
| `HP-REP-XGTO-IMPORT-FN-01` | Remove overlap, transfer, and ladder dependencies from owned edit paths. |
| `HP-REP-XGTO-PROTECT-SIDECAR-FN-01` | Remove overlap and ladder dependencies from owned edit paths. |
| `HP-RES-01` | Classify the live rejection boundary as canonical, not history. |
| `HP-RG-CUTOFF-TEST-01` | Restore the historical `HP-RG-CUTOFF-FN-01` dependency. |
| `HP-RG-CUTOFF-TEST-02` | Restore the active `HP-RG-CUTOFF-FN-02` dependency. |
| `HP-RG-FN-01` | Replace the numerical-complete misc-test evidence with the ordinary H2 endpoint evidence. |
| `HP-RG-IDTOL-FN-01` | Restore the ORTHO dependency and preserve the superseded-by relationship explicitly. |
| `HP-RG-IDTOL-TEST-01` | Restore the historical IDTOL function dependency. |
| `HP-RG-OBJ-01` | Restore the RG file dependency. |
| `HP-RG-ORTHO-TEST-01` | Restore the ORTHO function dependency. |
| `HP-RG-WIRE-01` | Restore dependencies on `HP-RG-FN-01` through `HP-RG-FN-04`. |
| `HP-WLTERM-FILE-01` | Add the existing module-wiring source path. |

The XGTO, protected-ladder, and packet findings expose an important modeling
rule: an implementation dependency is not automatically an owned edit path.
The candidate must keep ownership fail-closed and leave non-owned dependencies
in canonical contracts or an explicitly non-authorizing dependency field.

The IDTOL finding also requires an explicit decision on structural
`superseded_by`/`supersedes` metadata. It must not be hidden inside an ordinary
dependency edge if the relationship has different semantics.

## Confirmed Tooling Reconciliation

Before another authority rehearsal:

1. Parse the marked whitelist in full-document Markdown context so a fence or
   HTML comment cannot hide an apparently valid block.
2. Bind rehearsal manifests to the transition snapshot, prose registry, exact
   marked block, candidate/checker bytes, and Git commit.
3. Apply a conservative grammar or correct escaping to generated Markdown link
   destinations.
4. Put a conspicuous non-authoritative warning in the standalone whitelist
   preview while retaining a separately rendered exact block for parity.
5. Describe transition parity accurately: current automated comparison covers
   identity, owner-document paths, candidate bytes, and execution membership;
   semantic field parity still requires reviewed evidence.
6. Pin the reviewed candidate digest. Regenerating a transition snapshot is not
   itself semantic approval.

The current external preview writer is not an atomic authority-cutover tool and
must not be promoted as one. A later approved cutover needs generation from an
immutable checkout, complete output verification, one canonical machine-to-view
direction, and a required canonical-authority CI gate.

## Decision

Authority remains with the prose registry and marked `AGENTS.md` whitelist.
The schema-v3 candidate remains non-authoritative and
authorization-incomplete.

Next sequence:

1. perform one bounded candidate/tooling reconciliation pass;
2. regenerate transition metadata and external previews;
3. repeat an independent semantic/tooling rehearsal;
4. only then request a separate atomic cutover decision.
