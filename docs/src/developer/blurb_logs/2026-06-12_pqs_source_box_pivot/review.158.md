Pass 158 review - Be2/PQS Hamiltonian consumer audit

Accepted. The audit stayed read-only and clarified the next boundary.

Current Be2/PQS handoff already provides enough for a private inspect-only
consumer contract:

- final-space H1 matrix with dimension 221;
- pre-final density interaction representation;
- pre-final pair matrix and final-to-pre-final coefficients;
- support weights, raw pair numerator, and raw pair-factor terms;
- density gauge `:pre_final_localized_positive_weight`;
- raw pair convention `:raw_numerator`;
- nuclear charges/coordinates/repulsion;
- electron count and spin sector.

The audit also confirmed what this is not:

- not HFDMRG density-density ready, because HFDMRG expects a reviewed final-space
  `H, V` contract;
- not HFDMRG sliced-integral ready, because no `V6` or `Vblocks` representation
  is materialized;
- not CR2 production ready, because there is no bundle/JLD2/HamV6 handoff
  format or seed-orbital contract;
- not public/export/artifact behavior.

Next decision:

Add a private route-owned Hamiltonian consumer contract payload, but only as a
replacement for the current blocker:

```text
:missing_diatomic_hamiltonian_consumer_contract
```

The pass is valid only if the probe-enabled Be2 readiness no longer claims the
consumer contract is missing. It should instead report that a private inspector
contract is available while downstream-specific contracts remain blocked, such
as HFDMRG density-density `H,V`, sliced integrals/HamV6, and CR2 handoff format.

This is an adoption-boundary object, not another payload layer to preserve
indefinitely. Its expiration condition is explicit: once a real downstream
format/consumer is chosen, this inspect-only contract must either be replaced
by that format-specific contract or collapsed into its summary.

Deletion/shrinkage accounting:

- deleted: none by pass 158
- simplified: downstream boundary is now specific enough for one replacement
  implementation pass
- quarantined: current handoff remains private inspect-only
- not deleted because: audit only
- exact remaining caller/blocker: readiness still blocks on
  `:missing_diatomic_hamiltonian_consumer_contract` until the next pass replaces
  that blocker with downstream-specific missing contracts

-- repo-manager@macmini
