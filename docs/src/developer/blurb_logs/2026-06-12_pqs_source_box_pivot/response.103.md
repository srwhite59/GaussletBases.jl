Pass 103 response

No-edit fixture-rule design audit completed.

Files changed:

- none before this response/log

Read-only inspection used:

- `test/nested/pqs_direct_retained_final_h1_runtests.jl`
- `src/pqs_source_box_route_driver_helpers.jl`
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/successor_handoff_2026-06-12_pqs_source_box_pivot.md`
- `/Users/srw/Dropbox/chatarchive/handoff/software_packets/gausslet-methods-fundamentals.md`
- recent pass logs/reviews around 095-102
- existing q-ladder / side13 / H1/J probe notes surfaced by read-only search
- `git status --short --branch`

Fixture-rule proposal

- Treat compact fixtures as route-smoke/convention fixtures unless a separate
  fixture card promotes them to physics endpoints.
- Before RHF implementation, define a small fixture policy card with:
  - fixture role: route smoke, oracle/debug, convention diagnostic, or physics
    endpoint;
  - physical target: `Z`, electron count, charge/spin/closed-shell rule;
  - geometry/resolution family: spacing/distortion/radius/q/source-box/shell
    depth;
  - expected numerical claims and nonclaims;
  - timing category: routine, slow/manual, or ignored probe only.
- Do not let `final_dimension == 223`, H1/J materialization, or the 0.21 s
  warm route runtime turn the compact fixture into a scientific gate by inertia.

Current compact fixtures

- Tracked focused H1 fixture in
  `test/nested/pqs_direct_retained_final_h1_runtests.jl`:
  - `Z = 1.0` via `center.charge = 1.0`
  - `q = 5`
  - `n_s = 5` / core side `5`
  - parent axis counts: `7 x 7 x 7`
  - current/parent box: `(1:7, 1:7, 1:7)`
  - core/source box: `(2:6, 2:6, 2:6)` = `5 x 5 x 5`
  - support count: `343`
  - core support count: `125`
  - shell support count: `218`
  - shell retained count: `98`
  - final dimension: `223`
  - H1 energy: `-0.48047934800387226`
  - fixed-block H1 oracle agreement: `< 1.0e-6`
  - self-Coulomb/density gauge: not asserted in this tracked test.

- Same compact direct structured H1/J convention probe from pass 097:
  - final dimension: `223`
  - H1 energy: `-0.48047934800387126`
  - self-Coulomb: `0.6397851751855723`
  - density gauge: `pre_final_localized_positive_weight`
  - role: density-gauge/H1/J convention diagnostic over the compact H1 seam.

- Accepted one-center source-box driver H1/J dry-run from pass 101:
  - one-center route metadata uses `atom_symbols = ("Be",)` and
    `nuclear_charges = (4,)`, so `Z = 4`
  - `q = 5`, `n_s = 5`
  - parent axis counts: `7 x 7 x 7`
  - compact source/core geometry remains the `5 x 5 x 5` source-box family
  - final dimension: `223`
  - H1 energy: `-5.6629907690725245`
  - self-Coulomb: `1.8691288063594704`
  - density gauge: `pre_final_localized_positive_weight`
  - role: private route-owned H1/J smoke; not a neutral-Be physics endpoint.

Route-smoke facts

- Source-box-first route reaches shellification/lowering-backed region plan.
- Complete core/shell source plan is available.
- Complete core/shell final basis materializes with final dimension `223`.
- Complete H1 payload materializes.
- Route-owned density inputs feed H1/J without signed final-weight division,
  raw no-division density, or density-normalized pair-term authority.
- Driver H1/J status materializes with RHF/GTO/export/artifact flags false.
- Warm repeated route runtime is short (`0.210283` s), while cold first-call
  latency is slow (`120.102067` s).

Scientific endpoint / physics facts

- The tracked Z=1 H1 value and fixed-block agreement are physically meaningful
  one-electron seam checks for this compact geometry, but not convergence or
  endpoint acceptance.
- The q-ladder/side13 He probe history is closer to physics endpoint evidence:
  q=11 reached final dimension `1933`, Z=2 H1 `-1.9949555260353655`, RHF total
  `-2.8559475204289022`, and about `0.00573 Ha` error vs the He HF reference.
- Those side13/q-ladder data are ignored probe evidence, not permanent gates
  yet. A manager fixture decision is still needed before promoting q=9, q=11,
  or any side13 fixture.

Oracle/debug facts

- Fixed-block H1 comparison in the tracked test is oracle/debug.
- Explicit-box final-basis equivalence is compatibility/debug.
- Shell/support-row contraction remains oracle/debug, not the PQS algorithm.
- The H1/J self-Coulomb scalar is a convention diagnostic for the pre-final
  density gauge, not a final IDA/quadrature result.
- Retained diagnostic weights remain diagnostic only.

Parameter families that must move together

- `Z`, electron count, charge, and closed-shell/spin convention.
- Spacing parameters: `d`/core spacing, reference spacing, tail spacing, and
  distortion/mapping parameter `s`.
- `q`, `n_s`, source/core side length, shell depth, and retained-rule family.
- Parent radius/axis count and source-box placement.
- Coulomb expansion/factor inputs and density gauge.
- Compact `7^3 / 5^3 / dim 223` fixtures and side13/q-ladder fixtures should
  not be compared as if only `q` changed; they are different resolution and
  physics-evidence regimes unless a fixture family explicitly ties them.

What should not become an acceptance gate by inertia

- The compact driver H1/J dry-run with `Z=4`, final dimension `223`, and H1
  energy `-5.6629907690725245`.
- The tracked compact Z=1 H1 fixed-oracle seam test.
- The self-Coulomb scalar alone.
- The accepted pre-final positive density gauge alone.
- The side13/q-ladder ignored probes as a group; choose one reviewed fixture if
  a physics gate is needed.
- Any shell/support-row or fixed-block oracle path as route authority.

Fixture policy needed before RHF

- Decide whether the first RHF work is:
  - route smoke only over the compact driver fixture;
  - a closed-shell He physics endpoint probe;
  - or both, separated into different fixtures.
- If physics endpoint: choose the target (`He`, `Z=2`, electron count `2`,
  closed-shell RHF), fixture size (`q=9`, `q=11`, or side13 family), and
  accepted reference/error/timing thresholds before implementation.
- If route smoke: keep it compact and assert only object/status/non-promotion
  facts, not HF accuracy.
- Record explicitly that compact H1/J route materialization does not imply RHF
  readiness or production route adoption.

Smallest next pass recommendation

- Next design pass, if manager wants one more before code: write a short
  tracked fixture policy note or update the PQS framework with the fixture role
  table above.
- Smallest implementation pass after policy approval: add the private RHF input
  contract payload from pass 102, blocked/available only, with explicit
  electron-count and fixture-role fields. Do not run SCF in that first pass.

Git status:

```text
## main...origin/main
```

Deletion/shrinkage forecast:

- deleted: none in this no-edit audit.
- simplified: future tracked H1 cleanup can shrink explicit-box and duplicate
  nuclear convention assertions after a reviewed route smoke or fixture policy
  replaces that pressure.
- quarantined: side13/q-ladder RHF probes and compact H1/J timing probes remain
  ignored/historical evidence until a fixture is explicitly promoted.
- not deleted because: current tracked H1 fixture is still the only compact
  tracked seam/oracle gate; no replacement policy or smoke is approved yet.
- exact remaining caller/blocker: blocker is an explicit fixture-role policy
  for RHF that separates route smoke from physics endpoint acceptance.

-- repo-doer@macmini
