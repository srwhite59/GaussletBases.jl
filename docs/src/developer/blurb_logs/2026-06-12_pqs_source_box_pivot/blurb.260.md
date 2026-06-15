Pass 260 - independent H2 PQS supplement support-tiling audit

Context:
- Current HEAD should include
  `4ec57e95 Record independent H2 PQS provider seam audit`.
- Pass 259 identified the key provider-block hazard: existing CPB provider
  functions operate on rectangular CPBs, while independent H2 shared-shell
  support is outer-box minus inner-box. Provider calls must therefore be preceded
  by explicit support CPB tiling/row ownership or an equivalent row map.

Task:
Do a focused no-implementation audit of support tiling/row ownership for the
future independent H2 PQS supplement provider-block payload.

Questions to answer:
1. Where in the current independent H2 source-plan/final-basis payloads are the
   support rows, unit keys, source boxes, and parent row indices available?
2. Can atom-contact core support be decomposed into rectangular CPB/provider
   tiles directly?
   - Expected conceptual pieces are atom-local core supports plus midpoint/contact
     slab, but report the actual code objects and counts.
3. Can `shared_shell_1` and `shared_shell_2` support be decomposed into
   rectangular CPB/provider tiles without using filled outer source CPBs as if
   they were the support?
   - Identify outer box, inner box, shell support rows, and any existing
     terminal shellification or region objects that already encode this.
4. What exact row-map/fingerprint should a first support-partition payload
   expose?
   - unit keys;
   - tile count;
   - per-tile support count;
   - total parent-row coverage;
   - duplicate/missing/outside row counts;
   - source unit row ranges or parent row ids;
   - retained range association if available.
5. Is a support-partition summary implementable before provider blocks?
   - If yes, give the exact helper/object name and next-pass scope.
   - If no, state the exact missing source-plan data.
6. What is the smallest validation smoke for the support partition?
   - Prefer no full H2 route if existing payload fixtures can expose the facts.
   - If the full route is the only way, explain why and keep it focused.

Strict exclusions:
- Do not implement support partition code in this pass unless the audit shows
  it is a trivial deletion-only/no-risk extraction. Prefer audit only.
- Do not call CPB provider functions.
- Do not implement provider blocks, mixed matrices, residual MWG representation,
  combined density readiness, supplemented values, CR2/export, or public API.
- Do not use fake-PQS/WL source-backed data as independent-PQS evidence.

Validation:
- `git diff --check`.
- No Julia run required for no-edit audit.
- If you make a tiny cleanup offset, run parse smoke for the touched test.

Line budget:
- Scoped `src + test + bin` impact `0` is acceptable for this audit.
- If any source/test/bin file is touched, keep scoped net-negative.

Report:
- Direct answers to the six audit questions.
- Recommended exact next implementation pass.
- Whether a support-partition summary is ready to implement.
- Validation command(s).
- Scoped line count for `src + test + bin`.
- Deletion/shrinkage result:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
