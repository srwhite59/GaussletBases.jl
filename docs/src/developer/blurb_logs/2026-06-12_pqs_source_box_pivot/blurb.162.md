Pass 162 - choose the next downstream Hamiltonian contract

Role: repo-doer@macmini

Task type: no-edit design/decision audit. Do not edit files, do not commit.

Purpose:

Choose the next downstream Hamiltonian contract after the private Be2/PQS
consumer contract became available and readiness moved to:

```text
:missing_hfdmrg_density_density_contract
```

The next implementation should be chosen deliberately. Do not default to adding
another payload layer.

Line-budget rule for future implementation:

Any subsequent pass that edits `src/` or `test/` in this lane remains subject to
the corrective line-budget rule:

```text
git diff --numstat -- src test
sum(deleted) > sum(added)
```

If implementing the chosen downstream contract cannot be line-negative, the
next implementation pass must either get an explicit manager exception in
advance or write `.agent_handoffs/ATTENTION.md` and stop.

Read first:

- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/review.161.md`
- `docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.158.md`
- `src/pqs_source_box_diatomic_complete_core_shell.jl`
- `test/nested/pqs_source_box_route_driver_be2_ham_payload_fingerprint_runtests.jl`
- minimal relevant HFDMRG/CR2 files from the pass-158 audit, if needed

Decision options:

1. HFDMRG density-density `H,V` contract.

   This would define how to convert or expose the current Be2/PQS handoff as a
   reviewed final-space `H::N x N`, `V::N x N` density-density object suitable
   for HFDMRG's density-density backend. This is probably the closest route to
   a useful downstream comparison, but it requires reviewing the final-space
   two-body convention before claiming readiness.

2. CR2 read-only inspector handoff.

   This would produce no solver-ready `V`, but would expose the current private
   Be2/PQS handoff in a compact inspectable form for the CR2 agent to compare
   dimensions, H1 spectrum, nuclear metadata, electron/spin convention, and
   missing downstream blockers against the WL route. This may be useful sooner
   and avoids premature `V` claims.

3. HamV6/sliced-integral/export contract.

   This is likely too broad now because no `V6`/`Vblocks` or dense four-index
   object is materialized.

4. More cleanup before downstream work.

   If the current diatomic source file still has obvious stale surfaces or
   repeated nonclaim flags that can be line-negatively removed, recommend that
   before downstream contract implementation.

Audit questions:

- Which option most directly advances the medium-term goal: CR2 comparing Be2
  WL and PQS?
- What exact new contract would be implemented?
- What existing source/test surface would it replace or shrink to satisfy the
  line-budget rule?
- What physics/numerical convention must be verified before saying `ready`?
- What should remain explicitly false or blocked?
- What is the smallest focused validation for that next pass?

Forbidden:

- No source edits.
- No test edits.
- No Julia tests.
- No CR2/HFDMRG runs.
- No dense `Vee`, final-space `V`, `V6`, `Vblocks`, export/HamV6/JLD2,
  artifacts, public API, H1/J, RHF/SCF, WL run, or fixture promotion.

Validation:

Run only:

```text
git status --short --branch
```

Use `rg`, `sed`, and small file reads as needed.

Report back:

- recommended next contract option and why;
- options rejected/deferred and why;
- exact implementation seam for the next pass;
- line-budget strategy for the next implementation pass;
- expected status/blocker vocabulary;
- expected validation;
- git status;
- deletion/shrinkage report:
  - deleted: none in this audit
  - simplified: none in this audit
  - quarantined: none in this audit
  - not deleted because: audit only
  - exact remaining caller/blocker:

-- repo-manager@macmini
