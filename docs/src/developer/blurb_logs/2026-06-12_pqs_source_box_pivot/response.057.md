Pass 057 complete: fixture-quality audit done; q=9/q=11 gate promotion deferred.

Files edited:

- `docs/src/developer/pqs_near_term_final_basis_realization_plan.md`
  - added a concise status note that the q ladder is route/scaling validation, not accepted physical convergence;
  - recorded that one-shell q=9/q=11 gate promotion is deferred pending fixture review.

Ignored audit artifacts:

- `tmp/work/pqs_complete_core_shell_fixture_quality_audit.jl`
- `tmp/work/pqs_complete_core_shell_fixture_quality_audit_summary.txt`

No production source changed. No permanent test was added.

Fixture-quality table:

```text
q   box        inner      raw dims    radius  central dx  core dx range       shells  final dim  H1(Z=2)       J(Z=2)       RHF total        HF error
5   (1:7)^3    (2:6)^3    (5,5,5)     8.0     0.703230    0.703..2.363       1       223       -1.8476619225  1.0254002949 -2.7213372829  +0.1403427128
7   (1:9)^3    (2:8)^3    (7,7,7)     8.0     0.399447    0.399..2.382       1       561       -1.9599740172  1.2431154890 -2.8100680501  +0.0516119455
9   (1:11)^3   (2:10)^3   (9,9,9)     8.0     0.279568    0.280..2.295       1       1115      -1.9869007472  1.2260752734 -2.8499091618  +0.0117708338
11  (1:13)^3   (2:12)^3   (11,11,11)  8.0     0.216221    0.216..2.182       1       1933      -1.9949555260  1.2455645981 -2.8559475204  +0.0057324752
```

Common fixture facts:

- All four points use `AsinhMapping(a = 0.25, tail_spacing = 10.0)`.
- The mapping `s` changes with q:
  - q=5: `1.8905123339209908`
  - q=7: `1.2997272295706812`
  - q=9: `0.9902683653871858`
  - q=11: `0.7998321412742654`
- All four mapped parent axes have physical endpoints approximately `(-8.0, 8.0)` bohr.
- All four fixtures have exactly one surrounding shell layer.

Direct answers:

- Enough shells for a physical convergence claim: no.
  The ladder changes q and central spacing while keeping exactly one surrounding shell. That validates route scaling and the final/pre-final density convention, but it is not a full fixture-convergence study.

- Promote q=9 or q=11 now: no.
  q=9 is cheap enough to be a candidate later and q=11 is the stronger reference point, but neither should become a permanent gate before fixture review.

- Recommended next fixture direction:
  Vary shell count, mapping/distortion parameters, box radius, and central spacing together. Do not keep increasing q alone and call that physical convergence.

Interpretation:

- The ladder is strong route/scaling validation.
- It supports the accepted pre-final positive-weight density gauge and final-orbital consumption through the combined Lowdin cleanup.
- It does not establish accepted He RHF accuracy because one-shell topology may be the limiting fixture assumption.

Validation:

```text
julia --project=. tmp/work/pqs_complete_core_shell_fixture_quality_audit.jl
  PASS: audit completed and wrote summary artifact.

julia --project=. -e 'using GaussletBases; println("load ok")'
  PASS: load ok

git diff --check
  PASS
```

Deletion/shrinkage report:

- Old/fallback/oracle surface made less necessary: the q ladder should no longer be treated as a pending test-promotion ladder. It is now explicitly route/scaling history.
- Deleted/simplified: no code or test was deleted. The pass was audit/doc-only by design.
- No deletion happened because the manager explicitly requested no permanent test promotion and no production source change.
- No new test was added. The audit is an ignored `tmp/work` probe and does not replace permanent coverage.
- Remaining stale/duplicate surfaces to retire next: any future docs or blurbs that call q=9/q=11 accepted physical gates should be corrected; if a permanent gate is later approved, choose one fixture and avoid preserving the whole q ladder as tests.

-- repo-doer@macmini
