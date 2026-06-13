Pass 192 review

Accepted.

This pass answered the H2 supplement/comparability question without adding H2
code, and it kept the line-budget rule real by deleting one stale metadata
scaffold.

Accepted audit conclusions:

- The visible driver can express H2 as raw atom symbols, charges, and
  coordinates, but it does not yet expose the old H2 WL/QW reference controls
  as first-class driver inputs: `bond_axis`, `bond_length`, `xmax_parallel`,
  and `xmax_transverse`.
- The PQS diatomic route has private route-owned pieces for source plan, final
  basis, H1, and Hamiltonian handoff, but the current visible driver path does
  not yet provide a clean H2 physics endpoint artifact.
- The old WL/QW H2 HF/ED references include the H/cc-pVTZ S/P residual
  supplement. The current PQS driver/source-box path does not carry that
  supplement/residual policy.
- Therefore the correct next implementation choice is Option B: first implement
  a gausslet-only H2 PQS diagnostic with `comparison_ready = false`, and do not
  compare its HF value directly to the old supplemented WL/QW HF total.

This is the right guardrail. It avoids repeating the earlier Be2 mistake where
two available routes were treated as comparable before their basis and
construction policies actually matched.

Accepted future shape:

- A future `test/driver_inputs/h2_pqs_q5_wlmap.jl` can be introduced only after
  the manager chooses the gausslet-only target explicitly.
- It should record `supplement_policy = :none`, `comparison_ready = false`, and
  no direct `wl_rhf_total` comparison.
- A future supplemented H2 endpoint needs explicit driver/PQS support for the
  residual supplement policy before comparing to the documented WL/QW HF/ED
  references.
- Final-basis self-overlap remains a scalar identity-error diagnostic only, not
  downstream working data.

Deletion accepted:

- Deleted `test/nested/cartesian_pair_block_one_body_block_set_summary_runtests.jl`.
- The deleted file was a metadata-only summary shape scaffold for future mixed
  one-body block-set consumption.
- It was not included in the nested/default/integration runners, was not a
  scientific endpoint, and had no source caller.
- Active compact smoke coverage remains elsewhere for the live one-body
  consumer path.

Validation reviewed:

- deleted-file/helper search has no live `src/test/bin` hits;
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed;
- `git diff --check` passed;
- scoped line budget is `0 added / 237 deleted`.

No pass 193 is being issued now because the user requested review, commit,
push, and stop.

-- repo-manager@macmini
