Accepted with an important performance/design caveat.

Pass 095 does the intended functional work: the one-center `:pqs_source_box`
driver path now builds the complete core/shell final basis and H1 payload from
the route-owned source plan. The H1/J diagnostic now carries a final dimension
and H1 value while remaining blocked only on the density/J inputs:
`axis_weights` and `raw_pair_factor_terms`.

The implementation stayed inside the requested boundary:

- no RHF, SCF, Fock construction, density iteration, GTO, exports, artifacts,
  fixture promotion, explicit-box authority, or fixed-block authority;
- no permanent tests were added;
- the public report still exposes the existing compact H1/J summary fields,
  not dense final-basis/H1 objects.

The main caveat is cost shape. A normal one-center PQS driver dry-run now
performs real final-basis and H1 construction, including Lowdin/final-basis
work and electron-nuclear H1 construction. That explains the longer validation
times. This is acceptable as a temporary diagnostic seam, but it should not
become ordinary always-on driver behavior without an explicit diagnostic or
materialization request.

Follow-up rule:

- Use a tiny smoke while editing route seams.
- After the seam works, run one modest representative check.
- Do not jump straight to broad route-driver/report/materialization tests or
  a large physics fixture.

Next design cleanup should turn the current local assembly work into a compact
route object or explicitly requested diagnostic payload. `pqs_source_box_route_driver_helpers.jl`
should call that object and expose summaries; it should not keep accumulating
complete-core-shell route logic and scalar report fields.

Validation accepted:

- doer's focused one-center PQS source-box dry-run;
- doer's `test/nested/pqs_direct_retained_final_h1_runtests.jl` run;
- manager `git diff --check`;
- manager `/Users/srw/Dropbox/codexhome/cjulia -e 'using GaussletBases; println("load ok")'`.

The first manager plain-`julia` load check hit the known Juliaup lockfile
permission issue, not a code failure; the `cjulia` load check passed.

Deletion/shrinkage review:

- No old code became removable yet; this is still a connection seam.
- No tests were added, which is correct.
- The next stale pressure is fixture-local complete core/shell H1 construction:
  once the driver-owned seam is stable and optionally gated, keep fixture code
  as oracle coverage only or shrink it.

Do not resume the baton loop automatically after this review.

-- repo-manager@macmini
