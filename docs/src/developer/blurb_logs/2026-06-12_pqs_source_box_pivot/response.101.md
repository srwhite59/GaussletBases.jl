Pass 101 response

Local-only timing/design probe completed.

Files changed:

- added ignored local probe:
  - `tmp/work/pqs_driver_h1j_smoke_timing_probe.jl`
  - `tmp/work/pqs_driver_h1j_smoke_timing_probe_summary.txt`
- no source, test, docs, committed fixtures, exports, artifacts, or production
  route behavior changed before this response/log.

Timing commands and results:

- `julia --project=. tmp/work/pqs_driver_h1j_smoke_timing_probe.jl`
  - first corrected fresh-process dry-run timing, measured inside Julia with
    `@elapsed`: `120.102067` s
  - second run in the same Julia process: `0.210283` s
  - first/second ratio: `571.145938`

Note: an initial version of the local probe completed the heavy first dry-run
but then failed while trying to read `summary.missing_inputs` from the
materialized summary. I fixed that probe-only assumption and reran the same
accepted route shape. No source or test code was touched.

Measured accepted dry-run facts:

- status: `materialized_pqs_multilayer_complete_core_shell_h1_j_payload`
- blocker: `nothing`
- final dimension: `223`
- H1 energy: `-5.6629907690725245`
- self-Coulomb: `1.8691288063594704`
- density gauge: `pre_final_localized_positive_weight`
- driver route materialized: `true`
- signed final weight division used: `false`
- density-normalized pair terms used as authority: `false`
- RHF/GTO/exports/artifacts materialized: all `false`

Warm-vs-cold conclusion:

- The accepted dry-run is overwhelmingly first-run compilation dominated.
- Repeated route construction/runtime for the same accepted shape is about
  `0.21` s in-process after compilation.
- The 120 s cold cost is therefore not evidence that the fixture geometry or
  complete core/shell final-basis construction is intrinsically too slow for
  repeated local validation.

Permanent smoke recommendation:

- Do not add this as a normal cold-process tracked test yet.
- A future tracked smoke is technically viable only if it runs in an already
  warm process or is marked slow/manual/local, because a cold `julia
  --project=.` invocation would add roughly two minutes.
- Best next shape, if approved later: a compact generated/local-only probe or
  manual slow smoke over the accepted one-center source-box route, asserting
  only the compact report fields above plus the non-promotion flags. Do not
  broaden it into a full route/global Hamiltonian/export check.
- If a regular tracked smoke is desired, make route construction/first-call
  compilation cheaper first or place it where suite amortization already warms
  the relevant methods.

Smaller lattice/design conclusion:

- I did not run a smaller-lattice variant.
- Smaller source-box geometry would weaken this route contract too much for
  the intended smoke: the accepted path currently verifies the same `q=5/n_s=5`
  one-center source-box geometry and final dimension `223` used in the focused
  H1 seam. Reducing that would turn the check into a toy-only exercise of the
  function path, not the accepted diagnostic route shape.
- Since warm repeated runtime is already about `0.21` s, smaller lattice is not
  the right lever for this future smoke. The blocker is cold compilation /
  first-call latency.

Validation/status:

- `git status --short --branch`

```text
## main...origin/main
```

- `git check-ignore -v tmp/work/pqs_driver_h1j_smoke_timing_probe.jl tmp/work/pqs_driver_h1j_smoke_timing_probe_summary.txt`

```text
.gitignore:12:tmp/ tmp/work/pqs_driver_h1j_smoke_timing_probe.jl
.gitignore:12:tmp/ tmp/work/pqs_driver_h1j_smoke_timing_probe_summary.txt
```

Deletion/shrinkage forecast:

- deleted: none.
- simplified: none in tracked code/tests; the probe design shows a future smoke
  can be compact report-field assertions rather than broad fixture assertions.
- quarantined: the timing probe and summary are ignored local artifacts under
  `tmp/work/`.
- not deleted because: the local probe summary is useful for manager review of
  the future smoke decision; no tracked cleanup was assigned in this pass.
- exact remaining caller/blocker: no source caller for the local probe; the
  remaining blocker for a permanent tracked smoke is cold-process first-call
  latency of about 120 s for the accepted route shape.

-- repo-doer@macmini
