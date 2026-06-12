Review result:

Accepted as a developer-only H2 `R = 4.0` old nested/QW restricted HF
diagnostic. The probe stayed inside `tmp/work`, changed no production code,
and added no tests. It used H/cc-pVTZ with `lmax = 1`, so the old molecular
S+P supplement path was exercised as part of the reference route.

Main result:

- default complete rectangular route:
  - final dimension `481`
  - residual count `18`
  - RHF total `-0.9109382643524664`
  - documented total `-0.910938264352`
  - difference `-4.664046926450283e-13`
- endcap/panel `q = 4`, `L = 4` route:
  - final dimension `461`
  - residual count `18`
  - RHF total `-0.9109773150033322`
  - documented total `-0.910977315003`
  - difference `-3.3217872896784684e-13`

Scientific interpretation:

This is enough to treat the old nested/QW diatomic HF path as a trusted oracle
for this line. We do not need a broader old-route H2 sweep right now. The more
useful next question is whether the newer decomposed/final-basis machinery can
handle atomic S+P GTO residuals, since the active decomposed supplement
acceptance has so far mostly used S-only fixtures.

Validation reviewed:

- `julia --project=. tmp/work/h2_r4_qw_hf_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- artifact: `tmp/work/h2_r4_qw_hf_summary.txt`
- artifact: `tmp/work/h2_r4_qw_hf.tsv`

Deletion/shrinkage review:

No old surface became obsolete from this probe. No test was added. The old
nested/QW route remains useful as oracle/comparison infrastructure, but it
should not become another routine proof target unless a specific mismatch needs
localization.

Commit/push:

Pending manager commit/push of this tracked response/review log and the next
published blurb.

Next target:

Use a Be atom with S+P GTO supplement as the next physics probe. The old
nested/QW atomic route should be used as a trusted oracle; the live target is
the newer decomposed/final-basis S+P residual behavior and any precise blocker
to valid Be RHF.
