Review result:

Accepted as a developer-only `n_s = 7` He RHF accuracy diagnostic. The pass
used the old nested/QW MWG path as an oracle/comparison route, loaded AHGBS-9
S-only from the external GaussletModules basis file, and did not add tests or
change production code.

Main result:

All three `n_s = 7` points are sub-mH relative to their Fig. 8 rows:

- `d = 0.15`: RHF `-2.8616625064826398`, Fig. 8 error
  `+1.684323533845955e-5 Ha`
- `d = 0.10`: RHF `-2.861673961528321`, Fig. 8 error
  `+1.716095156645281e-6 Ha`
- `d = 0.20`: RHF `-2.861639584273142`, Fig. 8 error
  `+4.379386004993435e-5 Ha`

The `d = 0.10` point is the closest in this repo probe, both to Fig. 8 and to
the He HF reference. The trend does not exactly reproduce the Fig. 8 ordering,
where `d = 0.15` is the smallest positive-error point, but the discrepancies
are now in the microhartree-to-tens-of-microhartree range rather than mH.

Scientific interpretation:

`n_s = 7` is accurate enough to close the immediate atomic He accuracy check.
The useful next question is no longer whether the old nested/QW path can reach
the paper-like He accuracy; it can. The next physics target should move to a
diatomic HF problem rather than continuing He reproduction details.

Validation reviewed:

- `julia --project=. tmp/work/fig8_he_ns7_rhf_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- artifact: `tmp/work/fig8_he_ns7_rhf_summary.txt`
- artifact: `tmp/work/fig8_he_ns7_rhf.tsv`

Deletion/shrinkage review:

No old surface became obsolete. No test was added. The old nested/QW MWG route
still remains useful as an oracle comparator.

Commit/push:

Pending manager commit/push of the tracked response/review log.

Next target:

Move to diatomic H2 at the restricted closed-shell HF level. Use the documented
old nested/QW H2 `R = 4.0` route as the first reference diagnostic before
trying to generalize the decomposed route.
