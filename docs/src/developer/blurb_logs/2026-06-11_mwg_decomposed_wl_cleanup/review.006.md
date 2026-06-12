Review result:

Accepted as a developer-only oracle comparison. The pass constructed a matched
old nested fixed-block QW/MWG comparator for the side13 He + GTO diagnostic and
kept the old route in its proper role: oracle/reference only, not route
authority.

Main result:

The matched old nested QW/MWG result agrees with the current decomposed
side13 He + GTO final-basis diagnostic to roundoff:

- old nested RHF total: `-2.864394277839952`
- current decomposed RHF total: `-2.864394277839975`
- old-minus-decomposed total: `2.3092638912203256e-14`
- one-electron difference: `1.8207657603852567e-14`
- electron-electron difference: `4.440892098500626e-15`

This is strong evidence that the below-HF side13 result is not a new
decomposed-route convention mismatch. It is consistent with the old QW/MWG
Hamiltonian for the matched fixture.

Paper/Fig. 8 status:

No repo-local Fig. 8 numeric table was found. The `d = 0.3`, `s = 1.0`, side-11
proxy is useful as a constructible old-QW diagnostic, but it is not a
paper-value validation because it has final dimension `422`, not the manager's
recollected `447` basis functions.

Scientific interpretation:

The decomposed route now has the right local oracle evidence. The remaining
question is not "is the new decomposed route wrong?" but "what external or
paper-like target should define He + GTO/MWG acceptance?" The side13 result
should remain diagnostic until that target is chosen.

Validation reviewed:

- `julia --project=. tmp/work/side13_he_gto_old_qw_mwg_oracle_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- artifact: `tmp/work/side13_he_gto_old_qw_mwg_oracle_summary.txt`
- artifact: `tmp/work/side13_he_gto_old_qw_mwg_oracle.tsv`

Deletion/shrinkage review:

No old surface became obsolete. No test was added. The old nested fixed-block
QW/MWG path remains useful as an oracle comparator. The next retirement
decision should wait until He + GTO acceptance semantics are settled.

Commit/push:

Pending manager commit/push of the tracked response/review log.

Next target:

No automatic implementation blurb is queued from this review. The next step is
a manager decision: either locate/enter the Fig. 8 target data, define a
paper-like fixture that reaches the intended 447-function S-only comparison, or
keep He + GTO/MWG diagnostic-only while moving to another physics target.
