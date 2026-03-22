# Cartesian Nested Shell-sequence Coverage Fix

This note records the coverage fix that had to be made before the shell
sequence could be judged fairly.

The parent-orbital capture diagnostic was the right next test, but it also
exposed a concrete source-space bug in the current multi-shell construction:

- the sequence builder enforced disjointness
- but it did not enforce complete coverage of the intended parent working cube

On the previous count-17 three-shell test, the decomposition was effectively:

- shell 1: `13^3 - 11^3`
- shell 2: `11^3 - 9^3`
- shell 3: `9^3 - 7^3`
- direct core: `5^3`

so the whole `7^3 - 5^3` annulus was being dropped.

That violated the intended nesting rule:

- every parent fixed-space row in the working cube must belong to exactly one
  contraction object

before any retained-content policy is evaluated.

## Enforced Coverage Rule

The shell-sequence builder now checks both:

- disjointness
- completeness

For the inferred rectangular working cube of the sequence, every parent row
must belong to exactly one of:

- direct core
- shell face piece
- shell edge piece
- shell corner piece
- or some later contraction object in the same fixed-block source language

The current implementation enforces this by:

- forming the disjoint support union of the core and shell layers
- inferring the bounding working cube from those retained rows
- checking that the support union equals the full parent-row set of that cube

If rows overlap, the sequence still throws as before. If rows are missing, it
now also throws explicitly.

The current corrected complete sequence has:

- working cube `(3:15, 3:15, 3:15)`

and therefore full coverage of the intended `13^3` parent working cube.

## Corrected Multi-shell Test

The corrected complete-shell-sequence test now adds the missing shell layer:

- shell 1: `13^3 - 11^3`
- shell 2: `11^3 - 9^3`
- shell 3: `9^3 - 7^3`
- shell 4: `7^3 - 5^3`
- direct core: `5^3`

So there is no uncovered annulus left between the outer shells and the
retained core.

For reference, the face-only shell sequence is still kept as a deliberately
incomplete comparison object, but it is now marked as such and built with
coverage enforcement disabled.

## Updated He nearest/GGT Scalar Comparison

Same stabilized He fixed-`a` count-17 case:

- `a = 1/4`
- `xmax = 10`
- `s = 0.626026121152214`

Results:

- baseline unnested:
  - fixed dim `4913`
  - overlap error `4.096500916261903e-12`
  - `E1 = -1.998588440029027`
  - `⟨Vee⟩ = 1.248965557065494`

- shell-plus-core:
  - fixed dim `1403`
  - overlap error `2.2354340600827527e-12`
  - `E1 = -1.9985629071304167`
  - `⟨Vee⟩ = 1.2489197930296585`

- face-only four-shell sequence:
  - fixed dim `413`
  - overlap error `1.208411649553688e-14`
  - `E1 = -1.815211525572233`
  - `⟨Vee⟩ = 2.296897646259769`

- corrected complete four-shell sequence:
  - fixed dim `589`
  - overlap error `2.1921484365440793e-12`
  - `E1 = -1.9981842264017424`
  - `⟨Vee⟩ = 1.816303754935915`

So the coverage fix materially improves the complete shell sequence:

- `E1` is now back in the right one-body regime
- much closer to the shell-plus-core / baseline result

but the two-electron scalar is still poor.

## Updated Parent-orbital Capture Comparison

With the corrected complete sequence, the parent low-energy one-body capture
changes drastically.

Ground parent mode, energy `-1.9984765858974547`:

- shell-plus-core:
  - retained `99.999926149879%`
  - projected energy shift `+3.7100e-6`
- face-only sequence:
  - retained `89.100774497625%`
  - shift `+2.345254901576858`
- corrected complete sequence:
  - retained `99.998674577069%`
  - shift `+4.339694670905114e-4`

Next three low-energy modes near `-0.499730664992...`:

- shell-plus-core:
  - retained about `98.92%`
  - shifts about `+0.01`
- face-only sequence:
  - retained about `70.7%`
  - shifts about `+3.1`
- corrected complete sequence:
  - retained about `99.53%`
  - shifts about `+0.008`

Average retained fraction for the first four parent modes:

- shell-plus-core: `99.196065252963%`
- face-only sequence: `75.296206837706%`
- corrected complete sequence: `99.665509173334%`

## Interpretation

The coverage bug was real, and fixing it changes the diagnosis materially.

After the coverage fix:

- the corrected complete shell sequence is much better than the broken result
- it now captures the parent low-energy one-body content extremely well
- and its `E1` is correspondingly restored to the right physical regime

So the remaining problem is no longer well described as a failure to retain the
parent low-energy one-body span.

The remaining problem is now:

- still not consumer plumbing
- no longer simple source incompleteness
- and no longer an obvious one-body low-energy capture failure

The strongest remaining mismatch is in `⟨Vee⟩`.

That means the next diagnostic should probably move to the interaction side:

- pair-term transfer
- Coulomb / interaction packet quality
- or a two-electron retained-subspace diagnostic

rather than immediately changing the retained-count policy again.
That follow-up is now recorded in
[cartesian_nested_interaction_transfer_diagnostic.md](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/docs/cartesian_nested_interaction_transfer_diagnostic.md).

## Conclusion

Yes, the coverage fix was necessary first.

After it is removed:

- the corrected complete shell sequence is materially better than the broken
  result
- the remaining problem is still real
- but it no longer looks like a one-body orbital-capture problem

So retained-content policy may still matter later, but the next diagnostic
target should now be the interaction representation, not another immediate
change to the shell one-body retained counts.
