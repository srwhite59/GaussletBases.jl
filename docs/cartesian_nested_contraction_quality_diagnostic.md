# Cartesian Nested Contraction Quality Diagnostic

This note records the first internal approximation-quality diagnostic on the
current complete nonrecursive shell language.

This local diagnostic was run before the later shell-sequence coverage fix.
Its whole-sequence numbers therefore belong to the earlier three-shell test,
not the corrected four-shell sequence in
[cartesian_nested_sequence_coverage_fix.md](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/docs/cartesian_nested_sequence_coverage_fix.md).
The point of the note is still the same: local smooth-content loss did not look
like the main failure mode.

The complete shell layer fixed the representation-completeness issue:

- faces
- edges
- corners
- interior/core

are all now present as explicit contraction objects.

But the resulting nearest/GGT physics is still poor. So the next question is
not support completeness and not consumer plumbing. The question is whether the
local retained subspaces are already underresolving the smooth low-order
content.

## Diagnostic Definition

For one local contraction piece:

- let `C_p` be the parent-to-retained coefficient matrix restricted to that
  piece support
- let `S_p` be the parent overlap restricted to that same support

The retained-space projector in the local parent metric is

```math
P_p = C_p (C_p^T S_p C_p)^{-1} C_p^T S_p
```

The first target is the exact local constant-like parent-space vector:

- in the current orthonormalized PGDG parent basis, this is the vector of
  parent integral weights on that piece
- for one 3D parent row `(i_x, i_y, i_z)`, the coefficient is
  `w_{i_x} w_{i_y} w_{i_z}`

The retained constant fraction is

```math
\eta_{\mathrm{const}} =
\frac{v_{\mathrm{const}}^T S_p P_p v_{\mathrm{const}}}
     {v_{\mathrm{const}}^T S_p v_{\mathrm{const}}}
```

and the constant lost fraction is `1 - η_const`.

## Test Case

The diagnostic was run on the current stabilized He fixed-`a` count-17
complete-shell sequence:

- shell layer 1: annulus `13^3 - 11^3`
- shell layer 2: annulus `11^3 - 9^3`
- shell layer 3: annulus `9^3 - 7^3`
- retained direct core: `5^3`

The scratch script is
[nested_contraction_quality_diagnostic.jl](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/tmp/work/nested_contraction_quality_diagnostic.jl).

## Constant-like Retained Percentages

The constant-like result is essentially exact for every local piece.

Layer 1:

- faces: `100.000000000000%` retained, `~0%` lost
- edges: `100.000000000000%` retained, `~0%` lost
- corners: `100.000000000000%` retained, `0%` lost
- whole shell layer: `100.000000000000%` retained

Layer 2:

- faces: `100.000000000000%` retained, `~0%` lost
- edges: `100.000000000000%` retained, `~0%` lost
- corners: `100.000000000000%` retained, `0%` lost
- whole shell layer: `100.000000000000%` retained

Layer 3:

- faces: `100.000000000000%` retained, `~0%` lost
- edges: `100.000000000000%` retained, `~0%` lost
- corners: `100.000000000000%` retained, `0%` lost
- whole shell layer: `100.000000000000%` retained

Core:

- core block: `100.000000000000%` retained

Whole complete-shell sequence:

- `100.000000000000%` retained

So the constant-like direction does not reveal the problem clearly. It is
retained essentially exactly by construction.

## Minimal Second Smooth Probe

Because the constant result was unexpectedly good, I added one very smooth
second probe:

- a broad sampled Gaussian proxy
- on each parent row `(i_x, i_y, i_z)`,
  `v_g = w_{i_x} w_{i_y} w_{i_z} exp(-(x^2 + y^2 + z^2)/(2 sigma^2))`
- with `sigma = 4`

This is not an exact coefficient vector for a Gaussian function; it is a smooth
weight-times-value proxy. It is only meant as a very mild second check.

Even this probe is retained at a very high level:

Layer 1:

- faces: `99.999557146991%` retained
- edges: `99.999778573250%` retained
- corners: `100%` retained
- whole shell layer: `99.999610142419%` retained

Layer 2:

- faces: `99.999995656507%` retained
- edges: `99.999997828253%` retained
- corners: `100%` retained
- whole shell layer: `99.999996614559%` retained

Layer 3:

- faces: `99.999999970679%` retained
- edges: `99.999999985340%` retained
- corners: `100%` retained
- whole shell layer: `99.999999978894%` retained

Core:

- core block: `100%` retained

Whole complete-shell sequence:

- `99.999798314502%` retained

## Interpretation

The constant-like diagnostic does not show an obvious local contraction
failure. The second smooth probe also does not show one.

So the current face/edge/core retained counts are not obviously too small for
the very smooth low-order local content.

What this means in practice is:

- the poor physics is not explained by losing the local constant direction
- it is not even explained by a very mild broad smooth target
- the first obvious offenders are not visible in these local diagnostics

If one still wants a weak ranking, the outermost shell faces are the least good
on the broad Gaussian proxy, followed by the outermost shell edges, but their
losses are still only at the `10^-6` to `10^-5` percent level.

## Conclusion

This diagnostic is still useful because it narrows the target:

- the current complete shell language is not obviously failing on the lowest
  smooth local content
- so the next policy change should probably be guided by a stronger diagnostic
  than the constant direction

The most natural next diagnostics are:

- exact low-order moment transfer beyond the constant direction
- or, more directly, projection of the parent-space low-energy orbital content
  onto the nested retained blocks

That is likely where the real approximation loss now lives.
