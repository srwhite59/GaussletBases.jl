# Cartesian Nested Fixed-nside Compression Policy

This note records the first explicit nonrecursive replacement policy on top of
the working shell-sequence fixed-block interface.

This pass sits on the earlier face-only shell-sequence line. It is still a
useful negative result about the first compression policy, but it should not be
read as the later coverage-fixed complete-shell state.

For the current supported non-atomic recursive nested routes, the old
compressed complete-shell retain default `(4,3)` / `3` is no longer the
canonical source-layer behavior. The default route now derives complete-shell
retention directly from `nside`, and any future non-`nside` shell-resolution
mode should remain explicit rather than hidden behind the default path.

The shell-sequence interface by itself is structurally valid, but it is not yet
a useful hierarchy if each added shell only grows the fixed block. The missing
ingredient is an explicit rule for what direct/core content gets replaced once a
new shell layer is added.

## First Policy

The first policy now follows the legacy fixed-`nside` interpretation:

- choose `nside = 5`
- add shell layers so the remaining raw interior is reduced by one parent site
  on each side after each shell
- continue until the remaining raw interior is a `5 x 5 x 5` cube
- at that point, do not keep a direct `5^3` core as the final retained block
- instead, apply the same local `doside`/COMX contraction on the final
  one-dimensional length-`5` intervals
- this gives a retained inner block of length `3` per axis, i.e. a `3 x 3 x 3`
  contracted core inside the faces of the raw `5 x 5 x 5` cube

So for `nside = 5`, the intended stopping rule is:

- stop shrinking once the raw interior first reaches size `5`
- then replace that raw `5^3` interior by the contracted `3^3` inner block
- no further shrinking beyond that point

This is different from the earlier incorrect interpretation where the direct
core itself was kept at size `5^3`.

## Current First Test

On the stabilized He fixed-`a` count-17 line:

- baseline unnested fixed block: dimension `4913`
- grow-only three-shell sequence with raw `5^3` core: dimension `341`
- fixed-`nside = 5` sequence with contracted `3^3` core: dimension `243`

The overlap quality stays clean:

- grow-only overlap error: about `1.21e-14`
- fixed-`nside` overlap error: about `3.55e-15`

But the first physics result is not yet acceptable:

- baseline: `E1 = -1.998588440029027`, `⟨Vee⟩ = 1.248965557065494`
- grow-only three-shell: `E1 = -1.6437270930603496`, `⟨Vee⟩ = 4.4067758554156065`
- fixed-`nside = 5`: `E1 = -1.6148380896530632`, `⟨Vee⟩ = 24.694908673599514`

The timing is favorable but not physically meaningful yet:

- baseline constructor time: about `8.18 s`
- grow-only three-shell constructor time: about `1.17 s`
- fixed-`nside = 5` constructor time: about `1.08 s`

## Conclusion

This pass establishes the first explicit compression policy in code, and it now
matches the intended legacy `nside = 5` meaning. But it is not yet the first
genuinely useful hierarchy policy.

The reason is not structural failure:

- the shell-sequence interface still works
- the fixed-block consumer is unchanged
- the compressed blocks remain orthonormal and end-to-end runnable

The reason is physical failure:

- even the grow-only `5^3` core version loses the interior physics badly
- the contracted `3^3` core version loses it even more

So the next step is not another consumer rewrite. It is a better compression
policy for what is kept in the innermost block before true recursion begins.
