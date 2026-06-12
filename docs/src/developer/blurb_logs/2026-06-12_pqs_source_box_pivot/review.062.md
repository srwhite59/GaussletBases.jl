Review 062: accepted as a probe result, not an acceptance gate.

The side-13 multi-layer PQS He RHF result is coherent and useful:

- final dimension: 1549
- final overlap identity error: about `5.51e-13`
- Z = 2 H1: `-1.9755618232013417`
- self-Coulomb J: `1.2169264388860319`
- RHF total: `-2.8372556463894707`
- error vs He HF reference: about `+0.02442434922276826` Hartree
- delta vs WL side-13 RHF: about `-0.0007576466885570454` Hartree

The density convention is still framed correctly: final one-electron H, pre-final
positive-weight density interaction, final orbitals pulled back through
`combined_lowdin_cleanup`, no signed-final-weight route, no raw no-division route,
and no fixed-block pair authority.

This should not be promoted as a permanent gate yet. The fixture relation among
box radius, core spacing `d`, mapping/distortion `s`, direct core size, and shell
depth remains unsettled. The next useful step is a probe-only fixed-parent ladder:
hold the side-13, `R ~= 8`, `d = 0.1`, `s = 1.0` parent fixed and vary direct core
size/shell depth together. That will expose whether the current q/shell choice is
special without spending a full design pass on a general `Z,d,s,ns` rule.

Deletion/shrinkage: no source/test cleanup was expected in pass 062. The result
does reduce the need to treat old fixed-block pair data as route authority, but
those surfaces remain useful oracle/reference material until a compact PQS gate
is selected.

-- repo-manager@macmini
