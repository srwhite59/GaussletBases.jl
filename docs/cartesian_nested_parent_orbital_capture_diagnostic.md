# Cartesian Nested Parent-Orbital Capture Diagnostic

This note records the next diagnostic after the local smooth-content test.

This first orbital-capture pass was run before the later shell-sequence
coverage fix. So its complete-shell-sequence numbers are now superseded by
[cartesian_nested_sequence_coverage_fix.md](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/docs/cartesian_nested_sequence_coverage_fix.md), which corrects the missing `7^3 - 5^3` annulus and reruns the same comparison.

The previous pass ruled out the simplest local-failure story:

- the constant-like direction is retained essentially exactly
- even a broad smooth Gaussian proxy is retained at about `99.9998%`
  globally

So the poor He nearest/GGT scalar behavior is not explained by losing the
lowest smooth local content.

That makes the next useful diagnostic more global and more physical:

- build the parent fixed-block one-body problem
- take its lowest one-body orbitals
- project those parent orbitals into the nested retained spaces
- measure how much physically relevant parent content survives

## Diagnostic Definition

On the parent fixed block, let:

- `S_parent` be the parent overlap matrix
- `H_parent` be the parent fixed-block one-body Hamiltonian

The parent low-energy orbitals are taken from the generalized eigenproblem

```math
H_{\mathrm{parent}} u_i = \varepsilon_i S_{\mathrm{parent}} u_i
```

with the usual `S_parent`-orthonormalization.

For one retained nested space with parent-to-retained map `C` and retained
metric

```math
G = C^T S_{\mathrm{parent}} C
```

the whole-space parent-metric projector is

```math
P = C G^{-1} C^T S_{\mathrm{parent}}.
```

For one parent orbital `u_i`, the retained norm fraction is

```math
\eta_i = \frac{u_i^T S_{\mathrm{parent}} P u_i}
               {u_i^T S_{\mathrm{parent}} u_i}.
```

The lost fraction is `1 - η_i`.

I also computed the Rayleigh quotient of the projected orbital:

```math
\varepsilon_i^{\mathrm{proj}} =
\frac{(P u_i)^T H_{\mathrm{parent}} (P u_i)}
     {(P u_i)^T S_{\mathrm{parent}} (P u_i)}
```

and reported the shift

```math
\Delta \varepsilon_i = \varepsilon_i^{\mathrm{proj}} - \varepsilon_i.
```

The scratch driver is
[nested_parent_orbital_capture_diagnostic.jl](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/tmp/work/nested_parent_orbital_capture_diagnostic.jl).

## Test Case

The same stabilized He fixed-`a` count-17 case was used:

- `a = 1/4`
- `xmax = 10`
- `s = 0.626026121152214`
- base QW-PGDG fixed line
- nearest/GGT nested comparisons

The parent fixed block has dimension `4913`, and its overlap is already
identity-level:

- `||S_parent - I||_∞ = 1.2084116495536886e-14`

The retained spaces compared here are:

- shell-plus-core
- face-only three-shell sequence
- complete three-shell sequence

## First Five Parent One-body Orbitals

Mode 1, parent energy `-1.9984765858974547`:

- shell-plus-core:
  - retained `99.999926149879%`
  - lost `7.3850e-7`
  - projected energy `-1.9984728758646715`
  - shift `+3.7100e-6`
- face-only sequence:
  - retained `61.750482431713%`
  - lost `38.249517568287%`
  - projected energy `3.9536495609102995`
  - shift `+5.952126146807754`
- complete sequence:
  - retained `62.963911660281%`
  - lost `37.036088339719%`
  - projected energy `3.678074537518059`
  - shift `+5.676551123415514`

Mode 2, parent energy `-0.49973066499238866`:

- shell-plus-core:
  - retained `98.916503544796%`
  - lost `1.083496455204%`
  - shift `+0.011072502552376406`
- face-only sequence:
  - retained `65.733919317610%`
  - lost `34.266080682390%`
  - shift `+3.0850700743409796`
- complete sequence:
  - retained `89.063750187034%`
  - lost `10.936249812966%`
  - shift `+0.9302010694930936`

Mode 3, parent energy `-0.4997306649923789`:

- shell-plus-core:
  - retained `98.928111620657%`
  - lost `1.071888379343%`
  - shift `+0.01050009538362584`
- face-only sequence:
  - retained `65.760147664123%`
  - lost `34.239852335877%`
  - shift `+3.0648769192355445`
- complete sequence:
  - retained `89.089978533547%`
  - lost `10.910021466453%`
  - shift `+0.9277435463799071`

Mode 4, parent energy `-0.4997306649923544`:

- shell-plus-core:
  - retained `98.939719696518%`
  - lost `1.060280303482%`
  - shift `+0.00992782252990726`
- face-only sequence:
  - retained `65.786376010635%`
  - lost `34.213623989365%`
  - shift `+3.0446998657340942`
- complete sequence:
  - retained `89.116206880059%`
  - lost `10.883793119941%`
  - shift `+0.9252874698446781`

Mode 5, parent energy `-0.4994888533519622`:

- shell-plus-core:
  - retained `97.792083720871%`
  - lost `2.207916279129%`
  - shift `+0.01828706608583014`
- face-only sequence:
  - retained `71.735922785365%`
  - lost `28.264077214635%`
  - shift `+1.769956345992213`
- complete sequence:
  - retained `98.424393995627%`
  - lost `1.575606004373%`
  - shift `+0.14992334200871443`

## Low-energy Subspace Summary

Average retained fraction for the first four parent one-body modes:

- shell-plus-core: `99.196065252963%`
- face-only sequence: `64.757731356020%`
- complete sequence: `82.558461815230%`

## Interpretation

This diagnostic does identify the real problem clearly.

The bad nested He scalars are not coming from an obvious failure on local
constant-like or very broad smooth content. They are coming from failure to
retain the parent low-energy one-body orbital content of the full fixed block.

The ranking is clear:

- shell-plus-core is physically good because it keeps the parent low-energy
  content almost intact
- face-only sequence is the worst by far
- complete sequence is better than face-only, especially on the excited
  low-energy modes, but it still loses far too much of the parent ground
  orbital and enough of the next low-energy modes to spoil the physics

The most important single number is probably the parent ground-orbital capture:

- shell-plus-core: `99.9999%`
- face-only sequence: `61.75%`
- complete sequence: `62.96%`

So the current complete shell sequence is still not a viable replacement for
the physical low-energy fixed space, even though it is geometrically complete.

## Policy Consequence

This points to the next fix much more directly than the earlier smooth-content
diagnostic:

- the next issue is not support completeness
- and it is not consumer plumbing
- it is the retained-content policy for the nested shell sequence as a whole

In practical terms, the next improvement should target the global low-energy
parent content:

- either more retained content on the inner shell hierarchy
- or a better shell/core replacement rule that preserves the parent low-energy
  orbital span
- or a direct low-energy-informed contraction policy

The current shell-plus-core result is the proof of principle:

- the consumer model is fine
- the fixed-block language is fine
- but the grow-and-replace shell sequence is still discarding too much of the
  physical low-energy fixed-space content

So this pass does identify where the bad scalar behavior is coming from.
