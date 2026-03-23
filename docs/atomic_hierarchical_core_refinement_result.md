# Atomic Hierarchical Core-only Refinement Result

This pass implements and tests the first hierarchical core-only refinement
inside the trusted corrected complete-shell atomic anchor.

## Scope

Atomic-only, on the same active hybrid footing:

- `legacy_atomic_gaussian_supplement("He", "cc-pVTZ"; lmax = 0)`

The outer trusted complete-shell object stays fixed. Only the retained direct
`5^3` core is changed.

The fixed-block consumer stays unchanged.

## What Was Refined

Starting trusted reduced atomic anchor:

- four complete outer shell layers
- direct retained `5^3` core
- total fixed dimension `589`

The new core-only refinement replaces only that direct `5^3` core by one local
hierarchical object built from the original parent-space functions assigned to
the `5^3` region:

- one inner complete shell layer on the `5^3 - 3^3` annulus
- face retained counts `(2, 2)` on each tangential pair
- edge retained count `2` on each free axis
- corners kept direct
- direct retained sub-core `3^3`

So the local refined core has:

- shell columns `6 * 4 + 12 * 2 + 8 = 56`
- direct sub-core columns `27`
- total refined-core dimension `83`

and the whole outer-plus-core candidate becomes:

- `464 + 83 = 547`

instead of the current trusted:

- `464 + 125 = 589`

## Geometry Reporting

All ranges below use the repo’s current positive parent-index convention.

### Count-17 Primary Case

He, fixed-`a`, `a = 1/4`, `xmax = 10`, `s = 0.626026121152214`.

Trusted outer complete-shell layout:

- shell 1:
  - outer parent box: `ix = iy = iz = 3:15`
  - inner retained box: `ix = iy = iz = 4:14`
  - reference-center range: `(-6.0, 6.0)` outer, `(-5.0, 5.0)` inner
  - physical-center range: `(-4.1266299562, 4.1266299562)` outer,
    `(-2.4470817965, 2.4470817965)` inner
- shell 2:
  - outer parent box: `4:14`
  - inner retained box: `5:13`
  - reference-center range: `(-5.0, 5.0)` outer, `(-4.0, 4.0)` inner
  - physical-center range: `(-2.4470817965, 2.4470817965)` outer,
    `(-1.3904553657, 1.3904553657)` inner
- shell 3:
  - outer parent box: `5:13`
  - inner retained box: `6:12`
  - reference-center range: `(-4.0, 4.0)` outer, `(-3.0, 3.0)` inner
  - physical-center range: `(-1.3904553657, 1.3904553657)` outer,
    `(-0.7596040019, 0.7596040019)` inner
- shell 4:
  - outer parent box: `6:12`
  - inner retained box: `7:11`
  - reference-center range: `(-3.0, 3.0)` outer, `(-2.0, 2.0)` inner
  - physical-center range: `(-0.7596040019, 0.7596040019)` outer,
    `(-0.3900208287, 0.3900208287)` inner
- retained direct core:
  - parent box: `ix = iy = iz = 7:11`
  - reference-center range: `(-2.0, 2.0)`
  - physical-center range: `(-0.3900208287, 0.3900208287)`

New local refinement inside the `5^3` core:

- local outer box: `7:11`
- local inner retained box / direct sub-core: `8:10`
- reference-center range: `(-2.0, 2.0)` outer, `(-1.0, 1.0)` inner
- physical-center range: `(-0.3900208287, 0.3900208287)` outer,
  `(-0.1638565613, 0.1638565613)` inner
- trusted complete-shell working box: `(3:15, 3:15, 3:15)`
- local refined-core working box: `(7:11, 7:11, 7:11)`
- combined refined-sequence working box: `(3:15, 3:15, 3:15)`

### Count-15 Robustness Case

He, fixed-`a`, `a = 1/4`, `xmax = 10`, `s = 0.7303638080109164`.

Trusted outer complete-shell layout:

- shell 1:
  - outer parent box: `2:14`
  - inner retained box: `3:13`
  - reference-center range: `(-6.0, 6.0)` outer, `(-5.0, 5.0)` inner
  - physical-center range: `(-6.3072080399, 6.3072080399)` outer,
    `(-3.6786608731, 3.6786608731)` inner
- shell 2:
  - outer parent box: `3:13`
  - inner retained box: `4:12`
  - reference-center range: `(-5.0, 5.0)` outer, `(-4.0, 4.0)` inner
  - physical-center range: `(-3.6786608731, 3.6786608731)` outer,
    `(-1.9980860684, 1.9980860684)` inner
- shell 3:
  - outer parent box: `4:12`
  - inner retained box: `5:11`
  - reference-center range: `(-4.0, 4.0)` outer, `(-3.0, 3.0)` inner
  - physical-center range: `(-1.9980860684, 1.9980860684)` outer,
    `(-1.0225967528, 1.0225967528)` inner
- shell 4:
  - outer parent box: `5:11`
  - inner retained box: `6:10`
  - reference-center range: `(-3.0, 3.0)` outer, `(-2.0, 2.0)` inner
  - physical-center range: `(-1.0225967528, 1.0225967528)` outer,
    `(-0.4896496380, 0.4896496380)` inner
- retained direct core:
  - parent box: `ix = iy = iz = 6:10`
  - reference-center range: `(-2.0, 2.0)`
  - physical-center range: `(-0.4896496380, 0.4896496380)`

New local refinement inside the `5^3` core:

- local outer box: `6:10`
- local inner retained box / direct sub-core: `7:9`
- reference-center range: `(-2.0, 2.0)` outer, `(-1.0, 1.0)` inner
- physical-center range: `(-0.4896496380, 0.4896496380)` outer,
  `(-0.1947357748, 0.1947357748)` inner
- trusted complete-shell working box: `(2:14, 2:14, 2:14)`
- local refined-core working box: `(6:10, 6:10, 6:10)`
- combined refined-sequence working box: `(2:14, 2:14, 2:14)`

## Numerical Comparison

### Count-17

Current trusted corrected complete-shell anchor:

- fixed dim `589`
- overlap error `4.04e-12`
- `E1 = -1.9981842264017804`
- `⟨Vee⟩ = 1.24893460807307`
- projected fixed-only interaction shift `+1.33e-4`
- ground capture `99.99867457706937%`
- average first-four capture `99.66550917333363%`
- ground projected one-body shift `+4.34e-4`
- representative warmed time `1.29 s`

Hierarchical core-only refinement candidate:

- fixed dim `547`
- overlap error `1.61e-12`
- `E1 = -1.9974974043000224`
- `⟨Vee⟩ = 1.2473506228103564`
- projected fixed-only interaction shift `-1.14e-3`
- ground capture `99.95744298612602%`
- average first-four capture `99.65503769492917%`
- ground projected one-body shift `+5.51e-2`
- representative warmed time `1.25 s`

So at count `17` the refinement saves only:

- `589 -> 547`

while materially worsening:

- `ΔE1 ≈ +6.87e-4`
- `Δ⟨Vee⟩ ≈ -1.58e-3`

### Count-15

Current trusted corrected complete-shell anchor:

- fixed dim `589`
- overlap error `3.35e-12`
- `E1 = -1.9970883634888783`
- `⟨Vee⟩ = 1.247898324597102`
- projected fixed-only interaction shift `+6.69e-5`
- ground capture `99.99958917940271%`
- average first-four capture `99.9517256202187%`
- ground projected one-body shift `+6.88e-5`
- representative warmed time `0.82 s`

Hierarchical core-only refinement candidate:

- fixed dim `547`
- overlap error `1.95e-12`
- `E1 = -1.9969277398210528`
- `⟨Vee⟩ = 1.2472215233758854`
- projected fixed-only interaction shift `-1.67e-3`
- ground capture `99.9341878327217%`
- average first-four capture `99.93485338519488%`
- ground projected one-body shift `+5.96e-2`
- representative warmed time `0.76 s`

So at count `15` the same pattern remains:

- small reduction
- but clear degradation of the trusted physical regime

## Interpretation

This is a clean negative result.

What worked:

- the new core-only refinement is geometrically well defined
- it uses only the original parent-space functions assigned to the retained
  `5^3` core
- the outer complete-shell object stays fixed
- the fixed-block consumer stays unchanged
- overlap and transformed fixed-block weights remain clean

What did not work:

- the reduction from `589` to `547` is too small to matter much
- the trusted complete-shell anchor quality does not survive the refinement
- the main degradation is not orthogonality or weight pathology
- it shows up in the physical diagnostics:
  - `E1`
  - `⟨Vee⟩`
  - fixed-only interaction transfer
  - projected one-body energy shift

So this first hierarchical core-only refinement does **not** beat the current
trusted `589`-dimensional corrected complete-shell anchor enough to matter.

## Decision

The current atomic roadmap should **stop atomic hierarchy work here** rather
than pushing farther down the same core-only complete-shell line.

Practical conclusion:

- keep the current corrected complete-shell hybrid as the trusted reduced
  atomic anchor
- keep shell-plus-core as the conservative nested anchor
- do not continue with further atomic hierarchy work on the present
  core-only refinement pattern

The next move should therefore be:

- leave this negative result recorded
- and move on to the next broader phase only when ready, rather than spending
  more time on deeper atomic hierarchy compression along the same line

In particular, this result does **not** justify:

- broad shell-policy redesign in the same pass
- strong full-cube recursion
- diatomic work on an unsettled atomic hierarchy line

The atomic hierarchy branch is now informative enough to pause.
