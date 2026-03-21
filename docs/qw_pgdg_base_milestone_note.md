# Base QW-PGDG Milestone Note

## Status

The current base QW-PGDG path is now the starting point for the next phase.

The active PGDG side is quadrature-free, the final COMX step is in the right
place for the basis-realization branch, the legacy He `s` supplement is now
contracted, and the `:ggt_nearest` QW-PGDG validation path is running
robustly.

The immediate goal of this note is a larger-count convergence/scaling check
before any nesting work starts.

## Fixed path for this milestone

This check keeps the path fixed:

- base PGDG intermediate only, `refinement_levels = 0`
- He `cc-pVTZ` contracted `s` supplement
- `Z = 2`
- `ordinary_cartesian_qiu_white_operators(...; interaction_treatment = :ggt_nearest)`
- no refinement hierarchy rollout
- no sliced-basis logic
- no nesting logic

## Convergence sequence

```text
count = 9:
  contracted supplement functions = 3
  kept residual modes             = 3
  E1                              = -1.9653729039662613
  <Vee>                           = 1.2571992738145692
  Δ(<Vee> - 1.25)                 = +0.007199273814569196
  constructor time                = 3.600558206 s

count = 11:
  contracted supplement functions = 3
  kept residual modes             = 2
  E1                              = -1.9936548803801342
  <Vee>                           = 1.2462429132727513
  Δ(<Vee> - 1.25)                 = -0.0037570867272487263
  constructor time                = 0.869931379 s

count = 13:
  contracted supplement functions = 3
  kept residual modes             = 2
  E1                              = -1.9971817407657657
  <Vee>                           = 1.247306574897531
  Δ(<Vee> - 1.25)                 = -0.0026934251024690603
  constructor time                = 2.08453247 s

count = 15:
  contracted supplement functions = 3
  kept residual modes             = 2
  E1                              = -1.9969656012766186
  <Vee>                           = 1.2472140786487742
  Δ(<Vee> - 1.25)                 = -0.002785921351225795
  constructor time                = 7.794032719 s

count = 17:
  contracted supplement functions = 3
  kept residual modes             = 2
  E1                              = -1.9963315101634858
  <Vee>                           = 1.2467757151316037
  Δ(<Vee> - 1.25)                 = -0.003224284868396321
  constructor time                = 9.513510435 s
```

## Interpretation

This looks like a solid base milestone.

- `count = 9` is the coarse outlier
- from `count = 11` onward, `⟨Vee⟩` is already clustered within about
  `0.0038` of the `1.25` target
- the `11/13/15/17` sequence looks like clean finite-resolution convergence,
  not a blocked or unstable path
- runtime is still reasonable for this light He validation problem
- the residual space is stable at `2` kept modes for the larger counts

So the current base QW-PGDG route looks ready to serve as the clean starting
point for the next nesting phase.

## Proposed milestone commit message

```text
Stabilize base QW-PGDG path before nesting
```
