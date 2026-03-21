# QW-PGDG Outer Radius Scan Note

## Goal

Check whether the remaining small `1s^2` scalar offset in the stabilized base
QW-PGDG `:ggt_nearest` path is partly an outer-radius effect.

The path was held fixed:

- base PGDG intermediate only
- contracted He `cc-pVTZ` `s` supplement
- `Z = 2`
- `interaction_treatment = :ggt_nearest`
- `s = 0.8`

Only `xmax` and `count` were varied.

## Scan

### `xmax = 6`

```text
count =  9: E1 = -1.965372903966, <Vee> = 1.257199273815, Δ = +0.007199273815
count = 11: E1 = -1.993654880380, <Vee> = 1.246242913273, Δ = -0.003757086727
count = 13: E1 = -1.997181740766, <Vee> = 1.247306574898, Δ = -0.002693425102
count = 15: E1 = -1.996965601277, <Vee> = 1.247214078649, Δ = -0.002785921351
count = 17: E1 = -1.996331510163, <Vee> = 1.246775715132, Δ = -0.003224284868
```

This is acceptable, but the high-count sequence is fairly flat rather than
cleanly directional.

### `xmax = 8`

```text
count =  9: E1 = -1.945636997844, <Vee> = 1.335455333250, Δ = +0.085455333250
count = 11: E1 = -1.980524992203, <Vee> = 1.244174880366, Δ = -0.005825119634
count = 13: E1 = -1.995930350378, <Vee> = 1.247960340222, Δ = -0.002039659778
count = 15: E1 = -1.996734324111, <Vee> = 1.247970805578, Δ = -0.002029194422
count = 17: E1 = -1.996008115203, <Vee> = 1.247298636121, Δ = -0.002701363879
```

The coarse counts get worse. High counts improve modestly, but the sequence is
still not especially clean.

### `xmax = 10`

```text
count =  9: E1 = -1.941063559335, <Vee> = 1.384291248220, Δ = +0.134291248220
count = 11: E1 = -1.963793871622, <Vee> = 1.259458649839, Δ = +0.009458649839
count = 13: E1 = -1.994398017744, <Vee> = 1.246493323754, Δ = -0.003506676246
count = 15: E1 = -1.998188014161, <Vee> = 1.248287568120, Δ = -0.001712431880
count = 17: E1 = -1.998635811139, <Vee> = 1.248595270324, Δ = -0.001404729676
```

Here the pattern changes in the useful way:

- `count = 9` and `11` are clearly too coarse for the larger box
- but `count = 13/15/17` becomes a more systematic high-count sequence
- `E1` moves closer to `-2`
- `⟨Vee⟩` approaches `1.25` more cleanly from below for the larger counts

## Interpretation

The outer radius does matter.

But the effect is coupled to count:

- increasing `xmax` without enough count makes the basis too coarse overall
- once count is high enough, the larger box helps the asymptotic trend look
  more physical and more systematic

So the current best practical reading is:

- `xmax = 6` is safer for small-to-medium counts
- `xmax = 10` looks better for the higher-count convergence study
  (`count >= 13`)

That suggests the remaining few-millihartree-level `⟨Vee⟩` offset at
`xmax = 6` is at least partly a finite-extent effect rather than a deeper
structural problem.
