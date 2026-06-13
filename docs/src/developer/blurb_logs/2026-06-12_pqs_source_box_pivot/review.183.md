Accepted pass 183.

The matching old WL gausslet-only 419 baseline was not already present as a
complete scalar result, so doer computed it locally in an ignored probe:

```text
WL source = old_nested_fixed_block_wl_gausslet_only
d = 0.3
q = n_s = 5
core = 125
shell retained per layer = (98, 98, 98)
final dimension = 419
AHGBS residuals excluded = true
```

Manager reran:

```text
julia --project=. tmp/work/wl_he_419_gausslet_only_probe.jl
```

and reproduced:

```text
WL H1 = -1.991344469963435
WL H1-orbital self-Coulomb = 1.2420473874925473
WL RHF total = -2.85080350301779
WL RHF one-electron = -3.871408908674227
WL RHF electron-electron = 1.0206054056564373
```

Diagnostic deltas on matched 419 gausslet-only spaces:

```text
PQS RHF total - WL RHF total = +0.0014163805874733981 Ha
PQS H1 - WL H1 = +0.004662494788541416 Ha
PQS H1-J self-Coulomb - WL H1 self-Coulomb = -0.015884787180628912 Ha
PQS RHF two-body - WL RHF two-body = -0.002427934232285267 Ha
```

So, at this q=5/419 diagnostic point, current PQS is not better than the
matching WL gausslet-only baseline. The dominant sign points to the one-body
side: PQS has a less favorable H1 by about 4.66 mHa, partly offset by a lower
two-body contribution.

The supplemented WL 447 Fig.8-style result remains a separate baseline and
must not be treated as this gausslet-only comparison.

No tracked source/test changes were made.

-- repo-manager@macmini
