# Cr2 Staged Subspace Filter Audit (870498b54)

## Plain-Language Result

Staged subspace filtering is the first protected-original broad-injection
geometry that looks numerically healthy after the aspect-aware complete-shell
fix.

The successful ordering is:

1. start from the full non-narrow broad original subspace `W`;
2. filter by collective representability in `M = [G, R_compact]`;
3. optionally localize and apply shape filters;
4. diagonalize fake-RDM in the surviving subspace;
5. keep fake-RDM eigenvectors above an occupancy threshold.

The best tradeoff in this probe was:

```text
s_cut       = 0.95
shape       = none
occ_cut     = 0.003
broad dim   = 87
Z dim       = 117
fake trace  = 93.3726973285
B_min       = 0.9934658245
B < 0.99    = 0
```

This improves on both the old strict protected-original geometry and the failed
scalar relaxed geometry:

```text
old nonlocalized strict: broad 94,  B_min 0.990572
localized strict:        broad 80,  B_min 0.972017
failed scalar relaxed:   broad 100, B_min 0.869302
staged best:             broad 87,  B_min 0.993466
```

The important caveat is that the best no-shape-filter variant is not an
atom-local-only policy. After fake-RDM diagonalization, the retained final
directions reassemble into mostly balanced bond-axis and transverse
combinations. That is acceptable only if injection is understood as replacement
inside a well-represented main space, not as residual-MWG broad channels.

The conservative shape-filter variants keep the final directions more local,
but lose a large amount of fake-RDM trace. For example, `s_cut=0.95`,
conservative shape filtering, and `occ_cut=0.03` gives a clean `B_min` of about
`0.992856`, but retains only `82.9517` fake-RDM trace rather than `93.3727`.

## Decision Interpretation

- Subspace-first representability filtering is viable.
- Scalar per-mode fake/representability cuts are not viable.
- Combined fake/representability score prefixes were not viable.
- A production design should treat broad retained directions as injected
  replacement directions only, never as residual Gaussian/MWG channels.
- If source work is approved, the design should expose `s_cut` and `occ_cut`
  only as private/internal measurement parameters at first.
- No Cr2 HF, artifact writing, public driver input, or provenance/schema work
  is justified by this report alone.

## Method

The probe rebuilt the current aspect-shell Cr2 `ns=7`, `lmax=2`, `cc-pV5Z`
protected-original setup from the old artifact recipe. It used the compact
ordered-MGS residual selector to build:

```text
G
R_compact
M = [G, R_compact]
protected originals = 30
W = 108 broad original directions after protected-original projection and
    Gaussian Gram cleanup
```

For each representability cutoff, it diagonalized:

```text
B = M' S W
rho_rep = B'B
```

and retained singular directions with:

```text
s = sqrt(eig(rho_rep)) >= s_cut
```

It then optionally localized with the `z` operator and applied either no shape
filter or the conservative shape filter from the localization report. Finally,
it diagonalized fake-RDM in the surviving subspace and retained fake-RDM
eigenvectors above `occ_cut`.

## Best Variant

```text
variant:                       s0.95_shape_none_occ_0.003
Z dimension:                   117
broad retained dimension:       87
fake-RDM trace retained:        93.37269732851348
fake-RDM trace dropped:          0.05661708046312697
B min / median / max:            0.9934658245058705 / 0.9999186883626252 / 1.0
B counts <0.999/<0.99/<0.98:    26 / 0 / 0
Z' S M Qperp max:                8.729e-16
sampled Qperp identity max:      9.992e-16
F' S F - I block max:            1.164e-9
protected span singular min:     1.0
```

Retained channel summary for this variant:

```text
dxx:6, dxy:10, dxz:6, dyy:4, dyz:6,
px:13, py:13, pz:2, s:27
```

No dropped direction above the probe's "important dropped direction" reporting
threshold was recorded for this best variant.

## Files

- `stage_table.tsv`: stage-by-stage dimensions, fake-RDM trace, singular
  values, and shape counts for the tested ladder.
- `final_geometry.tsv`: final geometry diagnostics for each completed variant.
- `dropped_important_directions.tsv`: bounded details of important directions
  dropped at each stage. The full retained-direction table was intentionally
  not committed because it is about 1900 rows; the raw probe output remains in
  `/Users/srw/dmrgtmp/cr2_staged_subspace_filter_870498b54/`.

## Validation

- Package load: `0.455547708 s`.
- Probe outer elapsed: `70.178911167 s`.
- `git diff --check`: clean in the doer pass.
- No tracked source edits were made by the measurement pass.
