# Ordinary Mapped Hydrogen: Coulomb Calibration and PGDG Direction

This note records the current decision point for the ordinary mapped branch.

## 1. What the repo is using now

The current repo uses the internal `ShortGaucoulomb`-derived construction in
[`src/ordinary_coulomb.jl`](../src/ordinary_coulomb.jl) as its deterministic
Gaussian expansion of bare `1/r`.

With the current default high-accuracy choice, that means the sinh-mapping
parameters

- `del = 1.0`
- `s = 0.16`
- `c = 0.01`
- `maxu = 135.0`

which gives about `135` Gaussian terms.

This is a clean, deterministic source of separable Gaussianized Coulomb
factors, and it is not the present bottleneck.

## 2. How that differs from the historical plain-Coulomb default

The historical plain-Coulomb branch in `Gaucoulomb.newgaucoulomb(doacc=true)`
uses the same sinh-mapping idea and the same high-accuracy parameters

- `del = 1.0`
- `s = 0.16`
- `c = 0.01`

but stops earlier, at

- `maxu = 115.0`

so the historical plain-Coulomb default is about `115` terms rather than the
current repo default of about `135`.

This means the present repo is using a slightly tighter short-range plain-
Coulomb choice than the older historical default, not a different method.

## 3. Why the 29-term `Gaucoulomb` fits are not the right comparison

The very short `29`-term branches in `Gaucoulomb.jl` are not alternative bare
`1/r` defaults. They are fits for the Takeshi/Yanai range-separated form with a
finite screening parameter.

So they are not the right direct comparison for the current ordinary hydrogen
path, which is still using bare Coulomb attraction.

The meaningful plain-Coulomb comparison is therefore:

- historical sinh-mapped plain `1/r`: about `115` terms
- current repo sinh-mapped plain `1/r`: about `135` terms

## 4. Why the mapped-integral cost now matters more than the expansion length

In the current mapped ordinary branch, the expensive step is not the final
dense one-electron diagonalization, and it is not mainly the difference between
`115` and `135` Gaussian terms.

The expensive step is building the mapped one-dimensional operator factors,
especially when each Gaussianized factor is evaluated through numerical
integration on explicitly distorted primitives.

So the main performance question is no longer:

- how many Coulomb Gaussians should the bare `1/r` fit use?

It is now:

- should the mapped ordinary one-body path continue to rely on numerical
  mapped-primitive integration?

## 5. Where COMX and localization would enter historically

Historically, the PGDG line is not just:

- take the same final basis
- swap numerical integrals for analytic ones

The fuller historical path is closer to:

1. build a plain Gaussian primitive layer in physical space
2. form analytic one-body matrices on that primitive layer
3. orthogonalize and clean up the primitive-backed working space
4. apply COMX/localization/projected diagonalization steps
5. keep the final working basis after those basis-construction steps

So COMX/localization belongs to the basis-construction side of the historical
PGDG story, not just to the operator-evaluation side.

## 6. What the present prototype is and is not

The present pass stops earlier on purpose.

It is only a:

- primitive/contraction-level analytic prototype

and it is **not yet**:

- the full historical PGDG basis-construction path

In particular, this first prototype does **not** recreate the later
COMX/localization/projected-diagonalization stages. Instead, it keeps the
current mapped ordinary working basis and asks a narrower question:

- if the primitive layer is replaced by a locally matched plain-Gaussian proxy,
  do the one-body matrices and hydrogen energies stay close to the current
  mapped numerical reference?

That makes the comparison honest:

- current numerical mapped path: same working basis, numerical mapped primitive
  integrals
- present PGDG-style prototype: same working basis, analytic plain-Gaussian
  primitive proxy

So this is a comparison at the primitive/operator layer, not yet at the full
historical basis-construction layer.

## 7. Why the next prototype should be PGDG-style

The old ordinary branch suggests a cleaner direction:

- replace the explicitly distorted primitive layer by a locally matched plain
  Gaussian layer in physical space
- reuse analytic one-dimensional Gaussian matrix elements
- keep the current contraction layer from primitives to working basis

That is the point of the present `MappedPGDGPrototype1D` pass.

It is not a full port of the old PGDG workflows. It is only a one-body
decision prototype meant to answer one question:

**does the mapped ordinary hydrogen path look better when the primitive layer is
treated analytically rather than through repeated numerical mapped integrals?**

For now, the current numerical mapped ordinary path remains the reference and
validation route. The PGDG-style path is the candidate fast route.
