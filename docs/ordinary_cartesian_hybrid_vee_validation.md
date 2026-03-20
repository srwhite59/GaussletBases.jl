> **Status:** supporting development note. For the current ordinary branch,
> read `docs/current_ordinary_branch.md` first.

# Ordinary Cartesian IDA: Hybrid/Core-Gaussian `Vee` Validation

This note records the next narrow step after the pure ordinary Cartesian
`1s^2` interaction check.

## 1. Why this is the next pass

The pure ordinary check made the Cartesian IDA interaction layer concrete on a
simple orthonormal basis and gave one meaningful scalar target.

The next practical ordinary-branch question is different:

- not whether the pure ordinary branch is already enough for He
- but whether the friendlier hybrid/core-supported regime improves the same
  scalar interaction check in the direction the practical ordinary branch is
  supposed to work

That comes before any real He solver.

## 2. What this pass is testing

The target remains the same doubly occupied noninteracting `1s` reference
state:

1. build the one-body ordinary Cartesian Hamiltonian
2. take its lowest orbital
3. occupy that orbital once with spin up and once with spin down
4. evaluate the resulting density-density / IDA interaction expectation

But now do it on the current hybrid one-dimensional basis:

- mild mapped ordinary backbone
- explicit centered core Gaussians
- final overlap cleanup / localization already included in the basis object

## 3. Scope of the current hybrid/Gaussian route

This pass uses the current combined hybrid basis machinery that is already in
the package.

That means the static `Vee` matrix is built on the final hybrid working basis,
so the resulting matrix already includes:

- backbone-backbone contributions
- backbone-core Gaussian contributions
- core-core Gaussian contributions

There is **not yet** a separate residual-Gaussian transfer or DTA-style
treatment in this pass. That remains a later follow-on if the practical
ordinary branch needs it.

## 4. What counts as success here

The first meaningful question is simple:

**does the hybrid/core-supported ordinary basis improve the `1s^2` scalar
interaction check relative to the corresponding pure mapped localized basis?**

The exact hydrogenic reference is still

```text
⟨Vee⟩ = (5 / 8) Z,
```

so for `Z = 2` the scalar target remains `1.25 Eh`.

## 5. What this pass is not doing

This is still narrow.

It is not:

- a full He solve
- a residual-Gaussian DTA layer
- a broad hybrid basis-library system

It is only the first hybrid/core-Gaussian validation of the ordinary
Cartesian `Vee` layer on the same `1s^2` scalar used in the pure ordinary
check.
