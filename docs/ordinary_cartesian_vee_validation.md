> **Status:** supporting development note. For the current ordinary branch,
> read `docs/current_ordinary_branch.md` first.

# Ordinary Cartesian IDA: First `Vee` Validation

This note records the next ordinary-branch milestone after the static
Cartesian IDA object itself: make the ordinary Cartesian `Vee` layer explicit
enough that it can be checked on one simple two-electron scalar.

## 1. Why this comes before a He solver

The package already has the static ordinary Cartesian ingredients:

- a one-body Hamiltonian on the Cartesian product basis
- a separable Coulomb-expansion assembly
- a dense two-index density-density / IDA interaction matrix

The right next step is not yet a full He solver. The interaction layer should
first be validated on the simplest meaningful expectation value.

## 2. Validation target

The first scalar target is the doubly occupied noninteracting `1s` state:

1. solve the one-body He-like problem on the pure ordinary Cartesian basis
2. take the lowest orbital
3. occupy that spatial orbital once with spin up and once with spin down
4. evaluate the resulting IDA interaction expectation value

This is narrow but physically meaningful.

## 3. Reference value

For an exact hydrogenic `1s` orbital,

```text
J(1s, 1s) = (5 / 8) Z.
```

So in the He-like `Z = 2` case the scalar target is

```text
⟨Vee⟩ = 1.25 Eh.
```

That gives the present ordinary Cartesian IDA layer a concrete number to aim
at before any broader He workflow is opened.

## 4. Scope of this pass

This pass is intentionally pure ordinary gausslets first:

- no hybrid/core Gaussian augmentation yet
- no residual Gaussian IDA treatment yet
- no He solver yet

The point is only to make the static ordinary Cartesian `Vee` object explicit
and check it on the simplest `1s^2` expectation value.

## 5. Expected follow-on after this

Once the pure ordinary `Vee` check is in place, the next natural questions are:

- hybrid/core-Gaussian augmentation
- residual Gaussian IDA treatment
- and only then whether the ordinary branch is ready for a first real He solve
