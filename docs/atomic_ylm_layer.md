# The first angular `(l,m)` layer

This note describes the first explicit atomic angular layer in GaussletBases.

It sits on top of the radial substrate that is already in the package:

- a radial basis
- an explicit radial quadrature grid
- radial one-body operators
- the primitive/contraction story underneath those radial operators

The point of this layer is simple:

**turn the radial one-electron hydrogen problem into an explicit atomic `(l,m)` problem**

That is the natural next step before later work on helium, IDA-style angular coupling, or broader atomic structure.

## 1. What the radial layer already gives you

The radial line already provides the reduced-radial matrices

- overlap
- kinetic
- nuclear attraction
- centrifugal

on a chosen radial basis and quadrature grid.

For example:

```julia
rb = build_basis(spec)
grid = radial_quadrature(rb)
radial_ops = atomic_operators(rb, grid; Z = Z, lmax = 2)
```

At that point, the radial one-body problem is already well defined.

## 2. What the angular layer adds

The new angular layer adds:

- explicit channels labeled by `(l,m)`
- an ordered list of those channels up to some `lmax`
- an atomic one-body Hamiltonian assembled from the radial blocks

The narrow public objects are:

- `YlmChannel`
- `YlmChannelSet`
- `AtomicOneBodyOperators`

with the main construction calls:

```julia
channels = ylm_channels(2)
atom = atomic_one_body_operators(radial_ops, channels)
```

or equivalently:

```julia
atom = atomic_one_body_operators(radial_ops; lmax = 2)
```

## 3. Why this is simple for hydrogen

For a central one-electron Hamiltonian such as hydrogen,

```math
H = -\frac12 \nabla^2 - \frac{Z}{r},
```

the angular dependence is carried by the spherical harmonics `Y_{lm}`.

At this stage:

- the Hamiltonian is diagonal in the angular channels
- the dependence on `l` enters through the centrifugal term
- the dependence on `m` is only through degeneracy

So for one channel `(l,m)`, the one-body block is simply

```math
T + V_{\mathrm{nuc}} + C_l.
```

That is why the first Ylm layer can stay narrow and explicit.

## 4. How the radial and angular layers combine

The radial layer and the angular layer play different roles.

The radial layer carries:

- the radial basis functions
- the quadrature
- the radial matrix elements

The angular layer carries:

- the channel list
- the block structure
- the repeated `m` degeneracy

So the atomic one-body object is not a replacement for the radial layer. It is a clean way of organizing the radial data into an atomic basis labeled by `(l,m)`.

## 5. Why this comes before He / IDA

This is the right step to take before the later interacting atomic path.

The later He / IDA story should sit on top of:

1. the radial one-body substrate
2. the explicit `(l,m)` channel structure
3. the radial multipole operators
4. later Gaunt or related angular-coupling factors

In other words, the current one-electron Ylm layer is the place to settle:

- channel ordering
- block conventions
- how radial and angular data fit together

before introducing electron-electron structure.

## 6. What this layer does not do yet

This first angular slice does **not** yet provide:

- a two-electron Hamiltonian
- Gaunt-coupled electron-electron assembly
- a full atomic workflow beyond the one-electron block structure

It is intentionally only the first explicit one-electron atomic layer.

## 7. Example

The repository includes:

- `examples/15_atomic_hydrogen_ylm.jl`

That example:

- builds a recommended hydrogen radial basis
- builds the radial one-body operator bundle
- constructs the `(l,m)` channel list up to `lmax`
- assembles the block-diagonal one-electron atomic Hamiltonian
- diagonalizes it and reports low-lying energies

That is the cleanest way to see how the present radial substrate and the new angular layer fit together.
