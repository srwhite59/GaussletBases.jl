# Sectorized Angular Preparation for Atomic IDA

This note describes the next angular step after introducing the sparse/block
`GauntTable` backend.

## Why the new Gaunt backend changes the situation

The earlier atomic IDA layer already had the right mathematical ingredients:

- radial multipole operators
- explicit `(l,m)` channels
- Gaunt coefficients
- angular kernels built from those Gaunt coefficients

But it still paid the main angular cost too early. After building the Gaunt
data, it immediately expanded back out to dense full-channel objects:

- `gaunt_tensor`
- `angular_kernel`

That is acceptable for very small `lmax`, but it hides the main conserved
structure and scales badly.

The sparse/block `GauntTable` backend changes this because the angular data are
no longer ad hoc. Once the Gaunt coefficients are organized by legal
`(L,l_1,l_2)` blocks, the next natural step is to organize channel pairs by the
conserved angular quantum number that matters in the IDA kernels.

## Conserved pair sectors

For the complex `Y_{lm}` basis used here, the angular kernel couples channel
pairs

```math
(\alpha,\beta) \longleftrightarrow (\alpha',\beta')
```

only when the magnetic quantum numbers satisfy

```math
m_\alpha + m_\beta = m_{\alpha'} + m_{\beta'}.
```

So the full pair space naturally breaks into independent sectors labeled by the
conserved pair sum

```math
m_{\mathrm{sum}} = m_\alpha + m_\beta.
```

This is the same basic idea that appeared in the older
`atombasisYlmopt.jl` line: group pair indices by a conserved `m` quantity
before doing the expensive angular work.

## What is built now

The new internal preparation layer does three things:

1. enumerates channel pairs in the current repo channel ordering
2. groups those pairs by conserved `m` sum
3. builds one angular kernel matrix per sector and per multipole `L`

So the atomic IDA layer now keeps a sectorized angular representation
internally, rather than storing the full dense four-index kernel eagerly.

## What stays public

The public atomic story does not need to change for this pass.

The following public objects and accessors stay conceptually the same:

- `AtomicIDAOperators`
- `gaunt_tensor`
- `gaunt_coefficient`
- `angular_kernel`

The difference is internal:

- Gaunt coefficients come from the sparse/block `GauntTable`
- angular kernels are prepared sector by sector
- dense full-channel tensors are reconstructed only when explicitly requested

That keeps the present user-facing interface stable while giving later
application code a cleaner and more scalable angular foundation.

## Why this is the right middle step

This pass is intentionally narrower than porting the older `VeeYlm` or
`DiagYlm` application layers wholesale.

The goal here is only:

- exploit the exact angular selection rules already present
- make conserved pair sectors explicit
- prepare a better internal structure for later atomic IDA work

It is **not** yet:

- a Hartree/Fock application layer
- a Coulomb-application engine
- a packed external storage format

## Relation to later steps

Once the angular data are sectorized, later code can work with those sectors
directly instead of repeatedly materializing dense `n_{\mathrm{chan}}^4`
objects.

That should make the next atomic steps much cleaner:

- improved two-electron assembly paths
- later Hartree or exchange application code
- larger `lmax` studies without immediately paying the full dense angular cost
