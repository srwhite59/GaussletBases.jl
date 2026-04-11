# Ordinary Hybrid / Residual-Gaussian Clarification

This note records an important correction to the current ordinary-branch
interpretation.

It is not a solver note. It is a formulation note.

The point is to make clear:

- what the current code actually does
- why that is not yet the paper-faithful residual-Gaussian (RG) formulation
- where COMX/localization does and does not belong
- what the next correction should be before any ordinary He-style solve

## 1. Main correction

The current ordinary hybrid / residual-Gaussian implementation is **not yet**
the Qiu–White hybrid formulation.

The key correction is:

- the full three-dimensional Cartesian gausslet basis should be built first
- the added Gaussian orbitals should then be treated as genuine **3D**
  orbitals
- the residual Gaussians should be defined by orthogonalizing those 3D
  Gaussian orbitals against the **full 3D gausslet basis**
- the final basis for the residual-Gaussian interaction story should be:
  - the original 3D gausslet basis
  - plus the orthogonalized residual Gaussians
- there should be **no COMX/localization or other further rotation** on that
  combined gausslet-plus-RG set

That means the present 1D residual-Gaussian construction is only a provisional
simplification. It is not the right target for the paper-faithful hybrid/RG
route.

## 2. What the current code actually does

At present, the ordinary hybrid branch does something different.

### 2.1 Legacy combined hybrid basis construction

The legacy/internal hybrid basis constructor is:

- [`hybrid_mapped_ordinary_basis(...)`](../src/ordinary_hybrid.jl)

In the current implementation, it:

1. chooses a 1D ordinary backbone layer from the current mapped ordinary
   backend
2. appends the added core Gaussians as explicit extra 1D primitives
3. builds a combined seed basis from those pieces
4. applies `_cleanup_comx_transform(...)` to that combined set

The relevant code is in
[src/ordinary_hybrid.jl](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/src/ordinary_hybrid.jl#L109).

So the current final `HybridMappedOrdinaryBasis1D` is:

- not the original backbone gausslets
- not a gausslet-plus-RG basis
- but a COMX-cleaned/localized rotation of backbone plus raw added Gaussians

### 2.2 Residual-Gaussian interaction path

The current residual-Gaussian interaction work is in
[src/ordinary_cartesian_ida.jl](/Users/srw/Library/CloudStorage/Dropbox/codexhome/repositories/GaussletBases/src/ordinary_cartesian_ida.jl#L446).

That code currently:

1. reconstructs a 1D backbone block
2. orthogonalizes the added 1D Gaussian channel against that 1D backbone
3. forms 1D residual-Gaussian directions
4. builds either:
   - a nearest-center residual interaction ansatz
   - or an MWG-style residual interaction ansatz
5. transfers that interaction back into the final COMX-cleaned hybrid basis

So the current residual interaction path is:

- defined in **1D**
- then lifted separably into the 3D Cartesian IDA assembly
- then transferred into a final basis that has already been rotated by the
  hybrid COMX cleanup

That is not the paper-faithful residual-Gaussian construction.

## 3. Why this is a real mistake, not just a stylistic difference

This is not only a matter of presentation. It changes the scientific meaning
of the approximation.

The Qiu–White residual-Gaussian story is built around the fact that:

- the final working gausslet basis is already fixed
- the added Gaussian channel is orthogonalized against that full gausslet
  space
- the residual Gaussians then live in the orthogonal complement
- diagonal interaction approximations for those residual Gaussians are
  justified because their occupancies are very small

If the code instead:

- builds a 1D residual channel first
- or rotates the combined gausslet-plus-Gaussian set afterward with COMX
- or transfers residual interactions into a different final basis

then the residual channel is no longer the same object that the paper is
talking about.

So one can no longer interpret:

- nearest-center/GGT behavior
- MWG behavior
- residual occupancies
- or their effect on `⟨Vee⟩`

as direct evidence about the paper-style residual-Gaussian approximation.

## 4. What the recent occupancy check still tells us

Even with that caveat, the recent occupancy check was still informative.

In the current calibrated friendly regime (`count = 11`, `s = 0.6`) with the
legacy He `s` supplements:

- `cc-pVTZ`: residual-channel weight in the lowest orbital was about
  `4.600159e-4`
- `cc-pVQZ`: residual-channel weight in the lowest orbital was about
  `4.600030e-4`

So the residual-channel occupancy is indeed small, around `10^-4` to
`10^-3`-level once doubled for the `1s^2` reference state.

That supports the expected physical story:

- the residual channel should be low occupancy
- its interaction treatment should not strongly control the final answer

The fact that the present residual interaction approximations still change the
`1s^2` scalar noticeably is therefore another sign that the current
implementation path is not the right formulation yet.

## 5. Where COMX *does* belong later

COMX/localization is not being rejected in general.

The important distinction is stage and purpose.

For later work on:

- nested gausslets
- contraction hierarchies
- leaf/local contraction frameworks
- or later 1D compression/localization pipelines

COMX can be entirely appropriate.

But that is a different stage of development.

For the current paper-faithful hybrid/RG interaction story:

- COMX on the combined gausslet-plus-RG set is a mistake
- using a COMX-cleaned hybrid basis as the final basis for the RG/IDA path is
  a mistake

So COMX is a later nested/contraction tool, not part of the present
Qiu–White-style hybrid/RG formulation.

## 6. Consequences for the current MWG pass

This changes how the current MWG work should be interpreted.

The present MWG implementation is still useful as a coding experiment, but it
should be labeled honestly:

- it is **not yet** a faithful implementation of the paper's MWG route
- it is an MWG-like interaction model built on top of the current 1D residual
  construction and transferred into a COMX-cleaned hybrid basis

So the present nearest-center vs MWG comparison should not be treated as the
final scientific decision on GGT vs MWG for the ordinary branch.

## 7. Correct next implementation target

The next correct implementation target should be:

1. build the full 3D Cartesian gausslet product basis
2. build the actual 3D Gaussian supplement orbitals
3. orthogonalize those 3D Gaussian orbitals against the full 3D gausslet
   basis
4. keep the final basis explicitly as:
   - gausslets
   - plus orthogonalized 3D residual Gaussians
5. define nearest/GGT and MWG directly for that 3D residual channel
6. evaluate the same `1s^2` scalar and residual occupancy checks there

Only after that can the repo make a paper-faithful scientific comparison
between:

- combined hybrid treatment
- nearest/GGT residual treatment
- MWG residual treatment

## 8. What should *not* be done next

Before that formulation fix:

- do **not** open an ordinary He solver
- do **not** treat the current MWG numbers as final evidence
- do **not** tune residual interaction approximations further inside the
  current COMX-based hybrid path

The correct next milestone is to fix the basis/formulation mismatch first.

## 9. Short version

The short version is:

- the current hybrid/RG implementation is not yet Qiu–White-faithful
- the main reason is that residual Gaussians are currently constructed in 1D
  and then transferred into a COMX-cleaned hybrid basis
- the correct formulation is full 3D gausslets first, then 3D residual
  Gaussians, with no COMX on the combined set
- COMX belongs later for nested/contraction work, not here

That clarification should be treated as the current authoritative reading of
the ordinary hybrid/RG branch until the code is corrected.
