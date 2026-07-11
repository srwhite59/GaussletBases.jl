# Numerical-Complete Residual Gaussian Basis

Status: approved internal opt-in source and validation authority under
`HP-RG-NUMCOMP-FN-01` and `HP-RG-NUMCOMP-TEST-01`; implementation pending.

This page is the canonical contract for retaining the numerical Gaussian
complement of an explicit supplement while preserving the terminal gausslet
basis. The registry owns permission and lifecycle. This page owns the
numerical policy, additive-reference composition, validation, and exclusions.

## Purpose

For the first screened additive-reference comparison, use the simpler
augmented basis

```text
M = [G, R_num]
```

where `G` is the existing orthonormal terminal basis and `R_num` is the
owner-local numerical complement of the complete selected Gaussian supplement
`A`. This path appends true residual directions. It does not replace any
direction in `G`, inject original or occupied orbitals, or require a final
localization.

The ordinary production residual policy remains the `1e-6` occupation cutoff
owned by `HP-RG-CUTOFF-FN-02`. The `1e-10` value below is a separate explicit
internal numerical-complete policy, not a new default or a conditioning repair.

## Evidence And Decision

The accepted Cr2 fixed-density decomposition isolated a large protected
replacement penalty while the compact residual extension itself was small:

```text
M = G + R_compact                         +0.252 mHa H1 error
mandatory occupied replacement            +1.082 mHa
compact-source original replacement      +119.751 mHa
optional broad repair                      +28.510 mHa
F -> L localization change             effectively zero
```

Bounded H2 and padded Be2 measurements found low broad-residual RHF
occupations, `4.78e-5 e` and `2.07e-3 e`, respectively. Same-span owner-blocked
Boys rotations changed MWG/RHF endpoints by at most `0.00210 mHa` and
`0.03492 mHa`. These measurements support retaining the numerical residual
span without a compactness prefilter. They do not establish a correlated-method
claim or justify Gaussian-array enrichment.

## Numerical Contract

Inputs are:

```text
G       orthonormal terminal basis
A       explicit Gaussian supplement, partitioned by physical owner
X       <G|A>
S_AA    <A|A>
```

For each owner block `a`, form the residual metric

```text
M_a = S_AA[a,a] - X[:,a]' * X[:,a]
```

Validate the metric with the existing negative-eigenvalue tolerance. Retain
exactly the eigenmodes satisfying

```text
lambda > eta_num
eta_num = 1e-10.
```

The inequality is strict. Do not floor, clamp upward, or otherwise manufacture
rank. A materially negative eigenvalue or an unstable near-threshold direction
is a construction failure.

Use the existing owner-local normalization, concatenate the owner blocks, and
perform the existing single global inverse-square-root/Lowdin merge. The
result must satisfy

```text
R_num' S G       = 0
R_num' S R_num   = I
M                = [G, R_num].
```

The initial source path keeps the existing unlocalized final orientation. It
does not run Boys localization or any other same-span rotation.

## Existing Builder Is The Algorithm Owner

`build_residual_gaussian_basis(...)` already implements this algorithm when
called with:

```text
residual_occupation_cutoff = 1e-10
residual_injection_cutoff  = 0.0
residual_compactness       = nothing
```

Implementation must reuse that builder and the existing
`CartesianResidualGaussianBasis`. It must not duplicate the owner-local
eigensolve, define another residual object, or create a second merge path. A
narrow internal caller may name the numerical-complete policy, but all
construction must continue through `build_residual_gaussian_basis(...)`.

There is no compactness, width, zeta, center-distance, or tail prefilter in
this policy. The complete supplement selected by the caller is presented to
the owner-local residual metric.

## Packet Occupied-Space Contract

Any atomic reference packet used by a consumer must match the explicit
supplement basis and its occupied orbital span must lie in the selected
supplement span. Packet occupied orbitals do not select, append, replace, or
rotate `R_num`.

After constructing `M`, represent every original packet occupied block by its
exact cross overlap into `[G,R_num]`. Validate each packet separately for:

```text
occupied singular-value capture
occupied Gram/orthogonality error
electron trace
owner, center, basis, order, and overlap identity
```

Failure at the existing strict `1e-10` capture tolerance means that the
selected supplement or retained numerical rank is insufficient. Stop. Do not
append the occupied directions, weaken capture, reinterpret fitted
density/potential Gaussians as orbitals, or invoke occupied-first/protected
injection.

## Additive-Reference Consumer

The approved internal consumer builds one native unlocalized augmented member
in `M = [G,R_num]` order and reuses existing facilities:

- exact augmented kinetic, by-center unit-nuclear, position, and moment
  transforms;
- placed fitted-potential `GG/GA/AA` reference-Hartree blocks;
- the existing augmented-operator transform for `J0_M`;
- packet occupied-block representation for additive `P0_M/q0_M`;
- density-fit self and cross energies for additive `E0`;
- the existing residual-containing MWG/IDA interaction in the final merged
  residual orientation;
- the existing in-memory `ScreenedHartreeCorrection` algebra.

The correction remains separate from the member:

```text
Delta_J0 = J0_M - Diagonal(Vee_M * q0_M)
C        = 0.5*q0_M' * Vee_M * q0_M - 0.5*E0.
```

Screening off/on uses the same `H1_M` and `Vee_M`; only the separately returned
`Delta_J0/C` is applied. The one-body correction is direct electron-electron
accounting, not a replacement kinetic or nuclear operator.

The implementation may reuse private input, packet-placement, raw-field, and
energy helpers near the protected ladder owner. It must not construct
protected geometry, perform `F -> L` localization, write a protected artifact,
or alter protected-localized ladder and EGOI facilities.

## Approved Source Surfaces

`HP-RG-NUMCOMP-FN-01` approves only:

- `src/cartesian_residual_gaussians/residual_basis.jl`, only if a narrow policy
  assertion or diagnostic is needed around the existing builder;
- `src/cartesian_residual_gaussians/augmented_operators.jl` for packet
  occupied-block representation in native `[G,R_num]` order;
- `src/cartesian_protected_ladder_bundle.jl` for narrow private in-memory
  composition reusing existing inputs, placements, fields, energies, augmented
  operators, and MWG;
- `src/cartesian_reference_density/atomic_hf_reference_packets.jl` and
  `src/cartesian_reference_density/screened_hartree_correction.jl` only for
  narrow reuse when no contract or persistent result shape changes.

Target at most `180` added source lines. If the path requires another builder,
module, persistent struct, artifact field, or metadata cloud, stop and request
new authority.

## Validation

`HP-RG-NUMCOMP-TEST-01` approves:

1. A compact synthetic contract in `test/misc/runtests.jl` that distinguishes
   a retained residual eigenvalue above `1e-10` from a true numerical null and
   verifies malformed/negative metric failure.
2. Narrow explicit numerical-complete coverage in the existing
   `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`; ordinary
   default assertions must remain unchanged.
3. Ignored physically padded H2/Be2 construction and endpoint probes under
   `tmp/work`, with durable output only under `/Users/srw/dmrgtmp`.
4. Only after the small-system gates pass, one ignored Cr2 fixed imported-density
   comparison against the accepted protected `+28.510 mHa` result.

Every endpoint must inspect terminal due diligence. Report:

```text
candidate and retained counts by owner
owner residual-metric spectra and distance from eta_num
final merge spectrum/condition
max |G' S R_num|
max |R_num' S R_num - I|
packet occupied capture and trace
P_RR spectrum and broad residual occupation
P_GR coupling
low H1 and Fock modes with owner/channel makeup
MWG finiteness and symmetry
screened and unscreened endpoint deltas where run
```

Stop before Cr2 if H2 or Be2 develops a bad low mode, failed capture, unstable
metric, or material MWG-sensitive endpoint change. A Cr2 measurement is not a
production energy claim.

## Explicit Non-Goals

This authority does not approve:

- changing the production residual cutoff or any public/default behavior;
- public API, canonical driver, artifact, solver, or corrected-Hamiltonian
  workflow;
- protected-original or occupied-first replacement, residual injection, or
  automatic occupied-direction append;
- Boys or other localization policy;
- eigenvalue flooring or tolerance weakening;
- `C' V C`, an interaction rotation, or a source-Hamiltonian transform;
- Gaussian-array odd partners without a later measured occupation,
  instability, correlated-natural-occupation, or MWG-sensitivity trigger;
- EGOI changes, correlated-method claims, Cr2-specific source behavior, or a
  Cr2 production claim;
- changes to protected-localized artifacts, ladders, EGOI, or the existing
  protected additive-reference path.
