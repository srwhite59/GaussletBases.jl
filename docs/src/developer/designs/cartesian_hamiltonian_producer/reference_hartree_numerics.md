# Reference Hartree Numerics

Status: implemented neutral internal infrastructure under
`HP-RHO0-MIXH-GG-FN-01` / `HP-RHO0-MIXH-GG-TEST-01`,
`HP-RHO0-MIXH-GAAA-FN-01` / `HP-RHO0-MIXH-GAAA-TEST-01`, and
`HP-RHO0-MIXH-FEXACT-FN-01` / `HP-RHO0-MIXH-FEXACT-TEST-01`.

This page is the canonical contract for exact reference-density Hartree raw
blocks and their protected fixed/localized one-body transforms. These
facilities are neutral numerical infrastructure. They do not define `rho0`, a
screened-Hartree correction, an affine anchor, exchange policy, or corrected
Hamiltonian behavior.

Live correction ownership is separate:

- [Screened Hartree residual density](screened_hartree_residual_density.md)
  owns the `Delta_J0/C` physics;
- [Screened Hartree correction assembly](screened_hartree_correction_assembly.md)
  owns its implemented API;
- [Protected additive reference correction](protected_additive_reference_correction.md)
  owns additive molecular `P0/J0/E0` composition.

## Source Ownership

Raw exact Hartree construction is owned by:

- `src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl`;
- `src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl` for
  internal module wiring.

The existing dense oracle and reusable pair-term algebra remain in
`src/gaussian_coulomb_reference.jl`. Under these IDs that file is an
oracle/debug and narrow algebra-reuse surface only; it is not a second raw
block production owner.

Protected fixed/localized transformation is owned by:

- `src/cartesian_residual_gaussians/augmented_operators.jl`.

The primary implemented entry points are:

- `atomic_reference_hartree_gg_block`;
- `atomic_reference_hartree_ga_aa_blocks`;
- `transform_protected_original_fixed_sector_exact_hartree`;
- `transform_protected_original_localized_exact_hartree`;
- `atomic_reference_protected_original_fixed_sector_exact_hartree`.

## Reference Density And Potential Terms

Each raw construction consumes one finite, symmetric, same-center atomic
reference density matrix `P_A` in a reference Gaussian supplement. Its shape
must match the reference overlap, and normalization is reported through

```text
N_A = Tr(P_A * S_AA).
```

The source forms a vector-backed same-center reference pair-density term
stream. Angular and off-diagonal reference pairs are part of the contract;
this is not a diagonal or `s*s` approximation. The terms are convolved with
the supplied Coulomb Gaussian expansion to produce separable one-body
potential factors. Dense Gaussian Coulomb construction is an oracle/debug
path only.

The helpers accept an explicit Coulomb expansion. Their private default is the
existing compact `coulomb_gaussian_expansion(doacc = false)` convention;
callers selecting another approved Hamiltonian-wide expansion pass it
explicitly.

## Exact Raw Blocks

For one atom's reference density `P_A`, terminal/base functions `G`, and
supplement Gaussian functions `A`, the raw block contract is:

```text
GG = <G_i | v_P_A | G_j>
GA = <G_i | v_P_A | A_j>
AA = <A_i | v_P_A | A_j>
```

`atomic_reference_hartree_gg_block` constructs the symmetric terminal/final
`GG` block through terminal Gaussian-sum accumulation and returns compact
pair-term, packet-count, normalization, trace, symmetry, and finite-value
diagnostics.

`atomic_reference_hartree_ga_aa_blocks` reuses the same validated `P_A`, pair
term stream, and Coulomb-expanded potential terms. It returns parent/proxy
`GA`, symmetric `AA`, and compact diagnostics. This wrapper remains a neutral
oracle seam even when no correction workflow calls it directly.

A molecular reference Hartree field is a caller-side sum of one-center atomic
contributions. These helpers do not introduce cross-atom reference-density
pair products.

## Protected Fixed And Localized Transforms

The exact-side transform consumes the complete raw triple and existing
protected geometry:

```text
raw blocks:  GG, GA, AA
geometry:    F = [Z, M Q_perp]
output:      F_exact_Hartree[P0] in the F sector
```

`transform_protected_original_fixed_sector_exact_hartree` validates raw block
dimensions and finite values, then calls the existing exact one-body
`transform_protected_original_fixed_sector_operator`. It returns the dense
symmetric fixed-sector Hartree operator and compact geometry/symmetry
diagnostics. It does not select or modify protected geometry.

When a raw `GA` block is still in parent/proxy rows,
`atomic_reference_protected_original_fixed_sector_exact_hartree` projects it
to terminal-final `G-A` rows before applying the same fixed-sector transform.

For an existing protected localization transform `W`,
`transform_protected_original_localized_exact_hartree` defines

```text
J0_L = sym(W' * J0_F * W).
```

This is an exact one-body congruence. It is independent of the inherited-site
two-index `Vee_L` convention and does not transform an IDA/MWG interaction
matrix.

## Validation Contract

The compact accepted gates are:

- `GG`: H direct-core/scalar-screen parity, angular and off-diagonal
  same-center pair checks, bounded Be output, finite/symmetric output, and
  dense primitive-oracle spot checks;
- `GA/AA`: prior `GG` parity, angular reference pairs and angular supplement
  rows, finite dimensions, symmetric `AA`, and dense Gaussian-oracle spots;
- fixed/localized transforms: H/Be/Be2 only, complete raw-block dimension
  checks, finite/symmetric transformed output, dense-oracle spots, and compact
  protected-geometry diagnostics.

Accepted implementation evidence is recorded in manager running-log Passes
284, 286, and 288 and commits `efaee93f6`, `daac231d0`, and `40a6f7e99`.
Those passes used bounded ignored probes; the probe paths are historical
working files, not committed test ownership.

`test/nested/cartesian_screened_hartree_correction_runtests.jl` now exercises
the surviving `GG` construction as a consumer check. The accepted MIXH and
FEXACT implementation gates otherwise used the ignored bounded probes above;
there is no dedicated committed MIXH test file.

## Failure And Exclusions

Fail on invalid reference density shape, nonfinite data, raw block dimension
mismatch, or nonfinite transformed output. Keep symmetry, normalization, and
protected-geometry diagnostics available to callers.

These IDs do not authorize or own:

- row-gauge or scalar-rho0 substitutions;
- `F_app[P0]`, `Delta_F0`, `Delta_J0`, `C0`, or `C` policy;
- approximate IDA/MWG energy or Fock helpers;
- exchange or exchange/direct pairing;
- corrected Hamiltonian assembly;
- artifact, public driver, solver, or endpoint workflow;
- residual selection or protected-geometry changes;
- dense final four-index ERI production;
- Cr/Cr2 production work.

Removing or superseding a correction experiment does not retire these neutral
Hartree kernels or transforms. Their authority remains exactly the six
implemented MIXH/FEXACT IDs named at the top of this page.
