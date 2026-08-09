# Sliced Hydrogen-Chain Producer

Status: approved for bounded implementation; implementation preflight is
complete, but no source is implemented yet.

Authority IDs:

- `HP-SLICE-HCHAIN-FN-01`
- `HP-SLICE-HCHAIN-TEST-01`

This page is the canonical contract for one expert producer of minimal sliced
Hamiltonians for long, equally spaced hydrogen chains. It is a separate
one-dimensional-plus-transverse construction, not a Cartesian/PQS route, a
solver, or a general sliced-basis framework.

## Purpose

The producer represents one spatial orbital per longitudinal site,

```text
Psi_i(x,y,z) = g_i(z) phi_i(x,y),
```

where `g_i` is a uniform G10 gausslet and `phi_i` is the leading transverse
function obtained by the paper's local projected-density construction. The
result must support chains from bounded validation fixtures through H1000 and
H10000 without storing dense one- or two-electron matrices.

The bare operator depends on the nuclei and basis only. Electron count, spin
sector, HF, DMRG, and other solver choices remain consumer data.

## Paper Construction

### Longitudinal basis

Use the existing `UniformBasisSpec(:G10; xmin, xmax, spacing)` and `build_basis`
path as the authoritative G10 donor. Because the returned `UniformBasis` owns a
dense primitive-coefficient matrix, build only a bounded prototype sufficient
to validate its centers, integral weights, overlap, boundary convention, and
translated local primitive stencil. The full chain retains only that validated
stencil together with realized centers, weights, and compact boundary data. It
must not construct or retain a full-chain `UniformBasis`.

This is a compact realization of the existing G10 family, not a new gausslet
family or an independently fitted stencil. The producer must fail if the donor
prototype does not exhibit the expected translated interior stencil or if the
compact realization does not reproduce the bounded donor within the declared
tolerance. The consumer controls:

- atom count and equal nuclear spacing `R`;
- gausslet spacing, or a commensurate integer `sites_per_atom`;
- lattice phase/alignment;
- left/right padding and boundary policy; and
- numerical one-body and interaction tolerances.

The realized basis must report its centers, spacing, bounds, integral weights,
prototype identity and stencil parity, overlap error, atom/site commensurability,
and boundary facts. The construction must reject nonfinite geometry,
nonpositive spacing, inconsistent `sites_per_atom`, duplicate centers,
materially nonorthogonal longitudinal functions, or a boundary that does not
contain the requested chain and padding.

`R = 3.6`, spacing `0.2`, `sites_per_atom = 18`, and atom-centered alignment
belong only to the first paper-validation fixture. They are not defaults.

### Transverse function

Start from all relevant atom-centered contracted functions in the requested
minimal hydrogen basis. For a longitudinal center `z_i`, restrict each 3D
function to that plane,

```text
f_Ai(x,y) = chi_A(x,y,z_i).
```

In an orthonormal two-dimensional primitive span `xi_mu,i`, form

```text
eta_A,mu(i) = <xi_mu,i | f_Ai>
rho_mu,nu(i) = sum_A eta_A,mu(i) * eta_A,nu(i).
```

`phi_i` is the normalized leading eigenvector of `rho(i)`. Its phase is fixed
by making the largest-magnitude primitive coefficient positive, with the
lowest primitive index breaking ties. Fail on material negative spectrum,
lost numerical rank, an unresolved leading eigenspace at the declared
tolerance, nonfinite coefficients, or failed two-dimensional normalization.
Do not silently pick an arbitrary vector from a degenerate leading subspace.

Finite-chain mode uses the actual finite set of atom-centered contracted
functions. A periodic-template mode may reuse site classes only after proving
the translated local density matrices and phase-canonical transverse functions
agree within the declared tolerance. Any finite projection window used to
construct a periodic template must be enlarged until the omitted projected
norm is below that tolerance.

The initial producer supports exactly one transverse function per site and the
vendored contracted H/STO-6G fixture. Other minimal hydrogen bases may be
accepted only when the same one-function numerical-rank contract passes; there
is no multi-orbital slice policy in this lane.

## One-Body Operator

The one-body matrix is full analytic Galerkin:

```text
h_ij = <g_i phi_i |
          -1/2 nabla^2 - sum_A Z_A / |r-R_A|
        | g_j phi_j>.
```

It must include, without a one-body diagonal approximation:

```text
Tz_ij * Sperp_ij
+ Sz_ij * Tperp_ij
+ sum_A U_A,ij.
```

The implementation must retain the full analytic transverse overlap and
transverse kinetic factors, including off-diagonal values, even where exact
longitudinal orthogonality makes a product numerically small. Every nuclear
term uses the physical nucleus and the complete three-dimensional product
functions. Allowed evaluation paths are analytic Gaussian formulas and a
declared, validated Gaussian representation of `1/r`. A fixed Gaussian
expansion may be used only over its validated distance range; distant-nucleus
contributions require a checked analytic or asymptotic continuation whose
transition error and bound satisfy the requested one-body tolerance. The fixed
expansion's eventual false decay must never truncate the nuclear attraction of
a long chain. Finite differences, sampled slice energies, point evaluation,
diagonal slice-energy reduction, and one-body IDA are forbidden.

Nuclear repulsion is a separately reported scalar. No neutrality or electron
sector is required to construct `h`, the interaction resource, or nuclear
repulsion.

## Integral-Diagonal Interaction

Only the electron-electron interaction receives a longitudinal integral
diagonal approximation. For normalized transverse density
`rho_i_perp = |phi_i|^2`, define

```text
Kperp_ij(z,z') = integral dx dy dx' dy'
    rho_i_perp(x,y) rho_j_perp(x',y') /
    sqrt((x-x')^2 + (y-y')^2 + (z-z')^2).
```

With `w_i = integral g_i(z) dz`, the retained density-density matrix is

```text
V_ij = integral dz dz' g_i(z) Kperp_ij(z,z') g_j(z') / (w_i w_j),

V_ijkl -> delta_il delta_jk V_ij.
```

This mirrors the existing `IntegralDiagonal` convention. It is neither point
IDA nor the exact retained integral `(ii|jj)`. The implementation must reject
nonfinite or numerically zero longitudinal integral weights and must preserve
symmetry and positive Coulomb self-values within tolerance.

The accuracy contract covers the largest realized site separation. A fixed
Coulomb Gaussian expansion may be used only over its validated range. At
larger separation, use a validated analytic evaluation or an asymptotic
multipole/`1/R` tail with a checked transition and an error bound below the
requested interaction tolerance. This long-range rule applies independently to
both Vee and the distant-nucleus H1 terms above. Silently allowing a Gaussian
expansion to decay to zero at H10000 distances in either operator is forbidden.

## Compact Representation And Expert Interface

One exported concrete owner is approved:

```text
SlicedHydrogenChain
```

It is built by the exported expert constructor

```text
sliced_hydrogen_chain(...)
```

and exposes two lightweight operator views:

```text
sliced_h1(chain)
sliced_vee(chain)
```

The views support `size`, scalar `getindex`, `mul!`, and the exported
caller-buffer operation

```text
sliced_row!(destination, view, row)
```

After construction, scalar access must allocate nothing and row/action access
must allocate nothing beyond caller-supplied output and explicitly documented
reusable workspace. Neither view may own an `N x N` matrix.

The producer object may own only compact vector-backed numerical state: the
validated translated G10 stencil, realized centers and weights, compact
boundary data, nuclei, normalized transverse primitive data, validated
periodic-class data when requested, finite-edge data, Coulomb/tail resources,
tolerances, and a small fixed diagnostics summary. It must not retain the
full-chain `UniformBasis`, electron-sector state, a solver, an MPO, dense
matrices, all-pairs metadata, runtime-keyed `NamedTuple` inventories, or
recursive construction stages.

For a translated periodic transverse template with `b` sites per atom,

```text
V[(cell,alpha),(cell+d,beta)] = K[alpha,beta,d].
```

The optimized interaction state is one `b x b` block row over realized cell
separations, equivalent to `b` scalar rows and linear storage at fixed `b`.
Translated scalar and row access must agree within tolerance. Finite open-chain
mode uses the same public interface but makes no translation-invariance claim
and still may not store an `N x N` matrix. A periodic-interior/finite-edge
optimization is optional and may not change matrix elements.

## Provenance And Diagnostics

The in-memory object reports the resolved chain geometry, G10 family and
spacing, phase/alignment, padding/boundary policy, minimal-basis identity and
fingerprint, transverse construction mode, local density ranks and leading
eigenvalue gaps, phase pivots, normalization errors, periodic-class errors,
Coulomb evaluator identity, analytic/asymptotic transition facts, maximum
validated separation, and requested/achieved tolerances.

These are fixed object diagnostics, not an artifact schema or a new metadata
family. The producer writes no files.

## Implementation Ownership And Budget

Approved source surfaces are exactly:

- new `src/sliced_hydrogen_chain.jl`;
- `src/GaussletBases.jl`, only for one include and the five exported names
  listed above.

The implementation may call existing `UniformBasisSpec`, G10 primitive/stencil
data, `IntegralDiagonal` convention, legacy contracted-basis loading,
`GaussianAnalyticIntegrals`, and `ordinary_coulomb` kernels. It must not broaden
their contracts or copy their coefficient tables.

The implementation-preflight-adjusted budget is:

- preferred: at most `480` added source lines in total;
- hard: at most `520` added source lines in total;
- tests: at most `240` added lines in existing `test/ida/runtests.jl`;
- no other added source, test, tool, driver, helper, or module file.

The complete producer, both operator views, finite/periodic behavior, and the
validated long-range interaction policy must land together. Do not split
incomplete scaffolding to evade the budget. If a readable implementation
cannot fit the hard limit, make no source commit and report the smallest
specific reusable-kernel extraction or revised budget required.

## Validation Contract

The committed bounded validation in `test/ida/runtests.jl` must cover:

1. contracted H/STO-6G transverse functions with deterministic phases,
   normalization, finite spectra, and explicit rank/gap diagnostics;
2. bounded donor-prototype parity for the translated G10 stencil, followed by
   compact full-chain overlap, centers, integral weights, and boundary facts;
3. finite, Hermitian analytic `h` against an independently assembled bounded
   dense primitive-contraction oracle, including selected distant-nucleus
   transition checks;
4. selected onsite, near, translated-cell, and far `V_ij` entries against
   direct quadrature at the declared tolerance;
5. scalar/row/action parity and zero-allocation steady-state calls;
6. translated periodic rows and finite open-chain behavior;
7. the `R=3.6`, spacing `0.2`, `b=18`, atom-aligned H10 dense oracle;
8. H1000 construction with measured subquadratic storage, no full-chain
   `UniformBasis`, and no `N x N` allocation; and
9. malformed geometry, rank, phase, boundary, weight, tail-range, and
   tolerance rejection.

H10000 is a transient scaling gate, not a committed test. Report construction
time, steady-state scalar/row/action allocation, retained storage by category,
maximum separation, tail transition, and achieved errors. No HF, DMRG, MPO,
or publication energy is an acceptance requirement.

## Explicit Non-Goals

This authority does not permit:

- a provider hierarchy or general sliced-basis framework;
- a driver, artifact, sidecar, payload, schema, or persistent report;
- a solver, HF/DMRG loop, MPO/SVD compressor, or electron-sector policy;
- finite differences, a soft-Coulomb model, one-body IDA, or dense H1/Vee;
- multiple transverse functions per site;
- Cartesian/PQS shellification, PRFs, residual Gaussians, screening, or EGOI;
- use or expansion of `bond_aligned_homonuclear_chain_qw_basis`;
- use or expansion of `sliced_ham_export.jl`; or
- an H1000/H10000 energy or publication claim.

## References

- E. M. Stoudenmire and S. R. White,
  [Sliced Basis Density Matrix Renormalization Group for Electronic Structure](https://arxiv.org/abs/1702.03650).
- S. R. White,
  [Hybrid Grid/Basis Set Discretizations of the Schrodinger Equation](https://arxiv.org/abs/1709.08059).
- R. C. Sawaya and S. R. White,
  [Constructing Hubbard Models for the Hydrogen Chain Using Sliced-Basis DMRG](https://arxiv.org/abs/2109.05129).
