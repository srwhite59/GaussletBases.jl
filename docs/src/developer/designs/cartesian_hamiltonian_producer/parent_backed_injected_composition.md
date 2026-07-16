# Parent-Backed Injected Composition

Status: implemented internal facility with completed bounded validation.

Authority IDs:

- `HP-PQS-PRF-INJECT-COMP-FN-01`
- `HP-PQS-PRF-INJECT-COMP-TEST-01`
- `HP-PQS-PRF-INJECT-INTERACT-FN-01`
- `HP-PQS-PRF-INJECT-INTERACT-TEST-01`

This page is the canonical contract for composing consumer-selected
parent-backed injection with the numerical-complete Gaussian supplement and
the existing screened-Hartree correction. It does not select a physical state,
shell, source order, injection target, residual orientation, or solver method.

## Purpose And Evidence

The repo owns support-local parent residual functions (PRFs), exact
parent-backed one-body blocks, onsite-calibrated Gaussian direct blocks, and
the category-correct composition in which:

- injected localized rows with positive linear IDA weights remain terminal
  gausslet-like sites;
- their exact parent-backed complement remains a separate PRF sector;
- a Gaussian supplement is residualized against the complete parent-backed
  span; and
- interaction blocks are assembled by their appropriate existing convention.

The accepted Cr2 measurement used a fixed q7 parent span with eight PRFs per
owner. Moving the occupied-important content into injected terminal rows
reduced the fixed-state density-interaction error from about `4.081 mHa` to
`1.463 mHa` without changing the parent-backed span or dimension. After adding
the numerical-complete cc-pV5Z complement and the existing additive screened
Hartree correction, a bounded relaxed UHF result remained in the intended spin
basin and lay about `1.105 mHa` above matched PySCF.

This is evidence for reusable mechanics, not a q7, shell-9, eight-mode,
`1e-5`, Cr2, or energy policy. The corresponding q6 interaction error remained
large, so no automatic contraction or injection choice is approved.

## Basis Objects And Native Order

Let `G` be the current orthonormal localized terminal basis and let `R` be one
or more consumer-supplied, validated PRF blocks in the same mapped parent
construction. Define the original parent-backed span

```text
S_parent = span(G, R).
```

The consumer supplies an already-selected orthonormal target block `Y` inside
`S_parent`. Selection may have used an external RDM or other consumer
analysis, but the repo helper receives columns only. It must not receive or
interpret the RDM, occupations, shell-selection policy, or cutoff that produced
them.

The final native order is

```text
B = [G_inj, R_new, RG_external].
```

`G_inj` and `R_new` remain distinct categories even though both have explicit
parent coefficients. PRFs must not be appended to
`CartesianTerminalBasisRealization` and then passed through terminal IDA merely
to reuse terminal operator code. The implemented internal
`CartesianParentBackedInjectedComposition` holds the terminal realization,
PRF blocks, column ranges, dimension, and compact validation summary. It is the
sole composition object added by this lane; no public type, module, source
file, or persistent metadata cloud was added.

## Fixed-Span Injection

For each consumer-selected source block or symmetry-related block set, require
full target rank in

```text
C = G' S Y.
```

Let `Q_perp` be an orthonormal basis for `null(C')`. Construct

```text
Q      = G Q_perp
F      = [Y, Q]
Gbar   = P_F G = F (F' S G)
G_inj  = Gbar (Gbar' S Gbar)^(-1/2).
```

The projection of every old localized `G` seed is mandatory. Returning an
arbitrary orthonormal basis of `F`, including `[Y,Q]`, proves only span
replacement and is not the angular-style injected localized basis.

Construct the exact complement inside the original parent-backed span:

```text
R_new = S_parent intersect span(G_inj)_perp.
```

The result must preserve dimension and span:

```text
dim(G_inj)          = dim(G)
dim(R_new)          = dim(R)
span(G_inj, R_new)  = span(G, R).
```

The helper must validate target recovery, old/new principal singular values,
metric identity, cross orthogonality, exact shell support, symmetry-related
counts, and unchanged unaffected blocks. Every `G_inj` row must have a finite
positive terminal IDA weight. Rank loss, nonpositive weight, support leakage,
or a material span mismatch is a hard failure.

The consumer owns target orientation and any same-span orientation of `R_new`.
The repo may provide a deterministic metric-orthonormal complement and phase
convention, but it must not optimize an RDM, localize the complement, discard a
column, or choose a physical cutoff.

## Numerical-Complete External Residual

Let `A` be the explicit Gaussian supplement. Form exact mixed overlap against
the entire parent-backed span, not only the injected terminal rows:

```text
X_BA  = <[G_inj,R_new] | A>
M_a   = S_AA[a,a] - X_BA[:,a]' X_BA[:,a].
```

Reuse `build_residual_gaussian_basis(...)` and the numerical-complete policy:

```text
residual_occupation_cutoff = 1e-10
residual_injection_cutoff  = 0
residual_compactness       = nothing.
```

Owner-local rank selection and the one final global Lowdin merge remain the
existing algorithm. No second residual builder or object is approved. The
native unlocalized merge remains the repo construction default. A consumer may
supply an otherwise-authorized same-span residual rotation, but this lane does
not choose or optimize one and must rebuild all moment-derived interaction
blocks after any such rotation.

Packet occupied spaces validate capture only after the complete `B` basis is
built. They do not select or append residual directions. Capture failure is a
supplement/rank failure and must not trigger occupied injection or a relaxed
threshold.

## Exact One-Body Composition

Build exact kinetic, per-center unit-nuclear, physical `H1`, position, and
second-moment blocks in native `B` order through the existing factorized
parent and augmented-operator owners. Reuse the implemented PRF `G-R`/`R-R`
one-body blocks and the existing Residual-Gaussian transforms.

The source must not form a dense parent operator, transform a source
Hamiltonian, use a generalized final overlap, or infer interaction semantics
from an exact one-body transform. Omitted composition must preserve the current
terminal and numerical-complete `H1` paths exactly.

## Interaction Block Contract

Assemble the unscreened density interaction explicitly by category:

| Block | Required convention |
| --- | --- |
| `G_inj-G_inj` | Rebuild ordinary terminal/site IDA from the injected coefficients and their positive final weights |
| `G_inj-R_new` | Existing onsite-calibrated parent-Gaussian direct block |
| `R_new-R_new` | Existing onsite-calibrated parent-Gaussian direct block |
| `G_inj-RG_external` | Existing terminal-to-residual MWG block |
| `R_new-RG_external` | Existing moment-derived residual MWG block |
| `RG_external-RG_external` | Existing residual-to-residual MWG block |

`G_inj-G_inj` is reassembled through terminal IDA. It is neither copied from
the old `G-G` matrix nor obtained by `C' V C`. `R_new` has near-zero linear
parent weight by design and must bypass the positive terminal-weight division.

For `R_new-RG_external`, derive centers and widths from the exact parent-backed
position and second-moment matrices and the existing residual moments. Require
finite centers, finite positive widths, symmetry, and finiteness of the final
block. This is a moment-matched density-density approximation, not an exact
PRF-GTO integral.

The parent-Gaussian resource remains direct-only. A density-density
Hamiltonian may use these values in an exchange-like contraction, but the
facility does not construct transition products `(ai|ia)` and makes no
continuum-exact exchange claim.

## Screened-Hartree Delegation

Construct screening only after `B`, exact `H1_B`, and unscreened `Vee_B` are
fixed. Reuse separately placed atomic reference packets; do not perform a
molecular density or potential refit. Represent each original occupied packet
block separately in `B`, then use the existing additive field and energy
construction:

```text
Delta_J0 = J0_B - Diagonal(Vee_B * q0_B)
C        = 0.5*q0_B' * Vee_B * q0_B - 0.5*E0.
```

Return the existing `ScreenedHartreeCorrection` separately. Screening off/on
uses identical `H1_B` and `Vee_B`. Any basis, PRF orientation, residual
rotation, interaction, packet placement, or Coulomb-policy change invalidates
the old correction and requires rebuilding it.

## Implemented Source Boundary

Commits `cdd2c27af` and `006432e9d` implemented the facility in these existing
owners:

- `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
  owns the composition object, fixed-span terminal replacement, exact
  complement, and compact validation;
- `src/cartesian_final_basis_realization/pqs_terminal_one_body.jl` owns exact
  native parent-backed one-body operators;
- `src/cartesian_final_basis_realization/pqs_terminal_ida.jl` owns rebuilt
  terminal IDA and parent-Gaussian base interaction blocks;
- `src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl` owns
  parent-backed-to-GTO overlap, numerical-complete augmentation, native exact
  operators, separated interaction composition, and fitted-field
  representation;
- `src/cartesian_residual_gaussians/mwg_interaction.jl` owns the separated
  terminal, parent-residual, and external-residual MWG blocks; and
- `src/cartesian_protected_ladder_bundle.jl` owns the private in-memory
  composition and additive screened-Hartree delegation seam.

The implementation consumes the existing numerical-complete residual builder,
parent-Gaussian direct resource, augmented transforms, and screened-Hartree API
without moving their ownership into these IDs. No root export, public facade,
driver input, artifact field, restart format, module, or source file was added.

## Validation Contract

Committed validation ownership is limited to
`test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl`. Its bounded H2
gate owns construction, malformed and stale-input failures, exact one-body
oracles, all separated interaction blocks, and omitted-path parity.

Required gates are:

1. A synthetic fixture distinguishes complete old-seed projection followed by
   symmetric Lowdin from merely returning an orthonormal basis of `F`.
2. Malformed, rank-deficient, support-leaking, dimension-changing,
   nonorthogonal, nonpositive-weight, and span-changing inputs fail without
   dropping or repairing columns.
3. Omitted composition preserves the existing terminal, numerical-complete,
   additive-reference, `H1`, and `Vee` paths exactly.
4. Bounded H2 validates target recovery, positive injected weights, unchanged
   unaffected blocks, old/new span singular values, exact one-body oracle
   parity, and finite symmetric separated interaction blocks.
5. Padded Be2 residualizes the complete supplement against
   `[G_inj,R_new]`, validates packet capture and native ordering, recomputes the
   correction after interaction assembly, and inspects terminal due diligence.
6. Every `R_new-RG_external` center/width is finite with positive width; block
   orientation, dimensions, symmetry, and finiteness agree with a bounded
   explicit existing-MWG construction.

Manager rerun passed the committed H2 gate `434/434` and the unchanged facade
regression `69/69`. The padded Be2 construction is external manager-reviewed
evidence rather than a committed slow fixture. It used native dimensions
`1729 + 2 + 42 = 1773`; packet capture/accounting, additive decomposition,
correction recomputation, derivative anchor, unchanged `H`/`V`, and terminal
due diligence passed. The facade regression and external Be2 gate do not
expand these IDs' test ownership.

No committed Be2 or Cr2 fixture, HF convergence, endpoint energy, or exchange
assertion is approved.

## Failure And Stop Rules

Stop and report if implementation requires:

- automatic target, shell, source-order, RDM, cutoff, count, or orientation
  selection;
- representing `R_new` as positive-weight terminal sites;
- changing the existing numerical-complete rank or merge algorithm;
- rotating an old interaction matrix instead of rebuilding category blocks;
- weakening metric, span, support, capture, weight, or symmetry tolerances;
- a dense parent/final operator, new source file/module, second residual type,
  or persistent metadata field cloud; or
- source edits outside the listed owners.

Do not make a failed construction pass by dropping targets or PRFs, flooring
eigenvalues, renormalizing malformed charges, inventing PRF IDA weights, or
reusing a correction from another basis or interaction.

## Explicit Non-Goals

These IDs do not approve:

- any q7, q6, shell-9, eight-PRF, `1e-5`, element, molecule, or geometry
  default;
- consumer RDM construction, spin averaging, occupation analysis, automatic
  selection, or localization policy;
- public API, facade, canonical driver, artifact, sidecar, restart, solver, or
  HF workflow;
- `C' V C`, generalized-overlap interaction, source-Hamiltonian transform, or
  PRF-as-terminal compatibility adapter;
- exact PRF-GTO direct integrals, transition-density exchange, exact exchange,
  or a new exchange correction;
- screened-Hartree formula, packet-fit, Coulomb-policy, EGOI, or protected
  artifact changes; or
- Cr2 production energies, correlated-method claims, paper claims, or a
  complete continuum-accurate Hamiltonian claim.
