# Reference-Density Rho0 Implementation Plan

Status: planning memo for review. This is not source authority, artifact
authority, public API authority, solver authority, or Cr2 production authority.

This memo records the intended implementation shape for the
reference-density-matrix correction lane before source work begins. It should
be reviewed before approving any source ID.

Revision note: this version incorporates the initial design-manager review and
the external ChatGPT review saved on 2026-07-03. The main change is that the
first source target is now explicitly `GG`-only. The full `GG` / `GA` / `AA`,
protected-transform, approximate-Fock, and constant-assembly path remains the
longer plan, not the first implementation slice.

Governing authority today remains:

- `HP-RHO0-REFDENS-AUDIT-01` - measurement-only fixed `P0` audit;
- `HP-RHO0-REFDENS-MIXH-AUDIT-01` - measurement-only mixed Hartree seam audit.

Candidate future IDs such as `HP-RHO0-REFDENS-FN-01` and
`HP-RHO0-REFDENS-ERI-01` are not approved.

## Current Conclusion

The rho0 correction should be built around a fixed reference density matrix
`P0`, not a scalar row-gauge test. The correction is:

```text
Delta_F0_sigma = F_exact0_sigma[P0] - F_app0_sigma[P0]

C0 =
    E_exact0[P0]
  - E_app0[P0]
  - sum_sigma Tr(P0_sigma * Delta_F0_sigma)
```

The corrected model must satisfy:

```text
E_new[P0] = E_exact0[P0]
dE_new/dP_sigma at P0 = F_exact0_sigma[P0]
```

The old row-action diagnostics are retired as acceptance criteria:

```text
J*w = u*w      not a rho0 gauge requirement
diag(J) = u    not a rho0 gauge requirement
C' V C         invalid for two-index IDA/MWG interactions
```

## One-Center Atomic Reference Density

The intended reference density is a sum of atomic densities:

```text
rho0(r) = sum_A rho_A(r - R_A)
```

Each atomic density is built from same-center atomic GTO products:

```text
rho_A(r) = sum_ab P_A[ab] chi_Aa(r) chi_Ab(r)
```

This means:

- no cross-atom reference density blocks `chi_Aa chi_Bb` are needed;
- no molecular orthogonalization of atomic reference orbitals is implied;
- `P_A` is a density matrix for an atomic density, not a molecular Slater
  determinant block.

However, Coulomb interactions between atomic densities are still required:

```text
E_H[rho0] =
  1/2 sum_A   (rho_A | rho_A)
  +   sum_A<B (rho_A | rho_B)
```

The exact Hartree operator on the final basis is:

```text
F_exact[P0]_{ij} =
    sum_A sum_ab P_A[ab] (phi_i phi_j | chi_Aa chi_Ab)
```

So the reference pair products are one-center, but their potential is felt
everywhere in the molecule.

## Trace And Density Convention

Every audit and future source caller must state the density convention.

If `P_A` is in a nonorthogonal atomic GTO basis, electron count is normally:

```text
N_A = Tr(P_A S_A)
```

If `P_A` is already in an orthonormal atomic orbital basis, the code must say
how it was converted back to AO pair-product coefficients before forming
`chi_Aa chi_Ab`.

Required metadata for each atomic `P_A`:

- atom / nuclear center;
- reference electron count;
- spin convention: spin-summed or spin-resolved;
- full atom, core-only, or valence-screened;
- GTO basis, `lmax`, contraction flag, and basisfile if relevant;
- whether `P_A` came from internal or external atomic calculation;
- normalization check, including `Tr(P_A S_A)` when applicable.

### Two Distinct `P0` Objects

Keep these objects separate in every source design and audit:

```text
P_A       atomic reference density in the atomic AO/GTO basis
P0_final  represented reference density in the final/protected basis
```

`P_A` is used to build the exact reference density:

```text
rho_A(r) = sum_ab P_A[ab] chi_Aa(r) chi_Ab(r)
```

`P0_final` is the object that can appear in final-basis traces such as:

```text
Tr(P0_final * Delta_F0)
```

Do not put the atomic AO matrix `P_A` directly into a trace with a final-basis
operator. If the correction constant or derivative anchor is evaluated in the
final basis, the represented/projected density must be in that same basis, or
the probe must state an explicitly equivalent representation.

## Production-Shaped Math

The production-shaped seam should not build a dense four-index final-basis ERI.
It should use the fast separable Coulomb path:

```text
1/r ~= sum_t c_t exp(-alpha_t |r-r'|^2)
```

For each one-center pair density `chi_a(r') chi_b(r')`, and each Coulomb
Gaussian term, analytically integrate the reference side over `r'`. The result
is a separable polynomial-Gaussian one-body potential in `r`:

```text
chi_a(r') chi_b(r') * exp(-alpha_t |r-r'|^2)
    integrated over r'
  -> product over axes of polynomial-Gaussian factors in r
```

Those one-body factors should then be contracted through the existing
Cartesian product machinery:

```text
reference atomic P0
  -> weighted one-center pair-density terms
  -> Coulomb-expanded separable one-body terms
  -> GG / GA / AA exact Hartree blocks
  -> existing protected/final transform
```

Dense Gaussian ERIs remain useful only as oracle/debug for small systems.

## Existing Ingredients

### Pure Gaussian Coulomb Oracle

`src/gaussian_coulomb_reference.jl`

Useful functions and private ingredients:

- `gaussian_coulomb_pair_matrix(...)`;
- `_gaussian_coulomb_pair_terms(...)`;
- `_gaussian_coulomb_pair_integral(...)`;
- `_gaussian_coulomb_pair_matrix_compressed_checked(...)`;
- `_gaussian_coulomb_global_term_kernel(...)`.

This code already knows how to build polynomial Gaussian pair terms and exact
Coulomb couplings between Gaussian pair products. It is the right oracle and
the right source of pair-term algebra.

It is not the production final-basis Hartree path because dense pair matrices
scale like `N^4`.

### Cartesian Gaussian Raw Blocks

`src/cartesian_gaussian_raw_blocks/non_nuclear_blocks.jl`

- `gaussian_non_nuclear_raw_blocks(...)`;
- `_non_nuclear_ga_block_matrices(...)`;
- `_non_nuclear_aa_block_matrices(...)`.

`src/cartesian_gaussian_raw_blocks/nuclear_blocks.jl`

- `gaussian_nuclear_raw_blocks_by_center(...)`;
- `_fill_axis_factor_table!(...)`.

The nuclear block owner already builds `G-A` / `A-A` one-body blocks from a
Coulomb Gaussian expansion and centered factor terms. That is close in shape,
but the current path is center-factor based. The new seam must generalize the
factor side from a point/center screen to arbitrary one-center GTO pair-density
terms, including angular and off-diagonal products.

### Terminal Gaussian-Sum Accumulation

`src/cartesian_final_basis_realization/pqs_terminal_one_body.jl`

- `_accumulate_terminal_gaussian_sum!(...)`;
- `_terminal_gaussian_sum_action(...)`.

This is the fast terminal contraction shape: term coefficients and three
axis-factor packets are contracted over terminal support states and final
terminal coefficients. The reference-density Hartree seam should reuse this
style rather than materializing dense parent matrices.

### Protected/Residual Final Transform

`src/cartesian_residual_gaussians/augmented_operators.jl`

- `transform_augmented_operator(...)`;
- `protected_original_fixed_sector_components(...)`;
- `transform_protected_original_fixed_sector_operator(...)`.

The mixed Hartree owner should produce exact `GG`, `GA`, and `AA` one-body
blocks. Residual/protected code should consume those blocks and transform them
into the protected-localized final sector. It should not own the reference
GTO-pair math.

## Proposed Source Shape

The likely future source owner should be neutral, adjacent to:

```text
src/gaussian_coulomb_reference.jl
src/cartesian_gaussian_raw_blocks/
```

It should not live in residual Gaussian code.

The first reviewed source target should be narrow:

```text
reference_hartree_GG_block(
    terminal_basis_or_proxy,
    reference_supplements_by_atom,
    P0_by_atom;
    expansion = coulomb_gaussian_expansion(doacc = false),
) -> (; GG, diagnostics)
```

or an equivalent internal name under a neutral raw-block owner such as:

```text
src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl
```

The first source slice should prove only:

```text
one-center atomic P_A
  -> same-center reference pair-density term stream
  -> Coulomb-expanded separable one-body terms
  -> terminal/base GG Hartree block
```

Do not include `GA`, `AA`, protected-localized transforms, `F_app[P0]`, `C0`,
artifacts, public workflow, Cr/Cr2 diagnostics, or HF exchange in the first
source slice.

A later full mixed-block API may look like:

```text
mixed_hartree_blocks(
    terminal_basis_or_proxy,
    optional_supplement,
    reference_supplements_by_atom,
    P0_by_atom;
    expansion = coulomb_gaussian_expansion(doacc = false),
) -> (; GG, GA, AA, reference_self_energy, diagnostics)
```

That is a future shape, not the first implementation target.

For the later full path, two reasonable options remain:

1. Return raw `GG`, `GA`, `AA` blocks and let existing augmented/protected
   transform code build final-sector `F_exact`.
2. Return final-sector `F_exact` only for a higher-level convenience wrapper,
   backed internally by the raw blocks.

The first option is cleaner for reuse and validation.

## Proposed Algorithm

### 1. Normalize Atomic Reference Density

For each atom `A`:

1. Read or construct same-center atomic GTO supplement `chi_A`.
2. Read or construct `P_A`.
3. Validate symmetry, finite entries, spin convention, and electron count.
4. Compute or consume the atomic overlap `S_A`.
5. Report `Tr(P_A S_A)` and any trace discrepancy.
6. Keep only same-center pair products `(a,b)` for atom `A`.

Do not add cross-atom reference products unless a future design explicitly
switches from atomic density to molecular reference density.

### 2. Build Reference Pair-Density Terms

For each atom `A` and GTO pair `(a,b)` with non-negligible `P_A[ab]`:

1. Expand contracted `chi_Aa chi_Ab` into primitive pair terms.
2. Combine coefficients with `P_A[ab]`, including symmetric off-diagonal
   accounting.
3. Represent each pair term as a separable polynomial Gaussian density.
4. Compress repeated axis descriptors where possible.

This should reuse or adapt the pair-term logic in
`gaussian_coulomb_reference.jl`. The source seam should make the pair-density
term stream explicit enough to share between:

- reference self-energy oracle;
- mixed Hartree one-body operator;
- diagnostics.

### 3. Convert Pair Density To One-Body Factor Terms

For each Coulomb expansion term and each reference pair-density term, integrate
the reference coordinate analytically. The output is a separable one-body
factor packet on the final-basis coordinate.

This is the key fast path. It must support:

- spherical `s*s` products;
- angular products such as `p*p`, `s*d`, `d*d`;
- off-diagonal same-center products;
- off-center atomic centers;
- positive and negative density-matrix coefficients.

The factor packet should be compatible with terminal/product contraction:

```text
coefficients[t]
factor_terms_x[t, ix, jx]
factor_terms_y[t, iy, jy]
factor_terms_z[t, iz, jz]
```

or an equivalent blocked/vector-backed representation that avoids huge
runtime-keyed tuple types.

### 4. Build Exact Raw Hartree Blocks

For the first source slice, build only:

```text
GG = <G_i | v_P0 | G_j>
```

For later slices, extend the same term stream and factor packets to optional
Gaussian supplement `A`:

```text
GA = <G_i | v_P0 | A_j>
AA = <A_i | v_P0 | A_j>
```

Requirements:

- use the same Coulomb expansion convention for all blocks;
- preserve symmetry for `GG` and `AA`;
- report finite/symmetry errors;
- block over terms and support states so memory does not explode;
- expose enough diagnostics to compare against dense Gaussian oracle on small
  systems.

If the first `GG` slice cannot be implemented without a dense final-basis
four-index tensor, stop. Do not expand scope to hide that failure in later
`GA` / `AA` or protected-transform code.

### 5. Build Reference Self-Energy

Compute:

```text
E_exact0_H =
    1/2 sum_A,B sum_ab,cd
        P_A[ab] P_B[cd] (chi_Aa chi_Ab | chi_Bc chi_Bd)
```

For small systems, `gaussian_coulomb_pair_matrix(...)` can serve as the
oracle. A future fast path can share the same pair-term compressed kernel if
needed, but the first source seam may use the dense oracle for reference
self-energy if bounded to small atomic `P0` sizes.

### 6. Transform To Protected/Final Sector

Given `GG`, `GA`, and `AA`, use existing augmented/protected one-body
transform code to form final-sector `F_exact[P0]`.

This keeps reference-density math out of residual Gaussian code and preserves
the current protected-localized convention:

- compact RG rows are part of the final basis and receive the correction;
- injected/protected rows are main-basis replacement directions, not residual
  MWG channels;
- rejected broad directions remain out of the Hamiltonian.

### 7. Build Approximate Reference Side

`F_app[P0]` must come from the actual IDA/MWG Hamiltonian convention used by
the solver for the represented `P0`.

Do not replace it with:

- `(J*w)/w`;
- `diag(J)`;
- `q0` fitting;
- center metadata;
- row-vector subtraction;
- `C' V C`.

If the approximate Fock from arbitrary represented `P0` is not accessible,
that is a separate missing seam and must be reported.

### 8. Assemble Correction

After exact and approximate derivative checks pass:

```text
Delta_F0 = F_exact0 - F_app0
C0 = E_exact0 - E_app0 - Tr(P0 * Delta_F0)
```

For spin-resolved later:

```text
C0 = E_exact0 - E_app0 - sum_sigma Tr(P0_sigma * Delta_F0_sigma)
```

## Validation Ladder

### Slice A - Dense Oracle And Pair Terms

Purpose: prove same-center GTO pair-density term representation.

Checks:

- one-center H/Be `s` reference reproduces existing limited measurement;
- angular same-center pairs agree with `gaussian_coulomb_pair_matrix`;
- off-diagonal same-center pairs agree with dense oracle;
- self-energy finite/symmetric and trace convention reported.

### Slice B - Fast Terminal `GG` Hartree Blocks

Purpose: build exact `GG` one-body blocks from atomic `P0` using fast
separable contraction.

Checks:

- H direct-core and H q5 compare to limited `s` reference measurement;
- at least one same-center angular/off-diagonal pair (`s*p`, `p*p`, or
  similar) agrees with dense Gaussian oracle;
- small Be/Be2 terminal `GG` block finite/symmetric;
- dense oracle spot checks for small support subsets;
- no dense final-basis four-index tensor.

### Slice C - `GA` / `AA` Supplement Blocks

Purpose: support compact RG/protected-original final bases.

Checks:

- `GA` and `AA` finite/symmetric where applicable;
- compare supplement-only `AA` against dense Gaussian oracle;
- compare `GA` spot rows against direct analytic/probe reference;
- no residual/MWG channel reinterpretation.

### Slice D - Protected-Localized `F_exact[P0]`

Purpose: transform exact raw blocks into final/protected sector.

Checks:

- Be/Be2 protected-localized `P0` audit;
- representability loss reported before correction use;
- exact and approximate finite-difference checks;
- `E_new[P0]` and derivative anchor pass;
- `Delta_F0` spectra and occupied expectations by sector.

### Slice E - Cr Atom / Cr2 Measurement

Only after H/Be/Be2 pass.

Purpose: measurement-only Cr diagnostics, not production.

Checks:

- reference density definition, likely atomic/core variants;
- representability and lost channels;
- `Delta_F0` broad/compact sector spectra;
- bounded behavior only after static diagnostics look sane.

## Performance Design Notes

The source seam must be shaped for speed from the start.

Do:

- keep Coulomb expansion separable by axis;
- compress repeated reference pair terms and axis descriptors;
- reuse terminal Gaussian-sum contraction patterns;
- build `GG` first, then extend to `GA/AA` only after the terminal seam is
  trusted;
- block over terms/support states;
- keep dense Gaussian ERI only as small oracle/debug;
- use vector-backed term lists and summaries, not runtime-size `NamedTuple`
  inventories.

Do not:

- materialize `(final final | reference reference)` as a dense four-index
  final tensor;
- hide source kernels in ignored probes after the design is approved;
- place neutral GTO-pair Hartree ownership in residual Gaussian code;
- create artifact/provenance fields before the Hamiltonian convention is
  stable.

## Open Design Questions

These should be answered before source authority:

1. Source owner and filename:
   - new neutral module?
   - extension of `cartesian_gaussian_raw_blocks`?
   - companion to `gaussian_coulomb_reference.jl`?

2. Public/internal API shape:
   - first source slice returns `(; GG, diagnostics)` only?
   - later full mixed-block slice returns raw `GG/GA/AA`?
   - final-sector convenience wrapper remains later?

3. Pair-density term representation:
   - can existing `_GaussianCoulombPairTerm3D` be generalized without
     exporting private internals?
   - should a new compact reference-pair term struct be introduced?

4. Polynomial-Gaussian one-body factor packets:
   - reuse terminal Gaussian-sum packet shape?
   - add a generalized polynomial factor packet?
   - how to keep memory bounded?

5. Reference self-energy:
   - dense oracle acceptable for first source slice?
   - or should self-energy use the same compressed term stream immediately?
   - if included in the first slice, is it explicitly bounded to H/Be-size
     oracle validation only?

6. Approximate `F_app[P0]` seam:
   - is the solver convention already callable for a provided represented
     density?
   - if not, does that need its own source authority?

7. Representability:
   - where should projection/loss diagnostics live?
   - how to classify loss by owner/angular channel without adding field
     clouds?

8. Spin:
   - first source slice Hartree spin-summed only?
   - when to introduce spin-resolved exchange candidate?

## Stop Rules

Stop and request design review if:

- exact mixed Hartree construction would require dense final four-index
  tensors;
- the first source slice starts growing past terminal `GG` before the term
  stream and terminal contraction are validated;
- `F_app[P0]` cannot be obtained from the actual approximate Hamiltonian
  convention;
- one-center atomic density assumptions are violated by the desired `P0`;
- reference density is not represented well by the final/protected basis;
- broad rejected directions would need to become MWG residual channels;
- implementation pressure points toward artifact/public workflow before the
  Hamiltonian convention is stable.

## Review Target

Before source work, reviewers should decide:

```text
Is the first source target a neutral fast atomic-reference Hartree GG block
builder from one-center atomic P0?
```

If yes, the source design should be written around a candidate
`HP-RHO0-MIXH-GG-FN-01` lane:

```text
P_A in same-center atomic GTO basis
  -> reference pair-density term stream
  -> Coulomb-expanded separable one-body factor packets
  -> terminal/base GG Hartree block
```

Dense Gaussian ERIs remain oracle/debug only. `GA` / `AA`,
protected-localized transforms, `F_app[P0]`, `C0`, artifacts, public workflow,
Cr/Cr2 diagnostics, and HF exchange remain later lanes.

## External ChatGPT Review Incorporation - 2026-07-03

The external review agreed with the fixed-`P0` direction and the retirement of
`J*w`, `diag(J)`, and `C' V C` as governing tests. It also tightened the first
implementation target. Accepted changes in this revision:

- first source slice is `GG`-only, not full `GG/GA/AA`;
- likely owner is a neutral raw-block file such as
  `src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl`;
- the pair-density term stream should be a small vector-backed internal record,
  not a runtime-size `NamedTuple` inventory;
- at least one angular/off-diagonal same-center pair must be checked against a
  dense oracle before trusting the term stream;
- `P_A` in the atomic AO/GTO basis and `P0_final` in the final basis are
  distinct objects and must not be mixed in traces;
- dense self-energy oracle is acceptable for the first source slice only as a
  bounded H/Be-size validation aid;
- `F_app[P0]` and `C0` are not part of the first exact Hartree source slice.

## Initial Design-Manager Review - 2026-07-03

Status: initial planning review. This review does not approve source edits,
new modules, public API, artifacts, solver workflow, committed tests, or Cr2
production work.

### Verdict

The plan is pointed in the right direction. The correct first source target is
not another row-gauge probe and not a dense final-basis ERI path. It is a
neutral exact Hartree mixed-block seam that takes one-center atomic reference
density matrices and produces exact one-body Hartree blocks in the existing
`GG` / `GA` / `AA` block language.

The most important design choice in the memo is sound:

```text
atomic P0 pair-density terms
  -> Coulomb-expanded separable one-body terms
  -> GG / GA / AA exact Hartree blocks
  -> existing protected/final one-body transforms
```

That keeps reference-density math out of Residual Gaussian ownership, preserves
the protected-localized injection convention, and avoids treating broad
rejected Gaussian directions as MWG residual channels.

### Strong Points

- The fixed-`P0` algebra is correctly separated from row-action diagnostics.
  The source path is designed around `Delta_F0` and `C0`, so the corrected
  model can match both energy and first derivative at the reference density.
- The one-center atomic-density assumption is a good first production-shaped
  simplification. It avoids molecular reference density products while still
  preserving cross-atom Coulomb interactions in the reference energy and
  potential.
- The memo correctly treats dense Gaussian Coulomb matrices as oracle/debug,
  not as the production path.
- The proposed `GG` / `GA` / `AA` output is the right interoperability seam:
  existing augmented/protected one-body transforms can consume those blocks.
- The validation ladder has the right ordering. H and Be/Be2 must settle the
  exact Hartree seam before Cr atom or Cr2 diagnostics make sense.
- The stop rules are mostly the right ones: no dense final four-index tensor,
  no row-action replacement for `F_exact[P0]`, no `C'VC`, and no public or
  artifact workflow before the Hamiltonian convention is stable.

### Decisions Needed Before Source Authority

1. **Split the first implementation into source slices.**

   Do not approve one broad `HP-RHO0-REFDENS-FN-01` that attempts all of the
   memo at once. The first source authority should likely be a small exact
   mixed Hartree kernel slice, with later slices for `GA`/`AA`, reference
   self-energy, approximate `F_app[P0]`, and correction assembly.

2. **Choose the neutral owner carefully.**

   The owner should be neutral exact Gaussian/Hartree infrastructure, not
   `cartesian_residual_gaussians`. Two plausible shapes are:

   ```text
   src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl
   ```

   or a small new neutral reference-density owner that depends on
   `gaussian_coulomb_reference.jl` and Cartesian raw-block machinery.

   The owner decision should be based on the lowest-level reusable object. If
   the first source slice only builds mixed Hartree raw blocks, it belongs near
   raw blocks. If it also owns `P0` construction and correction constants, that
   is a higher-level reference-density owner and should be a later slice.

3. **Define the reference pair-density term object before coding.**

   The plan should become more explicit about the term stream produced from
   `P_A[ab] chi_Aa chi_Ab`:

   - primitive coefficient convention;
   - symmetric off-diagonal handling;
   - center and angular-power representation;
   - compression key;
   - sign handling for density-matrix coefficients;
   - normalization audit against `Tr(P_A S_A)`.

   This object is likely the durable seam shared by mixed Hartree blocks and
   reference self-energy. It should not be a runtime-keyed `NamedTuple` or a
   large route-style metadata cloud.

4. **Keep `F_exact[P0]` and `F_app[P0]` as separate design seams.**

   The memo correctly says `F_app[P0]` must come from the actual solver
   convention. That is a separate missing seam risk. Do not let the first exact
   Hartree source pass invent an approximate-side evaluator. If `F_app[P0]` is
   not callable, record a separate authority request.

5. **Treat reference self-energy as an audit requirement, not necessarily the
   first fast kernel.**

   It is acceptable for the first source slice to use a dense Gaussian oracle
   for small `E_exact0_H` self-energy checks, provided the limitation is
   explicit. The fast self-energy path can share the reference pair-density
   stream later if needed.

6. **Make the block shape exact.**

   Before source approval, name whether the first source output is:

   ```text
   (; GG, GA, AA, reference_self_energy, diagnostics)
   ```

   or only:

   ```text
   (; GG, diagnostics)
   ```

   For the first slice, a `GG`-only source pass may be more reviewable if it
   can validate the term stream and terminal contraction against dense oracles.
   `GA`/`AA` can follow once the mixed operator is trusted.

7. **Preserve protected-localized ownership boundaries.**

   The mixed Hartree owner should not know about broad rejected directions,
   injection selection, MWG channels, or final solver policy. It should build
   exact one-body blocks. Residual/protected code can transform those blocks
   after the fact.

### Suggested First Source Authority Shape

If chat review agrees, the first future source lane should be narrower than the
whole memo:

```text
HP-RHO0-MIXH-GG-FN-01
```

Possible scope:

- neutral owner file only;
- construct same-center atomic reference pair-density term stream;
- integrate Coulomb-expanded reference side into separable one-body factor
  packets;
- build exact `GG` Hartree block for terminal/base `G`;
- compare H and small Be/Be2 `GG` subsets against dense Gaussian oracle;
- report symmetry, finite checks, term counts, compression counts, and
  reference electron-count normalization.

Explicitly out of that first source slice:

- `GA` / `AA`;
- final protected-localized transform;
- `F_app[P0]`;
- `C0` assembly;
- artifact/provenance;
- public driver/API;
- Cr/Cr2 diagnostics;
- HF exchange.

That slice would prove the hard numerical seam without carrying the whole
reference-density correction framework at once.

### Review Questions For Chat

- Is `GG`-first too narrow, or should the first source lane include `GA` and
  `AA` because protected-localized Be/Be2 needs the full block set to be
  meaningful?
- Should the neutral owner live under `cartesian_gaussian_raw_blocks`, or
  should reference-density get a new module from the beginning?
- Can the pair-density term stream reuse private
  `gaussian_coulomb_reference.jl` structures safely, or should it define a new
  compact public-to-the-package internal record?
- Is dense oracle self-energy acceptable for the first source slice, or must
  reference self-energy share the compressed term stream immediately?
- What is the smallest validation that catches angular/off-diagonal pair-term
  mistakes without committing slow development-era tests?

### Bottom Line

The plan should be kept. It correctly identifies the missing source seam and
prevents the main category mistakes: row-gauge substitution, `C'VC` revival,
dense final ERIs, and RG-owned reference-density math.

Before implementation, tighten the first source lane around one reusable exact
mixed Hartree block builder. Keep the first approved source slice small enough
that a doer can prove the pair-density term stream and `GG` contraction before
adding the rest of the correction machinery.

## Design-Manager Source Authority - 2026-07-03

The revised plan is accepted as a review basis, and the first source lane is
approved only for the `GG` mixed Hartree seam. This approval deliberately does
not promote the whole reference-density correction framework.

Approved IDs:

- `HP-RHO0-MIXH-GG-FN-01`
- `HP-RHO0-MIXH-GG-TEST-01`

Approved source owner:

- `src/cartesian_gaussian_raw_blocks/mixed_hartree_blocks.jl`
- `src/cartesian_gaussian_raw_blocks/CartesianGaussianRawBlocks.jl` only for
  internal include/import wiring;
- optional narrow reuse or extraction from `src/gaussian_coulomb_reference.jl`
  only if needed to avoid duplicating existing Gaussian Coulomb oracle or
  pair-term algebra.

Approved behavior:

- consume one-center atomic reference density matrices `P_A`;
- build a vector-backed same-center reference pair-density term stream,
  including angular and off-diagonal same-center pairs;
- validate finite/symmetric `P_A`, dimensions, and electron-count
  normalization diagnostics;
- use the existing Coulomb Gaussian expansion convention to build separable
  one-body factor packets;
- build only the terminal/base exact Hartree `GG` block;
- return compact internal diagnostics such as pair-term counts, expansion
  packet counts, reference electron count, symmetry, and dense-oracle deltas;
- keep dense Gaussian Coulomb matrices as bounded oracle/debug validation only.

Explicitly out of scope:

- `GA` and `AA`;
- protected-localized, injected, residual, or final-basis transforms;
- `F_app[P0]`, `Delta_F0`, `C0`, reference self-energy production, or
  corrected Hamiltonian assembly;
- artifact/provenance/schema/writer/reader/manifest changes;
- public driver/API/export/defaults or solver workflow;
- Cr atom, Cr2, Cr2 HF, or Cr2 production diagnostics;
- HF exchange and `(final GTO | GTO final)` kernels;
- row-action, `diag(J)`, `q0`, center metadata, direct `C' V C`, or IDA proxy
  shortcuts;
- residual/MWG default changes, basis-fate changes, broad rejected directions
  as MWG residuals, dense final four-index ERIs, runtime-keyed construction
  inventories, committed tests, or fixtures.

Validation required:

- `git diff --check`;
- package load;
- ignored H direct-core or one-center terminal `GG` replay;
- bounded one-center Be or Be2 terminal/base `GG` smoke with finite/symmetric
  output;
- dense Gaussian Coulomb oracle spot checks for small `GG` subsets;
- at least one angular/off-diagonal same-center reference-pair oracle check;
- report pair-term counts, Coulomb expansion packet counts, reference
  electron-count normalization, dense-oracle deltas, and whether dense paths
  were validation-only;
- no `GA`/`AA`, no protected transform, no artifacts, no public workflow, and
  no Cr/Cr2 run.

Failure rule: if the source pass needs dense final four-index ERIs, files
outside the approved surface, `GA`/`AA`, protected transforms, artifact/schema
work, or broader reference-density correction machinery, stop and report the
exact missing owner or kernel. Do not widen the lane in source.
