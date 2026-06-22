# PQS Shell Construction

This page is the normative algorithm for localized PQS shell construction. It
defines where Lowdin orthogonalization is allowed and, more importantly, where
it is forbidden.

## Spaces and Dimensions

Let `X_s` be the boundary-selected product-mode columns generated from the full
source box for shell step `s`. Let `I_s` be the authoritative owned terminal
support rows for that shell, represented by `support.support_indices` and
`support.support_states`.

The full source box is used only to generate candidate product-mode columns.
The realized shell block lives only on rows `I_s`. If shell `s` retains `r_s`
functions, the only Lowdin matrix for that shell is `r_s x r_s` and is formed
from the shell-local Gram matrix on rows `I_s`.

## Inputs

- full source-box product modes;
- boundary product-mode column selection for shell `s`;
- authoritative owned shell support rows `I_s`;
- retained shell count and keep policy;
- tolerances for shell-local rank and overlap checks.

## Outputs

- shell-local realized basis block `B_s`;
- shell-local Lowdin transform `L_s`;
- shell Gram eigenvalues and rank diagnostics;
- appended terminal block record with unchanged owned support rows `I_s`.

## Pseudocode

The source-box shell stage order is:

```text
full source-box product modes
-> boundary product-mode column selection
-> restrict rows to authoritative owned shell support
-> shell-local Gram matrix
-> shell-local Lowdin
-> final sign canonicalization
-> append block with unchanged owned support
```

1. Generate full source-box product modes.

2. Select boundary product-mode columns for shell `s`.

3. Restrict candidate rows to authoritative owned shell support `I_s`:

   ```math
   X_s^{I} = R_{I_s} X_s.
   ```

4. Build the shell-local Gram matrix on the owned support:

   ```math
   G_s = (X_s^{I})^T S_{I_s I_s} X_s^{I}.
   ```

5. Apply a Lowdin transform only to this shell-local Gram matrix:

   ```math
   L_s = G_s^{-1/2}.
   ```

6. Realize the shell block:

   ```math
   B_s = X_s^{I} L_s.
   ```

7. Canonicalize final signs according to the current PQS final-weight gauge.

8. Append a terminal block whose support rows remain exactly `I_s`.

## Diatomic Atom-Contact Core Rule

For a bond-aligned diatomic, first define one odd `q`-side atom seed box around
each snapped nuclear grid index. These seed boxes describe direct gausslet
support near each nucleus; they are not shell-projected PQS sectors.

If the two seed boxes overlap, touch, or have a bond-axis gap shorter than `q`,
the direct molecular core is the discrete hull of the two seed boxes:

```text
atom_contact_core = hull(left_q_seed_box, right_q_seed_box)
```

This hull rule is the entire size rule. Do not force a double-core volume and
do not force the bond-axis length to be odd. For `q = 5`, coincident nuclei
therefore give a `5 x 5 x 5` direct core, a small sub-core nuclear separation
gives a slightly elongated direct core such as `5 x 5 x 6`, and the current
H2 q5 fixture gives `5 x 5 x 11`.

When the seed-box gap is at least `q`, keep the seed boxes as separate
atom-local direct cores, grow atom-local projected shells outside them, and
represent the central region by midpoint slabs or a distorted product box
according to the terminal shellification geometry. In all cases, separated
unit nuclear attraction matrices remain center-specific even when both nuclei
lie inside the same direct atom-contact core.

## Linear Algebra

Parent gausslet rows are orthonormal to machine precision. Terminal regions own
disjoint parent rows. Therefore block-local terminal basis supports are
structurally orthogonal across blocks.

The shell-local Lowdin step is responsible only for orthonormalizing retained
directions within the new shell-owned support. There is no projection against
previous terminal blocks and no corrective projection step.

If a cross-block overlap is nonzero after owned-support restriction, the
problem is duplicated support rows, incorrect row restriction, wrong support
ownership, or an indexing error. It is not a physical residual and must not be
repaired by mixing coefficients into previous supports.

## Allowed Orthogonalizations

- Shell-local Lowdin on `G_s`.
- Optional rank truncation inside `G_s` according to the shell keep policy.

## Forbidden Operations

- No previous-block projection or recursive projection.
- No Lowdin over all core and shell functions together.
- No rotation of previously completed core or shell blocks.
- No global orthogonalization to conceal support ownership or indexing errors.
- No growth of a shell block onto previous terminal regions.
- No signed-final-weight division for the localized IDA interaction.
- No promotion of a source-backed WL/QW retained transform as independent PQS
  shell authority.

## Numerical Invariants

Each completed shell must satisfy:

```math
B_s^T S_{I_s I_s} B_s \approx I.
```

The terminal realization must also satisfy:

- every block support equals its authoritative terminal support;
- terminal support sets are pairwise disjoint;
- structural cross-block overlap is zero by construction.

Failure of any structural check is a construction error. It must not be hidden
by a global Lowdin cleanup or previous-block projection.

## Operator and Gauge Conventions

The localized PQS basis produced by this algorithm is the IDA working basis.
The IDA pair matrix and the one-body matrices should be constructed in this
same localized basis for the base all-electron Hamiltonian.

Structural cross-block overlap zero does not make operators block diagonal.
Cross-block kinetic, nuclear-attraction, and localized IDA interactions may be
nonzero and remain assembled over terminal block pairs.

## Code Map

- `src/cartesian_final_basis_realization/pqs_terminal_basis_realization.jl`
  owns the current terminal block-local realization contract.
- `src/cartesian_final_basis_realization/pqs_source_shell_final_basis.jl`
  is a donor/reference path only where it still preserves older shell
  projection vocabulary; it is not current terminal support authority.
- `src/cartesian_shellification/terminal_geometry.jl` implements the
  diatomic atom-contact core seed-box hull rule and terminal region coverage
  checks before lowering.
- `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`
  is oracle/reference for legacy complete core/shell basis assembly, not the
  current block-local terminal support authority.
- `src/pqs_multilayer_complete_core_shell_h1.jl` consumes the completed basis
  in legacy/operator-reference paths.

## Current Implementation Deviations

The active H2 PQS route has removed the forbidden combined core/shell Lowdin.
The private H2 residual-GTO route now uses an internal one-basis IDA object.
The public `CartesianIDAHamiltonian` type and minimal writer/reader now exist,
and the shell construction itself must remain source-box-first, owned-support
restricted, and shell-local.

Production source still needs cleanup to replace cross-overlap audit /
`max_cross_overlap` plumbing with structural support checks. That cleanup is a
separate implementation handoff, not approved by this algorithm note.
