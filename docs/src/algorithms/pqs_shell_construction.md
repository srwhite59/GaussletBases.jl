# PQS Shell Construction

This page is the normative algorithm for localized PQS shell construction. It
defines where Lowdin orthogonalization is allowed and, more importantly, where
it is forbidden.

## Spaces and Dimensions

Let `B_<s` be the matrix of all already-retained inner/core and previous-shell
functions for shell step `s`. Let `X_s` be the candidate source-box functions
for shell `s`, represented in the same support metric `S`.

If shell `s` retains `r_s` functions, the only Lowdin matrix for that shell is
`r_s x r_s` after projection onto the shell complement.

## Inputs

- support overlap metric `S`;
- already-retained localized functions `B_<s`;
- shell candidate functions `X_s`;
- retained shell count and keep policy;
- tolerances for projection rank and shell overlap checks.

## Outputs

- shell-local realized basis block `B_s`;
- shell-local Lowdin transform `L_s`;
- shell Gram eigenvalues and rank diagnostics;
- appended localized basis `B_<=s = [B_<s, B_s]`.

## Pseudocode

The source-box shell stage order is:

```text
raw product-box source modes
-> boundary-mode selection
-> box-to-shell projection
-> shell-local Gram matrix
-> shell-local Lowdin
-> append
```

1. Start with accepted inner functions `B_<s`.

2. Form the projector onto the already-retained span:

   ```math
   P_{<s}
   =
   B_{<s}
   (B_{<s}^{T} S B_{<s})^{-1}
   B_{<s}^{T} S.
   ```

3. Project the shell candidates out of that span:

   ```math
   \widetilde X_s = (I - P_{<s}) X_s.
   ```

4. Build the shell-local Gram matrix:

   ```math
   G_s = \widetilde X_s^T S \widetilde X_s.
   ```

5. Apply a Lowdin transform only to this shell-local Gram matrix:

   ```math
   L_s = G_s^{-1/2}.
   ```

6. Realize the shell block:

   ```math
   B_s = \widetilde X_s L_s.
   ```

7. Append the block without rotating any earlier block:

   ```math
   B_{\leq s} = [B_{<s}, B_s].
   ```

8. Check, but do not globally repair, the full localized basis overlap.

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

The projection step is responsible for shell-to-inner orthogonality. The
shell-local Lowdin step is responsible only for orthonormalizing retained
directions within the new projected shell block.

If `B_<s^T S B_s` is not small, the projection step is wrong or insufficient.
The remedy is to repair projection, not to rotate `B_<s` and `B_s` together.

## Allowed Orthogonalizations

- Shell-local Lowdin on `G_s`.
- Optional rank truncation inside `G_s` according to the shell keep policy.

## Forbidden Operations

- No Lowdin over all core and shell functions together.
- No rotation of previously completed core or shell blocks.
- No global orthogonalization to conceal shell-projection error.
- No signed-final-weight division for the localized IDA interaction.
- No promotion of a source-backed WL/QW retained transform as independent PQS
  shell authority.

## Numerical Invariants

Each completed shell must satisfy:

```math
B_s^T S B_s \approx I,
\qquad
B_{<s}^T S B_s \approx 0.
```

The concatenated localized PQS basis should satisfy:

```math
B_{\leq s}^T S B_{\leq s} \approx I.
```

Failure of the concatenated check is a construction error. It must not be
hidden by a global Lowdin cleanup.

## Operator and Gauge Conventions

The localized PQS basis produced by this algorithm is the IDA working basis.
The IDA pair matrix and the one-body matrices should be constructed in this
same localized basis for the base all-electron Hamiltonian.

## Code Map

- `src/cartesian_final_basis_realization/pqs_source_shell_final_basis.jl`
  implements shell projection, shell-local Lowdin cleanup, and shell block
  realization.
- `src/cartesian_shellification/terminal_geometry.jl` implements the
  diatomic atom-contact core seed-box hull rule and terminal region coverage
  checks before lowering.
- `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`
  concatenates core and shell sectors and checks the completed overlap.
- `src/pqs_multilayer_complete_core_shell_h1.jl` consumes the completed basis
  for common H1 and density-interaction construction.

## Current Implementation Deviations

The active H2 PQS route has removed the forbidden combined core/shell Lowdin.
The private H2 residual-GTO route now uses an internal one-basis IDA object.
The public `CartesianIDAHamiltonian` type and minimal writer/reader now exist,
and the shell construction itself must remain source-box-first and shell-local.
