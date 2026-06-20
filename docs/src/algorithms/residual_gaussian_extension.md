# Residual-Gaussian Extension

This page defines the residual-Gaussian extension for a fixed localized
gausslet/PQS basis. It is the normative rule for adding Gaussian supplement
orbitals without changing the already-built gausslet/PQS functions.

## Spaces and Dimensions

Let `G` be the fixed localized gausslet/PQS basis, orthonormal in metric `S`.
Let `A` be the matrix of Gaussian supplement candidates. If `k` residual
directions survive the cutoff, the residual block `R` has `k` columns.

## Inputs

- fixed localized basis `G`;
- overlap metric `S`;
- Gaussian supplement candidates `A`;
- supplement self-overlap `A^T S A`;
- cross overlap `G^T S A`;
- residual eigenvalue cutoff.

## Outputs

- residual Gaussian basis block `R`;
- residual transform `L_R`;
- retained residual eigenvalues and rank;
- overlap diagnostics for `G`, `R`, and their cross block.

## Pseudocode

1. Keep the localized gausslet/PQS basis fixed.

2. Project Gaussian candidates out of the fixed basis. If `G^T S G = I`:

   ```math
   A_\perp = A - G (G^T S A).
   ```

   In the general form:

   ```math
   A_\perp =
   \left[I - G(G^T S G)^{-1}G^T S\right]A.
   ```

3. Form only the small residual overlap:

   ```math
   S_R = A_\perp^T S A_\perp.
   ```

   Equivalently, in the orthonormal fixed-basis case:

   ```math
   S_R = A^T S A - (G^T S A)^T(G^T S A).
   ```

4. Diagonalize the residual overlap:

   ```math
   S_R = U \Lambda U^T.
   ```

5. Keep eigenvectors with eigenvalues above cutoff and define:

   ```math
   L_R = U_{\rm keep}\Lambda_{\rm keep}^{-1/2},
   \qquad
   R = A_\perp L_R.
   ```

6. Append `R` to the fixed basis without reorthogonalizing `[G, R]`.

## Linear Algebra

Residualization is one-sided. The fixed gausslet/PQS basis is not changed. The
only Lowdin transform is the small residual-block transform `L_R`.

For an implementation that carries residual coefficients in the raw carrier
`[G, A]`, the residual carrier is:

```math
\begin{bmatrix}
-X \\
I
\end{bmatrix}
L_R,
\qquad
X = G^T S A.
```

There is no multiplication by a global PQS cleanup matrix.

## Allowed Orthogonalizations

- Lowdin of the residual overlap `S_R` after projection out of the fixed
  gausslet/PQS basis.
- Eigenvalue truncation inside the residual Gaussian block.

## Forbidden Operations

- No Lowdin of `[G, R]`.
- No rotation of the fixed gausslet/PQS basis.
- No use of a residual Gaussian as if it were an unprojected supplement
  Gaussian.
- No construction of residual two-electron data by blindly applying the
  one-body residual transform to unrelated pair data.

## Numerical Invariants

The completed residual extension must satisfy:

```math
G^T S G \approx I,
\qquad
G^T S R \approx 0,
\qquad
R^T S R \approx I.
```

The retained residual eigenvalues should be positive and above the selected
cutoff.

## Operator and Gauge Conventions

One-body matrices involving residual Gaussians can be built exactly from raw
Gaussian/PQS blocks and then transformed with `L_R`.

Two-electron terms involving residual Gaussians remain in the IDA model. The
current route uses matched-width-Gaussian (MWG) residual descriptors for
residual-containing IDA interaction terms; the Qiu-White reference route also
documents GGT and MWG variants.

## Code Map

- `src/ordinary_qw_residuals.jl` contains the Qiu-White donor/reference
  residual-space construction.
- `src/ordinary_qw_operator_assembly.jl` contains donor/reference GGT and MWG
  residual interaction logic.
- [Qiu-White residual-Gaussian route](qiu_white_residual_gaussian_route.md)
  remains the paper-faithful reference description for the old QW line.

## Current Implementation Deviations

The former private H2 residual-GTO handoff helper has been retired. A public
residual-Gaussian producer still needs to be rebuilt on the generic terminal
PQS route before this algorithm is available as a promoted construction path.
