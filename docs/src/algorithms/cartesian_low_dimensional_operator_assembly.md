# Cartesian Low-Dimensional Operator Assembly

This page defines how Cartesian product operators are assembled from
one-dimensional data. It is the low-dimensional operator contract used by the
PQS shell and IDA Hamiltonian pages.

## Spaces and Dimensions

For each Cartesian axis `a in (x, y, z)`, the one-dimensional data include
overlap, kinetic, position, second-moment, and Gaussian-factor tables. A
three-dimensional product index is ordered consistently with the route's
support-row convention. If the axis dimensions are `n_x`, `n_y`, and `n_z`, the
product support dimension is `n_x n_y n_z`.

## Inputs

- axis overlap matrices `S_x`, `S_y`, `S_z`;
- axis kinetic matrices `K_x`, `K_y`, `K_z`;
- optional position and second-moment matrices;
- Gaussian-expansion coefficients and exponents for `1/r`;
- per-center locations for nuclear attraction;
- pair-factor data for IDA electron-electron assembly.

## Outputs

- Cartesian overlap `S_3D`;
- Cartesian kinetic `K_3D`;
- separated, uncharged nuclear attraction matrices `U_A`;
- raw and density-normalized pair factors used to build two-index `Vee`;
- moment data used by residual-Gaussian MWG approximations.

## Pseudocode

1. Build or read the one-dimensional axis matrices and factor tables.

2. Assemble product overlap:

   ```math
   S_{3D} = S_x \otimes S_y \otimes S_z.
   ```

3. Assemble product kinetic:

   ```math
   K_{3D}
   =
   K_x\otimes S_y\otimes S_z
   +
   S_x\otimes K_y\otimes S_z
   +
   S_x\otimes S_y\otimes K_z.
   ```

4. For a center `A`, use the Gaussian expansion of `1/r` to assemble the
   uncharged attraction matrix:

   ```math
   U_A
   =
   -\sum_\mu c_\mu
   F^x_{A\mu}\otimes
   F^y_{A\mu}\otimes
   F^z_{A\mu}.
   ```

5. Build raw electron-electron pair numerators from product pair factors and
   support weights.

6. Convert raw pair numerators to density-normalized pair factors only at
   explicit boundaries that require that convention. Do not mix a
   density-normalized factor with a raw weighted contraction.

7. Contract the product data through the localized basis coefficients to obtain
   `K`, separated `{U_A}`, and `Vee` in the localized IDA basis.

## Linear Algebra

The Cartesian assembly is a product-form contraction. The efficient path keeps
the short Gaussian-expansion sum and axis factorizations visible for as long as
possible, then applies small retained/final coefficient contractions only after
the source or support operator is built.

## Allowed Orthogonalizations

None. This page assembles operators and factors. Orthogonalization belongs only
to PQS shell construction and residual-Gaussian residualization.

## Forbidden Operations

- Do not perform PQS shell projection or Lowdin cleanup inside raw product
  operator assembly.
- Do not sum nuclear centers into one charged matrix when the caller needs
  counterpoise.
- Do not identify raw pair numerators with density-normalized pair factors.
- Do not build four-index ERI tensors as the public Cartesian IDA route.

## Numerical Invariants

- `S_3D`, `K_3D`, each `U_A`, and `Vee` are symmetric up to numerical
  tolerance.
- Gaussian expansion terms use the same coefficient/exponent convention across
  nuclear attraction and electron-electron pair factors.
- Raw weighted and density-normalized pair factors are not mixed in one
  contraction without the explicit support-weight conversion.

## Operator and Gauge Conventions

`U_A` is always uncharged: it represents `-1/r_A`, not `-Z_A/r_A`.
Charges are applied later by the Hamiltonian assembly:

```math
H_1 = K + \sum_A Z_A U_A.
```

`Vee` is the two-index IDA interaction in the localized IDA basis. IDA is a
model choice, not a failed attempt to store full four-index ERIs.

## Code Map

- `src/cartesian_pair_block_materialization/pqs_source_safe_terms.jl`
  contains source-space safe one-body term assembly.
- `src/pqs_multilayer_complete_core_shell_h1.jl` consumes support operators
  for the common complete core/shell H1 path.
- `src/cartesian_gaussian_axis_integrals.jl` owns shared Gaussian axis integral
  kernels where present.
- `src/ordinary_qw_raw_blocks.jl` still contains donor/reference raw-block
  surfaces used by surviving QW routes.

## Current Implementation Deviations

Some older donor routes still expose QW/WL-specific wrapper names around shared
Gaussian and weighted-Hadamard kernels. Those wrappers are migration surfaces,
not the long-term public Cartesian operator contract.
