# Reference-Density Hartree Screening

Reference-density Hartree screening rewrites the direct electron-electron
energy around a declared fixed reference density. The public interface accepts
an already evaluated reference Hartree field in the same orthonormal
one-particle basis as the interaction and reference orbitals. It returns a
one-body correction and a separate scalar; it does not alter the physical
kinetic-plus-nuclear one-body operator.

## Supplied Reference Fields

Use one of two typed field routes:

```julia
exact_field = ExactRepresentedHartreeField(
    J0,
    reference_coulomb_self_integral;
    provenance = "declared represented-density evaluator",
)

fitted_field = FittedReferenceHartreeField(
    J0_fit,
    density_coulomb_self_integral;
    provenance = "declared fitted-potential evaluator",
)
```

The exact scalar is

```text
reference_coulomb_self_integral = (rho0|rho0) = Tr(P0 * J0),
```

which is twice the reference Hartree energy. A fitted potential supplies the
self integral of its underlying density representation separately; the
potential fit does not redefine that scalar. The distinct types prevent a
fitted field from silently entering the exact represented route.

The public API does not construct atomic density fits, evaluate fitted
potentials on a general terminal-plus-supplement basis, place atomic packets,
or form additive molecular references. Those remain internal facilities.

## Correction Assembly

Assemble the correction with:

```julia
correction = screened_hartree_correction(
    V_IDA,
    exact_field,
    reference_coefficients,
    occupations,
)
```

`V_IDA`, the supplied field matrix, and `reference_coefficients` must use the
same orthonormal basis and ordering. The occupations are finite, nonnegative,
spin-summed occupations and may be fractional. They define

```text
P0 = reference_coefficients * Diagonal(occupations) *
     reference_coefficients'
q0 = diag(P0).
```

No neutrality condition is imposed. `representation_atol` controls orbital,
trace, and occupation agreement. `density_nonnegativity_atol` only allows
small negative roundoff in `q0`.

The correction is

```text
Delta_J0 = J0 - Diagonal(V_IDA * q0)
C        = 0.5 * q0' * V_IDA * q0 - 0.5 * (rho0|rho0).
```

Read it through:

```julia
screened_hartree_delta_one_body(correction)
screened_hartree_energy_constant(correction)
screened_hartree_consistency_error(correction)
screened_hartree_field_kind(correction)
```

The signed consistency convention is

```text
Tr(P0 * J0) - (rho0|rho0).
```

The exact route enforces closure. The fitted route reports a finite nonzero
consistency error rather than rejecting the field solely for that
approximation.

## Example And Boundary

[`examples/40_screened_hartree_fixed_density.jl`](https://github.com/srwhite59/GaussletBases.jl/blob/main/examples/40_screened_hartree_fixed_density.jl)
builds a two-function `CartesianIDAHamiltonian`, supplies a bounded analytic
Gaussian reference field, and checks energy and occupied-action closure. It is
an assembly example, not reproduction of historical He.

This surface corrects only the direct Hartree accounting. It does not construct
an exchange field, run SCF, expose packet/fitting machinery, or provide the
paper-local Ximg/XHF finite-reference corrections.
