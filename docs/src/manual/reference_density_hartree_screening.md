# Reference-Density Hartree Screening

Status: approved for the proposed GaussletBases v0.2 public surface;
implementation is pending on the current development branch.

Reference-density Hartree screening rewrites the direct electron-electron
energy around a declared fixed reference density. The approved public surface
will assemble a one-body correction and a separate scalar from a supplied
same-basis reference Hartree field. It will not construct atomic density fits,
evaluate fitted potentials on general mixed bases, or run a solver.

The complete argument, sign, scalar, charge, exact-versus-fitted, validation,
and compatibility contract is tracked in the developer design
[PQS v0.2 public surface](../developer/designs/cartesian_hamiltonian_producer/pqs_public_surface.md).

The approved names are not yet implemented or supported on the current
package version. This page must not be used as executable API documentation
until the implementation lifecycle is accepted.

Paper-local Ximg/XHF finite-reference corrections are outside this public
surface. They are provenance for particular paper states, not a general
mixed-basis interaction, exchange correction, or correlation method.
