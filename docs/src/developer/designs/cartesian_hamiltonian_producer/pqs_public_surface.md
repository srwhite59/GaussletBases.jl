# PQS v0.2 Public Surface

Status: implemented and validated under:

- `HP-PQS-PUBLIC-DOC-01`;
- `HP-PQS-PUBLIC-MATCHED-FN-01` and
  `HP-PQS-PUBLIC-MATCHED-TEST-01`;
- `HP-PQS-PUBLIC-SCREEN-FN-01` and
  `HP-PQS-PUBLIC-SCREEN-TEST-01`.

Implementation commit: `058ee54f45c759949f70b54a699ccc318476f8ac`.

This page owns the narrow supported interface intended for GaussletBases
v0.2.0. It does not authorize a release tag. The existing numerical
contracts remain in [R1 public base producer](r1_public_base_producer.md),
[matched PQS/White-Lindsey aspect modes](pqs_complete_shell_aspect_source_modes.md),
and [screened-Hartree correction assembly](screened_hartree_correction_assembly.md).

## Purpose And Boundary

The public surface has two independent parts:

1. one bounded H2+ comparison that reconstructs the common full parent and
   matched PQS/White-Lindsey terminal spaces used for Table I; and
2. correction assembly from a consumer-supplied same-basis reference Hartree
   field and represented fixed density.

Neither part exposes the staged producer, atomic packet construction, general
mixed-basis field evaluation, a solver, or paper-local state exactification.
The existing `cartesian_base_hamiltonian` remains the general supported
PQS/White-Lindsey constructor. `cartesian_base_working_basis` remains an
expert/unstable construction surface whose result fields are not public
schema.

## Matched H2+ Public Comparison

The approved exported names are:

```julia
PQSH2PlusRow
PQSH2PlusComparison
pqs_h2plus_comparison
```

The single constructor is:

```julia
pqs_h2plus_comparison(;
    independent_reference_total_energy_hartree::Real,
)::PQSH2PlusComparison
```

The independent reference is required and is used only to form the reported
total errors. The Table I example supplies
`-0.6026342144949465` hartree. The producer does not relabel that external
reference as a repository-computed value.

`PQSH2PlusRow` is a small immutable result row with stable fields:

```text
basis
dimension
capture
electronic_energy_hartree
nuclear_repulsion_energy_hartree
total_energy_hartree
contraction_error_hartree
total_error_hartree
```

The basis values and row order are fixed as `:parent`, `:pqs`, and
`:white_lindsey`. Parent contraction error is zero by definition.

`PQSH2PlusComparison` contains the fixed three rows, the declared independent
reference, and compact common-construction facts:

```text
rows
independent_reference_total_energy_hartree
parent_axis_counts
parent_dimension
direct_core_columns
complete_shell_count
complete_shell_columns
slab_count
slab_columns
```

No internal parent, shell, block, transform, coefficient, due-diligence row,
or fingerprint object is returned. The constructor must fail before returning
unless its live objects establish one common physical construction:

- H2+ at `R = 2` bohr;
- public `ns = 5`, `core_spacing = 0.30` bohr, `s_factor = 1`;
- padding `10` bohr, tail spacing `2.8` bohr, ordinary source span, and
  high135 Coulomb representation;
- route-local PQS `q = 5` and White-Lindsey `q = 3` on the same physical
  parent;
- parent axes `21 x 21 x 29`, dimension `12789`;
- identical parent centers, weights, mapped-axis operators, physical bounds,
  and same-run numerical fingerprints;
- identical physical region keys, shell ownership/supports, direct core, and
  slabs;
- eight complete shells, two slabs, and equal aggregate column accounting
  `275 + 960 + 50 = 1285` for both terminal methods; and
- the current axis-specific White-Lindsey inner construction derived from the
  shared `(ns,ns,L)` shape, never the superseded scalar shell policy.

The full parent solve remains matrix-free. The implementation may call the
accepted staged and apply-based numerical owners internally, but those names
do not become public. It constructs no `Vee`, artifact, or solver state.

The supported example writes a TSV with exactly the stable row fields above.
Additional timing, allocation, Git, path, shell-ledger, or report fields from
the private paper driver are not part of this API.

## Screened-Hartree Public Assembly

The approved exported names are:

```julia
ExactRepresentedHartreeField
FittedReferenceHartreeField
ScreenedHartreeCorrection
screened_hartree_correction
screened_hartree_delta_one_body
screened_hartree_energy_constant
screened_hartree_consistency_error
screened_hartree_field_kind
```

The field constructors are:

```julia
ExactRepresentedHartreeField(
    reference_hartree_field,
    reference_coulomb_self_integral;
    provenance,
)

FittedReferenceHartreeField(
    fitted_reference_hartree_field,
    density_coulomb_self_integral;
    provenance,
)
```

Both matrices are copied finite symmetric matrices. Provenance is a required,
nonempty descriptive string. The two concrete types are intentionally
distinct: a fitted potential field cannot silently enter the exact represented
route.

The public assembly entry point is:

```julia
screened_hartree_correction(
    V_IDA,
    field,
    reference_coefficients,
    occupations;
    representation_atol = 1.0e-8,
    density_nonnegativity_atol = 1.0e-12,
    symmetry_atol = 1.0e-10,
    closure_atol = 1.0e-8,
)::ScreenedHartreeCorrection
```

`V_IDA`, the reference-field matrix, and the rows of
`reference_coefficients` all use the same orthonormal one-particle basis and
ordering. The columns of `reference_coefficients` are represented reference
orbitals. `occupations` are finite, nonnegative, spin-summed occupations;
fractional values are allowed.

Define the spin-summed reference density and localized charges by:

```text
P0 = reference_coefficients * Diagonal(occupations) *
     reference_coefficients'
q0 = diag(P0)
```

`representation_atol` owns orbital orthonormality, trace, and occupation-sum
agreement. `density_nonnegativity_atol` is only the numerical allowance for a
small negative entry of `q0`; it is not a neutrality or charge-state test.
The reference particle number is `sum(occupations)`. No relationship to the
sum of nuclear charges is required.

Public documentation names the exact constructor's scalar input
`reference_coulomb_self_integral` and defines it as:

```text
(rho0|rho0) = Tr(P0 * J0)
```

It is twice the reference Hartree energy. The fitted constructor's
`density_coulomb_self_integral` is the corresponding density-representation
quantity; the fitted potential does not redefine it.

The returned correction is:

```text
Delta_J0 = J0 - Diagonal(V_IDA * q0)
C        = 0.5 * q0' * V_IDA * q0 - 0.5 * (rho0|rho0)
```

`Delta_J0` and `C` are added to direct electron-electron/Hartree energy
accounting. They do not replace or mutate the physical kinetic-plus-nuclear
one-body operator. The scalar remains a separate total-energy constant. A
consumer may represent it as `(C/N)I` only inside a fixed nonzero
particle-number sector `N`; that optional rewrite is not provided by v0.2 and
is not valid across particle sectors.

The stable accessors return the correction matrix, scalar, signed consistency
error, and field kind. The signed convention is always:

```text
screened_hartree_consistency_error =
    Tr(P0 * J0) - (rho0|rho0)
```

For `ExactRepresentedHartreeField`, its magnitude must not exceed
`closure_atol`. Exact means exact for the declared represented density and
Coulomb representation, not continuum exactness. For
`FittedReferenceHartreeField`, the same signed value is a reported potential-
fit approximation and is not by itself a rejection condition. Finiteness,
symmetry, dimensions, represented-density validity, and the derivative/action
identity remain hard failures for both routes.

The result type is supported through the four accessors above. Its internal
field layout, packet summary, and diagnostic `NamedTuple` are not public
schema and receive no serialization promise.

This v0.2 surface assembles a correction from a supplied same-basis field. It
does not publicly construct atomic density fits, fit a potential, evaluate a
fitted potential on a terminal-plus-supplement basis, place atomic packets, or
form additive molecular self/cross energies. Those capabilities remain
internal under their existing contracts.

No exchange field or exchange correction is returned. The words "exact" and
"fitted" classify the supplied direct Hartree field only.

## User Documentation And Layering

The user documentation must distinguish:

1. stable PQS/White-Lindsey construction through
   `cartesian_base_hamiltonian` and the bounded matched H2+ comparison;
2. stable correction assembly from a supplied same-basis reference Hartree
   field; and
3. paper-local Ximg/XHF finite-reference exactification.

Ximg/XHF remains provenance-only. It is neither reconstruction of a physical
mixed-basis four-index interaction nor transferable exchange/correlation
theory. No Ximg/XHF name, constructor, result, or accessor is exported.

The implemented user surfaces are:

```text
docs/src/manual/projected_q_shells.md
docs/src/manual/reference_density_hartree_screening.md
docs/src/manual/index.md
docs/src/examples/index.md
docs/src/howto/example_guide.md
examples/40_screened_hartree_fixed_density.jl
examples/41_pqs_h2plus_table1.jl
```

The screening manual page and both examples are now source-backed and may be
rendered as executable public documentation. A candidate release still
requires one clean replay of both examples; this lifecycle acceptance does not
itself create a tag.

The screening example uses two bounded pure-GTO orbitals and an explicitly
declared one-center, two-electron fixed density. It obtains the direct IDA
matrix from a small `CartesianIDAHamiltonian`, keeps
`one_body_hamiltonian(ham)` separate and unchanged, and obtains the accurate
represented reference field from the existing bounded pure-GTO Coulomb
oracle. It demonstrates energy and occupied-action closure with a nontrivial
correction. It is an executable assembly example, not reproduction of the
historical He calculation or a new He endpoint.

## Validation Contract

Ordinary pull-request coverage uses small dimensions and requires:

- exported-name and constructor smoke coverage;
- exact structural dimensions and spin-summed reference charge;
- matrix symmetry at `1e-10`;
- energy and occupied-action closure at `1e-11`;
- a nonzero `Delta_J0` in the two-function fixture;
- the signed exact consistency error within `1e-11`;
- a fitted-field control that reports rather than rejects a finite nonzero
  signed consistency error;
- unchanged physical one-body data; and
- malformed dimensions, asymmetry, nonfinite data, negative occupations,
  representation loss, and exact-route closure failure rejected with
  `ArgumentError` or `DimensionMismatch` as appropriate.

The existing small `examples/39_pqs_h2plus.jl` remains the ordinary PQS/WL
construction smoke. The new fixed-density screening example may join the
existing quick-example set if its measured runtime remains bounded.

A separate release/nightly test owns the full Table I comparison. It requires:

- exact parent/PQS/White-Lindsey dimensions and structural counts;
- exact same-run parent and physical-shell identity between PQS and
  White-Lindsey;
- capture absolute tolerance `2e-8`, relative tolerance zero;
- electronic, total, contraction-error, and total-error absolute tolerance
  `1e-7` hartree, relative tolerance zero;
- same-run parent residual and capture closure at `1e-9`;
- same-run symmetry at `1e-10`; and
- no cross-platform raw-byte or matrix-fingerprint equality requirement.

The slower test is not part of every local edit gate. It is run by a focused
read-only release/nightly workflow and manually before the candidate release
is frozen.

## Compatibility

Historical Table I records remain immutable provenance. The public replay must
agree semantically and within the declared tolerances, not reproduce the
composite TSV byte-for-byte. No old staged object, paper-driver report, or
internal packet object becomes a public serialized format.

The accepted historical He screening result remains data-only. The public
assembly example validates the stable algebra but does not claim its basis,
field, or energy reproduces that calculation. Retired moment-polish packet
provenance remains rejected. Existing internal packet readers retain only
their current internal compatibility contract.

HFDMRG and its solver-dependent paper rows remain separately versioned and do
not block this GaussletBases surface.

## Implementation Surface And Budgets

Approved source ownership:

```text
src/GaussletBases.jl
src/pqs_matched_h2plus.jl
src/cartesian_reference_density/CartesianReferenceDensity.jl
src/cartesian_reference_density/screened_hartree_correction.jl
```

Approved test ownership:

```text
test/ida/runtests.jl
test/runtests.jl
test/pqs_h2plus_table1_release_runtests.jl
```

The focused release/nightly workflow and two example files are allowed only
for the contracts above. No general example runner or release framework is
approved.

Stop-and-report source budgets, including root wiring, are:

- matched H2+ wrapper: preferred at most `300`, hard at most `360` added
  source lines;
- public screened-Hartree layer: preferred at most `90`, hard at most `130`
  added source lines;
- examples together: hard at most `180` lines;
- tests together: hard at most `220` lines.

These are review and stop conditions, not instructions to compress equations,
validation, or readable implementation. If a coherent implementation exceeds
a hard limit, make no source commit and return the smallest justified revision
to this design manager.

The private paper driver may remain because it owns broader historical report,
supplemented-H2, and fixed-state operations. Do not copy its input parser,
report schema, timing fields, or unrelated interaction path into the public
wrapper. Reuse accepted numerical owners directly and account for any
remaining duplicate H2+ mechanics in the handoff.

## Explicit Non-Goals And Failure Rule

This authority adds no public staged parent/shell/transform object, packet
builder or reader, mixed-basis Hartree producer, density fitter, potential
fitter/evaluator, solver, artifact, checkpoint, Ximg/XHF surface, HFDMRG
surface, interaction/exchange correction, general four-index engine, or
release tag.

If either example requires a new scientific approximation, public stage
schema, general field-construction engine, or paper-local data dependency,
stop without a source commit and report the exact missing operation to this
same design manager.
