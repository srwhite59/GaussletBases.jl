# High-Order q4 Recipe CR2 Validation Handoff

## Status

This is a request note for CR2-facing validation of the internal q4
atom-growth/endcap-panel recipe path. It does not claim CR2 validation,
PySCF agreement, same-density validation, energy validation, or production
support.

Repo-side evidence currently covers:

- opt-in construction of the default q4 recipe into a nested fixed block
- PGDG QW construction smoke through the existing nested fixed-block route
- one dimensionally matched repo-native parent-reference capture/H1 sanity
  check

The CR2 external-target diagnostic used duplicated one-center PySCF Cr ccECP
occupied spaces with source-metric Lowdin in the two-center GTO AO metric. It
verified the exact q4 fixture, grid order, and centers for alpha `(735, 20)`,
beta `(735, 8)`, and combined `(735, 28)` target matrices. The q4 parent
capture reported by CR2 was alpha `0.9558876918` and beta `0.9811637481`.
Conditional on that q4 parent target, the opt-in fixed block captured alpha
`0.9994161870`, beta `0.9998427016`, and combined spin-sum `0.9995403327`.
Worst-column capture was `0.9987695368` for `alpha_right_mo7_col18` and
`0.9998415239` for `beta_left_mo0_col1`; max absolute H1 delta was
`7.249921e-03`. The diagnostic used `:pgdg_localized_experimental`, warning
log count `0`, and no numerical-reference fallback.

This updates the request note with one completed capture/H1 check only. It
does not validate relaxed Cr2 energy, two-electron terms, same-density
energies, production defaults, or public route readiness. The q label remains a
convergence/control parameter, not a validated accuracy tier.

The existing side-29 Cr occupied artifacts are not compatible with this q4
fixture and must not be reused for this validation.

## Exact Fixture Required

CR2 must produce or project the target into the exact parent basis used by the
small q4 smoke fixture:

- basis constructor family: `bond_aligned_homonuclear_qw_basis`
- `family = :G10`
- `bond_length = 5.0`
- `core_spacing = 0.7`
- `xmax_parallel = 8.0`
- `xmax_transverse = 4.0`
- `bond_axis = :z`
- parent axis counts: `(7, 7, 15)`
- parent dimension: `735`
- opt-in fixed dimension for comparison: `469`
- parent coefficient ordering:
  `flat = (ix - 1) * ny * nz + (iy - 1) * nz + iz`

The target rows must correspond exactly to that parent grid and ordering. A
different parent size, different spacing, different molecule, different
axis order, or different coefficient ordering is a blocker.

## Required CR2 Artifacts

The minimum useful handoff is:

- metadata file recording all fixture settings above
- parent grid file with `flat`, `ix`, `iy`, `iz`, and center coordinates for
  all `735` rows
- raw parent-space occupied coefficient matrix, shape `735 x nocc`, for the
  target columns
- occupied-label file mapping each column to spin, source orbital label, and
  occupation or intended electron count
- provenance file or report explaining how the parent-space coefficients were
  generated or projected

If spin-polarized targets are used, provide separate alpha and beta matrices
or a single matrix with unambiguous spin labels. If Lowdin or otherwise
orthonormalized variants are included, label them as secondary diagnostics;
the raw parent-space coefficients should remain the primary input.

## Accepted Provenance

Accepted target provenance includes:

- CR2/PySCF orbitals projected into the exact q4 repo parent basis and grid
  order above
- another external target only if its projection into the same `735`-row
  parent basis is fully documented
- a repo-side target only if it is explicitly labeled as an internal sanity
  target, not CR2/PySCF evidence

The current repo-native sanity target is
`repo_native_parent_one_body_generalized_eigenvectors`. It is useful only as a
dimension/order sanity check. It is not a substitute for the requested CR2
occupied target.

## Forbidden Mismatches

Do not mix this validation with:

- side-29 Cr artifacts with parent dimension `24389`
- side-27 or other Cr fixtures
- count-7 high-order smoke routes that do not share the exact q4 parent
- coefficient matrices whose row count is not `735`
- coefficient matrices whose grid/order contract is not explicitly verified
- QW/operator routes that silently select `:numerical_reference`

For repo-side QW checks, the expected backend is
`:pgdg_localized_experimental`. Numerical-reference routes are allowed only if
a task explicitly asks for a labeled reference/debug comparison; they are not
valid as the default q4 validation path.

## Expected Repo-Side Metrics

Once a matching target exists, repo should report:

- metadata and grid/order verification
- parent dimension and fixed dimension
- backend used, with explicit confirmation that no numerical-reference
  fallback occurred
- per-column and per-spin retained-space capture
- total captured electrons per spin, if occupations are supplied
- worst captured orbital label and value
- parent H1 expectation and fixed-space projected H1 expectation for each
  target column, if the target is compatible with the repo one-body operator
- alpha, beta, and total H1 deltas when spin data are supplied
- receipt/carried-space diagnostics for the PGDG QW route, if operators are
  built

These are capture/H1 validation metrics only. They do not validate
two-electron terms, same-density energies, relaxed HF/UHF, ED, or CR2
production readiness.

## Pass Or Blocker Criteria

The handoff is usable if:

- metadata exactly matches the q4 fixture
- parent grid has `735` rows and agrees with the repo flat-index convention
- coefficient matrices have `735` rows
- parent-space normalization is documented well enough to compute capture
- PGDG is the active repo backend for QW checks
- no numerical-reference fallback is accepted silently

The handoff is blocked if:

- any fixture or ordering field differs from the q4 smoke parent
- only side-29/side-27 Cr data are available
- coefficient provenance cannot explain how external orbitals were projected
  into the repo parent basis
- row count or grid order is ambiguous
- repo-side checks would require changing source builders, defaults, QW
  kernels, or backend policy

## Current Repo-Native Sanity Reference

The repo-native parent-reference artifact used to prepare this request lives
under `tmp/work/`:

- `tmp/work/q4_recipe_parent_reference_capture_report.txt`
- `tmp/work/q4_recipe_parent_reference_capture.tsv`

It records:

- parent dimension: `735`
- fixed dimension: `469`
- target vector count: `8`
- primary capture: `0.9999809353`
- primary H1 delta: `4.69e-5`
- worst capture: `0.9876930450`
- max absolute H1 delta: `0.0140`

Those numbers show the opt-in q4 fixed block can capture a compatible
repo-native parent one-body target. They do not count as CR2 validation.
