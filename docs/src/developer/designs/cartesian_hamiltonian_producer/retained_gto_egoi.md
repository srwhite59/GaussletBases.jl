# Retained-GTO Local-Product EGOI

EGOI means Exact Gaussian Orbital Integrals. This page is the durable record
for the protected-localized retained-original-GTO experiment; it is not active
implementation authority.

## Status

- `HP-RG-PROTECT-EGOI-AUDIT-01` is completed historical measurement.
- `HP-RG-PROTECT-EGOI-FN-01` is deferred with no source grant or owned path.
- `HP-RG-PROTECT-EGOI-TEST-01` is deferred with no test grant or planned path.

The generic matrix-level EGOI routines remain implemented in
`src/cartesian/hamiltonian_corrections.jl`:

- `egoi_target_product_matrix`;
- `egoi_target_coulomb_matrix`;
- `egoi_density_density_correction`;
- `egoi_stationary_hamiltonian_correction`.

The specialized adapter was never committed to `main`. Its exact `501`-line
source patch is preserved on remote branch
`archive/retained-gto-egoi-wip-2026-08-03` at commit `e8dee9ea3`. That commit
is archival evidence, not a merge candidate or execution grant. The planned
focused test was never created.

## Measured Convention

The completed experiment represented retained original supplement GTOs in the
native protected-localized basis and targeted one `s1+s2` pair per owner. It
used the six owner-local symmetric products, retained `AA-AA`, `BB-BB`, and
`AA-BB` Coulomb residuals in the acceptance metric, and excluded cross-owner
orbital-overlap products.

The `M2` correction variables were same-owner native-row pairs satisfying

```text
r_ij <= 1.75 * max(ell_i, ell_j)
ell_i = max(nearest-neighbor scale_i, core_spacing)
```

Disallowed entries were intended to remain exactly zero. The correction was
in-memory only; no artifact, solver, screening, HF, or public workflow was
part of the experiment.

H, Be, and Be2 probes gave clean target projection, strong training-residual
reduction, small interaction changes, and benign low-Fock shifts. The ignored
Cr2 replay used four targets, six products, and `120435` allowed variables;
it reduced the measured target residual by about `99.6%`, with a relative
interaction Frobenius change of order `1e-4` and no new bad low-Fock mode.
These are positive measurement results, not endpoint or transfer validation.

## Why The Adapter Is Deferred

The archived implementation had no committed caller or test and added `501`
source lines without replacing code. Review also found:

- incomplete duplicate-source, target-Gram, represented-rank, and material
  projection-loss gates;
- corrected-matrix validity reported rather than enforced;
- `max_disallowed_delta_v` reported as a constant instead of independently
  checked;
- an underdetermined fit with `21` independent symmetric target conditions and
  `120435` allowed variables, making regularization and clipping policy a
  material part of the result;
- no holdout or physical endpoint demonstrating transfer beyond the fitted
  products.

This does not reject EGOI as a scientific idea. It rejects carrying the old
specialized adapter as an approved pending implementation.

## Reactivation Boundary

Any specialized retained-GTO EGOI work requires a new docs-only authority
amendment tied to a current physics target. A fresh implementation must reuse
the generic matrix algebra, bind target and row geometry to one construction,
hard-fail provenance/rank/locality/finiteness violations, expose its actual
regularization, include meaningful holdouts, and satisfy current line-budget
and test-cost rules.

No corrected artifact, solver integration, screened-Hartree change, broader
target class, or Cr2 production claim follows from this historical record.
