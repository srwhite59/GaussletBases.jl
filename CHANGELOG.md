# Changelog

## Unreleased

### Fixed
- Radial `IntegralDiagonal` multipoles (`multipole_matrix`, `atomic_operators`) now avoid
  artificial global-scaling underflow. The previous `r^L` / `r^-(L+1)` scaling returned partially
  (L = 16-17) and then completely (L >= 18) zero matrices on the examined erf grids that start
  at r ~ 1e-23. The kernel now uses a running
  normalization `(r_q / r_p)^L` recursion. Independent extended-precision regressions cover
  the same discrete quadrature through L = 132 at a 1e-12 relative limit. This does not certify
  continuum quadrature accuracy or remove all Float64 range and cancellation limitations.

### Changed (numerics of default angular builds)
- The shell-local angular IDA interaction (`build_atomic_injected_angular_hfdmrg_payload`,
  `..._hfdmrg_hf_adapter`, `..._hf_style_benchmark`, `..._small_ed_benchmark`,
  `build_atomic_fixed_radial_angular_sequence`) now sums multipoles through the shell interaction
  moment `lcap` (new keyword `interaction_lmax = :auto`), as documented in
  `docs/angular_interaction_moment_span_note.md`. Previously the sum silently stopped at the
  multipoles stored by `atomic_operators`, `L <= 2 * lmax`, e.g. L <= 12 for the common
  `lmax = 6`, while the moment tables reach L = 34 / 46 / 52 / 66 for NΩ = 50 / 98 / 130 / 200.
  Reviewed packet cases: the Be2+ pair energy at s = 0.25 and NΩ = 98 drops by 0.41 µHa, and
  by 2.0 µHa at NΩ = 130; tested s-shell HF shifts are below 1e-12 Ha, not a general energy bound.
  `interaction_lmax = :stored` reproduces the old cap, not the old arithmetic bitwise;
  an integer sets an explicit cap bounded by the moment span.
- `atomic_operators` takes `multipole_lmax` (default `2 * lmax`, unchanged). It also keeps its
  quadrature samples, so the angular builders evaluate missing multipoles on demand.
- Payload, HF-style, and sequence-export metadata record the cap used: `interaction_lmax`,
  `interaction_lmax_required`, `multipole_lmax_stored`, `interaction_truncated`, and the
  `manifest/interaction/*` keys.
- Radial multipole checksum changes can change radial, sequence and level identities even
  in stored mode. Cap requests also enter sequence/level identities. SHA-pinned receipts
  and old/new level mixtures require deliberate downstream-owner review.

## v0.2.1

### Fixed

- Stabilized centered and displaced Gaussian Coulomb arithmetic, including
  determinant cancellation, damping, and product-local polynomial moments.
- Restored three v0.2.0 compatibility exports and the original timing carrier.
- Isolated annotated-tag verification from checkout's local tag references.

### Changed

- Allow five radial construction attempts, retaining early exit and explicit
  no-refinement behavior. Construction cost roughly doubles for affected setups.
- Report unmet construction and quadrature convergence criteria on exhaustion.
  Quadrature schedules and tolerances remain unchanged; returned grids are best
  effort, not a guarantee of whole-matrix convergence or energy accuracy.
- Expanded public documentation and validation. The near-origin inverse-radius
  limitation remains; source relocations do not introduce new numerical methods.

## v0.2.0

### Changed

- Reduced PQS construction cost through support-local shell-seed construction,
  safe call-local workspace reuse, and terminal-buffer reuse.
- Accelerated terminal Gaussian sums with an exact-order four-element path.
- Made path-aware CI fail closed: candidate and code changes retain the full
  numerical gates; proven docs-only main pushes use lightweight package/docs checks.
- Example 41 output and all unchanged release assertions now share one matched H2+
  comparison instead of executing the comparison twice.

### Scope

- These changes preserve accepted public numerics. v0.2.0 is the supported public
  package version closest to software used in the separate PQS and reference-density
  Hartree-screening work, not an exact archive of either paper's complete computational history.

## v0.2.0-rc2

### Added

- A public residual-GTO/MWG working-system constructor through
  `cartesian_residual_gto_mwg_system`.
- Version-1 external Cartesian-GTO transfer with a strict reader, a
  checkpoint-only PySCF exporter, and caller-thresholded closest-determinant
  preparation.

### Changed

- Split public CI into separate Supported-floor, PQS, and Screening gates.

### Fixed

- Removed invalid package exports that did not name usable bindings.

## v0.2.0-rc1

### Added

- A bounded matched H2+ comparison through `PQSH2PlusRow`,
  `PQSH2PlusComparison`, and `pqs_h2plus_comparison`.
- Supplied-field screened-Hartree assembly with explicit exact/fitted field
  identity, a typed correction result, and public correction accessors.
- Public examples 39-41 for PQS/White-Lindsey H2+, fixed-density screening,
  and the focused H2+ comparison table, with bounded release validation.
- Versioned documentation deployment for exact semantic-version tags.

### Changed

- Declared compatibility ranges for all six direct dependencies while retaining
  Julia 1.10 as the supported minimum.
- Expanded reader documentation for the current Cartesian and nested methods,
  including a curated API reference.
- Kept pull-request documentation build-only, deployed `main` to `/dev/`, and
  confined deployment credentials to the deployment step.

### Fixed

- Repaired legal unsplit-H2 packet construction, rectangular/endcap provenance,
  and forwarding of existing endcap policy, `q`, and `L` diagnostics.

### Public-surface reduction

- Removed six PRF-specific root exports while retaining the parent residual
  function implementations as private diagnostic and provenance machinery.

### Scope

- GaussletBases provides basis and operator construction plus bounded
  supplied-field screened-Hartree assembly.
- Self-consistent-field and correlated solvers, and paper-scale calculation
  campaigns, remain external to the package.
