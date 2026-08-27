# Changelog

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
