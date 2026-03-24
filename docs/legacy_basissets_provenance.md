# Legacy BasisSets Provenance

GaussletBases now vendors a curated narrow legacy-style `BasisSets` file at:

- `data/legacy/BasisSets`

## What this file is

This vendored file is not a general basis-library subsystem. It is a small
repo-local compatibility file for the current atomic named-basis supplement
line.

The file is intentionally narrow. It currently contains only:

- `H cc-pVTZ`
- `H cc-pVQZ`
- `He cc-pVTZ`
- `He cc-pVQZ`

Those are the atom/basis blocks currently exercised by the active atomic
supplement path, the first bond-aligned diatomic supplement path, and the
current test suite.

## Provenance

The numerical basis data came from Basis Set Exchange. They were then
reformatted into the legacy `ReadBasis.jl` / `BasisSets` block style used by
the historical GaussletModules tooling and by the loader in this repo.

Minor formatting normalization was applied for compatibility with that loader,
for example keeping the plain `e`-style numeric format expected by the legacy
parser.

The vendored file should therefore be understood as:

- Basis Set Exchange data
- curated to the narrow subset used here
- normalized into the legacy block format required by this code path

## Loader policy

The loader search order is:

1. explicit `basisfile = ...`
2. `GAUSSLETBASES_BASISSETS_PATH`
3. vendored repo copy `data/legacy/BasisSets`
4. legacy fallback `~/BasisSets`

So users who need unsupported atoms or basis names still have a clear extension
path without changing repo code.

## Extension path

Users needing a broader basis selection should supply either:

- `basisfile = ...`
- or `GAUSSLETBASES_BASISSETS_PATH`

This pass does not add a general import/conversion tool. If that becomes
necessary later, it should be a separate small utility rather than a broad
built-in basis-library subsystem.
