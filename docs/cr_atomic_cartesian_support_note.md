# Cr Atomic Cartesian Support Note

## Scope

This is a narrow repo-support note for the chromium / `Cr2` line.

It records the current ordinary/hybrid Cartesian capability relevant to an
atomic Cr comparison, the actual blocker that was present in the live repo, and
the minimum viable extension that makes the atomic test scientifically honest
without overbuilding for the later two-atom case.

## Present capability

There are three distinct pieces to keep separate.

1. Legacy named-basis loading:

- `src/legacy_basis_adapter.jl`
- public entry point: `legacy_atomic_gaussian_supplement(atom, basis_name; lmax, basisfile, ...)`

The loader itself is broader than the vendored data file. It can read arbitrary
atoms and basis names from a legacy `BasisSets`-style file supplied through:

- `basisfile = ...`
- or `GAUSSLETBASES_BASISSETS_PATH`

The repo-vendored basis file is intentionally tiny:

- `H cc-pVTZ`
- `H cc-pVQZ`
- `He cc-pVTZ`
- `He cc-pVQZ`

So current repo-local *vendored coverage* is minimal, but the loader contract
is already externally extensible and is not hardwired to H/He.

2. Legacy one-dimensional hybrid builder:

- `src/ordinary_hybrid.jl`
- legacy/internal constructor: `hybrid_mapped_ordinary_basis(...)`

This branch remains intentionally centered, separable, and `s`-only for active
legacy atomic supplements. It rejects any `LegacyAtomicGaussianSupplement` with
non-`s` shell content. That restriction is real, but it belongs only to this
older 1D hybrid builder, which is now quarantined from the supported public
surface.

3. Atomic ordinary-QW Cartesian route:

- `src/ordinary_qw_operator_assembly.jl`
- `src/ordinary_qw_raw_blocks.jl`
- public entry point: `ordinary_cartesian_qiu_white_operators(...)`

This is the scientifically relevant atomic Cartesian comparison path. It now
uses the explicit atomic-centered 3D Cartesian shell representation for all
active atomic supplement content up to `lmax = 2`, including pure `s` shells.

After the present pass, the atomic QW path now supports:

- `lmax = 0`: `s`
- `lmax = 1`: `s, p_x, p_y, p_z`
- `lmax = 2`: `s, p_x, p_y, p_z, dxx, dyy, dzz, dxy, dxz, dyz`

The separate two-center molecular shell route remains intentionally narrower
for now and still stops at `lmax <= 1`.

## Actual blocker

For an atomic Cr Cartesian comparison, the live blocker was **not** the legacy
loader itself.

The real blocker was that the active atomic ordinary-QW shell route stopped at
`lmax <= 1`, even though the current Cr atom scientific contract is already
centered on the checked `lmax = 2` atomic referee.

The concrete stop points were:

- `src/legacy_basis_adapter.jl`
  - `_atomic_cartesian_shell_labels(...)`
  - `_atomic_cartesian_shell_supplement_3d(...)`
- `src/ordinary_qw_raw_blocks.jl`
  - `_qwrg_atomic_derivative_terms(...)`

So the atomic QW route had true `SP` support, but not yet true `D` support.

That meant:

- the older 1D hybrid branch looked `s`-only because it is
- the active atomic QW branch was actually `sp`-only before this pass
- neither state was enough for a scientifically meaningful atomic Cr test on
  the current `lmax = 2` contract

## Minimum viable extension for Cr atom testing

The minimum viable repo extension is:

- keep the vendored basis library narrow
- rely on `basisfile = ...` or `GAUSSLETBASES_BASISSETS_PATH` for the actual Cr
  basis block
- extend only the **atomic** explicit 3D Cartesian shell route from `lmax <= 1`
  to `lmax <= 2`
- do **not** broaden the separate one-dimensional hybrid builder
- do **not** broaden the two-center molecular shell route yet

That is the correct first move because the first chromium question is atomic:

- can the repo build an honest atomic Cartesian QW comparison on the same
  low-`l` footing as the checked Cr atom `Ylm` referee?

For that question, `d` support is the decisive missing ingredient.

## What was changed in this pass

Atomic-only `d` support was added on the explicit QW shell route:

- `src/legacy_basis_adapter.jl`
  - atomic Cartesian shell labels now include the six Cartesian `d` functions
  - atomic explicit shell supplement now allows `lmax <= 2`
- `src/ordinary_qw_raw_blocks.jl`
  - atomic shell derivative helper now supports powers `0, 1, 2`

Scope boundary kept explicit:

- atomic ordinary-QW / nested-QW: now through `lmax <= 2`
- one-dimensional hybrid builder: still `s`-only
- bond-aligned molecular shell supplement route: still `lmax <= 1`

## What should be done first for Cr atom

Do **not** start with the older `hybrid_mapped_ordinary_basis(...)` branch.

The correct first repo-side Cr Cartesian test is:

- atomic ordinary-QW Cartesian route
- external Cr basis block supplied through `basisfile = ...`
- `legacy_atomic_gaussian_supplement("Cr", basis_name; lmax = 2, basisfile = ...)`
- `ordinary_cartesian_qiu_white_operators(...)`

That is the first scientifically meaningful atomic Cartesian comparison because
it matches the current Cr atom low-`l` contract much more closely than the
older centered separable hybrid branch.

## Atom versus later two-atom scope

This extension is enough for the **atom**.

It is only partially on the path to the later two-atom case:

- it proves the repo can now carry true atomic `d` shell content on the active
  Cartesian QW route
- but it does **not** yet solve the later two-center/molecular placement
  problem for `d` shells
- the current bond-aligned molecular supplement route is still intentionally
  narrower and remains a separate future extension

So this pass should be read as:

- atomic Cr comparison support: yes
- later Cr2 molecular shell support: not yet

## Validation

Focused atomic smoke run after the extension:

- built a temporary atomic `S/P/D` legacy basis block
- verified atomic `lmax = 2` supplement expansion into explicit Cartesian `d`
  labels
- verified `ordinary_cartesian_qiu_white_operators(...)` builds and produces a
  finite atomic QW object
- verified the separate two-center molecular shell route still rejects
  `lmax = 2`

The milestone-specific smoke printed:

- `cr-atomic-cartesian-spd-ok`
