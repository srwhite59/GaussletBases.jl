> **Status:** supporting development note. For the current ordinary branch,
> read `docs/current_ordinary_branch.md` first. For this note chain, use
> `docs/ordinary_pgdg_supporting_notes.md`.

# Ordinary Mapped PGDG: Backend Pivot

This note records the first controlled architectural pivot for the mapped
ordinary one-body branch.

## 1. What remains the reference route

The mapped numerical ordinary path remains the validation route.

That is still the trusted comparison path for:

- overlap
- kinetic
- Gaussianized one-body factors from the Coulomb expansion
- end-to-end mapped Cartesian hydrogen checks

It is not being removed.

## 2. What changes in this pass

For one-body mapped ordinary work in the mild-to-moderate distortion regime,
the experimental PGDG-style analytic path becomes the preferred implementation
route.

The distinction is explicit:

- `:numerical_reference`
- `:pgdg_experimental`

So the package can now lean on the analytic route while still keeping the
numerical mapped path alive as the reference and validation layer.

## 3. Why this is justified

The distortion-regime study showed that the relevant White–Lindsey standard is
already met well enough in the physically relevant regime:

- nearly identical span/subspace
- almost identical projection behavior
- acceptable mapped hydrogen energies

Strong-distortion cases remain useful stress tests, but they are not the main
target regime for this experimental ordinary one-body backend.

## 4. What this pass is not doing

This is still narrow.

It is not:

- a broad public PGDG surface
- a two-electron mapped ordinary layer
- a He workflow
- a claim that the strong-distortion stress-test regime is already finished

It is only the point where the one-body mapped ordinary branch can start to
prefer the analytic PGDG-style implementation while continuing to validate
against the numerical reference route.
