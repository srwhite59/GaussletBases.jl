# Cartesian R0 Baseline - June 2026

Status: measurement baseline, not implementation authority.

This document records the R0 quantitative baseline after the internal base PQS
Hamiltonian path reached `CartesianIDAHamiltonian` materialization. It is a
measurement snapshot for future carrying-cost and performance comparisons.

## Repository And Machine

- Commit: `2979514492cd4311426bd8f6e19c8be61c3e5a66`
- Branch at measurement start: `main...origin/main`
- Worktree at measurement start: clean
- Julia: `julia version 1.12.6`
- Host: `rh310l.ps.uci.edu`
- Machine/OS:
  `Darwin rh310l.ps.uci.edu 25.5.0 Darwin Kernel Version 25.5.0: Mon Apr 27 20:41:15 PDT 2026; root:xnu-12377.121.6~2/RELEASE_ARM64_T6041 arm64`

## Line And File Counts

Counts were computed with `find . -path "./PATTERN" -type f | sort` and
summed with `wc -l`.

| Pattern | Files | Lines |
| --- | ---: | ---: |
| `src/cartesian*` | 89 | 47713 |
| `src/pqs*` | 11 | 6592 |
| `src/*source_box*` | 6 | 6403 |
| `tools/cartesian*` | 2 | 471 |
| `tools/*pqs*` | 1 | 81 |
| `test/driver_inputs/*pqs*` | 8 | 80 |
| `docs/src/developer/designs/cartesian_hamiltonian_producer/*` | 23 | 3488 |
| `docs/src/developer/roadmaps/*` | 2 | 285 |

## H2 Base Hamiltonian Smoke

Command:

```sh
julia --project=. tools/h2_pqs_base_hamiltonian_smoke.jl
```

Cold fresh-process result:

- elapsed: `30.078453875 s`
- final dimension: `471`
- materialization type: `CartesianIDAHamiltonian{Float64}`
- `K` symmetry error: `3.9968028886505635e-15`
- `V` symmetry error: `1.6653345369377348e-15`
- `U_1` symmetry error: `1.1102230246251565e-16`
- `U_2` symmetry error: `1.0061396160665481e-16`
- H1 lowest delta from reviewed baseline
  `-0.79460371733658908`: `0.0`
- self-Coulomb delta from reviewed baseline
  `0.4569117646737212`: `1.3322676295501878e-15`
- artifact readback one-body delta: `0.0`

Cheap stage timings from the same cold run:

- `cartesian_pair_terms`: `0.662161 s`, `1.083 MiB`
- `cartesian_assembly`: `0.763517 s`, `1.074 MiB`

## H2 Warm Same-Process Measurement

Command:

```sh
julia --project=. tmp/work/r0_h2_base_hamiltonian_warm_measure.jl
```

Generated ignored script:

- `tmp/work/r0_h2_base_hamiltonian_warm_measure.jl`

The script included `tools/h2_pqs_base_hamiltonian_smoke.jl` once as warmup,
then measured a second include with Julia-level `@elapsed` and `@allocated`.

Warm measured result:

- measured include elapsed: `0.665858167 s`
- measured include allocation: `550.6315002441406 MiB`
- smoke internal harness elapsed: `0.369739875 s`
- final dimension: `471`
- H1 lowest delta: `0.0`
- self-Coulomb delta: `1.3322676295501878e-15`
- readback one-body delta: `0.0`

Cheap stage timings from the warm measured include:

- `cartesian_transforms`: `0.022864 s`, `117.606 MiB`
- `cartesian_pair_terms`: `0.000030 s`, `14.359 KiB`
- `cartesian_assembly`: `0.000047 s`, `14.484 KiB`

## Light Separated-Diatomic Status

Command:

```sh
julia --project=. tmp/work/light_separated_diatomic_one_body_validation.jl
```

This existing ignored validation script still runs. It tried several small
diatomic fixtures and selected `n2_separated`.

Selected system:

- atom symbols: `("N", "N")`
- nuclear charges: `(7, 7)`
- bond length: `8.0 bohr`
- `q = 5`
- `core_side = 5`
- `core_spacing = 0.15`
- `xmax_parallel = 10.0`
- `xmax_transverse = 4.0`
- parent axis counts: `(x = 11, y = 11, z = 23)`
- terminal roles:
  `(:atom_local_core, :atom_local_core, :atom_local_shell, :atom_local_shell, :atom_local_shell, :atom_local_shell, :midpoint_slab, :shared_molecular_shell, :z_low_outer_mismatch_slab, :z_high_outer_mismatch_slab)`
- support counts:
  `(125, 125, 218, 218, 386, 386, 81, 1002, 121, 121)`
- coverage counts: `(duplicates = 0, missing = 0, outside = 0)`
- final retained dimension: `1063`
- terminal basis max cross overlap: `1.0957649991595716e-14`
- largest local workspace: `49.26056671142578 MiB`

Stage timings:

- `cartesian_system`: `0.000000 s`
- `cartesian_recipe`: `0.000000 s`
- `cartesian_parent`: `0.169267 s`
- `cartesian_shells`: `1.659348 s`
- `cartesian_units`: `19.849568 s`
- `cartesian_transforms`: `9.838328 s`

One-body result:

- `K` shape: `(1063, 1063)`
- `K` finite: `true`
- `K` symmetry error: `2.842170943040401e-14`
- `U1` shape: `(1063, 1063)`
- `U1` finite: `true`
- `U1` symmetry error: `2.636779683484747e-16`
- `U2` shape: `(1063, 1063)`
- `U2` finite: `true`
- `U2` symmetry error: `3.3306690738754696e-16`
- H1 lowest diagnostic: `-24.938572197228144`
- one-body elapsed: `1.384113 s`
- one-body allocation: `1187.859 MiB`

## Cr2 Status

Cr2 was not rerun in this baseline pass. Current intended status:

- terminal basis previously realized during Slice A work;
- full Cr2 Hamiltonian stress/performance is not closed;
- Cr2-scale Hamiltonian validation remains a later roadmap gate and was not
  run here because the pass requested no full Cr2 Hamiltonian unless explicitly
  approved.

## Commands Run

```sh
git rev-parse HEAD
git status --short --branch
hostname
julia --version
git diff --check
zsh -lc '... find/wc line-count loop ...'
find docs/src/developer/designs/cartesian_hamiltonian_producer -type f | sort
find docs/src/developer/roadmaps -type f | sort
ls tools | rg "h2_pqs_base|h2_pqs|pqs|cartesian"
find tmp/work -maxdepth 1 -type f | sort | rg "h2|n2|terminal|wire|one_body|baseline|r0"
julia --project=. -e 'using GaussletBases; println("load ok")'
julia --project=. tools/h2_pqs_base_hamiltonian_smoke.jl
sed -n '1,220p' tools/h2_pqs_base_hamiltonian_smoke.jl
julia --project=. tmp/work/r0_h2_base_hamiltonian_warm_measure.jl
julia --project=. tmp/work/light_separated_diatomic_one_body_validation.jl
uname -a
sed -n '1,220p' docs/src/developer/roadmaps/README.md
sed -n '1,220p' docs/src/developer/roadmaps/cartesian_long_range_roadmap.md
```

## Not Measured

- Full Cr2 Hamiltonian stress/performance was not run; it is intentionally
  deferred.
- OS-level maximum RSS was not measured; this pass used Julia-level timing and
  allocation measurements per repo policy.
- No new committed test or implementation helper was added.

## Validation

- `git diff --check`: passed for the final docs diff.
- Package load: passed with `load ok`.
- H2 base Hamiltonian smoke: passed.
- Existing ignored light separated-diatomic validation: passed.
