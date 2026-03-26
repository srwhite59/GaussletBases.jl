# Angular Research Track

This page records the smallest current repo-side statement of the angular
gausslet line.

It is an active research track, not a frozen public workflow branch.

## Current interpretation

The present angular direction in `GaussletBases.jl` should be read as:

- manuscript-facing and experimental
- centered on shell-local injected angular basis work that has not yet been
  imported into the package
- intentionally narrower than the mature radial/atomic producer line

The repo now contains the first controlled scaffold needed to start that
import, plus the first shell-local experimental construction on top of it:

- the existing atomic `(l,m)` and Gaunt/sectorized IDA foundation
- a read-only vendored full sphere-point-set access layer
- an explicit curated fixture subset for tiny tests and pinned examples
- a shell-local injected angular basis constructor on one curated sphere-point
  set
- a first shell-to-atom angular assembly layer over shell radii
- a first atom-side one-electron angular benchmark path built on that assembly

It does **not** yet contain:

- the full `sphgatomps`-style workflow
- optimizer/sweep drivers for sphere-point generation
- manuscript figure scripts
- a full mixed-basis HamIO/HamV6 export contract
- a broader downstream migration beyond the current HF-facing bridge surfaces

## Near-term target

The near-term scientific target is an **atomic angular benchmark ladder**:

1. HF
2. small ED
3. one HamIO / HFDMRG-facing HF bridge

That is the first clean way to prove the angular line through the same
producer/consumer boundary now established for the current atomic IDA branch.

## Sphere-point backing store

The normal/default angular backing store is now the full vendored optimized
sphere-point collection under:

- `data/angular/SpherePoints.jld2`
- `sphere_point_set_orders()`
- `sphere_point_set(order)`

This full JLD2 collection is the default pool used by order-based shell-local
and shell-assembly helpers.

The smaller curated subset remains available through:

- `curated_sphere_point_set_orders()`
- `curated_sphere_point_set(order)`

Those curated helpers are experimental and read-only fixture surfaces. They are
there for tiny tests, pinned examples, and paper-stable demo cases, not as the
normal angular order pool. The first shell-local basis layer is now available
through:

- `build_shell_local_injected_angular_basis(order; ...)`
- `build_shell_local_injected_angular_basis(point_set; ...)`
- `shell_local_injected_angular_diagnostics(shell)`

This remains an experimental research-track surface, not a frozen public API.

The first shell-to-atom assembly layer is now available through:

- `assign_atomic_angular_shell_orders(shell_radii; ...)`
- `build_atomic_shell_local_angular_assembly(shell_radii; ...)`
- `atomic_shell_local_angular_diagnostics(assembly)`

This assembly layer is still angular-only. It stops before Coulomb assembly,
HF/ED workflow wiring, or end-to-end atomic angular drivers.

The first narrow atom-side benchmark path is now available through:

- `build_atomic_injected_angular_one_body_benchmark(radial_ops; ...)`
- `atomic_injected_angular_one_body_diagnostics(benchmark)`

This benchmark stays on the one-electron central-potential line. It is there to
prove that the shell-local angular assembly supports a real atom-side Galerkin
calculation before the later benchmark ladder stages are imported.

The first narrow angular HF-style benchmark is now available through:

- `build_atomic_injected_angular_hf_style_benchmark(radial_ops; ...)`
- `atomic_injected_angular_hf_style_diagnostics(benchmark)`

This remains benchmark-oriented and density-density / manuscript-facing. It is
there to prove that the angular line can support a real atom-side mean-field
solve before the later small-ED and bridge-facing stages are imported.

The first narrow angular small-ED benchmark is now available through:

- `build_atomic_injected_angular_small_ed_benchmark(radial_ops; ...)`
- `atomic_injected_angular_small_ed_diagnostics(benchmark)`

This remains a tiny `1 up, 1 down` benchmark on the same He-sized line. It is
there to prove that the shell-local angular assembly supports a real
interacting benchmark beyond the HF stage without importing the whole old
workflow.

The first narrow HamIO / HFDMRG-facing HF bridge surface is now available through:

- `angular_benchmark_exact_hamv6_payload(benchmark; ...)`
- `write_angular_benchmark_exact_hamv6_jld2(path, benchmark; ...)`

This bridge is intentionally honest about the current contract boundary. It
exports the benchmark line's exact common low-`l` reference in the proven
HamV6 / `HamIO` language for the current HF consumer path, but it does **not**
yet claim that the full mixed shell-local angular basis is directly
representable in that consumer language, and it is not yet a true many-body
DMRG bridge.

The first narrow in-memory HFDMRG-facing HF adapter is now available through:

- `build_atomic_injected_angular_hfdmrg_hf_seeds(one_body; nup, ndn)`
- `build_atomic_injected_angular_hfdmrg_hf_adapter(radial_ops; ...)`
- `build_atomic_injected_angular_hfdmrg_hf_adapter(benchmark; ...)`
- `run_atomic_injected_angular_hfdmrg_hf(adapter; ...)`

This adapter reuses the current dense `H, V, psiup0, psidn0` handshake and
lets the mixed-basis angular benchmark line call `HFDMRG.solve_hfdmrg(...)`
without a file round-trip. It now supports explicit open-shell control through
`nup`, `ndn`, `psiup0`, and `psidn0`, while the narrow default seed helper
builds first practical orbital guesses from the assembled one-body orbital
frame if explicit seeds are not supplied. The direct `radial_ops` adapter
entrypoint assembles the Hamiltonian and interaction without running the repo's
internal HF-style benchmark first. It does **not** solve the separate
mixed-basis HamIO/HamV6 export problem.

The intended post-whitening/injection working basis remains orthonormal.
Accordingly, any residual nonidentity part of the final overlap matrix is to be
treated as a conditioning/construction diagnostic, not as a physically
meaningful generalized-overlap model. The internal angular HF-style benchmark
therefore uses the same conventional orthonormal-basis density-density SCF
model as the current HFnn / HFDMRG references, rather than a generalized
overlap reinterpretation.

The durable repo-local boundary note for this state is:

- `docs/angular_consumer_contract_boundary.md`

The vendored full collection currently lives in:

- `data/angular/SpherePoints.jld2`
- `data/angular/SpherePoints_manifest.toml`

The explicit curated fixture subset currently lives in:

- `data/angular/curated_sphere_points.toml`

and carries:

- point-set cardinality/order
- Cartesian coordinates on `S^2`
- nearest-neighbor spread ratio `nn_ratio`
- source tag / source project / source artifact note

## What is deferred

The following items remain explicitly deferred:

- Hooke as its own later dedicated line
- heteronuclear/angular coupling beyond the first atomic benchmark ladder
- any claim that the present angular data API is stable enough to freeze

Hooke remains important, but it should come later as its own workflow/paper
line after the first atomic angular benchmark ladder is real.

## Next import boundary

If this scaffold stays sound, the next real scientific import should be:

- the first cleaner angular consumer-contract formalization beyond the current
  exact-common-reference bridge

That is the step that turns the current one-electron / HF-style / small-ED /
bridge-reference branch into a broader downstream-facing angular line inside
`GaussletBases`.
