## Angular Interaction Moment Span Note

This note records the current interaction-span policy for the experimental
shell-local injected angular line.

### Policy

For exact angular content through `lmax`, a density-density interaction needs
product moments through at least `Lmax = 2*lmax`.

For the mixed shell-local injected basis, `2*l_inject` is therefore only a
lower bound. It is not the final working interaction rule.

The repo now follows the same design principle as the legacy
`sphgatomps.jl` line:

- keep the bare injected-sector shell moments through `L <= l_inject` as the
  exact low-`l` lower-bound surface
- grow a separate per-shell interaction Y-moment table adaptively beyond that
  cutoff
- record per-shell interaction `lcap`, `lexpand`, and tail diagnostics
- assemble the interaction from those expanded shell moment tables rather than
  truncating at the injected cutoff

### Why this matters

Using only `L <= l_inject` is enough for the shell-local exact injected
subspace, but it is too small for the full mixed shell basis.

That truncation was the remaining low-order Ne failure mode:

- one-body assembly was already on the correct branch
- HFnn and HFDMRG agreed on the same repo payload
- the remaining mismatch was that the interaction assembly omitted the first
  nontrivial product angular content at low orders

The current repo policy now matches the legacy principle more closely:

- exact injected span and interaction span are distinct concepts
- one-body kinetic span and interaction span are also distinct concepts

This note is intentionally narrow. It records the current repo-side contract
without claiming a broader frozen angular API.

### Radial multipole cap (2026-10)

The assembled interaction is

    V = sum_L 4pi/(2L+1) R_L(a,b) mt_L(a)' mt_L(b).

Its multipole cap is now set explicitly by `interaction_lmax`. Before this
change it was implicitly `length(radial_ops.multipole_data) - 1`, which is
`2 * lmax` of `atomic_operators`. `lmax` there is the one-electron (centrifugal)
range, so the common `lmax = 6` cut the sum at `L = 12` although the expanded
shell tables above extend to `lcap = 34 / 46 / 52 / 66` for
`NΩ = 50 / 98 / 130 / 200`. That contradicted the policy in this note. The
default `:auto` now uses every `L` the tables carry. Missing radial multipoles
are evaluated on demand from the quadrature samples kept by `atomic_operators`.
`:stored` reproduces the old cap.

Raising the cap also needed the stable radial kernel described in
`radial_multipole_stabilization_milestone.md`, because the previous kernel
returned zero matrices for `L >= 18` on production grids.

Pass 653's bounded arithmetic/full-moment repair contract is in
[numerical contracts](src/developer/numerical_contracts.md#Radial-Multipole-and-Angular-Interaction-Repair).
The policy above describes the intended span; the pre-repair assembler still
stops at stored radial multipoles. Historical evidence does not certify it.
