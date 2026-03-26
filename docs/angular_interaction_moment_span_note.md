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
