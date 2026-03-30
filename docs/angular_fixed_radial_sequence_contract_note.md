# Fixed-Radial Angular Sequence Contract

This note records the current narrow producer contract for the fixed-radial
`N_sph` angular ladder.

What stays fixed across levels:

- one prebuilt radial substrate
- radial shell ids
- physical radial shell centers
- radial shell count
- one fixed angular-profile setting line:
  `beta`, `l_inject` request, `tau`, `whiten`, and gauge version

What changes across levels:

- `N_sph`
- the cached shell-local angular profile selected for that `N_sph`
- the resulting shell dimension
- the dense per-level bridge Hamiltonian built on that profile

Under the current contract, the shell-local construction is radius-independent
for fixed profile settings. That means every source-target shell-local profile
overlap inside one fixed-radial sequence is common across shells and can be
exported once as a shell-independent sidecar. The current sequence object keeps:

- `adjacent_overlaps` for neighboring `N_sph[k] -> N_sph[k+1]` pairs
- `direct_overlaps` for the full non-adjacent upper triangle inside the same
  sequence

What is not included yet:

- no common-target embedding export
- no smaller-to-max lift matrices
- no DMRG-side Givens or restart logic
- no many-body consumer contract beyond the existing per-level dense bridge
