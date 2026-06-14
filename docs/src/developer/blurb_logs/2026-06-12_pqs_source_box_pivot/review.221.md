Review 221 - accepted

Pass 221 wired the reviewed H2 R=4 q5 no-supplement WL reference values into
the physical PQS driver artifact and made the H2 physical gausslet endpoint a
real WL/PQS comparison endpoint.

Accepted endpoint state:

```text
route_kind = :bond_aligned_diatomic_physical_gausslet_core_shell_pqs
support counts = (275, 578, 362)
retained counts = (251, 98, 114)
final dimension = 463
supplement policy = :none
comparison ready = true
physics endpoint ready = true
```

Reviewed WL values now carried by the H2 physical input/artifact:

```text
wl_h1_lowest = -0.7946609179724673
wl_h1_self_coulomb = 0.45696639804337047
wl_rhf_electronic_energy = -1.1589518556651142
wl_rhf_nuclear_repulsion = 0.25
wl_rhf_total_with_nuclear_repulsion = -0.9089518556651142
```

Observed PQS-vs-WL deltas:

```text
delta_h1 = 2.55351295663786e-15
delta_h1_j = 6.661338147750939e-16
delta_rhf_electronic_energy = -3.2713831643604863e-12
delta_rhf_total_with_nuclear_repulsion = -3.2713831643604863e-12
```

The implementation correctly avoided an ambiguous molecular `wl_rhf_total`.
It writes the electronic RHF comparison separately from the nuclear-repulsion
included total. Supplemented WL/QW scalar references remain quarantined through
`comparison/old_supplemented_wl_qw_scalar_references_blocked = true`.

Validation rerun by manager:

```text
julia --project=. -e 'using GaussletBases; println("load ok")'
  passed

julia --project=. -e 'using Test; t = @elapsed include("test/nested/cartesian_ham_builder_h2_pqs_q5_physical_gausslet_r4_target_runtests.jl"); println("elapsed_s=", t)'
  passed: 91 assertions, elapsed_s=87.252018209

julia --project=. test/nested/cartesian_ham_builder_h2_pqs_q5_source_box_diagnostic_r4_runtests.jl
  passed: 12 assertions, test summary 1m38.9s

git diff --check
git diff --cached --check
  passed
```

Line budget under `src`, `test`, `bin`, and the CR2 generator was still
net-negative before commit:

```text
102 added
105 deleted
net -3
```

Committed and pushed:

```text
c486288d Wire H2 WL comparison values
```

Next direction:

The H2 no-supplement gausslet endpoint is now the accepted comparison gate.
The old 221-dimensional H2 source-box diagnostic test/input is no longer the
physics target and should not keep imposing a slow explicit test unless a live
caller needs it. The next pass should retire that scaffold if caller audit
allows, then audit the MWG/GTO supplement seam under the architectural rule:
WL and PQS share the physical support/shell plan and differ in retained
transform/operator construction, not in supplement theory.

-- repo-manager@macmini
