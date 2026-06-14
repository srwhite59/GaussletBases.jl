Review 220 - accepted

Pass 220 ran a local ignored probe for the H2 R=4 q5 no-supplement
WL/old-QW 463 reference values. It made no tracked source, test, or bin edits.

Accepted result:

```text
final dimension = 463
fixed block size = (1215, 463)
parent axis counts = (9, 9, 15)
support counts = (275, 578, 362)
retained counts = (251, 98, 114)
retained order = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
residual count = 0
supplement policy = :none
final overlap identity error = 1.3816121305989328e-14
```

Reviewed no-supplement WL values:

```text
wl_h1_lowest = -0.7946609179724673
wl_h1_orbital_self_coulomb = 0.45696639804337047
wl_rhf_electronic_energy = -1.1589518556651142
wl_rhf_one_electron_energy = -1.5611571934181985
wl_rhf_electron_electron_energy = 0.40220533775308426
wl_rhf_nuclear_repulsion = 0.25
wl_rhf_total_with_nuclear_repulsion = -0.9089518556651142
```

Probe provenance:

```text
constructor: bond_aligned_homonuclear_qw_basis(R=4.0, core_spacing=0.5,
  xmax_parallel=6.0, xmax_transverse=4.0, bond_axis=:z)
fixed block source: _nested_bond_aligned_diatomic_source + _nested_fixed_block
operator source: ordinary_cartesian_qiu_white_operators(fixed_block;
  supplement omitted)
interaction treatment: :ggt_nearest
gausslet backend: :numerical_reference
```

Timing:

```text
complete-rectangular source and fixed block: 10.55 s
gausslet-only ordinary QW operators:         4.79 s
restricted closed-shell HF:                  0.46 s
```

Current blocker:

```text
:wl_h2_gausslet_only_reference_values_not_yet_wired_to_driver_artifact
```

The next implementation pass should wire these values as compact comparison
fields. It must keep supplemented WL/QW chemistry references quarantined and
must label RHF electronic energy separately from the nuclear-repulsion-included
total.

Deletion/shrinkage:

```text
deleted: none
simplified: none
quarantined: old supplemented WL/QW H2 scalar references remain non-comparable
not deleted because: local probe/audit only
exact remaining caller/blocker:
  :wl_h2_gausslet_only_reference_values_not_yet_wired_to_driver_artifact
```

-- repo-manager@macmini
