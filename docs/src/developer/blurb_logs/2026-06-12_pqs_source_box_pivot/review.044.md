Review 044: accepted; the independent complete-shell final basis is real.

Pass 044 implemented the missing final-basis realization in the right place,
`CartesianFinalBasisRealization`, without exporting a new broad public surface.
The helper:

```text
pqs_complete_core_shell_final_basis
```

materializes the intended route-owned one-center basis:

```text
direct core support count:       125
surrounding shell support count: 218
shell retained count:             98
total final dimension:           223
core/shell support disjoint:     true
pre-final overlap rank:          223
final overlap identity error:    ~7.1e-14
```

This resolves the pass 043 rank deficiency. The shell sector is now the
surrounding shell outside the core, not boundary columns inside the same
`5 x 5 x 5` source space.

The implementation is appropriately scoped: final-basis realization only, no
H1 operator placement, no IDA/RHF/GTO/driver work, and no permanent test burden.
The old White-Lindsey fixed-block packet remains an oracle for geometry/counts
only.

Next pass can now add the one-body operator transfer for this final basis:
build overlap, kinetic, and uncharged by-center nuclear blocks over the
route-owned combined core/shell support layout, transform through the new final
coefficients, assemble the Z=1 Hamiltonian, and run the ordinary symmetric H1
probe. Keep this as a focused probe first; do not add a permanent test unless
it replaces older oracle pressure.

-- repo-manager@macmini
