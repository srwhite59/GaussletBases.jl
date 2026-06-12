Review 042: accepted; the smaller blocker is correctly framed.

The pass did not force the private White-Lindsey low-order seed/fixed-block
packet into the active PQS route, which is the important outcome. That packet
has the desired `223` shape:

```text
core retained range:   1:125
shell retained range:  126:223
total dimension:       223
```

but its operator source is `:nested_fixed_block` and its inventory status is
`:private_development_seed`. It may be an oracle/reference, not the active H1
authority.

The precise next blocker is:

```text
:missing_route_owned_combined_core_shell_retained_operator_blocks
```

That is narrower than the pass 041 blocker. The current route-owned final-basis
seam can handle the `98`-function PQS boundary shell, but it does not have a
combined retained layout or operator-block placement for:

```text
core-core
core-shell
shell-shell
```

for overlap, kinetic, and uncharged by-center nuclear blocks.

Next pass should implement the smallest combined core/shell retained-basis
surface needed to assemble those three block families and run the one-center
H1 probe. Keep the private fixed-block packet as an oracle only, and keep IDA,
density-density, RHF, GTO, and driver work out of scope.

-- repo-manager@macmini
