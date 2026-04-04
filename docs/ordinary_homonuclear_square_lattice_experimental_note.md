## Experimental ordinary homonuclear square-lattice line

This repo now includes an experimental ordinary Qiu-White square-lattice basis
constructor for axis-aligned homonuclear `n x n` lattices in the `xy` plane:

- `axis_aligned_homonuclear_square_lattice_qw_basis(...)`

Milestone A lands only the ordinary basis, mapping, and geometry layer:

- one combined multi-center inverse-sqrt map on `x` from the unique lattice
  `x` coordinates
- one combined multi-center inverse-sqrt map on `y` from the unique lattice
  `y` coordinates
- one shared single-center inverse-sqrt map on `z` at the common lattice plane
- cheap geometry-first validation on small H square lattices

Milestone B now also lands:

- `ordinary_cartesian_qiu_white_operators(::AxisAlignedHomonuclearSquareLatticeQWBasis3D; ...)`
- cheap H `2 x 2` and `3 x 3` ordinary-operator smoke coverage on the same
  experimental line
- no nested planar splitting policy or planar nested operator claim yet

Milestone C now begins:

- an experimental planar nested-geometry and split-tree diagnostics line for
  tiny square lattices
- a first geometry-backed experimental planar fixed-block probe
- explicit planar candidate reporting rather than an implicit split policy
- no production planar nested-operator claim yet

What is explicitly not claimed here:

- no nested planar splitting policy
- no nested planar operator milestone
- no molecular-shell supplement route
- no HF/DMRG workflow or application contract
- no settled square-lattice policy inherited from the diatomic or chain notes

The current diatomic and chain policy notes should be treated as guidance for
related experimental lines, not as a frozen square-lattice implementation
contract.
