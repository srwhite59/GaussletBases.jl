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

What is explicitly not claimed here:

- no nested planar splitting policy
- no square-lattice operator milestone yet
- no molecular-shell supplement route
- no HF/DMRG workflow or application contract
- no settled square-lattice policy inherited from the diatomic or chain notes

The current diatomic and chain policy notes should be treated as guidance for
related experimental lines, not as a frozen square-lattice implementation
contract.
