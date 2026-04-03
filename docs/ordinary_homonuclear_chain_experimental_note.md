## Experimental ordinary homonuclear chain line

This repo now includes an experimental ordinary Qiu-White chain constructor for
homonuclear linear chains:

- `bond_aligned_homonuclear_chain_qw_basis(...)`

This milestone is intentionally narrow:

- ordinary QW product basis only
- combined inverse-sqrt map on the chain axis
- one shared transverse projection map
- cheap geometry-first validation on small H-chain cases

What is explicitly deferred here:

- nested chain splitting
- molecular-shell supplement routes
- any HF/DMRG study contract for chains

Current chain-policy notes should be treated as guidance for this experimental
line, not as a frozen implementation contract.
