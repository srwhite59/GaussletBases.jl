## Experimental ordinary homonuclear chain line

This repo now includes an experimental ordinary Qiu-White chain constructor for
homonuclear linear chains:

- `bond_aligned_homonuclear_chain_qw_basis(...)`

Milestone A landed:

- ordinary QW product basis only
- combined inverse-sqrt map on the chain axis
- one shared transverse projection map
- cheap geometry-first validation on small H-chain cases

Milestone B now also lands:

- `ordinary_cartesian_qiu_white_operators(::BondAlignedHomonuclearChainQWBasis3D; ...)`
- cheap H-chain ordinary-operator smoke coverage on the same experimental line

Milestone C now begins:

- experimental nested-chain geometry and split-tree diagnostics
- a first geometry-backed nested fixed-block object for small chain tests
- exploratory repeated-midpoint chain split policy only
- no production nested-chain operator claim yet

Current experimental nested-chain status:

- even-chain central binary splitting looks coherent on cheap H-chain cases
- odd-chain nested policy is still under active experimental development
- the repo now exposes explicit odd-chain policy diagnostics rather than hiding
  the current no-split outcome behind one implicit heuristic
- the repo is not yet claiming a settled odd-chain split/local-resolution rule

Milestone E now also lands:

- a first experimental nested-chain ordinary-QW operator path
- odd chains on that operator path use `odd_chain_policy = :central_ternary_relaxed`
  explicitly
- `:strict_current` remains the conservative default/reference policy for the
  lower-level geometry path
- the nested-chain operator line is still experimental and not yet a final
  chain contract

What is explicitly deferred here:

- any production nested-chain operator route
- any frozen split/local-resolution or endcap policy
- molecular-shell supplement routes
- any HF/DMRG study contract for chains

Current chain-policy notes should be treated as guidance for this experimental
line, not as a frozen implementation contract.
