"""
Internal owner for exact Cartesian Gaussian raw blocks.

This module owns neutral raw numerical kernels that consume parent proxy axes,
Cartesian Gaussian supplements, Coulomb expansions, and nuclear centers, then
return raw matrices for downstream route/domain consumers. It does not own
terminal projection, residual Gaussian transforms, Hamiltonian assembly,
artifact writing, route reports, or public API.

File map:
- `nuclear_blocks.jl`: uncharged by-center Gaussian nuclear `G-A`/`A-A`
  raw blocks.
"""
module CartesianGaussianRawBlocks

include("nuclear_blocks.jl")

end # module CartesianGaussianRawBlocks
