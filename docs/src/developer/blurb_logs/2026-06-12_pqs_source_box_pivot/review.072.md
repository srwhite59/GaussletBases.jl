Review 072: accepted.

The convention is now clear enough to implement the support electron-nuclear
helper:

- one separated matrix record per center;
- support ordering is core rows followed by shell rows;
- uncharged by-center matrices carry the negative unit-charge attraction
  `-1/r_A`;
- Hamiltonian assembly applies charge later as `Z_A * V_A`;
- no center summation in the helper;
- origin PGDG factors are only valid for the matching origin center;
- off-origin support should use centered Gaussian factor construction from
  axis layers and Coulomb exponents;
- fixed-block/WL/raw-support probes remain oracle only.

Proceed with the helper, but do not widen into H1/RHF/driver work.

-- repo-manager@macmini
