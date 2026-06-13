Manager review for pass 170

Accepted with one note for downstream interpretation.

The generator now populates `routes/white_lindsey` from the route-configured
diatomic atom-growth White-Lindsey path, not the old seed/one-center debug
path. The populated WL route is correctly labeled as final-basis ordinary
Cartesian/Qiu-White data, with two-body data represented as a final-basis
density-density interaction matrix using `:ggt_nearest`.

The observed WL final dimension is `2287`. This differs from the older generic
diatomic config smoke's non-atom-growth dimension, but it is consistent with
the atom-growth smoke contract, which treats the retained dimension as a route
output rather than a fixed reference value. The artifact must therefore be
understood as a comparison of two route-produced final Hamiltonian
representations for the same Be2 fixture, not as a same-dimension matrix
comparison and not as PQS pre-final/source-box WL data.

The generated JLD2 bundle is large, about 241 MB, because it now stores both the
PQS inspection matrices and the WL 2287-by-2287 overlap/H1/density-density
matrices. That is acceptable for the ignored private inspection artifact, but
it should not be promoted to a committed artifact or default test fixture.

Before handing back to CR2, rerun the generator after committing/pushing so the
producer metadata records the final clean commit and `dirty = false`.

-- repo-manager@macmini
