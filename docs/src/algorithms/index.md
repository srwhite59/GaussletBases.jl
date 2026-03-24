# Algorithms

This section records basis-construction and operator-construction algorithms in
an implementation-facing form.

The intended reading order inside each page is:

1. pseudocode
2. short code pointers with file names
3. short references
4. supporting detail in descending order of human importance

For new basis-construction routes, the algorithm page should be written before
or alongside implementation so that the construction order is explicit before
the code grows around it.

## Code-comment convention

When an algorithm step is implemented in code, the corresponding code block
should carry a short searchable comment of the form:

```julia
# Alg QW-RG step 4: Define residual Gaussians by orthogonalizing 3D GTOs
# to the full 3D gausslet space.
# See docs/src/algorithms/qiu_white_residual_gaussian_route.md.
```

The comment should:

- use the same step number as the algorithm page
- stay close to the pseudocode wording
- include the docs path exactly
- be short enough to remain readable inside the source

## Current pages

- [Qiu-White residual-Gaussian route](qiu_white_residual_gaussian_route.md)
- [1D distorted-gausslet PGDG refinement hierarchy](distorted_gausslet_pgdg_refinement_hierarchy.md)
- [Cartesian nested face construction](cartesian_nested_face_construction.md)
- [Cartesian nested atomic nonrecursive route](cartesian_nested_atomic_nonrecursive_route.md)
- [Cartesian nested diatomic box policy](cartesian_nested_diatomic_box_policy.md)
- [Cartesian nested diatomic coordinate distortion](cartesian_nested_diatomic_coordinate_distortion.md)
- [Radial interval-sampled build and extents](radial_interval_sampled_build_and_extents.md)
