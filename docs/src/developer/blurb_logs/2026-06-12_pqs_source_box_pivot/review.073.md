Review 073: needs correction before commit.

The explicit origin-factor path appears to work and the H1 gate now exercises
the new support electron-nuclear helper. However, the implementation also claims
an off-origin `axis_layers` path, and that path looks inconsistent:

```julia
_pqs_multilayer_centered_factor_terms(...)
    factors = ntuple(axis -> gaussian_factor_matrices(...), 3)
    return (factors[1], factors[2], factors[3], :centered_axis_layers)

_pqs_multilayer_validate_factor_terms(axis_terms, term_count)
    ndims(axis_terms[axis]) == 3 || throw(...)
```

`gaussian_factor_matrices(...)` returns a collection of 2D matrices, not a
term-first 3D array. So the off-origin axis-layer path is likely blocked or
wrong even though the response says it is implemented.

Do not commit this pass as accepted until one of these is true:

1. the axis-layer path is converted into term-first 3D arrays and covered by a
   focused test, or
2. the helper/docs/response are narrowed to say only the explicit-factor path is
   implemented and the centered axis-layer path remains blocked.

Preferred: fix and test the axis-layer path, because the convention has already
been documented and this is a local shape issue.

-- repo-manager@macmini
