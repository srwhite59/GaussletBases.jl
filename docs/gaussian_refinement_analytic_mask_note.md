# Gaussian Refinement Analytic Mask Note

## Purpose

This note records the explicit analytic Gaussian quasi-refinement formula that
is now relevant to the proposed **1D distorted-gausslet PGDG refinement
hierarchy**.

The immediate repo question is narrow:

can the first tabulated hierarchy mask be taken directly from the analytic
Gaussian refinement formula, rather than only from a fitted numerical table?

## Core formula

Let

```math
\phi_\sigma(x)=\frac{1}{\sqrt{2\pi}\sigma}\exp\!\left(-\frac{x^2}{2\sigma^2}\right)
```

be a normalized one-dimensional Gaussian. If `Σ > σ`, define

```math
\tau^2=\Sigma^2-\sigma^2.
```

Then there is an exact continuous convolution identity

```math
\phi_\Sigma(x-s)=\int_{\mathbb R}\phi_\tau(t-s)\,\phi_\sigma(x-t)\,dt.
```

Sampling that identity on a uniform lattice of spacing `h` gives the explicit
Gaussian quasi-refinement

```math
\phi_\Sigma(x-s)\approx h\sum_{k\in\mathbb Z}\phi_\tau(kh-s)\,\phi_\sigma(x-kh).
```

So the refinement-mask coefficients are Gaussian:

```math
c_k \propto \phi_\tau(kh-s)
      \propto \exp\!\left(-\frac{(kh-s)^2}{2\tau^2}\right).
```

This is the main practical formula.

## Ternary specialization

For the ternary hierarchy, take coarse spacing `H = 3 h` and a fixed
shape ratio

```math
\Sigma = \rho H,\qquad \sigma = \rho h.
```

Then the analytic refinement mask for the `1 -> 1/3` step is a centered
Gaussian mask with width set by

```math
\tau^2=\Sigma^2-\sigma^2
=
\rho^2(H^2-h^2).
```

With `H = 1` and `h = 1/3`, this becomes

```math
\tau^2=\rho^2\left(1-\frac{1}{9}\right)
=
\frac{8}{9}\rho^2.
```

So a direct analytic candidate for the first hierarchy mask is

```math
c_k \propto \exp\!\left(-\frac{(k/3)^2}{2(8\rho^2/9)}\right)
=
\exp\!\left(-\frac{k^2}{16\rho^2}\right),
```

up to normalization and finite-window truncation.

## Important interpretation

This formula is a **quasi-refinement** formula for the Gaussian proxy layer.
It is not an exact discrete refinement identity in the usual wavelet sense.

That distinction matters:

- for the current practical mask question, quasi-refinement is enough
- for later exact distorted-gausslet convergence questions, some further
  correction or reprojection may still be desirable

## Direct fit versus recursive refinement

Two different uses should be kept separate.

### Direct coarse-to-current-level fit

Fit the original coarse Gaussian directly by the current finer Gaussian line.

- This improves with finer current level.
- At fixed `rho`, it has a nonzero saturation floor.

### Recursive level-to-level refinement

Apply the same `1 -> 1/3` transfer repeatedly:

```text
1 -> 1/3 -> 1/9 -> 1/27 -> ...
```

This gives a practical repeated local hierarchy, but by itself it is still an
approximate refinement ladder rather than an exact discrete nesting theorem.

## Current practical question

The current repo question is deliberately narrower than the full hierarchy
theory:

- build the analytic `1 -> 1/3` Gaussian refinement mask from this formula
- compare it directly to the fitted tabulated mask already studied numerically
- decide whether the analytic mask can serve as the default source of the
  first stored hierarchy mask

At this stage, the key comparison is:

- coefficient-by-coefficient agreement
- constant-reproduction behavior
- single-level Gaussian refinement error
- repeated-refinement behavior

not yet full transferability to all distorted-gausslet integral data.

## Practical reading

The local numerical studies already show that the fitted masks are:

- symmetric
- positive
- local
- non-cancelling
- and robust under Float64 tabulation

So if the explicit analytic mask closely matches those fitted tables, it is a
strong candidate for the first stored hierarchy mask.

## Related literature cue

This note was prompted by external discussion pointing to the Gaussian
approximate-refinement / approximate-MRA literature, especially work by
Beylkin, Monzon, and Satkauskas.

For current repo work, the most important usable result is the explicit
Gaussian quasi-refinement mask above.
