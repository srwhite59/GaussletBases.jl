Accepted.

The audit found the right convention boundary. The existing Coulomb expansion
stores positive coefficients/exponents for `1/r`; the active WL/QW by-center
storage uses `-coefficients` to build a negative unit-charge attraction matrix,
records the nuclear charge, and applies charge only later during Hamiltonian or
diagnostic assembly.

The old CCPM source-box lane is a useful oracle, but it is the wrong production
shape for the new PQS retained source-mode route because it builds a positive
Gaussian sum and then applies `-Z` in a physical wrapper. The new CPBM path
should not import that wrapper as authority.

The audit also identified the useful implementation split:

```text
first: source-space 3D nuclear contraction from caller-supplied term-first
       Gaussian factor matrices
later: source-axis helper that produces those matrices from actual PQS source
       axis facts
```

That split avoids guessing the source-axis object while still making progress
on sign, by-center metadata, and retained contraction.

Manager validation:

- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`

Next target:

Implement the minimal CPBM PQS source-space electron-nuclear by-center block
from supplied term-first 1D Gaussian factor matrices. Keep it uncharged,
by-center, and source-box-first. Do not generate Gaussian factors from parent
axis objects yet.

-- repo-manager@macmini
