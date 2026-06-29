# Residual Gaussian Orthogonality Robustness

Status: approved narrow source authority under `HP-RG-ORTHO-FN-01`,
`HP-RG-ORTHO-TEST-01`, `HP-RG-IDTOL-FN-01`, and
`HP-RG-IDTOL-TEST-01`, plus the superseding cutoff/tolerance policy
`HP-RG-CUTOFF-FN-01` and `HP-RG-CUTOFF-TEST-01`.

## Reason

Strict N2 q5 p10 at `core_spacing = 0.042857` failed only the final residual
identity validation:

```text
max |G' S R|       = 1.776e-14
max |R' S R - I|   = 1.673e-10
current tolerance  = 1.0e-10
owner retained     = 9, 9
merge eigen range  = 7.232e-2 .. 1.928
merge condition    = 26.65
```

The owner-local residual metrics are full rank and positive, the retained
counts are unchanged, and the final merge spectrum is healthy. Comparison N2
cases at `core_spacing = 0.075` and `0.05` pass with smaller final identity
errors. This is a final floating-point validation overshoot, not evidence for
rank loss, wrong owner grouping, severe conditioning, or a residual-selection
algorithm change.

Be atom cc-pV5Z `lmax = 1` later reproduced the same issue class with a
higher-zeta supplement:

```text
fixture                  = Be atom, ns = 5, core_spacing = 0.075, radius = 20
supplement               = cc-pV5Z, lmax = 1, contracted
retained residual count  = 21
minimum occupation       = 6.151e-6
merge condition          = 1.0
max |G' S R|             = 1.776e-14
max |R' S R - I|         = 2.183e-10
current allowed error    = about 2.000e-10
```

The base path is healthy, owner-local occupations are all above the
legacy/current `1.0e-8` residual occupation cutoff, and the final merge is
well-conditioned. Raising the residual occupation cutoff to `1.0e-5` would drop
a real retained direction and is a stronger basis-selection change than this
evidence supports.

## Approved IDs

- `HP-RG-ORTHO-FN-01` - robust final residual orthogonalization and validation.
- `HP-RG-ORTHO-TEST-01` - validation gates for the robustness pass.
- `HP-RG-IDTOL-FN-01` - final residual identity tolerance default.
- `HP-RG-IDTOL-TEST-01` - validation gates for the tolerance-default pass.
- `HP-RG-CUTOFF-FN-01` - prior residual cutoff/tolerance default.
- `HP-RG-CUTOFF-TEST-01` - validation gates for the prior cutoff policy.
- `HP-RG-CUTOFF-FN-02` - production residual occupation cutoff tightening.
- `HP-RG-CUTOFF-TEST-02` - residual-only validation for the tightened cutoff.

## Approved Source Surface

Approved files:

```text
src/cartesian_residual_gaussians/residual_basis.jl
src/cartesian_final_basis_realization/pqs_terminal_residual_gto.jl
```

`residual_basis.jl` is the primary owner. The terminal residual file is allowed
only for narrow existing internal keyword plumbing if the current compatibility
entry point needs to pass an approved tolerance/check option through to the RG
domain function.

No public export or public API change is approved.

## Approved Behavior

`HP-RG-ORTHO-FN-01` may make final residual normalization and identity
validation robust for small floating-point overshoots when the owner-local
selection and final merge are otherwise healthy.

Allowed approaches:

- use a numerically stable symmetric final merge normalization;
- recompute or validate the final residual overlap through an explicitly
  symmetrized matrix before measuring `R' S R - I`;
- evaluate final identity with a combined absolute/relative check:

```text
err_RR <= tau_identity_abs + tau_identity_rel * max(1, scale_RR)
```

where `err_RR = max(abs.(symmetrized(R' S R) - I))`,
`tau_identity_abs = 1.0e-10`, `tau_identity_rel = 1.0e-10`, and `scale_RR` is
the maximum absolute entry or an equivalent infinity-norm scale of the
symmetrized final residual overlap;
- keep the existing `G' S R` validation strict, with the same absolute
orthogonality tolerance used by the approved RG path unless a later amendment
changes it;
- report both absolute errors (`G' S R`, `R' S R - I`), retained owner counts,
  and final merge spectrum in the implementation handoff.

This is not a broad tolerance relaxation. The relative component is approved
only for final residual identity validation after owner-local selection has
succeeded and the final merge metric is positive and not near-singular under
the existing merge failure rule.

`HP-RG-IDTOL-FN-01` additionally approves changing the default final residual
identity validation tolerance to:

```text
identity_atol = 1.0e-8
```

This supersedes only the default final `R' S R` identity acceptance threshold.
It does not change the residual basis selection algorithm. This older default
policy is superseded for production by `HP-RG-CUTOFF-FN-01`; legacy
`PureGaussianGausslet.jl:addinGaussians(...; RGcut = 1.0e-8)` remains
historical context only, and width/zeta filtering remains explicit and
user-controlled.

The identity tolerance is a final validation/cleanup tolerance. It may accept a
small `R' S R - I` overshoot only when owner-local metric checks, final merge
metric checks, and `G' S R` orthogonality checks remain healthy. It is not a
direction-retention cutoff, an automatic zeta/width filter, eigenvalue
flooring, or a reason to alter the final merge failure rule.

`HP-RG-CUTOFF-FN-01` superseded the earlier production defaults:

```text
residual_occupation_cutoff = 5.0e-8
identity_atol = 5.0e-8
```

Evidence: Cr atom PQS `basis_ns = 9`, `map_ns = 11`, `lmax = 1` retained a
marginal residual direction at occupation `3.637e-8`. The production policy is
to discard such marginal residual directions by default, rather than retain
them because the older cutoff happened to be lower. `identity_atol = 5.0e-8`
keeps the final validation tolerance aligned with that stricter retained-space
policy.

This update changes only the default owner-local residual occupation cutoff and
the default final `R' S R` identity validation tolerance. Owner-local grouping,
negative-eigenvalue tolerances, final merge metric checks, `G' S R`
validation, width/zeta filtering, MWG/IDA, artifacts, driver workflow, public
API, and source ownership remain unchanged.

`HP-RG-CUTOFF-FN-02` supersedes only the residual occupation cutoff for
production:

```text
residual_occupation_cutoff = 1.0e-6
identity_atol = 5.0e-8
```

Evidence: Cr2 residual spectra showed the worst low-H1 modes were built from
marginal owner-local residual directions with occupations around `1.27e-7` to
`8.98e-7`; a `1.0e-6` cutoff drops `6` directions per owner. Broad residual
widths remain diagnostic evidence, but width filtering is not the first
production rule because one-center atoms can have broad candidates without the
same bad `H1_RR` sector.

This update keeps final identity validation at `5.0e-8` and changes only the
default owner-local residual occupation cutoff. It does not approve kinetic or
`H1_RR` spectral guards, width-filtering defaults, full HF, Cr2 artifact
workflow, or residual-selection algorithm changes.

## Forbidden

This amendment does not approve:

- residual selection semantic changes;
- global raw-candidate residual selection;
- global raw-column pivoted-Cholesky residual selection;
- occupation-cutoff changes beyond the explicit `HP-RG-CUTOFF-FN-01` default
  and `HP-RG-CUTOFF-FN-02` default updates;
- kinetic or `H1_RR` spectral guards without a separate amendment;
- numerical negative-eigenvalue tolerance changes;
- final merge eigenvalue flooring to retain low-occupation directions;
- width filtering as a conditioning repair;
- width/zeta filtering default changes;
- owner grouping changes;
- final merge metric failure-rule changes;
- MWG/IDA, nuclear, raw-block, exact-operator, artifact, driver, public API, or
  route-stage changes;
- status/report fields or payload objects;
- committed fixtures/tests unless a later amendment names the exact file.

If passing strict N2 requires changing residual selection, supplement
construction, or final-basis construction, implementation must stop and report
the blocker instead of widening this lane.

## Validation

`HP-RG-ORTHO-TEST-01` approves only:

- existing H2 residual-GTO/MWG endpoint;
- H2 base/supplemented readback only if the compatibility/facade path is
  touched;
- ignored strict N2 q5 p10 residual audit or artifact smoke at
  `core_spacing = 0.042857`;
- one passing N2 comparison at `core_spacing = 0.05` or `0.075`;
- reporting `max |G' S R|`, `max |R' S R - I|`, owner retained counts, and the
  final merge eigenvalue range/condition.
- Be atom cc-pV5Z `lmax = 1` residual audit/artifact path passes with the same
  `21` retained residual directions;
- Be atom cc-pVDZ `lmax = 1` still passes;
- H2 residual-GTO/MWG endpoint remains unchanged;
- reporting for the Be tolerance pass includes `max |R' S R - I|`, allowed
  tolerance, retained count, minimum retained occupation, final merge
  condition, and `max |G' S R|`.
- Cr atom `basis_ns = 9`, `map_ns = 11`, `lmax = 1` residual construction
  passes or cleanly drops the marginal `s4` direction at occupation `3.637e-8`
  as intended;
- Be atom cc-pV5Z still passes;
- H2 residual-GTO/MWG endpoint remains unchanged.
- The existing H2 endpoint test
  `test/nested/cartesian_r3a_h2_augmented_one_body_runtests.jl` may update only
  its two cutoff assertions from `1.0e-8` to `5.0e-8`:
  `residual.occupation_cutoff` and `values[:occupation_cutoff]`.
- Cr2 residual-only audit after `HP-RG-CUTOFF-FN-02`, not full HF or a Cr2
  artifact/workflow run;
- Cr2 owner retained counts should drop from `68 + 68` to `62 + 62`;
- report `min eig(K_RR)`, `min eig(H1_RR)`, and low-mode candidate
  composition;
- if low-H1 ghost modes remain, stop and request separate kinetic/`H1_RR`
  spectral-guard authority;
- the same H2 endpoint test may update only its two cutoff assertions from
  `5.0e-8` to `1.0e-6`.

No Cr2 full Hamiltonian, Cr2 artifact, Cr2 facade support, new committed test
file, other committed test/fixture edit, driver workflow, artifact schema
change, solver/RHF, ECP, or EGOI work is approved by this robustness lane.

If Be cc-pV5Z cannot pass by changing only the final identity tolerance
default, implementation must stop and report the exact blocker. Do not change
residual selection, drop directions, alter width/zeta filtering defaults,
change merge rules, alter MWG/IDA conventions, or expand the source surface.

If the Cr atom case cannot pass or cleanly drop the marginal direction by
changing only the two approved defaults, make no source commit in the later
implementation pass and report the exact blocker.
