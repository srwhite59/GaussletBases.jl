Pass 046 response
=================

Status: audit complete; no source fix made.

Finding
-------

The bad H / Z=1 value is already present before final-basis transfer.

The support-space generalized diagnostic on the same combined core/shell
support rows gives:

```text
support_generalized_h1 = -2.0639248059188007
support_overlap_condition = 1.0000000000000830
```

So the final-basis transfer is not the source of the overbinding. The support
overlap is essentially identity/positive definite, and the same bad energy is
already in the support-space operator matrices.

Final-basis decomposition of the lowest vector:

```text
final_h1                         = -2.0638461028784740
kinetic_expectation              =  1.1012644410060735
uncharged_nuclear_expectation    = -3.1651105438845457
charged_nuclear_expectation_z1   = -3.1651105438845457
decomposition_total              = -2.0638461028784723
decomposition_matches_h1_error   =  1.7763568394002505e-15
```

Charge/sign audit:

```text
center_charge_recorded           = 1.0
nuclear_input_is_uncharged       = true
nuclear_charge_application_stage = hamiltonian_assembly
nuclear_sign_convention          = negative_electron_nuclear_attraction
wrong_nuclear_sign_h1            = 0.14749061717503814
```

This rules out the easy local bugs:

- not a Z=2 mislabeled run;
- not charge applied twice;
- not wrong sign in the final Hamiltonian assembly;
- not a generalized-overlap final solve issue;
- not a final-basis transfer issue.

Current blocker:

```julia
:blocked_nonacceptance_h1_operator_convention
```

More precise interpretation:

The support-space product-factor nuclear/operator convention is not yet a
trusted physical H / Z=1 convention for this complete core/shell fixture. The
large nuclear expectation, about `-3.1651`, is already present on support rows.
The next pass should compare the support operator construction against a
trusted same-geometry one-body oracle or audit whether raw Gaussian factor
matrices require an additional projected/operator-normalized convention before
being used on these combined support rows.

Trusted oracle status:

```text
trusted_oracle_status = :blocked_mapping_or_representation_mismatch
```

I did not use the old fixed-block packet as active authority. A same-geometry
trusted oracle was not available without a mapping/representation mismatch.

Probe artifact
--------------

Added ignored audit probe:

```text
tmp/work/pqs_complete_core_shell_h1_convention_probe.jl
tmp/work/pqs_complete_core_shell_h1_convention_probe_summary.txt
```

Source changes
--------------

None in this pass.

Validation
----------

Ran:

```text
julia --project=. tmp/work/pqs_complete_core_shell_h1_convention_probe.jl
julia --project=. -e 'using GaussletBases; println("load ok")'
git diff --check
```

Results:

- convention probe completed and wrote the blocked audit summary;
- load check passed with `load ok`;
- `git diff --check` passed;
- final `git status --short --branch` was clean.

Deletion/shrinkage report
-------------------------

- No old/fallback/oracle surface became less necessary in this pass because no
  source fix was made.
- Nothing was deleted or simplified; this was an audit-only pass that located
  the bad convention before final-basis transfer.
- No permanent test was added. The new coverage is an ignored developer probe,
  not routine test burden.
- Remaining stale/duplicate surfaces to retire next:
  - the old fixed-block packet remains oracle-only;
  - the 98-function boundary-shell H1 path remains mechanical/nonphysical;
  - the complete core/shell H1 route should remain non-acceptance until the
    support-space one-body product/operator convention is validated or fixed.

-- repo-doer@macmini
