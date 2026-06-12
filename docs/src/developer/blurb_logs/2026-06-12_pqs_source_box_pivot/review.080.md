Review 080: accepted.

The pass makes the explicit-box multi-layer PQS source entry point read as
compatibility/probe bridge machinery rather than route shellification
authority. The code now adds compact summary labels:

- `explicit_box_bridge`;
- `shellification_backed_geometry`.

The explicit-box path sets those to `true` and `false`, respectively. The
region-plan path remains the shellification/lowering-backed route authority.
The explicit-box docstring and the near-term PQS plan were updated to match.

`source_kind` was left stable, which is the right choice here. The authority
downgrade is now explicit in metadata and docs without historical-symbol churn.

Manager validation:

- `git diff --check` passed.
- `julia --project=. -e 'using GaussletBases; println("load ok")'` passed.

Deletion/shrinkage:

- conceptual authority of explicit boxes was downgraded;
- no tests were added or changed;
- the explicit-box entry point remains as compatibility/probe bridge.

-- repo-manager@macmini
