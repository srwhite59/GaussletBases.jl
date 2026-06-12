Review result:

Accepted as a useful specialization-shape audit. It did not change code, but
it identified a concrete first cleanup target: the hot decomposed WL inventory
result still encodes route size in tuple-valued summary fields.

Key evidence:

```text
field                         side-7 synthetic        side-15 Be shape
parent axis counts             (7,7,7)                 (15,15,15)
retained dimension             223                     615
unit count                     27                      131
pair count                     378                     8646
coefficient matrix size        (343,223)               (3375,615)
unit pair storage              UnitPairIndexTable      UnitPairIndexTable
retained units storage         Vector                  Vector
unit_keys type                 NTuple{27,Symbol}       NTuple{131,Symbol}
unit_summaries type            27-element tuple        131-element tuple
pair_summaries type            NTuple{378,...}         Symbol
```

Interpretation:

The production data structures are already vector/table based:
`Vector{RetainedUnitRecord}` and `UnitPairIndexTable`. The remaining
shape-specialization pressure comes from report-style fields on the hot
inventory result:

- `unit_keys = Tuple(unit.unit_key for unit in units)`
- `unit_summaries = Tuple(...)`
- small-inventory `pair_summaries = Tuple(...)`

That explains why the side-7 synthetic precompile workload does not fully cover
side-15 Be. It compiles a 27-unit inventory result and a 378-pair summary tuple,
while Be needs a 131-unit inventory result and omits pair summaries.

Stable concepts review:

- `ParentAxisContext3D` is still a good route-level struct candidate, but it is
  not the first cleanup to implement.
- `DecomposedWLInventorySummary` or a smaller count/status summary is the first
  practical target because it directly removes type-level unit-count summaries
  from the hot result.
- `DecomposedWLFactorizedBasis3D`, `ResidualMomentMatrices`, and
  `AtomGTOFinalBasisMatrices` remain good future concepts, but implementing
  them now would broaden the pass.

Deletion/shrinkage review:

No production code, tests, metadata, or compatibility path became unnecessary
in this audit-only pass. No tests were added. The next cleanup should shrink
hot inventory summary shape, not add coverage.

Validation reviewed:

- `julia --project=. tmp/work/atom_gto_specialization_shape_audit.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- artifact: `tmp/work/atom_gto_specialization_shape_audit_summary.txt`

Next target:

Change the hot decomposed WL inventory result so `unit_keys`,
`unit_summaries`, and small `pair_summaries` no longer encode inventory size in
their concrete tuple types. Preserve compact access for probes/tests, but move
detailed tuple-shaped summaries out of the production route shape or make them
vector/count based.
