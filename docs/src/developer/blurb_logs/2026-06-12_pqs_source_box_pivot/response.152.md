Pass 152 response - audit Be2 PQS Ham handoff contract

No source edits, no Julia commands, no commits.

Recommended next object:

```text
_PQSDiatomicCompleteCoreShellHamiltonianHandoffPayload
```

I recommend a private factor-first Hamiltonian handoff object, not a public
export/HamV6 adapter yet.

Minimal fields:

```text
status
blocker
route_family
source_plan
final_basis
h1_payload
ham_input_payload
overlap
one_body_hamiltonian
density_interaction
pre_final_pair_matrix
final_to_pre_final_coefficients
pre_final_weights
support_weights
support_pair_raw_numerator
raw_pair_factor_terms
center_records
nuclear_repulsion
electron_count
spin_sector
ordering
conventions
available_objects
missing_objects
summary
metadata
```

Shape decision:

- Prefer the private factor handoff object first.
- Include a dense inspect view only if it is explicitly labeled
  `private_inspect_only` and tied to the `pre_final_localized_positive_weight`
  density gauge.
- Do not create a compatibility adapter toward HamV6 or the cartesian bundle
  writer in the first pass.
- Do not claim CR2/HFDMRG readiness until electron count, nuclear repulsion,
  and solver-basis/order conventions are explicit.

What the current Be2/PQS payloads already satisfy:

- Source plan:
  - object kind: `:pqs_diatomic_complete_core_shell_source_plan`
  - support order: `(:product, :pqs_left, :pqs_right)`
  - route retained order: `(:pqs_left, :pqs_right, :product)`
  - source-plan support row order:
    `:core_product_then_shell_left_right_pqs`
  - core support count `25`, shell support count `250`
  - precleanup/final dimension `221`
- Final basis:
  - status: `:available_pqs_complete_core_shell_final_basis`
  - support row order: `:core_then_shell`
  - final dimension `221`
  - pre-final coefficients, final coefficients, and combined Lowdin cleanup are
    available through the final-basis object.
- H1:
  - status: `:available_diatomic_complete_core_shell_h1_payload`
  - final kinetic, by-center electron-nuclear matrices, final one-electron
    Hamiltonian, and H1 solve are materialized.
  - center count is `2`; H1 summary records support/pre-final ordering.
- Ham input/electron-electron:
  - status: `:available_diatomic_complete_core_shell_ham_input_payload`
  - density provenance available from source-plan bundles
  - support weights materialized
  - raw pair factor terms available
  - support raw-pair numerator materialized in diatomic support-row order
  - pre-final density interaction materialized
  - density gauge: `:pre_final_localized_positive_weight`
  - raw-pair convention: `:raw_numerator`
  - support weight count `275`
  - support raw-pair shape `(275, 275)`
  - pre-final pair matrix shape `(221, 221)`
  - final-to-pre-final coefficients and pre-final weights are available inside
    the density interaction object.
- Nuclear metadata:
  - current readiness center summary has center count `2` and nuclear charges
    `(4, 4)`.
  - H1 helper builds by-center electron-nuclear records, so center records are
    available internally, but no handoff-owned center/nuclear-repulsion package
    exists yet.

What remains missing before CR2 should consume it:

- A named private handoff object that declares the solver-basis order and which
  matrices/vectors are in that order.
- Nuclear repulsion / center-center term for Be2 total-energy comparisons.
- Electron count, charge, and spin/closed-shell convention. For neutral Be2 the
  likely electron count is `8`, but the handoff must state it rather than make a
  CR2-side guess.
- A decision whether CR2 consumes:
  - pre-final density interaction plus `final_to_pre_final_coefficients`, or
  - a dense density-density matrix in the final solver basis.
- If a dense `Vee` is provided, it needs an explicit derivation/convention name.
  Current downstream scripts often read dense `ham/interaction_matrix`, while
  HFDMRG accepts dense `H,V` density-density input. That is useful compatibility
  pressure, but it should not silently collapse the PQS pre-final gauge into a
  public export.
- Basis centers/orbital labels or a deliberate statement that they are
  unavailable/non-authoritative for this private route.
- Explicit nonclaims: no H1-J scalar diagnostic, no RHF/SCF, no HamV6, no
  public export/artifact, no WL comparison payload.

Recommended next implementation pass:

Add a private helper:

```text
_pqs_source_box_route_driver_diatomic_complete_core_shell_hamiltonian_handoff_payload
```

Wire it after `diatomic_complete_core_shell_ham_input_payload` in
`cartesian_assembly(...)`.

First pass should:

- require available source plan, final basis, H1 payload, and Ham-input payload;
- carry references to those payloads rather than copying their internals;
- expose compact inspect fields:
  - final dimension
  - overlap/final-basis orthonormality convention
  - one-body Hamiltonian reference
  - density interaction reference
  - pre-final pair matrix shape
  - final-to-pre-final coefficient shape
  - density gauge
  - raw-pair convention
  - support/pre-final/final ordering labels
  - center count, nuclear charges, nuclear coordinates if available
  - nuclear repulsion status
  - electron-count/spin status
- leave public/export/artifact/HamV6/CR2/HFDMRG flags false;
- keep readiness blocked on export/consumer convention until the consumer shape
  is explicitly chosen.

Proposed status/blocker labels:

Available private handoff:

```text
:available_diatomic_complete_core_shell_hamiltonian_handoff_payload
```

Blocked generic:

```text
:blocked_diatomic_complete_core_shell_hamiltonian_handoff_payload
```

Input blockers:

```text
:missing_diatomic_complete_core_shell_source_plan
:missing_diatomic_complete_core_shell_final_basis
:missing_diatomic_complete_core_shell_h1_payload
:missing_diatomic_complete_core_shell_ham_input_payload
```

Convention/data blockers:

```text
:missing_diatomic_hamiltonian_handoff_ordering_convention
:missing_diatomic_nuclear_repulsion
:missing_diatomic_electron_count_convention
:missing_diatomic_spin_sector_convention
:missing_diatomic_dense_interaction_convention
```

Next readiness blocker after private handoff exists:

```text
:missing_diatomic_hamiltonian_consumer_contract
```

or, if manager wants to stay CR2-specific:

```text
:missing_diatomic_cr2_inspection_contract
```

Evaluation of options:

- Private dense handoff with `H1` + dense/factorized electron-electron:
  useful for CR2/HFDMRG shape pressure, but too easy to overclaim if first pass
  calls it an export or solver-ready dense `Vee`.
- Private factor handoff:
  best first step. It preserves the PQS pre-final density interaction and
  ordering metadata while still giving CR2 enough structure to inspect.
- HamV6/export compatibility adapter:
  premature. Current export docs are producer-side public contracts for atomic,
  sliced, angular, and experimental chain surfaces, not Be2/PQS.
- Sharper blocked payload only:
  acceptable as fallback, but less useful now because all private ingredients
  except the handoff conventions are already materialized.

Validation:

```text
git status --short --branch
```

Result:

```text
## main...origin/main
```

Deletion/shrinkage forecast:

- deleted: none in the next pass unless manager explicitly retires old report
  aliases
- simplified: readiness can move from `:missing_diatomic_ham_consumer_contract`
  to a compact handoff payload status plus one sharper consumer/export blocker
- quarantined: dense `Vee`, HamV6/export, CR2/HFDMRG execution, RHF/SCF, H1-J,
  and WL comparison should remain nonclaimed
- not deleted because: current compatibility/report fields and old cartesian
  bundle/export surfaces are still separate live contracts
- exact remaining caller/blocker:
  `cartesian_assembly(...)` currently exposes
  `diatomic_complete_core_shell_ham_input_payload`; the next blocker is the
  missing private Hamiltonian handoff/consumer convention,
  `:missing_diatomic_hamiltonian_consumer_contract`

-- repo-doer@macmini
