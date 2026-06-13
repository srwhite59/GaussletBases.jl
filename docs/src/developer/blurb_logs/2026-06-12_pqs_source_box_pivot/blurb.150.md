Pass 150 - audit Be2 PQS post-H1 Hamiltonian seam

Purpose:

Do a no-edit audit of the next Be2/PQS step after private H1 materialization.
The current readiness blocker is:

```text
:missing_diatomic_complete_core_shell_h1_j_consumer
```

but the medium-term user target is a Hamiltonian constructor usable by a
downstream CR2/Be2 comparison of WL and PQS. H1-J is a private diagnostic
validator; it may or may not be the next useful Hamiltonian-constructor seam.

Audit question:

Should the next implementation pass build:

1. a private Be2/PQS H1-J diagnostic payload,
2. a private electron-electron/Ham-input payload,
3. a blocked Ham-readiness payload with sharper missing electron-electron facts,
4. or a small compatibility handoff object for downstream Hamiltonian consumers?

Read first:

- `src/pqs_source_box_route_driver_helpers.jl`
- `src/pqs_multilayer_complete_core_shell_h1.jl`
- `src/pqs_multilayer_support_density.jl`
- `src/cartesian_final_basis_realization/pqs_complete_core_shell_final_basis.jl`
- existing Ham/H1-J payload code reachable from `rg "complete_core_shell_ham|h1_j|density_interaction|raw_pair_factor"`
- `docs/src/developer/pqs_source_box_operator_framework.md`
- `docs/src/developer/pqs_source_box_fixture_policy.md`
- `~/Dropbox/chatarchive/handoff/software_packets/gausslet-methods-fundamentals.md`

Audit requirements:

1. Identify the current one-center H1-J diagnostic contract:
   - required inputs
   - density gauge
   - axis weights source
   - raw pair-factor source
   - what scalar it reports
   - why it is diagnostic/private rather than a full two-electron Hamiltonian

2. Identify what a downstream Be2 WL/PQS Hamiltonian comparison will actually
   need from GaussletBases:
   - final basis / retained basis description
   - H1 matrix
   - electron-electron representation or factor payload
   - nuclear metadata
   - ordering/convention labels
   - export/handoff constraints, if visible from current code

3. Compare those needs with what Be2/PQS now has:
   - private diatomic source plan
   - private final basis
   - private H1 payload
   - parent axis bundle / source-box raw-box data
   - Coulomb expansion access

4. Recommend the smallest safe next implementation pass:
   - H1-J diagnostic first, if it directly validates required two-electron
     conventions without adding much surface;
   - or electron-electron/Ham-input payload first, if H1-J would be a detour;
   - or a sharper blocked payload if required two-electron facts are not yet
     structured.

5. Propose blocker/status vocabulary for the next pass.

Trust boundary:

- No source edits.
- No commits.
- No H1-J/Ham/electron-electron materialization.
- No RHF/SCF/DIIS work.
- No WL payload, public API, exports, artifacts, hfdmrg execution, or CR2
  execution.
- Reading local hfdmrg/CR2 docs or code is allowed only if necessary for
  understanding downstream handoff shape; do not run downstream jobs.
- Do not promote shell/support-row contraction, raw product-box probes, or old
  WL adapter paths to route authority.
- Do not reinterpret retained diagnostic/self-integral weights as
  IDA/quadrature weights.
- Do not ask for interactive approval during unattended baton work. If approval
  would be required, write `.agent_handoffs/ATTENTION.md` and stop.

Validation:

- Read-only audit only.
- `git status --short --branch`
- No Julia commands are required unless you need a short local introspection
  probe. If you do run one, report why and keep it local/ignored.

Report back:

- H1-J diagnostic contract summary with file/line references.
- Downstream Ham-constructor needs as far as current code shows.
- What Be2/PQS already has versus what is missing.
- Recommendation for the next implementation pass.
- Proposed status/blocker labels.
- Git status.
- Deletion/shrinkage forecast:
  - deleted:
  - simplified:
  - quarantined:
  - not deleted because:
  - exact remaining caller/blocker:

-- repo-manager@macmini
