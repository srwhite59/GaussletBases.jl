Pass 225 - audit H2 PQS MWG/GTO supplement provider seam

Role:
You are `repo-doer@macmini` doing one bounded no-edit audit for GaussletBases.
Follow `AGENTS.md`, `JuliaStyle.md`, `BlurbStyle.md`, and the unattended baton
rules in `.agent_handoffs/RUN.md`.

Loop/approval rule:
- Unattended baton mode is active.
- Do not ask the user for permission through an escalation prompt.
- If a needed command requires approval, write `.agent_handoffs/ATTENTION.md`
  with the exact command, why it is needed, and the blocking condition, then
  stop.

Current state:
- Head before this pass should be the commit that accepts pass 224.
- H2 R=4 q=n_s=5 gausslet-only PQS/WL physical endpoint is accepted and
  comparison-ready at final dimension 463.
- Pass 224 added a private supplement preflight boundary. For
  `supplement_policy = :mwg_residual_gto`, the first blocker is:

```text
:missing_provider_gto_supplement_blocks
```

Architectural guardrail:
WL and PQS should share the same physical support/intermediate gausslet plan:

```text
support counts  = (275, 578, 362)
retained counts = (251, 98, 114)
retained order  = (:atom_contact_core, :shared_shell_1, :shared_shell_2)
gausslet final dimension = 463
```

WL and PQS differ in retained transform and fast operator construction. GTO/MWG
supplement policy should sit above that distinction, not become a PQS-specific
physics fork.

Task:
No source/test/bin edits. Audit the existing WL/MWG/GTO supplement machinery and
report the smallest seam that can make the first PQS supplement preflight fact
available:

```text
:provider_gto_supplement_blocks
```

Read/search for existing code that provides or consumes:

```text
raw GTO supplement basis/request objects
gausslet/GTO mixed overlap and one-body blocks
GTO/GTO overlap and one-body blocks
MWG residualization / residual Gaussian orthogonalization
combined gausslet+GTO transform construction
density-density or two-body supplement blocks
old H2 WL/QW supplemented reference path
```

Likely search terms:

```text
residual
MWG
GTO
supplement
cc-pVTZ
bond_aligned
mixed
gaussian
```

Questions to answer:
1. Which existing files/functions own the WL/MWG supplement provider logic?
2. Which parts are route-independent enough to reuse for PQS?
3. What exact data does the H2 PQS physical route already have that those
   providers need: parent basis, axis bundles, final-basis transform, center
   metadata, nuclear charges, support plan, retained order?
4. What is still missing before `:provider_gto_supplement_blocks` can become
   available?
5. Is the first implementation pass likely to be:
   - a supplement request payload,
   - a provider-block adapter,
   - a mixed-block convention wrapper,
   - or a blocked payload with sharper missing facts?
6. Which stale source/test surfaces could be deleted in the first coding pass to
   keep line count negative?

Do not:
- implement provider blocks;
- build GTO/GTO, mixed gausslet/GTO, or MWG residual matrices;
- add H2 supplemented scalar values;
- change the accepted no-supplement H2 endpoint;
- touch CR2/HFDMRG/export/HamV6/public API behavior;
- edit source, tests, bin scripts, or docs besides the response file.

Validation:
No Julia validation is required for this no-edit audit. If you run only read-only
commands such as `rg`, `sed`, `git status`, or `git log`, report them.

Response file:
Write `.agent_handoffs/response.225.md` and also copy it to:

```text
docs/src/developer/blurb_logs/2026-06-12_pqs_source_box_pivot/response.225.md
```

Report:
- files/functions audited;
- existing WL/MWG/GTO supplement provider surfaces;
- reusable pieces versus WL-specific pieces;
- exact PQS data already available;
- exact missing facts for `:provider_gto_supplement_blocks`;
- smallest recommended implementation pass;
- candidate deletion/shrink surfaces for that implementation pass.

Continue polling after writing the response. Manager will review, commit/push,
and publish the next blurb unless a real blocker appears.

-- repo-manager@macmini
