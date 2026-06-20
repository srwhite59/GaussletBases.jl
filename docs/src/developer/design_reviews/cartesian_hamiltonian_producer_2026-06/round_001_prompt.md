# Round 001 Review Prompt

Review only:

- `docs/src/developer/cartesian_hamiltonian_producer_design.md`
- files under this review directory

Do not edit:

- `src/`
- `test/`
- `bin/`
- `tools/`
- `AGENTS.md`

The current design document is a draft and is not implementation authority. The
goal of this round is to improve the design before any source implementation
starts.

## Review Hats

Review with three explicit hats.

### Numerical Correctness

Check:

- shell projection and shell-local Lowdin sequence;
- no global core/shell or PQS/RG Lowdin;
- direct terminal sectors as identity blocks;
- cross-block overlap acceptance and stop condition;
- center-specific unit nuclear attraction convention;
- localized IDA gauge and final-column sign handling;
- whether pseudocode is scientifically wrong or underspecified.

### Performance And Memory

Check:

- asymptotic cost by slice;
- Cr2-scale temporary dense matrices;
- whether the design avoids global support operators and unnecessary global
  coefficient matrices;
- whether blockwise contraction strategy is realistic;
- whether line budgets and validation gates are plausible.

### Minimality And Deletion

Check:

- whether each proposed HP-* surface has a live physics consumer;
- which proposed objects/functions can be removed or merged;
- which fields are redundant or derived;
- exact old blocked preflight/payload/report/test surfaces each slice should
  delete;
- whether any item creates a new framework instead of replacing an old one.

## Required Output

Write a review in `round_001_consolidated_review.md` using the template in
`round_001_review_template.md`.

The review must include:

- blocking issues;
- candidate document edits, preferably as patch-style snippets or precise
  replacement text;
- surfaces to remove from the registry;
- surfaces that should stay candidate rather than approved;
- missing invariants or forbidden paths;
- questions that must be settled before ChatGPT-Pro review.

Do not request implementation work. If a source change seems necessary, express
it as a future design-slice requirement, not a code task.
