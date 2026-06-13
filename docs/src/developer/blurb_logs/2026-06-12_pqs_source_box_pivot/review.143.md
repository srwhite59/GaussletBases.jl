Pass 143 manager review

Accepted.

Plain-language state:

- The Be2 raw-box payload has enough raw ingredients to define a plausible
  source realization:
  product/doside unit as the core/body sector, and left/right PQS raw plans as
  shell/source sectors.
- Counts line up cleanly: product 25 plus left/right 98 retained modes gives
  retained dimension 221.
- The mapping is not the old one-center `:pqs_multilayer_shell_source_plan`
  vocabulary. It needs an explicit diatomic source-realization payload and a
  later adapter/new-object decision.

Decision:

- Add a private diatomic complete-core/shell source-realization payload next.
- It should record support order, retained order, range/permutation facts,
  support counts, retained counts, expected shell coefficient block structure,
  and the explicit nonclaim that it is not a
  `:pqs_multilayer_shell_source_plan`.
- Do not build the actual source plan, final basis, H1, H1/J, Ham data, RHF,
  exports, or artifacts.

Deletion/shrinkage assessment:

- deleted: none.
- simplified: the next pass can move source-realization ordering/count/shape
  facts into one private payload instead of keeping them spread across raw-box,
  support-window, and source-plan summaries.
- quarantined: returning `:pqs_multilayer_shell_source_plan`, final-basis/H1
  consumption, support-row contraction authority, retained diagnostic weights as
  IDA weights, RHF/SCF, WL payloads, exports, artifacts, hfdmrg, and CR2
  execution remain out of scope.
- not deleted because: existing one-center final-basis/H1 consumers remain the
  validated consumer shape; Be2 does not yet have an approved adapter or new
  consumer.
- exact remaining caller/blocker: Be2/PQS needs a private source-realization
  payload plus adapter/new-object contract before existing complete-core/shell
  final-basis/H1 callers can consume it honestly.

-- repo-manager@macmini
