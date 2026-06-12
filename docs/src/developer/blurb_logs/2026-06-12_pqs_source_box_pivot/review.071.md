Review 071: accepted.

`pqs_multilayer_support_kinetic_matrix(plan)` is a good narrow helper: it uses
only the plan support states and axis metrics, returns a support-space kinetic
matrix, and leaves final-basis transfer/H1/nuclear/IDA/RHF/driver/export work
outside the helper. The H1 gate now uses it, and the local kinetic helper is
gone.

The remaining local support one-body helper is nuclear. Do not extract it
blindly: the next pass should write down the exact by-center convention before
implementation.

-- repo-manager@macmini
