# Numerical Contracts

This page records a small internal engineering policy that is easy for machine
generated code to miss.

It is intentionally developer-facing rather than part of the main user manual.

## Orthonormal blocks

When a block is constructed to be orthonormal, the intended contract is:

1. build it to be orthonormal in the relevant metric
2. check that the resulting overlap is `I +` small unavoidable Float64 noise
3. then treat the block as orthonormal

Examples include:

- finalized PGDG / COMX-cleaned blocks
- nested fixed blocks
- any other internal block whose construction is explicitly meant to produce an
  orthonormal basis

## What not to do

Do not keep propagating a near-identity overlap matrix as if it were meaningful
mathematical data.

In particular:

- do not store `S = I + eps` by default when `eps` is just Float64 residue
- do not build downstream logic that keeps consulting such an `S`
- do not interpret tiny nonorthogonality at the `1e-12` to `1e-14` level as a
  real structural feature

## What to do instead

Use overlap matrices in this regime only for:

- construction
- validation
- assertions
- or one final cleanup step if the deviation from identity is not small enough

For transfers between final orthonormal working bases:

- use only the cross overlap between the two final bases
- treat self-overlaps of those final bases as diagnostic-only
- do not turn the final working path into a generalized-overlap formulation

For decomposed White--Lindsey plus GTO supplement acceptance, the raw combined
gausslet+GTO Galerkin generalized solve is diagnostic-only. The active final
path first projects the GTO supplement into the residual space against the
decomposed WL gausslet sector, orthonormalizes the retained residual directions,
and then solves an ordinary Hermitian problem in the resulting final working
basis. The raw combined overlap and supplement self-overlap are construction
and validation data, not normal downstream working-basis data.

After that, the working representation should simply regard the block as having
identity overlap.

This is a coding and design policy, not a user-facing scientific statement.

## Nested fixed-block kinetic

`_NestedFixedBlock3D.kinetic` follows the nested packet contract.

That means:

- it is the kinetic matrix carried by the assembled nested packet
- it is the kinetic payload that downstream nested operator routes should use
- it is not automatically interchangeable with "contract the ordinary parent
  kinetic later and call that the same thing"

For current nested diatomic routes, this distinction is real. A nested
fixed-block kinetic can differ measurably from a later contraction of a
separately assembled ordinary parent one-body path even when both live on the
same final basis dimension.

## One-body reassembly

`assembled_one_body_hamiltonian(...)` reassembles from the operator payload that
was actually stored:

- stored `kinetic_one_body`
- stored `nuclear_one_body_by_center`
- requested `nuclear_charges`

So the meaningful contract is:

- compare reassembled one-body matrices against other operators built from the
  same stored kinetic contract

Do not use `assembled_one_body_hamiltonian(...)` to assert equality between two
routes that already disagree about what the kinetic payload is supposed to be.

In particular, for nested fixed-block routes:

- a by-center payload should be compared against another by-center payload on
  the same nested kinetic contract
- not against a separate total-only route that rebuilds one-body terms through a
  different parent-space contraction path

## Decomposed WL Gausslet-Only Acceptance Boundary

The intended gausslet-only scientific acceptance path is a true decomposed
White-Lindsey calculation with `q = 5` and `ns = 5`, not a single CPB covering
the full parent product window. The active readiness check is
`test/nested/cartesian_wl_gausslet_h_atom_acceptance_runtests.jl`.

Current status:

- decomposed White-Lindsey overlap and kinetic local pair-block paths exist
- global pilots exist for those supported safe one-body terms
- CPB-local electron-nuclear by-center provider blocks exist
- decomposed White-Lindsey electron-nuclear by-center local block support now
  exists for one center
- placement-plan and global by-center matrix pilots now accept
  `:electron_nuclear_by_center` records
- the route-global electron-nuclear by-center adapter can now consume the real
  decomposed low-order WL seed unit-pair inventory and produce separated
  retained/global by-center matrices
- by-center records keep center identity separated and defer physical nuclear
  charge application to acceptance/Hamiltonian assembly
- decomposed route-global overlap and kinetic matrices can now be materialized
  from the same real unit-pair inventory source
- a narrow decomposed WL one-electron Hamiltonian assembly helper now combines
  route-global kinetic with separated unit-charge nuclear-attraction matrices;
  it applies recorded nuclear charges and sums centers only at Hamiltonian
  assembly
- a decomposed WL unit-pair inventory source is now exposed from the
  materialized low-order seed retained ranges, including the direct-core
  retained operator inventory
- therefore active H and H2+ scientific acceptance through the decomposed WL
  path now reach real one-electron solves without a full-window CPB or direct
  Cartesian fallback

The current q = 5, ns = 5 route metadata exposes terminal shellification unit
inventory at terminal-region granularity, and the local White-Lindsey adapter can
materialize overlap, kinetic, and one-center electron-nuclear by-center blocks
for a supplied decomposed unit pair. The terminal shellification summary still
marks its own pair inventory as deferred, but the materialized low-order seed now
provides a narrow decomposed inventory source through its direct-core, face,
edge, and corner retained ranges. A compact
`white_lindsey_decomposed_unit_pair_inventory` validator now exists in the
pair-block materialization layer; it accepts a `UnitPairPlan` or unit-pair
records when they carry retained dimensions and column ranges, and reports
compact pair/range/global-dimension metadata. The active readiness audit now
validates the seed-backed source as 27 decomposed units, 378 upper-triangular
unit pairs, retained/global dimension 223, and retained column coverage
`1:223`. The direct-core unit covers columns `1:125`; boundary units cover the
shell range `126:223`. The route-global by-center
nuclear adapter uses that inventory plus the existing local
`electron_nuclear_by_center` block path to materialize one uncharged
retained/global matrix per supplied center. A focused one-center fingerprint
currently materializes all 378 decomposed local pair blocks into a 223 by 223
retained matrix. Centers are not summed and nuclear charges are recorded but not
applied in the by-center matrix path. The decomposed Hamiltonian helper consumes
the unit-charge nuclear-attraction convention already carried by those matrices
and multiplies by the recorded nuclear charge at Hamiltonian assembly.

The current H atom audit materializes decomposed route-global overlap, kinetic,
one-center electron-nuclear by-center, and the one-electron Hamiltonian for
`q = 5`, `ns = 5`, retained dimension 223. The decomposed unit inventory spans
retained columns `1:223`, so the assembled overlap is full rank with no missing
prefix columns. The current diagnostic values are: minimum overlap eigenvalue
about `0.999999999999839`, maximum about `1.000000000000165`, condition
estimate about `1.000000000000327`, symmetry error about `2.7e-17`, zero
near-zero eigenvalues, and rank estimate 223. The one-electron H solve uses the
ordinary symmetric path and gives `-0.4788666674548281` Hartree, which is
variational relative to the exact `-0.5` Hartree value.

The current H2+ audit uses the same decomposed route at `R = 2.0` bohr with
protons at `z = +/-1.0`. It materializes two separated uncharged by-center
nuclear matrices, applies both unit charges only in Hamiltonian assembly, and
uses the same full-rank retained overlap metric. The electronic energy is
`-1.033841728044377` Hartree. With nuclear repulsion `1/R = 0.5`, the
Born-Oppenheimer total energy is `-0.533841728044377` Hartree. This is
variational relative to the exact total reference near `-0.6026342144949465`.

The H2+ plus GTO supplement acceptance uses the same decomposed WL route and
final-basis GTO residual projection as the H plus GTO fixture. For `R = 2.0`
bohr with H cc-pVTZ `lmax = 0` supplements on both centers, the combined
dimension is 229, all six raw supplement residual directions are retained, and
the final overlap identity error is about `3.0953017926549364e-13`. The final
ordinary symmetric electronic energy is `-1.0988624733888488` Hartree and the
Born-Oppenheimer total energy is `-0.5988624733888488` Hartree. The raw
generalized combined solve is diagnostic-only and is not the active final-basis
acceptance result.

The first He atom decomposed WL acceptance uses the q/ns = 5/5 gausslet-only
fixture with one center at the origin and `Z = 2`. The active fixture follows
the one-center WL atomic spacing rule: `d = 0.2`, `s = sqrt(d Z) =
0.6324555320336759`, and `tail_spacing = 10.0`. An earlier readiness checkpoint
mixed the standard Z = 2 decomposed seed inventory with an H-style shared
adapter axis source; that violated the Z-dependent spacing contract and is not
the accepted He baseline.

The corrected active fixture is still a deliberately tiny low-order box. The
mapped one-dimensional parent has 7 centers with reference endpoints `(-3, 3)`
and physical endpoints about `(-0.9666200087560217, 0.9666200087560217)` bohr.
The minimum adjacent physical spacing is about `0.20858672857920835`, the
maximum adjacent spacing is about `0.4698383534446874`, the direct-core retained
range is `1:125`, and the shell retained range is `126:223`. The current
low-order decomposed seed inventory only accepts the one-shell fixture; a
`parent_side_count = 9` exploratory probe is blocked by the existing
single-shell inventory contract rather than by the RHF convention.

The corrected decomposed route materializes overlap, kinetic, one separated
uncharged electron-nuclear by-center matrix, the charge-applied one-electron
Hamiltonian, and the full retained density-density/IDA electron-electron
interaction matrix in the 223-column retained basis. The full interaction
matrix has shape `(223, 223)`, uses the existing WL pair-factor-term convention
with integral weights deferred to the IDA/HF density interpretation stage, and
is the object later correlation work should consume. Restricted closed-shell HF
with one alpha and one beta electron converges in 17 iterations. The bare
closed-shell one-electron value from the lowest one-electron orbital is
`-3.7575402050745312` Hartree. The self-consistent RHF one-electron contribution
is `-3.7316519035708953` Hartree, the electron-electron contribution is
`1.6861351364925603` Hartree, and the accepted total HF energy is
`-2.045516767078335` Hartree. The converged retained density has trace 1 for
the occupied spatial orbital, electron count 2 after closed-shell occupation,
peak retained-column density about `0.016237877162231747` at column 63, direct
core fraction about `0.7661258457949129`, shell/boundary fraction about
`0.23387415420508706`, and direct-core RMS radius about `0.544451699989865`
bohr. The converged-density Coulomb contribution is positive and equals the RHF
electron-electron contribution, `1.6861351364925603` Hartree, under the current
full retained two-index density-density convention.

One supported exploratory probe with the same one-shell decomposed topology and
finer Z = 2 spacing, `d = 0.15`, shrinks the physical endpoints to about
`(-0.6557127550383339, 0.6557127550383339)` bohr and worsens the RHF energy to
about `0.0633231599983839` Hartree. That probe took about 79 seconds and is not
promoted into the active test. The current evidence points first to low-order
box/basis quality, especially the very compact one-shell physical extent, not a
failure of the closed-shell density-density scalar convention. A larger-box He
acceptance fixture needs a reviewed decomposed inventory that supports more than
one shell before it should replace this baseline.

The He RHF energy is above the He HF reference near `-2.861679995612234`
Hartree. The acceptance path does not use a full-parent CPB, direct Cartesian
fallback, `ordinary_cartesian_ida_operators`, a generalized final solve, GTO
supplements, PQS transforms, exports, or artifacts.

Do not use the existing nested fixed-block operator matrices as the acceptance
path. They remain useful historical/oracle material, but they bypass the
decomposed retained-unit pair inventory that the active WL scientific tests are
intended to protect.

The retired transition helpers built one full-parent CPB with role
`:wl_cpb_acceptance_full_parent_window`. That path exercised CPB-local
operators but not White-Lindsey boundary-unit decomposition, so it is not an
active acceptance contract. Its observed values are retained only as transition
notes:

- H full-window CPB energy `-0.4832079279118124` Hartree
- H2+ full-window CPB electronic energy `-1.0971828374927926` Hartree
- H2+ full-window CPB total energy `-0.5971828374927926` Hartree

Older direct-route transition baselines remain historical comparison points
only:

- H small fitted direct route `-0.4706400351534759` Hartree
- H coarse/distorted direct route `-0.4966106635473884` Hartree
- H2+ R = 2.0 direct electronic energy `-1.0654839328172023` Hartree
- H2+ R = 2.0 direct total energy `-0.5654839328172023` Hartree

The next implementation after these gausslet-only H/H2+ acceptance baselines is
not a full-window fallback; it is either a reviewed accuracy improvement for the
decomposed WL fixture or a separate, explicitly scoped GTO/PQS extension. Do not
reintroduce the full-parent CPB helper or a direct Cartesian product fallback as
the active route.
