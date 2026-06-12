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
  decomposed WL unit-pair inventory and produce separated retained/global
  by-center matrices
- by-center records keep center identity separated and defer physical nuclear
  charge application to acceptance/Hamiltonian assembly
- decomposed route-global overlap and kinetic matrices can now be materialized
  from the same real unit-pair inventory source
- a narrow decomposed WL one-electron Hamiltonian assembly helper now combines
  route-global kinetic with separated unit-charge nuclear-attraction matrices;
  it applies recorded nuclear charges and sums centers only at Hamiltonian
  assembly
- a shellification-backed decomposed WL unit-pair inventory source now assigns
  retained dimensions and column ranges from retained-unit records, including
  the direct-core retained operator inventory
- therefore active H and H2+ scientific acceptance through the decomposed WL
  path now reach real one-electron solves without a full-window CPB or direct
  Cartesian fallback

The current q = 5, ns = 5 route metadata exposes terminal shellification unit
inventory at terminal-region granularity, and the local White-Lindsey adapter can
materialize overlap, kinetic, and one-center electron-nuclear by-center blocks
for a supplied decomposed unit pair. A compact
`white_lindsey_decomposed_unit_pair_inventory` validator now exists in the
pair-block materialization layer; for shellification-derived retained units it
assigns retained dimensions and column ranges from the White-Lindsey
complete-shell retained-count policy and exposes upper-triangular unit pairs
through a lightweight index table. The active He readiness audit now validates
the shellification-backed source as 27 decomposed units, 378 upper-triangular
unit pairs, retained/global dimension 223, and retained column coverage `1:223`.
The direct-core unit covers columns `1:125`; boundary units cover the shell
range `126:223`. The driver-facing helper
`white_lindsey_shellification_decomposed_unit_pair_inventory` now owns the
shellification -> lowering -> retained-unit -> lightweight unit-pair index table
-> decomposed-inventory handoff for this path. The route-global by-center
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
low-order decomposed seed inventory remains a one-shell bridge/oracle source; it
is not the production path for larger parent boxes.

The larger-box He readiness path now starts from shellification-derived units
instead of extending the low-order seed. With `AsinhMapping(c = 0.1, s = 1.0,
tail_spacing = 10.0)` and `parent_side_count = 13`, the reference endpoints are
`(-6, 6)`, the physical endpoints are about
`(-8.565228460168399, 8.565228460168399)` bohr, and the parent covers the target
radius `R = 6`. Shellification produces one direct-core region plus four
complete-shell regions. White-Lindsey lowering selects source CPB counts
`(1, 26, 26, 26, 26)`, yielding 105 retained units, 5,565 upper-triangular unit
pairs, and a 517-column retained inventory. The inventory source is
`:cartesian_shellification_retained_unit_pair_plan`; the old
`:white_lindsey_low_order_materialized_seed_ranges` source is kept only as the
compact one-shell acceptance bridge/oracle.

The corrected decomposed route materializes overlap, kinetic, one separated
uncharged electron-nuclear by-center matrix, the charge-applied one-electron
Hamiltonian, and the full retained density-density/IDA electron-electron
interaction matrix in the 223-column retained basis. The full interaction
matrix has shape `(223, 223)`. The density-density route follows the legacy IDA
transfer rule: contract `pgdg_intermediate.pair_factor_terms_raw` as the raw
pair numerator, project retained density weights through the same retained unit
coefficient maps, and divide the assembled retained numerator by the final
retained density-weight outer product. The weight-divided
`pgdg_intermediate.pair_factor_terms` table remains a PGDG storage relation and
local-oracle convenience, not the authority for transformed retained IDA
interactions. The resulting matrix is the object later correlation work should
consume. As a compact physics convention check before RHF, the lowest
one-electron He+ orbital for the same Z = 2 center has H1 eigenvalue about
`-1.878770102537269` Hartree versus the hydrogenic 1s reference `-2.0`. Using
that retained orbital, the IDA Coulomb self term is about
`1.4202542835594492` Hartree versus the hydrogenic 1s reference
`5Z/8 = 1.25`, an error of about `0.17025428355944916` on this compact
one-shell fixture. This self-Coulomb diagnostic is not accepted physics yet; it
is a convention probe before RHF interpretation.

The same ownership rule applies to projected-shell/PQS IDA paths: raw pair
numerators must be carried through the projection and Lowdin/final-basis
transformation first. Division by density weights belongs only at the final
retained/final-basis density-interaction boundary, using weights in that final
basis. Source-level or parent PGDG weight division is not a substitute for the
post-projection/PQS weight boundary.

For shellification-backed decomposed White-Lindsey inventories, the factorized
retained-basis backend is the fast production path for route-global overlap,
kinetic, electron-nuclear by-center, and density-density construction. The
pair-streaming route remains a reference/fallback path and is used for compact
equivalence checks. The current factorized retained-basis sidecar is derived
from shellification retained-unit coefficient maps through a temporary dense
retained coefficient-matrix bridge into the existing factorized extraction
kernel. That bridge is implementation scaffolding, not a change of authority:
old nested fixed-block matrices remain timing or oracle comparators only and
are not the source of the decomposed route matrices.

Restricted closed-shell HF with one alpha and one beta electron converges in
16 iterations. The bare closed-shell one-electron value from the lowest
one-electron orbital is `-3.757540205074538` Hartree. The self-consistent RHF
one-electron contribution is `-3.705023012192527` Hartree, the electron-electron
contribution is `1.3106054775285387` Hartree, and the accepted total HF energy
is `-2.3944175346639884` Hartree. The converged retained density has trace 1 for
the occupied spatial orbital, electron count 2 after closed-shell occupation,
peak retained-column density about `0.014711725584653477` at column 63, direct
core fraction about `0.7281451450072138`, shell/boundary fraction about
`0.27185485499278617`, and direct-core RMS radius about `0.5509019498518761`
bohr. The converged-density Coulomb contribution is positive and equals the RHF
electron-electron contribution, `1.3106054775285387` Hartree, under the current
full retained two-index density-density convention.

The larger-box shellification path is currently exploratory rather than a
routine acceptance gate. The side-13 probe with
`AsinhMapping(c = 0.1, s = 1.0, tail_spacing = 10.0)` reports source kind
`:cartesian_shellification_retained_unit_pair_plan`, retained dimension 517, 105
units, 5,565 upper-triangular unit pairs represented by
`CartesianUnitPairs.UnitPairIndexTable`, retained column coverage from `1:125`
through `517:517`, omitted large pair keys/summaries, and no low-order seed,
full-parent CPB, direct Cartesian, or ordinary IDA fallback. The same probe
materialized route-global overlap, kinetic, electron-nuclear by-center, and the
full retained density-density matrix. The lowest one-electron Z = 2 orbital
improved to about `-1.9748150892830352` Hartree versus the hydrogenic reference
`-2.0`. With raw pair numerators contracted and final retained density weights
applied after retained assembly, the retained IDA self-Coulomb term is about
`1.2158294767735702` Hartree versus `1.25`. The self-Coulomb error moves from
about `0.1703` in the compact fixture to about `0.0342` in the side-13 fixture,
so the corrected retained-weight boundary is the active convention. A follow-up
factorized-backend side-13 He RHF probe converged in 19 iterations with
one-electron contribution about `-3.8457802974697652` Hartree,
electron-electron contribution about `1.0092822977688516` Hartree, and total
RHF energy about `-2.8364979997009137` Hartree. The retained density trace was
1.0000000000000004 for the occupied spatial orbital, giving electron count
2.000000000000001 under the restricted closed-shell convention. This is an
exploratory larger-box physics result, not yet a permanent acceptance baseline.

A short fixed-`ns = 5` gausslet-only He ladder then found the current best
exploratory point at `d = 0.075`, `s = 0.75`, side count 17, retained dimension
713, and 157 retained units / 12,403 unit pairs. That point gives H1 error about
`0.002753`, IDA self-Coulomb error about `0.003389`, and RHF total energy about
`-2.858531351214` Hartree, about `+0.003149` Hartree above the He HF reference.
A larger diagnostic point at `d = 0.05`, `s = 0.5`, side count 23, retained
dimension 1007, improved the H1 and self-Coulomb diagnostics further but moved
the RHF total slightly away from the HF reference. The side-17 fixture is the
current best gausslet-only He exploratory point, but it is not a routine
acceptance test because cold runtime is still too long. Further gausslet-only He
convergence work should vary `ns`, `s`, and `d` together while keeping an
adequate box, roughly `R >= 8`, rather than overfitting fixed-`ns = 5` spacing.
The side-23 result is diagnostic only, not the chosen next fixture.

The Fig. 8 He RHF comparison line is now anchored to the archive-local
White-Lindsey table and provenance extract under
`/Users/srw/Library/CloudStorage/Dropbox/chatarchive/references/papers/canonical/`.
The relevant setup is He RHF with AHGBS-9 S-only supplement (`lmaxadd = 0`),
`doside = n_s`, `corespacing = d`, `dwidth = 10.0`, and
`gscalefac = sqrt(2.0)`. A repo old nested/QW MWG reproduction of the
`n_s = 5`, `d = 0.3` point loads AHGBS-9 from the external
`GaussletModules/BasisSets` file and matches the expected 447-function
structure: 419 gausslet functions plus 28 residual S directions. Its RHF energy
is about `-2.862102144533723` Hartree, while the Fig. 8 table gives
`-2.861543784624258` Hartree, a difference of about `-0.558 mHa`. This is close
enough for the current reproduction audit and is much smaller than the earlier
cc-pVTZ side13 diagnostic discrepancy. The `n_s = 5` Fig. 8 reproduction is
therefore recorded as structurally matched and scientifically useful, but not a
permanent acceptance test. The follow-up old nested/QW `n_s = 7` He RHF
diagnostic ran `d = 0.15`, `0.10`, and `0.20` with AHGBS-9 S-only. All three
were sub-mH relative to their Fig. 8 rows. In this repo probe, `d = 0.10` was
closest: final dimension 1897, gausslet/residual counts 1869/28, RHF total
`-2.861673961528321` Hartree, error `+1.716095156645281e-6` Hartree versus the
Fig. 8 row, and error `+6.034083917860755e-6` Hartree versus the He HF
reference. This is accurate enough to close the immediate atomic He accuracy
check.

The next atomic supplement target is Be with S+P GTO residuals. The old
nested/QW Be S+P oracle uses the external GaussletModules `BasisSets` file,
`legacy_atomic_gaussian_supplement("Be", "cc-pV5Z"; lmax = 1)`, `Z = 4`,
q/ns `5 / 5`, `interaction_treatment = :mwg`, and
`residual_keep_policy = :near_null_only`. It gives fixed dimension `615`,
`21` raw S+P supplement orbitals, `21` retained residual directions, final
dimension `636`, final overlap identity error about `1.01e-11`, and closed-
shell RHF total `-14.574514244574694` Hartree. This old route is now an oracle
comparator, not a target to keep re-proving.

The current decomposed/final-basis route has a private one-center
atom+supplement seam,
`_white_lindsey_decomposed_atom_gto_final_basis_route(...)`, that wires mapped
parent axes, shellification-backed decomposed WL inventory, combined-GTO
one-electron assembly, final-basis projection, and optional residual-MWG
density-density materialization. On the Be S+P q/ns `5 / 5` probe, the
one-electron final-basis route materializes with retained gausslet dimension
`615`, `21` retained supplement directions, final dimension `636`, final
overlap identity error about `1.01e-10`, and no full-parent CPB, direct
Cartesian, ordinary IDA, raw-GTO-final-density, or generalized-final-solve
fallback. Final density-density/RHF for Be S+P is not yet accepted; the next
step is a phase-attributed density-density/RHF run through the same seam rather
than a blind long run.

The current coarse timing split for the active tiny-box He RHF acceptance is
reported by the test as diagnostics, not asserted as performance thresholds.
Before the electron-nuclear cache and precompile workload, a representative
cold run was:

| phase | elapsed seconds |
| --- | ---: |
| parent seed report | 4.413 |
| parent-axis setup | 0.025 |
| decomposed inventory | 4.753 |
| route-global overlap | 26.606 |
| route-global kinetic | 0.659 |
| route-global electron-nuclear by center | 34.453 |
| one-electron Hamiltonian assembly | 0.067 |
| density-density matrix build | 0.713 |
| route operator and interaction build total | 62.498 |
| Hamiltonian/interactions build total | 67.276 |
| RHF solve | 0.689 |
| total acceptance elapsed | 74.117 |

The production route now carries TimeG scopes for decomposed WL inventory,
overlap, kinetic, one-body local batching, one-body global placement,
electron-nuclear by-center local batching, and electron-nuclear by-center
global placement. In a two-pass same-process probe after threading the
already-built decomposed inventory into the He route-global operator calls, the
cold pass took about 79.8 seconds total and the warm pass took about 15.0
seconds. The warm route-global overlap, kinetic, and density-density phases were
about `0.013`, `0.017`, and `0.040` seconds respectively. The warm by-center
electron-nuclear phase was still about `14.675` seconds, with TimeG attributing
that cost to `decomposed_wl.electron_nuclear_by_center.local_batch`; global
placement was about `0.29` seconds in the cold pass and below the live threshold
in the warm pass. Threading the inventory did not materially change the cold
one-electron timings, so the first redundant inventory rebuild was not the main
cost center. The remaining slowness is local by-center nuclear block
construction; the large cold overlap cost is primarily compilation.

The older flat/non-module WL path in this repo is the right comparison target
for the next optimization pass. The corresponding optimized functions are
`_qwrg_diatomic_overlap_matrix`,
`_qwrg_diatomic_kinetic_matrix`, `_qwrg_diatomic_nuclear_one_body_by_center`,
and the direct/staged by-center nuclear helpers in `ordinary_qw_raw_blocks.jl`,
with orchestration in `_ordinary_cartesian_qiu_white_operators_pure_bond_aligned_direct`.
That code batches work by reusing parent 1D bundles, filling product matrices
directly for overlap/kinetic, precomputing per-axis Gaussian nuclear term
tables once per unique center, and then contracting/filling blocks without
reconstructing route inventories or local placement plans per operator. The
next concrete optimization target is to make the decomposed/module
electron-nuclear by-center route reuse per-axis Gaussian term tables and local
unit-pair coefficient work across all by-center block construction, following
the older flat WL batching pattern, rather than reconstructing those local
inputs per decomposed unit pair.

The decomposed WL electron-nuclear by-center route now caches centered PGDG
Gaussian axis-term tables once per center/axis/expansion/source context and
reuses them across all decomposed unit pairs. A post-cache compile-attribution
probe measured the active electron-nuclear local batch at about `0.072` seconds
warm and the full route-global electron-nuclear by-center call at about `0.084`
seconds warm. The avoided diagnostic repeated-axis setup remains about 14 to 17
seconds on the q/ns = 5/5 He fixture, confirming that the cache targets the
previous warm bottleneck.

A narrow `PrecompileTools` workload now compiles only the existing decomposed WL
one-body production route calls for the q/ns = 5/5 seed: route-global overlap,
kinetic, and electron-nuclear by-center. No route behavior, fallback,
Hamiltonian solve, acceptance assertion, GTO/PQS path, export, or artifact
logic lives in the precompile workload. After a content-changing edit to the
workload, package precompilation took about `39.5` seconds. A cached fresh
process `using GaussletBases` measured about `0.67` seconds. Fresh-process route
calls then landed at warm-scale timings for the seed-backed one-shell
route-global operator calls: route-global overlap about `0.012` seconds cold,
route-global kinetic about `0.016` seconds cold, and route-global
electron-nuclear by-center about `0.089` seconds cold. After redirecting the
active He acceptance route to shellification-derived retained units, the He RHF
energy remained within the same regression window at about
`-2.045516767078339` Hartree. A representative fresh test run spent about
`9.0` seconds in shellification retained-unit/lightweight-pair inventory
construction, about `6.6` seconds in one-electron operator build, about `0.70`
seconds in density-density matrix build, and about `0.78` seconds in the RHF
solve. The side-13 inventory probe improved from about `77.9` seconds to about
`45.9` seconds after replacing the stored rich `UnitPairRecord` tuple with the
module-owned `CartesianUnitPairs.UnitPairIndexTable` and omitting duplicate pair
summary materialization. After the WL-only selected-contract path and one-pass
retained-unit planning cleanup, the same standalone side-13 inventory probe
reported about `27.5` seconds while preserving 105 units, 5,565 pairs, and
retained dimension 517. A same-process phase probe showed the warm planning path
is no longer the main algorithmic cost; remaining cold-process time is dominated
by first-call compilation and retained-range wrapping rather than rich pair
tuple storage. The shellification-backed one-body assembly path now streams over
`CartesianUnitPairs.UnitPairIndexTable`, materializes each local WL block only
when needed, and inserts it immediately into the retained/global matrix. On the
side-13 fixture, a cold probe materialized 5,565 local blocks for each of
overlap, kinetic, and one by-center electron-nuclear matrix with same-process
warm route-global times of about `0.092`, `0.218`, and `1.179` seconds,
respectively, for the 517-column retained basis.

The shellification-backed density-density assembly path now also streams over
`CartesianUnitPairs.UnitPairIndexTable`. It caches retained-unit coefficient maps
once per unit, builds one retained pair block at a time, and inserts directly
into the retained two-index interaction matrix instead of materializing a full
pair-coefficients batch first. On the accepted one-shell He fixture, the
streaming and preserved rich-pair batch paths agreed exactly with maximum matrix
difference `0.0`. The active He RHF energy stayed at about
`-2.045516767078339` Hartree, with the density-density build about `0.75`
seconds in a fresh acceptance run. On the side-13 probe, the streaming
density-density route materialized 5,565 unit pairs for the 517-column retained
basis in about `17.5` seconds in a cold ordered probe, with the scoped local
pair-stream phase about `1.08` seconds. A direct side-13 comparison against the
preserved rich-pair batch route agreed exactly with maximum matrix difference
`0.0`; forcing the old rich-pair inventory view cost about `10.4` seconds, and
second same-process density-density calls were essentially tied at about `1.04`
seconds for streaming versus `1.03` seconds for the preserved batch route. The
streaming path removes the full pair-coefficients result array and rich-pair
conversion pressure, but the local support-block contraction is now the shared
warm cost center rather than a placement-sidecar problem.

One supported exploratory probe with the same one-shell decomposed topology and
finer Z = 2 spacing, `d = 0.15`, shrinks the physical endpoints to about
`(-0.6557127550383339, 0.6557127550383339)` bohr and worsens the RHF energy to
about `0.0633231599983839` Hartree. That probe took about 79 seconds and is not
promoted into the active test. The current evidence points first to low-order
box/basis quality, especially the very compact one-shell physical extent, not a
failure of the closed-shell density-density scalar convention. The larger-box
readiness probe now exercises shellification-derived decomposed unit inventory;
the next larger-box milestone is a complete RHF run through the now-streaming
one-body and density-density route, not another seed inventory extension.

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
