# Current ordinary branch

PQS is the standard current Cartesian construction and default base route.
Supported interfaces and experimental backend choices should be distinguished
from scientific convergence or general molecular solver support.

## What the ordinary branch is today

The [PQS manual](../manual/projected_q_shells.md) is the public starting point.
`cartesian_base_hamiltonian` builds an unsupplemented, uncorrected all-electron
localized IDA Hamiltonian for origin-centered positive-charge atoms and
equal-label/equal-charge homonuclear z-axis diatomics. Both `q` and `ns` are
supported: PQS uses `q = ns`; White-Lindsey uses `q = ns - 2`. Supplying both
requires consistency. White-Lindsey remains a matched comparison/alternative.

Separate released guides describe:

- [Reference-density Hartree screening](../manual/reference_density_hartree_screening.md):
  supplied-field correction assembly in one orthonormal basis/order with a
  separate energy scalar, not field generation, fitting, exchange, or SCF.
- [External Cartesian GTO transfer](../manual/external_cartesian_gto_transfer.md):
  validated version-1 packets, raw projection with capture loss, and optional
  caller-thresholded determinant cleanup. The checkpoint-only PySCF/NumPy
  exporter is optional and supplies neither a Hamiltonian nor a full restart;
  it does not broaden destination-basis geometry.

The lower-level ordinary/Qiu-White routes also provide:

- Coulomb-expansion first, not 3D grid first
- one-dimensional mapped ordinary bases on each Cartesian axis
- exact Cartesian workflow primitives:
  - `cross_overlap`
  - `basis_projector`
  - `transfer_orbitals`
- an explicit backend split between validation and experimental PGDG-style
  analytic construction
- a separate paper-faithful Qiu-White residual-Gaussian route
- real one-center nested Cartesian fixed-block / consumer support
- real bond-aligned diatomic nested source / fixed-block / diagnostics /
  geometry payload support
- the old COMX-cleaned 1D hybrid route quarantined as legacy/internal code,
  not current workflow

The exported bond-aligned diatomic source, fixed-block, and diagnostics entry
points remain supported within their contracts. Experimental backend and
geometry policies below do not broaden the public base facade.

## Start here

Within the rendered docs, the best supporting pages are:

- [Example guide](../howto/example_guide.md)
- [Atomic and ordinary workflows reference](../reference/atomic_and_ordinary.md)
- [Qiu-White residual-Gaussian route](../algorithms/qiu_white_residual_gaussian_route.md)
- [Cartesian nested atomic nonrecursive route](../algorithms/cartesian_nested_atomic_nonrecursive_route.md)
- [Cartesian nested diatomic box policy](../algorithms/cartesian_nested_diatomic_box_policy.md)
- [Cartesian nested endcap/panel shared-shell route](../algorithms/cartesian_nested_endcap_panel_shared_shell.md)

The radial tutorial remains useful context because the package’s operator and
quadrature story was established there first, but it is no longer the whole
ordinary-branch story.

## Documentation authority for this branch

For the current ordinary / Cartesian story, use these pages in this order:

- this page for branch status and scope
- [Reference](../reference/index.md) for exported API surface
- [Algorithms](../algorithms/index.md) for path-specific construction recipes
- [Developer Notes](../developer/index.md) for lower-priority architecture and
  supporting material

Older flat `docs/*.md` notes remain available as supporting or historical
material, but they should not override this page as the current ordinary branch
summary.

## Current interpretation

The wording discipline for the ordinary line remains:

- `:numerical_reference` is the validation route
- the PGDG-style analytic route is good enough on the mapped ordinary backbone
- on that PGDG lane, distorted 1D primitives are replaced by plain Gaussian
  proxies before matrix assembly and numerical primitive quadrature is outside
  the intended production contract
- the public mapped ordinary backbone is the PGDG-capable surface today:
  `mapped_ordinary_one_body_operators` and the mapped Cartesian IDA path may
  use `:pgdg_experimental` or `:pgdg_localized_experimental`
- the pure bond-aligned ordinary direct-product Qiu-White route is also
  PGDG-capable on `gausslet_backend = :pgdg_localized_experimental`
- the pure prebuilt nested fixed-block ordinary/Qiu-White operator route is
  also PGDG-capable on `gausslet_backend = :pgdg_localized_experimental` when
  the nested fixed block still represents a pure Cartesian parent space
- the bond-aligned diatomic molecular supplement direct-product and prebuilt
  nested fixed-block Qiu-White routes are now also PGDG-capable on
  `gausslet_backend = :pgdg_localized_experimental`
- on those widened pure PGDG-capable routes, molecular nuclear one-body terms
  are intended to contract directly into the exposed parent/final space rather
  than first materializing avoidable dense parent 3D nuclear matrices
- on the widened diatomic molecular supplement route, only the carried GG/fixed
  nuclear one-body backbone is widened; GA/AA supplement blocks and
  residual-Gaussian closure remain on the existing formulation
- bond-aligned diatomic nested frontends also expose an opt-in experimental
  `shared_shell_layer_policy = :endcap_panel_owned` route; it is mainline
  code, but remains disabled by default and should be compared by consumers
  before being treated as production policy; this opt-in path may use
  `gausslet_backend = :pgdg_localized_experimental`
- on that same diatomic molecular supplement route, the final one-body mix is
  now intended to use the carried-plus-residual `raw_to_final` structure
  directly rather than a generic dense final-space congruence
- the Cr-facing one-center atomic full-parent fixed-block and atomic
  supplement-bearing Qiu-White routes now resolve omitted
  `gausslet_backend` through `:auto` to `:pgdg_localized_experimental`; the
  nested operator overload reuses recorded fixed-block backend provenance, and
  refuses `:auto` for older/hand-built fixed blocks whose backend provenance is
  `:unknown`
- remaining nested source-building front doors and other supplement-bearing
  Qiu-White routes remain `gausslet_backend = :numerical_reference` unless a
  route is explicitly documented otherwise; that includes the default
  complete-rectangular diatomic nested source builder and the experimental
  chain/square nested QW source wrappers
- exact overlap / projector / transfer are now first-class workflow primitives
- Gaussian-supplement comparisons should go through the separate paper-faithful
  3D Qiu-White route, not the old 1D COMX-cleaned surrogate path
- the public one-center White-Lindsey-style mapping language is now the
  `d`-driven helper documented in [Bases and mappings](../reference/bases_and_mappings.md),
  not the older fixed-`a` / `count -> s` historical note language
- hard pure mapped small-`c` cases remain stress tests
- `AsinhMapping` is the current working map, not final truth

The ordinary branch also includes narrow experimental producer-side chain and
square-lattice lines. Those are real code, but they should still be read as
experimental rather than as settled broad workflows.
