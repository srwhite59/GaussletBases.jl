# One-Center Cartesian Nested Fixed-Block Atomic Contract Correction Note

This note supersedes the earlier `ns7/ns9` one-center localization diagnosis.

The earlier pass used the wrong object for the atomic case:

- it built a central inferred working box
- it did not cover the full parent lattice
- it therefore diagnosed a contaminated fixed-block backbone

For the intended one-center atomic line, the correct contract is:

- start from the original largest parent lattice
- keep full parent coverage
- use `working_box = (1:n, 1:n, 1:n)`
- peel complete shells until the direct inner cube is about `ns^3`
- transfer that final inner cube with the identity

That full-parent shell-sequence route is the right backbone for the atomic
anchor cases.

The repo now hardens that contract through the canonical helpers:

- `build_one_center_atomic_full_parent_shell_sequence(...)`
- `one_center_atomic_full_parent_fixed_block(...)`

## Scope

This corrected comparison stays on the same narrow atomic anchor pair:

- `ns7_family`, `Z = 2`, `rmax = 10`
- `ns9_family`, `Z = 2`, `rmax = 10`

and keeps the intended narrow execution line:

- repo-local depot only
- no helper include from the old heavy pass
- no KrylovKit
- no parent eigensolve
- no `$HOME/.julia` fallback

So this correction re-evaluates:

- fixed-block dimension
- coverage / ownership structure
- low fixed one-body `H1`

It does not recompute `1s` / `2p` capture in this pass, because that would
reintroduce the parent eigensolve path the narrowed contract was meant to
avoid.

## Corrected Atomic Fixture

The corrected one-center atomic fixture is:

- `mapping = white_lindsey_atomic_mapping(Z = Z, d = d, tail_spacing = 10.0)`
- resolve the parent mapped 1D count so `rmax = 10` is actually covered
- build the mapped 1D gausslet bundle
- wrap it into the 3-axis nested bundle
- build the shell sequence on the full parent box
  - `working_box = (1:n, 1:n, 1:n)`
  - `nside = ns`
- build the fixed block from that full-parent shell sequence

For this correction pass, the reproducible helper is:

- `tmp/work/diagnose_one_center_atomic_full_parent_contract.jl`

That local helper now keeps the old central-box construction only as an
explicitly named wrong comparison object:

- `wrong_central_box_atomic_fixture`

and it should not be reused for one-center atomic contraction diagnosis.

## Anchor Results

### `ns7_family`

Old contaminated central-box fixture:

- `count = 19`
- `fixed_dim = 517`
- `nshells = 4`
- `direct_core_len = 5`
- `support_count = 2197`
- `expected_support_count = 6859`
- `missing_row_count = 4662`
- `working_box = (4:16, 4:16, 4:16)`
- `full_parent_working_box = false`
- `low_fixed_h1 = [-1.998074204, -0.470998929, -0.470998929, -0.470998929, -0.463481685]`

Corrected full-parent atomic contract:

- `count = 19`
- `fixed_dim = 931`
- `nshells = 6`
- `direct_core_len = 7`
- `support_count = 6859`
- `expected_support_count = 6859`
- `missing_row_count = 0`
- `working_box = (1:19, 1:19, 1:19)`
- `full_parent_working_box = true`
- `ownership_group_count_min = 1`
- `ownership_group_count_max = 1`
- `low_fixed_h1 = [-1.998226506, -0.499382797, -0.490167263, -0.490167263, -0.490167263]`

### `ns9_family`

Old contaminated central-box fixture:

- `count = 25`
- `fixed_dim = 517`
- `nshells = 4`
- `direct_core_len = 5`
- `support_count = 2197`
- `expected_support_count = 15625`
- `missing_row_count = 13428`
- `working_box = (7:19, 7:19, 7:19)`
- `full_parent_working_box = false`
- `low_fixed_h1 = [-1.982604385, 0.050419602, 0.050419602, 0.050419602, 0.525979157]`

Corrected full-parent atomic contract:

- `count = 25`
- `fixed_dim = 1513`
- `nshells = 8`
- `direct_core_len = 9`
- `support_count = 15625`
- `expected_support_count = 15625`
- `missing_row_count = 0`
- `working_box = (1:25, 1:25, 1:25)`
- `full_parent_working_box = true`
- `ownership_group_count_min = 1`
- `ownership_group_count_max = 1`
- `low_fixed_h1 = [-1.997931605, -0.499137420, -0.485603493, -0.485603493, -0.485603493]`

## Diagnosis Boundary After Rebasing

The earlier bug hunt was aimed at the wrong object.

The contaminated central-box fixture made `ns9` look catastrophically bad:

- the likely triplet was already positive
- most parent rows were not even covered

That failure does not survive the contract correction.

On the corrected full-parent atomic contract:

- `ns7` stays sane
- `ns9` also stays qualitatively sane
- both anchor cases have full parent coverage
- both anchor cases pass the one-piece row-ownership audit on the retained
  support
- the low fixed `H1` triplet is negative again for `ns9`

So the earlier “shell-2-and-inward balance problem” diagnosis is superseded,
not merely weakened. It was downstream of the wrong coverage fixture.

## Conclusion

The durable conclusion for the one-center atomic anchor pair is:

- the corrected atomic contract is full parent coverage
- the canonical repo-native path is now
  `build_one_center_atomic_full_parent_shell_sequence(...)` /
  `one_center_atomic_full_parent_fixed_block(...)`
- the old central-box fixture should not be used for one-center atomic
  localization diagnosis
- on the corrected full-parent backbone, `ns9_family, Z = 2, rmax = 10` is no
  longer exhibiting the earlier qualitative fixed-block failure

If contraction-subspace debugging is revisited later, it should start from the
full-parent atomic shell-sequence contract only.
