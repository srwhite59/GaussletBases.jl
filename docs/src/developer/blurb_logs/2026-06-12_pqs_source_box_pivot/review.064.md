Review 064: accepted with one manager-side correction.

The implementation added the intended direct retained PQS position/x2 moment
blocks and validated them against the raw-source-block plus retained-selector
oracle. The shared helper in `pqs_source_safe_terms.jl` is the right shape:
choose the requested axis, use the operator factor on that axis and overlaps on
the others, then fill directly over retained mode tuples.

I made one correction before committing: the response said the six new wrappers
were exported through `CartesianPairBlockMaterialization`, but the blurb
explicitly said not to add exports. The tests call them as qualified CPBM
functions, so export entries were unnecessary. I removed the six new export
entries and left the functions available as module-qualified internal helpers.

Validation status:

- doer ran the focused CPBM contract file, 1140 tests, elapsed about 71s;
- doer ran the load check and `git diff --check`;
- after the manager-side export removal, I reran `julia --project=. -e 'using
  GaussletBases; println("load ok")'` and `git diff --check`.

Remaining small cleanup: the named direct retained moment wrappers exist, but
the generic retained one-body selector/matrix path still appears to expose only
overlap and kinetic. It may now be reasonable to wire position/x2 into that
existing selector path, as long as doing so does not widen public exports or
turn into driver/physics work.

-- repo-manager@macmini
