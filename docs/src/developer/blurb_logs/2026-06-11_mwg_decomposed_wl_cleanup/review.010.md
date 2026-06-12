Review result:

Accepted as a correctly blocked readiness pass. The doer followed the blurb
trust boundary and did not substitute a different Be basis source when
`/Users/srw/BasisSets` was missing.

Main result:

- status: `blocked`
- blocker: `:missing_required_be_ccpv5z_basisfile`
- required basisfile: `/Users/srw/BasisSets`
- `isfile("/Users/srw/BasisSets") == false`
- no old nested/QW Be S+P oracle was built
- no decomposed/final-basis Be S+P audit was attempted

Manager follow-up:

I checked the local basis sources and found Be `cc-pV5Z` in the existing
GaussletModules basis files:

- `/Users/srw/Dropbox/GaussletModules/BasisSets`
- `/Users/srw/Library/CloudStorage/Dropbox/GaussletModules/BasisSets`
- `/Users/srw/Library/CloudStorage/Dropbox/oldgaussletmodules/BasisSets`

The active follow-up should authorize the CloudStorage GaussletModules file
explicitly instead of asking doer to rediscover or guess a replacement.

Validation reviewed:

- `julia --project=. tmp/work/be_atom_sp_qw_oracle_probe.jl`
- `julia --project=. -e 'using GaussletBases; println("load ok")'`
- `git diff --check`
- artifact: `tmp/work/be_atom_sp_qw_oracle_summary.txt`

Deletion/shrinkage review:

No old surface became obsolete. No test was added. The pass stopped before any
route construction because the required machine-local basisfile was missing.

Commit/push:

Pending manager commit/push of this tracked response/review log and the next
published blurb.

Next target:

Rerun the Be S+P old nested/QW oracle using the explicitly authorized
GaussletModules `BasisSets` file, then continue to the decomposed/final-basis
audit only after the old oracle is clean.
