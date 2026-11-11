#!/usr/bin/env bash
# Re-aggregate one or more recordings after re-curating them in suite2p, then refresh the database.
# Usage: ./reaggregate.sh AnimalRB24_261109_1506 [AnimalRB24_261106_1317 ...]
# Before: in suite2p GUI curate, then File > "Save to Matlab" (re-exports Fall.mat; aggregation reads Fall.mat, NOT iscell.npy).
set -euo pipefail
[ $# -ge 1 ] || { echo "usage: $0 SESSION [SESSION ...]"; exit 1; }
REPO="C:/Users/rbondarenko/projects/2P_Scripts"
INTER="C:/Users/rbondarenko/projects/IntermediateDir"
XS="C:/Users/rbondarenko/2P_cache/xs"
LIST=""
for s in "$@"; do
  animal="${s%%_*}"                        # AnimalRB24_261109_1506 -> AnimalRB24
  rm -fv "$INTER/$animal/${s}_S2Presult.mat"      # cached copy of the OLD Fall.mat: must go, else the curation is ignored
  rm -fv "$XS/${s}_xs.mat"                # cross-session cache (FOV colours); rebuilt by the database build
  LIST="$LIST'$s',"
done
LIST="{${LIST%,}}"
# 1) aggregate (deletes the old *_preprocessed.mat on Z: first; 30-75 min per session)
matlab -batch "addpath('$REPO/newAggregation'); run_aggregation_batch($LIST, 'overwrite', true)"
# 2) refresh the database (~1-2 min)
matlab -batch "cd('$REPO/axon analysis/recordings_db'); build_recdb"
echo "done. Check sessions.csv: stale must be blank for: $*"
