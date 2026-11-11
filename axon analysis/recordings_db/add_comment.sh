#!/usr/bin/env bash
# Add a comment (or set any notes.csv column) for a session, then refresh the database tables.
# Usage: ./add_comment.sh AnimalRB24_261109_1506 "dead 2nd ball sensor, ignore flow asymmetry"
#        ./add_comment.sh AnimalRB24_261109_1506 "no" pool            # set another column (pool, flag_reason, curation, tags, quality, visible_activity)
# Comments are appended (dated), other columns are overwritten. Needs: Z: mounted. Edit notes.csv in Excel instead if you prefer (close it before building).
set -euo pipefail
[ $# -ge 2 ] || { echo "usage: $0 SESSION \"text\" [column=comments]"; exit 1; }
PY="C:/Users/rbondarenko/miniforge3/envs/suite2p_venv/python.exe"
NOTES="${RECDB_NOTES:-Z:/Group Members/Rima/Recordings_DB/notes.csv}"
"$PY" - "$NOTES" "$1" "$2" "${3:-comments}" <<'PY'
import csv, sys, datetime, shutil
path, sess, text, col = sys.argv[1:5]
with open(path, newline='', encoding='utf-8-sig') as f:
    r = csv.DictReader(f); rows = list(r); head = r.fieldnames
if col not in head: sys.exit(f"no column '{col}'; columns: {head}")
row = next((x for x in rows if x['session'] == sess), None)
if row is None: sys.exit(f"session '{sess}' not in notes.csv")
shutil.copy(path, path + '.bak')
if col == 'comments':
    new = f"[{datetime.date.today()}] {text}"
    row[col] = (row[col] + ' | ' + new) if row[col] else new
else:
    row[col] = text
with open(path, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=head); w.writeheader(); w.writerows(rows)
print(f"{sess}.{col} = {row[col]}")
PY
if [ "${RECDB_NOBUILD:-0}" != "1" ]; then
  matlab -batch "cd('C:/Users/rbondarenko/projects/2P_Scripts/axon analysis/recordings_db'); build_recdb"
fi
