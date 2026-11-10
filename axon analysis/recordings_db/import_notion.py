"""Merge an exported Notion table into notes.csv (the `notion_*` columns). Never touches your other columns, never overwrites a non-empty notion_* cell
unless --overwrite is given. Dry run by default.

    miniforge3\\envs\\suite2p_venv\\python.exe import_notion.py "C:\\path\\to\\notion_export.csv"            # shows what would be written
    miniforge3\\envs\\suite2p_venv\\python.exe import_notion.py "C:\\path\\to\\notion_export.csv" --apply

Notion: open the table -> "..." menu -> Export -> format CSV (or Markdown & CSV, then take the .csv). The row is matched to a session through a recording id
anywhere in the row (e.g. "RB24_261021_1630" or "AnimalRB24_261021_1630" in the Name / title / date columns). Columns are matched by name (case, spaces and
punctuation ignored): planes, quality, visible activity, comments, suite2p, stimuli / camera. Anything else is appended to notion_comments as "column: value".
Unmatched rows and unmatched sessions are listed so nothing is silently lost.
"""
import os, re, sys, csv, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_recdb as B

ID = re.compile(r'(?:Animal)?(RB2[234])[_\- ]?(\d{6})[_\- ]?(\d{4})')
TARGETS = {'planes': 'notion_planes', 'quality': 'notion_quality', 'visibleactivity': 'notion_visible_activity', 'activity': 'notion_visible_activity', 'comments': 'notion_comments',
           'comment': 'notion_comments', 'notes': 'notion_comments', 'suite2p': 'notion_suite2p', 'stimulicamera': 'notion_stimuli_camera', 'stimuli': 'notion_stimuli_camera',
           'stimulicameras': 'notion_stimuli_camera'}


def norm(h): return re.sub(r'[^a-z0-9]', '', h.lower())


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('csv'); ap.add_argument('--apply', action='store_true'); ap.add_argument('--overwrite', action='store_true'); a = ap.parse_args()
    with open(a.csv, newline='', encoding='utf-8-sig') as f: rows = list(csv.DictReader(f))
    notes_p = os.path.join(B.DB, 'notes.csv'); notes = B.read_csv(notes_p)
    if not notes: sys.exit('notes.csv does not exist yet: run build_recdb first')
    by = {r['session']: r for r in notes}
    matched, unmatched, extra_cols = {}, [], set()
    for r in rows:
        sid = None
        for v in r.values():
            m = ID.search(v or '')
            if m: sid = 'Animal%s_%s_%s' % m.groups(); break
        if sid is None or sid not in by: unmatched.append((sid, next(iter(r.values()), ''))); continue
        upd = {}
        for h, v in r.items():
            if v is None or not v.strip(): continue
            t = TARGETS.get(norm(h))
            if t: upd[t] = (upd.get(t, '') + ' | ' + v.strip()).strip(' |')
            elif not ID.search(v) and norm(h) not in ('name', 'title', 'session', 'date', 'recording'):
                extra_cols.add(h); upd['notion_comments'] = (upd.get('notion_comments', '') + ' | %s: %s' % (h, v.strip())).strip(' |')
        matched[sid] = upd
    n_w = 0
    for sid, upd in matched.items():
        for k, v in upd.items():
            if by[sid].get(k) and not a.overwrite: continue
            if by[sid].get(k) != v: by[sid][k] = v; n_w += 1
    print('%d Notion rows, %d matched to sessions, %d unmatched, %d cells %s' % (len(rows), len(matched), len(unmatched), n_w, 'written' if a.apply else 'would be written (dry run)'))
    if extra_cols: print('columns without a notion_* target (appended to notion_comments):', sorted(extra_cols))
    for sid, first in unmatched: print('  unmatched Notion row: id=%s  first cell=%r' % (sid, first[:60]))
    miss = [s for s in by if s not in matched]
    print('sessions with no Notion row: %d %s' % (len(miss), ', '.join(m.replace('Animal', '') for m in miss[:50])))
    if a.apply:
        cols = list(notes[0].keys()); B.write_csv(notes_p, [by[r['session']] for r in notes], cols); print('notes.csv updated; re-run build_recdb to refresh sessions.csv / index.html')


if __name__ == '__main__':
    main()
