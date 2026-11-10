"""recdb - read and index the RB22-24 recordings database (stdlib only, works in every env; numpy not required).

    import sys; sys.path.insert(0, r'C:\\Users\\rbondarenko\\projects\\2P_Scripts\\axon analysis\\recordings_db')
    import recdb
    sel = recdb.select(animal='RB24', has_block='loom_expand', min_rois=10, pool=True)       # list of session rows (dicts), sorted by animal / date
    names = recdb.names(sel)                                                                    # ['AnimalRB24_261021_1630', ...]
    recdb.show(sel, ['session', 'stimuli', 'n_roi_accepted', 'flag_reason'])                    # quick table in the console
    rois = recdb.rois(names[0], accepted=True, cls='signal')                                    # ROI rows of one session
    blocks = recdb.blocks(names[0], key='loom_expand')                                          # stimulus blocks of one session

select() filters (all optional, combined with AND):
    animal='RB24' | ['RB23','RB24']      status='analysed' | [...]      task_type='lick task'      pool=True (pool == 'yes') | False (pool == 'no')
    tags=['loom','task']  all of these in all_tags      any_tags=[...]  at least one      exclude_tags=['stale','interrupted']
    has_block='loom_expand' (key of blocks.csv; with min_trials=10: at least that many completed trials in total)
    min_rois=30 (accepted ROIs)      min_minutes=10      date_from='2026-10-27'  date_to='2026-11-06'      aggregated=True      where=lambda row: ...
The database folder is Z:\\Group Members\\Rima\\Recordings_DB (override with the environment variable RECDB_DIR; falls back to the git snapshot when Z: is not mounted).
"""
import os, csv

_HERE = os.path.dirname(os.path.abspath(__file__))
_BOOL = {'has_stim', 'has_tif', 'has_suite2p', 'has_fall', 'has_camera', 'has_dlc', 'has_aggregated', 'in_session_qc', 'stale', 'name_date_mismatch', 'accepted', 'multi_plane'}
_cache = {}


def db_dir():
    cands = [os.environ.get('RECDB_DIR', ''), r'Z:\joeschgrp\Group Members\Rima\Recordings_DB', r'Z:\Group Members\Rima\Recordings_DB', os.path.join(_HERE, 'snapshot')]
    for d in cands:
        if d and os.path.isfile(os.path.join(d, 'sessions.csv')): return d
    raise FileNotFoundError('no recordings database found (tried %s)' % [c for c in cands if c])


def _isnum(v):
    try: float(v); return v.lower() not in ('nan', 'inf', '-inf', 'infinity')
    except ValueError: return False


def table(name, reload=False):
    """rows of <name>.csv (sessions | blocks | rois | axons | animals | notes) as a list of dicts.
    Per column: yes/no columns -> bool; columns whose every non-blank value is a number -> int / float (blank = None); everything else stays text ('' when blank)."""
    if reload or name not in _cache:
        with open(os.path.join(db_dir(), name + '.csv'), newline='', encoding='utf-8') as f: raw = list(csv.DictReader(f))
        cols = list(raw[0].keys()) if raw else []
        kind = {}
        for c in cols:
            vals = [r[c] for r in raw if r[c] != '']
            if c in _BOOL: kind[c] = 'bool'
            elif vals and all(_isnum(v) for v in vals): kind[c] = 'int' if all(v.lstrip('-').isdigit() for v in vals) else 'float'
            else: kind[c] = 'text'
        def cv(c, v):
            k = kind[c]
            if k == 'text': return v
            if v == '': return None
            if k == 'bool': return {'yes': True, 'no': False}.get(v)
            return int(v) if k == 'int' else float(v)
        _cache[name] = [{c: cv(c, r[c]) for c in cols} for r in raw]
    return _cache[name]


def _as_list(x): return [] if x is None else ([x] if isinstance(x, str) else list(x))


def select(animal=None, status=None, task_type=None, pool=None, tags=None, any_tags=None, exclude_tags=None, has_block=None, min_trials=None, min_rois=None, min_minutes=None,
           date_from=None, date_to=None, aggregated=None, where=None):
    rows = table('sessions'); blk = table('blocks') if (has_block or min_trials) else []
    an, st = ['Animal' + a if not a.startswith('Animal') else a for a in _as_list(animal)], _as_list(status)
    out = []
    for r in rows:
        tg = set(t for t in (r['all_tags'] or '').split(';') if t)
        if an and r['animal'] not in [a.replace('Animal', '') for a in an]: continue
        if st and r['status'] not in st: continue
        if task_type and r['task_type'] != task_type: continue
        if pool is True and r['pool'] != 'yes': continue
        if pool is False and r['pool'] != 'no': continue
        if tags and not set(_as_list(tags)) <= tg: continue
        if any_tags and not (set(_as_list(any_tags)) & tg): continue
        if exclude_tags and (set(_as_list(exclude_tags)) & tg): continue
        if has_block:
            bs = [b for b in blk if b['session'] == r['session'] and b['key'] == has_block]
            if not bs: continue
            if min_trials and sum((b['n_trials'] or 0) for b in bs) < min_trials: continue
        if min_rois is not None and (r['n_roi_accepted'] or 0) < min_rois: continue
        if min_minutes is not None and (r['minutes'] or 0) < min_minutes: continue
        if date_from and r['date'] < date_from: continue
        if date_to and r['date'] > date_to: continue
        if aggregated is not None and bool(r['has_aggregated']) != aggregated: continue
        if where and not where(r): continue
        out.append(r)
    return sorted(out, key=lambda r: (r['animal'], r['date'], r['time']))


def names(sel): return [r['session'] for r in sel]


def rois(session, accepted=None, cls=None):
    return [r for r in table('rois') if r['session'] == session and (accepted is None or r['accepted'] == accepted) and (cls is None or r['class'] == cls)]


def blocks(session, key=None):
    return [b for b in table('blocks') if b['session'] == session and (key is None or b['key'] == key)]


def show(sel, cols=('session', 'status', 'task_type', 'stimuli', 'minutes', 'n_roi_accepted', 'pool')):
    w = [max(len(c), *(len(str(r[c]) if r[c] is not None else '') for r in sel)) if sel else len(c) for c in cols]
    print('  '.join(c.ljust(n) for c, n in zip(cols, w)))
    for r in sel: print('  '.join((str(r[c]) if r[c] is not None else '').ljust(n) for c, n in zip(cols, w)))
    print('%d session(s)' % len(sel))


if __name__ == '__main__':
    show(select())
