"""Figure 8 - results ledger: every tested hypothesis with its effect size, uncertainty, n and verdict (collected from the figure scripts: cache/ledger.json).
usage: python fig8_ledger.py -> plots/summary_figures/Fig8_ledger.pdf/.png + results_ledger.csv
Verdict (never colour only): SUPPORTED = dark chip, AMBIGUOUS = grey chip, NOT SEEN = open chip, NOT AVAILABLE = hatched.
Verdict rule used by the figure scripts: 'supported' needs the 95% CI (animal-clustered bootstrap) away from 0 and agreement across sessions; any result carried by a single animal is downgraded to 'ambiguous' automatically.
"""
import csv, json, textwrap
import numpy as np
from fig_common import *
from matplotlib.patches import FancyBboxPatch, Rectangle

SECTIONS = [('Fig2', 'Fig. 2  Coupling of the axon mean to behaviour'), ('Fig3', 'Fig. 3  Arousal events and state tuning'), ('Fig4', 'Fig. 4  ROI structure'), ('Fig5', 'Fig. 5  Visual stimuli'), ('Fig6', 'Fig. 6  Lick task and grey screens'),
            ('Fig7', 'Fig. 7  Controls and limits of the evidence'), ('Data', 'Data limits (no statistics)'),
            ('Fig9', 'Fig. 9  Whisking vs axon calcium (shuffle tests)'), ('Fig10', 'Fig. 10  Looming discs: heatmaps and per-session responses'), ('Fig11', 'Fig. 11  Lick task vs passive square, engagement (behaviours regressed)'), ('Fig12', 'Fig. 12  What the rejected ROIs are (control population)')]
PAGES = [('Fig8_ledger', ['Fig2', 'Fig3', 'Fig4', 'Fig5', 'Fig6', 'Fig7', 'Data'], 'Results ledger: what was tested and what came out (exploratory, 3 mice, n = unflagged sessions unless stated)'),
         ('Fig8b_ledger_figs9-12', ['Fig9', 'Fig10', 'Fig11', 'Fig12'], 'Results ledger, continued: Figs 9-12 (whisking shuffle tests, loom heatmaps, lick task vs passive, rejected ROIs)')]
CHIP = {'supported': ('SUPPORTED', INK, 'white', INK), 'ambiguous': ('AMBIGUOUS', '#B8BDC6', INK, '#B8BDC6'), 'not seen': ('NOT SEEN', 'white', MUTED, MUTED), 'not available': ('NOT AVAILABLE', '#F3F4F6', MUTED, FAINT)}


def static_rows():
    n_ball = sum(1 for n in ALL if QC[n]['ball'].startswith('used')); n_ball_u = sum(1 for n in UNFLAGGED if QC[n]['ball'].startswith('used'))
    return [dict(figure='Data', key='n_mice', label='animals with analysable axon sessions', est=None, lo=None, hi=None, n_sessions=len(ALL), n_animals=3, verdict='not available', note='RB24 %d, RB23 %d, RB22 %d sessions: 3 mice, session is not an independent sample' % tuple(sum(1 for n in ALL if animal(n) == a) for a in ('RB24', 'RB23', 'RB22')), unit=''),
            dict(figure='Data', key='ball', label='locomotion (ball) usable', est=None, lo=None, hi=None, n_sessions=n_ball, n_animals=1, verdict='not available', note='%d of %d sessions (%d unflagged): open-loop logs truncated, implausible readings on 11-06, 2nd ball sensor dead' % (n_ball, len(ALL), n_ball_u), unit=''),
            dict(figure='Data', key='ref', label='reference channel (tdTomato) / axial-motion control', est=None, lo=None, hi=None, n_sessions=0, n_animals=0, verdict='not available', note='only GCaMP imaged; 5-plane z-registration (Suite3D) failed at 3 um plane spacing', unit=''),
            dict(figure='Data', key='unaggregated', label='recorded sessions not aggregated (RB22-24)', est=None, lo=None, hi=None, n_sessions=19, n_animals=3, verdict='not available', note='see Fig. 1g: several have suite2p ROIs (not curated here); aggregating them would add RB22 / RB23 sessions', unit='')]


UNIT_SHORT = {'net sign (+1 = all positive)': 'net sign', 'dF/F (high-passed)': 'dF/F', '% of subsession 1': '%', 'SD of metric': 'SD', 'fraction': 'frac.', 'fraction of sessions': 'of sessions', 'fraction of ROIs': 'of ROIs', 'dF/F per hit-vs-miss': 'dF/F', 'dF/F per SD of hit rate': 'dF/F per SD'}


def fmt_est(r):
    if r['est'] is None: return 'n/a'
    big = abs(r['est']) >= 10; small = 0 < abs(r['est']) < 0.1
    s = ('%+.3f' if small else '%+.2f' if not big else '%.0f') % r['est']
    if r['lo'] is not None and r['hi'] is not None: s += '  [' + (('%.3f, %.3f' if small else '%.2f, %.2f' if not big else '%.0f, %.0f') % (r['lo'], r['hi'])) + ']'
    return s


def main():
    L = json.load(open(LEDGER_PATH)); L = [r for r in L if r['figure'] in dict(SECTIONS)] + static_rows()
    with open(os.path.join(OUT, 'results_ledger.csv'), 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f); w.writerow(['figure', 'key', 'what was tested', 'estimate', 'ci_low', 'ci_high', 'unit', 'n_sessions', 'n_animals', 'verdict', 'note'])
        for r in sorted(L, key=lambda r: [s[0] for s in SECTIONS].index(r['figure'])): w.writerow([r['figure'], r['key'], r['label'], r['est'], r['lo'], r['hi'], r['unit'], r['n_sessions'], r['n_animals'], r['verdict'], r['note']])
    for name, secs, ptitle in PAGES: page(L, secs, name, ptitle)


def page(L, secs, name, ptitle):
    rows = []
    for sec, title in SECTIONS:
        if sec not in secs: continue
        rs = [r for r in L if r['figure'] == sec]
        if not rs: continue
        rows.append(('H', title))
        for r in rs: rows.append(('R', r))
    H = 4.2 * len(rows) + 22; fig = plt.figure(figsize=(W2, max(H, 60) * MM)); ax = fig.add_axes([0.01, 0.035, 0.98, 0.905]); ax.set_xlim(0, 100); ax.set_ylim(len(rows), -1.6); ax.axis('off')
    fig.text(0.01, 0.985, ptitle, fontsize=8, fontweight='bold', va='top')
    X = dict(label=1.0, est=61.0, n=80.0, chip=88.0)
    for h, x in (('what was tested', X['label']), ('effect  [95% CI]  unit', X['est']), ('sessions / mice', X['n']), ('verdict', X['chip'])): ax.text(x, -1.0, h, fontsize=5.8, fontweight='bold', va='center', color=MUTED)
    fig.text(0.01, 0.006, textwrap.fill('Verdict rule: SUPPORTED = 95% CI (hierarchical bootstrap over animals / sessions) excludes 0 and most sessions agree; any result carried by one animal is downgraded to AMBIGUOUS automatically; NOT SEEN = CI includes 0 and no session-level support. '
             'Effects are in different units and are not comparable with each other. "Unit" columns: r = correlation, z = session z, dF/F = high-passed, SD = SD of the metric over the session\'s windows. n = 3 mice at most: every number describes these animals.', 200), fontsize=4.8, color=MUTED, va='bottom')
    ax.plot([0, 100], [-0.45, -0.45], color=MUTED, lw=0.5)
    for i, (kind, r) in enumerate(rows):
        if kind == 'H':
            ax.add_patch(Rectangle((0, i - 0.45), 100, 0.9, color='#EEF0F3', lw=0)); ax.text(0.6, i, r, fontsize=6.4, fontweight='bold', va='center'); continue
        lab = r['label']; ax.text(X['label'], i, lab if len(lab) < 112 else lab[:109] + '…', fontsize=5.4, va='center')
        ax.text(X['est'], i, fmt_est(r) + ('  ' + UNIT_SHORT.get(r['unit'], r['unit']) if r['unit'] and r['est'] is not None else ''), fontsize=5.4, va='center', color=INK)
        ax.text(X['n'], i, '%d / %d' % (r['n_sessions'], r['n_animals']) if r['n_sessions'] else '-', fontsize=5.4, va='center', color=MUTED)
        txt, fc, tc, ec = CHIP[r['verdict']]; ax.add_patch(FancyBboxPatch((X['chip'], i - 0.36), 11.5, 0.72, boxstyle='round,pad=0.0,rounding_size=0.3', fc=fc, ec=ec, lw=0.7, hatch='////' if r['verdict'] == 'not available' else None)); ax.text(X['chip'] + 5.75, i, txt, fontsize=4.6, fontweight='bold', color=tc, ha='center', va='center')
    savefig(fig, name)
    print('ledger rows', sum(1 for k, _ in rows if k == 'R'))


if __name__ == '__main__':
    main()
