"""Rebuild the whole paper-style figure set: python make_all.py [--clean]
--clean  deletes the derived caches first (cache/events_*, roi_*, excerpt_*, whisk_*, loom_roi_*, ctx_*, ledger.json); the Z:-drive inventory cache (cache/inventory.pkl) is kept unless REFRESH_INV=1.
Order: analysis modules (event-triggered, ROI heterogeneity) -> figures 1-7 and 9-10 -> context_state.py (analysis of Fig. 11) -> fig11 (each writes its entries to the results ledger) -> ledger figure (2 pages) + results_ledger.csv.
Needs: C:\\Users\\rbondarenko\\2P_cache\\xs\\*.mat (export_xs.m), Z:\\ (suite2p folders for Fig. 1 / 4, rejected ROIs for Fig. 2 excerpts), the cross_session pickles in plots/cross_session.
Interpreter used for development: C:\\Users\\rbondarenko\\miniforge3\\envs\\suite2p_venv\\python.exe
"""
import os, sys, glob, subprocess, time

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, '..', 'plots', 'summary_figures', 'cache')
STEPS = ['event_triggered.py', 'roi_heterogeneity.py', 'fig1_overview.py', 'fig2_state_coupling.py', 'fig3_events.py', 'fig4_heterogeneity.py', 'fig5_stimuli.py', 'fig6_context.py', 'fig7_negative.py',
         'fig9_whisking.py', 'fig10_loom_heatmaps.py', 'context_state.py', 'fig11_task_engagement.py', 'rejected_rois.py', 'fig12_rejected_rois.py', 'fig8_ledger.py']          # fig8 (the ledger, 2 pages) must stay last

if __name__ == '__main__':
    if '--clean' in sys.argv:
        for pat in ('events_*.pkl', 'roi_*.pkl', 'excerpt_*.pkl', 'whisk_ex_*.pkl', 'whisk_partial_*.pkl', 'loom_roi_*.pkl', 'ctx_*.pkl', 'rejroi_*.pkl', 'ledger.json'):
            for f in glob.glob(os.path.join(CACHE, pat)): os.remove(f)
        print('caches cleaned')
    for s in STEPS:
        t0 = time.time(); r = subprocess.run([sys.executable, os.path.join(HERE, s)], cwd=HERE, capture_output=True, text=True)
        tail = [l for l in r.stdout.splitlines() if l.startswith(('saved', 'ledger rows'))]
        print('%-26s %s  %.0f s  %s' % (s, 'ok' if r.returncode == 0 else 'FAILED', time.time() - t0, '; '.join(tail)[:80]))
        if r.returncode != 0: print(r.stderr[-1500:]); sys.exit(1)
