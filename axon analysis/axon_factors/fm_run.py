"""python fm_run.py [--domain band|dff] [session ...]  -> <2P_cache>/axon_factors/<session>_fm_<domain>.pkl and the printed summary"""
import os, sys, pickle
import fm_fit as fm
DEFAULT = ['AnimalRB24_261021_1630', 'AnimalRB22_261027_1138', 'AnimalRB23_261023_1642', 'AnimalRB24_261105_1423', 'AnimalRB24_261105_1520']
if __name__ == '__main__':
    args = sys.argv[1:]; domain = 'band'
    if '--domain' in args: domain = args[args.index('--domain') + 1]; args = [a for a in args if a not in ('--domain', domain)]
    os.makedirs(fm.CACHE, exist_ok=True)
    for nm in (args or DEFAULT):
        o = fm.analyse(nm, domain=domain)
        with open(os.path.join(fm.CACHE, '%s_fm_%s.pkl' % (nm, domain)), 'wb') as fh: pickle.dump(o, fh)
