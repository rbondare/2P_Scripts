import h5py, numpy as np, os, time
here = os.path.dirname(os.path.abspath(__file__))
p = r'Z:\Group Members\Rima\Stimulus_Axons\AnimalRB24\AnimalRB24_261109_1506\stim_001__training_20261109_150656.mat'
f = h5py.File(p, 'r')
b = f['BallCamData']
n = b.shape[0]
out = np.zeros((n, 5), np.float64)
t0 = time.time()
for i in range(n):
    out[i] = np.array(f[b[i, 0]]).ravel()[:5]
    if i % 100000 == 0:
        print(i, time.time() - t0, flush=True)
np.save(os.path.join(here, 'ball_AnimalRB24_261109_1506.npy'), out)
print('done', out.shape, 'nonzero fraction per column', (out != 0).mean(0).round(3), flush=True)
