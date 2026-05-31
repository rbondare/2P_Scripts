# Suite2P Axon Detection Guide

## Quick Start

```bash
conda activate suite2p_venv
cd C:\Users\rbondarenko\projects\2P_Scripts\suite2p_analysis
python run_suite2p_RB_axons_V1.py
```

---

## Key Differences: Axons vs Soma

| Parameter | Soma Detection | Axon Detection | Reason |
|-----------|---|---|---|
| **Diameter** | 15-20 μm | 3-5 μm | Axons are much thinner |
| **Cellprob threshold** | 0.5 | 0.1 | Axons are fainter |
| **Threshold scaling** | 1.0 | 0.5 | Need lower sensitivity threshold |
| **Spatial HP detect** | 25 | 10 | Preserve small structures |
| **Block size** | 128×128 | 64×64 | Finer spatial detail |
| **Inner neuropil radius** | 2-3 | 1 | Axons have minimal neuropil |
| **Min neuropil pixels** | 350 | 50 | Much smaller surround area |
| **Max pixel norm** | 200+ | 100 | Axons are smaller |
| **Spatial scale** | 0 (default) | 2 | Better for fine structures |

---

## Parameter Tuning Guide

### If you detect TOO MANY axons (including noise):
**Increase these:**
```python
'cellprob_threshold': 0.3  # instead of 0.1
'threshold_scaling': 0.8   # instead of 0.5
'flow_threshold': 0.6      # instead of 0.4
```

### If you're missing faint axons:
**Decrease these:**
```python
'cellprob_threshold': 0.05  # even lower
'threshold_scaling': 0.3    # even lower
'spatial_hp_detect': 5.0    # preserve more detail
'diameter': [2.0, 2.0]      # even smaller diameter
```

### If ROI masks are too large (including surrounding noise):
**Adjust these:**
```python
'npix_norm_max': 50.0       # smaller cap
'inner_neuropil_radius': 0  # no neuropil
'allow_overlap': True       # allow overlapping axons
```

### If processing is too slow:
**Change these:**
```python
'block_size': [32.0, 32.0]  # smaller blocks
'batch_size': 50            # smaller batches
'max_iterations': 10        # fewer iterations
```

---

## Expected Output

After running, you'll find:
```
suite2p/
├── plane0/
│   ├── stat.npy           # Detected axon locations and masks
│   ├── iscell.npy         # Classification (cell vs non-cell)
│   ├── F.npy              # Fluorescence traces
│   ├── Fneu.npy           # Neuropil fluorescence
│   └── spks.npy           # Inferred spikes
├── plane1/
├── plane2/
└── plane3/
```

---

## Validation Checklist

After detection:

1. **Open in Suite2P GUI** to visualize detected axons
2. **Check for false positives** (neuropil, blood vessels)
3. **Check for missed axons** (very faint structures)
4. **Inspect iscell classifications** (should separate axons from garbage)
5. **Check ROI masks** (should be thin, linear for axons)

---

## Example Adjustments for Your Data

Based on the faint, sparse axon image you showed:

**Conservative (fewer false positives):**
```python
'cellprob_threshold': 0.2
'threshold_scaling': 0.6
'diameter': [2.5, 2.5]
'spatial_hp_detect': 8.0
```

**Aggressive (more sensitivity):**
```python
'cellprob_threshold': 0.05
'threshold_scaling': 0.3
'diameter': [2.0, 2.0]
'spatial_hp_detect': 5.0
```

---

## Troubleshooting

### Error: "overlap key missing"
- Delete `stat.npy` and `iscell.npy` from plane directories
- Re-run script

### Very slow processing
- Reduce `block_size` to `[64, 64]` or smaller
- Reduce `batch_size` to 50 or 25
- Disable `nonrigid`: set to `False`

### Memory errors
- Set `do_regmetrics: False` (already done in script)
- Reduce `nimg_init` to 200 or 100
- Use smaller `batch_size`

### Too many false detections
- Increase `cellprob_threshold` (currently 0.1)
- Increase `threshold_scaling` (currently 0.5)
- Increase `spatial_hp_detect` (currently 10.0)

---

## Comparing to Your Soma Detection

Run both:
1. **Soma**: `run_suite2p_RB_soma_V1.py` → somas only
2. **Axons**: `run_suite2p_RB_axons_V1.py` → axons only

Then compare the outputs to see differences in:
- Number of detected ROIs
- ROI sizes and shapes
- Signal quality

---

## References

- **Suite2P Detection**: https://suite2p.readthedocs.io/en/latest/
- **Cellpose models**: https://www.cellpose.org/
- **For axon-specific analysis**: Consider post-processing with morphological filtering
