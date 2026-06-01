# Pupil-Calcium Activity Distribution Visualization Guide

## 📊 Quick Comparison of Visualization Methods

| Visualization | Best For | Pros | Cons | Interpreting |
|---|---|---|---|---|
| **Scatter + Marginals** | Overall relationship & outliers | See all data points; shows density hotspots | Can be dense with large datasets | Look for linear/nonlinear trends and density clusters |
| **Violin Plot** | Distribution shapes across bins | Beautiful; shows full distribution; easy to compare | Harder to see individual points; requires binning | Wider = more values at that calcium level; symmetry indicates distribution type |
| **2D Density Heatmap** | Joint distribution & correlation | Clear hotspots; compact; good for papers | Less detail on tails; requires interpolation | Red regions = most common combinations |
| **Ridge Plot** | Distribution changes aesthetically | Professional appearance; shows transitions clearly | Just visualizing; harder to extract stats | Left-skewed = more low activity at small pupils? |
| **Box/Violin by Bin** | Statistical comparison | Shows quartiles, means, outliers, sample sizes | Multiple sub-plots can be confusing | Increasing/decreasing trend tells you about state dependence |
| **Scatter with Time** | Temporal dynamics | See when high/low pupil occurs; behavioral context | Can be cluttered if many ROIs | Colored regions show pupil states over time |

---

## 🎯 Recommended Analysis Workflow

### Step 1: Start with **Scatter + 2D Density Heatmap**
- Get the overall picture
- Check for correlation (positive = calcium increases with pupil; negative = opposite)
- Identify nonlinearity

### Step 2: Use **Violin/Box Plots by Pupil Bins**
- Quantify changes across pupil states
- Compare means, medians, and variability
- Test for statistical differences between bins

### Step 3: Add **Ridge Plot** or **Time Series**
- For presentations/papers
- Show temporal patterns
- Demonstrate state-dependent modulation

---

## 📈 What Different Patterns Mean

### Pattern: Scatter Points in Upper-Right (High Pupil → High Calcium)
**Interpretation**: Arousal state correlates with increased neural activity
- Likely explanation: Pupil dilation = sympathetic activation = heightened alertness
- Behavioral context: Active exploration, attention

### Pattern: Scatter Points in Lower-Right (High Pupil → Low Calcium)
**Interpretation**: Arousal state correlates with decreased neural activity
- Likely explanation: Inhibitory neurons active during arousal
- Behavioral context: Might reflect inhibitory circuits controlling sensory processing

### Pattern: Vertical Scatter (Pupil → Variable Calcium)
**Interpretation**: Calcium activity independent of pupil
- Likely explanation: These neurons driven by other factors (stimulus, internal state)
- Behavioral context: Task-specific activity that doesn't correlate with arousal

### Pattern: Horizontal Scatter (Calcium → Relatively Fixed Pupil)
**Interpretation**: Pupil relatively stable despite calcium fluctuations
- Likely explanation: Activity in these neurons doesn't affect autonomic state
- Behavioral context: Local processing without systemic effect

### Pattern: V-shaped or U-shaped
**Interpretation**: Calcium activity extreme at both small and large pupils
- Likely explanation: Complex arousal-activity relationship
- Behavioral context: Activity at both rest and high arousal, but low at medium arousal

---

## 💡 Statistical Measures to Include

### Correlation Analysis
```
Pearson r: Linear relationship strength (-1 to 1)
  - r > 0.3: Weak positive correlation
  - r > 0.5: Moderate positive correlation
  - r > 0.7: Strong positive correlation

Spearman ρ: Rank correlation (less sensitive to outliers)
  - Use if relationship is nonlinear
```

### By-Bin Statistics
```
For each pupil size bin, report:
  - Mean ± Std (for normally distributed data)
  - Median ± IQR (for skewed data)
  - Sample size (n)
  - p-value (if comparing bins statistically)
```

### Effect Size
```
Cohen's d = (mean1 - mean2) / pooled_std
  - d < 0.2: Negligible
  - d = 0.2-0.5: Small effect
  - d = 0.5-0.8: Medium effect
  - d > 0.8: Large effect
```

---

## 🔍 Practical Analysis Tips

### 1. **Preprocessing**
- Remove outliers (e.g., > 3σ from mean)
- Handle NaNs carefully (don't just interpolate behavioral data)
- Consider temporal smoothing (but disclose window size)

### 2. **Binning Strategy**
- Too few bins (2-3): Loss of detail
- Too many bins (>10): Sparse data, noisy distributions
- **Recommended**: 5 bins for most analyses
- Consider **adaptive binning**: Equal-count bins instead of equal-width

### 3. **Multiple ROIs**
If you have multiple ROIs/axons:
- **Option A**: Average across ROIs, then visualize (simpler)
- **Option B**: Visualize top N responsive ROIs separately (more detail)
- **Option C**: Create distribution of all correlations (meta-analysis)

```matlab
% Correlation for each ROI
correlations = zeros(size(Ca_selected, 1), 1);
for roi = 1:size(Ca_selected, 1)
    correlations(roi) = corr(Beh_selected', Ca_selected(roi, :)');
end
histogram(correlations);
```

### 4. **Control for Confounds**
- Locomotion speed: Does pupil correlate with movement?
  ```
  corr_residual = corr(Beh_selected - predicted_from_locomotion, Ca_selected)
  ```
- Stimulus timing: Separate spontaneous from stimulus-driven activity

### 5. **Statistical Testing**
```matlab
% ANOVA across pupil bins
p_anova = anova1(calcium_by_bin, []);

% Regression: Does pupil predict calcium?
mdl = fitlm(Beh_selected, Ca_selected(1,:));
r2 = mdl.Rsquared.Ordinary;  % Variance explained

% Control for confounds: partial correlation
```

---

## 📊 Code Examples

### Python (Minimal Example)
```python
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Your data
pupil = data['Behav']['PupilArea']
calcium = data['CaData'][0]['Ca_dFF'].mean(axis=0)  # Average ROIs

# Quick visualization
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Scatter + density
axes[0].hexbin(pupil, calcium, gridsize=30, cmap='YlOrRd')
axes[0].set_xlabel('Pupil')
axes[0].set_ylabel('Calcium')

# By bins
pupil_bins = pd.cut(pupil, 5)
sns.violinplot(x=pupil_bins, y=calcium, ax=axes[1])

plt.tight_layout()
plt.show()
```

### MATLAB (Minimal Example)
```matlab
% Quick check
figure;
subplot(1,2,1);
scatter(Beh_selected, Ca_selected(1,:), 20, 'filled', 'o');
xlabel('Pupil'); ylabel('Calcium');

subplot(1,2,2);
boxplot(Ca_selected(1,:), discretize(Beh_selected, 5));
xlabel('Pupil Bin'); ylabel('Calcium');
```

---

## 🎨 Visualization Aesthetics & Publication

### For Presentations
1. **Violin/Box plots**: Professional, clear message
2. **Ridge plots**: Visually striking, memorable
3. **Time series with colored background**: Shows temporal dynamics

### For Papers
1. **Scatter + heatmap**: Shows raw data & density
2. **Statistical plot**: Means ± SEM by bin with p-values
3. **Supplementary**: Individual ROI analyses

### Color Schemes
- **Diverging** (correlation): Blues ↔ Reds (negative ↔ positive)
- **Sequential** (density): Yellow → Red (low → high count)
- **Qualitative** (categories): Distinct colors for each bin

---

## ❓ Troubleshooting

### Q: My scatter plot looks like noise
**A**: 
- Increase bin size for heatmap
- Check if relationship is nonlinear (polynomial fit)
- Look at individual ROIs vs. average
- Consider removing high-noise calcium periods

### Q: Violin plots don't look smooth
**A**: 
- You need enough samples (>50 per bin minimum)
- Increase bin smoothing bandwidth
- Consider wider pupil size bins

### Q: I see strong correlation but no clear trend in violin plots
**A**: 
- Relationship might be nonlinear (U-shaped, V-shaped)
- Try polynomial or spline fits instead of linear
- Correlations can be driven by outliers → check Spearman

### Q: How do I compare between animals?
**A**: 
- **Normalize** each animal's data first: `(x - mean) / std`
- Then combine and re-analyze
- Or: Report correlations for each animal separately → compare distributions

---

## 📚 Recommended Reading & Tools

### Python
- `seaborn.violinplot()` / `seaborn.stripplot()`
- `pandas.cut()` for binning
- `scipy.stats.gaussian_kde()` for KDE

### MATLAB
- `violinplot()` (or use `boxplot()` as fallback)
- `discretize()` for binning
- `ksdensity()` for KDE
- `histogram2()` for 2D histograms

### R (if you want to use it)
```r
ggplot(data, aes(x=pupil, y=calcium)) +
  geom_hex(bins=30) +
  geom_smooth(method="loess", se=TRUE) +
  theme_minimal()
```

---

## 🚀 Next Steps

1. **Load your actual data** and run one visualization
2. **Look for the pattern** (correlation, nonlinearity, etc.)
3. **Quantify** with statistics (correlation, ANOVA by bins)
4. **Interpret** considering your axon type and stimulus
5. **Validate** by checking individual ROIs
