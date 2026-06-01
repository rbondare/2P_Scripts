#!/usr/bin/env python3
"""
Pupil-Calcium Activity Relationship Visualization
Multiple plotting approaches to explore how calcium activity changes across pupil sizes
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Assuming you have loaded the data:
# data = scipy.io.loadmat('your_file.mat')
# pupil_area = data['Behav']['PupilArea']
# ca_activity = data['CaData'][0]['Ca_dFF']  # or dFF_z for normalized

def prepare_data(pupil_area, ca_activity, time_vector=None):
    """
    Prepare and align pupil and calcium data
    
    Parameters:
    -----------
    pupil_area : 1D array of pupil sizes
    ca_activity : 2D array (n_rois, n_timepoints) or 1D array (n_timepoints,)
    time_vector : 1D array of timestamps (optional)
    
    Returns:
    --------
    pupil_clean, ca_clean : aligned and cleaned data
    """
    # Flatten if 2D
    if len(ca_activity.shape) > 1:
        ca_activity = np.nanmean(ca_activity, axis=0)  # Average across ROIs
    
    # Remove NaNs and align
    valid_idx = ~(np.isnan(pupil_area) | np.isnan(ca_activity))
    pupil_clean = pupil_area[valid_idx]
    ca_clean = ca_activity[valid_idx]
    
    return pupil_clean, ca_clean


def plot_scatter_with_marginals(pupil, ca_activity, title="Pupil-Calcium Relationship"):
    """
    Scatter plot with marginal distributions (KDE on margins)
    Good for: Seeing individual data points and distributions
    """
    fig = plt.figure(figsize=(10, 10))
    gs = fig.add_gridspec(3, 3, hspace=0.05, wspace=0.05)
    
    # Main scatter plot
    ax_main = fig.add_subplot(gs[1:, :-1])
    ax_top = fig.add_subplot(gs[0, :-1], sharex=ax_main)
    ax_right = fig.add_subplot(gs[1:, -1], sharey=ax_main)
    
    # Main scatter with density coloring
    hb = ax_main.hexbin(pupil, ca_activity, gridsize=30, cmap='YlOrRd', mincnt=1, edgecolors='face')
    ax_main.set_xlabel('Pupil Area (px²)', fontsize=12, fontweight='bold')
    ax_main.set_ylabel('Calcium Activity (ΔF/F)', fontsize=12, fontweight='bold')
    
    # Correlation
    corr = np.corrcoef(pupil, ca_activity)[0, 1]
    ax_main.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax_main.transAxes, 
                 fontsize=11, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Top histogram (pupil)
    ax_top.hist(pupil, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    ax_top.set_ylabel('Count', fontsize=10)
    ax_top.tick_params(labelbottom=False)
    ax_top.set_title(title, fontsize=14, fontweight='bold')
    
    # Right histogram (calcium)
    ax_right.hist(ca_activity, bins=50, orientation='horizontal', color='coral', alpha=0.7, edgecolor='black')
    ax_right.set_xlabel('Count', fontsize=10)
    ax_right.tick_params(labelleft=False)
    
    # Add colorbar
    cbar = plt.colorbar(hb, ax=ax_main)
    cbar.set_label('Count', fontsize=10)
    
    return fig


def plot_violin_distribution(pupil, ca_activity, n_bins=5, title="Calcium Distribution by Pupil Size"):
    """
    Violin plot: shows distribution of calcium activity for different pupil sizes
    Good for: Seeing how distributions change across pupil ranges
    """
    # Bin pupil size
    pupil_bins = pd.cut(pupil, bins=n_bins)
    data_for_violin = pd.DataFrame({
        'Calcium Activity': ca_activity,
        'Pupil Size': pupil_bins
    })
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    sns.violinplot(data=data_for_violin, x='Pupil Size', y='Calcium Activity', 
                   palette='Set2', ax=ax, inner='quartile')
    ax.set_xlabel('Pupil Size Range', fontsize=12, fontweight='bold')
    ax.set_ylabel('Calcium Activity (ΔF/F)', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Add sample size
    for i, label in enumerate(ax.get_xticklabels()):
        bin_label = data_for_violin['Pupil Size'].unique()[i]
        n = len(data_for_violin[data_for_violin['Pupil Size'] == bin_label])
        ax.text(i, ax.get_ylim()[0], f'n={n}', ha='center', va='top', fontsize=9)
    
    plt.tight_layout()
    return fig


def plot_2d_density_heatmap(pupil, ca_activity, title="2D Density: Pupil vs Calcium"):
    """
    2D density plot using hexbin
    Good for: Seeing joint distribution and hotspots
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    
    hb = ax.hexbin(pupil, ca_activity, gridsize=40, cmap='YlOrRd', mincnt=1, edgecolors='face')
    ax.set_xlabel('Pupil Area (px²)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Calcium Activity (ΔF/F)', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    cbar = plt.colorbar(hb, ax=ax)
    cbar.set_label('Count', fontsize=11)
    
    # Add trend line
    z = np.polyfit(pupil, ca_activity, 2)
    p = np.poly1d(z)
    pupil_sorted = np.sort(pupil)
    ax.plot(pupil_sorted, p(pupil_sorted), 'r--', linewidth=2, label='Polynomial fit (degree 2)', alpha=0.8)
    ax.legend(fontsize=10)
    
    # Correlation
    corr = np.corrcoef(pupil, ca_activity)[0, 1]
    ax.text(0.05, 0.95, f'Correlation: r = {corr:.3f}', transform=ax.transAxes, 
            fontsize=11, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
    
    return fig


def plot_ridge_distribution(pupil, ca_activity, n_bins=5, title="Ridge Plot: Calcium Distribution Across Pupil Sizes"):
    """
    Ridge plot: stacked density plots for each pupil size bin
    Good for: Beautiful visualization of distribution changes
    """
    from scipy.stats import gaussian_kde
    
    # Bin pupil size
    pupil_min, pupil_max = pupil.min(), pupil.max()
    bin_edges = np.linspace(pupil_min, pupil_max, n_bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    colors = plt.cm.viridis(np.linspace(0, 1, n_bins))
    
    for i in range(n_bins):
        # Get data for this pupil bin
        mask = (pupil >= bin_edges[i]) & (pupil < bin_edges[i+1])
        ca_bin = ca_activity[mask]
        
        if len(ca_bin) > 5:  # Need enough points for KDE
            # Calculate KDE
            kde = gaussian_kde(ca_bin)
            x_range = np.linspace(ca_activity.min(), ca_activity.max(), 200)
            density = kde(x_range)
            
            # Normalize and offset for ridge effect
            density_normalized = density / density.max() * 0.8
            offset = i * 0.5
            
            # Plot
            ax.fill_between(x_range, offset, density_normalized + offset, alpha=0.7, color=colors[i], 
                            label=f'Pupil: {bin_edges[i]:.0f}-{bin_edges[i+1]:.0f}')
            ax.plot(x_range, density_normalized + offset, color=colors[i], linewidth=2)
            
            # Add sample size
            ax.text(ca_activity.min(), offset, f'n={len(ca_bin)}', fontsize=9, va='center')
    
    ax.set_xlabel('Calcium Activity (ΔF/F)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Pupil Size Range', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=10)
    ax.set_ylim(-0.2, n_bins * 0.5)
    
    return fig


def plot_heatmap_binned_statistics(pupil, ca_activity, n_pupil_bins=10, n_ca_bins=10, 
                                   title="Mean Activity Heatmap"):
    """
    Heatmap of mean/median calcium activity binned by pupil size
    Good for: Simple, interpretable overview
    """
    # Create 2D bins
    pupil_bins = np.linspace(pupil.min(), pupil.max(), n_pupil_bins + 1)
    ca_bins = np.linspace(ca_activity.min(), ca_activity.max(), n_ca_bins + 1)
    
    # Create 2D histogram
    H, xedges, yedges = np.histogram2d(pupil, ca_activity, bins=[pupil_bins, ca_bins])
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    im = ax.imshow(H.T, origin='lower', aspect='auto', cmap='hot', 
                   extent=[pupil_bins[0], pupil_bins[-1], ca_bins[0], ca_bins[-1]],
                   interpolation='nearest')
    
    ax.set_xlabel('Pupil Area (px²)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Calcium Activity (ΔF/F)', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Count', fontsize=11)
    
    return fig


def plot_calcium_by_pupil_bins(pupil, ca_activity, n_bins=5, 
                               title="Calcium Activity Statistics by Pupil Size"):
    """
    Box/violin plot with statistics for each pupil bin
    Good for: Comparing statistics across pupil ranges
    """
    pupil_min, pupil_max = pupil.min(), pupil.max()
    bin_edges = np.linspace(pupil_min, pupil_max, n_bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_labels = [f'{bin_edges[i]:.0f}-{bin_edges[i+1]:.0f}' for i in range(n_bins)]
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(title, fontsize=14, fontweight='bold')
    
    # Prepare data for each bin
    ca_by_bin = []
    counts = []
    
    for i in range(n_bins):
        mask = (pupil >= bin_edges[i]) & (pupil < bin_edges[i+1])
        ca_bin = ca_activity[mask]
        ca_by_bin.append(ca_bin)
        counts.append(len(ca_bin))
    
    # Plot 1: Violin plot
    ax = axes[0, 0]
    parts = ax.violinplot([ca_bin for ca_bin in ca_by_bin], positions=range(n_bins), 
                          showmeans=True, showmedians=True)
    ax.set_xticks(range(n_bins))
    ax.set_xticklabels(bin_labels, rotation=45, ha='right')
    ax.set_ylabel('Calcium Activity', fontsize=11, fontweight='bold')
    ax.set_title('Distribution by Pupil Size', fontsize=12, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
    
    # Plot 2: Mean and std
    ax = axes[0, 1]
    means = [np.mean(ca_bin) for ca_bin in ca_by_bin]
    stds = [np.std(ca_bin) for ca_bin in ca_by_bin]
    ax.errorbar(bin_centers, means, yerr=stds, marker='o', markersize=8, capsize=5, 
                capthick=2, linewidth=2, color='steelblue', ecolor='coral')
    ax.set_xlabel('Pupil Size (center)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Mean Calcium Activity', fontsize=11, fontweight='bold')
    ax.set_title('Mean ± Std', fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3)
    
    # Plot 3: Median and IQR
    ax = axes[1, 0]
    medians = [np.median(ca_bin) for ca_bin in ca_by_bin]
    q1 = [np.percentile(ca_bin, 25) for ca_bin in ca_by_bin]
    q3 = [np.percentile(ca_bin, 75) for ca_bin in ca_by_bin]
    ax.errorbar(bin_centers, medians, yerr=[np.array(medians)-np.array(q1), np.array(q3)-np.array(medians)],
                marker='s', markersize=8, capsize=5, capthick=2, linewidth=2, color='green', ecolor='orange')
    ax.set_xlabel('Pupil Size (center)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Median Calcium Activity', fontsize=11, fontweight='bold')
    ax.set_title('Median and IQR', fontsize=12, fontweight='bold')
    ax.grid(alpha=0.3)
    
    # Plot 4: Sample sizes
    ax = axes[1, 1]
    ax.bar(range(n_bins), counts, color='steelblue', alpha=0.7, edgecolor='black')
    ax.set_xticks(range(n_bins))
    ax.set_xticklabels(bin_labels, rotation=45, ha='right')
    ax.set_ylabel('Sample Size (n)', fontsize=11, fontweight='bold')
    ax.set_title('Samples per Bin', fontsize=12, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
    
    # Add sample counts on bars
    for i, count in enumerate(counts):
        ax.text(i, count, str(count), ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    return fig


def print_statistics(pupil, ca_activity, n_bins=5):
    """
    Print summary statistics
    """
    print("\n" + "="*70)
    print("PUPIL-CALCIUM ACTIVITY ANALYSIS")
    print("="*70)
    
    print(f"\nGlobal Statistics:")
    print(f"  Pupil Area: μ={np.mean(pupil):.2f}, σ={np.std(pupil):.2f}")
    print(f"  Calcium Activity: μ={np.mean(ca_activity):.2f}, σ={np.std(ca_activity):.2f}")
    
    # Correlation
    corr, pval = stats.pearsonr(pupil, ca_activity)
    print(f"\nPearson Correlation: r={corr:.4f}, p={pval:.2e}")
    
    # Spearman correlation (non-parametric)
    spearman_r, spearman_p = stats.spearmanr(pupil, ca_activity)
    print(f"Spearman Correlation: ρ={spearman_r:.4f}, p={spearman_p:.2e}")
    
    # By bins
    print(f"\nStatistics by Pupil Size Bins (n={n_bins}):")
    print("-"*70)
    
    pupil_min, pupil_max = pupil.min(), pupil.max()
    bin_edges = np.linspace(pupil_min, pupil_max, n_bins + 1)
    
    for i in range(n_bins):
        mask = (pupil >= bin_edges[i]) & (pupil < bin_edges[i+1])
        ca_bin = ca_activity[mask]
        
        if len(ca_bin) > 0:
            print(f"\n  Bin {i+1}: Pupil {bin_edges[i]:.0f}-{bin_edges[i+1]:.0f} (n={len(ca_bin)})")
            print(f"    Calcium: μ={np.mean(ca_bin):.3f}, σ={np.std(ca_bin):.3f}, " + 
                  f"median={np.median(ca_bin):.3f}, IQR={np.percentile(ca_bin, 75)-np.percentile(ca_bin, 25):.3f}")


# ============================================================================
# MAIN: How to use
# ============================================================================

if __name__ == "__main__":
    print("Pupil-Calcium Activity Visualization Script")
    print("\nUsage:")
    print("------")
    print("1. Load your data (MATLAB or NPY file)")
    print("2. Extract pupil area and calcium activity")
    print("3. Run the plotting functions")
    print("\nExample:")
    print("--------")
    print("""
    import scipy.io
    import pandas as pd
    
    # Load data
    data = scipy.io.loadmat('your_file.mat')
    pupil_area = data['Behav']['PupilArea'][0]
    ca_dff = data['CaData'][0]['Ca_dFF']  # Shape: (n_rois, n_timepoints)
    
    # Prepare data
    pupil_clean, ca_clean = prepare_data(pupil_area, ca_dff)
    
    # Create visualizations
    fig1 = plot_scatter_with_marginals(pupil_clean, ca_clean)
    fig2 = plot_violin_distribution(pupil_clean, ca_clean, n_bins=5)
    fig3 = plot_2d_density_heatmap(pupil_clean, ca_clean)
    fig4 = plot_ridge_distribution(pupil_clean, ca_clean, n_bins=5)
    fig5 = plot_calcium_by_pupil_bins(pupil_clean, ca_clean, n_bins=5)
    
    # Print statistics
    print_statistics(pupil_clean, ca_clean, n_bins=5)
    
    plt.show()
    """)
