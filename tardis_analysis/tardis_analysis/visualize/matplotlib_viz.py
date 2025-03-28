# tardis_analysis/visualize/matplotlib_viz.py

import matplotlib.pyplot as plt
import os
from ..config import MATPLOTLIB_COLORS


def plot_all_commits_matplotlib(all_data, spectrum_keys, output_dir):
    fig, axes = plt.subplots(2, 2, figsize=(20, 15), sharex=True, sharey=True)
    axes = axes.flatten()

    for idx, key in enumerate(spectrum_keys):
        ax = axes[idx]
        for commit_idx, data in enumerate(all_data, 1):
            if key in data:
                wavelength = data[key]['wavelength']
                luminosity = data[key]['luminosity']
                color = MATPLOTLIB_COLORS[(commit_idx - 1) % len(MATPLOTLIB_COLORS)]
                ax.plot(wavelength, luminosity, label=f'Commit {commit_idx}', color=color)
        ax.set_title(f'Luminosity for {key}')
        ax.set_xlabel('Wavelength')
        ax.set_ylabel('Luminosity')
        ax.legend()
        ax.grid(True)

    plt.suptitle('Comparison of Spectrum Solvers Across All Commits', fontsize=16)
    plt.tight_layout()
    file_name = os.path.join(output_dir, "all_commits_spectrum_comparison.pdf")
    try:
        plt.savefig(file_name)
        print(f"Saved Matplotlib plot as {file_name}")
    except Exception as e:
        print(f"Failed to save {file_name}: {e}")
    plt.close()