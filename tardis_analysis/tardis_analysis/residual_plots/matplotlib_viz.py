import matplotlib.pyplot as plt
import os
from ..config import MATPLOTLIB_COLORS
from .utils import calculate_residuals

def plot_residuals_matplotlib(all_data, spectrum_keys, output_dir, reference_index=0):
    if not all_data:
        print("No data to plot.")
        return

    fig, axes = plt.subplots(2, 2, figsize=(20, 15), sharex=True, sharey=True)
    axes = axes.flatten()

    for idx, key in enumerate(spectrum_keys):
        ax = axes[idx]
        
        if reference_index >= len(all_data) or key not in all_data[reference_index]:
            print(f"Reference commit does not have data for {key}. Skipping.")
            continue

        ref_wavelength = all_data[reference_index][key]['wavelength']
        ref_luminosity = all_data[reference_index][key]['luminosity']


        for commit_idx, data in enumerate(all_data):
            if key in data:
                wavelength = data[key]['wavelength']
                luminosity = data[key]['luminosity']
                
                wavelength, residuals, is_valid = calculate_residuals(
                    wavelength, luminosity, ref_wavelength, ref_luminosity
                )
                
                if not is_valid:
                    print(f"Wavelengths differ for commit {commit_idx} and reference for {key}. Skipping.")
                    continue
                    
                color = MATPLOTLIB_COLORS[commit_idx % len(MATPLOTLIB_COLORS)]
                ax.plot(wavelength, residuals, label=f'Commit {commit_idx + 1}', color=color)


        ax.axhline(0, color='black', linestyle='--', linewidth=0.8)
        ax.set_title(f'Fractional Residuals for {key}')
        ax.set_xlabel('Wavelength')
        ax.set_ylabel('Fractional Residual')
        ax.legend()
        ax.grid(True)

    plt.suptitle('Fractional Residuals Across All Commits', fontsize=16)
    plt.tight_layout()
    file_name = os.path.join(output_dir, "residuals_all_commits.pdf")
    try:
        plt.savefig(file_name)
        print(f"Saved Matplotlib residual plot as {file_name}")
    except Exception as e:
        print(f"Failed to save {file_name}: {e}")
    plt.close()