import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
from ..config import PLOTLY_COLORS
from .utils import calculate_residuals

def plot_residuals_plotly(all_data, spectrum_keys, output_dir, commit_hashes=None, reference_index=0):
    
    if not all_data:
        print("No data to plot.")
        return

   
    fig = make_subplots(rows=2, cols=2, subplot_titles=[f'Fractional Residuals for {key}' for key in spectrum_keys])

    for idx, key in enumerate(spectrum_keys):
        row = idx // 2 + 1
        col = idx % 2 + 1
       
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
                    
                color = PLOTLY_COLORS[commit_idx % len(PLOTLY_COLORS)]
                
               
                commit_label = f'{commit_hashes[commit_idx][:6]}' if commit_hashes else f'Commit {commit_idx + 1}'
                
                fig.add_trace(
                    go.Scatter(
                        x=wavelength,
                        y=residuals,
                        mode='lines',
                        name=commit_label,
                        legendgroup=commit_label,
                        showlegend=(idx == 0), 
                        line=dict(color=color)
                    ),
                    row=row,
                    col=col
                )
        
        fig.add_hline(y=0, line=dict(color='black', dash='dash', width=0.8), row=row, col=col)

    fig.update_layout(
        title='Fractional Residuals Across All Commits',
        height=800,
        showlegend=True
    )

    file_name = os.path.join(output_dir, "residuals_all_commits.html")
    try:
        fig.write_html(file_name)
        print(f"Saved Plotly residual plot as {file_name}")
    except Exception as e:
        print(f"Failed to save {file_name}: {e}")
    return fig