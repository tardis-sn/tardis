from pathlib import Path
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from .config import PLOTLY_COLORS, SPECTRUM_KEYS
from .data_processing import calculate_residuals

def plot_combined_analysis_plotly(all_data, spectrum_keys, output_dir, commit_hashes=None, reference_index=0):
    if not all_data:
        print("No data to plot.")
        return
    fig = make_subplots(
        rows=4, 
        cols=2, 
        subplot_titles=[
            'spectrum_integrated', 'spectrum_real_packets',
            'Residuals - spectrum_integrated', 'Residuals - spectrum_real_packets',
            'spectrum_real_packets_reabsorbed', 'spectrum_virtual_packets',
            'Residuals - spectrum_real_packets_reabsorbed', 'Residuals - spectrum_virtual_packets'
        ]
    )

    for idx, key in enumerate(spectrum_keys):
        if idx < 2:
            spectrum_row = 1  
            residual_row = 2  
        else:
            spectrum_row = 3  
            residual_row = 4  
        
        spectrum_col = idx % 2 + 1 
        for commit_idx, data in enumerate(all_data):
            if key not in data:
                continue

            wavelength = data[key]['wavelength']
            luminosity = data[key]['luminosity']
            color = PLOTLY_COLORS[commit_idx % len(PLOTLY_COLORS)]
            commit_label = f'{commit_hashes[commit_idx][:6]}' if commit_hashes else f'Commit {commit_idx + 1}'

            fig.add_trace(
                go.Scatter(
                    x=wavelength,
                    y=luminosity,
                    mode='lines',
                    name=commit_label,
                    legendgroup=commit_label,
                    showlegend=(idx == 0), 
                    line=dict(color=color)
                ),
                row=spectrum_row,
                col=spectrum_col
            )
            if commit_idx != reference_index: 
                ref_wavelength = all_data[reference_index][key]['wavelength']
                ref_luminosity = all_data[reference_index][key]['luminosity']
                
                wavelength, residuals, is_valid = calculate_residuals(
                    wavelength, luminosity, ref_wavelength, ref_luminosity
                )
                if is_valid:
                    fig.add_trace(
                        go.Scatter(
                            x=wavelength,
                            y=residuals,
                            mode='lines',
                            name=commit_label,
                            legendgroup=commit_label,
                            showlegend=False, 
                            line=dict(color=color)
                        ),
                        row=residual_row,
                        col=spectrum_col
                    )

        fig.add_hline(
            y=0, 
            line=dict(color='black', dash='dash', width=0.8),
            row=residual_row,
            col=spectrum_col
        )

    fig.update_layout(
        height=1200,  
        showlegend=True,
        title='Spectrum Analysis with Residuals'
    )

    for i in range(4):
        row = i + 1
        for col in [1, 2]:
            if row % 2 == 1:  
                fig.update_yaxes(title_text="Luminosity (erg/s/Å)", row=row, col=col)
            else:  
                fig.update_yaxes(title_text="Residual (%)", row=row, col=col)
            fig.update_xaxes(title_text="Wavelength (Å)", row=row, col=col)

    output_path = Path(output_dir)
    file_name = output_path / "combined_analysis.html"
    
    try:
        file_name.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(str(file_name))
        print(f"Saved combined analysis plot as {file_name}")
    except Exception as e:
        print(f"Failed to save {file_name}: {e}")
    
    return fig 