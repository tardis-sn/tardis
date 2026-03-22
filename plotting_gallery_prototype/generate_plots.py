import argparse
import os
import matplotlib.pyplot as plt

import tardis
from tardis import run_tardis

def generate_plots(yml_path, atom_data=None):
    print(f"Running TARDIS v{tardis.__version__} with {yml_path}")
    if atom_data:
        sim = run_tardis(yml_path, atom_data=atom_data)
    else:
        sim = run_tardis(yml_path)
    
    # SDec Plot
    try:
        from tardis.visualization import SDECPlotter
        plotter = SDECPlotter.from_simulation(sim)
        ax_sdec = plotter.generate_plot_mpl(packets_mode='real')
        sdec_out = "sdec_plot.pdf"
        ax_sdec.figure.savefig(sdec_out, bbox_inches='tight')
        print(f"SDec plot successfully saved to {sdec_out}")
    except Exception as e:
        print("SDecPlotter failed:", type(e).__name__, e)
        print("SDecPlotter failed:", type(e).__name__, e)
        # Fallback
        fig, ax = plt.subplots(figsize=(10, 6))
        try:
            spectrum = sim.spectrum_solver.spectrum
        except AttributeError:
            try:
                spectrum = sim.transport.spectrum
            except AttributeError:
                spectrum = sim.runner.spectrum
        ax.plot(spectrum.wavelength, spectrum.luminosity_density_lambda)
        ax.set_xlabel("Wavelength [Å]")
        ax.set_ylabel("Luminosity Density [erg/s/Å]")
        ax.set_title("TARDIS Spectrum")
        fig.savefig("sdec_fallback_plot.pdf")
        print("Saved fallback spectrum plot.")

    # LIV Plot
    try:
        if hasattr(tardis.visualization, 'LIVPlotter'):
            from tardis.visualization import LIVPlotter
            liv_plotter = LIVPlotter.from_simulation(sim)
            fig_liv = liv_plotter.generate_plot()
        else:
            # support legacy LineInfoWidget
            from tardis.visualization import LineInfoWidget
            liv_plotter = LineInfoWidget.from_simulation(sim)
            fig_liv = liv_plotter.figure_widget
            
        liv_out = "liv_plot.pdf"
        if hasattr(fig_liv, 'write_image'):
            try:
                fig_liv.write_image(liv_out)
                print(f"Saved LIV plot to {liv_out}")
            except Exception as e:
                print(f"PDF export failed ({e}). Saving as HTML.")
                html_out = "liv_plot.html"
                fig_liv.write_html(html_out)
                print(f"Saved interactive map to {html_out}")
        else:
            fig_liv.savefig(liv_out)
            print(f"Saved LIV plot to {liv_out}")
    except Exception as e:
        print(f"Could not generate LIV Plot: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate TARDIS SDec and LIV plots.")
    parser.add_argument("config_yaml", help="Path to the TARDIS YAML config file")
    parser.add_argument("atom_data", help="Path to the atomic data file (optional)", nargs="?", default=None)
    args = parser.parse_args()
    
    generate_plots(args.config_yaml, args.atom_data)
