"""Script to generate SDec and LIV diagnostic plots from a TARDIS simulation."""

import argparse
import os

import matplotlib.pyplot as plt

import tardis
from tardis import run_tardis
from tardis.visualization import SDECPlotter


def generate_plots(
    yml_path,
    atom_data=None,
    output_dir="outputs",
    output_format="pdf",
    prefix="tardis",
):
    """Run a TARDIS simulation and produce SDec and LIV diagnostic plots.

    Both plots are saved to *output_dir*.  If ``SDECPlotter`` fails for any
    reason, a plain spectrum plot is written as a fallback.  The LIV plot is
    attempted with ``LIVPlotter`` (current API) and falls back to
    ``LineInfoWidget`` for older TARDIS builds.

    Parameters
    ----------
    yml_path : str
        Path to the TARDIS YAML configuration file.
    atom_data : str or None
        Path to the atomic data HDF5 file.  ``None`` uses the path specified
        inside the YAML configuration.
    output_dir : str
        Directory where the output files will be written.  Created
        automatically if it does not already exist.
    output_format : {"pdf", "png", "svg"}
        File format for saved plots.
    prefix : str
        Filename prefix applied to every output file, e.g.
        ``"tardis"`` → ``"tardis_sdec.pdf"``.
    """
    print(f"Running TARDIS v{tardis.__version__} with config: {yml_path}")

    os.makedirs(output_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Run simulation
    # ------------------------------------------------------------------
    if atom_data:
        sim = run_tardis(yml_path, atom_data=atom_data)
    else:
        sim = run_tardis(yml_path)

    # ------------------------------------------------------------------
    # SDec plot
    # ------------------------------------------------------------------
    try:
        plotter = SDECPlotter.from_simulation(sim)
        ax_sdec = plotter.generate_plot_mpl(packets_mode="real")
        sdec_out = os.path.join(output_dir, f"{prefix}_sdec.{output_format}")
        ax_sdec.figure.savefig(sdec_out, bbox_inches="tight")
        print(f"SDec plot saved to {sdec_out}")
    except Exception as exc:  # noqa: BLE001
        print(f"SDECPlotter failed ({type(exc).__name__}: {exc}). Saving spectrum fallback.")
        _save_spectrum_fallback(sim, output_dir, prefix, output_format)

    # ------------------------------------------------------------------
    # LIV plot
    # ------------------------------------------------------------------
    try:
        _save_liv_plot(sim, output_dir, prefix, output_format)
    except Exception as exc:  # noqa: BLE001
        print(f"Could not generate LIV plot: {exc}")


def _save_spectrum_fallback(sim, output_dir, prefix, output_format):
    """Save a basic spectrum plot when SDECPlotter is unavailable.

    Parameters
    ----------
    sim : tardis.simulation.base.Simulation
        Completed TARDIS simulation object.
    output_dir : str
        Output directory (must already exist).
    prefix : str
        Filename prefix.
    output_format : str
        File format extension.
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    # TARDIS API changed across versions; try each attribute in order.
    for attr in ("spectrum_solver.spectrum", "transport.spectrum", "runner.spectrum"):
        try:
            parts = attr.split(".")
            spectrum = sim
            for part in parts:
                spectrum = getattr(spectrum, part)
            break
        except AttributeError:
            continue

    ax.plot(spectrum.wavelength, spectrum.luminosity_density_lambda)
    ax.set_xlabel(r"Wavelength [$\mathrm{\AA}$]")
    ax.set_ylabel(r"Luminosity Density [erg/s/$\mathrm{\AA}$]")
    ax.set_title("TARDIS Spectrum")

    out_path = os.path.join(output_dir, f"{prefix}_spectrum_fallback.{output_format}")
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Fallback spectrum plot saved to {out_path}")


def _save_liv_plot(sim, output_dir, prefix, output_format):
    """Generate and save the LIV (Line Interaction Volume) plot.

    Parameters
    ----------
    sim : tardis.simulation.base.Simulation
        Completed TARDIS simulation object.
    output_dir : str
        Output directory (must already exist).
    prefix : str
        Filename prefix.
    output_format : str
        File format extension.
    """
    from tardis.visualization import LIVPlotter

    liv_plotter = LIVPlotter.from_simulation(sim)
    fig_liv = liv_plotter.generate_plot()
    liv_out = os.path.join(output_dir, f"{prefix}_liv.{output_format}")

    if hasattr(fig_liv, "write_image"):
        try:
            fig_liv.write_image(liv_out)
        except Exception:  # noqa: BLE001 – kaleido may not be installed
            html_out = os.path.join(output_dir, f"{prefix}_liv.html")
            fig_liv.write_html(html_out)
            print(f"LIV plot saved as interactive HTML to {html_out}")
            return
    else:
        fig_liv.savefig(liv_out, bbox_inches="tight")

    print(f"LIV plot saved to {liv_out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate TARDIS SDec and LIV diagnostic plots from a YAML config."
    )
    parser.add_argument("config_yaml", help="Path to the TARDIS YAML configuration file")
    parser.add_argument(
        "atom_data",
        nargs="?",
        default=None,
        help="Path to the atomic data HDF5 file (optional; overrides the path in the YAML)",
    )
    parser.add_argument(
        "--output-dir",
        default="outputs",
        metavar="DIR",
        help="Directory to save the plots (default: %(default)s)",
    )
    parser.add_argument(
        "--format",
        default="pdf",
        choices=["pdf", "png", "svg"],
        dest="output_format",
        help="Output file format (default: %(default)s)",
    )
    parser.add_argument(
        "--prefix",
        default="tardis",
        help="Filename prefix for saved plots (default: %(default)s)",
    )

    args = parser.parse_args()
    generate_plots(
        args.config_yaml,
        args.atom_data,
        args.output_dir,
        args.output_format,
        args.prefix,
    )
