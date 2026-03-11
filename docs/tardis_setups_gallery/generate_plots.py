"""
Run a TARDIS simulation from a YAML configuration file,
generate Spectrum, SDEC, and LIV plots, and save them as PNG files.

Usage:
    python generate_plots.py --yml /path/to/config.yml --output /output/directory
"""

import argparse
import os
import sys

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt

from tardis import run_tardis
from tardis.visualization import LIVPlotter, SDECPlotter


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate TARDIS plots from a YAML configuration file."
    )
    parser.add_argument(
        "--yml",
        type=str,
        required=True,
        help="Path to the TARDIS YAML configuration file.",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="outputs",
        help="Output folder for plots.",
    )
    return parser.parse_args()


def run_simulation(yml_path):
    print(f"\n[1/4] Running TARDIS simulation from: {yml_path}")
    try:
        sim = run_tardis(
            yml_path,
            virtual_packet_logging=True,
        )
    except Exception as e:
        print(f"Error: Simulation failed: {e}")
        sys.exit(1)
    print("      Simulation complete!")
    return sim


def generate_spectrum_plot(sim, output_dir):
    print("\n[2/4] Generating Spectrum plot...")
    try:
        fig, ax = plt.subplots(figsize=(10, 6))

        spectrum = sim.spectrum_solver.spectrum_real_packets
        spectrum_virtual = sim.spectrum_solver.spectrum_virtual_packets

        ax.plot(
            spectrum.wavelength.value,
            spectrum.luminosity_density_lambda.value,
            label="Real packets",
            color="steelblue",
        )
        ax.plot(
            spectrum_virtual.wavelength.value,
            spectrum_virtual.luminosity_density_lambda.value,
            label="Virtual packets",
            color="darkorange",
            linestyle="--",
        )
        ax.set_xlabel("Wavelength (A)")
        ax.set_ylabel("Luminosity Density (erg/s/A)")
        ax.set_title("TARDIS Synthetic Spectrum")
        ax.legend()
        ax.grid(True, alpha=0.3)
        output_path = os.path.join(output_dir, "spectrum_plot.png")
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"      Saved -> {output_path}")
    except Exception as e:
        print(f"      Warning: Spectrum plot failed: {e}")
        plt.close("all")


def generate_sdec_plot(sim, output_dir):
    print("\n[3/4] Generating SDEC plot...")
    try:
        plotter = SDECPlotter.from_simulation(sim)
        ax = plotter.generate_plot_mpl(packets_mode="real")
        output_path = os.path.join(output_dir, "sdec_plot.png")
        ax.get_figure().savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close("all")
        print(f"      Saved -> {output_path}")
    except Exception as e:
        print(f"      Warning: SDEC plot failed: {e}")
        plt.close("all")


def generate_liv_plot(sim, output_dir):
    print("\n[4/4] Generating LIV plot...")
    try:
        plotter = LIVPlotter.from_simulation(sim)
        ax = plotter.generate_plot_mpl(packets_mode="real")
        output_path = os.path.join(output_dir, "liv_plot.png")
        ax.get_figure().savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close("all")
        print(f"      Saved -> {output_path}")
    except Exception as e:
        print(f"      Warning: LIV plot failed: {e}")
        plt.close("all")


def main():
    args = parse_arguments()
    if not os.path.exists(args.yml):
        print(f"Error: YAML file not found: {args.yml}")
        sys.exit(1)
    os.makedirs(args.output, exist_ok=True)
    sim = run_simulation(args.yml)
    generate_spectrum_plot(sim, args.output)
    generate_sdec_plot(sim, args.output)
    generate_liv_plot(sim, args.output)
    print("\nAll plots generated successfully!")
    print(f"Find your plots in: {args.output}/")


if __name__ == "__main__":
    main()
    