"""
Run a TARDIS simulation from a YAML configuration file,
generate Spectrum, SDEC, and LIV plots, and save them as PNG or PDF files.

Usage:
    python generate_plots.py config.yml
    python generate_plots.py config.yml --atom-data kurucz.h5 --format pdf --output my_plots
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
        "config",
        type=str,
        help="Path to the TARDIS YAML configuration file.",
    )
    parser.add_argument(
        "--atom-data",
        type=str,
        default=None,
        help="Path to the atomic dataset (.h5).",
    )
    parser.add_argument(
        "--format",
        type=str,
        choices=["png", "pdf"],
        default="png",
        help="Output format: png or pdf.",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="outputs",
        help="Output folder for plots.",
    )
    return parser.parse_args()


def run_simulation(config_path, atom_data=None):
    print(f"\n[1/4] Running TARDIS simulation from: {config_path}")
    try:
        kwargs = {}
        if atom_data is not None:
            kwargs["atom_data"] = atom_data
        sim = run_tardis(
            config_path,
            virtual_packet_logging=True,
            **kwargs,
        )
    except Exception as e:
        print(f"Error: Simulation failed: {e}")
        sys.exit(1)
    print("      Simulation complete!")
    return sim


def generate_spectrum_plot(sim, output_dir, fmt="png"):
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
        output_path = os.path.join(output_dir, f"spectrum_plot.{fmt}")
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"      Saved -> {output_path}")
    except Exception as e:
        print(f"      Warning: Spectrum plot failed: {e}")
        plt.close("all")


def generate_sdec_plot(sim, output_dir, fmt="png"):
    print("\n[3/4] Generating SDEC plots...")
    try:
        plotter = SDECPlotter.from_simulation(sim)
        for mode in ["real", "virtual"]:
            print(f"\n  Generating SDEC plot [packets_mode={mode}]...")
            ax = plotter.generate_plot_mpl(packets_mode=mode)
            output_path = os.path.join(output_dir, f"sdec_{mode}.{fmt}")
            ax.get_figure().savefig(output_path, dpi=150, bbox_inches="tight")
            plt.close("all")
            print(f"      Saved -> {output_path}")
    except Exception as e:
        print(f"      Warning: SDEC plot failed: {e}")
        plt.close("all")


def generate_liv_plot(sim, output_dir, fmt="png"):
    print("\n[4/4] Generating LIV plots...")
    try:
        plotter = LIVPlotter.from_simulation(sim)
        for mode in ["real", "virtual"]:
            print(f"\n  Generating LIV plot [packets_mode={mode}]...")
            ax = plotter.generate_plot_mpl(packets_mode=mode)
            output_path = os.path.join(output_dir, f"liv_{mode}.{fmt}")
            ax.get_figure().savefig(output_path, dpi=150, bbox_inches="tight")
            plt.close("all")
            print(f"      Saved -> {output_path}")
    except Exception as e:
        print(f"      Warning: LIV plot failed: {e}")
        plt.close("all")


def main():
    args = parse_arguments()
    if not os.path.exists(args.config):
        print(f"Error: YAML file not found: {args.config}")
        return
    os.makedirs(args.output, exist_ok=True)
    sim = run_simulation(args.config, atom_data=args.atom_data)
    generate_spectrum_plot(sim, args.output, fmt=args.format)
    generate_sdec_plot(sim, args.output, fmt=args.format)
    generate_liv_plot(sim, args.output, fmt=args.format)
    print("\nAll plots generated successfully!")
    print(f"Find your plots in: {args.output}/")


if __name__ == "__main__":
    main()
    