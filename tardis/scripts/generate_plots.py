"""
Run a TARDIS simulation from a YAML configuration and an atom-data file,
generate the SDEC and LIV plots using the built-in plotter classes, and
save them as PNG files.

Usage:
    python generate_plots.py --config /path/to/config.yml --atom_data /path/to/atom_data.h5
        --output_path /output/directory
"""

import argparse
from pathlib import Path

from tardis import run_tardis
from tardis.visualization import SDECPlotter, LIVPlotter


def save_sdec_plot(sim, output_path):
    sdec = SDECPlotter.from_simulation(sim)
    plot = sdec.generate_plot_mpl(packets_mode="real").figure
    plot.savefig(Path(output_path)/"SDEC.png", dpi=300, bbox_inches="tight")


def save_liv_plot(sim, output_path):
    liv = LIVPlotter.from_simulation(sim)
    plot = liv.generate_plot_mpl(packets_mode="real").figure
    plot.savefig(Path(output_path)/"LIV.png", dpi=300, bbox_inches="tight")


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="Run TARDIS simulation from YAML and atom-data and save SDEC & LIV plots as PNG."
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to TARDIS YAML configuration file."
    )
    parser.add_argument(
        "--atom_data",
        default=None,
        help="Path to TARDIS atom-data HDF5 file (optional, can be specified in config file)."
    )
    parser.add_argument(
        "--output_path",
        default=None,
        help="Output directory where PNG files will be saved. Default is a new directory with same name as config file.",
    )

    return parser.parse_args(argv)


def main():
    args = parse_arguments()

    if args.output_path is None:
        output_path = Path(args.config).parent / Path(args.config).stem
    else:
        output_path = Path(args.output_path)
    output_path.mkdir(exist_ok=True)

    sim = run_tardis(
        config=args.config,
        atom_data=args.atom_data,
        virtual_packet_logging=True,
    )

    save_sdec_plot(sim, output_path)
    save_liv_plot(sim, output_path)


if __name__ == "__main__":
    main()