# TARDIS Plot Generator

This contribution introduces a small utility script that runs a TARDIS simulation 
from a YAML configuration file and automatically generates commonly used diagnostic plots.

The goal of this tool is to provide a simple way to explore simulation outputs and 
visualize important properties of the synthetic spectrum.

## Features

The script generates three plots from a TARDIS simulation:

- **Synthetic Spectrum** — Shows luminosity density as a function of wavelength for real and virtual packets.
- **SDEC Plot** — Visualizes packet interaction processes such as electron scattering and line interactions.
- **LIV Plot** — Displays the last interaction velocity distribution of packets by element.

All plots are saved automatically as PNG files.

## Files Included

- `generate_plots.py` — Command-line script that runs the simulation and produces plots.
- `tardis_example.yml` — Example configuration file used to run a minimal TARDIS simulation.

## Requirements

- Python
- TARDIS
- matplotlib

Make sure TARDIS and its dependencies are installed before running the script:
```bash
pip install tardis matplotlib
```

## Usage

Run the script from the command line:
```bash
python generate_plots.py --yml tardis_example.yml
```

Optional: specify a custom output directory (default is `outputs/`):
```bash
python generate_plots.py --yml tardis_example.yml --output my_outputs
```

## Output

The script generates the following files inside the output directory:
```
outputs/
├── spectrum_plot.png
├── sdec_plot.png
└── liv_plot.png
```

These plots provide a quick visual overview of the simulation results.

## Example Workflow

1. Provide a TARDIS YAML configuration file.
2. Run the script.
3. Inspect the generated plots in the output directory.

This utility is intended as a lightweight exploration tool for users working with 
TARDIS simulations and may serve as a starting point for building automated plot galleries.