from tardis.io.model.readers.snec.xg_files import read_xg_file, XGData
import numpy as np

SNEC_XG_OUTPUT_QUANTITIES = {
    "radius": {
        "long_name": "Radius",
        "units": "cm",
        "description": "Radial coordinate at each point.",
    },
    "mass": {
        "long_name": "Mass",
        "units": "g",
        "description": "Mass coordinate in grams.",
    },
    "vel": {
        "long_name": "Velocity",
        "units": "cm / s",
        "description": "Fluid velocity.",
    },
    "rho": {"long_name": "Density", "units": "g / cm3", "description": "Mass density."},
    "temp": {
        "long_name": "Temperature",
        "units": "K",
        "description": "Gas temperature.",
    },
    "logT": {
        "long_name": "Log Temperature",
        "units": "dex",
        "description": "Base-10 logarithm of temperature.",
    },
    "tau": {
        "long_name": "Optical Depth",
        "units": "",
        "description": "Integrated optical depth from this point to the outer boundary.",
    },
    "lum": {
        "long_name": "Luminosity",
        "units": "erg / s",
        "description": "Local radiative luminosity (see Eq. 4 of the notes).",
    },
    "p_rad": {
        "long_name": "Radiation Pressure",
        "units": "dyn / cm2",
        "description": "Pressure due to radiation field.",
    },
    "press": {
        "long_name": "Gas Pressure",
        "units": "dyn / cm2",
        "description": "Thermal gas pressure.",
    },
    "E_shell": {
        "long_name": "Shell Energy",
        "units": "erg",
        "description": "Auxiliary shell energy quantity (see Sec. 7 of the notes).",
    },
    "Ni_deposit_function": {
        "long_name": "Ni Deposition Function",
        "units": "",
        "description": "Fractional gamma-ray energy deposition function (Eq. 46 of the notes).",
    },
    "ye": {
        "long_name": "Electron Fraction",
        "units": "",
        "description": "Ratio of free electrons to baryons.",
    },
    "free_electron_frac": {
        "long_name": "Free Electron Fraction",
        "units": "",
        "description": "Local fraction of electrons not bound in atoms.",
    },
    "photosphere_tracer": {
        "long_name": "Photosphere Tracer",
        "units": "",
        "description": "1 at the photosphere position, 0 elsewhere.",
    },
    "time_diff": {
        "long_name": "diffusion timescale",
        "units": "s",
        "description": "Auxiliary time difference quantity (used in analysis; Sec. 7 of the notes).",
    },
    "delta_time": {
        "long_name": "Timestep",
        "units": "s",
        "description": "Local timestep as computed in timestep.F90 (Eq. 34 of the notes).",
    },
    "time_exp": {
        "long_name": "Time Since Explosion",
        "units": "s",
        "description": "Time since explosion at each grid point (used in analysis; Sec. 7 of the notes).",
    },
    "Q": {
        "long_name": "Artificial Viscosity",
        "units": "dyn / cm2",
        "description": "Artificial viscosity term added in the momentum equation.",
    },
    "kappa": {
        "long_name": "Opacity (floored)",
        "units": "cm2 / g",
        "description": "Rosseland-mean opacity after applying the floor (see Sec. 3.3).",
    },
    "kappa_table": {
        "long_name": "Opacity (raw)",
        "units": "cm2 / g",
        "description": "Rosseland-mean opacity before applying the floor.",
    },
    "eps": {
        "long_name": "Internal Energy",
        "units": "erg / g",
        "description": "Specific internal energy.",
    },
    "logR_op": {
        "long_name": "Log R",
        "units": "dex",
        "description": "Logarithm of R as defined in opacity calculations (see Sec. 3.3).",
    },
    "cs2": {
        "long_name": "Sound Speed Squared",
        "units": "cm2 / s2",
        "description": "Square of the adiabatic sound speed.",
    },
    "H_1": {
        "long_name": "Neutral Hydrogen Fraction",
        "units": "",
        "description": "Fraction of hydrogen in the neutral state.",
    },
    "H_2": {
        "long_name": "Ionized Hydrogen Fraction",
        "units": "",
        "description": "Fraction of hydrogen in the first ionized state.",
    },
    "He_1": {
        "long_name": "Neutral Helium Fraction",
        "units": "",
        "description": "Fraction of helium in the neutral state.",
    },
    "He_2": {
        "long_name": "Singly Ionized Helium Fraction",
        "units": "",
        "description": "Fraction of helium in the first ionized state.",
    },
    "He_3": {
        "long_name": "Doubly Ionized Helium Fraction",
        "units": "",
        "description": "Fraction of helium in the second ionized state.",
    },
}

def read_snec_output_xg(snec_dir, show_progress=False):
    all_xg_data = {}
    mass_xg_data = read_xg_file(
        snec_dir / "output" / "mass.xg", column_names=["radius", "enclosed_mass"]
    )

    for output_quantity, metadata in SNEC_XG_OUTPUT_QUANTITIES.items():
        xg_file = snec_dir / "output" / f"{output_quantity}.xg"
        if not xg_file.exists():
            raise FileNotFoundError(f"File {xg_file} does not exist.")

        xg_data = read_xg_file(
            xg_file,
            column_names=["enclosed_mass", output_quantity],
            show_progress=show_progress,
        )

        # Ensure timestamps match
        assert np.all(np.isclose(mass_xg_data.timestamps, xg_data.timestamps)), (
            f"Time stamps do not match for {output_quantity} and mass."
        )

        # Add metadata to XGData
        xg_data.metadata = metadata
        all_xg_data[output_quantity] = xg_data

    merged_data_blocks = []
    for i, time_stamp in enumerate(mass_xg_data.timestamps):
        merged_df = mass_xg_data.data_blocks[i][["radius"]].copy()
        merged_df["enclosed_mass"] = mass_xg_data.data_blocks[i]["enclosed_mass"]

        for quantity_name, xg_data in all_xg_data.items():
            merged_df[quantity_name] = xg_data.data_blocks[i][quantity_name]

        merged_data_blocks.append(merged_df)

    # Create a single XGData object with all merged data
    merged_xg_data = XGData(
        timestamps=mass_xg_data.timestamps,
        data_blocks=merged_data_blocks,
        metadata=SNEC_XG_OUTPUT_QUANTITIES  # Use the full SNEC_XG_OUTPUT_QUANTITIES as metadata
    )
    return merged_xg_data
