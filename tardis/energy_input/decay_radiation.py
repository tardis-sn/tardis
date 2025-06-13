import pandas as pd
import xarray as xr

from tardis.energy_input.energy_source import get_radioactive_isotopes


def get_decay_radiation_data(decay_radiation_data, isotopic_mass_fraction_index):
    """
    Extracts and processes radiation data for gamma rays and beta particles from decay radiation data.

    Parameters
    ----------
    decay_radiation_data : pd.DataFrame
        DataFrame containing decay radiation data with columns including 'Z', 'A', 'Radiation',
        'Rad Energy', 'Rad subtype', and 'Rad Intensity'.
    isotopic_mass_fraction_index : pd.Index
        Index of isotopic mass fractions used to filter relevant radioactive isotopes.

    Returns
    -------
    em_radiation_data : pd.DataFrame
        DataFrame containing processed electromagnetic radiation data for gamma rays with columns:
        'radiation_energy_kev', 'energy_per_decay_kev', and 'radiation_type'.
    bp_radiation_data : pd.DataFrame
        DataFrame containing processed radiation data for beta particles with columns:
        'radiation_energy_kev', 'energy_per_decay_kev', and 'radiation_type'.
    """
    decay_radiation_data = decay_radiation_data.set_index(["Z", "A"])
    decay_radiation_data.index.names = ["atomic_number", "mass_number"]
    relevant_decay_radiation_data = decay_radiation_data.loc[
        get_radioactive_isotopes(isotopic_mass_fraction_index)
    ]
    relevant_decay_radiation_data = relevant_decay_radiation_data.rename(
        columns={"Rad Energy": "radiation_energy_kev", "Rad subtype": "radiation_type"}
    )
    relevant_decay_radiation_data[
        "energy_per_decay_kev"
    ] = relevant_decay_radiation_data["radiation_energy_kev"] * (
        relevant_decay_radiation_data["Rad Intensity"] / 100
    )  # given per 100 decays

    # Add channel_id as an additional multi_index
    relevant_decay_radiation_data = relevant_decay_radiation_data.reset_index()
    relevant_decay_radiation_data["channel_id"] = relevant_decay_radiation_data.groupby(
        ["atomic_number", "mass_number"]
    ).cumcount()
    relevant_decay_radiation_data = relevant_decay_radiation_data.set_index(
        ["atomic_number", "mass_number", "channel_id"]
    )

    em_radiation_data = relevant_decay_radiation_data[
        relevant_decay_radiation_data.Radiation == "g"
    ]
    # move 'radiation_type' to the last column
    em_radiation_data = em_radiation_data[
        ["radiation_energy_kev", "energy_per_decay_kev", "radiation_type"]
    ]

    bp_radiation_data = relevant_decay_radiation_data[
        relevant_decay_radiation_data.Radiation == "bp"
    ]
    bp_radiation_data = bp_radiation_data[
        ["radiation_energy_kev", "energy_per_decay_kev", "radiation_type"]
    ]

    return em_radiation_data, bp_radiation_data


def convert_decay_radiation_to_ds(radiation_data: pd.DataFrame) -> xr.Dataset:
    """
    Convert decay radiation DataFrame to xarray Dataset.

    Parameters
    ----------
    radiation_data : pd.DataFrame
        DataFrame with MultiIndex (atomic_number, mass_number, channel_id)
        and columns like 'energy_per_decay_kev', 'radiation_energy_kev', 'radiation_type'

    Returns
    -------
    xr.Dataset
        Dataset with decay_radiation_channel_id as main dimension
    """
    dataset = xr.Dataset(
        {
            "energy_per_decay_kev": (
                "decay_radiation_channel_id",
                radiation_data["energy_per_decay_kev"].values,
            ),
            "radiation_energy_kev": (
                "decay_radiation_channel_id",
                radiation_data["radiation_energy_kev"].values,
            ),
            "radiation_type": (
                "decay_radiation_channel_id",
                radiation_data["radiation_type"].values,
            ),
        },
        coords={
            "atomic_number": (
                "decay_radiation_channel_id",
                pd.Index(radiation_data.index.get_level_values("atomic_number")),
            ),
            "mass_number": (
                "decay_radiation_channel_id",
                pd.Index(radiation_data.index.get_level_values("mass_number")),
            ),
            "decay_radiation_channel_id": (
                "decay_radiation_channel_id",
                pd.Index(radiation_data.index.get_level_values("channel_id")),
            ),
        },
    )
    return dataset
