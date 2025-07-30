import pandas as pd

from tardis.energy_input.energy_source import get_radioactive_isotopes


def process_decay_radiation_data(
    decay_radiation_data: pd.DataFrame, isotopic_mass_fraction_index: pd.Index
) -> pd.DataFrame:
    """
    Extract and process radiation data from decay radiation data.

    Extracts and processes radiation data for gamma rays and beta particles from
    decay radiation data, filtering based on relevant radioactive isotopes and
    calculating energy per decay.

    Parameters
    ----------
    decay_radiation_data : pd.DataFrame
        DataFrame containing decay radiation data with columns including 'Z', 'A',
        'Radiation', 'Rad Energy', 'Rad subtype', and 'Rad Intensity'.
    isotopic_mass_fraction_index : pd.Index
        Index of isotopic mass fractions used to filter relevant radioactive isotopes.

    Returns
    -------
    pd.DataFrame
        DataFrame containing processed radiation data with columns:
        'radiation_energy_kev', 'energy_per_decay_kev', and 'radiation_type'.
        Index includes 'atomic_number', 'mass_number', and 'channel_id'.

    Notes
    -----
    The function performs the following operations:
    - Filters decay radiation data for relevant radioactive isotopes
    - Renames columns to standardized format
    - Calculates energy per decay from radiation intensity
    - Adds channel_id as additional multi-index level
    """
    decay_radiation_data = decay_radiation_data.set_index(["Z", "A"])
    decay_radiation_data.index.names = ["atomic_number", "mass_number"]
    relevant_decay_radiation_data = decay_radiation_data.loc[
        get_radioactive_isotopes(isotopic_mass_fraction_index)
    ]
    relevant_decay_radiation_data = relevant_decay_radiation_data.rename(
        columns={
            "Radiation": "radiation_type",
            "Rad Energy": "radiation_energy_kev",
            "Rad subtype": "radiation_sub_type",
        }
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

    return relevant_decay_radiation_data[
        ["radiation_energy_kev", "energy_per_decay_kev", "radiation_type"]
    ]
