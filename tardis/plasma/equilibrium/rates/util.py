import pandas as pd


def get_ground_state_multi_index(multi_index_full: pd.MultiIndex):
    """Return the ground-state index of the next ionization stage."""
    atomic_number = multi_index_full.get_level_values(0)
    ion_number = multi_index_full.get_level_values(1) + 1
    level_number = ion_number * 0
    return pd.MultiIndex.from_arrays([atomic_number, ion_number, level_number])


def get_ion_multi_index(
    multi_index_full: pd.MultiIndex, next_higher: bool = True
):
    """Convert a level index to its corresponding ion index."""
    atomic_number = multi_index_full.get_level_values(0)
    ion_number = multi_index_full.get_level_values(1)
    if next_higher:
        ion_number = ion_number + 1
    return pd.MultiIndex.from_arrays([atomic_number, ion_number])


def reindex_ionization_rate_dataframe(
    rate_dataframe: pd.DataFrame, recombination=False
):
    """Add ion source and destination numbers the
    ionization rate dataframe.

    Parameters
    ----------
    rate_dataframe : pd.DataFrame
        Dataframe of ionization rates, indexed by atomic number, ion number,
        level source number and level destination number.
    recombination : bool, optional
        If the rates are recombination rates, by default False

    Returns
    -------
    pd.DataFrame
        Dataframe with additional columns for ion source and destination
    """
    rate_dataframe.index.names = [
        "atomic_number",
        "ion_number",
        "level_number_source",
    ]

    rate_dataframe = rate_dataframe.reset_index()

    if recombination:
        rate_dataframe["ion_number_destination"] = rate_dataframe["ion_number"]
        rate_dataframe["ion_number_source"] = rate_dataframe["ion_number"] + 1
    else:
        rate_dataframe["ion_number_source"] = rate_dataframe["ion_number"]
        rate_dataframe["ion_number_destination"] = (
            rate_dataframe["ion_number"] + 1
        )

    # ionized electrons are assumed to leave the ion in the ground state for now
    rate_dataframe["level_number_destination"] = 0

    not_fully_ionized_mask = (
        rate_dataframe["atomic_number"] != rate_dataframe["ion_number"]
    )

    rate_dataframe = rate_dataframe[not_fully_ionized_mask]

    rate_dataframe = rate_dataframe.set_index(
        [
            "atomic_number",
            "ion_number",
            "ion_number_source",
            "ion_number_destination",
            "level_number_source",
            "level_number_destination",
        ]
    )

    return rate_dataframe
