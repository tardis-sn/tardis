import pandas as pd


def reindex_ion_number_density_to_level_number_density(
    ion_number_density: pd.DataFrame,
    level_number_density: pd.DataFrame,
    next_higher: bool = True,
) -> pd.DataFrame:
    """Align ion number density indexes with level number densities by ionization stage.

    Parameters
    ----------
    ion_number_density : pd.DataFrame
        Ion-resolved values indexed by atomic number and ion number.
    level_number_density : pd.DataFrame
        Level-resolved values indexed by atomic number, ion number, and level
        number.
    next_higher : bool, optional
        Select the continuum ion for each level when ``True``; otherwise select
        the ion containing the level.

    Returns
    -------
    pd.DataFrame
        Ion-resolved values with the level number-density index and columns.
    """
    if ion_number_density.index.nlevels == level_number_density.index.nlevels:
        aligned_ion_number_density = ion_number_density.reindex(
            index=level_number_density.index,
            columns=level_number_density.columns,
        )
    else:
        ion_number = level_number_density.index.get_level_values("ion_number")
        if next_higher:
            ion_number = ion_number + 1
        ion_index = pd.MultiIndex.from_arrays(
            [
                level_number_density.index.get_level_values("atomic_number"),
                ion_number,
            ],
            names=ion_number_density.index.names,
        )
        aligned_ion_number_density = ion_number_density.reindex(
            index=ion_index, columns=level_number_density.columns
        )
    return pd.DataFrame(
        aligned_ion_number_density.to_numpy(),
        index=level_number_density.index,
        columns=level_number_density.columns,
    )


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
