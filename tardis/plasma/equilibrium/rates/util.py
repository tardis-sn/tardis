from collections.abc import Collection

import pandas as pd


def suppress_ground_recombination(
    rates: pd.DataFrame, case_b_species: Collection[tuple[int, int]]
) -> pd.DataFrame:
    """Copy and mask recombination into selected bound-ion ground levels.

    Parameters
    ----------
    rates : pandas.DataFrame
        Level-resolved recombination rates or thermal contributions. The bound
        species is identified by ``atomic_number`` and ``ion_number``; the
        original bound level is ``level_number_source`` in rate frames.
    case_b_species : collection of tuple[int, int]
        Bound species whose level-zero recombination is suppressed.

    Returns
    -------
    pandas.DataFrame
        Masked copy, or the unmodified input when no policy is configured.
    """
    if not case_b_species:
        return rates
    index = rates.index
    species_index = pd.MultiIndex.from_arrays(
        [
            index.get_level_values("atomic_number"),
            index.get_level_values("ion_number"),
        ]
    )
    level_name = (
        "level_number_source"
        if "level_number_source" in index.names
        else "level_number"
    )
    ground_mask = species_index.isin(case_b_species) & (
        index.get_level_values(level_name) == 0
    )
    result = rates.copy(deep=True)
    result.loc[ground_mask] = 0.0
    return result


def reindex_ion_population_to_level_population(
    ion_population: pd.DataFrame,
    level_population: pd.DataFrame,
    next_higher: bool = True,
) -> pd.DataFrame:
    """Align ion population indexes with level populations by ionization stage.

    Parameters
    ----------
    ion_population : pd.DataFrame
        Ion-resolved values indexed by atomic number and ion number.
    level_population : pd.DataFrame
        Level-resolved values indexed by atomic number, ion number, and level
        number.
    next_higher : bool, optional
        Select the continuum ion for each level when ``True``; otherwise select
        the ion containing the level.

    Returns
    -------
    pd.DataFrame
        Ion-resolved values with the level-population index and columns.
    """
    if ion_population.index.nlevels == level_population.index.nlevels:
        aligned_ion_population = ion_population.reindex(
            index=level_population.index,
            columns=level_population.columns,
        )
    else:
        ion_number = level_population.index.get_level_values("ion_number")
        if next_higher:
            ion_number = ion_number + 1
        ion_index = pd.MultiIndex.from_arrays(
            [
                level_population.index.get_level_values("atomic_number"),
                ion_number,
            ],
            names=ion_population.index.names,
        )
        aligned_ion_population = ion_population.reindex(
            index=ion_index, columns=level_population.columns
        )
    return pd.DataFrame(
        aligned_ion_population.to_numpy(),
        index=level_population.index,
        columns=level_population.columns,
    )


def reindex_ionization_rate_dataframe(
    rate_dataframe: pd.DataFrame, recombination: bool = False
) -> pd.DataFrame:
    """Add ion source and destination labels without mutating the input.

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
    rate_dataframe = rate_dataframe.rename_axis(
        [
            "atomic_number",
            "ion_number",
            "level_number_source",
        ]
    )

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
