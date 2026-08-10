import pandas as pd


def reindex_ionization_rate_dataframe(
    rate_dataframe: pd.DataFrame, recombination: bool = False
) -> pd.DataFrame:
    """Add source and destination states to an ionization-rate table.

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

    # Legacy ion-population solvers represent the destination ion by its
    # total-population slot for both directions.
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


def reindex_level_resolved_recombination_rate_dataframe(
    rate_dataframe: pd.DataFrame,
) -> pd.DataFrame:
    """Index recombination into the originating bound-level destination."""
    indexed_rate = reindex_ionization_rate_dataframe(
        rate_dataframe, recombination=True
    ).reset_index()
    indexed_rate["level_number_destination"] = indexed_rate[
        "level_number_source"
    ]
    indexed_rate["level_number_source"] = 0
    return indexed_rate.set_index(
        [
            "atomic_number",
            "ion_number",
            "ion_number_source",
            "ion_number_destination",
            "level_number_source",
            "level_number_destination",
        ]
    )
