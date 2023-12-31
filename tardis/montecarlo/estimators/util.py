import pandas as pd


def bound_free_estimator_array2frame(
    bound_free_estimator_array, level2continuum_idx
):
    """
    Transform a bound-free estimator array to a DataFrame.

    This function transforms a bound-free estimator array with entries
    sorted by frequency to a multi-indexed DataFrame sorted by level.

    Parameters
    ----------
    bf_estimator_array : numpy.ndarray, dtype float
        Array of bound-free estimators (e.g., for the stimulated recombination rate)
        with entries sorted by the threshold frequency of the bound-free continuum.
    level2continuum_idx : pandas.Series, dtype int
        Maps a level MultiIndex (atomic_number, ion_number, level_number) to
        the continuum_idx of the corresponding bound-free continuum (which are
        sorted by decreasing frequency).

    Returns
    -------
    pandas.DataFrame, dtype float
        Bound-free estimators indexed by (atomic_number, ion_number, level_number).
    """
    bf_estimator_frame = pd.DataFrame(
        bound_free_estimator_array, index=level2continuum_idx.index
    ).sort_index()
    bf_estimator_frame.columns.name = "Shell No."
    return bf_estimator_frame
