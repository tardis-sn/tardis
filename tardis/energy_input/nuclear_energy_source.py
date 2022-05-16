from nuclear.ejecta import Ejecta
from nuclear.io.nndc import get_decay_radiation_database, store_decay_radiation


def decay_nuclides(shell_mass, initial_composition, epoch):
    """Decay model

    Parameters
    ----------
    shell_mass : float
        Mass of the shell in grams
    initial_composition : DataFrame
        Initial ejecta composition
    epoch : float
        Time in days

    Returns
    -------
    DataFrame
        New composition at time epoch
    """
    decay_model = Ejecta(shell_mass, initial_composition)

    new_fractions = decay_model.decay(epoch)
    return new_fractions


def get_decay_database(
    isotope_abundance,
):
    """Gets the decay radiation database for a set
    of isotopes

    Parameters
    ----------
    isotope_abundance : DataFrame
        DataFrame of simulation isotope masses per shell

    Returns
    -------
    DataFrame
        Decay radiation database
    DataFrame
        Metadata for the decay radiation database
    """
    for column in isotope_abundance:
        if column == "Fe56":
            continue
        store_decay_radiation(column, force_update=False)

    decay_rad_db, meta = get_decay_radiation_database()

    return decay_rad_db, meta
