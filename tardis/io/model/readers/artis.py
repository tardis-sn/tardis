import numpy as np
from astropy import units as u
from numpy import recfromtxt


def read_artis_density(fname):
    """
    Reading a density file of the following structure (example; lines starting with a hash will be ignored):
    The first density describes the mean density in the center of the model and is not used.
    5
    #index velocity [km/s] log10(density) [log10(g/cm^3)]
    0 1.1e4 1.6e8
    1 1.2e4 1.7e8

    Parameters
    ----------
    fname : str
        filename or path with filename

    Returns
    -------
    time_of_model : astropy.units.Quantity
        time at which the model is valid
    data : pandas.DataFrame
        data frame containing index, velocity (in km/s) and density
    """

    with open(fname) as fh:
        for i, line in enumerate(open(fname)):
            if i == 0:
                no_of_shells = np.int64(line.strip())
            elif i == 1:
                time_of_model = u.Quantity(float(line.strip()), "day").to("s")
            elif i == 2:
                break

    artis_model_columns = [
        "index",
        "velocities",
        "mean_densities_0",
        "ni56_fraction",
        "co56_fraction",
        "fe52_fraction",
        "cr48_fraction",
    ]
    artis_model = recfromtxt(
        fname,
        skip_header=2,
        usecols=(0, 1, 2, 4, 5, 6, 7),
        unpack=True,
        dtype=[(item, np.float64) for item in artis_model_columns],
    )

    velocity = u.Quantity(artis_model["velocities"], "km/s").to("cm/s")
    mean_density = u.Quantity(10 ** artis_model["mean_densities_0"], "g/cm^3")[
        1:
    ]

    return time_of_model, velocity, mean_density
