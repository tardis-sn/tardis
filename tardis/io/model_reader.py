#reading different model files
from numpy import recfromtxt, genfromtxt
import pandas as pd
from astropy import units as u

from tardis.io import config_reader
from tardis.util import parse_quantity


def read_simple_ascii_density(fname):
    """
    Reading a density file of the following structure (example; lines starting with a hash will be ignored):
    The first density describes the mean density in the center of the model and is not used.
    5 s
    #index velocity [km/s] density [g/cm^3]
    0 1.1e4 1.6e8
    1 1.2e4 1.7e8

    Parameters
    ----------

    fname: str
        filename or path with filename


    Returns
    -------

    time_of_model: ~astropy.units.Quantity
        time at which the model is valid

    data: ~pandas.DataFrame
        data frame containing index, velocity (in km/s) and density
    """

    with open(fname) as fh:
        time_of_model_string = fh.readline().strip()
        time_of_model = parse_quantity(time_of_model_string)

    data = recfromtxt(fname, skip_header=1, names=('index', 'velocity', 'density'), dtype=(int, float, float))
    velocity = data['velocity'] * u.km / u.s
    v_inner, v_outer = velocity[:-1], velocity[1:]
    mean_density = (data['density'] * u.Unit('g/cm^3'))[1:]

    return time_of_model, data['index'], v_inner, v_outer, mean_density


def read_simple_ascii_abundances(fname):
    """
    Reading an abundance file of the following structure (example; lines starting with hash will be ignored):

    #index element1, element2, ..., element30
    0 0.4 0.3, .. 0.2

    Parameters
    ----------

    fname: str
        filename or path with filename

    Returns
    -------

    data: ~pandas.DataFrame
        data frame containing index, element1 - element30
    """
    data = pd.DataFrame(recfromtxt(fname))
    return data