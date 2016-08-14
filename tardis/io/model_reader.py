#reading different model files

import numpy as np
from numpy import recfromtxt, genfromtxt
import pandas as pd
from astropy import units as u

import logging
# Adding logging support
logger = logging.getLogger(__name__)

from tardis.util import parse_quantity


class ConfigurationError(Exception):
    pass


def read_density_file(filename, filetype):
    """
    read different density file formats

    Parameters
    ----------

    filename: ~str
        filename or path of the density file

    filetype: ~str
        type of the density file

    Returns
    -------
    time_of_model: ~astropy.units.Quantity
        time at which the model is valid

    velocity: ~np.ndarray
        the array containing the velocities

    unscaled_mean_densities: ~np.ndarray
        the array containing the densities

    """
    file_parsers = {'artis': read_artis_density,
                    'simple_ascii': read_simple_ascii_density}

    time_of_model,velocity, unscaled_mean_densities = file_parsers[filetype](filename)
    return time_of_model, velocity, unscaled_mean_densities

def read_abundances_file(abundance_filename, abundance_filetype,
                         inner_boundary_index=None, outer_boundary_index=None):
    """
    read different density file formats

    Parameters
    ----------

    abundance_filename: ~str
        filename or path of the density file

    abundance_filetype: ~str
        type of the density file

    inner_boundary_index: int
        index of the inner shell, default None

    outer_boundary_index: int
        index of the outer shell, default None


    """

    file_parsers = {'simple_ascii': read_simple_ascii_abundances,
                    'artis': read_simple_ascii_abundances}

    index, abundances = file_parsers[abundance_filetype](abundance_filename)
    if outer_boundary_index is not None:
        outer_boundary_index_m1 = outer_boundary_index - 1
    else:
        outer_boundary_index_m1 = None
    index = index[inner_boundary_index:outer_boundary_index]
    abundances = abundances.ix[:, slice(inner_boundary_index, outer_boundary_index_m1)]
    abundances.columns = np.arange(len(abundances.columns))
    return index, abundances


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

    data = recfromtxt(fname, skip_header=1,
                      names=('velocity', 'density'), dtype=(float, float))
    velocity = (data['velocity'] * u.km / u.s).to('cm/s')
    mean_density = (data['density'] * u.Unit('g/cm^3'))[1:]

    return time_of_model, velocity, mean_density

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
        for i, line in enumerate(file(fname)):
            if i == 0:
                no_of_shells = np.int64(line.strip())
            elif i == 1:
                time_of_model = u.Quantity(float(line.strip()), 'day').to('s')
            elif i == 2:
                break

    artis_model_columns = ['index', 'velocities', 'mean_densities_0', 'ni56_fraction', 'co56_fraction', 'fe52_fraction',
                           'cr48_fraction']
    artis_model = recfromtxt(fname, skip_header=2, usecols=(0, 1, 2, 4, 5, 6, 7), unpack=True,
                                dtype=[(item, np.float64) for item in artis_model_columns])


    velocity = u.Quantity(artis_model['velocities'], 'km/s').to('cm/s')
    mean_density = u.Quantity(10 ** artis_model['mean_densities_0'], 'g/cm^3')[1:]

    return time_of_model, velocity, mean_density



def read_simple_ascii_abundances(fname):
    """
    Reading an abundance file of the following structure (example; lines starting with hash will be ignored):
    The first line of abundances describe the abundances in the center of the model and are not used.
    #index element1, element2, ..., element30
    0 0.4 0.3, .. 0.2

    Parameters
    ----------

    fname: str
        filename or path with filename

    Returns
    -------

    index: ~np.ndarray
        containing the indices

    abundances: ~pandas.DataFrame
        data frame containing index, element1 - element30 and columns according to the shells
    """
    data = np.loadtxt(fname)

    index = data[1:,0].astype(int)
    abundances = pd.DataFrame(data[1:,1:].transpose(), index=np.arange(1, data.shape[1]))

    return index, abundances
