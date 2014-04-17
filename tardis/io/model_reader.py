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


def read_density_file(density_filename, density_filetype, time_explosion, v_inner_boundary=0.0, v_outer_boundary=np.inf):
    """
    read different density file formats

    Parameters
    ----------

    density_filename: ~str
        filename or path of the density file

    density_filetype: ~str
        type of the density file

    time_explosion: ~astropy.units.Quantity
        time since explosion used to scale the density

    """
    file_parsers = {'artis': read_artis_density,
                    'simple_ascii': read_simple_ascii_density}

    time_of_model, index, v_inner, v_outer, unscaled_mean_densities = file_parsers[density_filetype](density_filename)
    mean_densities = calculate_density_after_time(unscaled_mean_densities, time_of_model, time_explosion)

    if v_inner_boundary > v_outer_boundary:
        raise ConfigurationError('v_inner_boundary > v_outer_boundary ({0:s} > {1:s}). unphysical!'.format(v_inner_boundary, v_outer_boundary))
    
    if not np.isclose(v_inner_boundary, 0.0) and v_inner_boundary > v_inner[0]:

        if v_inner_boundary > v_outer[-1]:
            raise ConfigurationError('Inner boundary selected outside of model')

        inner_boundary_index = v_inner.searchsorted(v_inner_boundary) - 1

    else:
        inner_boundary_index = None
        v_inner_boundary = v_inner[0]
        logger.warning("v_inner_boundary requested too small for readin file. Boundary shifted to match file.")

    if not np.isinf(v_outer_boundary) and v_outer_boundary < v_outer[-1]:
        outer_boundary_index = v_outer.searchsorted(v_outer_boundary) + 1
    else:
        outer_boundary_index = None
        v_outer_boundary = v_outer[-1]
        logger.warning("v_outer_boundary requested too large for readin file. Boundary shifted to match file.")

        
    v_inner = v_inner[inner_boundary_index:outer_boundary_index]
    v_inner[0] = v_inner_boundary

    v_outer = v_outer[inner_boundary_index:outer_boundary_index]
    v_outer[-1] = v_outer_boundary

    mean_densities = mean_densities[inner_boundary_index:outer_boundary_index]


    return v_inner, v_outer, mean_densities, inner_boundary_index, outer_boundary_index

def read_abundances_file(abundance_filename, abundance_filetype, inner_boundary_index=None, outer_boundary_index=None):
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

    data = recfromtxt(fname, skip_header=1, names=('index', 'velocity', 'density'), dtype=(int, float, float))
    velocity = (data['velocity'] * u.km / u.s).to('cm/s')
    v_inner, v_outer = velocity[:-1], velocity[1:]
    mean_density = (data['density'] * u.Unit('g/cm^3'))[1:]

    return time_of_model, data['index'], v_inner, v_outer, mean_density

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
    mean_density = u.Quantity(10 ** artis_model['mean_densities_0'], 'g/cm^3')
    v_inner, v_outer = velocity[:-1], velocity[1:]

    return time_of_model, artis_model['index'], v_inner, v_outer, mean_density



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





def calculate_density_after_time(densities, time_0, time_explosion):
    """
    scale the density from an initial time of the model to the time of the explosion by ^-3

    Parameters:
    -----------

    densities: ~astropy.units.Quantity
        densities

    time_0: ~astropy.units.Quantity
        time of the model

    time_explosion: ~astropy.units.Quantity
        time to be scaled to

    Returns:
    --------

    scaled_density
    """

    return densities * (time_explosion / time_0) ** -3


