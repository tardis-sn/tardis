#initializing values; reading config files and so on.
import ConfigParser
import constants
from collections import OrderedDict
import numpy as np
import logging
import atomic
import os

#TODO Fix use pkgutil


default_lucy99_general_fname = os.path.abspath(os.path.join(os.path.dirname(__file__),
    'data',
    'lucy99_general_default.ini'))

logger = logging.getLogger(__name__)

default_lucy99_abundance = dict(C=0.01,
    O=0.01,
    Ne=0.01,
    Mg=0.01,
    Si=0.45,
    S=0.35,
    Ar=0.04,
    Ca=0.03,
    Fe=0.07,
    Co=0.01,
    Ni=0.01)


def read_simple_tardis_config(fname, default_general_fname=default_lucy99_general_fname,
                              default_abundances=default_lucy99_abundance):
    """
    Read standard tardis config file
    :param fname:
    :return:
    """

    tardis_config = {}
    config = ConfigParser.ConfigParser()

    config.read([default_general_fname, fname])

    tardis_config['max_atom'] = config.getint('general', 'max_atom')
    if tardis_config['max_atom'] != 30:
        logger.warn('max_atom is not 30; normal atomic models do not work with this setting')

    tardis_config['max_ion'] = config.getint('general', 'max_ion')
    if tardis_config['max_ion'] != 30:
        logger.warn('max_ion is not 30; normal atomic models do not work with this setting')

    tardis_config['time_exp'] = config.getfloat('general', 'time_exp') * constants.days2seconds
    logger.debug('Time since explosion %.2f s', tardis_config['time_exp'])

    tardis_config['v_inner'] = config.getfloat('general', 'v_inner') * 1e5 # converting from km/s to cm/s
    tardis_config['r_inner'] = tardis_config['v_inner'] * tardis_config['time_exp']
    logger.debug('Boundary velocity %.2f cm/s (at %.2f cm)', tardis_config['v_inner'], tardis_config['r_inner'])

    #TODO add luminosity check to get something else than bolometric
    log_l_lsun = config.getfloat('general', 'log_l_lsun')

    tardis_config['luminosity_outer'] = 10 ** (log_l_lsun + constants.log_lsun)
    tardis_config['t_outer'] = (tardis_config['luminosity_outer'] / (
    4 * np.pi * constants.sigma_sb * tardis_config['r_inner'] ** 2)) ** .25
    tardis_config['time_of_simulation'] = 1 / tardis_config['luminosity_outer']

    logger.debug('Required output luminosity is %s ergs/s => outer temperature %.2f K; duration of simulation %s s',
        tardis_config['luminosity_outer'],
        tardis_config['t_outer'],
        tardis_config['time_of_simulation'])

    tardis_config['no_of_calibration_packets'] = int(config.getfloat('general', 'calibration_packets'))
    tardis_config['no_of_spectrum_packets'] = int(config.getfloat('general', 'spectrum_packets'))
    tardis_config['iterations'] = config.getint('general', 'iterations')

    if 'abund' not in ''.join(config.sections()):
        logger.warn('No abundance section specified. Using defaults')
        named_abundances = default_abundances
        oxygen_buffer = False
    elif config.has_section('uniform_abundances'):
        named_abundances = {}
        if config.has_option('uniform_abundances', 'oxygen_buffer'):
            oxygen_buffer = config.getboolean('uniform_abundances', 'oxygen_buffer')
        else:
            logger.warn('oxygen_buffer keyword not found setting oxygen_buffer to True')
            oxygen_buffer = True

        for atom, abundance in config.items('uniform_abundances'):
            named_abundances[atom.lower()] = abundance
    else:
        raise ValueError('only uniform abundances supported at the moment. ')
        #TODO add reading of different abundance specifications



    abundances = named2array_abundances(named_abundances, tardis_config['max_atom'], oxygen_buffer=oxygen_buffer)
    logger.info('Included elements: %s', ','.join(named_abundances.keys()))
    tardis_config['abundances'] = abundances

    return tardis_config


def named2array_abundances(named_abundances, max_atom, oxygen_buffer=True):
    symbol2z = read_symbol2z()

    #converting the abundances dictionary dict to an array
    abundances = np.zeros(max_atom)

    for symbol in named_abundances:
        if symbol.lower() == 'o' and oxygen_buffer:
            logger.warn('Requested oxygen as buffer but specifying an abundance (O=%.2f)' % named_abundances[symbol])
            continue
        abundances[symbol2z[symbol.lower()] - 1] = named_abundances[symbol]

    if oxygen_buffer:
        if max_atom < 8:
            raise ValueError('Can\'t balance elements with oxygen if max_atom less than 8 (max_atom=%d)' % max_atom)
        abundance_sum = np.sum(abundances)
        if abundance_sum >= 1:
            raise ValueError(
                'Requesting oxygen as buffer but specifying relative abundances adding up to more than 1. (sum=%.2f)' % abundance_sum)
        oxygen_abundance = 1 - abundance_sum
        abundances[7] = oxygen_abundance
    return abundances


def read_symbol2z(fname=None):
    #a lookup dictionary converting between atomic symbols and number
    atom_data = atomic.read_atomic_data(fname)
    return OrderedDict([(symbol.lower(), atom) for atom, symbol in zip(atom_data['atom'], atom_data['symbol'])])


def read_z2symbol(fname=None):
    #a lookup dictionary converting between atomic symbols and number
    atom_data = atomic.read_atomic_data(fname)
    return OrderedDict([(atom, symbol.lower()) for atom, symbol in zip(atom_data['atom'], atom_data['symbol'])])