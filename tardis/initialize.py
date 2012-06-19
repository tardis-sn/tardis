#initializing values; reading config files and so on.
import ConfigParser
import constants
import pkgutil
from collections import OrderedDict
import StringIO
import numpy as np
import logging

logger = logging.getLogger(__name__)


def read_atomic_data(fname=None):
    if fname is None:
        data = np.recfromtxt(StringIO.StringIO(
            pkgutil.get_data('tardis', 'data/atoms.dat')),
            names=('atom', 'symbol', 'mass'))
    else:
        data = np.recfromtxt(fname,
            names=('atom', 'symbol', 'mass'))
    return data


def read_ionization_data(fname=None):
    if fname is None:
        data = np.recfromtxt(StringIO.StringIO(
            pkgutil.get_data('tardis', 'data/ionization.dat')),
            names=('atom', 'ion', 'energy'))
    else:
        data = recfromtxt(fname,
            names=('atom', 'ion', 'energy'))
    return data


def read_simple_config(fname):
    config = ConfigParser.ConfigParser()
    config.read(fname)
    #trad = config.getfloat('general', 'trad')
    time_exp = config.getfloat('general', 'time_exp') * 86400.
    v_inner = config.getfloat('general', 'v_inner') * 1e5
    v_outer = config.getfloat('general', 'v_outer') * 1e5
    density = config.getfloat('general', 'density')
    log_l_lsun = config.getfloat('general', 'log_l_lsun')
    calibration_no_of_packets = int(config.getfloat('general', 'calibration_packets'))
    spectrum_no_of_packets = int(config.getfloat('general', 'spectrum_packets'))
    iterations = config.getint('general', 'iterations')
    # packet energies sum up to 1.


    luminosity_outer = 10 ** (log_l_lsun + constants.log_lsun) # in cgs

    r_inner = v_inner * time_exp # in cm
    r_outer = v_outer * time_exp # in cm

    t_outer = (luminosity_outer / (4 * np.pi * constants.sigma_sb * r_inner ** 2)) ** .25

    time_of_simulation = 1 / luminosity_outer

    oxygen_buffer = config.get('general', 'oxygen_buffer').lower() in ('1', 'true', 'yes')
    #using oxygen as buffer without checking
    named_abundances = OrderedDict([(key, float(value)) for key, value in config.items('abundance')])

    if oxygen_buffer:
        abundance_sum = sum([value for key, value in named_abundances.items() if key != 'o'])
        assert abundance_sum <= 1
        named_abundances['o'] = 1. - abundance_sum
    else:
        abundance_sum = sum(named_abundances.values())

        for key in named_abundances:
            named_abundances[key] /= abundance_sum
    logger.info('Included elements: %s', ','.join(named_abundances.keys()))

    return {'time_exp': time_exp,
            'v_inner': v_inner,
            'density': density,
            'luminosity_outer': luminosity_outer,
            'r_inner': r_inner,
            'r_outer': r_outer,
            't_outer': t_outer,
            'time_simulation': time_of_simulation,
            'abundances': named_abundances,
            'iterations': iterations,
            'no_of_calibration_packets': calibration_no_of_packets,
            'no_of_spectrum_packets': spectrum_no_of_packets
    }


def read_normal_config(fname):
    config = ConfigParser.ConfigParser()
    config.read(fname)
    #trad = config.getfloat('general', 'trad')
    time_exp = config.getfloat('general', 'time_exp') * constants.days2seconds
    v_inner = config.getfloat('general', 'v_inner') * 1e5

    log_l_lsun = config.getfloat('general', 'log_l_lsun')
    calibration_no_of_packets = int(config.getfloat('general', 'calibration_packets'))
    spectrum_no_of_packets = int(config.getfloat('general', 'spectrum_packets'))
    iterations = config.getint('general', 'iterations')
    # packet energies sum up to 1.


    luminosity_outer = 10 ** (log_l_lsun + constants.log_lsun) # in cgs

    r_inner = v_inner * time_exp # in cm

    t_outer = (luminosity_outer / (4 * np.pi * constants.sigma_sb * r_inner ** 2)) ** .25

    time_of_simulation = 1 / luminosity_outer

    oxygen_buffer = config.get('general', 'oxygen_buffer').lower() in ('1', 'true', 'yes')
    #using oxygen as buffer without checking
    named_abundances = OrderedDict([(key, float(value)) for key, value in config.items('abundance')])

    if oxygen_buffer:
        abundance_sum = sum([value for key, value in named_abundances.items() if key != 'o'])
        assert abundance_sum <= 1
        named_abundances['o'] = 1. - abundance_sum
    else:
        abundance_sum = sum(named_abundances.values())

        for key in named_abundances:
            named_abundances[key] /= abundance_sum
    logger.info('Included elements: %s', ','.join(named_abundances.keys()))

    return {'time_exp': time_exp,
            'v_inner': v_inner,
            'luminosity_outer': luminosity_outer,
            'r_inner': r_inner,
            't_outer': t_outer,
            'time_simulation': time_of_simulation,
            'abundances': named_abundances,
            'iterations': iterations,
            'no_of_calibration_packets': calibration_no_of_packets,
            'no_of_spectrum_packets': spectrum_no_of_packets
    }


def read_symbol2z(fname=None):
    #a lookup dictionary converting between atomic symbols and number
    atom_data = read_atomic_data(fname)
    return OrderedDict([(symbol.lower(), atom) for atom, symbol in zip(atom_data['atom'], atom_data['symbol'])])


def read_z2symbol(fname=None):
    #a lookup dictionary converting between atomic symbols and number
    atom_data = read_atomic_data(fname)
    return OrderedDict([(atom, symbol.lower()) for atom, symbol in zip(atom_data['atom'], atom_data['symbol'])])