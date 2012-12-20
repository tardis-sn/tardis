#initializing values; reading config files and so on.
import ConfigParser
import constants
from collections import OrderedDict
import numpy as np
import logging
import atomic
import os
import re



#TODO Fix use pkgutil


default_lucy99_general_fname = os.path.abspath(os.path.join(os.path.dirname(__file__),
    'data',
    'lucy99_general_default.ini'))

logger = logging.getLogger(__name__)

default_lucy99_abundance = OrderedDict()
default_lucy99_abundance['C'] = 0.01
default_lucy99_abundance['O'] = 0.01
default_lucy99_abundance['Ne'] = 0.01
default_lucy99_abundance['Mg'] = 0.01
default_lucy99_abundance['Si'] = 0.45
default_lucy99_abundance['S'] = 0.35
default_lucy99_abundance['Ar'] = 0.04
default_lucy99_abundance['Ca'] = 0.03
default_lucy99_abundance['Fe'] = 0.07
default_lucy99_abundance['Co'] = 0.01
default_lucy99_abundance['Ni'] = 0.01

def write_uniform_abundance_config(fname, default_general_fname=default_lucy99_general_fname,
                                   default_abundances=default_lucy99_abundance):
    general_section = file(default_general_fname).read()
    atomic_data = atomic.read_atomic_data()
    with file(fname, 'w') as fh:
        fh.write(general_section)
        fh.write('\n\n[uniform_abundances]\n')
        fh.write('oxygen_buffer=False\n')
        for line in atomic_data:
            if line['symbol'] not in default_abundances:
                continue
            fh.write('%s=%.2f\n' % (line['symbol'], default_abundances[line['symbol']]))


def read_generic_general_config(config_object, config_dict):
    config_dict['max_atom'] = config_object.getint('general', 'max_atom')
    if config_dict['max_atom'] != 30:
        logger.warn('max_atom is not 30; normal atomic models do not work with this setting (yet)')

    config_dict['max_ion'] = config_object.getint('general', 'max_ion')
    if config_dict['max_ion'] != 30:
        logger.warn('max_ion is not 30; normal atomic models do not work with this setting (yet)')

    config_dict['time_exp'] = config_object.getfloat('general', 'time_exp') * constants.days2seconds
    logger.debug('Time since explosion %.2f s', config_dict['time_exp'])

    if config_object.has_option('uniform_abundances', 'oxygen_buffer'):
        config_dict['oxygen_buffer'] = config_object.getboolean('uniform_abundances', 'oxygen_buffer')

    else:
        logger.warn('oxygen_buffer keyword not found setting oxygen_buffer to True')
        config_dict['oxygen_buffer'] = True

    if config_object.has_option('general', 'line_interaction_type'):
        config_dict['line_interaction_type'] = config_object.get('general', 'line_interaction_type')
        if config_dict['line_interaction_type'] not in ('macro', 'scatter'):
            raise ValueError('line_interaction_type keyword can only be macro or scatter')
    else:
        logger.info('Line interaction type keyword (line_interaction_type) not set. Setting "macro"')
        config_dict['line_interaction_type'] = 'macro'

    if config_object.has_option('general', 'exclude_ions'):
        config_dict['exclude_ions'] = map(int, config_object.get('general', 'exclude_ions').split(','))
    else:
        config_dict['exclude_ions'] = None

    return config_dict


def read_uniform_w7_config(fname, default_general_fname=default_lucy99_general_fname,
                           default_abundances=default_lucy99_abundance):
    """
    Read standard tardis config file
    :param fname:
    :return:
    """

    tardis_config = {'config_type': 'uniform_w7'}
    config = ConfigParser.ConfigParser()

    config.read([default_general_fname, fname])

    tardis_config = read_generic_general_config(config, tardis_config)
    tardis_config['no_of_shells'] = config.getint('general', 'no_of_shells')
    tardis_config['v_inner'] = config.getfloat('general', 'v_inner') * 1e5 # converting from km/s to cm/s
    tardis_config['r_inner'] = tardis_config['v_inner'] * tardis_config['time_exp']
    logger.debug('Inner boundary velocity %.2f cm/s (at %.2f cm)', tardis_config['v_inner'], tardis_config['r_inner'])

    #TODO add luminosity check to get something else than bolometric
    log_l_lsun = config.getfloat('general', 'log_l_lsun')

    tardis_config['luminosity_outer'] = 10 ** (log_l_lsun + constants.log_lsun)
    tardis_config['t_outer'] = (tardis_config['luminosity_outer'] / (
        4 * np.pi * constants.sigma_sb * tardis_config['r_inner'] ** 2)) ** 0.25

    tardis_config['time_of_simulation'] = 1 / tardis_config['luminosity_outer']

    logger.debug('Required output luminosity is %s ergs/s => outer temperature %.2f K; duration of simulation %s s',
        tardis_config['luminosity_outer'],
        tardis_config['t_outer'],
        tardis_config['time_of_simulation'])
    tardis_config['t_rads'] = np.ones(tardis_config['no_of_shells']) * tardis_config['t_outer']
    tardis_config['no_of_calibration_packets'] = int(config.getfloat('general', 'calibration_packets'))
    tardis_config['no_of_spectrum_packets'] = int(config.getfloat('general', 'spectrum_packets'))
    tardis_config['iterations'] = config.getint('general', 'iterations')

    if 'abund' not in ''.join(config.sections()):
        logger.warn('No abundance section specified. Using defaults')
        named_abundances = default_abundances

    elif config.has_section('uniform_abundances'):
        named_abundances = {}

        for atom, abundance in config.items('uniform_abundances'):
            if atom == 'oxygen_buffer':
                continue
            named_abundances[atom.lower()] = config.getfloat('uniform_abundances', atom)
    else:
        raise ValueError('only uniform abundances supported at the moment. ')
        #TODO add reading of different abundance specifications

    abundances = named2array_abundances(named_abundances, tardis_config['max_atom'],
        oxygen_buffer=tardis_config['oxygen_buffer'])
    logger.info('Included elements: %s', ','.join(named_abundances.keys()))
    tardis_config['abundances'] = abundances

    return tardis_config


def read_shell_config(fname, default_general_fname=default_lucy99_general_fname):
    tardis_config = {'config_type': 'shell'}
    config = ConfigParser.ConfigParser()

    config.read([default_general_fname, fname])

    tardis_config = read_generic_general_config(config, tardis_config)

    #TODO add luminosity check to get something else than bolometric


    tardis_config['no_of_calibration_packets'] = int(config.getfloat('general', 'calibration_packets'))
    tardis_config['no_of_spectrum_packets'] = int(config.getfloat('general', 'spectrum_packets'))
    tardis_config['iterations'] = config.getint('general', 'iterations')

    shell_numbers = []
    shell_names = []
    for section in config.sections():
        if section == 'general': continue
        shell_match = re.match('shell(\d+)', section, flags=re.IGNORECASE)
        if shell_match is None:
            raise ValueError('Only the "general" section and shellxx'
                             '(where xx is the shell number) are allowed')
        else:
            shell_numbers.append(int(shell_match.groups()[0]))
            shell_names.append(section)

    symbol2z = read_symbol2z()
    atom_symbols = [item.lower() for item in symbol2z]


    #TODO add check and warning to see if all shells are there (not only shell 1 and 4 and 7).
    densities = np.zeros(len(shell_numbers) + 1)
    velocities = np.empty(len(shell_numbers) + 1)
    abundances = np.zeros((len(shell_numbers), tardis_config['max_atom']))
    ws = -np.ones(len(shell_numbers))
    t_rads = -np.ones(len(shell_numbers))
    for i, idx in enumerate(np.argsort(shell_numbers)):
        current_abundances_dict = {}
        for keyword, value in config.items(shell_names[idx]):
            if keyword.lower() == 'density':
                densities[i + 1] = value
            elif keyword.lower() == 'v_inner':
                velocities[i] = float(value) * 1e5
            elif keyword.lower() == 'v_outer':
                if idx == np.argmax(shell_numbers):
                    velocities[i + 1] = float(value) * 1e5
                else:
                    raise ValueError(
                        "Can't specify v_outer in any of the inner layers. Only allowed on the outermost shell")
            elif keyword.lower() == 'w':
                ws[i] = value
            elif keyword.lower() == 't_rad':
                t_rads[i] = value

            elif keyword.lower() in atom_symbols:
                current_abundances_dict[keyword] = value
            else:
                raise ValueError('Keyword %s not allowed. Allowed Keywords: atomic symbols, v_inner, v_outer & density')



        #TODO implement oxygen buffer for this shell config reader (change them all to be in the general config)
        current_abundances = named2array_abundances(current_abundances_dict, tardis_config['max_atom'],
            oxygen_buffer=tardis_config['oxygen_buffer'])
        if sum(current_abundances) == 0:
            raise ValueError('Either no abundances specified or abundances add up to 0.')
        abundances[i] = current_abundances
    tardis_config['log_l_lsun'] = config.getfloat('general', 'log_l_lsun')
    tardis_config['r_inner'] = velocities[0] * tardis_config['time_exp']
    tardis_config['luminosity_outer'] = 10 ** (tardis_config['log_l_lsun'] + constants.log_lsun)
    tardis_config['t_outer'] = (tardis_config['luminosity_outer'] /
                                (4 * np.pi * constants.sigma_sb * tardis_config['r_inner'] ** 2)) ** .25
    tardis_config['time_of_simulation'] = 1 / tardis_config['luminosity_outer']
    t_rads = np.array([tardis_config['t_outer'] if t_rad < 0 else t_rad for t_rad in t_rads])
    logger.debug('Required output luminosity is %s ergs/s => outer temperature %.2f K; duration of simulation %s s',
        tardis_config['luminosity_outer'],
        tardis_config['t_outer'],
        tardis_config['time_of_simulation'])

    tardis_config['abundances'] = abundances
    tardis_config['velocities'] = velocities
    tardis_config['densities'] = densities
    tardis_config['t_rads'] = t_rads
    tardis_config['ws'] = ws
    print tardis_config
    return tardis_config


def named2array_abundances(named_abundances, max_atom, oxygen_buffer=True):
    symbol2z = read_symbol2z()
    atom_symbols = [item.lower() for item in symbol2z]
    #converting the abundances dictionary dict to an array
    abundances = np.zeros(max_atom)

    for symbol in named_abundances:
        if symbol.lower() == 'o' and oxygen_buffer:
            print symbol, named_abundances
            logger.warn('Requested oxygen as buffer but specifying an abundance (O=%.2f)' % named_abundances[symbol])
            continue
        abundances[symbol2z[symbol.lower()] - 1] = named_abundances[symbol]

    if oxygen_buffer:
        if max_atom < 8:
            raise ValueError('Can\'t balance elements with oxygen if max_atom less than 8 (max_atom=%d)' % max_atom)
        abundance_sum = np.sum(abundances)
        if abundance_sum > 1:
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

