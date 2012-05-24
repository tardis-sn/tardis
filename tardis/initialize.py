#initializing values; reading config files and so on.
import ConfigParser
import constants
import pkgutil
from collections import OrderedDict
import StringIO
import numpy as np

def read_atomic_data(fname=None):
    if fname is None:
        data = np.recfromtxt(StringIO.StringIO(
                pkgutil.get_data('tardis', 'data/atoms.dat')),
                names=('atom', 'symbol', 'mass'))
    else:
        data = recfromtxt(fname,
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
    trad = config.getfloat('general', 'trad')
    t_exp = config.getfloat('general', 't_exp') *  86400.
    vph = config.getfloat('general', 'vph') * 1e5
    dens =  config.getfloat('general', 'dens')
    oxygen_buffer = config.get('general', 'oxygen_buffer').lower() in ('1','true', 'yes')
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
    print named_abundances
    return trad, t_exp, vph, dens, named_abundances
    
def read_symbol2z(fname=None):
    #a lookup dictionary converting between atomic symbols and number
    atom_data = read_atomic_data(fname)
    return OrderedDict([(symbol.lower(), atom) for atom, symbol in zip(atom_data['atom'], atom_data['symbol'])])

def read_z2symbol(fname=None):
    #a lookup dictionary converting between atomic symbols and number
    atom_data = read_atomic_data(fname)
    return OrderedDict([(atom, symbol.lower()) for atom, symbol in zip(atom_data['atom'], atom_data['symbol'])])