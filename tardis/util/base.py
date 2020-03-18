import logging
import os
import re
from collections import OrderedDict

import numexpr as ne
import numpy as np
import pandas as pd
import yaml
from tardis import constants
from astropy import units as u
from pyne import nucname

import tardis
from tardis.io.util import get_internal_data_path

k_B_cgs = constants.k_B.cgs.value
c_cgs = constants.c.cgs.value
h_cgs = constants.h.cgs.value
m_e_cgs = constants.m_e.cgs.value
e_charge_gauss = constants.e.gauss.value

logger = logging.getLogger(__name__)
tardis_dir = os.path.realpath(tardis.__path__[0])

ATOMIC_SYMBOLS_DATA = pd.read_csv(get_internal_data_path('atomic_symbols.dat'), delim_whitespace=True,
                                  names=['atomic_number', 'symbol']).set_index('atomic_number').squeeze()

ATOMIC_NUMBER2SYMBOL = OrderedDict(ATOMIC_SYMBOLS_DATA.to_dict())
SYMBOL2ATOMIC_NUMBER = OrderedDict((y, x) for x, y in ATOMIC_NUMBER2SYMBOL.items())

synpp_default_yaml_fname = get_internal_data_path('synpp_default.yaml')


NUMERAL_MAP = tuple(zip(
    (1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1),
    ('M', 'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I')
))

class MalformedError(Exception):
    pass


class MalformedSpeciesError(MalformedError):

    def __init__(self, malformed_element_symbol):
        self.malformed_element_symbol = malformed_element_symbol

    def __str__(self):
        return ('Expecting a species notation (e.g. "Si 2", "Si II", "Fe IV") '
                '- supplied {0}'.format(self.malformed_element_symbol))


class MalformedElementSymbolError(MalformedError):

    def __init__(self, malformed_element_symbol):
        self.malformed_element_symbol = malformed_element_symbol

    def __str__(self):
        return ('Expecting an atomic symbol (e.g. Fe) - supplied {0}').format(
            self.malformed_element_symbol)


class MalformedQuantityError(MalformedError):

    def __init__(self, malformed_quantity_string):
        self.malformed_quantity_string = malformed_quantity_string

    def __str__(self):
        return ('Expecting a quantity string(e.g. "5 km/s") for keyword '
                '- supplied {0}').format(self.malformed_quantity_string)


def int_to_roman(i):
    """
    Convert an integer into its roman numeral representation.

    Parameters
    ----------
    i : int
        Integer to be converted into roman numerals

    Returns
    -------
    str
        Returns roman numeral representation of i in str format.
    """
    result = []
    for integer, numeral in NUMERAL_MAP:
        count = i // integer
        result.append(numeral * count)
        i -= integer * count
    return ''.join(result)

def roman_to_int(roman_string):
    """
    Convert a roman numeral into its corresponding integer.

    Parameters
    ----------
    roman_string : str
        Roman numeral to be converted into an integer

    Returns
    -------
    int
        Returns integer representation of roman_string
    """

    NUMERALS_SET = set(list(zip(*NUMERAL_MAP))[1])
    roman_string = roman_string.upper()
    if len(set(list(roman_string.upper())) - NUMERALS_SET) != 0:
        raise ValueError('{0} does not seem to be a roman numeral'.format(
            roman_string))
    i = result = 0
    for integer, numeral in NUMERAL_MAP:
        while roman_string[i:i + len(numeral)] == numeral:
            result += integer
            i += len(numeral)
    if result < 1:
        raise ValueError('Can not interpret Roman Numeral {0}'.format(roman_string))
    return result


def calculate_luminosity(
        spec_fname, distance, wavelength_column=0,
        wavelength_unit=u.angstrom, flux_column=1,
        flux_unit=u.Unit('erg / (Angstrom cm2 s)')):
    """
    Calculates luminosity of star.

    Parameters
    ----------
    spec_fname : file or str
        File or file name to be read
    distance : float
        Distance to star
    wavelength_column : int, optional(default = 0)
        Column index in which the wavelength is stored
    wavelength_unit : float, optional(default = u.angstrom)
        Dictates units used for calculating wavelength.
    flux_column : int, optional(default = 1)
        Column index in which the flux is stored
    flux_unit : str, optional(default = u.Unit('erg / (Angstrom cm2 s)')
        Dictates units used for flux

    Returns
    -------
    luminosity.value : float
        Returned luminosity value of star.
    wavelength.min() : float
        Minimum value of wavelength of light
    wavelength.max() : float
        Maximum value of wavelength of light
    """
    #BAD STYLE change to parse quantity
    distance = u.Unit(distance)

    wavelength, flux = np.loadtxt(spec_fname, usecols=(wavelength_column, flux_column), unpack=True)

    flux_density = np.trapz(flux, wavelength) * (flux_unit * wavelength_unit)
    luminosity = (flux_density * 4 * np.pi * distance**2).to('erg/s')

    return luminosity.value, wavelength.min(), wavelength.max()


def create_synpp_yaml(radial1d_mdl, fname, shell_no=0, lines_db=None):
    """
    Create a yaml file that is readable from syn++

    Parameters
    ----------
    radial1d_mdl : Radial1DModel
        Inputted object that will be read into YAML file
    fname : str
        File name for the synpp yaml
    shell_no : int, optional(default = 0)
        Number of shells
    lines_db : file, optional(default = None)

    Raises
    ------
    ValueError
        If the current dataset does not contain necessary reference files
    """

    logger.warning('Currently only works with Si and a special setup')
    if radial1d_mdl.atom_data.synpp_refs is not None:
        raise ValueError(
            'The current atom dataset does not contain the '
            'necessary reference files (please contact the authors)')

    radial1d_mdl.atom_data.synpp_refs['ref_log_tau'] = -99.0
    for key, value in radial1d_mdl.atom_data.synpp_refs.iterrows():
        try:
            radial1d_mdl.atom_data.synpp_refs['ref_log_tau'].loc[key] = np.log10(
                radial1d_mdl.plasma.tau_sobolevs[0].loc[value['line_id']])
        except KeyError:
            pass


    relevant_synpp_refs = radial1d_mdl.atom_data.synpp_refs[
        radial1d_mdl.atom_data.synpp_refs['ref_log_tau'] > -50]

    with open(synpp_default_yaml_fname) as stream:
        yaml_reference = yaml.load(stream, Loader=yaml.CLoader)

    if lines_db is not None:
        yaml_reference['opacity']['line_dir'] = os.path.join(lines_db, 'lines')
        yaml_reference['opacity']['line_dir'] = os.path.join(lines_db, 'refs.dat')

    yaml_reference['output']['min_wl'] = float(
        radial1d_mdl.runner.spectrum.wavelength.to('angstrom').value.min())
    yaml_reference['output']['max_wl'] = float(
        radial1d_mdl.runner.spectrum.wavelength.to('angstrom').value.max())


    #raise Exception("there's a problem here with units what units does synpp expect?")
    yaml_reference['opacity']['v_ref'] = float(
        (radial1d_mdl.tardis_config.structure.v_inner[0].to('km/s') /
         (1000. * u.km / u.s)).value)
    yaml_reference['grid']['v_outer_max'] = float(
        (radial1d_mdl.tardis_config.structure.v_outer[-1].to('km/s') /
         (1000. * u.km / u.s)).value)

    #pdb.set_trace()

    yaml_setup = yaml_reference['setups'][0]
    yaml_setup['ions'] = []
    yaml_setup['log_tau'] = []
    yaml_setup['active'] = []
    yaml_setup['temp'] = []
    yaml_setup['v_min'] = []
    yaml_setup['v_max'] = []
    yaml_setup['aux'] = []

    for species, synpp_ref in relevant_synpp_refs.iterrows():
        yaml_setup['ions'].append(100 * species[0] + species[1])
        yaml_setup['log_tau'].append(float(synpp_ref['ref_log_tau']))
        yaml_setup['active'].append(True)
        yaml_setup['temp'].append(yaml_setup['t_phot'])
        yaml_setup['v_min'].append(yaml_reference['opacity']['v_ref'])
        yaml_setup['v_max'].append(yaml_reference['grid']['v_outer_max'])
        yaml_setup['aux'].append(1e200)

    with open(fname, 'w') as f:
        yaml.dump(yaml_reference, stream=f, explicit_start=True)


def intensity_black_body(nu, T):
    """
    Calculate the intensity of a black-body according to the following formula

    .. math::
        I(\\nu, T) = \\frac{2h\\nu^3}{c^2}\frac{1}
        {e^{h\\nu \\beta_\\textrm{rad}} - 1}

    Parameters
    ----------
    nu : float
        Frequency of light
    T : float
        Temperature in kelvin

    Returns
    -------
    Intensity : float
        Returns the intensity of the black body
    """
    beta_rad = 1 / (k_B_cgs * T)
    coefficient = 2 * h_cgs / c_cgs ** 2
    intensity = ne.evaluate('coefficient * nu**3 / '
                            '(exp(h_cgs * nu * beta_rad) -1 )')
    return intensity


def species_tuple_to_string(species_tuple, roman_numerals=True):
    """
    Convert a species tuple to its corresponding string representation.

    Parameters
    ----------
    species_tuple : tuple
        Tuple of 2 values indicated atomic number and number of
        electrons missing

    roman_numerals : bool, optional(default = TRUE)
        Indicates whether the returned ion number is in roman numerals

    Returns
    -------
    element_symbol, roman_ion_number : str
        Returns corresponding string representation of given tuple
    """
    atomic_number, ion_number = species_tuple
    element_symbol = ATOMIC_NUMBER2SYMBOL[atomic_number]
    if roman_numerals:
        roman_ion_number = int_to_roman(ion_number+1)
        return '{0} {1}'.format(str(element_symbol), roman_ion_number)
    else:
        return '{0} {1:d}'.format(element_symbol, ion_number)


def species_string_to_tuple(species_string):
    """
    Convert a species string to its corresponding tuple representation

    Parameters
    ----------
    species_string : str
        String containing species symbol (e.g. Si II, Fe III)

    Returns
    -------
    atomic_number, ion_number : tuple
        Returns tuple of length 2 indicating atomic number and ion number

    Raises
    ------
    MalformedSpeciesError
        If the inputted string does not match the species format
    """

    try:
        element_symbol, ion_number_string = re.match(r'^(\w+)\s*(\d+)',
                                                     species_string).groups()
    except AttributeError:
        try:
            element_symbol, ion_number_string = species_string.split()
        except ValueError:
            raise MalformedSpeciesError(
                'Species string "{0}" is not of format <element_symbol><number>'
                ' (e.g. Fe 2, Fe2, ..)'.format(species_string))

    atomic_number = element_symbol2atomic_number(element_symbol)

    try:
        ion_number = roman_to_int(ion_number_string)
    except ValueError:
        try:
            ion_number = int(ion_number_string)
        except ValueError:
            raise MalformedSpeciesError(
                "Given ion number ('{}') could not be parsed".format(
                    ion_number_string))

    if ion_number > atomic_number:
        raise ValueError(
            'Species given does not exist: ion number > atomic number')

    return atomic_number, ion_number - 1


def parse_quantity(quantity_string):
    """
    Changes a string into it's corresponding astropy.Quantity object.

    Parameters
    ----------
    quantity_string : str
        String to be converted into astropy.Quantity

    Returns
    -------
    q : ~u.Quantity
        Corresponding astropy.Quantity object for passed string

    Raises
    ------
    MalformedQuantityError
        If string is not properly formatted for Astropy Quantity
    """

    if not isinstance(quantity_string, str):
        raise MalformedQuantityError(quantity_string)

    try:
        value_string, unit_string = quantity_string.split()
    except ValueError:
        raise MalformedQuantityError(quantity_string)

    try:
        value = float(value_string)
    except ValueError:
        raise MalformedQuantityError(quantity_string)

    try:
        q = u.Quantity(value, unit_string)
    except ValueError:
        raise MalformedQuantityError(quantity_string)

    return q


def element_symbol2atomic_number(element_string):
    """
    Takes an element symbol and returns its corresponding atomic number

    Parameters
    ----------
    element_string : str
        Inputted element symbol

    Returns
    -------
    int
        Returned atomic number
    """
    reformatted_element_string = reformat_element_symbol(element_string)
    if reformatted_element_string not in SYMBOL2ATOMIC_NUMBER:
        raise MalformedElementSymbolError(element_string)
    return SYMBOL2ATOMIC_NUMBER[reformatted_element_string]


def atomic_number2element_symbol(atomic_number):
    """
    Convert atomic number to string

    Parameters
    ----------
    atomic_number : int
        Inputted atomic number

    Returns
    -------
    str
       Returned corresponding element symbol
    """
    return ATOMIC_NUMBER2SYMBOL[atomic_number]


def reformat_element_symbol(element_string):
    """
    Reformat the string so the first letter is uppercase and all subsequent
     letters lowercase.

    Parameters
    ----------
    element_string : str
        Inputted element symbol

    Returns
    -------
    str
        Returned reformatted element symbol
    """

    return element_string[0].upper() + element_string[1:].lower()


def quantity_linspace(start, stop, num, **kwargs):
    """
    Essentially the same input parameters as linspace, but
    calculated for an astropy quantity start and stop.

    Parameters
    ----------
    start : ~astropy.Quantity
        Starting value of the sequence
    stop : ~astropy.Quantity
        End value of the sequence
    num : int
        Number of samples to generate

    Returns
    -------
    ~astropy.Quantity
        Returns num evenly spaced characters of type astropy.Quantity

    Raises
    ------
    ValueError
        If start and stop values have no unit attribute.
    """
    if not (hasattr(start, 'unit') and hasattr(stop, 'unit')):
        raise ValueError('Both start and stop need to be quantities with a '
                         'unit attribute')

    return (np.linspace(start.value, stop.to(start.unit).value, num, **kwargs)
            * start.unit)


def convert_abundances_format(fname, delimiter=r'\s+'):
    """
    Changes format of file containing abundances into data frame

    Parameters
    ----------
    fname : file, str
        File or file name that contains abundance info
    delimiter : str, optional(default = '\\s+')
        Determines the separator for splitting file

    Returns
    -------
    DataFrame
        Corresponding data frame
    """
    df = pd.read_csv(fname, delimiter=delimiter, comment='#', header=None)
    # Drop shell index column
    df.drop(df.columns[0], axis=1, inplace=True)
    # Assign header row
    df.columns = [nucname.name(i)
                  for i in range(1, df.shape[1] + 1)]
    return df