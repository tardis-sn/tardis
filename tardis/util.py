# Utilities for TARDIS

from astropy import units as u, constants
import numpy as np
import os
import yaml

import logging

k_B_cgs = constants.k_B.cgs.value
c_cgs = constants.c.cgs.value
h_cgs = constants.h.cgs.value
m_e_cgs = constants.m_e.cgs.value
e_charge_gauss = constants.e.gauss.value


logger = logging.getLogger(__name__)

synpp_default_yaml_fname = os.path.join(os.path.dirname(__file__), 'data', 'synpp_default.yaml')

def int_to_roman(input):
   """
   from http://code.activestate.com/recipes/81611-roman-numerals/
   Convert an integer to Roman numerals.

   Examples:
   >>> int_to_roman(0)
   Traceback (most recent call last):
   ValueError: Argument must be between 1 and 3999

   >>> int_to_roman(-1)
   Traceback (most recent call last):
   ValueError: Argument must be between 1 and 3999

   >>> int_to_roman(1.5)
   Traceback (most recent call last):
   TypeError: expected integer, got <type 'float'>

   >>> for i in range(1, 21): print int_to_roman(i)
   ...
   I
   II
   III
   IV
   V
   VI
   VII
   VIII
   IX
   X
   XI
   XII
   XIII
   XIV
   XV
   XVI
   XVII
   XVIII
   XIX
   XX
   >>> print int_to_roman(2000)
   MM
   >>> print int_to_roman(1999)
   MCMXCIX
   """
   input = int(input)
   if type(input) != type(1):
      raise TypeError, "expected integer, got %s" % type(input)
   if not 0 < input < 4000:
      raise ValueError, "Argument must be between 1 and 3999"
   ints = (1000, 900,  500, 400, 100,  90, 50,  40, 10,  9,   5,  4,   1)
   nums = ('M',  'CM', 'D', 'CD','C', 'XC','L','XL','X','IX','V','IV','I')
   result = ""
   for i in range(len(ints)):
      count = int(input / ints[i])
      result += nums[i] * count
      input -= ints[i] * count
   return result

def roman_to_int(input):
   """
   from http://code.activestate.com/recipes/81611-roman-numerals/
   Convert a roman numeral to an integer.

   >>> r = range(1, 4000)
   >>> nums = [int_to_roman(i) for i in r]
   >>> ints = [roman_to_int(n) for n in nums]
   >>> print r == ints
   1

   >>> roman_to_int('VVVIV')
   Traceback (most recent call last):
    ...
   ValueError: input is not a valid roman numeral: VVVIV
   >>> roman_to_int(1)
   Traceback (most recent call last):
    ...
   TypeError: expected string, got <type 'int'>
   >>> roman_to_int('a')
   Traceback (most recent call last):
    ...
   ValueError: input is not a valid roman numeral: A
   >>> roman_to_int('IL')
   Traceback (most recent call last):
    ...
   ValueError: input is not a valid roman numeral: IL
   """
   if type(input) != type(""):
      raise TypeError, "expected string, got %s" % type(input)
   input = input.upper()
   nums = ['M', 'D', 'C', 'L', 'X', 'V', 'I']
   ints = [1000, 500, 100, 50,  10,  5,   1]
   places = []
   for c in input:
      if not c in nums:
         raise ValueError, "input is not a valid roman numeral: %s" % input
   for i in range(len(input)):
      c = input[i]
      value = ints[nums.index(c)]
      # If the next place holds a larger number, this value is negative.
      try:
         nextvalue = ints[nums.index(input[i +1])]
         if nextvalue > value:
            value *= -1
      except IndexError:
         # there is no next place.
         pass
      places.append(value)
   sum = 0
   for n in places: sum += n
   # Easiest test for validity...
   if int_to_roman(sum) == input:
      return sum
   else:
      raise ValueError, 'input is not a valid roman numeral: %s' % input


def calculate_luminosity(spec_fname, distance, wavelength_column=0, wavelength_unit=u.angstrom, flux_column=1,
                         flux_unit=u.Unit('erg / (Angstrom cm2 s)')):

    #BAD STYLE change to parse quantity
    distance = u.Unit(distance)

    wavelength, flux = np.loadtxt(spec_fname, usecols=(wavelength_column, flux_column), unpack=True)

    flux_density = np.trapz(flux, wavelength) * (flux_unit * wavelength_unit)
    luminosity = (flux_density * 4 * np.pi * distance**2).to('erg/s')

    return luminosity.value, wavelength.min(), wavelength.max()

def create_synpp_yaml(self, fname, lines_db=None):
    logger.warning('Currently only works with Si and a special setup')
    if not self.atom_data.has_synpp_refs:
        raise ValueError(
            'The current atom dataset does not contain the necesarry reference files (please contact the authors)')

    self.atom_data.synpp_refs['ref_log_tau'] = -99.0
    for key, value in self.atom_data.synpp_refs.iterrows():
        try:
            tau_sobolev_idx = self.atom_data.lines_index.ix[value['line_id']]
        except KeyError:
            continue

        self.atom_data.synpp_refs['ref_log_tau'].ix[key] = np.log10(self.plasmas[0].tau_sobolevs[tau_sobolev_idx])

    relevant_synpp_refs = self.atom_data.synpp_refs[self.atom_data.synpp_refs['ref_log_tau'] > -50]

    yaml_reference = yaml.load(file(synpp_default_yaml_fname))

    if lines_db is not None:
        yaml_reference['opacity']['line_dir'] = os.path.join(lines_db, 'lines')
        yaml_reference['opacity']['line_dir'] = os.path.join(lines_db, 'refs.dat')

    yaml_reference['output']['min_wl'] = float(self.spec_angstrom.min())
    yaml_reference['output']['max_wl'] = float(self.spec_angstrom.max())


    #raise Exception("there's a problem here with units what units does synpp expect?")
    yaml_reference['opacity']['v_ref'] = float(self.tardis_config.structure.v_inner.to('cm/s').value[0] / 1e8)
    yaml_reference['grid']['v_outer_max'] = float(self.tardis_config.structure.v_outer[-1] / 1e8)

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

    yaml.dump(yaml_reference, stream=file(fname, 'w'), explicit_start=True)


def intensity_black_body(nu, T):
    """
        Calculate the intensity of a black-body according to the following formula

        .. math::
            I(\\nu, T) = \\frac{2h\\nu^3}{c^2}\frac{1}{e^{h\\nu \\beta_\\textrm{rad}} - 1}

    """
    beta_rad = 1 / (k_B_cgs * T)

    return (2 * (h_cgs * nu ** 3) / (c_cgs ** 2)) / (
        np.exp(h_cgs * nu * beta_rad) - 1)