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

def create_synpp_yaml(radial1d_mdl, fname, shell_no=0, lines_db=None):
    logger.warning('Currently only works with Si and a special setup')
    if not radial1d_mdl.atom_data.has_synpp_refs:
        raise ValueError(
            'The current atom dataset does not contain the necesarry reference files (please contact the authors)')

    radial1d_mdl.atom_data.synpp_refs['ref_log_tau'] = -99.0
    for key, value in radial1d_mdl.atom_data.synpp_refs.iterrows():
        try:
            radial1d_mdl.atom_data.synpp_refs['ref_log_tau'].ix[key] = np.log10(
                radial1d_mdl.plasma_array.tau_sobolevs[0].ix[value['line_id']])
        except KeyError:
            pass


    relevant_synpp_refs = radial1d_mdl.atom_data.synpp_refs[radial1d_mdl.atom_data.synpp_refs['ref_log_tau'] > -50]

    yaml_reference = yaml.load(file(synpp_default_yaml_fname))

    if lines_db is not None:
        yaml_reference['opacity']['line_dir'] = os.path.join(lines_db, 'lines')
        yaml_reference['opacity']['line_dir'] = os.path.join(lines_db, 'refs.dat')

    yaml_reference['output']['min_wl'] = float(radial1d_mdl.spectrum.wavelength.to('angstrom').value.min())
    yaml_reference['output']['max_wl'] = float(radial1d_mdl.spectrum.wavelength.to('angstrom').value.max())


    #raise Exception("there's a problem here with units what units does synpp expect?")
    yaml_reference['opacity']['v_ref'] = float((radial1d_mdl.tardis_config.structure.v_inner[0].to('km/s') /
                                               (1000. * u.km / u.s)).value)
    yaml_reference['grid']['v_outer_max'] = float((radial1d_mdl.tardis_config.structure.v_outer[-1].to('km/s') /
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

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
