******
Models
******

The models in TARDIS encode the general information about the setup of the
star such as geometry (currently 1d radial with homologous expansion is the
only option), densities, abundances, etc. Here is a simple example of a radial
one-dimensional model (if no units are given the code **ALWAYS** assumes cgs)::

    >>> from tardis.models import HomologousRadial1D
    >>> from astropy import units as u, constants as const
    >>> import numpy as np
    >>> velocity=np.linspace(8000, 10000, 10) * u.km / u.s
    >>> time = 10 * u.day
    >>> time0 = 2 * u.s # time of the model
    >>> density0 = np.ones(10) * 1.4 * u.Msun / (0.6 * const.R_earth)**3  #density at the time of the model
    >>> my_model = HomologousRadial1D(velocity, time, time0, density0)
    >>> my_model.velocity
        <Quantity [  8.00000000e+08,  8.22222222e+08,  8.44444444e+08,
                     8.66666667e+08,  8.88888889e+08,  9.11111111e+08,
                     9.33333333e+08,  9.55555556e+08,  9.77777778e+08,
                     1.00000000e+09] cm / s>
    >>> my_model.density0 # the given density converted to cgs system
    <Quantity [ 10732560.04461544, 10732560.04461544, 10732560.04461544,
                10732560.04461544, 10732560.04461544, 10732560.04461544,
                10732560.04461544, 10732560.04461544, 10732560.04461544,
                10732560.04461544] g / cm3>
    >>> my_model.density # the density converted to from time0 to time [(time0/time)**3]
    <Quantity [  1.33122691e-10,  1.33122691e-10,  1.33122691e-10,
             1.33122691e-10,  1.33122691e-10,  1.33122691e-10,
             1.33122691e-10,  1.33122691e-10,  1.33122691e-10,
             1.33122691e-10] g / cm3>


.. automodapi:: tardis.models