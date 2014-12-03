******
Models
******

The models in TARDIS encode the general information about the setup of the
star such as geometry (currently 1d radial with homologous expansion is the
only option), densities, abundances, etc. Here is a simple example of a radial
one-dimensional model (if no units are given the code **ALWAYS** assumes cgs)::

    >>> from tardis.models import Radial1D
    >>> from astropy import units as u
    >>> import numpy as np
    >>> my_model = Radial1D(radius=np.linspace(8000, 10000, 10) * u.km)
    >>> my_model.radius
        <Quantity [  8.00000000e+08,  8.22222222e+08,  8.44444444e+08,
                     8.66666667e+08,  8.88888889e+08,  9.11111111e+08,
                     9.33333333e+08,  9.55555556e+08,  9.77777778e+08,
                     1.00000000e+09] cm>

