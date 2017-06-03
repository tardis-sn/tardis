import os
import numpy as np
from astropy import constants, units as u

from tardis.io.util import to_hdf


def require_field(name):
    def _require_field(f):
        def wrapper(self, *args, **kwargs):
            if not getattr(self, name, None):
                raise AttributeError(
                        '{} is required as attribute of {} to calculate {}'.format(
                            name,
                            self.__class__.__name__,
                            f.__name__
                            )
                        )
            else:
                return f(self, *args, **kwargs)
        return wrapper
    return _require_field


class TARDISSpectrum(object):
    """
    TARDIS Spectrum object
    """

    def __init__(self, frequency, luminosity, distance=None):
        self.distance = distance

        self._frequency = frequency
        self.delta_frequency = frequency[1] - frequency[0]
        self.wavelength = self.frequency.to('angstrom', u.spectral())

        self._luminosity = luminosity

        self.luminosity_density_nu = (
                luminosity / self.delta_frequency).to('erg / (s Hz)')
        self.luminosity_density_lambda = u.Quantity(
                self.f_nu_to_f_lambda(self.luminosity_density_nu.value),
                'erg / (s Angstrom)'
                )

        if self.distance is not None:
            self._flux_nu = (
                    self.luminosity_density_nu /
                    (4 * np.pi * self.distance.to('cm')**2))

            self._flux_lambda = u.Quantity(
                    self.f_nu_to_f_lambda(self.flux_nu.value),
                    'erg / (s Angstrom cm^2)'
                    )

    @property
    def frequency(self):
        return self._frequency[:-1]

    @property
    @require_field('distance')
    def flux_nu(self):
        return self._flux_nu

    @property
    @require_field('distance')
    def flux_lambda(self):
        return self._flux_lambda

    def f_nu_to_f_lambda(self, f_nu):
        return f_nu * self.frequency.value**2 / constants.c.cgs.value / 1e8

    def plot(self, ax, mode='wavelength'):
        if mode == 'wavelength':
            ax.plot(self.wavelength.value, self.flux_lambda.value)
            ax.set_xlabel('Wavelength [%s]' % self.wavelength.unit._repr_latex_())
            ax.set_ylabel('Flux [%s]' % self.flux_lambda.unit._repr_latex_())

    def to_ascii(self, fname, mode='luminosity_density'):
        if mode == 'luminosity_density':
            np.savetxt(fname, zip(self.wavelength.value, self.luminosity_density_lambda.value))
        elif mode == 'flux':
            np.savetxt(fname, zip(self.wavelength.value, self.flux_lambda.value))
        else:
            raise NotImplementedError('only mode "luminosity_density" and "flux" are implemented')

    def to_hdf(self, path_or_buf, path='', name=''):
        """
        Store the spectrum to an HDF structure.

        Parameters
        ----------
        path_or_buf
            Path or buffer to the HDF store
        path : str
            Path inside the HDF store to store the spectrum
        name : str, optional
            A different name than 'spectrum', if needed.

        Returns
        -------
        None

        """
        if not name:
            name = 'spectrum'
        spectrum_path = os.path.join(path, name)
        properties = ['luminosity_density_nu', 'delta_frequency', 'wavelength',
                      'luminosity_density_lambda']
        to_hdf(path_or_buf, spectrum_path, {name: getattr(self, name) for name
                                            in properties})
