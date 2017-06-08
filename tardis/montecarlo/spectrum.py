import warnings
import numpy as np
import os
from astropy import units as u

from tardis.io.util import HDFReaderWriter

class TARDISSpectrum(HDFReaderWriter, object):
    """
    TARDISSpectrum(_frequency, luminosity)

    _frequency: astropy.units.Quantity with unit 'Hz' or a length
        These are bin edges of frequency or wavelenght bins for the spectrum.

    luminosity: astropy.units.Quantity with unit Energy per second
        The luminosity in each bin of the spectrum.

    After manually adding a distance attribute, the properties 'flux_nu' and
    'flux_lambda' become available
    """
    hdf_properties = ['_frequency', 'luminosity']
    quantity_attrs = {'_frequency': 'Hz', 'luminosity': 'erg/s'}

    def __init__(self, _frequency, luminosity):

        # Check for correct inputs
        if not _frequency.shape[0] == luminosity.shape[0] + 1:
            raise ValueError(
                    "shape of '_frequency' and 'luminosity' are not compatible"
                    ": '{}' and '{}'".format(
                        _frequency.shape[0],
                        luminosity.shape[0])
                    )
        self._frequency = _frequency.to('Hz', u.spectral())
        self.luminosity = luminosity.to('erg / s')

        self.frequency = self._frequency[:-1]
        self.delta_frequency = self._frequency[1] - self._frequency[0]
        self.wavelength = self.frequency.to('angstrom', u.spectral())

        self.luminosity_density_nu = (
                self.luminosity / self.delta_frequency).to('erg / (s Hz)')
        self.luminosity_density_lambda = self.f_nu_to_f_lambda(
                self.luminosity_density_nu,
                )

    @property
    def flux_nu(self):
        warnings.simplefilter('always', DeprecationWarning)
        warnings.warn(
                "TARDISSpectrum.flux_nu is deprecated, "
                "please use TARDISSpectrum.luminosity_to_flux() in the "
                "future.",
                category=DeprecationWarning, stacklevel=2)
        warnings.simplefilter('default', DeprecationWarning)
        try:
            return self.luminosity_to_flux(
                    self.luminosity_density_nu,
                    self.distance)
        except AttributeError:
                raise AttributeError(
                        'distance is required as attribute of'
                        '{} to calculate "{}"'.format(
                            self.__class__.__name__,
                            'flux_nu'
                            )
                        )

    @property
    def flux_lambda(self):
        warnings.simplefilter('always', DeprecationWarning)
        warnings.warn(
                "TARDISSpectrum.flux_lambda is deprecated, "
                "please use TARDISSpectrum.luminosity_to_flux() in the "
                "future.",
                category=DeprecationWarning, stacklevel=2)
        warnings.simplefilter('default', DeprecationWarning)
        try:
            return self.luminosity_to_flux(
                    self.luminosity_density_lambda,
                    self.distance
                    )
        except AttributeError:
                raise AttributeError(
                        'distance is required as attribute of'
                        '{} to calculate "{}"'.format(
                            self.__class__.__name__,
                            'flux_lambda'
                            )
                        )

    @staticmethod
    def luminosity_to_flux(luminosity, distance):
        return luminosity / (4 * np.pi * distance.to('cm')**2)

    def f_nu_to_f_lambda(self, f_nu):
        return f_nu * self.frequency / self.wavelength

    def plot(self, ax, mode='wavelength'):
        if mode == 'wavelength':
            ax.plot(
                    self.wavelength.value,
                    self.luminosity_density_lambda.value)
            ax.set_xlabel('Wavelength [{}]'.format(
                self.wavelength.unit._repr_latex_())
                )
            ax.set_ylabel(
                    'Flux [{:s}]'.format(
                        self.luminosity_density_lambda.unit._repr_latex_())
                    )

    def to_ascii(self, fname, mode='luminosity_density'):
        if mode == 'luminosity_density':
            np.savetxt(
                    fname, zip(
                        self.wavelength.value,
                        self.luminosity_density_lambda.value))
        elif mode == 'flux':
            np.savetxt(
                    fname,
                    zip(self.wavelength.value, self.flux_lambda.value))
        else:
            raise NotImplementedError(
                    'only mode "luminosity_density"'
                    'and "flux" are implemented')

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
        self.to_hdf_util(path_or_buf, spectrum_path, {name: getattr(self, name) for name
                                         in self.hdf_properties})

    @classmethod
    def from_hdf(cls, file_path, path=''):
        buff_path = os.path.join(path, 'spectrum')
        data = cls.from_hdf_util(file_path, buff_path)
        return cls(data['_frequency'], data['luminosity'])