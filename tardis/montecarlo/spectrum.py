import os
import numpy as np
from astropy import constants, units as u

from tardis.io.util import HDFReaderWriter

class TARDISSpectrum(HDFReaderWriter, object):
    """
    TARDIS Spectrum object
    """
    hdf_properties = ['luminosity_density_nu', 'delta_frequency', 'wavelength',
                      'luminosity_density_lambda', 'frequency', 'distance']
    quantity_attrs = {'frequency': 'Hz'}

    def __init__(self, frequency, distance=None):
        self._frequency = frequency
        self.wavelength = self.frequency.to('angstrom', u.spectral())

        self.distance = distance



        self.delta_frequency = frequency[1] - frequency[0]

        self._flux_nu = np.zeros_like(frequency.value) * u.Unit('erg / (s Hz cm^2)')
        self._flux_lambda = np.zeros_like(frequency.value) * u.Unit('erg / (s Angstrom cm^2)')

        self.luminosity_density_nu = np.zeros_like(self.frequency) * u.Unit('erg / (s Hz)')
        self.luminosity_density_lambda = np.zeros_like(self.frequency) * u.Unit('erg / (s Angstrom)')

    @property
    def frequency(self):
        return self._frequency[:-1]


    @property
    def flux_nu(self):
        if self.distance is None:
            raise AttributeError('supernova distance not supplied - flux calculation impossible')
        else:
            return self._flux_nu

    @property
    def flux_lambda(self):
        if self.distance is None:
            raise AttributeError('supernova distance not supplied - flux calculation impossible')
        return self._flux_lambda

    def update_luminosity(self, spectrum_luminosity):
        self.luminosity_density_nu = (spectrum_luminosity / self.delta_frequency).to('erg / (s Hz)')
        self.luminosity_density_lambda = self.f_nu_to_f_lambda(self.luminosity_density_nu.value) \
                                         * u.Unit('erg / (s Angstrom)')

        if self.distance is not None:
            self._flux_nu = (self.luminosity_density_nu / (4 * np.pi * self.distance.to('cm')**2))

            self._flux_lambda = self.f_nu_to_f_lambda(self.flux_nu.value) * u.Unit('erg / (s Angstrom cm^2)')

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
        self.to_hdf_util(path_or_buf, spectrum_path, {name: getattr(self, name) for name
                                         in self.hdf_properties})

    @classmethod
    def from_hdf(cls, path, file_path):
        buff_path = path + '/spectrum'
        data = cls.from_hdf_util(buff_path, file_path)
        return cls(data['frequency'], data['distance'])
