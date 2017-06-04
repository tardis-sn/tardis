import numpy as np

from astropy import units as u


def require_field(name):
    def _require_field(f):
        def wrapper(self, *args, **kwargs):
            if not getattr(self, name, None):
                raise AttributeError(
                        '{} is required as attribute of'
                        '{} to calculate {}'.format(
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
    TARDISSpectrum(_frequency, luminosity)

    _frequency: astropy.units.Quantity with unit 'Hz' or a length
        These are bin edges of frequency or wavelenght bins for the spectrum.

    luminosity: astropy.units.Quantity with unit Energy per second
        The luminosity in each bin of the spectrum.

    After manually adding a distance attribute, the properties 'flux_nu' and
    'flux_lambda' become available
    """

    def __init__(self, _frequency, luminosity):
        if not _frequency.shape[0] == luminosity.shape[0] + 1:
            raise ValueError(
                    "shape of '_frequency' and 'luminosity' are not compatible"
                    ": '{}' and '{}'".format(
                        _frequency.shape[0],
                        luminosity.shape[0])
                    )
        self._frequency = _frequency.to('Hz', u.spectral())
        self.luminosity = luminosity.to('erg / s')

    @property
    def frequency(self):
        return self._frequency[:-1]

    @property
    def delta_frequency(self):
        return self.frequency[1] - self.frequency[0]

    @property
    def wavelength(self):
        return self.frequency.to('angstrom', u.spectral())

    @property
    def luminosity_density_nu(self):
        return (self.luminosity / self.delta_frequency).to('erg / (s Hz)')

    @property
    def luminosity_density_lambda(self):
        return self.f_nu_to_f_lambda(
                self.luminosity_density_nu,
                )

    @property
    @require_field('distance')
    def flux_nu(self):
        return self.luminosity_to_flux(
                self.luminosity_density_nu,
                self.distance)

    @property
    @require_field('distance')
    def flux_lambda(self):
        return self.luminosity_to_flux(
                self.luminosity_density_lambda,
                self.distance
                )

    @staticmethod
    def luminosity_to_flux(luminosity, distance):
        return luminosity / (4 * np.pi * distance.to('cm')**2)

    def f_nu_to_f_lambda(self, f_nu):
        return f_nu * self.frequency / self.wavelength

    def plot(self, ax, mode='wavelength'):
        if mode == 'wavelength':
            ax.plot(self.wavelength.value, self.flux_lambda.value)
            ax.set_xlabel('Wavelength [%s]' % self.wavelength.unit._repr_latex_())
            ax.set_ylabel('Flux [%s]' % self.flux_lambda.unit._repr_latex_())

    def to_ascii(self, fname, mode='luminosity_density'):
        if mode == 'luminosity_density':
            np.savetxt(fname, zip(self.wavelength.value, self.luminosity_density_lambda.value))
        elif mode == 'flux':
            np.savetxt(
                    fname,
                    zip(self.wavelength.value, self.flux_lambda.value))
        else:
            raise NotImplementedError('only mode "luminosity_density" and "flux" are implemented')

    def to_hdf(self, path_or_buf, path):
        pass

    @classmethod
    def from_hdf(cls, path_or_buf, path):
        pass
