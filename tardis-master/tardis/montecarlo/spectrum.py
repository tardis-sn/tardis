import warnings
import numpy as np
from astropy import units as u
from tardis.io.util import HDFWriterMixin


class TARDISSpectrum(HDFWriterMixin):
    """
    TARDISSpectrum(_frequency, luminosity)

    _frequency: astropy.units.Quantity with unit 'Hz' or a length
        These are bin edges of frequency or wavelenght bins for the spectrum.

    luminosity: astropy.units.Quantity with unit Energy per second
        The luminosity in each bin of the spectrum.

    After manually adding a distance attribute, the properties 'flux_nu' and
    'flux_lambda' become available
    """

    hdf_properties = ["_frequency", "luminosity"]

    def __init__(self, _frequency, luminosity):

        # Check for correct inputs
        if not _frequency.shape[0] == luminosity.shape[0] + 1:
            raise ValueError(
                "shape of '_frequency' and 'luminosity' are not compatible"
                ": '{}' and '{}'".format(
                    _frequency.shape[0], luminosity.shape[0]
                )
            )
        self._frequency = _frequency.to("Hz", u.spectral())
        self.luminosity = luminosity.to("erg / s")

        self.frequency = self._frequency[:-1]
        self.delta_frequency = self._frequency[1] - self._frequency[0]
        self.wavelength = self.frequency.to("angstrom", u.spectral())

        self.luminosity_density_nu = (
            self.luminosity / self.delta_frequency
        ).to("erg / (s Hz)")
        self.luminosity_density_lambda = self.f_nu_to_f_lambda(
            self.luminosity_density_nu,
        )

    @property
    def flux_nu(self):
        warnings.simplefilter("always", DeprecationWarning)
        warnings.warn(
            "TARDISSpectrum.flux_nu is deprecated, "
            "please use TARDISSpectrum.luminosity_to_flux() in the "
            "future.",
            category=DeprecationWarning,
            stacklevel=2,
        )
        warnings.simplefilter("default", DeprecationWarning)
        try:
            return self.luminosity_to_flux(
                self.luminosity_density_nu, self.distance
            )
        except AttributeError:
            raise AttributeError(
                "distance is required as attribute of"
                '{} to calculate "{}"'.format(
                    self.__class__.__name__, "flux_nu"
                )
            )

    @property
    def flux_lambda(self):
        warnings.simplefilter("always", DeprecationWarning)
        warnings.warn(
            "TARDISSpectrum.flux_lambda is deprecated, "
            "please use TARDISSpectrum.luminosity_to_flux() in the "
            "future.",
            category=DeprecationWarning,
            stacklevel=2,
        )
        warnings.simplefilter("default", DeprecationWarning)
        try:
            return self.luminosity_to_flux(
                self.luminosity_density_lambda, self.distance
            )
        except AttributeError:
            raise AttributeError(
                "distance is required as attribute of"
                '{} to calculate "{}"'.format(
                    self.__class__.__name__, "flux_lambda"
                )
            )

    @staticmethod
    def luminosity_to_flux(luminosity, distance):
        return luminosity / (4 * np.pi * distance.to("cm") ** 2)

    def f_nu_to_f_lambda(self, f_nu):
        return f_nu * self.frequency / self.wavelength

    def plot(self, ax=None, mode="wavelength"):
        if ax is None:
            from matplotlib.pyplot import gca

            ax = gca()
        if mode == "wavelength":
            ax.plot(self.wavelength.value, self.luminosity_density_lambda.value)
            ax.set_xlabel(
                "Wavelength [{}]".format(self.wavelength.unit._repr_latex_())
            )
            ax.set_ylabel(
                "Flux [{:s}]".format(
                    self.luminosity_density_lambda.unit._repr_latex_()
                )
            )
        else:
            warnings.warn(
                "Did not find plotting mode {}, doing nothing.".format(mode)
            )

    def to_ascii(self, fname, mode="luminosity_density"):
        if mode == "luminosity_density":
            np.savetxt(
                fname,
                list(
                    zip(
                        self.wavelength.value,
                        self.luminosity_density_lambda.value,
                    )
                ),
            )
        elif mode == "flux":
            np.savetxt(
                fname, zip(self.wavelength.value, self.flux_lambda.value)
            )
        else:
            raise NotImplementedError(
                'only mode "luminosity_density"' 'and "flux" are implemented'
            )
